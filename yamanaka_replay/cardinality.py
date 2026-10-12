"""
Inference when the *size* of the recipe is unknown.

Everything in ``inference.py`` assumes the recipe has exactly four factors,
because Yamanaka's did. Experiment 3 showed what that assumption costs when it
is wrong: searching for a four-factor recipe when the truth has five fails
**silently**, converging confidently on a wrong answer. Experiment 6 then
showed that this is the only genuinely structural failure mode -- the others
are calibration error.

A limitation that is the single most dangerous property of a method is worth
removing rather than documenting, which is what this module does. Recipe size
becomes a latent variable: the hypothesis space is the union of all subsets of
every plausible size, and the posterior over sizes is just a marginal of the
posterior over subsets.

The computational trick that makes the fixed-size version cheap still applies,
one dimension larger. The colony rate for hypothesis ``h`` under pool ``S``
depends on ``h`` only through the pair ``(|h|, |S n h|)``, so likelihood and
expected information gain collapse onto that small grid of classes rather than
over all hypotheses. Allowing sizes 3-6 over 24 candidates gives 189,750
hypotheses but only 22 reachable classes.

The price is paid in a different currency: with sizes 3-6 allowed, a pool that
yields colonies no longer pins the recipe size, so more experiments are needed
to separate "these four" from "these four plus a fifth that happened to be in
every positive pool". Experiment 10 measures that price.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np

from .inference import _MAX_COUNT, poisson_pmf
from .oracle import ColonyOracle


class HypothesisSpace:
    """
    All subsets of ``n_candidates`` whose size lies in ``sizes``.

    Attributes
    ----------
    hypotheses:
        The subsets, as frozensets, in a fixed order.
    mask:
        Boolean ``(n_hypotheses, n_candidates)`` membership matrix.
    size:
        Integer array of ``|h|`` per hypothesis.
    classes:
        For each ``(size, overlap)`` pair that can occur, an index into a flat
        class axis. Both inference routines work on that axis.
    """

    def __init__(self, n_candidates: int, sizes=(3, 4, 5, 6)):
        self.n_candidates = n_candidates
        self.sizes = tuple(sorted(sizes))
        subsets: list[frozenset[int]] = []
        for k in self.sizes:
            subsets.extend(frozenset(c) for c in combinations(range(n_candidates), k))
        self.hypotheses = subsets
        self.mask = np.zeros((len(subsets), n_candidates), dtype=bool)
        for i, h in enumerate(subsets):
            self.mask[i, list(h)] = True
        self.size = self.mask.sum(axis=1).astype(np.int16)

        self.max_size = max(self.sizes)
        # Flat class axis: class(size, overlap) = size * (max_size + 1) + overlap
        self._stride = self.max_size + 1
        self.n_classes = (self.max_size + 1) * self._stride

    def __len__(self) -> int:
        return len(self.hypotheses)

    def overlaps(self, pool) -> np.ndarray:
        """``|pool n h|`` for every hypothesis."""
        sel = np.zeros(self.n_candidates, dtype=bool)
        sel[list(pool)] = True
        return self.mask[:, sel].sum(axis=1).astype(np.int16)

    def class_index(self, pool) -> np.ndarray:
        """The ``(|h|, |pool n h|)`` class of every hypothesis, flattened."""
        return self.size * self._stride + self.overlaps(pool)

    def decode_class(self, index: int) -> tuple[int, int]:
        return divmod(int(index), self._stride)

    def uniform_prior(self) -> np.ndarray:
        """
        Flat over hypotheses, which is *not* flat over sizes: there are far more
        6-subsets than 3-subsets, so a uniform prior over subsets already
        believes the recipe is large. Kept deliberately, as the honest
        no-knowledge default; ``size_prior`` offers the alternative.
        """
        return np.full(len(self), 1.0 / len(self))

    def size_prior(self, weights: dict[int, float] | None = None) -> np.ndarray:
        """
        A prior that is uniform *within* each size and spreads mass *across*
        sizes according to ``weights`` (equal by default).

        This is usually what an experimenter means by "I have no idea how many
        factors there are": each plausible count is equally likely, not each
        subset.
        """
        weights = weights or {k: 1.0 for k in self.sizes}
        prior = np.zeros(len(self))
        for k in self.sizes:
            rows = self.size == k
            prior[rows] = weights.get(k, 0.0) / max(1, int(rows.sum()))
        total = prior.sum()
        if total <= 0:
            return self.uniform_prior()
        return prior / total

    def marginal_size(self, posterior: np.ndarray) -> dict[int, float]:
        """P(recipe size = k) under ``posterior``."""
        return {k: float(posterior[self.size == k].sum()) for k in self.sizes}

    def marginal_factor(self, posterior: np.ndarray) -> np.ndarray:
        """P(factor in recipe) per candidate."""
        return self.mask.T @ posterior


def expected_rate(h_size: int, overlap: int, pool_size: int) -> float:
    """
    The searcher's assumed colony rate, generalised to any recipe size.

    Same shape as the fixed-size model: full efficiency when the pool contains
    the whole recipe, a partial effect when it is one factor short, nothing
    otherwise, times the pool-size penalty.
    """
    if overlap >= h_size:
        efficiency = 1.0
    elif overlap == h_size - 1:
        efficiency = ColonyOracle.PARTIAL_3
    else:
        efficiency = 0.0
    penalty = np.exp(-ColonyOracle.DECAY * max(0, pool_size - h_size))
    return ColonyOracle.BASE * efficiency * penalty


def _class_pmfs(space: HypothesisSpace, pool_size: int) -> np.ndarray:
    """Count distribution for every ``(size, overlap)`` class."""
    pmfs = np.zeros((space.n_classes, _MAX_COUNT + 1))
    for k in space.sizes:
        for overlap in range(k + 1):
            idx = k * space._stride + overlap
            pmfs[idx] = poisson_pmf(expected_rate(k, overlap, pool_size))
    return pmfs


def _entropy_bits(p: np.ndarray) -> float:
    p = p[p > 1e-15]
    return float(-(p * np.log2(p)).sum())


def expected_information_gain(space: HypothesisSpace, posterior: np.ndarray, pool) -> float:
    """Mutual information in bits between the recipe and the outcome of ``pool``."""
    classes = space.class_index(pool)
    class_prob = np.bincount(classes, weights=posterior, minlength=space.n_classes)
    live = class_prob > 0
    if not live.any():
        return 0.0
    pmfs = _class_pmfs(space, len(pool))[live]
    probs = class_prob[live]
    marginal = probs @ pmfs
    conditional = float(sum(p * _entropy_bits(row) for p, row in zip(probs, pmfs)))
    return _entropy_bits(marginal) - conditional


def bayes_update(space: HypothesisSpace, posterior: np.ndarray, pool, observed: int) -> np.ndarray:
    """Posterior after seeing ``observed`` colonies from ``pool``."""
    classes = space.class_index(pool)
    column = _class_pmfs(space, len(pool))[:, min(observed, _MAX_COUNT)]
    updated = posterior * column[classes]
    total = updated.sum()
    if total <= 0:
        return space.uniform_prior()
    return updated / total


def adaptive_design_mixed(
    oracle: ColonyOracle,
    space: HypothesisSpace,
    seed: int = 0,
    max_experiments: int = 80,
    confidence: float = 0.95,
    prior: np.ndarray | None = None,
    n_random_pools: int = 200,
) -> tuple[set[int], int, float, dict[int, float]]:
    """
    Closed-loop design with recipe size treated as unknown.

    Returns ``(recipe, n_experiments, final_confidence, posterior over sizes)``.
    The size posterior is the point of the exercise: it is what tells the
    experimenter that their assumed cardinality was wrong, which the fixed-size
    designer can never do.
    """
    rng = np.random.default_rng(seed)
    posterior = space.size_prior() if prior is None else prior.copy()

    while oracle.n_experiments < max_experiments and posterior.max() < confidence:
        pools = [
            tuple(sorted(rng.choice(space.n_candidates,
                                    size=int(rng.integers(2, 13)),
                                    replace=False).tolist()))
            for _ in range(n_random_pools)
        ]
        ranked = np.argsort(-space.marginal_factor(posterior))
        for size in (3, 4, 5, 6, 7, 8):
            pools.append(tuple(sorted(int(x) for x in ranked[:size])))
            pools.append(tuple(sorted(int(x) for x in ranked[1 : size + 1])))
        pools = list(dict.fromkeys(pools))

        pool = max(pools, key=lambda p: expected_information_gain(space, posterior, p))
        observed = oracle.run(pool)
        posterior = bayes_update(space, posterior, pool, observed)

    best = space.hypotheses[int(np.argmax(posterior))]
    return set(best), oracle.n_experiments, float(posterior.max()), space.marginal_size(posterior)


class GeneralOracle(ColonyOracle):
    """
    A colony oracle whose hidden recipe may be of any size.

    The fixed-size :class:`ColonyOracle` hard-codes four, and
    :class:`MisspecifiedOracle` fakes a fifth factor by gating on it, which
    entangles recipe size with the partial-effect structure. This states it
    directly: full efficiency when the pool contains the whole recipe, the
    partial effect one factor short, nothing otherwise.
    """

    def rate(self, pool) -> float:
        pool = set(pool)
        if not pool:
            return 0.0
        hits = len(pool & self.truth)
        k = len(self.truth)
        if hits >= k:
            efficiency = 1.0
        elif hits == k - 1:
            efficiency = self.PARTIAL_3
        else:
            efficiency = 0.0
        penalty = np.exp(-self.DECAY * max(0, len(pool) - k))
        return self.BASE * efficiency * penalty
