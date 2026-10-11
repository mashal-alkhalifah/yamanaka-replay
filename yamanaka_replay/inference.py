"""
Exact Bayesian inference over candidate recipes, and expected information gain.

The hypothesis space is every 4-subset of the 24 candidates -- C(24,4) = 10,626
hypotheses -- which is small enough to carry an *exact* posterior. No sampling,
no variational approximation.

The key computational observation
---------------------------------
For a tested pool ``S``, the colony-count distribution depends on a hypothesis
``h`` only through the scalar ``|S n h|``, which takes one of five values
(0..4). So both the likelihood and the expected information gain collapse to a
five-way mixture, regardless of how many hypotheses are alive. Choosing the
next experiment therefore costs the same whether the posterior spans ten
hypotheses or ten thousand.

This is what makes the adaptive designer practical to scale to larger
candidate sets later.
"""

from __future__ import annotations

from itertools import combinations
from math import lgamma

import numpy as np

from .oracle import ColonyOracle, N_CANDIDATES, RECIPE_SIZE

# ---------------------------------------------------------------- hypotheses

HYPOTHESES: list[frozenset[int]] = [
    frozenset(c) for c in combinations(range(N_CANDIDATES), RECIPE_SIZE)
]

#: Boolean (n_hypotheses, n_candidates) membership matrix.
HYPOTHESIS_MASK = np.zeros((len(HYPOTHESES), N_CANDIDATES), dtype=bool)
for _i, _h in enumerate(HYPOTHESES):
    for _g in _h:
        HYPOTHESIS_MASK[_i, _g] = True

# ------------------------------------------------------------ count support

_MAX_COUNT = 90  # colony counts above this are clipped; the rate never approaches it
_COUNTS = np.arange(_MAX_COUNT + 1)
_LOG_FACTORIAL = np.array([lgamma(y + 1) for y in _COUNTS])


def poisson_pmf(rate: float) -> np.ndarray:
    """Poisson pmf over ``0.._MAX_COUNT``, renormalised after truncation."""
    if rate <= 1e-12:
        pmf = np.zeros(_MAX_COUNT + 1)
        pmf[0] = 1.0
        return pmf
    log_pmf = _COUNTS * np.log(rate) - rate - _LOG_FACTORIAL
    pmf = np.exp(log_pmf)
    total = pmf.sum()
    return pmf / total if total > 0 else pmf


def expected_rate(overlap: int, pool_size: int) -> float:
    """The searcher's *assumed* colony rate for a given overlap and pool size.

    This mirrors :meth:`ColonyOracle.rate` but is the model the searcher
    believes. When the two diverge, the searcher is misspecified -- which is
    precisely what ``03_misspecification.py`` probes.
    """
    if overlap == RECIPE_SIZE:
        efficiency = 1.0
    elif overlap == RECIPE_SIZE - 1:
        efficiency = ColonyOracle.PARTIAL_3
    else:
        efficiency = 0.0
    penalty = np.exp(-ColonyOracle.DECAY * max(0, pool_size - RECIPE_SIZE))
    return ColonyOracle.BASE * efficiency * penalty


def _overlaps(pool) -> np.ndarray:
    """``|pool n h|`` for every hypothesis ``h``."""
    mask = np.zeros(N_CANDIDATES, dtype=bool)
    mask[list(pool)] = True
    return HYPOTHESIS_MASK[:, mask].sum(axis=1)


def _entropy_bits(p: np.ndarray) -> float:
    p = p[p > 1e-15]
    return float(-(p * np.log2(p)).sum())


def expected_information_gain(posterior: np.ndarray, pool) -> float:
    """
    Mutual information in bits between the unknown recipe and the outcome of
    testing ``pool``, under the current ``posterior``.

        EIG(S) = H[ p(y | S) ] - E_h[ H[ p(y | S, h) ] ]

    Computed over the five overlap classes rather than all hypotheses, so the
    cost is independent of posterior size.
    """
    pool_size = len(pool)
    overlap = _overlaps(pool)
    class_prob = np.bincount(overlap, weights=posterior, minlength=RECIPE_SIZE + 1)
    class_prob = class_prob[: RECIPE_SIZE + 1]
    if class_prob.sum() <= 0:
        return 0.0

    per_class_pmf = np.vstack(
        [poisson_pmf(expected_rate(k, pool_size)) for k in range(RECIPE_SIZE + 1)]
    )
    marginal = class_prob @ per_class_pmf
    conditional = sum(
        class_prob[k] * _entropy_bits(per_class_pmf[k]) for k in range(RECIPE_SIZE + 1)
    )
    return _entropy_bits(marginal) - float(conditional)


def bayes_update(posterior: np.ndarray, pool, observed: int) -> np.ndarray:
    """Posterior after observing ``observed`` colonies from ``pool``."""
    pool_size = len(pool)
    overlap = _overlaps(pool)
    likelihood_by_class = np.array(
        [
            poisson_pmf(expected_rate(k, pool_size))[min(observed, _MAX_COUNT)]
            for k in range(RECIPE_SIZE + 1)
        ]
    )
    updated = posterior * likelihood_by_class[overlap]
    total = updated.sum()
    if total <= 0:  # observation impossible under every hypothesis -> reset
        return np.full_like(posterior, 1.0 / len(posterior))
    return updated / total


def uniform_prior() -> np.ndarray:
    """A flat prior: no knowledge of which factors matter."""
    return np.full(len(HYPOTHESES), 1.0 / len(HYPOTHESES))


def weighted_prior(factor_weights: np.ndarray) -> np.ndarray:
    """
    Build a hypothesis prior from per-factor plausibility weights.

    A hypothesis' weight is the product of its four factors' weights. This is
    how a literature-derived prior -- "these genes are reported pluripotency
    regulators" -- enters the search. Note it is a *weighting*, not a
    constraint: a factor the literature missed still carries nonzero mass and
    can be recovered from data, which is why the method degrades gracefully
    under a wrong prior (see ``02_prior_robustness.py``).
    """
    prior = np.array(
        [float(np.prod([factor_weights[g] for g in h])) for h in HYPOTHESES]
    )
    return prior / prior.sum()


def marginal_factor_probability(posterior: np.ndarray) -> np.ndarray:
    """P(factor is in the recipe) for each candidate, under ``posterior``."""
    return HYPOTHESIS_MASK.T @ posterior
