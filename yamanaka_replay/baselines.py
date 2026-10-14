"""
Baselines from the group-testing literature.

The main study compares an adaptive Bayesian designer against the 2006
leave-one-out protocol and a random-pool ablation. That is not a fair test of
novelty, because the problem has a precise and well-studied home:

    A pool is negative when it contains at most ``l`` essential factors,
    positive when it contains at least ``u``, and stochastic in between.

That is **threshold group testing** (Damaschke 2006), in its stochastic form
(Cheraghchi et al. 2013). The colony assay is the case ``u = |R|``,
``l = |R| - 2`` and gap ``g = u - l - 1 = 1``, with the gap response being the
graded three-of-four signal. Known bounds for that model are ``O(g^2 d log n)``
tests, which at ``g = 1``, ``d = 4``, ``n = 24`` is about 13 -- the same order
as the 10 the Bayesian designer achieves. So the right question is not "does
adaptive design beat the 2006 protocol" (it does, and theory says it must) but
"does the Bayesian machinery earn its keep against a classical algorithm".

This module provides two such comparators.

``half_removal_search``
    The natural algorithm for a conjunctive (AND) readout, and the structural
    analogue of binary splitting. Remove half the remaining candidates and
    retest: if the pool is still positive, every removed factor was
    inessential and the candidate set halves; if it goes negative, the removed
    half contains an essential factor, found by binary search. Costs
    ``O(k log(n/k))`` tests and carries **no model of the assay beyond a
    positivity threshold** -- in particular, it never assumes how many factors
    the recipe has, so the cardinality failure mode of the fixed-size Bayesian
    designer cannot arise.

``two_stage_cardinality``
    The standard group-testing recipe for unknown ``d``: estimate the number of
    defectives first, then run a fixed-size search at the estimate. This is the
    comparator for ``cardinality.adaptive_design_mixed``, which instead carries
    size and identity in one joint posterior. Without this baseline, the claim
    that joint inference helps is untested.

Fairness note on the threshold
------------------------------
``half_removal_search`` needs to call a count positive or negative. It derives
that threshold from the same assumed constants the Bayesian designer uses
(``BASE``, ``DECAY``, ``PARTIAL_3``), as the geometric mean of the assumed full
and partial rates at the pool size being tested. Both methods are therefore
exposed to the *same* calibration error, which is what makes the comparison
under misspecification meaningful rather than rigged in either direction.
"""

from __future__ import annotations

import numpy as np

from .inference import (
    HYPOTHESES,
    bayes_update,
    expected_information_gain,
    marginal_factor_probability,
    uniform_prior,
)
from .oracle import ColonyOracle


def positivity_threshold(pool_size: int, recipe_size: int = 4) -> float:
    """
    The count above which a pool is called positive.

    Placed at the geometric mean of the assumed full-recipe rate and the
    assumed one-short rate, so that it sits midway between them on the scale
    the Poisson noise actually lives on. Scales with pool size because the
    pool-size penalty lowers both rates together: at 60 candidates a full
    recipe yields about 6 colonies, not 45, and a fixed threshold of 12 would
    call every pool negative.
    """
    penalty = float(np.exp(-ColonyOracle.DECAY * max(0, pool_size - recipe_size)))
    full = ColonyOracle.BASE * penalty
    partial = ColonyOracle.BASE * ColonyOracle.PARTIAL_3 * penalty
    return float(np.sqrt(max(full, 1e-9) * max(partial, 1e-9)))


def _is_positive(oracle: ColonyOracle, pool, recipe_size: int = 4) -> bool:
    pool = sorted(set(pool))
    count = oracle.run(pool)
    return count >= positivity_threshold(len(pool), recipe_size)


def half_removal_search(
    oracle: ColonyOracle,
    n_candidates: int = 24,
    max_experiments: int = 150,
    seed: int = 0,
    recipe_size_hint: int = 4,
) -> tuple[set[int], int]:
    """
    Model-free conjunctive search: shrink a known-positive pool to its core.

    Returns ``(recovered set, experiments used)``. The returned set is whatever
    survived elimination, of whatever size -- the algorithm is told nothing
    about how many factors to expect. ``recipe_size_hint`` only enters the
    positivity threshold, never the search.

    The procedure:

    1.  Confirm the full candidate pool is positive. If it is not -- which
        happens once the pool-size penalty drags a full recipe's yield down to
        the noise floor -- there is nothing to shrink and the method reports
        failure rather than guessing. That ceiling is a real property of the
        assay, not of the algorithm.
    2.  Repeatedly remove half of the not-yet-confirmed factors and retest.
        Still positive: discard that half outright. Gone negative: it holds at
        least one essential factor, located by binary search within it.
    3.  Stop when every surviving factor has been confirmed essential.
    """
    rng = np.random.default_rng(seed)
    live = list(range(n_candidates))
    rng.shuffle(live)

    if not _is_positive(oracle, live, recipe_size_hint):
        # No positive control: the search has no foothold.
        return set(), oracle.n_experiments

    confirmed: set[int] = set()

    def still_positive(candidate_pool) -> bool:
        return _is_positive(oracle, candidate_pool, recipe_size_hint)

    while oracle.n_experiments < max_experiments:
        undecided = [x for x in live if x not in confirmed]
        if not undecided:
            break

        half = undecided[: max(1, len(undecided) // 2)]
        reduced = [x for x in live if x not in set(half)]

        if reduced and still_positive(reduced):
            live = reduced  # every factor in `half` was inessential
            continue

        # `half` contains at least one essential factor: binary search for one.
        block = list(half)
        current = list(live)
        while len(block) > 1 and oracle.n_experiments < max_experiments:
            first = block[: len(block) // 2]
            trial = [x for x in current if x not in set(first)]
            if trial and still_positive(trial):
                current = trial  # `first` inessential -- drop it for good
                block = block[len(block) // 2 :]
            else:
                block = first  # the essential factor is in `first`
        live = current
        confirmed.add(block[0])

    return set(live), oracle.n_experiments


def two_stage_cardinality(
    oracle: ColonyOracle,
    space,
    seed: int = 0,
    stage1_experiments: int = 6,
    max_experiments: int = 80,
    confidence: float = 0.95,
    n_random_pools: int = 200,
) -> tuple[set[int], int, float, int]:
    """
    The classical two-stage treatment of unknown cardinality.

    Stage 1 spends a fixed budget of information-gain experiments on the
    mixed-size hypothesis space purely to estimate the recipe size. Stage 2
    then **commits** to the most probable size and runs the ordinary
    fixed-size search within it.

    This is the baseline the joint latent-size designer has to beat. If
    committing early is as good as carrying the size posterior throughout,
    there is no case for the joint treatment.

    Returns ``(recipe, experiments, confidence, estimated size)``.
    """
    from .cardinality import adaptive_design_mixed, bayes_update as mixed_update
    from .cardinality import expected_information_gain as mixed_eig

    rng = np.random.default_rng(seed)
    posterior = space.size_prior()

    # ---- stage 1: estimate the size, do not try to finish the search
    while oracle.n_experiments < stage1_experiments:
        pools = [
            tuple(sorted(rng.choice(space.n_candidates,
                                    size=int(rng.integers(2, 13)),
                                    replace=False).tolist()))
            for _ in range(n_random_pools)
        ]
        ranked = np.argsort(-space.marginal_factor(posterior))
        for size in (3, 4, 5, 6):
            pools.append(tuple(sorted(int(x) for x in ranked[:size])))
        pools = list(dict.fromkeys(pools))
        pool = max(pools, key=lambda p: mixed_eig(space, posterior, p))
        posterior = mixed_update(space, posterior, pool, oracle.run(pool))

    size_posterior = space.marginal_size(posterior)
    estimated = max(size_posterior, key=size_posterior.get)

    # ---- stage 2: commit to that size and search within it
    restricted = np.where(space.size == estimated, posterior, 0.0)
    total = restricted.sum()
    restricted = (
        restricted / total
        if total > 0
        else np.where(space.size == estimated, 1.0, 0.0)
          / max(1, int((space.size == estimated).sum()))
    )

    while oracle.n_experiments < max_experiments and restricted.max() < confidence:
        pools = [
            tuple(sorted(rng.choice(space.n_candidates,
                                    size=int(rng.integers(2, 13)),
                                    replace=False).tolist()))
            for _ in range(n_random_pools)
        ]
        ranked = np.argsort(-space.marginal_factor(restricted))
        for size in (3, 4, 5, 6, 8):
            pools.append(tuple(sorted(int(x) for x in ranked[:size])))
            pools.append(tuple(sorted(int(x) for x in ranked[1 : size + 1])))
        pools = list(dict.fromkeys(pools))
        pool = max(pools, key=lambda p: mixed_eig(space, restricted, p))
        observed = oracle.run(pool)
        restricted = mixed_update(space, restricted, pool, observed)
        # Stage 2 cannot revisit the size decision: mass is held inside it.
        restricted = np.where(space.size == estimated, restricted, 0.0)
        s = restricted.sum()
        if s <= 0:  # the committed size became impossible -- stage 2 is stuck
            break
        restricted /= s

    best = space.hypotheses[int(np.argmax(restricted))]
    return set(best), oracle.n_experiments, float(restricted.max()), int(estimated)


def theoretical_bound(n_candidates: int, recipe_size: int = 4, gap: int = 1) -> float:
    """
    The ``O(g^2 d log n)`` scaling for stochastic threshold group testing,
    with the constant dropped -- useful only as a growth-rate reference curve,
    not as a prediction of absolute test counts.
    """
    return float(gap**2 * recipe_size * np.log(n_candidates))
