"""
Search strategies, all spending the same currency: wet-lab experiments.

Four strategies are compared:

``historical_protocol``
    A faithful replay of the 2006 procedure. Fixed, non-adaptive, 61
    experiments every time.

``adaptive_design``
    Maintains an exact posterior and picks each next pool to maximise expected
    information gain. This is the contribution being measured.

``adaptive_design(prior=...)``
    The same, warm-started from a literature-derived prior.

``adaptive_design(adaptive=False)``
    **The load-bearing ablation.** Identical Bayesian machinery, but pools are
    chosen at random. Isolating this is what shows the saving comes from
    *informative experiment selection* rather than from Bayesian updating.
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
from .oracle import ColonyOracle, N_CANDIDATES, RECIPE_SIZE


# --------------------------------------------------------------- historical


def historical_protocol(oracle: ColonyOracle) -> tuple[set[int], int]:
    """
    Replay Takahashi & Yamanaka (2006):

    1. every factor alone                      (24 experiments)
    2. the full 24-factor pool                 (1)
    3. leave-one-out from the 24               (24)
    4. the surviving 10-factor pool            (1)
    5. leave-one-out from the 10               (10)
    6. confirm the final four                  (1)

    Total: 61 experiments, regardless of what the data says -- the protocol is
    fixed in advance. That rigidity is exactly what adaptive design removes.
    """
    all_factors = list(range(N_CANDIDATES))

    for factor in all_factors:  # step 1
        oracle.run([factor])

    oracle.run(all_factors)  # step 2

    dropout_24 = {  # step 3
        f: oracle.run([x for x in all_factors if x != f]) for f in all_factors
    }
    shortlist = sorted(dropout_24, key=lambda f: dropout_24[f])[:10]

    oracle.run(shortlist)  # step 4

    dropout_10 = {  # step 5
        f: oracle.run([x for x in shortlist if x != f]) for f in shortlist
    }
    final = sorted(dropout_10, key=lambda f: dropout_10[f])[:RECIPE_SIZE]

    oracle.run(final)  # step 6
    return set(final), oracle.n_experiments


# ----------------------------------------------------------------- adaptive


def _candidate_pools(rng: np.random.Generator, posterior: np.ndarray, n_random: int):
    """
    The original proposal generator, kept for comparison. A mix of two sources:

    * random pools of size 2-12, which split the hypothesis space coarsely and
      are informative while the posterior is still broad;
    * pools built from the currently most-probable factors, which discriminate
      between the few surviving hypotheses late in the run.

    Its defect is the fixed size window. At 24 candidates a pool of 2-12 is a
    sizeable fraction of the space and stands a fair chance of containing most
    of the recipe; at 60 it is a sliver, and a sliver that misses the recipe
    returns the same uninformative zero whatever else is in it. Experiment 8
    measured the consequence: cost growing as about ``n^1.8`` where the
    classical algorithm grows as ``n^0.4``. See :func:`structured_pools`.
    """
    pools = [
        tuple(sorted(rng.choice(N_CANDIDATES, size=int(rng.integers(2, 13)),
                                replace=False).tolist()))
        for _ in range(n_random)
    ]
    ranked = np.argsort(-marginal_factor_probability(posterior))
    for size in (3, 4, 5, 6, 8):
        pools.append(tuple(sorted(int(x) for x in ranked[:size])))
        pools.append(tuple(sorted(int(x) for x in ranked[1 : size + 1])))
    return list(dict.fromkeys(pools))  # de-duplicate, preserve order


def structured_pools(rng: np.random.Generator, posterior: np.ndarray, n_random: int):
    """
    Proposals that include the queries a classical search would make.

    The Bayesian designer's scoring rule was never the problem -- expected
    information gain is exact here, over the full posterior. The problem was
    that it only ever got to score pools that could not be very informative at
    scale. This generator adds the three families the old one lacked, and the
    scoring rule is left untouched, so any change in performance is
    attributable to the proposals alone.

    **Removal pools.** ``half_removal_search`` works by holding a pool that is
    known positive and deleting half of what is left in it. Those are exactly
    the queries that halve the live set in one step, and they were absent from
    the old proposal set at every pool size. Here the live set is taken to be
    the factors carrying appreciable posterior mass, and pools are formed by
    deleting a random half, a random quarter, and a random eighth of it.

    **Ranked prefixes and their complements.** Cuts of the marginal ranking at
    a geometric ladder of depths, plus the complement of each cut. A prefix
    asks "is the recipe inside the top k?"; its complement asks the same
    question from the other side and is the more informative of the two
    whenever the answer is no.

    **Scale-aware random pools.** Sizes drawn log-uniformly between the recipe
    size and the whole candidate set rather than from a fixed 2-12 window, so
    that the coarse splits remain available when ``n`` is large.

    The budget is spent in that order of priority: structure first, random
    filler second. At the default ``n_random`` the structured families are a
    few dozen pools, so most of the budget is still random -- the point is that
    the informative few are now present to be found.
    """
    marginal = marginal_factor_probability(posterior)
    ranked = [int(x) for x in np.argsort(-marginal)]
    n = N_CANDIDATES

    pools: list[tuple[int, ...]] = []

    # -- ranked prefixes and complements, on a geometric ladder
    depths = sorted({3, 4, 5, 6, 8, n // 8, n // 4, n // 2, 3 * n // 4})
    for depth in depths:
        if 1 <= depth < n:
            prefix = ranked[:depth]
            pools.append(tuple(sorted(prefix)))
            pools.append(tuple(sorted(ranked[depth:])))
    pools.append(tuple(sorted(ranked[1:5])))  # drop the leader: a direct test of it

    # -- removal pools against the live set
    live = [f for f in ranked if marginal[f] > 1e-4] or ranked
    for fraction in (2, 4, 8):
        take = max(1, len(live) // fraction)
        for _ in range(3):
            drop = set(rng.choice(live, size=min(take, len(live)), replace=False).tolist())
            kept = tuple(sorted(f for f in live if f not in drop))
            if kept:
                pools.append(kept)

    # -- scale-aware random pools for the rest of the budget
    low, high = np.log(RECIPE_SIZE), np.log(n)
    for _ in range(max(0, n_random - len(pools))):
        size = int(round(np.exp(rng.uniform(low, high))))
        size = min(max(size, 2), n - 1)
        pools.append(
            tuple(sorted(rng.choice(n, size=size, replace=False).tolist()))
        )

    return list(dict.fromkeys(pools))


def adaptive_design(
    oracle: ColonyOracle,
    seed: int = 0,
    max_experiments: int = 60,
    confidence: float = 0.95,
    prior: np.ndarray | None = None,
    n_random_pools: int = 260,
    adaptive: bool = True,
    proposals: str = "original",
) -> tuple[set[int], int, float]:
    """
    Closed-loop experimental design.

    Each round: score candidate pools by expected information gain, run the
    best one, update the posterior. Stop when one hypothesis exceeds
    ``confidence`` or the experiment budget runs out.

    Parameters
    ----------
    prior:
        Optional hypothesis prior (see :func:`inference.weighted_prior`).
    adaptive:
        ``False`` selects pools at random instead -- the ablation.
    proposals:
        Which proposal generator supplies the pools to score: ``"original"``
        (default, so published results stay reproducible) or ``"structured"``.
        See :func:`structured_pools`.

    Returns
    -------
    (recipe, n_experiments, final_confidence)
    """
    rng = np.random.default_rng(seed)
    posterior = uniform_prior() if prior is None else prior.copy()
    propose = structured_pools if proposals == "structured" else _candidate_pools

    while oracle.n_experiments < max_experiments and posterior.max() < confidence:
        pools = propose(rng, posterior, n_random_pools)
        if adaptive:
            pool = max(pools, key=lambda p: expected_information_gain(posterior, p))
        else:
            pool = pools[int(rng.integers(len(pools)))]
        observed = oracle.run(pool)
        posterior = bayes_update(posterior, pool, observed)

    best = HYPOTHESES[int(np.argmax(posterior))]
    return set(best), oracle.n_experiments, float(posterior.max())
