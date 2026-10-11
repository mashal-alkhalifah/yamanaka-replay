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
    Propose pools to score. A mix of two sources:

    * random pools of size 2-12, which split the hypothesis space coarsely and
      are informative while the posterior is still broad;
    * pools built from the currently most-probable factors, which discriminate
      between the few surviving hypotheses late in the run.
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


def adaptive_design(
    oracle: ColonyOracle,
    seed: int = 0,
    max_experiments: int = 60,
    confidence: float = 0.95,
    prior: np.ndarray | None = None,
    n_random_pools: int = 260,
    adaptive: bool = True,
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

    Returns
    -------
    (recipe, n_experiments, final_confidence)
    """
    rng = np.random.default_rng(seed)
    posterior = uniform_prior() if prior is None else prior.copy()

    while oracle.n_experiments < max_experiments and posterior.max() < confidence:
        pools = _candidate_pools(rng, posterior, n_random_pools)
        if adaptive:
            pool = max(pools, key=lambda p: expected_information_gain(posterior, p))
        else:
            pool = pools[int(rng.integers(len(pools)))]
        observed = oracle.run(pool)
        posterior = bayes_update(posterior, pool, observed)

    best = HYPOTHESES[int(np.argmax(posterior))]
    return set(best), oracle.n_experiments, float(posterior.max())
