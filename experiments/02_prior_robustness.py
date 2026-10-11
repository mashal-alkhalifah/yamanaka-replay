#!/usr/bin/env python3
"""
Experiment 2 -- how wrong can the literature prior be?

A literature-derived prior is necessarily incomplete: published candidate
lists miss genes. This sweep drops 0, 1, 2 or 3 of the genuinely essential
factors out of the prior and replaces them with decoys, then measures the
cost.

The question is whether a wrong prior *misleads* the search or merely *slows*
it. Because the prior enters as a weighting rather than a constraint, data
should be able to overturn it.

Writes ``results/02_prior_robustness.csv``.
"""

from __future__ import annotations

import csv
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import (  # noqa: E402
    ColonyOracle,
    N_CANDIDATES,
    RECIPE_SIZE,
    adaptive_design,
    weighted_prior,
)

N_RUNS = 14
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"


def degraded_prior(truth, rng, n_dropped: int, strength: float = 3.0,
                   n_informed: int = 8):
    """Prior that omits ``n_dropped`` true factors, substituting decoys."""
    weights = np.ones(N_CANDIDATES)
    true_factors = list(truth)
    rng.shuffle(true_factors)
    retained = true_factors[n_dropped:]
    decoys = [f for f in range(N_CANDIDATES) if f not in truth]
    rng.shuffle(decoys)
    boosted = set(retained) | set(decoys[: n_informed - len(retained)])
    for f in boosted:
        weights[f] = strength
    return weighted_prior(weights)


CASES = [
    (0, "Prior contains all four true factors"),
    (1, "Prior drops one true factor"),
    (2, "Prior drops two true factors"),
    (3, "Prior drops three true factors"),
]


def main() -> None:
    print(f"\n{'=' * 72}")
    print("  Experiment 2 -- robustness to an incomplete literature prior")
    print(f"{'=' * 72}")
    print(f"  {'prior quality':<44}{'median':>9}{'success':>10}")
    print("  " + "-" * 64)

    rows = []
    for n_dropped, label in CASES:
        counts, correct = [], []
        for run in range(N_RUNS):
            rng = np.random.default_rng(7000 + run)
            truth = set(
                int(x) for x in rng.choice(N_CANDIDATES, size=RECIPE_SIZE, replace=False)
            )
            prior = degraded_prior(truth, np.random.default_rng(9000 + run), n_dropped)
            oracle = ColonyOracle(truth, seed=run)
            found, n, _ = adaptive_design(oracle, seed=run, prior=prior)
            counts.append(n)
            correct.append(found == truth)

        median = float(np.median(counts))
        success = 100.0 * float(np.mean(correct))
        print(f"  {label:<44}{median:>9.0f}{success:>9.0f}%")
        rows.append(
            {
                "prior_quality": label,
                "true_factors_dropped": n_dropped,
                "median_experiments": median,
                "success_rate_pct": round(success, 1),
            }
        )

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "02_prior_robustness.csv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  wrote {out.relative_to(RESULTS.parent)}\n")


if __name__ == "__main__":
    main()
