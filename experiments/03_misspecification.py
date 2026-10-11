#!/usr/bin/env python3
"""
Experiment 3 -- where the method breaks.

The adaptive designer assumes a likelihood model. Here the oracle is made to
disagree with it, one assumption at a time, to find the boundary of validity.

Two classes of error are distinguished, and they behave very differently:

*Calibration* errors -- the rates are off, but the shape of the landscape is
intact. The designer absorbs these.

*Structural* errors -- the landscape itself is different. Two of these are
fatal:

1. **No graded synergy.** If three-factor pools produce nothing at all, there
   is no gradient to climb and the search is reduced to blind enumeration.
2. **Wrong recipe cardinality.** Searching a fixed size of four can never find
   a five-factor truth. Success is exactly zero, and -- the dangerous part --
   the posterior still concentrates confidently on a wrong answer. The failure
   is silent.

The second finding is the reason recipe size must be an inferred latent
variable in any real deployment, not a fixed hyperparameter.

Writes ``results/03_misspecification.csv``.
"""

from __future__ import annotations

import csv
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import (  # noqa: E402
    MisspecifiedOracle,
    N_CANDIDATES,
    RECIPE_SIZE,
    adaptive_design,
)

N_RUNS = 14
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

CASES = [
    ("Oracle matches the designer's model", "calibration", {}),
    ("Pool-size penalty 3x harsher", "calibration", {"decay": 0.107}),
    ("Triplet effect stronger than assumed", "calibration", {"partial_3": 0.35}),
    ("Overall efficiency 3x lower", "calibration", {"base": 15.0}),
    ("Triplets have no effect at all", "structural", {"partial_3": 0.0}),
    ("Truth is five factors, not four", "structural", {"five_factors": True}),
]


def main() -> None:
    print(f"\n{'=' * 78}")
    print("  Experiment 3 -- robustness to model misspecification")
    print(f"{'=' * 78}")
    print(f"  {'oracle reality vs designer assumption':<46}{'kind':>13}"
          f"{'median':>8}{'success':>9}")
    print("  " + "-" * 74)

    rows = []
    for label, kind, overrides in CASES:
        counts, correct = [], []
        for run in range(N_RUNS):
            rng = np.random.default_rng(7000 + run)
            truth = set(
                int(x) for x in rng.choice(N_CANDIDATES, size=RECIPE_SIZE, replace=False)
            )
            kwargs = dict(overrides)
            if kwargs.pop("five_factors", False):
                # a fifth factor is also required, unknown to the designer
                kwargs["extra_required"] = next(
                    f for f in range(N_CANDIDATES) if f not in truth
                )
            oracle = MisspecifiedOracle(truth, seed=run, **kwargs)
            found, n, _ = adaptive_design(oracle, seed=run)
            counts.append(n)
            correct.append(found == truth)

        median = float(np.median(counts))
        success = 100.0 * float(np.mean(correct))
        print(f"  {label:<46}{kind:>13}{median:>8.0f}{success:>8.0f}%")
        rows.append(
            {
                "case": label,
                "error_kind": kind,
                "median_experiments": median,
                "success_rate_pct": round(success, 1),
            }
        )

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "03_misspecification.csv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  wrote {out.relative_to(RESULTS.parent)}\n")


if __name__ == "__main__":
    main()
