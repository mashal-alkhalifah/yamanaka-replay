#!/usr/bin/env python3
"""
Experiment 6 -- how much graded synergy does the method actually need?

Experiment 3 showed that removing the three-of-four partial signal entirely
breaks the search. That leaves the quantitative question unanswered: the
calibrated value is 0.10 (three of four gives ~10% of the full rate, ~4
colonies), and somewhere between 0.10 and 0 the method must fail. Where?

The sweep separates two effects that experiment 3 conflates:

``matched``
    The searcher's likelihood model is told the true partial efficiency. Any
    failure here is a **real loss of gradient** -- the information is not in
    the assay.
``misspecified``
    The searcher keeps believing 0.10 while the assay does something else.
    Failure here is a **modelling error**, repairable by calibration.

The gap between the two arms is the part of the failure that better modelling
could fix. Writes ``results/06_synergy_threshold.csv``.
"""

from __future__ import annotations

import time

from _harness import run_cells, write_rows

N_RUNS = 50
PARTIAL_GRID = (0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.25, 0.40)


def main() -> None:
    cells = []
    for p in PARTIAL_GRID:
        for run in range(N_RUNS):
            # Arm 1: the searcher is told the truth about the partial effect.
            cells.append(
                {
                    "strategy": "adaptive",
                    "run": run,
                    "partial_3": p,
                    "assumed_partial_3": max(p, 1e-6),
                    "arm": "matched",
                    "oracle_partial_3": p,
                }
            )
            # Arm 2: the searcher keeps the calibrated assumption of 0.10.
            cells.append(
                {
                    "strategy": "adaptive",
                    "run": run,
                    "partial_3": p,
                    "arm": "misspecified",
                    "oracle_partial_3": p,
                }
            )
    print("\n  Experiment 6 -- graded-synergy threshold")
    t0 = time.time()
    rows = run_cells(cells, label="synergy")
    write_rows(rows, "06_synergy_threshold.csv")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
