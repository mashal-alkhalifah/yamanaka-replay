#!/usr/bin/env python3
"""
Experiment 9 -- how good does the literature prior have to be?

Experiment 2 degrades the prior along one axis (how many true factors it
misses) at one strength. In a real deployment three things vary independently:

* **coverage** -- how many of the true factors the published candidate list
  actually contains (``prior_dropped``);
* **breadth** -- how many factors it up-weights in total, i.e. how diluted the
  signal is among decoys (``prior_informed``);
* **conviction** -- how hard the prior is pushed (``prior_strength``), from a
  gentle nudge to a near-hard constraint.

The last axis is the interesting one. A strong prior is the cheapest strategy
in the study *when it is right*; the question is what it costs when it is
wrong. A prior strong enough to be a constraint should be expected to turn a
missing factor into an outright failure, and the point of the grid is to find
the strength at which that happens.

Writes ``results/09_prior_quality.csv``.
"""

from __future__ import annotations

import time

from _harness import run_cells, write_rows

N_RUNS = 30
STRENGTHS = (1.5, 3.0, 10.0, 100.0)
DROPPED = (0, 1, 2, 3, 4)
INFORMED = (6, 8, 12)


def main() -> None:
    cells = [
        {
            "strategy": "informed",
            "run": run,
            "prior_strength": s,
            "prior_dropped": d,
            "prior_informed": k,
        }
        for s in STRENGTHS
        for d in DROPPED
        for k in INFORMED
        for run in range(N_RUNS)
    ]
    print("\n  Experiment 9 -- prior quality grid")
    t0 = time.time()
    rows = run_cells(cells, label="prior")
    write_rows(rows, "09_prior_quality.csv")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
