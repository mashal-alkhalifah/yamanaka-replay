#!/usr/bin/env python3
"""
Experiment 7 -- the cost of certainty, and whether the posterior is honest.

The designer stops when the leading hypothesis passes a confidence threshold.
That threshold is a free parameter nobody has justified, and it buys two
different things: fewer experiments when it is loose, fewer mistakes when it is
tight. This sweep prices that trade.

The second question is more important than the first. A Bayesian stopping rule
is only usable if the posterior is **calibrated** -- if "95% confident" means
wrong about one time in twenty. If the reported confidence systematically
exceeds the empirical success rate, the rule is dangerous no matter how cheap
it is, because in a real campaign there is nothing to check it against.

Writes ``results/07_stopping_rule.csv``.
"""

from __future__ import annotations

import time

from _harness import run_cells, write_rows

N_RUNS = 120
CONFIDENCE_GRID = (0.50, 0.75, 0.90, 0.95, 0.99, 0.999, 0.9999)


def main() -> None:
    cells = [
        {
            "strategy": strategy,
            "run": run,
            "confidence": c,
            "confidence_target": c,
        }
        for c in CONFIDENCE_GRID
        for strategy in ("adaptive", "informed")
        for run in range(N_RUNS)
    ]
    print("\n  Experiment 7 -- stopping threshold and posterior calibration")
    t0 = time.time()
    rows = run_cells(cells, label="stopping")
    write_rows(rows, "07_stopping_rule.csv")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
