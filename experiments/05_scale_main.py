#!/usr/bin/env python3
"""
Experiment 5 -- the main comparison at scale.

Experiment 1 runs 20 replicates, which is enough to show the effect but not
enough to put an interval on it or to characterise the tail. This runs 300
independent replicates of all four strategies on the same 300 hidden recipes,
so every comparison is **paired**: the same recipe and the same noise seed are
given to each strategy.

Writes one row per (strategy, run) to ``results/05_scale_main.csv``.
Aggregation and statistics are done separately in ``analysis/analyze.py``.
"""

from __future__ import annotations

import time

from _harness import run_cells, write_rows

N_RUNS = 300
STRATEGIES = ("historical", "adaptive", "informed", "random")


def main() -> None:
    cells = [
        {"strategy": s, "run": r} for r in range(N_RUNS) for s in STRATEGIES
    ]
    print("\n  Experiment 5 -- main comparison, 300 paired replicates")
    t0 = time.time()
    rows = run_cells(cells, label="main")
    write_rows(rows, "05_scale_main.csv")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
