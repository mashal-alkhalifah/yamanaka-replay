#!/usr/bin/env python3
"""
Experiment 8 -- does the saving survive a bigger candidate list?

Yamanaka's 24 candidates were themselves the product of a prior literature
filter. A method that only works once that filter has been applied is of
limited use, so this asks what happens as the candidate list grows: 24, 32, 40,
48 and 60 factors, i.e. hypothesis spaces from 10,626 up to 487,635.

Two quantities are tracked separately:

* **experiments** -- the wet-lab cost, which is what the study is about;
* **wall-clock seconds** -- the compute cost of *choosing* each experiment,
  which the five-overlap-class collapse in ``inference.py`` predicts should
  grow far more slowly than the hypothesis count.

The historical protocol is included at each size as the honest baseline: it
costs ``2n + 13`` experiments by construction, so it grows linearly while the
search is expected to grow logarithmically.

Writes ``results/08_scaling.csv``.
"""

from __future__ import annotations

import time

from _harness import run_cells, write_rows

N_RUNS = 20
SIZES = (24, 32, 40, 48, 60)

#: Held fixed across sizes so that cost comparisons are not confounded by the
#: number of proposals scored per round.
N_POOLS = 120


def main() -> None:
    cells = [
        {
            "strategy": strategy,
            "run": run,
            "n_candidates": n,
            "max_experiments": 150,
            "n_random_pools": N_POOLS,
        }
        for n in SIZES
        for strategy in ("historical", "adaptive", "informed")
        for run in range(N_RUNS)
    ]
    print("\n  Experiment 8 -- scaling in candidate-list size")
    t0 = time.time()
    rows = run_cells(cells, label="scaling")
    write_rows(rows, "08_scaling.csv")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
