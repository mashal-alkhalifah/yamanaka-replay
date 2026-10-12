#!/usr/bin/env python3
"""
Experiment 10 -- removing the method's one structural failure.

Experiments 3 and 6 together isolate a single structural weakness: if the true
recipe is not the size the searcher assumes, the search fails, and it fails
*silently* -- the posterior concentrates confidently on a wrong answer. Every
other perturbation tested is calibration error, which degrades gracefully.

A silent failure is not a limitation to document; it is a bug to fix. This
experiment compares two designers on hidden recipes of three, four and five
factors:

``fixed-4``
    The designer from the main study, which assumes four.
``latent``
    The designer in ``cardinality.py``, which carries a posterior over both the
    recipe and its size.

The quantities that matter are (a) does the latent designer recover recipes of
the wrong size, (b) what premium does it pay on the four-factor case it did not
need to generalise for, and (c) when it is wrong, does it *say* so -- which is
the property the fixed designer lacks entirely.

Writes ``results/10_unknown_cardinality.csv``.
"""

from __future__ import annotations

import csv
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import adaptive_design  # noqa: E402
from yamanaka_replay.cardinality import (  # noqa: E402
    GeneralOracle,
    HypothesisSpace,
    adaptive_design_mixed,
)

N_RUNS = 40
TRUE_SIZES = (3, 4, 5)
ALLOWED_SIZES = (3, 4, 5, 6)
N_CANDIDATES = 24
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

_SPACE: HypothesisSpace | None = None


def space() -> HypothesisSpace:
    """Built once per worker: 180,180 hypotheses is a second of setup."""
    global _SPACE
    if _SPACE is None:
        _SPACE = HypothesisSpace(N_CANDIDATES, ALLOWED_SIZES)
    return _SPACE


def run_one(cell: dict) -> dict:
    true_size, run, designer = cell["true_size"], cell["run"], cell["designer"]
    rng = np.random.default_rng(5_000_000 + 100 * run + true_size)
    truth = set(int(x) for x in rng.choice(N_CANDIDATES, size=true_size, replace=False))

    lab = GeneralOracle(truth, seed=run)
    t0 = time.perf_counter()
    if designer == "fixed-4":
        found, n, confidence = adaptive_design(lab, seed=run, max_experiments=80)
        size_posterior = {4: 1.0}
    else:
        found, n, confidence, size_posterior = adaptive_design_mixed(
            lab, space(), seed=run, max_experiments=80
        )
    elapsed = time.perf_counter() - t0

    return {
        "designer": designer,
        "true_size": true_size,
        "run": run,
        "n_experiments": n,
        "correct": int(found == truth),
        "found_size": len(found),
        "size_correct": int(len(found) == true_size),
        "n_recovered": len(found & truth),
        "n_spurious": len(found - truth),
        "final_confidence": round(confidence, 4),
        "p_true_size": round(float(size_posterior.get(true_size, 0.0)), 4),
        "hit_budget": int(n >= 80),
        "elapsed_s": round(elapsed, 3),
    }


def main() -> None:
    cells = [
        {"true_size": k, "run": r, "designer": d}
        for k in TRUE_SIZES
        for d in ("fixed-4", "latent")
        for r in range(N_RUNS)
    ]
    print("\n  Experiment 10 -- unknown recipe cardinality")
    print(f"  hypothesis space: {len(space()):,} subsets of sizes {ALLOWED_SIZES}")
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=2) as pool:
        step = max(1, len(cells) // 20)
        for i, rec in enumerate(pool.map(run_one, cells, chunksize=2), start=1):
            rows.append(rec)
            if i % step == 0 or i == len(cells):
                print(f"    {i:5d}/{len(cells)}  ({100 * i / len(cells):3.0f}%)", flush=True)

    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / "10_unknown_cardinality.csv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote {path.relative_to(RESULTS.parent)}  ({len(rows)} rows)")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
