#!/usr/bin/env python3
"""
Experiment 11 -- the Bayesian designer against a classical algorithm.

Experiment 8 found the adaptive designer's wet-lab cost growing as about
``n^1.8`` in the candidate count, while the published bounds for stochastic
threshold group testing are ``O(g^2 d log n)``. Either the theory does not
apply or our designer is leaving something on the table. This experiment
settles it by running the classical algorithm the theory describes.

``half_removal_search`` carries no model of the assay beyond a positivity
threshold and makes no assumption about recipe size. If it recovers
logarithmic growth where the Bayesian designer does not, the defect is in our
proposal distribution, not in the problem -- and the honest headline becomes a
*regime* claim: information-gain design wins at the scale Yamanaka actually
worked at, and loses beyond it.

The sweep is dense around the expected crossover. Writes
``results/11_classical_baseline.csv``.
"""

from __future__ import annotations

import csv
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import ColonyOracle, adaptive_design, historical_protocol  # noqa: E402
from yamanaka_replay.baselines import half_removal_search  # noqa: E402
from yamanaka_replay import inference, oracle as oracle_mod, strategies  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _harness import literature_prior, set_candidate_count  # noqa: E402

SIZES = (24, 28, 32, 36, 40, 48, 60)
#: Sizes experiment 8 did not cover; its runs are reused for the rest.
FILL_SIZES = (28, 36)
N_RUNS = 30
N_POOLS = 120
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"


def run_one(cell: dict) -> dict:
    n, run, method = cell["n_candidates"], cell["run"], cell["method"]
    if n != inference.N_CANDIDATES:
        set_candidate_count(n)

    rng = np.random.default_rng(7_000_000 + 1000 * run + n)
    truth = set(int(x) for x in rng.choice(n, size=oracle_mod.RECIPE_SIZE, replace=False))
    lab = ColonyOracle(truth, seed=run)

    t0 = time.perf_counter()
    if method == "historical":
        found, used = historical_protocol(lab)
    elif method == "half-removal":
        found, used = half_removal_search(
            lab, n_candidates=n, max_experiments=200, seed=run
        )
    elif method == "adaptive":
        found, used, _ = adaptive_design(
            lab, seed=run, max_experiments=200, n_random_pools=N_POOLS
        )
    else:  # adaptive with the literature prior held at constant quality
        prior = literature_prior(truth, np.random.default_rng(9_000_000 + run), n)
        found, used, _ = adaptive_design(
            lab, seed=run, max_experiments=200, prior=prior, n_random_pools=N_POOLS
        )
    elapsed = time.perf_counter() - t0

    return {
        "method": method,
        "n_candidates": n,
        "run": run,
        "n_experiments": used,
        "correct": int(found == truth),
        "n_recovered": len(found & truth),
        "found_size": len(found),
        "elapsed_s": round(elapsed, 3),
    }


def main() -> None:
    cells = [
        {"n_candidates": n, "run": r, "method": "half-removal"}
        for n in SIZES
        for r in range(N_RUNS)
    ] + [
        {"n_candidates": n, "run": r, "method": m}
        for n in FILL_SIZES
        for m in ("historical", "adaptive", "informed")
        for r in range(20)
    ]
    print("\n  Experiment 11 -- classical threshold-group-testing baseline")
    print(f"  sizes {SIZES}, {N_RUNS} runs, {len(cells)} cells")
    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=2) as pool:
        step = max(1, len(cells) // 20)
        for i, rec in enumerate(pool.map(run_one, cells, chunksize=4), start=1):
            rows.append(rec)
            if i % step == 0 or i == len(cells):
                print(f"    {i:5d}/{len(cells)}  ({100 * i / len(cells):3.0f}%)", flush=True)

    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / "11_classical_baseline.csv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote {path.relative_to(RESULTS.parent)}  ({len(rows)} rows)")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
