#!/usr/bin/env python3
"""
Experiment 15 -- is the scaling defect in the problem or in the proposals?

Experiment 8 measured the adaptive designer's cost growing as about ``n^1.8``
in the candidate count. Experiment 11 ran the classical algorithm on the same
runs and found ``n^0.4``, which forced the paper's headline down from "adaptive
design wins" to "adaptive design wins below about 48 candidates". Experiment 11
also named the suspect: if a model-free algorithm with no posterior at all
scales better, the Bayesian machinery is not what is failing -- the set of
pools it is allowed to score is.

This experiment tests that directly. The scoring rule, the posterior, the
oracle, the recipes and the seeds are all held fixed; the only thing that
changes is which pools the designer gets to choose among:

``adaptive-original``    the published generator: random pools of size 2-12
                         plus a few prefixes of the marginal ranking
``adaptive-structured``  the same budget of pools, but including removal pools
                         against the live set, complements of ranked prefixes,
                         and random sizes drawn across the whole range rather
                         than a fixed window

Both are given the same number of pools to score, so a difference cannot be
bought with a larger search. What the structured generator costs instead is
time per pool, since its pools are larger -- ``elapsed_s`` is recorded so that
price is visible rather than hidden.

``half-removal`` is not re-run: experiment 11 already ran it over these sizes
with the same seeds and recipes, and its rows are the comparison curve.

Writes ``results/15_proposal_ablation.csv``, flushing each row.
"""

from __future__ import annotations

import csv
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import ColonyOracle, adaptive_design  # noqa: E402
from yamanaka_replay import inference, oracle as oracle_mod  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _harness import set_candidate_count  # noqa: E402

#: Runs thin out as the hypothesis space grows: at 60 candidates there are
#: 487,635 hypotheses and a single run costs minutes.
RUNS_AT = {24: 20, 32: 20, 40: 20, 48: 20, 60: 12}
N_POOLS = 120
MAX_EXPERIMENTS = 150
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"


def run_one(cell: dict) -> dict:
    n, run, proposals = cell["n_candidates"], cell["run"], cell["proposals"]
    if n != inference.N_CANDIDATES:
        set_candidate_count(n)

    # Identical to experiments 8 and 11, so the recipes and noise are shared.
    rng = np.random.default_rng(7_000_000 + 1000 * run + n)
    truth = set(int(x) for x in rng.choice(n, size=oracle_mod.RECIPE_SIZE, replace=False))
    lab = ColonyOracle(truth, seed=run)

    t0 = time.perf_counter()
    found, used, confidence = adaptive_design(
        lab,
        seed=run,
        max_experiments=MAX_EXPERIMENTS,
        n_random_pools=N_POOLS,
        proposals=proposals,
    )
    elapsed = time.perf_counter() - t0

    pools = [p for p, _ in lab.log]
    return {
        "method": f"adaptive-{proposals}",
        "proposals": proposals,
        "n_candidates": n,
        "run": run,
        "n_experiments": used,
        "correct": int(found == truth),
        "n_recovered": len(found & truth),
        "final_confidence": round(confidence, 4),
        "hit_budget": int(used >= MAX_EXPERIMENTS),
        "mean_pool_size": round(float(np.mean([len(p) for p in pools])), 2),
        "max_pool_size": max((len(p) for p in pools), default=0),
        "elapsed_s": round(elapsed, 3),
    }


def main() -> None:
    # Small sizes first, so an interrupted campaign still leaves a usable curve.
    cells = [
        {"n_candidates": n, "run": r, "proposals": p}
        for n in sorted(RUNS_AT)
        for r in range(RUNS_AT[n])
        for p in ("original", "structured")
    ]
    print("\n  Experiment 15 -- proposal generator ablation")
    print(f"  sizes {sorted(RUNS_AT)}, {len(cells)} cells\n")

    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / "15_proposal_ablation.csv"
    t0 = time.time()
    with path.open("w", newline="") as fh, ProcessPoolExecutor(max_workers=2) as pool:
        writer = None
        for i, rec in enumerate(pool.map(run_one, cells, chunksize=1), start=1):
            if writer is None:
                writer = csv.DictWriter(fh, fieldnames=list(rec))
                writer.writeheader()
            writer.writerow(rec)
            fh.flush()
            if i % 10 == 0 or i == len(cells):
                print(f"    {i:4d}/{len(cells)}  ({100 * i / len(cells):3.0f}%)"
                      f"  [{time.time() - t0:.0f}s]", flush=True)

    print(f"\n  wrote {path.relative_to(RESULTS.parent)}")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
