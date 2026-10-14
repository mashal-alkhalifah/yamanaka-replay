#!/usr/bin/env python3
"""
Experiment 12 -- is joint inference over recipe size worth it?

Experiment 10 showed the latent-size designer recovering recipes of three,
four and five factors where the fixed-size designer scores 0%. That comparison
is against a straw man of our own making: nobody facing an unknown number of
defectives would fix it arbitrarily. Group testing's standard answer is to
**estimate the count first, then search** -- and a model-free conjunctive
search sidesteps the question entirely.

Four methods on the same hidden recipes:

``fixed-4``        assumes four factors (the straw man, kept for continuity)
``two-stage-6/12`` spends 6 or 12 experiments estimating the size, then commits
                   to it; the two budgets guard against the baseline failing
                   merely because its first stage was starved
``latent``         carries size and identity in one posterior throughout
``half-removal``   classical, model-free, never represents a size at all

If either ``two-stage`` matches ``latent``, the joint treatment is unnecessary
machinery. If ``half-removal`` matches both at lower cost, the Bayesian
treatment of cardinality is unnecessary altogether. Either outcome is worth
knowing before writing it up as a contribution.

Writes ``results/12_cardinality_baselines.csv``.
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
from yamanaka_replay.baselines import half_removal_search, two_stage_cardinality  # noqa: E402
from yamanaka_replay.cardinality import (  # noqa: E402
    GeneralOracle,
    HypothesisSpace,
    adaptive_design_mixed,
)

N_RUNS = 15
TRUE_SIZES = (3, 4, 5)
ALLOWED_SIZES = (3, 4, 5, 6)
N_CANDIDATES = 24
METHODS = ("fixed-4", "two-stage-6", "two-stage-12", "latent", "half-removal")
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

_SPACE: HypothesisSpace | None = None


def space() -> HypothesisSpace:
    global _SPACE
    if _SPACE is None:
        _SPACE = HypothesisSpace(N_CANDIDATES, ALLOWED_SIZES)
    return _SPACE


def run_one(cell: dict) -> dict:
    k, run, method = cell["true_size"], cell["run"], cell["method"]
    rng = np.random.default_rng(5_000_000 + 100 * run + k)
    truth = set(int(x) for x in rng.choice(N_CANDIDATES, size=k, replace=False))
    lab = GeneralOracle(truth, seed=run)

    t0 = time.perf_counter()
    estimated_size = ""
    if method == "fixed-4":
        found, used, confidence = adaptive_design(lab, seed=run, max_experiments=80)
        estimated_size = 4
    elif method.startswith("two-stage"):
        found, used, confidence, estimated_size = two_stage_cardinality(
            lab,
            space(),
            seed=run,
            stage1_experiments=int(method.split("-")[-1]),
            max_experiments=80,
        )
    elif method == "latent":
        found, used, confidence, size_post = adaptive_design_mixed(
            lab, space(), seed=run, max_experiments=80
        )
        estimated_size = max(size_post, key=size_post.get)
    else:
        found, used = half_removal_search(
            lab, n_candidates=N_CANDIDATES, max_experiments=80, seed=run,
            recipe_size_hint=k,
        )
        confidence = float("nan")
        estimated_size = len(found)
    elapsed = time.perf_counter() - t0

    return {
        "method": method,
        "true_size": k,
        "run": run,
        "n_experiments": used,
        "correct": int(found == truth),
        "n_recovered": len(found & truth),
        "n_spurious": len(found - truth),
        "found_size": len(found),
        "estimated_size": estimated_size,
        "size_correct": int(len(found) == k),
        "final_confidence": round(confidence, 4) if confidence == confidence else "",
        "elapsed_s": round(elapsed, 3),
    }


def main() -> None:
    cells = [
        {"true_size": k, "run": r, "method": m}
        for k in TRUE_SIZES
        for m in METHODS
        for r in range(N_RUNS)
    ]
    print("\n  Experiment 12 -- cardinality baselines")
    print(f"  hypothesis space for the Bayesian methods: {len(space()):,} subsets")
    t0 = time.time()
    rows = []
    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / "12_cardinality_baselines.csv"
    # Written row by row and flushed: a run this long is worth nothing if the
    # process is killed before the final write, which has already happened twice.
    fh = path.open("w", newline="")
    writer = None
    with ProcessPoolExecutor(max_workers=2) as pool:
        step = max(1, len(cells) // 20)
        for i, rec in enumerate(pool.map(run_one, cells, chunksize=1), start=1):
            if writer is None:
                writer = csv.DictWriter(fh, fieldnames=list(rec))
                writer.writeheader()
            writer.writerow(rec)
            fh.flush()
            rows.append(rec)
            if i % step == 0 or i == len(cells):
                print(f"    {i:5d}/{len(cells)}  ({100 * i / len(cells):3.0f}%)", flush=True)
    fh.close()
    print(f"  wrote {path.relative_to(RESULTS.parent)}  ({len(rows)} rows)")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
