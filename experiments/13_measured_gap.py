#!/usr/bin/env python3
"""
Experiment 13 -- the designer under the gap structure measured data shows.

Every number in experiments 01-12 is conditional on one assumed constant: that
a pool one factor short of the recipe performs at ``rho = 0.10`` of a complete
one. That figure was read off Takahashi & Yamanaka (2006), the same paper the
study replays, so replaying it cannot test it. ``analysis/real_data_audit.py``
tested it instead against measured combinatorial perturbation data and found
the analogous quantity has median 0.77 in Norman et al. (2019), never below
0.20 in 131 combinations, and at or above 1.0 in 29% of them -- cases where one
constituent moves the cell as far as the pair does.

This experiment stops treating that as a caveat and runs the method against it.

For each replicate a value of ``rho`` is drawn from the 131 measured ratios and
held fixed across every arm, so the arms are paired on the same landscape:

``calibrated``        the paper's own rho = 0.10, as a reference line
``adaptive-assumed``  information-gain design whose likelihood still says 0.10
                      -- what a reader applying this repository unmodified gets
``adaptive-matched``  the same designer told the true rho. The gap between this
                      and ``adaptive-assumed`` is the cost of misspecification;
                      the gap between this and ``calibrated`` is the cost of
                      the landscape itself, which no amount of correct
                      modelling removes
``half-removal``      the classical algorithm, threshold from assumed constants
``half-removal-matched``  the same with the threshold rebuilt at the true rho
``historical``        the 2006 protocol, which models nothing and so cannot be
                      misspecified

Two fall-off shapes are run. Pair data constrains only the one-short step, so
whether efficiency then drops to zero (``gap``) or decays smoothly
(``geometric``) is undetermined by the measurements; reporting one alone would
be choosing the convenient half of an open question.

Rows are written as they complete -- a long campaign in this container has been
lost twice to reclamation. Writes ``results/13_measured_gap.csv``.
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
from yamanaka_replay import inference  # noqa: E402
from yamanaka_replay.baselines import half_removal_search  # noqa: E402
from yamanaka_replay.norman_oracle import EmpiricalGapOracle, measured_one_short  # noqa: E402

N_CANDIDATES = 24
N_RUNS = 200
MAX_EXPERIMENTS = 120
N_POOLS = 200
SHAPES = ("gap", "geometric")
METHODS = (
    "calibrated",
    "adaptive-assumed",
    "adaptive-matched",
    "half-removal",
    "half-removal-matched",
    "historical",
)
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

#: One draw of rho per replicate, shared by every arm so comparisons are paired.
_SAMPLE = measured_one_short()
RHO = {
    r: float(np.random.default_rng(31_000 + r).choice(_SAMPLE)) for r in range(N_RUNS)
}


def run_one(cell: dict) -> dict:
    run, method, shape = cell["run"], cell["method"], cell["shape"]
    rho = 0.10 if method == "calibrated" else RHO[run]

    rng = np.random.default_rng(7_000_000 + 1000 * run + N_CANDIDATES)
    truth = set(int(x) for x in rng.choice(N_CANDIDATES, size=4, replace=False))
    lab = EmpiricalGapOracle(truth, seed=run, rho=rho, decay=shape)

    # What the searcher believes. inference.expected_rate and the half-removal
    # positivity threshold both read these, so they are set together.
    matched = method.endswith("matched") or method == "calibrated"
    ColonyOracle.PARTIAL_3 = rho if matched else 0.10
    inference.set_assumed_shape(shape if matched else "gap")

    t0 = time.perf_counter()
    if method == "historical":
        found, used = historical_protocol(lab)
        confidence = float("nan")
    elif method.startswith("half-removal"):
        found, used = half_removal_search(
            lab, n_candidates=N_CANDIDATES, max_experiments=MAX_EXPERIMENTS, seed=run
        )
        confidence = float("nan")
    else:
        found, used, confidence = adaptive_design(
            lab,
            seed=run,
            max_experiments=MAX_EXPERIMENTS,
            n_random_pools=N_POOLS,
        )
    elapsed = time.perf_counter() - t0

    return {
        "method": method,
        "shape": shape,
        "run": run,
        "rho": round(rho, 4),
        "n_experiments": used,
        "correct": int(found == truth),
        "n_recovered": len(found & truth),
        "found_size": len(found),
        "final_confidence": round(confidence, 4) if confidence == confidence else "",
        # A run that stops confident and wrong is a silent failure: the
        # practitioner has no signal that anything went astray. Counting these
        # separately from budget exhaustion is the point of the experiment.
        "silent_failure": int(
            method.startswith("adaptive")
            and used < MAX_EXPERIMENTS
            and found != truth
        ),
        "hit_budget": int(used >= MAX_EXPERIMENTS),
        "elapsed_s": round(elapsed, 3),
    }


def main() -> None:
    cells = [
        {"run": r, "method": m, "shape": s}
        for s in SHAPES
        for m in METHODS
        # The calibrated reference does not depend on the measured rho, so it is
        # run once under the gap shape only.
        if not (m == "calibrated" and s == "geometric")
        for r in range(N_RUNS)
    ]
    print("\n  Experiment 13 -- measured gap structure")
    print(f"  rho drawn from {_SAMPLE.size} Norman ratios "
          f"(median {np.median(_SAMPLE):.2f}, {100 * np.mean(_SAMPLE >= 1):.0f}% at or above 1)")
    print(f"  {len(cells)} cells\n")

    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / "13_measured_gap.csv"
    t0 = time.time()
    written = 0
    with path.open("w", newline="") as fh, ProcessPoolExecutor(max_workers=2) as pool:
        writer = None
        step = max(1, len(cells) // 20)
        for i, rec in enumerate(pool.map(run_one, cells, chunksize=4), start=1):
            if writer is None:
                writer = csv.DictWriter(fh, fieldnames=list(rec))
                writer.writeheader()
            writer.writerow(rec)
            fh.flush()  # survive container reclamation mid-campaign
            written += 1
            if i % step == 0 or i == len(cells):
                print(f"    {i:5d}/{len(cells)}  ({100 * i / len(cells):3.0f}%)", flush=True)

    print(f"\n  wrote {path.relative_to(RESULTS.parent)}  ({written} rows)")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
