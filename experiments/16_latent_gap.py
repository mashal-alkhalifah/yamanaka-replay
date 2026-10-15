#!/usr/bin/env python3
"""
Experiment 16 -- inferring the gap response, against calibrating it first.

Experiment 13 showed the designer's assumed gap response has to be right, and
that this study's assumed value is not: at the measured rho the published
designer recovers 2.5% of recipes while reporting 0.991 confidence. So rho must
come from somewhere. There are two places it can come from, and the choice is
not obvious.

A laboratory would **calibrate first**: run a recipe already known to work and
each of its one-short pools, take the ratio, then fix rho and search. It costs
whole experiments up front and ends in a point estimate carrying its own
Poisson noise.

The alternative is to **infer rho in the loop**, carrying it in the posterior
with the recipe. It spends no experiments up front, but every experiment now
has to inform a nuisance parameter as well as the answer, which could cost more
than the pilot it replaces.

Arms, all on the same drawn rho and the same recipes, with structured
proposals throughout so the proposal generator is not a confound:

``fixed-assumed``     rho = 0.10, the published calibration. The control.
``fixed-oracle``      told the true rho and the true shape. Not available in
                      practice; it bounds what any calibration could buy.
``pilot-then-fixed``  5 experiments on a known working recipe to estimate rho,
                      then the ordinary fixed-rho search at that estimate.
``latent-rho``        joint posterior over recipe and rho on a 9-point grid.
``latent-rho-shape``  the same, with the fall-off shape latent as well (18
                      models), so the designer assumes neither the depth of the
                      gap nor what lies beyond it.

The pilot arm is deliberately advantaged: in this benchmark only the hidden
recipe yields colonies, so the pilot has to be pointed at it, and it therefore
pays a real calibration's cost without its risk of calibrating on the wrong
system. If joint inference matches it anyway, that settles the comparison in
the harder direction.

Writes ``results/16_latent_gap.csv``, flushing each row.
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
from yamanaka_replay import inference  # noqa: E402
from yamanaka_replay.latent_gap import (  # noqa: E402
    GapGrid,
    adaptive_design_latent_gap,
    marginal_model,
    pilot_estimate_rho,
    posterior_rho,
)
from yamanaka_replay.norman_oracle import EmpiricalGapOracle, measured_one_short  # noqa: E402

N_CANDIDATES = 24
N_RUNS = 150
MAX_EXPERIMENTS = 120
N_POOLS = 120
SHAPES = ("gap", "geometric")
METHODS = (
    "fixed-assumed",
    "fixed-oracle",
    "pilot-then-fixed",
    "latent-rho",
    "latent-rho-shape",
)
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

_SAMPLE = measured_one_short()
#: Same draws as experiment 13, so the two campaigns sit on one landscape.
RHO = {
    r: float(np.random.default_rng(31_000 + r).choice(_SAMPLE)) for r in range(N_RUNS)
}


def run_one(cell: dict) -> dict:
    run, method, shape = cell["run"], cell["method"], cell["shape"]
    rho = RHO[run]

    rng = np.random.default_rng(7_000_000 + 1000 * run + N_CANDIDATES)
    truth = set(int(x) for x in rng.choice(N_CANDIDATES, size=4, replace=False))
    lab = EmpiricalGapOracle(truth, seed=run, rho=rho, decay=shape)

    record = {"method": method, "shape": shape, "run": run, "rho": round(rho, 4)}
    t0 = time.perf_counter()

    if method.startswith("latent"):
        shapes = ("gap", "geometric") if method.endswith("shape") else ("gap",)
        found, used, confidence, posterior, grid = adaptive_design_latent_gap(
            lab,
            grid=GapGrid(shapes=shapes),
            seed=run,
            max_experiments=MAX_EXPERIMENTS,
            n_pools=N_POOLS,
        )
        mean, lo, hi = posterior_rho(grid, posterior)
        record["rho_hat"] = round(mean, 4)
        record["rho_in_ci"] = int(lo <= rho <= hi)
        record["pilot_cost"] = 0
        if len(shapes) > 1:
            weights = marginal_model(posterior)
            correct = sum(
                w for w, s in zip(weights, grid.shape) if s == shape
            )
            record["p_true_shape"] = round(float(correct), 4)
    else:
        pilot_cost = 0
        if method == "fixed-oracle":
            assumed, assumed_shape = rho, shape
        elif method == "pilot-then-fixed":
            # A pilot measures the one-short step; it cannot see the shape
            # beyond it, so the search that follows still assumes "gap".
            assumed, pilot_cost = pilot_estimate_rho(lab, truth)
            assumed_shape = "gap"
        else:
            assumed, assumed_shape = 0.10, "gap"

        ColonyOracle.PARTIAL_3 = float(assumed)
        inference.set_assumed_shape(assumed_shape)
        found, used, confidence = adaptive_design(
            lab,
            seed=run,
            max_experiments=MAX_EXPERIMENTS,
            n_random_pools=N_POOLS,
            proposals="structured",
        )
        record["rho_hat"] = round(float(assumed), 4)
        record["rho_in_ci"] = ""
        record["pilot_cost"] = pilot_cost

    record.update(
        {
            "n_experiments": used,  # the pilot's experiments are inside this
            "correct": int(found == truth),
            "n_recovered": len(found & truth),
            "final_confidence": round(float(confidence), 4),
            "silent_failure": int(used < MAX_EXPERIMENTS and found != truth),
            "hit_budget": int(used >= MAX_EXPERIMENTS),
            "rho_abs_error": round(abs(float(record["rho_hat"]) - rho), 4),
            "elapsed_s": round(time.perf_counter() - t0, 3),
        }
    )
    return record


def main() -> None:
    cells = [
        {"run": r, "method": m, "shape": s}
        for s in SHAPES
        for m in METHODS
        for r in range(N_RUNS)
    ]
    print("\n  Experiment 16 -- the gap response as a latent variable")
    print(f"  rho drawn from {_SAMPLE.size} measured ratios; {len(cells)} cells\n")

    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / "16_latent_gap.csv"
    t0 = time.time()
    fields = [
        "method", "shape", "run", "rho", "rho_hat", "rho_in_ci", "p_true_shape",
        "pilot_cost", "n_experiments", "correct", "n_recovered",
        "final_confidence", "silent_failure", "hit_budget", "rho_abs_error",
        "elapsed_s",
    ]
    with path.open("w", newline="") as fh, ProcessPoolExecutor(max_workers=2) as pool:
        writer = csv.DictWriter(fh, fieldnames=fields, restval="")
        writer.writeheader()
        step = max(1, len(cells) // 20)
        for i, rec in enumerate(pool.map(run_one, cells, chunksize=2), start=1):
            writer.writerow(rec)
            fh.flush()
            if i % step == 0 or i == len(cells):
                print(f"    {i:5d}/{len(cells)}  ({100 * i / len(cells):3.0f}%)"
                      f"  [{time.time() - t0:.0f}s]", flush=True)

    print(f"\n  wrote {path.relative_to(RESULTS.parent)}")
    print(f"  elapsed {time.time() - t0:.0f}s\n")


if __name__ == "__main__":
    main()
