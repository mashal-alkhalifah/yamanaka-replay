#!/usr/bin/env python3
"""
Experiment 14 -- discovery against measured data, with no simulator anywhere.

    python experiments/14_norman_retrospective.py <pseudobulk.csv> <perturbations.csv>

Every other experiment in this repository queries a simulator. This one queries
Norman et al. (2019): each "experiment" returns the measured transcriptional
response of a condition that was actually profiled, perturbed only by the
dataset's own internal disagreement. Nothing is interpolated, and a condition
the authors did not profile simply cannot be asked about.

That restriction changes the task. With only singles and pairs in the library,
there are no large pools to design, so this is not threshold group testing; it
is sequential discovery within a measured library -- find the strongest
combination while profiling as few conditions as possible. It tests the half of
the claim that survives the restriction: that choosing the next experiment from
a model beats choosing it by protocol.

Arms
----
``exhaustive``        profile everything. Always correct, and the cost the
                      other arms are trying to avoid.
``random``            profile in random order. The floor.
``singles-first``     the protocol a bench scientist would actually run, and
                      the structural cousin of the 2006 leave-one-out design:
                      profile every single, rank them, then profile the pairs
                      built from the best singles.
``surrogate``         fit a cheap model of pair response from the singles seen
                      so far and profile the condition with the highest
                      upper-confidence score. This is the "virtual cell as
                      surrogate" arrangement in miniature -- a model standing
                      in for the expensive assay and spending its budget where
                      the model is both optimistic and unsure.

What to watch
-------------
The audit found measured combinations to be far more additive than the 2006
calibration assumes: a single constituent reaches a median 0.77 of its pair's
effect. If that holds, ``singles-first`` should be strong here -- the cheap
protocol is nearly sufficient *because* the landscape is shallow. A surrogate
that cannot beat it by much is reporting a property of the biology, not a
failure of the method, and the paper should say so.

Writes ``results/14_norman_retrospective.csv`` (one row per arm per replicate)
and ``results/14_norman_curves.csv`` (best-found-so-far against budget).
"""

from __future__ import annotations

import csv
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay.norman_oracle import (  # noqa: E402
    SPLIT_HALF_CV,
    NormanOracle,
    estimate_split_half_cv,
)

N_RUNS = 40
BUDGET = 120  # about half the 236-condition library
KAPPA = 1.5  # optimism in the surrogate's acquisition score
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"


# ------------------------------------------------------------------ surrogate


def _fit(singles: dict, pairs: dict) -> tuple[np.ndarray, float, float, float]:
    """
    Ridge fit of ``log v_pair`` on its two constituents' ``log v_single``.

    Returns ``(coefficients, residual sd, mean log single, sd log single)``.
    The constituent logs are entered sorted (larger first), because the pair is
    unordered and letting the fit see an arbitrary order would invent a
    distinction the data does not have.
    """
    usable = [
        (k, v) for k, v in pairs.items() if all(g in singles for g in k)
    ]
    logs = np.log([v for v in singles.values()]) if singles else np.array([0.0])
    mu_s, sd_s = float(logs.mean()), float(logs.std() or 0.3)

    if len(usable) < 5:
        return np.array([mu_s, 0.0, 0.0]), max(sd_s, 0.3), mu_s, sd_s

    rows, y = [], []
    for key, value in usable:
        a, b = sorted((np.log(singles[g]) for g in key), reverse=True)
        rows.append([1.0, a, b])
        y.append(np.log(value))
    X, y = np.asarray(rows), np.asarray(y)
    ridge = 1e-3 * np.eye(3)
    coef = np.linalg.solve(X.T @ X + ridge, X.T @ y)
    resid = float(np.std(y - X @ coef))
    return coef, max(resid, 1e-3), mu_s, sd_s


def _predict(key, singles, coef, resid, mu_s, sd_s) -> tuple[float, float]:
    """Predicted ``log`` value and its uncertainty for an unmeasured condition."""
    genes = sorted(key)
    if len(genes) == 1:
        # Nothing in the data predicts one single from another: the best the
        # surrogate can say is the population it was drawn from.
        return mu_s, sd_s
    known = [np.log(singles[g]) for g in genes if g in singles]
    if len(known) < 2:
        # One or both constituents unprofiled. Substituting the mean keeps the
        # prediction defined, but the inflated sd records that it is a guess --
        # which is what makes the acquisition explore the missing singles.
        filled = sorted(known + [mu_s] * (2 - len(known)), reverse=True)
        penalty = sd_s * (2 - len(known))
        return float(coef @ [1.0, *filled]), float(np.hypot(resid, penalty))
    a, b = sorted(known, reverse=True)
    return float(coef @ [1.0, a, b]), resid


# --------------------------------------------------------------------- arms


def _record(trace: list[float], target: float, arm: str, run: int) -> dict:
    """Summarise one replicate's trace of best-value-so-far."""
    hit = next((i + 1 for i, v in enumerate(trace) if v >= target - 1e-9), None)
    return {
        "arm": arm,
        "run": run,
        "n_to_best": hit if hit is not None else "",
        "found_best": int(hit is not None),
        "best_at_20": round(trace[19] / target, 4) if len(trace) >= 20 else "",
        "best_at_60": round(trace[59] / target, 4) if len(trace) >= 60 else "",
        "n_experiments": len(trace),
    }


def _play(order, oracle, target, arm, run):
    """Profile conditions in a fixed order, tracking the best true value seen."""
    best, trace = 0.0, []
    for key in order:
        oracle.run(key)
        best = max(best, oracle.truth(key))
        trace.append(best)
    return _record(trace, target, arm, run), trace


def run_surrogate(oracle: NormanOracle, target: float, run: int, budget: int):
    rng = np.random.default_rng(5_000 + run)
    pool = list(oracle.available())
    singles_obs: dict[str, float] = {}
    pairs_obs: dict[frozenset, float] = {}
    best, trace = 0.0, []

    # Seed: three conditions at random. Without any observation the surrogate
    # has nothing to rank by, and seeding at random rather than by a rule keeps
    # the arm from smuggling in prior knowledge of the answer.
    seed_idx = rng.choice(len(pool), size=3, replace=False)
    queue = [pool[i] for i in seed_idx]

    while len(trace) < budget:
        if queue:
            key = queue.pop()
        else:
            coef, resid, mu_s, sd_s = _fit(singles_obs, pairs_obs)
            scored = [
                (_predict(k, singles_obs, coef, resid, mu_s, sd_s), k)
                for k in pool
            ]
            key = max(scored, key=lambda t: t[0][0] + KAPPA * t[0][1])[1]
        pool.remove(key)
        value = oracle.run(key)
        if len(key) == 1:
            singles_obs[next(iter(key))] = value
        else:
            pairs_obs[key] = value
        best = max(best, oracle.truth(key))
        trace.append(best)
        if not pool:
            break
    return _record(trace, target, "surrogate", run), trace


def run_singles_first(oracle: NormanOracle, target: float, run: int, budget: int):
    rng = np.random.default_rng(6_000 + run)
    conditions = oracle.available()
    singles = [k for k in conditions if len(k) == 1]
    pairs = [k for k in conditions if len(k) == 2]
    rng.shuffle(singles)

    observed = {}
    best, trace = 0.0, []
    for key in singles:
        if len(trace) >= budget:
            break
        observed[next(iter(key))] = oracle.run(key)
        best = max(best, oracle.truth(key))
        trace.append(best)

    # Rank pairs by the sum of their constituents' measured singles -- the
    # additive guess, which is what ranking by "best single factors" means.
    def score(key):
        return sum(observed.get(g, 0.0) for g in key)

    for key in sorted(pairs, key=score, reverse=True):
        if len(trace) >= budget:
            break
        oracle.run(key)
        best = max(best, oracle.truth(key))
        trace.append(best)
    return _record(trace, target, "singles-first", run), trace


def main() -> None:
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    pseudobulk, perturbations = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])

    probe = NormanOracle(pseudobulk, perturbations, noiseless=True)
    target_key, target = probe.best()
    print("\n  Experiment 14 -- retrospective discovery on measured data")
    print(f"  {len(probe.available())} measured conditions, {len(probe.candidates)} genes")
    print(f"  target: {'+'.join(sorted(target_key))}  ||delta|| = {target:.2f}")
    measured_cv = estimate_split_half_cv(pseudobulk)
    print(f"  split-half noise: measured {measured_cv:.3f}, "
          f"oracle uses {SPLIT_HALF_CV:.3f}\n")

    rows, curves = [], []
    for run in range(N_RUNS):
        conditions = probe.available()

        oracle = NormanOracle(pseudobulk, perturbations, seed=run)
        order = list(conditions)
        np.random.default_rng(4_000 + run).shuffle(order)
        rec, trace = _play(order[:BUDGET], oracle, target, "random", run)
        rows.append(rec), curves.append((run, "random", trace))

        rec, trace = run_singles_first(
            NormanOracle(pseudobulk, perturbations, seed=run), target, run, BUDGET
        )
        rows.append(rec), curves.append((run, "singles-first", trace))

        rec, trace = run_surrogate(
            NormanOracle(pseudobulk, perturbations, seed=run), target, run, BUDGET
        )
        rows.append(rec), curves.append((run, "surrogate", trace))

        if (run + 1) % 10 == 0:
            print(f"    {run + 1:3d}/{N_RUNS}", flush=True)

    # Exhaustive is deterministic: one row, stating the cost of certainty.
    rows.append(
        {
            "arm": "exhaustive",
            "run": 0,
            "n_to_best": len(probe.available()),
            "found_best": 1,
            "best_at_20": "",
            "best_at_60": "",
            "n_experiments": len(probe.available()),
        }
    )

    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / "14_norman_retrospective.csv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    curve_path = RESULTS / "14_norman_curves.csv"
    with curve_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["arm", "run", "n_experiments", "best_fraction_of_target"])
        for run, arm, trace in curves:
            for i, value in enumerate(trace, start=1):
                writer.writerow([arm, run, i, round(value / target, 5)])

    for arm in ("random", "singles-first", "surrogate"):
        hits = [r["n_to_best"] for r in rows if r["arm"] == arm and r["found_best"]]
        found = sum(1 for r in rows if r["arm"] == arm and r["found_best"])
        at20 = [r["best_at_20"] for r in rows if r["arm"] == arm and r["best_at_20"] != ""]
        print(f"  {arm:15s} found the best in {found}/{N_RUNS} runs"
              + (f", median {int(np.median(hits))} experiments" if hits else "")
              + (f"; at 20 experiments reached {100 * np.median(at20):.0f}% of it" if at20 else ""))
    print(f"\n  wrote {path.name} and {curve_path.name}\n")


if __name__ == "__main__":
    main()
