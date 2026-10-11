#!/usr/bin/env python3
"""
Experiment 1 -- main comparison.

How many wet-lab experiments does each strategy need to identify the hidden
four-factor recipe? Twenty independent runs, each with a different hidden
recipe and a different noise seed.

Writes ``results/01_main_comparison.csv``.
"""

from __future__ import annotations

import csv
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import (  # noqa: E402
    ColonyOracle,
    N_CANDIDATES,
    RECIPE_SIZE,
    adaptive_design,
    historical_protocol,
    weighted_prior,
)

N_RUNS = 20
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

LABELS = {
    "historical": "Historical protocol (2006)",
    "adaptive": "Adaptive design (info gain)",
    "informed": "Adaptive + literature prior",
    "random": "Random pools + Bayes [ablation]",
}


def literature_prior(truth, rng, strength: float = 3.0, n_informed: int = 8):
    """
    A plausible-but-imperfect literature prior.

    Up-weights ``n_informed`` factors: the true ones plus some decoys, since
    published pluripotency candidate lists contain both. Experiment 2 degrades
    this prior deliberately.
    """
    weights = np.ones(N_CANDIDATES)
    informed = set(truth)
    others = [f for f in range(N_CANDIDATES) if f not in informed]
    rng.shuffle(others)
    informed |= set(others[: n_informed - len(informed)])
    for f in informed:
        weights[f] = strength
    return weighted_prior(weights)


def main() -> None:
    results = {k: {"experiments": [], "correct": []} for k in LABELS}

    print(f"\n{'=' * 72}")
    print("  Experiment 1 -- wet-lab experiments needed to find the recipe")
    print(f"{'=' * 72}")
    print(f"  candidates: {N_CANDIDATES}   recipe size: {RECIPE_SIZE}   runs: {N_RUNS}\n")
    print("  run │ historical │ adaptive │ informed │ random")
    print("  ────┼────────────┼──────────┼──────────┼────────")

    for run in range(N_RUNS):
        rng = np.random.default_rng(7000 + run)
        truth = set(
            int(x) for x in rng.choice(N_CANDIDATES, size=RECIPE_SIZE, replace=False)
        )

        oracle = ColonyOracle(truth, seed=run)
        found, n = historical_protocol(oracle)
        results["historical"]["experiments"].append(n)
        results["historical"]["correct"].append(found == truth)

        oracle = ColonyOracle(truth, seed=run)
        found, n, _ = adaptive_design(oracle, seed=run)
        results["adaptive"]["experiments"].append(n)
        results["adaptive"]["correct"].append(found == truth)

        oracle = ColonyOracle(truth, seed=run)
        prior = literature_prior(truth, np.random.default_rng(9000 + run))
        found, n, _ = adaptive_design(oracle, seed=run, prior=prior)
        results["informed"]["experiments"].append(n)
        results["informed"]["correct"].append(found == truth)

        oracle = ColonyOracle(truth, seed=run)
        found, n, _ = adaptive_design(oracle, seed=run, adaptive=False)
        results["random"]["experiments"].append(n)
        results["random"]["correct"].append(found == truth)

        print(
            f"  {run + 1:3d} │ {results['historical']['experiments'][-1]:10d} │"
            f" {results['adaptive']['experiments'][-1]:8d} │"
            f" {results['informed']['experiments'][-1]:8d} │"
            f" {results['random']['experiments'][-1]:6d}"
        )

    baseline = float(np.median(results["historical"]["experiments"]))

    print(f"\n{'=' * 72}")
    print(f"  {'strategy':<34}{'median':>8}{'mean':>8}{'success':>10}{'saving':>9}")
    print("  " + "-" * 68)
    rows = []
    for key in ("historical", "adaptive", "informed", "random"):
        exp = np.array(results[key]["experiments"])
        success = 100.0 * float(np.mean(results[key]["correct"]))
        median = float(np.median(exp))
        saving = "--" if key == "historical" else f"{100 * (1 - median / baseline):.0f}%"
        print(
            f"  {LABELS[key]:<34}{median:>8.0f}{exp.mean():>8.1f}"
            f"{success:>9.0f}%{saving:>9}"
        )
        rows.append(
            {
                "strategy": LABELS[key],
                "median_experiments": median,
                "mean_experiments": round(float(exp.mean()), 2),
                "success_rate_pct": round(success, 1),
                "saving_vs_historical_pct": (
                    0.0 if key == "historical" else round(100 * (1 - median / baseline), 1)
                ),
            }
        )

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "01_main_comparison.csv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  wrote {out.relative_to(RESULTS.parent)}\n")


if __name__ == "__main__":
    main()
