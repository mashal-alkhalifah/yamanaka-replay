#!/usr/bin/env python3
"""
Experiment 4 -- a single run, experiment by experiment.

Prints every pool the designer chooses, the colony count it observes, and how
its belief moves. Useful for seeing *why* the saving happens:

* zero-colony pools are not failures -- each one eliminates every hypothesis
  it contains, which is classical group testing;
* the partial three-factor signal is the thread the search pulls on;
* the designer never wastes an experiment on a single factor, which is where
  the historical protocol spent 24 of its 61.

Usage::

    python experiments/04_trace_single_run.py [run_index]
"""

from __future__ import annotations

import pathlib
import platform
import os
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import (  # noqa: E402
    HYPOTHESES,
    ColonyOracle,
    N_CANDIDATES,
    RECIPE_SIZE,
    bayes_update,
    expected_information_gain,
    marginal_factor_probability,
    uniform_prior,
)
from yamanaka_replay.strategies import _candidate_pools  # noqa: E402

RUN = int(sys.argv[1]) if len(sys.argv) > 1 else 3
CONFIDENCE = 0.95
MAX_STEPS = 30


def main() -> None:
    rng = np.random.default_rng(7000 + RUN)
    truth = set(
        int(x) for x in rng.choice(N_CANDIDATES, size=RECIPE_SIZE, replace=False)
    )
    names = dict(zip(sorted(truth), ["OCT4", "SOX2", "KLF4", "MYC"]))
    name_of = lambda g: names.get(g, f"F{g:02d}")  # noqa: E731

    print("=" * 78)
    print(f"  Trace of a single adaptive run (run index {RUN})")
    print("=" * 78)
    print(
        "  hidden recipe (never shown to the searcher): "
        + ", ".join(f"{name_of(g)}=F{g:02d}" for g in sorted(truth))
    )
    print(f"  hypothesis space: {len(HYPOTHESES):,} possible 4-subsets of {N_CANDIDATES}")
    print(f"  initial belief:   {1 / len(HYPOTHESES):.6f} per hypothesis (uniform)\n")

    oracle = ColonyOracle(truth, seed=RUN)
    posterior = uniform_prior()
    pool_rng = np.random.default_rng(RUN)
    started = time.time()

    print("  exp │ size │ pool tested                          │ colonies │ "
          "best hypothesis     │ conf")
    print("  ────┼──────┼──────────────────────────────────────┼──────────┼"
          "─────────────────────┼──────")

    step = 0
    while posterior.max() < CONFIDENCE and step < MAX_STEPS:
        step += 1
        pools = _candidate_pools(pool_rng, posterior, 260)
        pool = max(pools, key=lambda p: expected_information_gain(posterior, p))
        observed = oracle.run(pool)
        posterior = bayes_update(posterior, pool, observed)

        best = HYPOTHESES[int(np.argmax(posterior))]
        shown = ",".join(f"{g:02d}" for g in pool)
        if len(shown) > 36:
            shown = shown[:33] + "..."
        print(
            f"  {step:3d} │ {len(pool):4d} │ {shown:<36} │ {observed:5d}    │ "
            f"{'+'.join(name_of(g) for g in sorted(best)):<19} │ {posterior.max():.3f}"
        )

    elapsed = time.time() - started
    final = set(HYPOTHESES[int(np.argmax(posterior))])

    print("\n  " + "-" * 74)
    print(f"  recovered recipe  : {'+'.join(name_of(g) for g in sorted(final))}")
    print(f"  correct           : {'yes' if final == truth else 'NO'}")
    print(f"  wet-lab experiments: {oracle.n_experiments}   (historical protocol: 61)")
    print(f"  final confidence  : {posterior.max():.1%}")
    print(f"  compute time      : {elapsed:.2f} s")

    print("\n  marginal P(factor in recipe), top 8:")
    marginal = marginal_factor_probability(posterior)
    for g in np.argsort(-marginal)[:8]:
        bar = "#" * int(round(marginal[g] * 30))
        flag = "  <- truly in recipe" if g in truth else ""
        print(f"    F{g:02d} {name_of(g):<6} {marginal[g]:6.1%} {bar}{flag}")

    print("\n  full experiment log:")
    for i, (pool, count) in enumerate(oracle.log, 1):
        print(
            f"    {i:2d}. |pool|={len(pool):2d}  overlap with recipe="
            f"{len(set(pool) & truth)}  ->  {count} colonies"
        )

    print("\n" + "=" * 78)
    print("  Runtime environment")
    print("=" * 78)
    print(f"  os      : {platform.system()} {platform.release()}")
    print(f"  cpu     : {platform.machine()}, {os.cpu_count()} cores")
    print(f"  python  : {platform.python_version()}, numpy {np.__version__}")
    print("  gpu     : none -- CPU only")
    print("  network : none -- no data downloaded, no API calls\n")


if __name__ == "__main__":
    main()
