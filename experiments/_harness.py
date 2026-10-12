"""
Shared harness for the large-scale campaign (experiments 05-08).

Experiments 01-04 are the small, fast, reviewer-facing demonstrations. The
scripts in this family instead run hundreds of independent replicates and write
**one row per run** rather than a summary table, so that the aggregation and the
statistics live in ``analysis/`` and can be re-done without re-running anything.

Three things in here matter:

``run_one``
    Executes a single (strategy, condition, seed) cell and returns a flat
    record. Every source of randomness is derived from the cell's own seeds, so
    a row is reproducible in isolation.

``run_cells``
    Maps ``run_one`` over a list of cells across worker processes. Each worker
    re-imports the library, so no state is shared and no lock is needed.

``set_candidate_count``
    Rebuilds the hypothesis space for a different number of candidate factors.
    ``N_CANDIDATES`` is a module constant bound at import time in three places;
    changing it safely means rebinding all three plus the derived hypothesis
    tables. Used only by experiment 08.
"""

from __future__ import annotations

import csv
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from yamanaka_replay import (  # noqa: E402
    ColonyOracle,
    MisspecifiedOracle,
    adaptive_design,
    historical_protocol,
    weighted_prior,
)
from yamanaka_replay import inference, oracle as oracle_mod, strategies  # noqa: E402

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

#: (n_hypotheses, recipe_size) factor-index array per candidate count, kept so
#: that priors can be built with array arithmetic instead of a Python loop.
_HYPOTHESIS_INDEX: dict[int, np.ndarray] = {}


# ------------------------------------------------------------ hypothesis space


def set_candidate_count(n_candidates: int) -> None:
    """
    Rebind the candidate count across the library and rebuild the hypothesis
    tables. Call once per worker, before any run at that size.

    Exists because ``N_CANDIDATES`` is imported by value into ``inference`` and
    ``strategies``; setting only ``oracle.N_CANDIDATES`` would leave the other
    two stale and silently produce a mismatched hypothesis space.
    """
    recipe_size = oracle_mod.RECIPE_SIZE
    oracle_mod.N_CANDIDATES = n_candidates
    inference.N_CANDIDATES = n_candidates
    strategies.N_CANDIDATES = n_candidates

    # Vectorised: at 60 candidates there are 487,635 hypotheses, and filling the
    # membership matrix row by row in Python costs more than a whole search run.
    index = np.fromiter(
        (g for c in combinations(range(n_candidates), recipe_size) for g in c),
        dtype=np.int16,
    ).reshape(-1, recipe_size)
    mask = np.zeros((index.shape[0], n_candidates), dtype=bool)
    mask[np.arange(index.shape[0])[:, None], index] = True
    hypotheses = [frozenset(row.tolist()) for row in index]

    inference.HYPOTHESES = hypotheses
    inference.HYPOTHESIS_MASK = mask
    strategies.HYPOTHESES = hypotheses
    _HYPOTHESIS_INDEX[n_candidates] = index


# -------------------------------------------------------------------- priors


def literature_prior(
    truth,
    rng,
    n_candidates: int,
    strength: float = 3.0,
    n_informed: int = 8,
    n_dropped: int = 0,
):
    """
    A plausible-but-imperfect literature prior, with controllable quality.

    ``n_dropped`` true factors are deliberately left out of the up-weighted set,
    modelling a published candidate list that misses genes. The remaining slots
    are filled with decoys, because real candidate lists contain both.
    """
    weights = np.ones(n_candidates)
    truth = list(truth)
    rng.shuffle(truth)
    informed = set(truth[n_dropped:])
    others = [f for f in range(n_candidates) if f not in set(truth)]
    rng.shuffle(others)
    informed |= set(others[: max(0, n_informed - len(informed))])
    for f in informed:
        weights[f] = strength

    index = _HYPOTHESIS_INDEX.get(n_candidates)
    if index is None:  # 24 candidates: the library's own tables are still live
        return weighted_prior(weights)
    prior = weights[index].prod(axis=1)
    return prior / prior.sum()


# ----------------------------------------------------------------- one cell


def run_one(cell: dict) -> dict:
    """
    Execute one replicate.

    Recognised keys, all optional except ``strategy`` and ``run``:

    ===================  ======================================================
    ``strategy``         ``historical`` | ``adaptive`` | ``informed`` | ``random``
    ``run``              replicate index; seeds the hidden recipe and the noise
    ``n_candidates``     candidate-pool size (default 24)
    ``confidence``       stopping threshold (default 0.95)
    ``max_experiments``  budget (default 60)
    ``partial_3``        oracle's three-of-four efficiency
    ``assumed_partial_3``  what the *searcher* believes it to be; set equal to
                         ``partial_3`` for the matched arm, leave ``None`` to
                         keep the library default and so be misspecified
    ``base``, ``decay``, ``extra_required``  further oracle overrides
    ``prior_strength``, ``prior_informed``, ``prior_dropped``  prior quality
    ``n_random_pools``   proposals scored per round
    ===================  ======================================================
    """
    strategy = cell["strategy"]
    run = int(cell["run"])
    n_candidates = int(cell.get("n_candidates", 24))
    if n_candidates != inference.N_CANDIDATES:
        set_candidate_count(n_candidates)

    recipe_size = oracle_mod.RECIPE_SIZE
    rng = np.random.default_rng(7_000_000 + 1000 * run + n_candidates)
    truth = set(
        int(x) for x in rng.choice(n_candidates, size=recipe_size, replace=False)
    )

    # The searcher's assumed three-of-four efficiency is read off the class
    # attribute by inference.expected_rate, so the matched arm sets it here.
    assumed = cell.get("assumed_partial_3")
    default_partial = 0.10
    ColonyOracle.PARTIAL_3 = default_partial if assumed is None else float(assumed)

    overrides = {
        k: cell[k]
        for k in ("partial_3", "base", "decay", "extra_required")
        if cell.get(k) is not None
    }
    if "extra_required" in overrides:
        # A fifth essential factor, drawn from outside the four the searcher
        # could possibly name.
        outside = [f for f in range(n_candidates) if f not in truth]
        overrides["extra_required"] = int(outside[run % len(outside)])

    if overrides:
        lab = MisspecifiedOracle(truth, seed=run, **overrides)
    else:
        lab = ColonyOracle(truth, seed=run)
        lab.PARTIAL_3 = default_partial  # instance-level, immune to the class poke

    record = {
        **{k: v for k, v in cell.items() if k != "strategy"},
        "strategy": strategy,
        "n_candidates": n_candidates,
    }

    t0 = time.perf_counter()
    if strategy == "historical":
        found, n = historical_protocol(lab)
        confidence = float("nan")
    else:
        prior = None
        if strategy == "informed":
            prior = literature_prior(
                truth,
                np.random.default_rng(9_000_000 + 1000 * run + n_candidates),
                n_candidates,
                strength=float(cell.get("prior_strength", 3.0)),
                n_informed=int(cell.get("prior_informed", 8)),
                n_dropped=int(cell.get("prior_dropped", 0)),
            )
        found, n, confidence = adaptive_design(
            lab,
            seed=run,
            max_experiments=int(cell.get("max_experiments", 60)),
            confidence=float(cell.get("confidence", 0.95)),
            prior=prior,
            n_random_pools=int(cell.get("n_random_pools", 260)),
            adaptive=(strategy != "random"),
        )

    elapsed = time.perf_counter() - t0
    pools = [p for p, _ in lab.log]
    record.update(
        {
            "n_experiments": n,
            "correct": int(found == truth),
            "n_recovered": len(found & truth),
            "final_confidence": confidence,
            # The historical protocol has no budget: it always runs its 61 steps.
            "hit_budget": (
                0
                if strategy == "historical"
                else int(n >= int(cell.get("max_experiments", 60)))
            ),
            "mean_pool_size": round(float(np.mean([len(p) for p in pools])), 2),
            "n_zero_colony": sum(1 for _, c in lab.log if c == 0),
            "max_colonies": max((c for _, c in lab.log), default=0),
            "elapsed_s": round(elapsed, 4),
        }
    )
    return record


# ------------------------------------------------------------------ dispatch


def run_cells(cells: list[dict], workers: int = 2, label: str = "") -> list[dict]:
    """Map ``run_one`` across processes, printing progress every 5%."""
    total = len(cells)
    step = max(1, total // 20)
    out: list[dict] = []
    print(f"  {label}: {total} runs on {workers} workers")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, rec in enumerate(pool.map(run_one, cells, chunksize=4), start=1):
            out.append(rec)
            if i % step == 0 or i == total:
                print(f"    {i:5d}/{total}  ({100 * i / total:3.0f}%)", flush=True)
    return out


def write_rows(rows: list[dict], filename: str) -> pathlib.Path:
    """Write per-run records, with the union of all keys as the header."""
    RESULTS.mkdir(exist_ok=True)
    fields: list[str] = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    path = RESULTS / filename
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote {path.relative_to(RESULTS.parent)}  ({len(rows)} rows)")
    return path
