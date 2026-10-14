#!/usr/bin/env python3
"""
Does the gap structure the oracle assumes exist in measured data?

    python analysis/real_data_audit.py norman   <pseudobulk.csv> <perturbations.csv>
    python analysis/real_data_audit.py tfatlas  <pseudobulk.csv> <perturbations.csv>

The oracle in this repository models the assay as threshold group testing with a
gap of one: a pool holding the whole recipe acts fully, a pool one factor short
acts at a fraction rho of that, calibrated to 0.10 from the 2006 paper. Every
result in the study is conditional on gap structure of that shape. This script
checks it against measured combinatorial perturbation data.

The quantity measured, for every combination whose one-short constituents were
also profiled:

    one_short = max over constituents of  ||delta_constituent|| / ||delta_combination||

with delta taken against the unperturbed reference. The oracle's calibrated
value is 0.10. A value near 1 means one factor alone nearly reproduces the
combination, i.e. a far more gradual landscape than the oracle assumes.

Two datasets, two designs
-------------------------
``norman``   Norman et al. (2019), K562, CRISPRa. 105 single perturbations, 131
             pairs, one explicit control. Singles and the control are read
             directly.
``tfatlas``  Joung et al. (2023) combinatorial subset, human embryonic stem
             cells, TF overexpression. 10 TFs plus GFP as an inert partner, so
             ``X+GFP`` is the single for X, and 40 TF-TF pairs have both singles.
             **There is no GFP-only condition**, so the unperturbed reference is
             not measured. Rather than pick one substitute, the script reports
             three, including a least-squares fit that is deliberately biased
             toward additivity; a conclusion is only taken as real if it holds
             under all three.

Labels in that dataset are alphabetical (``GFP+HNF4A``, not ``HNF4A+GFP``), so
every lookup tries both orders. Writes ``results/real_data_audit_<dataset>.csv``.
"""

from __future__ import annotations

import csv
import itertools
import pathlib
import sys

import numpy as np

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"


def read_table(pseudobulk: pathlib.Path, meta: pathlib.Path):
    rows, index = [], []
    with pseudobulk.open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            index.append(row[0])
            rows.append([float(v) for v in row[1:]])
    X = np.asarray(rows, dtype=float)
    info = list(csv.DictReader(meta.open(newline="")))
    return X, {label: i for i, label in enumerate(index)}, info


def summarise(name: str, one_short: np.ndarray, synergy: np.ndarray, r2: np.ndarray) -> dict:
    print(f"\n  {name}")
    print(f"    one-short     median {np.median(one_short):.2f}"
          f"   IQR {np.percentile(one_short, 25):.2f}-{np.percentile(one_short, 75):.2f}"
          f"   below 0.20: {100 * np.mean(one_short < 0.20):.0f}%"
          f"   below 0.35: {100 * np.mean(one_short < 0.35):.0f}%")
    print(f"    |comb|/|sum|  median {np.median(synergy):.2f}"
          f"   super-additive {100 * np.mean(synergy > 1.1):.0f}%"
          f"   sub-additive {100 * np.mean(synergy < 0.9):.0f}%")
    print(f"    R^2 on its own constituents   median {np.median(r2):.2f}")
    return {
        "reference": name,
        "n": int(one_short.size),
        "one_short_median": round(float(np.median(one_short)), 3),
        "one_short_q25": round(float(np.percentile(one_short, 25)), 3),
        "one_short_q75": round(float(np.percentile(one_short, 75)), 3),
        "frac_below_0.20": round(float(np.mean(one_short < 0.20)), 3),
        "synergy_median": round(float(np.median(synergy)), 3),
        "r2_median": round(float(np.median(r2)), 3),
    }


def metrics(combination: np.ndarray, parts: list[np.ndarray]) -> tuple[float, float, float]:
    """one-short ratio, magnitude ratio against the additive sum, and R^2."""
    n_comb = float(np.linalg.norm(combination))
    if n_comb < 1e-9:
        return float("nan"), float("nan"), float("nan")
    one_short = max(float(np.linalg.norm(p)) for p in parts) / n_comb
    additive = np.sum(parts, axis=0)
    synergy = n_comb / max(float(np.linalg.norm(additive)), 1e-12)
    design = np.column_stack(parts)
    coef, *_ = np.linalg.lstsq(design, combination, rcond=None)
    residual = combination - design @ coef
    r2 = 1.0 - float(np.sum(residual**2)) / max(float(np.sum(combination**2)), 1e-12)
    return one_short, synergy, r2


# ------------------------------------------------------------------- norman


def audit_norman(X, pos, info) -> list[dict]:
    control = next((c for c in ("control", "ctrl", "non-targeting") if c in pos), None)
    if control is None:
        sys.exit("no control row found")
    x0 = X[pos[control]]

    one, syn, r2s, records = [], [], [], []
    for row in info:
        if row["type"] != "pair" or row["perturbation"] not in pos:
            continue
        a, b = row["genes"].split(";")
        if a not in pos or b not in pos:
            continue
        d_ab = X[pos[row["perturbation"]]] - x0
        parts = [X[pos[a]] - x0, X[pos[b]] - x0]
        o, s, r = metrics(d_ab, parts)
        one.append(o), syn.append(s), r2s.append(r)
        records.append({"combination": row["perturbation"], "size": 2,
                        "one_short": round(o, 4), "synergy": round(s, 4), "r2": round(r, 4)})
    print(f"  Norman et al. 2019 -- K562, CRISPRa -- {len(records)} pairs with both singles")
    summary = [summarise("measured control", np.array(one), np.array(syn), np.array(r2s))]
    return records, summary


# ------------------------------------------------------------------ tfatlas


def audit_tfatlas(X, pos, info) -> list[dict]:
    def find(*names):
        for order in itertools.permutations(names):
            key = "+".join(order)
            if key in pos:
                return pos[key]
        return None

    factors = sorted({g for row in info for g in row["genes"].split(";")} - {"GFP"})
    single = {t: find(t, "GFP") for t in factors}
    single = {t: i for t, i in single.items() if i is not None}

    pairs = [(row["perturbation"], *row["genes"].split(";")) for row in info
             if row["type"] == "pair" and "GFP" not in row["genes"].split(";")]
    pairs = [(p, a, b) for p, a, b in pairs
             if a in single and b in single and p in pos]
    print(f"  Joung et al. 2023 -- hESC, TF overexpression -- "
          f"{len(single)}/{len(factors)} singles, {len(pairs)} TF-TF pairs")

    P = np.array([X[pos[p]] for p, _, _ in pairs])
    A = np.array([X[single[a]] for _, a, _ in pairs])
    B = np.array([X[single[b]] for _, _, b in pairs])

    references = {
        "least-squares fit (biased toward additivity)": (A + B - P).mean(axis=0),
        "centroid of all conditions": X.mean(axis=0),
        "centroid of the single-TF conditions": X[[single[t] for t in single]].mean(axis=0),
    }

    records, summary = [], []
    for name, x0 in references.items():
        one, syn, r2s = [], [], []
        for i, (p, a, b) in enumerate(pairs):
            o, s, r = metrics(P[i] - x0, [A[i] - x0, B[i] - x0])
            one.append(o), syn.append(s), r2s.append(r)
            if name.startswith("least-squares"):
                records.append({"combination": p, "size": 2, "one_short": round(o, 4),
                                "synergy": round(s, 4), "r2": round(r, 4)})
        summary.append(summarise(name, np.array(one), np.array(syn), np.array(r2s)))

    # The two three-way combinations: one short of a triple is a measured pair.
    x0 = references["least-squares fit (biased toward additivity)"]
    triples = [row for row in info if row["type"] == "3-way" and row["perturbation"] in pos]
    if triples:
        print("\n  three-way combinations (one short of a triple is a measured pair)")
    for row in triples:
        genes = row["genes"].split(";")
        n_triple = float(np.linalg.norm(X[pos[row["perturbation"]]] - x0))
        got = [(f"{a}+{b}", float(np.linalg.norm(X[find(a, b)] - x0)) / n_triple)
               for a, b in itertools.combinations(genes, 2) if find(a, b) is not None]
        if not got:
            continue
        best = max(v for _, v in got)
        print(f"    {row['perturbation']}: "
              + "  ".join(f"{s} {v:.2f}" for s, v in got)
              + f"   -> one-short {best:.2f}")
        records.append({"combination": row["perturbation"], "size": 3,
                        "one_short": round(best, 4), "synergy": "", "r2": ""})
    return records, summary


def main() -> None:
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    dataset = sys.argv[1]
    X, pos, info = read_table(pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3]))
    print(f"  {X.shape[0]} conditions x {X.shape[1]} genes\n")

    if dataset == "norman":
        records, summary = audit_norman(X, pos, info)
    elif dataset == "tfatlas":
        records, summary = audit_tfatlas(X, pos, info)
    else:
        sys.exit("dataset must be 'norman' or 'tfatlas'")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / f"real_data_audit_{dataset}.csv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    out_sum = RESULTS / f"real_data_summary_{dataset}.csv"
    with out_sum.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(f"\n  wrote {out.name} and {out_sum.name}")
    print(f"  the oracle assumes one-short = 0.10")


if __name__ == "__main__":
    main()
