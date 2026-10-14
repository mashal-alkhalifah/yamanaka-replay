#!/usr/bin/env python3
"""
Does the gap structure the oracle assumes exist in measured data?

Run:
    python analysis/norman_audit.py <pseudobulk.csv> <perturbations.csv>

Background
----------
The simulator in this repository models the colony assay as threshold group
testing with a gap of one: a pool containing the whole recipe acts fully, a
pool one factor short acts at a fraction rho of that (calibrated to 0.10 from
the 2006 paper), and a pool two short does nothing. Every result in the study
is conditional on gap structure of roughly that shape existing.

That assumption has never been checked against measured combinatorial
perturbation data, which is the single largest weakness of the work. Norman
et al. (2019) measured 105 single CRISPRa perturbations and 131 pairs in K562
with single-cell readout, which makes the two-factor analogue of the question
answerable: given a measured pair, how much of its effect does one of its two
constituents reproduce alone?

What is measured
----------------
For a perturbation p with pseudobulk profile x_p and control x_0, write
delta_p = x_p - x_0. For each pair (A,B) whose two singles were also measured:

``one_short``   max(|delta_A|, |delta_B|) / |delta_AB|
    The direct empirical analogue of rho. The oracle's calibrated value is
    0.10; a value near 1 means one factor alone nearly reproduces the pair,
    i.e. a much more gradual landscape than the oracle assumes.

``synergy``     |delta_AB| / |delta_A + delta_B|
    Above 1, the pair exceeds the additive prediction; below 1, it falls short
    (buffering). Exactly 1 is additivity.

``residual``    |delta_AB - (delta_A + delta_B)| / |delta_AB|
    How much of the pair's effect the additive model fails to explain.

``r2``          R^2 of the least-squares fit of delta_AB on delta_A, delta_B
    Whether the pair's direction in expression space lies in the plane spanned
    by its constituents, regardless of magnitude.

Noise floor
-----------
Pseudobulk means carry sampling error that inflates every residual, and with
pseudobulk input the per-cell variance is gone. It is still recoverable: the
sampling variance of a residual scales as 1/n_AB + 1/n_A + 1/n_B in the cell
counts, so regressing squared residual norms on that quantity splits the
observed non-additivity into a sampling-noise slope and a biological
intercept. A residual is called real only if it exceeds the fitted noise
level. This avoids the mistake of comparing residuals against an absolute
constant, which conflates measurement error with biology.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np


def load(pseudobulk: pathlib.Path, meta: pathlib.Path):
    """Read both tables with the standard library; no pandas dependency."""
    import csv

    rows, index = [], []
    with pseudobulk.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for row in reader:
            index.append(row[0])
            rows.append([float(v) for v in row[1:]])
    X = np.asarray(rows, dtype=float)

    with meta.open(newline="") as fh:
        info = {r["perturbation"]: r for r in csv.DictReader(fh)}
    return X, index, header[1:], info


def main() -> None:
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    X, index, genes, info = load(pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]))
    pos = {label: i for i, label in enumerate(index)}
    print(f"  {X.shape[0]} perturbations x {X.shape[1]} genes")

    control = next((c for c in ("control", "ctrl", "non-targeting") if c in pos), None)
    if control is None:
        sys.exit("no control row found")
    x0 = X[pos[control]]
    n_cells = {k: float(v["n_cells"]) for k, v in info.items()}

    records = []
    skipped = 0
    for label, row in info.items():
        if row["type"] != "pair" or label not in pos:
            continue
        a, b = row["genes"].split(";")
        if a not in pos or b not in pos:
            skipped += 1
            continue

        d_ab = X[pos[label]] - x0
        d_a = X[pos[a]] - x0
        d_b = X[pos[b]] - x0

        n_ab, n_a, n_b = (np.linalg.norm(v) for v in (d_ab, d_a, d_b))
        if n_ab < 1e-9:
            skipped += 1
            continue

        additive = d_a + d_b
        residual = d_ab - additive
        design = np.column_stack([d_a, d_b])
        coef, *_ = np.linalg.lstsq(design, d_ab, rcond=None)
        fitted = design @ coef
        ss_res = float(np.sum((d_ab - fitted) ** 2))
        ss_tot = float(np.sum(d_ab**2))

        records.append({
            "pair": label,
            "a": a,
            "b": b,
            "one_short": max(n_a, n_b) / n_ab,
            "synergy": n_ab / max(np.linalg.norm(additive), 1e-12),
            "residual": float(np.linalg.norm(residual)) / n_ab,
            "r2": 1.0 - ss_res / max(ss_tot, 1e-12),
            "cos_additive": float(d_ab @ additive /
                                  max(n_ab * np.linalg.norm(additive), 1e-12)),
            "coef_sum": float(coef.sum()),
            "norm_pair": n_ab,
            "inv_n": 1.0 / n_cells[label] + 1.0 / n_cells[a] + 1.0 / n_cells[b],
            "resid_sq": float(np.sum(residual**2)),
        })

    print(f"  usable pairs: {len(records)}   skipped: {skipped}\n")
    if not records:
        sys.exit("nothing to analyse")

    def col(key):
        return np.array([r[key] for r in records], dtype=float)

    # ---- noise decomposition -------------------------------------------------
    inv_n, resid_sq = col("inv_n"), col("resid_sq")
    slope, intercept = np.polyfit(inv_n, resid_sq, 1)
    noise_sq = np.clip(slope, 0, None) * inv_n
    biological_sq = np.clip(resid_sq - noise_sq, 0, None)
    share = biological_sq.sum() / max(resid_sq.sum(), 1e-12)

    print("  ------------------------------------------------------------------")
    print("  Noise decomposition of the additivity residual")
    print("  ------------------------------------------------------------------")
    print(f"  sampling-noise slope      {slope:12.2f}  (per unit of 1/n)")
    print(f"  intercept                 {intercept:12.2f}")
    print(f"  share of residual that survives the noise model: {100 * share:5.1f}%")
    detected = int(np.sum(biological_sq > noise_sq))
    print(f"  pairs whose residual exceeds its own noise level: "
          f"{detected}/{len(records)}  ({100 * detected / len(records):.0f}%)")

    # ---- the quantity the oracle assumes ------------------------------------
    one_short = col("one_short")
    print("\n  ------------------------------------------------------------------")
    print("  How much of a pair's effect does ONE of its factors reproduce?")
    print("  (the empirical analogue of the oracle's rho = 0.10)")
    print("  ------------------------------------------------------------------")
    for q in (5, 25, 50, 75, 95):
        print(f"  {q:>3d}th percentile   {np.percentile(one_short, q):6.3f}")
    print(f"  mean              {one_short.mean():6.3f}")
    print(f"\n  fraction of pairs below 0.20 (oracle-like, sharp threshold): "
          f"{100 * np.mean(one_short < 0.20):5.1f}%")
    print(f"  fraction between 0.20 and 0.80 (graded)                     : "
          f"{100 * np.mean((one_short >= 0.20) & (one_short <= 0.80)):5.1f}%")
    print(f"  fraction above 0.80 (one factor nearly suffices)            : "
          f"{100 * np.mean(one_short > 0.80):5.1f}%")

    # ---- additivity ---------------------------------------------------------
    syn, res, r2 = col("synergy"), col("residual"), col("r2")
    print("\n  ------------------------------------------------------------------")
    print("  Departure from additivity")
    print("  ------------------------------------------------------------------")
    print(f"  |pair| / |sum of singles|   median {np.median(syn):.3f}   "
          f"IQR {np.percentile(syn, 25):.3f}-{np.percentile(syn, 75):.3f}")
    print(f"  residual / |pair|           median {np.median(res):.3f}   "
          f"IQR {np.percentile(res, 25):.3f}-{np.percentile(res, 75):.3f}")
    print(f"  R^2 of the two-single fit   median {np.median(r2):.3f}   "
          f"IQR {np.percentile(r2, 25):.3f}-{np.percentile(r2, 75):.3f}")
    print(f"  super-additive (ratio > 1.1): {100 * np.mean(syn > 1.1):5.1f}%")
    print(f"  near-additive  (0.9 - 1.1)  : "
          f"{100 * np.mean((syn >= 0.9) & (syn <= 1.1)):5.1f}%")
    print(f"  sub-additive   (ratio < 0.9): {100 * np.mean(syn < 0.9):5.1f}%")

    # ---- the strongest and most additive pairs ------------------------------
    order = np.argsort(-syn)
    print("\n  most super-additive pairs")
    for i in order[:8]:
        r = records[i]
        print(f"    {r['pair']:<22} synergy {r['synergy']:5.2f}   "
              f"one-short {r['one_short']:5.2f}   R2 {r['r2']:5.2f}")
    print("\n  most sub-additive pairs")
    for i in order[-8:]:
        r = records[i]
        print(f"    {r['pair']:<22} synergy {r['synergy']:5.2f}   "
              f"one-short {r['one_short']:5.2f}   R2 {r['r2']:5.2f}")

    out = pathlib.Path(__file__).resolve().parents[1] / "results" / "norman_audit.csv"
    out.parent.mkdir(exist_ok=True)
    import csv as _csv

    fields = [k for k in records[0] if k not in ("resid_sq",)]
    with out.open("w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)
    print(f"\n  wrote {out.relative_to(out.parents[1])}")


if __name__ == "__main__":
    main()
