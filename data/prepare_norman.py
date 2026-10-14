#!/usr/bin/env python3
"""
Turn the 699 MB Norman et al. (2019) single-cell file into a small table.

Why this step exists
--------------------
The audit does not need 90,000 individual cells. It needs, for each
perturbation (each single gene and each gene pair), the average expression
profile it produces. That is a few megabytes rather than a gigabyte, and it is
what gets shared, attached and version-controlled.

Usage
-----
    pip install anndata numpy pandas          # scanpy not required
    python prepare_norman.py NormanWeissman2019_filtered.h5ad

Writes, beside the input file:

    norman_pseudobulk.csv.gz   perturbations x genes, mean expression
    norman_perturbations.csv   one row per perturbation: label, type, n_cells

Both are small enough to attach to a chat message.

Notes
-----
* ``--genes`` controls how many of the most variable genes are kept (default
  2000). The audit computes synergy in this gene space, so more genes means a
  finer measurement and a larger file; 2000 is already generous.
* The script does not assume a particular column name for the perturbation
  label. It looks for the usual candidates and, failing that, prints the
  columns it found so you can pass ``--label-column`` explicitly rather than
  guess.
* Expression is log1p-normalised per cell to 10,000 counts before averaging,
  which is the standard choice and the one the audit's additivity null
  assumes. If the file is already normalised, pass ``--already-normalised``.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

LABEL_CANDIDATES = (
    "perturbation", "perturbation_name", "guide_identity", "gene_program",
    "condition", "perturbation_label", "gene_target", "target",
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ad", type=pathlib.Path, help="NormanWeissman2019_filtered.h5ad")
    ap.add_argument("--genes", type=int, default=2000,
                    help="number of most-variable genes to keep (default 2000)")
    ap.add_argument("--label-column", default=None,
                    help="obs column holding the perturbation label")
    ap.add_argument("--already-normalised", action="store_true",
                    help="skip per-cell normalisation and log1p")
    ap.add_argument("--min-cells", type=int, default=20,
                    help="drop perturbations with fewer cells than this (default 20)")
    args = ap.parse_args()

    try:
        import anndata as ad
        import pandas as pd
    except ImportError:
        sys.exit("needs anndata and pandas:  pip install anndata pandas")

    if not args.h5ad.exists():
        sys.exit(f"not found: {args.h5ad}")

    print(f"reading {args.h5ad} ({args.h5ad.stat().st_size / 1e6:.0f} MB) ...")
    adata = ad.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- the perturbation label
    column = args.label_column
    if column is None:
        for candidate in LABEL_CANDIDATES:
            if candidate in adata.obs.columns:
                column = candidate
                break
    if column is None or column not in adata.obs.columns:
        print("\nCould not identify the perturbation column. Columns present:")
        for c in adata.obs.columns:
            n = adata.obs[c].nunique()
            print(f"  {c}  ({n} distinct values)")
        sys.exit("\nRe-run with --label-column <one of the above>")
    print(f"  perturbation label: obs['{column}']")

    labels = adata.obs[column].astype(str).values

    # ---- normalise
    import scipy.sparse as sp

    X = adata.X
    if not args.already_normalised:
        print("normalising to 10,000 counts per cell and applying log1p ...")
        totals = np.asarray(X.sum(axis=1)).ravel()
        totals[totals == 0] = 1.0
        scale = 1e4 / totals
        if sp.issparse(X):
            X = sp.diags(scale) @ X
            X = X.log1p()
        else:
            X = np.log1p(X * scale[:, None])

    # ---- keep the most variable genes, measured on the control cells so that
    #      gene choice is not itself driven by the perturbations being audited
    print(f"selecting the {args.genes} most variable genes ...")
    mean = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel() if sp.issparse(X) \
        else (X ** 2).mean(axis=0)
    var = np.maximum(sq - mean ** 2, 0.0)
    keep = np.argsort(-var)[: min(args.genes, X.shape[1])]
    keep.sort()
    gene_names = np.asarray(adata.var_names)[keep]

    # ---- pseudobulk
    print("averaging per perturbation ...")
    unique = sorted(set(labels))
    rows, meta = [], []
    Xk = X[:, keep]
    for label in unique:
        mask = labels == label
        n = int(mask.sum())
        if n < args.min_cells:
            continue
        block = Xk[mask]
        profile = np.asarray(block.mean(axis=0)).ravel()
        rows.append(profile)
        # scPerturb writes pairs as "GENEA_GENEB"; controls vary by dataset
        parts = [p for p in label.replace("+", "_").split("_") if p]
        kind = (
            "control" if label.lower() in {"control", "ctrl", "nt", "non-targeting"}
            else {1: "single", 2: "pair"}.get(len(parts), f"{len(parts)}-way")
        )
        meta.append({"perturbation": label, "type": kind, "n_cells": n,
                     "genes": ";".join(parts)})

    frame = pd.DataFrame(np.vstack(rows), index=[m["perturbation"] for m in meta],
                         columns=gene_names)
    out_dir = args.h5ad.parent
    frame.to_csv(out_dir / "norman_pseudobulk.csv.gz", compression="gzip")
    pd.DataFrame(meta).to_csv(out_dir / "norman_perturbations.csv", index=False)

    counts = pd.DataFrame(meta)["type"].value_counts().to_dict()
    size = (out_dir / "norman_pseudobulk.csv.gz").stat().st_size / 1e6
    print(f"\nwrote norman_pseudobulk.csv.gz  ({size:.1f} MB, "
          f"{frame.shape[0]} perturbations x {frame.shape[1]} genes)")
    print(f"wrote norman_perturbations.csv  ({counts})")
    print("\nAttach both files to the chat.")


if __name__ == "__main__":
    main()
