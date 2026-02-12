#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
8_Task3_PCA_Baseline_fixed.py

Generate PCA baseline embeddings for Chem2Gen Task3 (K562 matched evaluation set).

This script is aligned to the SAME data-loading interface used by your uploaded FM scripts:
- It requires `Task3_Functions.Chem2GenEvalDataLoader`.
- It defaults to `Task3_Functions.EVAL_DATA_DIR` and `Task3_Functions.OUTPUT_DIR`.
- It will hard-fail if --eval_dir does not contain the expected Evaluation_Set_K562 files
  to prevent accidentally loading a different dataset (the most common reason for N mismatch).

Key choices
-----------
- CRISPR and Drug are concatenated and jointly fitted to ONE PCA basis (coordinates aligned).
- Preprocessing follows Scanpy-style pipeline on concatenated matrix:
    1) sc.pp.normalize_total(target_sum=1e4)
    2) sc.pp.log1p
    3) sc.pp.scale(max_value=10)   (disable by setting --scale_max_value <= 0)
- PCA is fitted once up to max(dims), then sliced to save multiple K (e.g., 50/100/200).

Outputs
-------
Saved to OUTPUT_DIR (or --out_dir) as:
- PCA{K}_CRISPR_emb.npy
- PCA{K}_Drug_emb.npy
Plus: PCA_baseline_provenance.json

Usage
-----
python 8_Task3_PCA_Baseline_fixed.py
python 8_Task3_PCA_Baseline_fixed.py --dims 50,100,200
python 8_Task3_PCA_Baseline_fixed.py --eval_dir /mnt/NAS_21T/.../Evaluation_Set_K562 --out_dir /mnt/NAS_21T/.../Model_Evaluation_Results
"""

from __future__ import annotations

import os
import argparse
from typing import List

import numpy as np
import torch

try:
    import anndata as ad
    import scanpy as sc
except Exception as e:  # pragma: no cover
    raise ImportError(
        "This script requires 'scanpy' and 'anndata'. Please install them in your environment.\n"
        "Example (conda): conda install -c conda-forge scanpy anndata\n"
        f"Original import error: {e}"
    )

# ---- Required: use the same loader interface as other FM scripts ----
try:
    from Task3_Functions import EVAL_DATA_DIR, OUTPUT_DIR, Chem2GenEvalDataLoader
except Exception as e:  # pragma: no cover
    raise ImportError(
        "Failed to import Task3_Functions (required for consistent data interface).\n"
        "Please run this script from the same project environment where Task3_Functions.py is available,\n"
        "or fix PYTHONPATH so that Task3_Functions can be imported.\n"
        f"Original error: {e}"
    )


def _parse_dims(dims_str: str) -> List[int]:
    dims: List[int] = []
    for tok in dims_str.split(","):
        tok = tok.strip()
        if not tok:
            continue
        v = int(tok)
        if v <= 0:
            raise ValueError(f"Invalid PCA dim: {v}")
        dims.append(v)
    if not dims:
        raise ValueError("No PCA dims provided.")
    return sorted(set(dims))


def _validate_eval_dir(eval_dir: str) -> None:
    """Hard validation to prevent accidentally loading the wrong dataset directory."""
    required = [
        "shared_var_names.csv",
        "CRISPR_counts.pt",
        "CRISPR_meta.csv",
        "Drug_counts.pt",
        "Drug_meta.csv",
        "Common_Targets_K562.csv",
    ]
    missing = [f for f in required if not os.path.exists(os.path.join(eval_dir, f))]
    if missing:
        raise FileNotFoundError(
            "eval_dir does not look like a valid Evaluation_Set_K562 directory.\n"
            f"  eval_dir = {eval_dir}\n"
            f"Missing files: {missing}\n"
            "Tip: Use the same eval_dir as your other FM scripts (typically Task3_Functions.EVAL_DATA_DIR)."
        )


def main() -> None:
    ap = argparse.ArgumentParser("Chem2Gen Task3 PCA baseline (Scanpy preprocessing + joint PCA).")

    ap.add_argument(
        "--eval_dir",
        type=str,
        default=EVAL_DATA_DIR,
        help="Evaluation set directory (default: Task3_Functions.EVAL_DATA_DIR).",
    )
    ap.add_argument(
        "--out_dir",
        type=str,
        default=OUTPUT_DIR,
        help="Output directory (default: Task3_Functions.OUTPUT_DIR).",
    )
    ap.add_argument(
        "--dims",
        type=str,
        default="50,100,200",
        help="Comma-separated PCA dimensions to save (e.g., '50,100,200').",
    )
    ap.add_argument(
        "--target_sum",
        type=float,
        default=1e4,
        help="sc.pp.normalize_total target_sum (default 1e4).",
    )
    ap.add_argument(
        "--scale_max_value",
        type=float,
        default=10.0,
        help="sc.pp.scale(max_value=...). Set <=0 to disable scaling.",
    )
    ap.add_argument(
        "--svd_solver",
        type=str,
        default="arpack",
        choices=["arpack", "randomized", "auto", "full"],
        help="PCA SVD solver (default arpack).",
    )
    ap.add_argument(
        "--random_state",
        type=int,
        default=0,
        help="Random seed (used by some solvers).",
    )

    args = ap.parse_args()
    dims = _parse_dims(args.dims)
    max_comps = max(dims)

    _validate_eval_dir(args.eval_dir)
    os.makedirs(args.out_dir, exist_ok=True)

    # ---- Load via the same loader as other FM scripts ----
    loader = Chem2GenEvalDataLoader(args.eval_dir)

    g_counts_t = loader.g_counts
    c_counts_t = loader.c_counts
    var_names = loader.var_names
    g_meta = loader.g_meta
    c_meta = loader.c_meta

    g_counts = (
        g_counts_t.detach().cpu().numpy().astype(np.float32, copy=False)
        if torch.is_tensor(g_counts_t)
        else np.asarray(g_counts_t, dtype=np.float32)
    )
    c_counts = (
        c_counts_t.detach().cpu().numpy().astype(np.float32, copy=False)
        if torch.is_tensor(c_counts_t)
        else np.asarray(c_counts_t, dtype=np.float32)
    )

    # Consistency checks
    if g_counts.shape[0] != len(g_meta):
        raise RuntimeError(f"CRISPR counts/meta mismatch: counts={g_counts.shape[0]} meta={len(g_meta)}")
    if c_counts.shape[0] != len(c_meta):
        raise RuntimeError(f"Drug counts/meta mismatch: counts={c_counts.shape[0]} meta={len(c_meta)}")
    if g_counts.shape[1] != c_counts.shape[1]:
        raise RuntimeError(f"Gene dim mismatch: CRISPR={g_counts.shape[1]} Drug={c_counts.shape[1]}")
    if g_counts.shape[1] != len(var_names):
        raise RuntimeError(f"var_names mismatch: X has {g_counts.shape[1]} genes but var_names has {len(var_names)}")

    n_g, n_genes = g_counts.shape
    n_c, _ = c_counts.shape

    print(">>> PCA baseline will use the SAME inputs as your FM scripts.")
    print(f">>> eval_dir: {args.eval_dir}")
    print(f">>> out_dir : {args.out_dir}")
    print(f">>> CRISPR  : {n_g} cells | Drug: {n_c} cells | genes: {n_genes}")
    print(f">>> dims to save: {dims} (fit max_comps={max_comps})")
    print(f">>> preprocess: normalize_total({args.target_sum}) + log1p + "
          f"{'scale(max_value='+str(args.scale_max_value)+')' if args.scale_max_value>0 else 'NO scale'}")
    print(f">>> PCA: solver={args.svd_solver}, random_state={args.random_state}")

    # ---- Concatenate and preprocess with Scanpy ----
    X = np.vstack([g_counts, c_counts])  # (n_g+n_c, n_genes)

    obs = {"modality": (["CRISPR"] * n_g) + (["Drug"] * n_c)}
    adata = ad.AnnData(X=X, obs=obs)
    adata.var["gene_symbol"] = [str(g) for g in var_names]

    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    if args.scale_max_value > 0:
        sc.pp.scale(adata, max_value=args.scale_max_value)

    sc.tl.pca(
        adata,
        n_comps=max_comps,
        svd_solver=args.svd_solver,
        random_state=args.random_state,
        use_highly_variable=False,
    )

    Z = np.asarray(adata.obsm["X_pca"], dtype=np.float32)
    if Z.shape != (n_g + n_c, max_comps):
        raise RuntimeError(f"Unexpected PCA shape: got {Z.shape}, expected {(n_g+n_c, max_comps)}")
    if not np.isfinite(Z).all():
        bad = int(np.sum(~np.isfinite(Z)))
        raise RuntimeError(f"PCA embedding contains non-finite values: {bad} entries")

    Zg = Z[:n_g, :]
    Zc = Z[n_g:, :]

    # ---- Save sliced embeddings ----
    for k in dims:
        g_out = os.path.join(args.out_dir, f"PCA{k}_CRISPR_emb.npy")
        c_out = os.path.join(args.out_dir, f"PCA{k}_Drug_emb.npy")
        np.save(g_out, Zg[:, :k].astype(np.float32, copy=False))
        np.save(c_out, Zc[:, :k].astype(np.float32, copy=False))
        print(f"✅ Saved PCA{k}: CRISPR {Zg[:, :k].shape} -> {g_out}")
        print(f"✅ Saved PCA{k}: Drug   {Zc[:, :k].shape} -> {c_out}")

    # ---- Provenance JSON ----
    try:
        import json
        prov = {
            "eval_dir": args.eval_dir,
            "out_dir": args.out_dir,
            "n_crispr": int(n_g),
            "n_drug": int(n_c),
            "n_genes": int(n_genes),
            "dims_saved": dims,
            "preprocess": {
                "normalize_total_target_sum": float(args.target_sum),
                "log1p": True,
                "scale": bool(args.scale_max_value > 0),
                "scale_max_value": float(args.scale_max_value) if args.scale_max_value > 0 else None,
            },
            "pca": {"n_comps_fit": int(max_comps), "svd_solver": args.svd_solver, "random_state": int(args.random_state)},
        }
        prov_path = os.path.join(args.out_dir, "PCA_baseline_provenance.json")
        with open(prov_path, "w", encoding="utf-8") as f:
            json.dump(prov, f, indent=2)
        print(f"🧾 Wrote provenance: {prov_path}")
    except Exception as e:
        print(f"⚠️  Failed to write provenance JSON: {e}")

    print(">>> Done.")


if __name__ == "__main__":
    main()
