#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Chem2Gen → UCE cell embedding extractor

Inputs (Chem2Gen Eval Set):
  - CRISPR_counts.pt, CRISPR_meta.csv
  - Drug_counts.pt, Drug_meta.csv
  - shared_var_names.csv (gene_symbol)

Outputs:
  - {out_dir}/UCE_CRISPR_embeddings.npy
  - {out_dir}/UCE_Drug_embeddings.npy
  - (optional) UCE embedded .h5ad files in {work_dir}

How it works:
  1) Build two AnnData (.h5ad) from Chem2Gen tensors + meta
  2) Run UCE embedding script/CLI on each .h5ad
  3) Extract adata.obsm["X_uce"] and save as .npy

UCE interface reference:
  - CLI: python eval_single_anndata.py --adata_path ... --dir ... --species ... --model_loc ... --batch_size ...
  - Output: dir/{dataset_name}.h5ad with embeddings in .obsm["X_uce"]
"""

import os
import sys
import json
import shutil
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from scipy import sparse

try:
    import scanpy as sc
except ImportError:
    sc = None
    import anndata as ad


def _read_var_names(eval_dir: str) -> list[str]:
    var_path = os.path.join(eval_dir, "shared_var_names.csv")
    df = pd.read_csv(var_path)
    if "gene_symbol" not in df.columns:
        raise ValueError(f"shared_var_names.csv must contain 'gene_symbol' column, got: {df.columns.tolist()}")
    genes = df["gene_symbol"].astype(str).tolist()
    if len(genes) == 0:
        raise ValueError("shared_var_names.csv has 0 genes.")
    return genes


def _load_counts_and_meta(eval_dir: str, source: str):
    """
    source: "CRISPR" or "Drug"
    returns: counts (torch.Tensor on CPU), meta (pd.DataFrame)
    """
    if source == "CRISPR":
        counts_path = os.path.join(eval_dir, "CRISPR_counts.pt")
        meta_path = os.path.join(eval_dir, "CRISPR_meta.csv")
    elif source == "Drug":
        counts_path = os.path.join(eval_dir, "Drug_counts.pt")
        meta_path = os.path.join(eval_dir, "Drug_meta.csv")
    else:
        raise ValueError("source must be 'CRISPR' or 'Drug'")

    counts = torch.load(counts_path, map_location="cpu")
    if not isinstance(counts, torch.Tensor):
        raise TypeError(f"{counts_path} did not load as torch.Tensor; got {type(counts)}")
    if counts.dim() != 2:
        raise ValueError(f"{counts_path} must be 2D (cells, genes); got shape {tuple(counts.shape)}")

    meta = pd.read_csv(meta_path, index_col=0)
    if len(meta) != counts.shape[0]:
        raise ValueError(
            f"Row mismatch for {source}: meta has {len(meta)} rows but counts has {counts.shape[0]} rows."
        )

    return counts, meta


def _build_anndata(counts: torch.Tensor,
                   genes: list[str],
                   meta: pd.DataFrame,
                   make_sparse: bool = True):
    """
    Build AnnData where:
      - X: counts (cells x genes)
      - obs: meta (preserve row order)
      - var_names: genes
    """
    if counts.shape[1] != len(genes):
        raise ValueError(f"Gene mismatch: counts has {counts.shape[1]} cols but shared_var_names has {len(genes)} genes.")

    # Torch -> numpy (shares memory if possible)
    x = counts.numpy()

    # Ensure float32 to keep files smaller; UCE expects "counts" in .X (numeric)
    if x.dtype != np.float32:
        x = x.astype(np.float32, copy=False)

    if make_sparse:
        # If data is very dense, CSR may not save memory; but for scRNA counts it's usually beneficial.
        X = sparse.csr_matrix(x, dtype=np.float32)
    else:
        X = x

    obs = meta.copy()
    # Important: keep original index as obs_names
    # anndata/scanpy will preserve obs.index as obs_names
    var = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))

    if sc is not None:
        adata = sc.AnnData(X=X, obs=obs, var=var)
    else:
        adata = ad.AnnData(X=X, obs=obs, var=var)

    return adata


def _write_h5ad(adata, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # Avoid heavy compression by default; you can toggle if needed.
    adata.write_h5ad(path)


def _run_cmd(cmd: list[str], cwd: str | None = None):
    print(f"\n[UCE] Running command:\n  {' '.join(cmd)}")
    p = subprocess.run(cmd, cwd=cwd)
    if p.returncode != 0:
        raise RuntimeError(f"UCE command failed with return code {p.returncode}.")


def _invoke_uce(adata_path: str,
                work_dir: str,
                species: str,
                batch_size: int,
                model_loc: str,
                nlayers: int,
                uce_repo_dir: str):
    """
    Prefer:
      1) running repo script: python eval_single_anndata.py ... (if uce_repo_dir provided)
      2) running console script: uce-eval-single-anndata ... (if installed from PyPI)
    """
    os.makedirs(work_dir, exist_ok=True)

    base_args = [
        "--adata_path", adata_path,
        "--dir", work_dir,
        "--species", species,
        "--batch_size", str(batch_size),
    ]
    if model_loc:
        base_args += ["--model_loc", model_loc]
    if nlayers is not None:
        base_args += ["--nlayers", str(nlayers)]

    # (1) Repo mode
    if uce_repo_dir:
        script_path = os.path.join(uce_repo_dir, "eval_single_anndata.py")
        if not os.path.exists(script_path):
            raise FileNotFoundError(f"--uce_repo_dir provided but eval_single_anndata.py not found: {script_path}")
        cmd = [sys.executable, script_path] + base_args
        _run_cmd(cmd, cwd=uce_repo_dir)
        return

    # (2) PyPI console script
    exe = shutil.which("uce-eval-single-anndata")
    if exe is not None:
        cmd = [exe] + base_args
        _run_cmd(cmd, cwd=None)
        return

    # Fallback hints
    raise RuntimeError(
        "Could not find UCE runner.\n"
        "Please either:\n"
        "  A) clone UCE repo and pass --uce_repo_dir /path/to/UCE\n"
        "  B) pip install uce-model (then ensure 'uce-eval-single-anndata' is on PATH)\n"
    )


def _read_h5ad(path: str):
    if sc is not None:
        return sc.read_h5ad(path)
    import anndata as ad
    return ad.read_h5ad(path)


def _extract_embeddings_preserve_order(embedded_h5ad: str,
                                      original_obs_names: list[str],
                                      obsm_key: str = "X_uce") -> np.ndarray:
    adata = _read_h5ad(embedded_h5ad)
    if obsm_key not in adata.obsm_keys():
        raise KeyError(f"Expected embeddings in adata.obsm['{obsm_key}'], but keys are: {list(adata.obsm_keys())}")

    emb = adata.obsm[obsm_key]
    # Ensure 2D numpy
    emb = np.asarray(emb)
    if emb.ndim != 2:
        raise ValueError(f"Embedding must be 2D; got shape {emb.shape}")

    # Reorder/align to original cell order if needed
    out_names = adata.obs_names.astype(str)
    orig_names = pd.Index([str(x) for x in original_obs_names])

    if len(out_names) == len(orig_names) and out_names.equals(orig_names):
        return emb

    # Build mapping via index
    emb_df = pd.DataFrame(emb, index=pd.Index(out_names, name="cell_id"))
    aligned = emb_df.reindex(orig_names)

    n_missing = int(aligned.isna().any(axis=1).sum())
    if n_missing > 0:
        print(f"[WARN] {n_missing} cells are missing embeddings after UCE processing. "
              f"These rows will be NaN in the saved matrix. "
              f"Consider checking UCE filtering behavior or input validity.")
    return aligned.to_numpy()


def main():
    parser = argparse.ArgumentParser("Chem2Gen → UCE embedding extraction")

    parser.add_argument("--eval_dir", type=str,
                        default="/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562",
                        help="Chem2Gen evaluation dataset directory.")
    parser.add_argument("--out_dir", type=str,
                        default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results",
                        help="Where to save final .npy embeddings.")
    parser.add_argument("--work_dir", type=str,
                        default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/UCE_workdir_K562/",
                        help="UCE working directory (intermediate + embedded .h5ad).")

    parser.add_argument("--species", type=str, default="human",
                        help="UCE species string (e.g., human, mouse, etc.).")
    parser.add_argument("--batch_size", type=int, default=24,
                        help="Per-GPU batch size for UCE.")
    parser.add_argument("--model_loc", type=str, default="/home/lenislin/Experiment/projects/Chem2Gen/benchmark/UCE/model_files/4layer_model.torch",
                        help="Path to UCE weights .torch/.pt if you want to pin it; otherwise UCE may download.")
    parser.add_argument("--nlayers", type=int, default=None,
                        help="Set to 33 if using large model; otherwise default is typically 4-layer.")
    parser.add_argument("--uce_repo_dir", type=str, default="/home/lenislin/Experiment/projects/Chem2Gen/benchmark/UCE",
                        help="If provided, run `python eval_single_anndata.py` from this repo path.")

    parser.add_argument("--make_sparse", action="store_true",
                        help="Store AnnData.X as CSR sparse matrix (recommended for count matrices).")

    args = parser.parse_args()

    eval_dir = args.eval_dir
    out_dir = args.out_dir
    work_dir = args.work_dir

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    genes = _read_var_names(eval_dir)
    print(f"[INFO] Loaded {len(genes)} genes from shared_var_names.csv")

    for source in ["CRISPR", "Drug"]:
        print(f"\n==============================")
        print(f"[STEP] Processing {source}")
        print(f"==============================")

        counts, meta = _load_counts_and_meta(eval_dir, source)
        print(f"[INFO] {source}: counts shape = {tuple(counts.shape)} | meta rows = {len(meta)}")

        # Build input h5ad under work_dir/inputs to avoid name collision with UCE output
        input_dir = os.path.join(work_dir, "inputs")
        os.makedirs(input_dir, exist_ok=True)
        dataset_name = f"Chem2Gen_K562_{source}"
        input_h5ad = os.path.join(input_dir, f"{dataset_name}.h5ad")

        if not os.path.exists(input_h5ad):
            print(f"[STEP] Writing input AnnData: {input_h5ad}")
            adata = _build_anndata(counts=counts, genes=genes, meta=meta, make_sparse=args.make_sparse)
            _write_h5ad(adata, input_h5ad)
        else:
            print(f"[INFO] Input AnnData exists, reuse: {input_h5ad}")

        # Run UCE → output expected at work_dir/{dataset_name}.h5ad (UCE standard behavior)
        _invoke_uce(
            adata_path=input_h5ad,
            work_dir=args.work_dir,
            species=args.species,
            batch_size=args.batch_size,
            model_loc=args.model_loc,
            nlayers=args.nlayers,
            uce_repo_dir=args.uce_repo_dir
        )

        embedded_h5ad = os.path.join(work_dir, f"{dataset_name}_uce_adata.h5ad")
        if not os.path.exists(embedded_h5ad):
            # Fallback: search for the newest h5ad in work_dir
            h5ads = sorted(Path(work_dir).glob("*.h5ad"), key=lambda p: p.stat().st_mtime, reverse=True)
            if len(h5ads) == 0:
                raise FileNotFoundError(f"No embedded .h5ad found in {work_dir} after running UCE.")
            embedded_h5ad = str(h5ads[0])
            print(f"[WARN] Expected {dataset_name}.h5ad not found. "
                  f"Using newest .h5ad instead: {embedded_h5ad}")
        else:
            print(f"[INFO] Found embedded AnnData: {embedded_h5ad}")

        # Extract embeddings
        original_obs_names = meta.index.astype(str).tolist()
        emb = _extract_embeddings_preserve_order(
            embedded_h5ad=embedded_h5ad,
            original_obs_names=original_obs_names,
            obsm_key="X_uce"
        )
        print(f"[INFO] Extracted embeddings: shape = {emb.shape}")

        save_path = os.path.join(out_dir, f"UCE_{source}_embeddings.npy")
        np.save(save_path, emb)
        print(f"[DONE] Saved: {save_path}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
