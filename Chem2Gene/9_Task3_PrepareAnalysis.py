#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import argparse
import numpy as np
import pandas as pd
import torch


# --------------------------
# IO helpers
# --------------------------
def read_common_targets(path: str) -> set:
    """
    Common_Targets_K562.csv may be an empty dataframe but has Index.
    Treat its index as the common target set.
    """
    df = pd.read_csv(path, index_col=0)
    return set(df.index.astype(str).tolist())


def detect_control_mask(meta: pd.DataFrame) -> np.ndarray:
    """
    Robust control detection across CRISPR/Drug meta.
    Priority:
      - benchmark_group == 'Control'
      - specificity_tier == 'Control'
      - perturbation_raw == 'control' / 'Control'
      - clean_target_mapped == 'control' / 'Control'
      - target contains 'non-targeting'
    """
    m = meta.copy()

    def col_eq(col, val):
        return (m[col].astype(str) == val) if col in m.columns else pd.Series(False, index=m.index)

    def col_lower_eq(col, val_lower):
        return (m[col].astype(str).str.lower() == val_lower) if col in m.columns else pd.Series(False, index=m.index)

    bm = col_eq("benchmark_group", "Control")
    st = col_eq("specificity_tier", "Control")
    pr = col_lower_eq("perturbation_raw", "control")
    ct = col_lower_eq("clean_target_mapped", "control")
    nt = (
        m["target"].astype(str).str.contains("non-targeting", case=False, regex=False)
        if "target" in m.columns else pd.Series(False, index=m.index)
    )

    mask = (bm | st | pr | ct | nt).values
    return mask


def split_target_list(x: str) -> list:
    """
    clean_target_mapped sometimes is like "CDK2;CDK9".
    """
    if x is None:
        return []
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return []
    return [t.strip() for t in s.split(";") if t.strip()]


def ensure_dir(p: str):
    if p is None or p == "":
        return
    os.makedirs(p, exist_ok=True)


def save_npy(path: str, arr: np.ndarray):
    ensure_dir(os.path.dirname(path))
    np.save(path, arr.astype(np.float32))


def save_parquet(path: str, df: pd.DataFrame):
    """
    Save parquet with fallback to CSV if parquet engine missing.
    """
    ensure_dir(os.path.dirname(path))
    df2 = df.reset_index(drop=True)
    try:
        df2.to_parquet(path)
    except Exception as e:
        csv_path = os.path.splitext(path)[0] + ".csv"
        print(f"[WARN] Parquet save failed ({e}). Fallback to CSV: {csv_path}")
        df2.to_csv(csv_path, index=False)


# --------------------------
# Normalization & sampling
# --------------------------
def sample_controls(ctrl_pool: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    if len(ctrl_pool) == 0:
        raise RuntimeError("No control cells detected; cannot sample controls.")
    replace = len(ctrl_pool) < n
    return rng.choice(ctrl_pool, size=n, replace=replace)


def log1p_cpm(x: torch.Tensor, scale: float = 1e4) -> torch.Tensor:
    """
    x: (B, G) non-negative counts (float/int), CPU tensor
    returns: (B, G) float32 log1p(CPM) on CPU
    """
    x = x.to(dtype=torch.float32)
    lib = torch.clamp(x.sum(dim=1, keepdim=True), min=1.0)
    x = x / lib * float(scale)
    return torch.log1p(x)


def streaming_mean_std_log1p_cpm(counts: torch.Tensor,
                                 row_idx: np.ndarray,
                                 scale: float = 1e4,
                                 batch_size: int = 4096) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute mean/std of log1p(CPM) over selected rows (streaming).
    Returns:
      mu: (G,) float32
      sigma: (G,) float32 (unbiased=False)
    """
    idx = row_idx.astype(np.int64)
    G = int(counts.shape[1])
    s1 = torch.zeros((G,), dtype=torch.float64)
    s2 = torch.zeros((G,), dtype=torch.float64)
    n = 0

    for i in range(0, len(idx), batch_size):
        b = idx[i:i + batch_size]
        xb = log1p_cpm(counts[b], scale=scale).to(dtype=torch.float64)  # (B,G)
        s1 += xb.sum(dim=0)
        s2 += (xb * xb).sum(dim=0)
        n += xb.shape[0]

    if n <= 0:
        raise RuntimeError("No rows for streaming mean/std.")

    mu = (s1 / n)
    var = (s2 / n) - (mu * mu)
    var = torch.clamp(var, min=1e-8)
    sigma = torch.sqrt(var)

    return mu.to(dtype=torch.float32).cpu().numpy(), sigma.to(dtype=torch.float32).cpu().numpy()


# --------------------------
# Cell-level delta builders (Task1-aligned)
# --------------------------
def compute_cell_deltas_gene(counts: torch.Tensor,
                             meta: pd.DataFrame,
                             ctrl_mask: np.ndarray,
                             pert_mask: np.ndarray,
                             rng: np.random.Generator,
                             k_controls: int = 1,
                             scale: float = 1e4,
                             use_zscore: bool = False,
                             z_eps: float = 1e-2,
                             z_clip: float = 10.0,
                             batch_size: int = 2048) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Cell-level deltas:
      For each perturbed cell i:
        baseline = mean( log1p_cpm(control_j) for j in sampled controls )
        delta = log1p_cpm(pert_i) - baseline
        optional zscore: delta / (sigma_ctrl + eps)

    Returns:
      meta_out: (N_pert, ...) metadata for each pert cell (with audit fields)
      delta_std: (N_pert, G) float32
    """
    m = meta.copy()
    # Keep original cell IDs as "unique_id"
    if "unique_id" not in m.columns:
        m["unique_id"] = m.index.astype(str)

    m["_row_idx"] = np.arange(len(m), dtype=np.int64)

    ctrl_pool = m.loc[ctrl_mask, "_row_idx"].values.astype(np.int64)
    pert_rows = m.loc[pert_mask].copy()
    if pert_rows.empty:
        raise RuntimeError("No perturbed cells detected under pert_mask.")

    pert_idx = pert_rows["_row_idx"].values.astype(np.int64)
    Np = len(pert_idx)
    G = int(counts.shape[1])

    # control stats for zscore (in normalized space)
    sigma_ctrl = None
    if use_zscore:
        mu_ctrl, sigma_ctrl = streaming_mean_std_log1p_cpm(
            counts=counts,
            row_idx=ctrl_pool,
            scale=scale,
            batch_size=max(1024, batch_size),
        )
        # avoid zeros
        sigma_ctrl = np.maximum(sigma_ctrl, 1e-6).astype(np.float32)

    deltas = []
    meta_rows = []

    k = int(max(1, k_controls))

    for i in range(0, Np, batch_size):
        b_idx = pert_idx[i:i + batch_size]
        B = len(b_idx)

        # sample control indices for each pert cell
        if len(ctrl_pool) >= 1:
            # sample B*k controls
            ctrl_sample = sample_controls(ctrl_pool, B * k, rng).astype(np.int64)
        else:
            raise RuntimeError("Control pool is empty after ctrl_mask filtering.")

        # normalize pert
        x_pert = log1p_cpm(counts[b_idx], scale=scale)  # (B,G)

        # normalize sampled controls
        x_ctrl = log1p_cpm(counts[ctrl_sample], scale=scale)  # (B*k,G)
        x_ctrl = x_ctrl.view(B, k, G).mean(dim=1)  # (B,G) baseline in normalized space

        d = (x_pert - x_ctrl).cpu().numpy().astype(np.float32)

        if use_zscore:
            d = d / (sigma_ctrl[None, :] + float(z_eps))
            if z_clip is not None and z_clip > 0:
                d = np.clip(d, -float(z_clip), float(z_clip)).astype(np.float32)

        deltas.append(d)

        # meta audit
        # map row_idx -> cell unique_id (index)
        ctrl_cell_ids = m.iloc[ctrl_sample]["unique_id"].astype(str).tolist()
        ctrl_cell_ids = np.array(ctrl_cell_ids, dtype=object).reshape(B, k)

        sub_meta = pert_rows.iloc[i:i + B].copy()
        sub_meta["n_ctrl_pool"] = int(len(ctrl_pool))
        sub_meta["n_ctrl_used"] = int(k)
        # store only a short head to keep files small
        sub_meta["ctrl_ids_head"] = ["|".join(ctrl_cell_ids[r, :min(10, k)].tolist()) for r in range(B)]
        sub_meta["ctrl_ids_n"] = int(k)

        # standardize fields used downstream (best-effort)
        if "cell_std" not in sub_meta.columns:
            if "cell_line" in sub_meta.columns:
                sub_meta["cell_std"] = sub_meta["cell_line"].astype(str)
            else:
                sub_meta["cell_std"] = "K562"
        if "target_std" not in sub_meta.columns:
            sub_meta["target_std"] = sub_meta["clean_target_mapped"].astype(str).str.upper()

        meta_rows.append(sub_meta.drop(columns=["_row_idx"], errors="ignore"))

    meta_out = pd.concat(meta_rows, ignore_index=True)
    delta_std = np.concatenate(deltas, axis=0).astype(np.float32)  # (N_pert,G)
    return meta_out, delta_std


def compute_cell_deltas_rep(X: np.ndarray,
                            meta: pd.DataFrame,
                            ctrl_mask: np.ndarray,
                            pert_mask: np.ndarray,
                            rng: np.random.Generator,
                            k_controls: int = 1,
                            batch_size: int = 8192) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Cell-level deltas in embedding space:
      delta_i = emb_pert_i - mean(emb_ctrl_samples)
    Returns:
      meta_out: (N_pert, ...)
      delta_std: (N_pert, D)
    """
    m = meta.copy()
    if "unique_id" not in m.columns:
        m["unique_id"] = m.index.astype(str)
    m["_row_idx"] = np.arange(len(m), dtype=np.int64)

    ctrl_pool = m.loc[ctrl_mask, "_row_idx"].values.astype(np.int64)
    pert_rows = m.loc[pert_mask].copy()
    if pert_rows.empty:
        raise RuntimeError("No perturbed cells detected under pert_mask.")

    pert_idx = pert_rows["_row_idx"].values.astype(np.int64)
    Np = len(pert_idx)
    D = int(X.shape[1])
    k = int(max(1, k_controls))

    deltas = []
    meta_rows = []

    for i in range(0, Np, batch_size):
        b_idx = pert_idx[i:i + batch_size]
        B = len(b_idx)

        ctrl_sample = sample_controls(ctrl_pool, B * k, rng).astype(np.int64)
        ctrl_sample = ctrl_sample.reshape(B, k)

        mu_p = X[b_idx]  # (B,D)
        mu_c = X[ctrl_sample].mean(axis=1)  # (B,D)
        d = (mu_p - mu_c).astype(np.float32)

        deltas.append(d)

        sub_meta = pert_rows.iloc[i:i + B].copy()
        sub_meta["n_ctrl_pool"] = int(len(ctrl_pool))
        sub_meta["n_ctrl_used"] = int(k)

        ctrl_cell_ids = m.iloc[ctrl_sample.reshape(-1)]["unique_id"].astype(str).tolist()
        ctrl_cell_ids = np.array(ctrl_cell_ids, dtype=object).reshape(B, k)
        sub_meta["ctrl_ids_head"] = ["|".join(ctrl_cell_ids[r, :min(10, k)].tolist()) for r in range(B)]
        sub_meta["ctrl_ids_n"] = int(k)

        if "cell_std" not in sub_meta.columns:
            if "cell_line" in sub_meta.columns:
                sub_meta["cell_std"] = sub_meta["cell_line"].astype(str)
            else:
                sub_meta["cell_std"] = "K562"
        if "target_std" not in sub_meta.columns:
            sub_meta["target_std"] = sub_meta["clean_target_mapped"].astype(str).str.upper()

        meta_rows.append(sub_meta.drop(columns=["_row_idx"], errors="ignore"))

    meta_out = pd.concat(meta_rows, ignore_index=True)
    delta_std = np.concatenate(deltas, axis=0).astype(np.float32)
    return meta_out, delta_std


def apply_systema_mean(delta_std: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Simple Systema-like correction (NaN-safe):
      bg = nanmean(delta_std)
      delta_sys = delta_std - bg

    Rationale:
      If any rows/dimensions contain NaN/Inf (e.g., upstream FM embedding issues),
      using mean() will propagate NaNs into bg and then into all rows.
      nanmean() keeps valid rows usable while preserving NaNs in the affected rows.
    """
    # Compute background with NaN-safe mean across cells
    bg = np.nanmean(delta_std, axis=0, keepdims=True).astype(np.float32)

    # If a dimension is entirely NaN/Inf, nanmean may still return NaN.
    # Set those bg entries to 0 so subtraction does not destroy otherwise-valid rows.
    bg[~np.isfinite(bg)] = 0.0

    delta_sys = (delta_std - bg).astype(np.float32)
    return delta_sys, bg.squeeze(0)



# --------------------------
# Embedding loaders
# --------------------------
def load_embeddings(base_path: str,
                    meta_crispr: pd.DataFrame, meta_drug: pd.DataFrame) -> dict:
    """
    Return:
      embs = {
        "Geneformer": {"crispr": np.ndarray, "drug": np.ndarray},
        ...
      }
    All arrays aligned to meta order.
    """
    embs = {}

    gf_c = pd.read_csv(os.path.join(base_path, "Geneformer_Chem2Gen_K562_CRISPR.csv"), index_col=0)
    gf_d = pd.read_csv(os.path.join(base_path, "Geneformer_Chem2Gen_K562_Drug.csv"), index_col=0)
    embs["Geneformer"] = {"crispr": gf_c.values.astype(np.float32), "drug": gf_d.values.astype(np.float32)}

    embs["UCE"] = {
        "crispr": np.load(os.path.join(base_path, "UCE_CRISPR_embeddings.npy")).astype(np.float32),
        "drug": np.load(os.path.join(base_path, "UCE_Drug_embeddings.npy")).astype(np.float32),
    }
    embs["scBERT"] = {
        "crispr": np.load(os.path.join(base_path, "scBERT_CRISPR_emb.npy")).astype(np.float32),
        "drug": np.load(os.path.join(base_path, "scBERT_Drug_emb.npy")).astype(np.float32),
    }
    embs["scFoundation"] = {
        "crispr": np.load(os.path.join(base_path, "scFoundation_CRISPR_embeddings.npy")).astype(np.float32),
        "drug": np.load(os.path.join(base_path, "scFoundation_Drug_embeddings.npy")).astype(np.float32),
    }
    embs["scGPT"] = {
        "crispr": np.load(os.path.join(base_path, "scGPT_CRISPR_emb.npy")).astype(np.float32),
        "drug": np.load(os.path.join(base_path, "scGPT_Drug_emb.npy")).astype(np.float32),
    }
    embs["STATE"] = {
        "crispr": np.load(os.path.join(base_path, "STATE_CRISPR_embeddings.npy")).astype(np.float32),
        "drug": np.load(os.path.join(base_path, "STATE_Drug_embeddings.npy")).astype(np.float32),
    }
    embs["TahoeX1_3b"] = {
        "crispr": np.load(os.path.join(base_path, "TahoeX1_3b_CRISPR_embeddings.npy")).astype(np.float32),
        "drug": np.load(os.path.join(base_path, "TahoeX1_3b_Drug_embeddings.npy")).astype(np.float32),
    }
    
    # --------------------------
    # PCA baselines (optional)
    # --------------------------
    for k in [50, 100, 200]:
        p_c = os.path.join(base_path, f"PCA{k}_CRISPR_emb.npy")
        p_d = os.path.join(base_path, f"PCA{k}_Drug_emb.npy")
        if os.path.exists(p_c) and os.path.exists(p_d):
            embs[f"PCA{k}"] = {
                "crispr": np.load(p_c).astype(np.float32),
                "drug": np.load(p_d).astype(np.float32),
            }
        else:
            print(f"[Info] PCA{k} embeddings not found; skip ({p_c}, {p_d})")


    n_c = len(meta_crispr)
    n_d = len(meta_drug)
    for name, dd in embs.items():
        if dd["crispr"].shape[0] != n_c:
            raise RuntimeError(f"[{name}] CRISPR N mismatch: emb={dd['crispr'].shape[0]} vs meta={n_c}")
        if dd["drug"].shape[0] != n_d:
            raise RuntimeError(f"[{name}] Drug N mismatch: emb={dd['drug'].shape[0]} vs meta={n_d}")

    return embs


# --------------------------
# Pathway matrix (edge list)
# --------------------------
def build_pathway_matrix_from_edges(csv_path: str, gene_symbols: list[str]):
    """
    Build row-normalized pathway mask matrix (P x G) from edge list CSV: columns ['source','target'].
    """
    gene_meta_df = pd.DataFrame({
        "gene_symbol": [str(g).upper() for g in gene_symbols],
        "local_idx": np.arange(len(gene_symbols), dtype=np.int64)
    })
    gene_to_idx = dict(zip(gene_meta_df["gene_symbol"], gene_meta_df["local_idx"]))

    df = pd.read_csv(csv_path)
    if not {"source", "target"}.issubset(df.columns):
        raise ValueError(f"Pathway CSV must contain columns ['source','target'], got {df.columns.tolist()}")

    pathways = sorted(df["source"].astype(str).unique().tolist())
    pathway_to_idx = {p: i for i, p in enumerate(pathways)}

    P = len(pathways)
    G = len(gene_symbols)
    mask = torch.zeros((P, G), dtype=torch.float32)

    for _, row in df.iterrows():
        p = str(row["source"])
        g = str(row["target"]).upper()
        if g in gene_to_idx:
            mask[pathway_to_idx[p], gene_to_idx[g]] = 1.0

    hit_counts = mask.sum(dim=1).cpu().numpy().astype(np.int64)

    row_sums = mask.sum(dim=1, keepdim=True)
    row_sums[row_sums == 0] = 1.0
    mask = mask / row_sums
    return mask, pathways, hit_counts


def project_gene_delta_to_pathway(delta_gene: np.ndarray, mask: torch.Tensor) -> np.ndarray:
    dg = torch.from_numpy(delta_gene.astype(np.float32))
    dp = dg @ mask.T
    return dp.cpu().numpy().astype(np.float32)


# --------------------------
# Main
# --------------------------
def main():
    ap = argparse.ArgumentParser("Task3 (Revised): Prepare cell-level delta (Gene/Pathway/FM) for K562 Chem2Gen evaluation set")
    ap.add_argument("--eval_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562")
    ap.add_argument("--base_path", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results")
    ap.add_argument("--out_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Deltas_K562_CellLevel")
    ap.add_argument("--seed", type=int, default=0)

    # kept for backward compatibility; will be forced to "cell"
    ap.add_argument("--group_level", type=str, default="cell",
                    choices=["cell", "target", "condition"],
                    help="DEPRECATED: this revised Task3 always outputs cell-level deltas; non-cell will be ignored.")

    # control sampling
    ap.add_argument("--k_controls", type=int, default=50,
                    help="Number of controls sampled per pert cell (default=1). Set 50 to mimic Task1's typical k.")
    ap.add_argument("--batch_size", type=int, default=2048)

    # gene delta options
    ap.add_argument("--pb_scale", type=float, default=1e4, help="CPM scale for log1p(CPM) (default 1e4)")
    ap.add_argument("--use_zscore", action="store_true",
                    help="Optional: z-score delta by control std per gene (Task1-like).")
    ap.add_argument("--z_eps", type=float, default=1e-2)
    ap.add_argument("--z_clip", type=float, default=10.0)

    # optional pathway (edge list CSV: source,target)
    ap.add_argument("--pathway_csv", type=str,
                    default="/mnt/NAS_21T/ProjectData/OSMOSIS/resource/OmniPath_annotations/HALLMARK_human.csv",
                    help="Optional pathway edge list CSV with columns [source,target]. If empty, skip pathway delta.")
    ap.add_argument("--pathway_name", type=str, default="HALLMARK")

    # source_db tags (useful if you later merge with Level3 pools)
    ap.add_argument("--source_db_crispr", type=str, default="EvalSet_K562_CRISPR")
    ap.add_argument("--source_db_drug", type=str, default="EvalSet_K562_Drug")

    args = ap.parse_args()
    ensure_dir(args.out_dir)

    if args.group_level != "cell":
        print(f"[WARN] --group_level={args.group_level} ignored. Revised Task3 outputs CELL-LEVEL deltas only.")
        args.group_level = "cell"

    # ---- load eval data ----
    common_targets = read_common_targets(os.path.join(args.eval_dir, "Common_Targets_K562.csv"))
    gene_symbols = pd.read_csv(os.path.join(args.eval_dir, "shared_var_names.csv"))["gene_symbol"].astype(str).tolist()

    pathway_mask = None
    pathway_names = None
    pathway_hit_counts = None
    if args.pathway_csv and os.path.exists(args.pathway_csv):
        pathway_mask, pathway_names, pathway_hit_counts = build_pathway_matrix_from_edges(
            csv_path=args.pathway_csv,
            gene_symbols=gene_symbols,
        )
        print(f"[Info] Pathway loaded: {args.pathway_csv} | P={len(pathway_names)}")
    else:
        print("[Info] No valid pathway CSV; skip pathway delta.")

    crispr_meta = pd.read_csv(os.path.join(args.eval_dir, "CRISPR_meta.csv"), index_col=0)
    drug_meta = pd.read_csv(os.path.join(args.eval_dir, "Drug_meta.csv"), index_col=0)

    crispr_counts = torch.load(os.path.join(args.eval_dir, "CRISPR_counts.pt"), map_location="cpu")
    drug_counts = torch.load(os.path.join(args.eval_dir, "Drug_counts.pt"), map_location="cpu")

    # align checks
    if crispr_counts.shape[0] != len(crispr_meta):
        raise RuntimeError("CRISPR counts rows != CRISPR_meta rows")
    if drug_counts.shape[0] != len(drug_meta):
        raise RuntimeError("Drug counts rows != Drug_meta rows")
    if crispr_counts.shape[1] != len(gene_symbols) or drug_counts.shape[1] != len(gene_symbols):
        raise RuntimeError("Gene columns != shared_var_names length")

    # ---- normalize meta schema ----
    for m in (crispr_meta, drug_meta):
        if "clean_target_mapped" not in m.columns:
            m["clean_target_mapped"] = m["target"].astype(str)
        if "specificity_tier" not in m.columns:
            m["specificity_tier"] = "NA"
        if "perturbation_raw" not in m.columns:
            m["perturbation_raw"] = m["target"].astype(str)
        if "benchmark_group" not in m.columns:
            m["benchmark_group"] = "NA"
        if "perturbation_type" not in m.columns:
            m["perturbation_type"] = "NA"
        if "time" not in m.columns:
            m["time"] = np.nan
        if "dose_value" not in m.columns:
            m["dose_value"] = np.nan
        if "unique_id" not in m.columns:
            m["unique_id"] = m.index.astype(str)

    drug_meta["target_list"] = drug_meta["clean_target_mapped"].apply(split_target_list)
    crispr_meta["target_list"] = crispr_meta["clean_target_mapped"].apply(split_target_list)

    # ---- control / pert masks ----
    c_ctrl = detect_control_mask(crispr_meta)
    d_ctrl = detect_control_mask(drug_meta)
    c_pert = ~c_ctrl
    d_pert = ~d_ctrl

    # ---- keep only perturbations with overlap to common targets (controls retained) ----
    c_keep = np.array([
        any([(g in common_targets) for g in split_target_list(t)])
        for t in crispr_meta["clean_target_mapped"].astype(str).tolist()
    ])
    d_keep = np.array([
        any([(g in common_targets) for g in lst])
        for lst in drug_meta["target_list"].tolist()
    ])

    c_pert = c_pert & c_keep
    d_pert = d_pert & d_keep

    print(f"[Info] CRISPR: N={len(crispr_meta)} | ctrl={int(c_ctrl.sum())} | pert(kept)={int(c_pert.sum())}")
    print(f"[Info] Drug:   N={len(drug_meta)}   | ctrl={int(d_ctrl.sum())} | pert(kept)={int(d_pert.sum())}")
    print(f"[Info] Common targets: {len(common_targets)}")

    # ---- load FM embeddings ----
    embs = load_embeddings(args.base_path, crispr_meta, drug_meta)

    rng = np.random.default_rng(args.seed)

    manifest = {
        "eval_dir": args.eval_dir,
        "out_dir": args.out_dir,
        "seed": args.seed,
        "group_level": "cell",
        "k_controls": int(args.k_controls),
        "tracks": []
    }

    def record(track, view, modality, meta_path, delta_path, shape, extra=None):
        item = {
            "track": track,
            "view": view,
            "modality": modality,
            "meta_path": meta_path,
            "delta_path": delta_path,
            "shape": list(shape),
        }
        if extra:
            item.update(extra)
        manifest["tracks"].append(item)

    # --------------------------
    # Gene / Pathway deltas (cell-level)
    # --------------------------
    for modality in ["CRISPR", "Drug"]:
        if modality == "CRISPR":
            meta = crispr_meta
            counts = crispr_counts
            ctrl_mask = c_ctrl
            pert_mask = c_pert
            source_db = args.source_db_crispr
            modality_std = "Genetic"
        else:
            meta = drug_meta
            counts = drug_counts
            ctrl_mask = d_ctrl
            pert_mask = d_pert
            source_db = args.source_db_drug
            modality_std = "Chemical"

        meta_std, delta_std = compute_cell_deltas_gene(
            counts=counts,
            meta=meta,
            ctrl_mask=ctrl_mask,
            pert_mask=pert_mask,
            rng=rng,
            k_controls=int(args.k_controls),
            scale=float(args.pb_scale),
            use_zscore=bool(args.use_zscore),
            z_eps=float(args.z_eps),
            z_clip=float(args.z_clip),
            batch_size=int(args.batch_size),
        )

        # enrich meta
        meta_std["Modality"] = modality
        meta_std["modality"] = modality_std
        meta_std["source_db"] = source_db
        meta_std["Track"] = "Gene"
        meta_std["View"] = "Standard"
        meta_std["target_list"] = meta_std["clean_target_mapped"].astype(str).apply(
            lambda x: ";".join(split_target_list(x))
        )

        delta_sys, bg = apply_systema_mean(delta_std)
        meta_sys = meta_std.copy()
        meta_sys["View"] = "Systema"

        subdir = os.path.join(args.out_dir, "group_level-cell", "Gene")
        meta_path_std = os.path.join(subdir, f"meta_Gene_Standard_{modality}.parquet")
        meta_path_sys = os.path.join(subdir, f"meta_Gene_Systema_{modality}.parquet")
        delta_path_std = os.path.join(subdir, f"delta_Gene_Standard_{modality}.npy")
        delta_path_sys = os.path.join(subdir, f"delta_Gene_Systema_{modality}.npy")
        bg_path = os.path.join(subdir, f"bg_Gene_mean_{modality}.npy")
        genes_path = os.path.join(subdir, "genes_shared_var_names.txt")

        save_parquet(meta_path_std, meta_std)
        save_parquet(meta_path_sys, meta_sys)
        save_npy(delta_path_std, delta_std)
        save_npy(delta_path_sys, delta_sys)
        save_npy(bg_path, bg[None, :])

        ensure_dir(os.path.dirname(genes_path))
        with open(genes_path, "w", encoding="utf-8") as f:
            for g in gene_symbols:
                f.write(g + "\n")

        record("Gene", "Standard", modality, meta_path_std, delta_path_std, delta_std.shape,
               extra={"feature_names": genes_path, "bg_path": bg_path})
        record("Gene", "Systema", modality, meta_path_sys, delta_path_sys, delta_sys.shape,
               extra={"feature_names": genes_path, "bg_path": bg_path})

        print(f"[Saved] Gene (cell-level) {modality}: std={delta_std.shape} sys={delta_sys.shape}")

        # ---- Pathway delta ----
        if pathway_mask is not None:
            delta_p_std = project_gene_delta_to_pathway(delta_std, pathway_mask)
            delta_p_sys, bgp = apply_systema_mean(delta_p_std)

            meta_p_std = meta_std.copy()
            meta_p_std["Track"] = f"Pathway_{args.pathway_name}"
            meta_p_std["View"] = "Standard"
            meta_p_sys = meta_std.copy()
            meta_p_sys["Track"] = f"Pathway_{args.pathway_name}"
            meta_p_sys["View"] = "Systema"

            subdirp = os.path.join(args.out_dir, "group_level-cell", f"Pathway_{args.pathway_name}")
            meta_path_p_std = os.path.join(subdirp, f"meta_Pathway_{args.pathway_name}_Standard_{modality}.parquet")
            meta_path_p_sys = os.path.join(subdirp, f"meta_Pathway_{args.pathway_name}_Systema_{modality}.parquet")
            delta_path_p_std = os.path.join(subdirp, f"delta_Pathway_{args.pathway_name}_Standard_{modality}.npy")
            delta_path_p_sys = os.path.join(subdirp, f"delta_Pathway_{args.pathway_name}_Systema_{modality}.npy")
            bgp_path = os.path.join(subdirp, f"bg_Pathway_{args.pathway_name}_mean_{modality}.npy")
            pathways_path = os.path.join(subdirp, f"pathways_{args.pathway_name}.txt")
            hitcount_path = os.path.join(subdirp, f"pathways_{args.pathway_name}_hit_counts.csv")

            save_parquet(meta_path_p_std, meta_p_std)
            save_parquet(meta_path_p_sys, meta_p_sys)
            save_npy(delta_path_p_std, delta_p_std)
            save_npy(delta_path_p_sys, delta_p_sys)
            save_npy(bgp_path, bgp[None, :])

            ensure_dir(os.path.dirname(pathways_path))
            with open(pathways_path, "w", encoding="utf-8") as f:
                for p in pathway_names:
                    f.write(p + "\n")
            pd.DataFrame({"pathway": pathway_names, "n_hit_genes": pathway_hit_counts}).to_csv(hitcount_path, index=False)

            record(f"Pathway_{args.pathway_name}", "Standard", modality, meta_path_p_std, delta_path_p_std, delta_p_std.shape,
                   extra={"feature_names": pathways_path, "bg_path": bgp_path, "pathway_csv": args.pathway_csv, "hit_counts": hitcount_path})
            record(f"Pathway_{args.pathway_name}", "Systema", modality, meta_path_p_sys, delta_path_p_sys, delta_p_sys.shape,
                   extra={"feature_names": pathways_path, "bg_path": bgp_path, "pathway_csv": args.pathway_csv, "hit_counts": hitcount_path})

            print(f"[Saved] Pathway({args.pathway_name}) (cell-level) {modality}: std={delta_p_std.shape} sys={delta_p_sys.shape}")

    # --------------------------
    # FM embedding deltas (cell-level)
    # --------------------------
    for model_name, dd in embs.items():
        for modality in ["CRISPR", "Drug"]:
            if modality == "CRISPR":
                meta = crispr_meta
                X = dd["crispr"]
                ctrl_mask = c_ctrl
                pert_mask = c_pert
                source_db = args.source_db_crispr
                modality_std = "Genetic"
            else:
                meta = drug_meta
                X = dd["drug"]
                ctrl_mask = d_ctrl
                pert_mask = d_pert
                source_db = args.source_db_drug
                modality_std = "Chemical"

            meta_std, delta_std = compute_cell_deltas_rep(
                X=X,
                meta=meta,
                ctrl_mask=ctrl_mask,
                pert_mask=pert_mask,
                rng=rng,
                k_controls=int(args.k_controls),
                batch_size=max(4096, int(args.batch_size)),
            )

            meta_std["Modality"] = modality
            meta_std["modality"] = modality_std
            meta_std["source_db"] = source_db
            meta_std["Track"] = model_name
            meta_std["View"] = "Standard"
            meta_std["target_list"] = meta_std["clean_target_mapped"].astype(str).apply(
                lambda x: ";".join(split_target_list(x))
            )

            delta_sys, bg = apply_systema_mean(delta_std)
            meta_sys = meta_std.copy()
            meta_sys["View"] = "Systema"

            subdir = os.path.join(args.out_dir, "group_level-cell", "FM", model_name)
            meta_path_std = os.path.join(subdir, f"meta_{model_name}_Standard_{modality}.parquet")
            meta_path_sys = os.path.join(subdir, f"meta_{model_name}_Systema_{modality}.parquet")
            delta_path_std = os.path.join(subdir, f"delta_{model_name}_Standard_{modality}.npy")
            delta_path_sys = os.path.join(subdir, f"delta_{model_name}_Systema_{modality}.npy")
            bg_path = os.path.join(subdir, f"bg_{model_name}_mean_{modality}.npy")

            save_parquet(meta_path_std, meta_std)
            save_parquet(meta_path_sys, meta_sys)
            save_npy(delta_path_std, delta_std)
            save_npy(delta_path_sys, delta_sys)
            save_npy(bg_path, bg[None, :])

            record(model_name, "Standard", modality, meta_path_std, delta_path_std, delta_std.shape,
                   extra={"bg_path": bg_path})
            record(model_name, "Systema", modality, meta_path_sys, delta_path_sys, delta_sys.shape,
                   extra={"bg_path": bg_path})

            print(f"[Saved] FM={model_name} (cell-level) {modality}: std={delta_std.shape} sys={delta_sys.shape}")

    # ---- save manifest ----
    manifest_path = os.path.join(args.out_dir, "manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    print(f"\n[Done] All CELL-LEVEL deltas prepared. Manifest: {manifest_path}")
    print("You can do target-level centroids / pseudobulk later during retrieval/statistics, without losing per-cell distribution.")


if __name__ == "__main__":
    main()
