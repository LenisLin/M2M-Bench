#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chem2Gen-Bench — Task3 Analysis (NaN-robust)

This script consumes Task3 "cell-level deltas" produced by 9_Task3_PrepareAnalysis.py and runs:
  A) Retrieval evaluation
  B) Target-level pairwise alignment (centroid cosine + Energy Distance)

Key improvement:
  - Robust handling of NaN/Inf rows in embeddings/deltas so that a small amount of NaN
    does not produce garbage metrics (e.g., NaN -> rank==1 artifact).
  - Writes QC summary: Task3_QC_NaNFiltering.csv

Inputs (expected in --task3_out_dir):
  - manifest.json
  - meta_*.parquet (or .csv fallback) + delta_*.npy per track/view/modality

Outputs (written to --out_dir):
  - Task3_Retrieval_PerQuery.csv
  - Task3_Retrieval_Summary.csv
  - Task3_Pairwise_PerTarget.csv
  - Task3_Pairwise_Summary.csv
  - Task3_QC_NaNFiltering.csv
"""
from __future__ import annotations

import os
import json
import argparse
from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional

import numpy as np
import pandas as pd

try:
    from scipy.spatial.distance import pdist, cdist
except Exception:
    pdist = None
    cdist = None


# --------------------------
# I/O helpers
# --------------------------
def read_json(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def ensure_dir(p: str):
    if p:
        os.makedirs(p, exist_ok=True)

def read_meta(path: str) -> pd.DataFrame:
    """
    Task3 meta is saved as parquet with CSV fallback in the Task3 script.
    Here we support both.
    """
    if os.path.exists(path):
        if path.endswith(".parquet"):
            try:
                return pd.read_parquet(path)
            except Exception:
                csv_path = os.path.splitext(path)[0] + ".csv"
                if os.path.exists(csv_path):
                    return pd.read_csv(csv_path)
                raise
        else:
            return pd.read_csv(path)

    csv_path = os.path.splitext(path)[0] + ".csv"
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    raise FileNotFoundError(f"Meta file not found: {path} (or CSV fallback)")

def load_npy(path: str) -> np.ndarray:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return np.load(path, mmap_mode=None)

def safe_upper(x) -> str:
    if x is None:
        return ""
    s = str(x).strip()
    return s.upper()

def get_target_series(meta: pd.DataFrame) -> pd.Series:
    """
    Robustly get target labels from meta.
    Priority:
      target_std > clean_target_mapped > target
    """
    if "target_std" in meta.columns:
        s = meta["target_std"]
    elif "clean_target_mapped" in meta.columns:
        s = meta["clean_target_mapped"]
    elif "target" in meta.columns:
        s = meta["target"]
    else:
        raise KeyError("Meta must contain one of: target_std / clean_target_mapped / target")
    return s

def get_query_ids(meta: pd.DataFrame) -> np.ndarray:
    """
    Stable Query_ID:
      unique_id > obs_names > index > fallback range
    """
    if "unique_id" in meta.columns:
        return meta["unique_id"].astype(str).to_numpy()
    if "obs_names" in meta.columns:
        return meta["obs_names"].astype(str).to_numpy()
    if meta.index is not None:
        try:
            return meta.index.astype(str).to_numpy()
        except Exception:
            pass
    return np.arange(len(meta)).astype(str)


# --------------------------
# NaN / Inf handling
# --------------------------
def finite_row_mask(X: np.ndarray) -> np.ndarray:
    return np.isfinite(X).all(axis=1)

def apply_nan_policy(
    meta: pd.DataFrame,
    X: np.ndarray,
    *,
    policy: str,
    target_std: Optional[pd.Series] = None,
    qc_rows: Optional[List[dict]] = None,
    qc_key: Optional[dict] = None,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    policy:
      - drop: drop rows with any non-finite values
      - zero: replace non-finite with 0.0 (keeps rows)
    Also drops empty target labels if target_std is provided.
    """
    if qc_rows is None:
        qc_rows = []
    if qc_key is None:
        qc_key = {}

    n0 = int(X.shape[0])
    nonfinite = (~finite_row_mask(X)).sum()

    # drop empty targets if provided
    empty_target = 0
    keep_label_mask = None
    if target_std is not None:
        t = target_std.astype(str).fillna("").str.strip()
        keep_label_mask = (t != "")
        empty_target = int((~keep_label_mask).sum())

    meta2 = meta.copy()
    X2 = X

    if policy.lower() == "drop":
        keep_mask = finite_row_mask(X2)
        if keep_label_mask is not None:
            keep_mask = keep_mask & keep_label_mask.to_numpy()
        meta2 = meta2.loc[keep_mask].reset_index(drop=True)
        X2 = X2[keep_mask]
    elif policy.lower() == "zero":
        # replace non-finite to 0, but still drop empty target if provided
        X2 = np.where(np.isfinite(X2), X2, 0.0).astype(np.float32, copy=False)
        if keep_label_mask is not None:
            meta2 = meta2.loc[keep_label_mask].reset_index(drop=True)
            X2 = X2[keep_label_mask.to_numpy()]
    else:
        raise ValueError(f"Unknown nan_policy: {policy}. Use 'drop' or 'zero'.")

    n1 = int(X2.shape[0])

    qc_rows.append({
        **qc_key,
        "n_total": n0,
        "n_nonfinite_rows": int(nonfinite),
        "n_empty_target": int(empty_target),
        "nan_policy": policy,
        "n_kept": n1,
        "pct_kept": (n1 / n0 * 100.0) if n0 > 0 else np.nan,
    })
    return meta2, X2


# --------------------------
# Math helpers (Task1-style)
# --------------------------
def robust_cosine_matrix(A: np.ndarray, B: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """
    Cosine similarity between row-vectors:
      sim(i,j) = <Ai,Bj> / (||Ai|| * ||Bj||)
    """
    A = A.astype(np.float32, copy=False)
    B = B.astype(np.float32, copy=False)

    # sanitize any remaining non-finite (should not happen under drop; helps under zero/edge cases)
    A = np.where(np.isfinite(A), A, 0.0)
    B = np.where(np.isfinite(B), B, 0.0)

    An = np.linalg.norm(A, axis=1, keepdims=True)
    Bn = np.linalg.norm(B, axis=1, keepdims=True)
    An = np.maximum(An, eps)
    Bn = np.maximum(Bn, eps)

    S = (A @ B.T) / (An @ Bn.T)
    # any numeric accidents -> set to very small similarity
    S = np.where(np.isfinite(S), S, -np.inf).astype(np.float32, copy=False)
    return S

def tie_aware_rank(sim_row: np.ndarray, true_idx: int) -> int:
    """
    1-based rank, tie-aware:
      if ties at true score, return best possible rank within tie block.
    Robustified against NaN/Inf: non-finite treated as -inf; true score non-finite -> worst rank.
    """
    scores = np.asarray(sim_row, dtype=np.float32)
    scores = np.where(np.isfinite(scores), scores, -np.inf)

    ts = float(scores[true_idx])
    if not np.isfinite(ts):
        return int(len(scores))  # worst
    higher = int(np.sum(scores > ts))
    return higher + 1

def norm_success(rank_1based: int, n_gallery: int) -> float:
    """
    Normalized success score in [0,1]:
      1.0 if rank==1, 0.0 if rank==n_gallery.
    """
    if n_gallery <= 1:
        return 1.0
    r = int(rank_1based)
    r = max(1, min(r, n_gallery))
    return float((n_gallery - r) / (n_gallery - 1))


# --------------------------
# Aggregation helpers
# --------------------------
def centroid_by_group(X: np.ndarray, labels: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute group means of X by label (finite-only assumed).
    Returns:
      centroids: (K,D)
      uniq_labels: (K,)
      counts: (K,)
    """
    labels = labels.astype(object)
    uniq = np.array(sorted(pd.unique(labels).tolist()), dtype=object)
    K = len(uniq)
    D = int(X.shape[1])
    cent = np.zeros((K, D), dtype=np.float32)
    cnt = np.zeros((K,), dtype=np.int64)
    for i, u in enumerate(uniq):
        idx = np.where(labels == u)[0]
        cnt[i] = len(idx)
        if cnt[i] > 0:
            cent[i] = X[idx].mean(axis=0).astype(np.float32)
    return cent, uniq, cnt


# --------------------------
# Energy distance (distributional alignment)
# --------------------------
def energy_distance_sqeuclidean(X: np.ndarray, Y: np.ndarray) -> float:
    """
    E-distance with squared Euclidean distances:
      E = 2*δXY - σX - σY
    Requires >=2 samples on each side (after filtering).
    """
    # filter any remaining non-finite rows (defensive)
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    X = X[np.isfinite(X).all(axis=1)]
    Y = Y[np.isfinite(Y).all(axis=1)]
    n, m = X.shape[0], Y.shape[0]
    if n < 2 or m < 2:
        return float("nan")

    if cdist is not None and pdist is not None:
        dxy = cdist(X, Y, metric="sqeuclidean")
        dxx = pdist(X, metric="sqeuclidean")
        dyy = pdist(Y, metric="sqeuclidean")
        delta_xy = float(dxy.mean())
        sigma_x = float(dxx.mean()) * 2.0
        sigma_y = float(dyy.mean()) * 2.0
        return 2.0 * delta_xy - sigma_x - sigma_y

    dxy = ((X[:, None, :] - Y[None, :, :]) ** 2).sum(axis=2)
    delta_xy = float(dxy.mean())
    dxx = ((X[:, None, :] - X[None, :, :]) ** 2).sum(axis=2)
    dyy = ((Y[:, None, :] - Y[None, :, :]) ** 2).sum(axis=2)
    sigma_x = float(dxx[~np.eye(n, dtype=bool)].mean())
    sigma_y = float(dyy[~np.eye(m, dtype=bool)].mean())
    return 2.0 * delta_xy - sigma_x - sigma_y

def energy_distance_bootstrap(
    X: np.ndarray,
    Y: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 50,
    max_n: int = 200
) -> Tuple[float, float]:
    """
    Bootstrap mean+sd of energy distance by subsampling up to max_n cells per side.
    Returns (nan,nan) if insufficient samples.
    """
    # filter defensive
    X = X[np.isfinite(X).all(axis=1)]
    Y = Y[np.isfinite(Y).all(axis=1)]
    if X.shape[0] < 2 or Y.shape[0] < 2:
        return float("nan"), float("nan")

    vals = []
    for _ in range(int(n_boot)):
        ix = rng.choice(X.shape[0], size=min(max_n, X.shape[0]), replace=X.shape[0] < max_n)
        iy = rng.choice(Y.shape[0], size=min(max_n, Y.shape[0]), replace=Y.shape[0] < max_n)
        ed = energy_distance_sqeuclidean(X[ix], Y[iy])
        vals.append(ed)
    vals = np.array(vals, dtype=np.float64)
    return float(np.nanmean(vals)), float(np.nanstd(vals))


# --------------------------
# Manifest parsing
# --------------------------
@dataclass
class TrackEntry:
    track: str
    view: str
    modality: str  # "CRISPR" or "Drug"
    meta_path: str
    delta_path: str
    shape: Tuple[int, int]
    extra: dict

def load_manifest(task3_out_dir: str) -> Tuple[dict, List[TrackEntry]]:
    man_path = os.path.join(task3_out_dir, "manifest.json")
    man = read_json(man_path)
    entries = []
    for item in man.get("tracks", []):
        e = TrackEntry(
            track=item["track"],
            view=item["view"],
            modality=item["modality"],
            meta_path=item["meta_path"],
            delta_path=item["delta_path"],
            shape=tuple(item["shape"]),
            extra={k: v for k, v in item.items() if k not in {"track", "view", "modality", "meta_path", "delta_path", "shape"}},
        )
        entries.append(e)
    return man, entries

def index_entries(entries: List[TrackEntry]) -> Dict[Tuple[str, str, str], TrackEntry]:
    out = {}
    for e in entries:
        out[(e.track, e.view, e.modality)] = e
    return out


# --------------------------
# Analysis A: Retrieval
# --------------------------
def run_retrieval(
    index: Dict[Tuple[str, str, str], TrackEntry],
    tracks: List[str],
    views: List[str],
    directions: List[Tuple[str, str]],
    out_dir: str,
    nan_policy: str,
    qc_rows: List[dict],
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    per_query_rows = []

    for track in tracks:
        for view in views:
            for q_mod, g_mod in directions:
                key_q = (track, view, q_mod)
                key_g = (track, view, g_mod)
                if key_q not in index or key_g not in index:
                    continue

                eq = index[key_q]
                eg = index[key_g]

                mq0 = read_meta(eq.meta_path)
                mg0 = read_meta(eg.meta_path)

                Xq0 = load_npy(eq.delta_path).astype(np.float32)
                Xg0 = load_npy(eg.delta_path).astype(np.float32)

                if Xq0.shape[0] != len(mq0) or Xg0.shape[0] != len(mg0):
                    raise RuntimeError(f"Meta/delta mismatch for {track}/{view} {q_mod} or {g_mod}")

                # targets
                tq = get_target_series(mq0).map(safe_upper)
                tg = get_target_series(mg0).map(safe_upper)
                mq0 = mq0.copy()
                mg0 = mg0.copy()
                mq0["target_std"] = tq
                mg0["target_std"] = tg

                # apply NaN policy + drop empty targets
                mq, Xq = apply_nan_policy(
                    mq0, Xq0,
                    policy=nan_policy,
                    target_std=mq0["target_std"],
                    qc_rows=qc_rows,
                    qc_key={"phase": "retrieval", "track": track, "view": view, "modality": q_mod, "role": "query"},
                )
                mg, Xg = apply_nan_policy(
                    mg0, Xg0,
                    policy=nan_policy,
                    target_std=mg0["target_std"],
                    qc_rows=qc_rows,
                    qc_key={"phase": "retrieval", "track": track, "view": view, "modality": g_mod, "role": "gallery"},
                )

                if len(mq) == 0 or len(mg) == 0:
                    qc_rows.append({
                        "phase": "retrieval",
                        "track": track, "view": view,
                        "direction": f"{q_mod}->{g_mod}",
                        "note": "empty_after_filtering",
                        "n_queries_kept": int(len(mq)),
                        "n_gallery_kept": int(len(mg)),
                    })
                    continue

                # Build gallery centroids by target
                G_cent, G_labels, G_counts = centroid_by_group(Xg, mg["target_std"].to_numpy())
                n_gallery = int(len(G_labels))
                label_to_idx = {lab: i for i, lab in enumerate(G_labels.tolist())}

                # Cosine similarity
                S = robust_cosine_matrix(Xq, G_cent)

                # stable query id
                q_ids = get_query_ids(mq)
                true_labels = mq["target_std"].to_numpy()

                n_skipped_missing = 0
                for i in range(S.shape[0]):
                    tl = true_labels[i]
                    if tl not in label_to_idx:
                        n_skipped_missing += 1
                        continue
                    j = label_to_idx[tl]
                    r = tie_aware_rank(S[i], j)
                    per_query_rows.append({
                        "Track": track,
                        "View": view,
                        "Direction": f"{q_mod}->{g_mod}",
                        "Query_Modality": q_mod,
                        "Gallery_Modality": g_mod,
                        "Query_ID": q_ids[i],
                        "True_Target": tl,
                        "True_Rank": int(r),
                        "MRR": float(1.0 / r),
                        "Success_Score": norm_success(r, n_gallery),
                        "N_Gallery": n_gallery,
                        "N_Cells_Target_in_Gallery": int(G_counts[j]),
                    })

                qc_rows.append({
                    "phase": "retrieval",
                    "track": track, "view": view,
                    "direction": f"{q_mod}->{g_mod}",
                    "note": "after_eval",
                    "n_queries_kept": int(len(mq)),
                    "n_gallery_targets": int(n_gallery),
                    "n_skipped_queries_true_target_not_in_gallery": int(n_skipped_missing),
                })

    per_query = pd.DataFrame(per_query_rows)
    if per_query.empty:
        raise RuntimeError("No retrieval results produced (after NaN filtering). Check manifest & filters.")

    summary = (
        per_query
        .groupby(["Track", "View", "Direction"], as_index=False)
        .agg(
            n_queries=("Query_ID", "count"),
            mean_mrr=("MRR", "mean"),
            mean_success=("Success_Score", "mean"),
            median_rank=("True_Rank", "median"),
        )
        .sort_values(["Track", "View", "Direction"])
    )

    ensure_dir(out_dir)
    per_query_path = os.path.join(out_dir, "Task3_Retrieval_PerQuery.csv")
    summary_path = os.path.join(out_dir, "Task3_Retrieval_Summary.csv")
    per_query.to_csv(per_query_path, index=False)
    summary.to_csv(summary_path, index=False)

    return per_query, summary


# --------------------------
# Analysis B: Pairwise per-target alignment
# --------------------------
def run_pairwise(
    index: Dict[Tuple[str, str, str], TrackEntry],
    tracks: List[str],
    views: List[str],
    out_dir: str,
    nan_policy: str,
    qc_rows: List[dict],
    n_boot: int = 30,
    max_n: int = 200,
    seed: int = 0
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    rng = np.random.default_rng(seed)
    rows = []

    for track in tracks:
        for view in views:
            key_c = (track, view, "Drug")
            key_g = (track, view, "CRISPR")
            if key_c not in index or key_g not in index:
                continue

            ec = index[key_c]
            eg = index[key_g]

            mc0 = read_meta(ec.meta_path)
            mg0 = read_meta(eg.meta_path)

            Xc0 = load_npy(ec.delta_path).astype(np.float32)
            Xg0 = load_npy(eg.delta_path).astype(np.float32)

            if Xc0.shape[0] != len(mc0) or Xg0.shape[0] != len(mg0):
                raise RuntimeError(f"Meta/delta mismatch for {track}/{view} Drug or CRISPR")

            tc = get_target_series(mc0).map(safe_upper)
            tg = get_target_series(mg0).map(safe_upper)
            mc0 = mc0.copy()
            mg0 = mg0.copy()
            mc0["target_std"] = tc
            mg0["target_std"] = tg

            mc, Xc = apply_nan_policy(
                mc0, Xc0,
                policy=nan_policy,
                target_std=mc0["target_std"],
                qc_rows=qc_rows,
                qc_key={"phase": "pairwise", "track": track, "view": view, "modality": "Drug"},
            )
            mg, Xg = apply_nan_policy(
                mg0, Xg0,
                policy=nan_policy,
                target_std=mg0["target_std"],
                qc_rows=qc_rows,
                qc_key={"phase": "pairwise", "track": track, "view": view, "modality": "CRISPR"},
            )

            if len(mc) == 0 or len(mg) == 0:
                qc_rows.append({
                    "phase": "pairwise",
                    "track": track, "view": view,
                    "note": "empty_after_filtering",
                    "n_drug_kept": int(len(mc)),
                    "n_crispr_kept": int(len(mg)),
                })
                continue

            common = sorted(set(mc["target_std"].unique()).intersection(set(mg["target_std"].unique())))
            n_targets_common = 0
            for t in common:
                idx_c = np.where(mc["target_std"].to_numpy() == t)[0]
                idx_g = np.where(mg["target_std"].to_numpy() == t)[0]
                if len(idx_c) < 1 or len(idx_g) < 1:
                    continue

                n_targets_common += 1
                cloud_c = Xc[idx_c]
                cloud_g = Xg[idx_g]

                # centroid cosine always if >=1 each side
                cent_c = cloud_c.mean(axis=0)
                cent_g = cloud_g.mean(axis=0)
                cos = float(robust_cosine_matrix(cent_c[None, :], cent_g[None, :])[0, 0])

                # edist only if >=2 each side; otherwise NaN
                if len(idx_c) >= 2 and len(idx_g) >= 2:
                    ed_mean, ed_sd = energy_distance_bootstrap(
                        cloud_c, cloud_g, rng=rng, n_boot=n_boot, max_n=max_n
                    )
                else:
                    ed_mean, ed_sd = float("nan"), float("nan")

                rows.append({
                    "Track": track,
                    "View": view,
                    "Target": t,
                    "n_chem": int(len(idx_c)),
                    "n_gene": int(len(idx_g)),
                    "centroid_cosine": cos,
                    "edist_mean": ed_mean,
                    "edist_sd": ed_sd,
                })

            qc_rows.append({
                "phase": "pairwise",
                "track": track, "view": view,
                "note": "after_eval",
                "n_common_targets": int(n_targets_common),
            })

    per_target = pd.DataFrame(rows)
    if per_target.empty:
        raise RuntimeError("No pairwise results produced (after NaN filtering). Check manifest & filters.")

    summary = (
        per_target
        .groupby(["Track", "View"], as_index=False)
        .agg(
            n_targets=("Target", "count"),
            mean_centroid_cos=("centroid_cosine", "mean"),
            median_centroid_cos=("centroid_cosine", "median"),
            mean_edist=("edist_mean", "mean"),          # pandas mean skips NaN by default
            median_edist=("edist_mean", "median"),      # pandas median skips NaN by default
        )
        .sort_values(["Track", "View"])
    )

    ensure_dir(out_dir)
    per_t_path = os.path.join(out_dir, "Task3_Pairwise_PerTarget.csv")
    sum_path = os.path.join(out_dir, "Task3_Pairwise_Summary.csv")
    per_target.to_csv(per_t_path, index=False)
    summary.to_csv(sum_path, index=False)

    return per_target, summary


# --------------------------
# Main
# --------------------------
def main():
    ap = argparse.ArgumentParser("Task3 Analysis — Retrieval + Pairwise using Task3 deltas (NaN-robust)")
    ap.add_argument("--task3_out_dir", default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Deltas_K562_CellLevel", type=str)
    ap.add_argument("--out_dir", default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Analysis", type=str)
    ap.add_argument("--views", type=str, default="Standard,Systema", help="Comma-separated: Standard,Systema")
    ap.add_argument("--tracks", type=str, default="ALL", help="Comma-separated track names, or ALL.")
    ap.add_argument("--do_retrieval", action="store_true", help="Run retrieval evaluation")
    ap.add_argument("--do_pairwise", action="store_true", help="Run per-target pairwise alignment")

    # NaN handling
    ap.add_argument("--nan_policy", type=str, default="drop", choices=["drop", "zero"],
                    help="How to handle non-finite embedding rows. drop=recommended; zero=replaces non-finite with 0.0.")

    # pairwise options
    ap.add_argument("--edist_boot", type=int, default=30)
    ap.add_argument("--edist_max_n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    man, entries = load_manifest(args.task3_out_dir)
    idx = index_entries(entries)

    views = [v.strip() for v in args.views.split(",") if v.strip()]
    all_tracks = sorted(set([e.track for e in entries]))
    if args.tracks.strip().upper() == "ALL":
        tracks = all_tracks
    else:
        tracks = [t.strip() for t in args.tracks.split(",") if t.strip()]
        tracks = [t for t in tracks if t in all_tracks]

    if not args.do_retrieval and not args.do_pairwise:
        args.do_retrieval = True
        args.do_pairwise = True

    ensure_dir(args.out_dir)

    qc_rows: List[dict] = []

    if args.do_retrieval:
        directions = [("Drug", "CRISPR"), ("CRISPR", "Drug")]
        run_retrieval(
            index=idx,
            tracks=tracks,
            views=views,
            directions=directions,
            out_dir=args.out_dir,
            nan_policy=args.nan_policy,
            qc_rows=qc_rows,
        )
        print(f"[OK] Retrieval outputs written to: {args.out_dir}")

    if args.do_pairwise:
        run_pairwise(
            index=idx,
            tracks=tracks,
            views=views,
            out_dir=args.out_dir,
            nan_policy=args.nan_policy,
            qc_rows=qc_rows,
            n_boot=int(args.edist_boot),
            max_n=int(args.edist_max_n),
            seed=int(args.seed),
        )
        print(f"[OK] Pairwise outputs written to: {args.out_dir}")

    # Write QC
    if qc_rows:
        qc_df = pd.DataFrame(qc_rows)
        qc_path = os.path.join(args.out_dir, "Task3_QC_NaNFiltering.csv")
        qc_df.to_csv(qc_path, index=False)
        print(f"[OK] QC written: {qc_path}")

    print("\n[Done] Task3 analysis complete.")

if __name__ == "__main__":
    main()
