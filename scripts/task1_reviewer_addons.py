#!/usr/bin/env python3
"""
Task1 reviewer-oriented add-on analyses for M2M-Bench.

This script does not rerun heavy Task1 pipelines. It consumes existing Task1 outputs
and produces additional tables/figures required by reviewer-facing robustness checks:

1) Retrieval dual-report (raw + balanced) with theoretical/random baselines.
2) Effect-size calibration (cross-modality vs within-modality replicate consistency).
3) Protocol mismatch continuous sensitivity (raw + deconfounded partial Spearman).
4) LINCS internal consistency (same cell/target, same cell/target/dose/time).
5) Composition-bias sensitivity (cell-stratified + leave-one-cell-out summaries).
6) Benchmarkability map (context comparability zones).
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, spearmanr


# ---------------------------------------------------------------------
# Section 0. Configuration
# ---------------------------------------------------------------------


@dataclass
class Task1ReviewerAddonsConfig:
    task0_unified_meta_path: str = "./outputs/task0_curated/metadata/unified_meta.csv"
    m1_candidates_path: str = "./outputs/task1/data/m1_candidates.csv"
    m1_matched_pairs_path: str = "./outputs/task1/data/m1_matched_pairs.csv"
    per_pair_path: str = "./outputs/task1/analysis/modality_gap_per_pair.csv"
    retrieval_summary_path: str = "./outputs/task1/retrieval/analysis/retrieval_summary.csv"
    retrieval_null_path: str = "./outputs/task1/retrieval/analysis/retrieval_null_summary.csv"
    retrieval_per_query_path: str = "./outputs/task1/retrieval/analysis/retrieval_per_query.csv"
    context_overlap_path: str = "./outputs/task1_audit/analysis/context_overlap_counts.csv"
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    lincs_tensor_name: str = "LINCS_Engine1_TrainData.pt"
    output_dir: str = "./outputs/task1_reviewer_fixes"
    max_within_pair_samples_per_context: int = 8000
    protocol_n_bins: int = 8
    random_seed: int = 42

    @property
    def lincs_data_path(self) -> str:
        return str(Path(self.processed_dir) / "LINCS_Processed" / self.lincs_tensor_name)

    @property
    def sc_dir(self) -> str:
        return str(Path(self.processed_dir) / "scPerturb_Processed")

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def qc_dir(self) -> str:
        return str(Path(self.output_dir) / "qc")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task1_reviewer_addons.json")


# ---------------------------------------------------------------------
# Section 1. Helpers
# ---------------------------------------------------------------------


def _ensure_dirs(paths: Iterable[str]) -> None:
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)


def _read_csv_or_parquet(path_like: str) -> pd.DataFrame:
    path = Path(path_like)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".parquet":
        try:
            return pd.read_parquet(path)
        except Exception:
            fallback = path.with_suffix(".csv")
            if fallback.exists():
                return pd.read_csv(fallback)
            raise
    return pd.read_csv(path)


def _read_unified_meta_minimal(path_like: str) -> pd.DataFrame:
    """
    Read only columns needed for LINCS internal consistency analysis.
    """
    path = Path(path_like)
    need = ["source_db", "modality", "global_idx", "cell_std", "target_std", "dose_val", "time_val"]
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path, columns=need)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, usecols=lambda c: c in set(need))
    return pd.read_csv(path, usecols=lambda c: c in set(need))


def _harmonic_number(n: int) -> float:
    return float(np.sum(1.0 / np.arange(1, n + 1)))


def _normalize_rows(x: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(x, axis=1, keepdims=True)
    denom = np.clip(denom, 1e-12, None)
    return x / denom


def _sample_within_cosine_mean(
    emb: np.ndarray,
    max_samples: int,
    rng: np.random.Generator,
) -> tuple[float, int]:
    """
    Estimate within-source pairwise cosine mean by random pair sampling.
    """
    n = int(emb.shape[0])
    if n < 2:
        return np.nan, 0

    sample_n = int(min(max_samples, n * (n - 1) // 2))
    x = _normalize_rows(emb.astype(np.float32, copy=False))

    # Random i!=j sampling (approximate without enumerating all O(n^2) pairs)
    i = rng.integers(0, n, size=sample_n * 3, endpoint=False)
    j = rng.integers(0, n, size=sample_n * 3, endpoint=False)
    mask = i != j
    i = i[mask][:sample_n]
    j = j[mask][:sample_n]
    if len(i) < sample_n:
        # Rare fallback for very small n
        needed = sample_n - len(i)
        ii = rng.integers(0, n, size=needed * 5, endpoint=False)
        jj = rng.integers(0, n, size=needed * 5, endpoint=False)
        mm = ii != jj
        i = np.concatenate([i, ii[mm][:needed]])
        j = np.concatenate([j, jj[mm][:needed]])

    cos = np.sum(x[i] * x[j], axis=1)
    return float(np.mean(cos)), int(len(cos))


def _sample_within_cosine_stats(
    emb: np.ndarray,
    max_samples: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    """
    Estimate within-source pairwise cosine summary by random pair sampling.
    Returns (mean, median, n_pairs_sampled).
    """
    n = int(emb.shape[0])
    if n < 2:
        return np.nan, np.nan, 0

    sample_n = int(min(max_samples, n * (n - 1) // 2))
    x = _normalize_rows(emb.astype(np.float32, copy=False))

    i = rng.integers(0, n, size=sample_n * 3, endpoint=False)
    j = rng.integers(0, n, size=sample_n * 3, endpoint=False)
    mask = i != j
    i = i[mask][:sample_n]
    j = j[mask][:sample_n]
    if len(i) < sample_n:
        needed = sample_n - len(i)
        ii = rng.integers(0, n, size=needed * 5, endpoint=False)
        jj = rng.integers(0, n, size=needed * 5, endpoint=False)
        mm = ii != jj
        i = np.concatenate([i, ii[mm][:needed]])
        j = np.concatenate([j, jj[mm][:needed]])

    cos = np.sum(x[i] * x[j], axis=1)
    return float(np.mean(cos)), float(np.median(cos)), int(len(cos))


def _cross_cosine_mean(a: np.ndarray, b: np.ndarray) -> tuple[float, int]:
    if len(a) == 0 or len(b) == 0:
        return np.nan, 0
    a_n = _normalize_rows(a.astype(np.float32, copy=False))
    b_n = _normalize_rows(b.astype(np.float32, copy=False))
    sims = a_n @ b_n.T
    return float(np.mean(sims)), int(sims.size)


def _weighted_mean(df: pd.DataFrame, value_col: str, weight_col: str) -> float:
    sub = df[[value_col, weight_col]].copy()
    sub = sub[np.isfinite(sub[value_col]) & (sub[weight_col] > 0)]
    if sub.empty:
        return np.nan
    return float(np.average(sub[value_col].to_numpy(dtype=float), weights=sub[weight_col].to_numpy(dtype=float)))


def _quantile_bins(series: pd.Series, n_bins: int) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    valid = values.dropna()
    if valid.nunique() < 2:
        return pd.Series([np.nan] * len(series), index=series.index, dtype="object")
    bins = int(min(n_bins, valid.nunique()))
    return pd.qcut(values, q=bins, duplicates="drop")


def _linear_residual(y: np.ndarray, x_design: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(x_design, y, rcond=None)
    return y - x_design @ beta


def _partial_spearman(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    covariate_numeric: list[str] | None = None,
    covariate_categorical: list[str] | None = None,
) -> dict:
    covariate_numeric = covariate_numeric or []
    covariate_categorical = covariate_categorical or []
    use_cols = [x_col, y_col] + covariate_numeric + covariate_categorical
    sub = df[use_cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if len(sub) < 6:
        return {"n_pairs": int(len(sub)), "rho_partial": np.nan, "p_value_partial": np.nan}

    x_rank = sub[x_col].rank(method="average").to_numpy(dtype=float)
    y_rank = sub[y_col].rank(method="average").to_numpy(dtype=float)

    x_parts = [np.ones((len(sub), 1), dtype=float)]
    for c in covariate_numeric:
        v = pd.to_numeric(sub[c], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(v).all() and np.nanstd(v) > 1e-12:
            vv = (v - np.nanmean(v)) / (np.nanstd(v) + 1e-12)
            x_parts.append(vv.reshape(-1, 1))
    if covariate_categorical:
        dummies = pd.get_dummies(sub[covariate_categorical].astype(str), drop_first=True, dtype=float)
        if not dummies.empty:
            x_parts.append(dummies.to_numpy(dtype=float))

    design = np.hstack(x_parts)
    if design.shape[1] <= 1:
        rho, pval = spearmanr(sub[x_col].to_numpy(dtype=float), sub[y_col].to_numpy(dtype=float))
        return {"n_pairs": int(len(sub)), "rho_partial": float(rho), "p_value_partial": float(pval)}

    x_res = _linear_residual(x_rank, design)
    y_res = _linear_residual(y_rank, design)
    if np.std(x_res) < 1e-12 or np.std(y_res) < 1e-12:
        return {"n_pairs": int(len(sub)), "rho_partial": np.nan, "p_value_partial": np.nan}
    rho, pval = pearsonr(x_res, y_res)
    return {"n_pairs": int(len(sub)), "rho_partial": float(rho), "p_value_partial": float(pval)}


def _fetch_embeddings_for_candidates(
    candidates: pd.DataFrame,
    cfg: Task1ReviewerAddonsConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fetch y_delta_gene and y_delta_pathway for all m1 candidates.
    Returned arrays align with candidates row order.
    """
    out_gene = [None] * len(candidates)
    out_path = [None] * len(candidates)

    lincs_mask = candidates["source_db"].eq("LINCS").to_numpy()
    if lincs_mask.any():
        lincs_data = torch.load(cfg.lincs_data_path, map_location="cpu")
        idx = torch.as_tensor(candidates.loc[lincs_mask, "global_idx"].to_numpy(dtype=int), dtype=torch.long)
        gene_sub = lincs_data["y_delta_gene"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
        path_sub = lincs_data["y_delta_pathway"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
        row_pos = np.where(lincs_mask)[0]
        for i, pos in enumerate(row_pos):
            out_gene[pos] = gene_sub[i]
            out_path[pos] = path_sub[i]

    sc_mask = candidates["source_db"].eq("scPerturb").to_numpy()
    if sc_mask.any():
        sc_rows = candidates.loc[sc_mask].copy()
        for chunk_file, grp in sc_rows.groupby("chunk_file", sort=False):
            chunk_path = Path(cfg.sc_dir) / str(chunk_file)
            if not chunk_path.exists():
                raise FileNotFoundError(chunk_path)
            chunk = torch.load(chunk_path, map_location="cpu")
            idx = torch.as_tensor(grp["chunk_idx"].to_numpy(dtype=int), dtype=torch.long)
            gene_sub = chunk["y_delta_gene"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
            path_sub = chunk["y_delta_pathway"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
            for i, pos in enumerate(grp.index.to_numpy(dtype=int)):
                out_gene[pos] = gene_sub[i]
                out_path[pos] = path_sub[i]

    if any(v is None for v in out_gene) or any(v is None for v in out_path):
        raise RuntimeError("Failed to fetch all candidate embeddings.")

    return np.stack(out_gene, axis=0), np.stack(out_path, axis=0)


# ---------------------------------------------------------------------
# Section 2. Analysis modules
# ---------------------------------------------------------------------


def build_retrieval_dual_report(
    retrieval_summary: pd.DataFrame,
    retrieval_null: pd.DataFrame,
    retrieval_per_query: pd.DataFrame,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    key = ["track", "direction", "modality"]
    merged = retrieval_summary.merge(
        retrieval_null[
            key
            + [
                "null_mrr_mean",
                "null_top1_mean",
                "null_top5_mean",
                "null_top10_mean",
                "null_mrr_balanced_mean",
                "null_top1_balanced_mean",
                "null_top5_balanced_mean",
                "null_top10_balanced_mean",
            ]
        ],
        on=key,
        how="left",
    )

    # Theoretical random baseline for balanced setting (uniform single-positive rank)
    q_bal = retrieval_per_query[retrieval_per_query["balanced_valid"].astype(bool)].copy()
    if q_bal.empty:
        balanced_n_gallery = np.nan
        balanced_n_true = np.nan
    else:
        balanced_n_gallery = int(pd.to_numeric(q_bal["balanced_n_gallery_eval"], errors="coerce").dropna().median())
        balanced_n_true = int(pd.to_numeric(q_bal["balanced_n_true_eval"], errors="coerce").dropna().median())

    if np.isfinite(balanced_n_gallery) and int(balanced_n_true) == 1:
        n = int(balanced_n_gallery)
        theo = {
            "theoretical_random_mrr_balanced": _harmonic_number(n) / n,
            "theoretical_random_top1_balanced": 1.0 / n,
            "theoretical_random_top5_balanced": min(5.0 / n, 1.0),
            "theoretical_random_top10_balanced": min(10.0 / n, 1.0),
        }
    else:
        theo = {
            "theoretical_random_mrr_balanced": np.nan,
            "theoretical_random_top1_balanced": np.nan,
            "theoretical_random_top5_balanced": np.nan,
            "theoretical_random_top10_balanced": np.nan,
        }

    for col, val in theo.items():
        merged[col] = val

    merged["mrr_balanced_over_theoretical"] = merged["mrr_balanced"] / merged["theoretical_random_mrr_balanced"]
    merged["top1_balanced_over_theoretical"] = merged["top1_balanced"] / merged["theoretical_random_top1_balanced"]

    merged.to_csv(out_dir / "retrieval_dual_report.csv", index=False)

    # Compact section-level summary
    compact = (
        merged.groupby(["modality"], as_index=False)
        .agg(
            mean_mrr_raw=("mrr", "mean"),
            mean_mrr_balanced=("mrr_balanced", "mean"),
            mean_top1_raw=("top1", "mean"),
            mean_top1_balanced=("top1_balanced", "mean"),
            mean_mrr_balanced_over_theoretical=("mrr_balanced_over_theoretical", "mean"),
        )
        .sort_values("modality")
    )
    compact.to_csv(out_dir / "retrieval_dual_compact.csv", index=False)
    return merged, compact


def build_effect_size_calibration(
    candidates: pd.DataFrame,
    y_gene: np.ndarray,
    y_path: np.ndarray,
    per_pair: pd.DataFrame,
    cfg: Task1ReviewerAddonsConfig,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build context-level within-source and cross-source cosine consistency estimates.
    """
    rng = np.random.default_rng(cfg.random_seed)
    key_cols = ["cell_std", "modality", "target_std"]
    candidates = candidates.copy().reset_index(drop=True)
    candidates["embed_idx"] = np.arange(len(candidates), dtype=np.int64)
    emb_map = {"gene": y_gene, "path": y_path}

    rows: list[dict] = []
    for keys, grp in candidates.groupby(key_cols, sort=False):
        cell, modality, target = keys
        lincs_idx = grp.loc[grp["source_db"] == "LINCS", "embed_idx"].to_numpy(dtype=int)
        sc_idx = grp.loc[grp["source_db"] == "scPerturb", "embed_idx"].to_numpy(dtype=int)
        for track, emb in emb_map.items():
            row = {
                "cell_std": str(cell),
                "modality": str(modality),
                "target_std": str(target),
                "track": str(track),
                "n_lincs": int(len(lincs_idx)),
                "n_sc": int(len(sc_idx)),
            }
            if len(lincs_idx) >= 2:
                m, n_samp = _sample_within_cosine_mean(
                    emb[lincs_idx],
                    max_samples=cfg.max_within_pair_samples_per_context,
                    rng=rng,
                )
                row["within_lincs_mean"] = m
                row["n_within_lincs_samples"] = n_samp
            else:
                row["within_lincs_mean"] = np.nan
                row["n_within_lincs_samples"] = 0

            if len(sc_idx) >= 2:
                m, n_samp = _sample_within_cosine_mean(
                    emb[sc_idx],
                    max_samples=cfg.max_within_pair_samples_per_context,
                    rng=rng,
                )
                row["within_sc_mean"] = m
                row["n_within_sc_samples"] = n_samp
            else:
                row["within_sc_mean"] = np.nan
                row["n_within_sc_samples"] = 0

            m_cross, n_cross = _cross_cosine_mean(emb[lincs_idx], emb[sc_idx])
            row["cross_allpairs_mean"] = m_cross
            row["n_cross_allpairs"] = n_cross
            rows.append(row)

    ctx = pd.DataFrame(rows)
    ctx.to_csv(out_dir / "effect_calibration_context_level.csv", index=False)

    # Match-aware cross metric from Task1 matched-pair outputs
    track_to_col = {"gene": "cosine_gene", "path": "cosine_path"}
    summary_rows: list[dict] = []
    for modality, mod_df in ctx.groupby("modality", sort=False):
        pp_mod = per_pair[per_pair["modality"] == modality].copy()
        for track in ["gene", "path"]:
            sub = mod_df[mod_df["track"] == track].copy()
            within_lincs = _weighted_mean(sub, "within_lincs_mean", "n_within_lincs_samples")
            within_sc = _weighted_mean(sub, "within_sc_mean", "n_within_sc_samples")
            cross_allpairs = _weighted_mean(sub, "cross_allpairs_mean", "n_cross_allpairs")
            col = track_to_col[track]
            cross_matched = float(pp_mod[col].mean()) if col in pp_mod.columns and len(pp_mod) > 0 else np.nan
            within_avg = np.nanmean([within_lincs, within_sc])
            summary_rows.append(
                {
                    "modality": str(modality),
                    "track": track,
                    "within_lincs_weighted_mean": within_lincs,
                    "within_sc_weighted_mean": within_sc,
                    "within_avg_weighted_mean": float(within_avg) if np.isfinite(within_avg) else np.nan,
                    "cross_allpairs_weighted_mean": cross_allpairs,
                    "cross_matched_mean": cross_matched,
                    "cross_matched_over_within_avg": float(cross_matched / within_avg) if np.isfinite(cross_matched) and np.isfinite(within_avg) and abs(within_avg) > 1e-12 else np.nan,
                    "cross_allpairs_over_within_avg": float(cross_allpairs / within_avg) if np.isfinite(cross_allpairs) and np.isfinite(within_avg) and abs(within_avg) > 1e-12 else np.nan,
                    "n_contexts": int(len(sub)),
                    "n_matched_pairs": int(len(pp_mod)),
                }
            )

    summary = pd.DataFrame(summary_rows).sort_values(["modality", "track"])
    summary.to_csv(out_dir / "effect_calibration_summary.csv", index=False)
    return ctx, summary


def build_protocol_continuous_sensitivity(
    per_pair: pd.DataFrame,
    cfg: Task1ReviewerAddonsConfig,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    chem = per_pair[per_pair["modality"] == "Chemical"].copy()
    chem["dose_logdiff"] = pd.to_numeric(chem["dose_logdiff"], errors="coerce")
    chem["time_absdiff"] = pd.to_numeric(chem["time_absdiff"], errors="coerce")
    chem["cosine_gene"] = pd.to_numeric(chem["cosine_gene"], errors="coerce")
    chem["cosine_path"] = pd.to_numeric(chem["cosine_path"], errors="coerce")

    # Spearman trend tests
    test_rows: list[dict] = []
    for x_col in ["dose_logdiff", "time_absdiff"]:
        for y_col in ["cosine_gene", "cosine_path"]:
            sub = chem[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(sub) >= 3:
                rho, pval = spearmanr(sub[x_col].to_numpy(), sub[y_col].to_numpy())
                test_rows.append(
                    {
                        "x_var": x_col,
                        "y_var": y_col,
                        "n_pairs": int(len(sub)),
                        "spearman_rho": float(rho),
                        "p_value": float(pval),
                    }
                )
    tests = pd.DataFrame(test_rows)
    tests.to_csv(out_dir / "protocol_continuous_spearman.csv", index=False)

    # Partial Spearman (deconfounded by cell/target; optionally plus other mismatch axis)
    partial_rows: list[dict] = []
    for x_col in ["dose_logdiff", "time_absdiff"]:
        for y_col in ["cosine_gene", "cosine_path"]:
            other = "time_absdiff" if x_col == "dose_logdiff" else "dose_logdiff"
            specs = [
                ("cell_target", []),
                ("cell_target_plus_other_mismatch", [other]),
            ]
            for control_spec, cov_num in specs:
                res = _partial_spearman(
                    df=chem,
                    x_col=x_col,
                    y_col=y_col,
                    covariate_numeric=cov_num,
                    covariate_categorical=["cell_std", "target_std"],
                )
                partial_rows.append(
                    {
                        "x_var": x_col,
                        "y_var": y_col,
                        "control_spec": control_spec,
                        "n_pairs": int(res["n_pairs"]),
                        "partial_spearman_rho": res["rho_partial"],
                        "p_value": res["p_value_partial"],
                    }
                )
    partial = pd.DataFrame(partial_rows)
    partial.to_csv(out_dir / "protocol_continuous_partial_spearman.csv", index=False)

    # 1D quantile-bin trend curves
    curve_rows: list[dict] = []
    for x_col in ["dose_logdiff", "time_absdiff"]:
        bins = _quantile_bins(chem[x_col], n_bins=cfg.protocol_n_bins)
        sub = chem.copy()
        sub["bin"] = bins
        for y_col in ["cosine_gene", "cosine_path"]:
            grouped = (
                sub.dropna(subset=["bin", y_col, x_col])
                .groupby("bin", observed=True)
                .agg(
                    n_pairs=(y_col, "size"),
                    x_center=(x_col, "median"),
                    y_mean=(y_col, "mean"),
                    y_median=(y_col, "median"),
                )
                .reset_index()
            )
            grouped["x_var"] = x_col
            grouped["y_var"] = y_col
            grouped["bin_label"] = grouped["bin"].astype(str)
            curve_rows.extend(grouped[["x_var", "y_var", "bin_label", "x_center", "n_pairs", "y_mean", "y_median"]].to_dict("records"))
    curves = pd.DataFrame(curve_rows)
    curves.to_csv(out_dir / "protocol_continuous_curves.csv", index=False)

    # 2D dose-time grid
    chem2 = chem.copy()
    chem2["dose_bin"] = _quantile_bins(chem2["dose_logdiff"], n_bins=5)
    chem2["time_bin"] = _quantile_bins(chem2["time_absdiff"], n_bins=5)
    heat_rows: list[dict] = []
    for y_col in ["cosine_gene", "cosine_path"]:
        g = (
            chem2.dropna(subset=["dose_bin", "time_bin", y_col])
            .groupby(["dose_bin", "time_bin"], observed=True)
            .agg(
                n_pairs=(y_col, "size"),
                mean_cosine=(y_col, "mean"),
            )
            .reset_index()
        )
        g["track"] = "gene" if y_col.endswith("gene") else "path"
        g["dose_bin"] = g["dose_bin"].astype(str)
        g["time_bin"] = g["time_bin"].astype(str)
        heat_rows.extend(g[["track", "dose_bin", "time_bin", "n_pairs", "mean_cosine"]].to_dict("records"))
    heat = pd.DataFrame(heat_rows)
    heat.to_csv(out_dir / "protocol_continuous_heatmap.csv", index=False)

    # Optional quick QC plots for manuscript prep
    try:
        import matplotlib.pyplot as plt

        for y_var in ["cosine_gene", "cosine_path"]:
            fig, ax = plt.subplots(figsize=(6, 4))
            sub = curves[curves["y_var"] == y_var].copy()
            for x_var in ["dose_logdiff", "time_absdiff"]:
                ss = sub[sub["x_var"] == x_var].sort_values("x_center")
                if ss.empty:
                    continue
                ax.plot(ss["x_center"], ss["y_mean"], marker="o", label=x_var)
            ax.axhline(0.0, color="gray", linewidth=1, linestyle="--")
            ax.set_xlabel("Mismatch (bin center)")
            ax.set_ylabel(f"Mean {y_var}")
            ax.set_title(f"Protocol mismatch sensitivity ({y_var})")
            ax.legend(frameon=False)
            fig.tight_layout()
            fig.savefig(out_dir / f"protocol_continuous_curve_{y_var}.png", dpi=160)
            plt.close(fig)
    except Exception:
        pass
    return tests, partial, curves, heat


def build_lincs_internal_consistency(
    unified_meta: pd.DataFrame,
    cfg: Task1ReviewerAddonsConfig,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Internal LINCS consistency over the full LINCS chemical pool:
      - same (cell, target)
      - same (cell, target, dose, time)
    """
    rng = np.random.default_rng(cfg.random_seed)
    meta = unified_meta.copy()
    meta["global_idx"] = pd.to_numeric(meta["global_idx"], errors="coerce")
    meta = meta[
        (meta["source_db"] == "LINCS")
        & (meta["modality"] == "Chemical")
        & meta["global_idx"].notna()
    ].copy()
    meta["global_idx"] = meta["global_idx"].astype(int)
    meta = meta.dropna(subset=["cell_std", "target_std", "dose_val", "time_val"])
    meta = meta.reset_index(drop=True)

    lincs_data = torch.load(cfg.lincs_data_path, map_location="cpu")
    idx = torch.as_tensor(meta["global_idx"].to_numpy(dtype=int), dtype=torch.long)
    y_gene = lincs_data["y_delta_gene"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
    y_path = lincs_data["y_delta_pathway"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
    meta["embed_idx"] = np.arange(len(meta), dtype=np.int64)

    def _compute(group_cols: list[str], level_name: str) -> pd.DataFrame:
        rows: list[dict] = []
        for keys, grp in meta.groupby(group_cols, sort=False):
            emb_idx = grp["embed_idx"].to_numpy(dtype=int)
            n_rows = int(len(emb_idx))
            if n_rows < 2:
                continue
            g_mean, g_median, g_n = _sample_within_cosine_stats(
                y_gene[emb_idx], max_samples=cfg.max_within_pair_samples_per_context, rng=rng
            )
            p_mean, p_median, p_n = _sample_within_cosine_stats(
                y_path[emb_idx], max_samples=cfg.max_within_pair_samples_per_context, rng=rng
            )
            row = {
                "level": level_name,
                "n_rows": n_rows,
                "n_pairs_sampled_gene": int(g_n),
                "n_pairs_sampled_path": int(p_n),
                "cosine_gene_mean": g_mean,
                "cosine_gene_median_sample": g_median,
                "cosine_path_mean": p_mean,
                "cosine_path_median_sample": p_median,
            }
            if isinstance(keys, tuple):
                for k, v in zip(group_cols, keys):
                    row[k] = v
            else:
                row[group_cols[0]] = keys
            rows.append(row)
        out = pd.DataFrame(rows)
        ordered = (
            ["level"]
            + group_cols
            + [
                "n_rows",
                "n_pairs_sampled_gene",
                "n_pairs_sampled_path",
                "cosine_gene_mean",
                "cosine_gene_median_sample",
                "cosine_path_mean",
                "cosine_path_median_sample",
            ]
        )
        if out.empty:
            return pd.DataFrame(columns=ordered)
        return out[ordered]

    by_cell_target = _compute(["cell_std", "target_std"], "cell_target")
    by_cell_target_dt = _compute(["cell_std", "target_std", "dose_val", "time_val"], "cell_target_dose_time")
    by_cell_target.to_csv(out_dir / "lincs_internal_consistency_cell_target.csv", index=False)
    by_cell_target_dt.to_csv(out_dir / "lincs_internal_consistency_cell_target_dose_time.csv", index=False)

    def _summary_row(name: str, df: pd.DataFrame) -> dict:
        if df.empty:
            return {
                "name": name,
                "n_groups": 0,
                "n_rows_total": 0,
                "n_pairs_sampled_total": 0,
                "weighted_cosine_gene": np.nan,
                "weighted_cosine_path": np.nan,
                "median_cosine_gene": np.nan,
                "median_cosine_path": np.nan,
                "q25_gene": np.nan,
                "q75_gene": np.nan,
                "q25_path": np.nan,
                "q75_path": np.nan,
            }
        return {
            "name": name,
            "n_groups": int(len(df)),
            "n_rows_total": int(df["n_rows"].sum()),
            "n_pairs_sampled_total": int(df["n_pairs_sampled_gene"].sum()),
            "weighted_cosine_gene": _weighted_mean(df, "cosine_gene_mean", "n_pairs_sampled_gene"),
            "weighted_cosine_path": _weighted_mean(df, "cosine_path_mean", "n_pairs_sampled_path"),
            "median_cosine_gene": float(df["cosine_gene_mean"].median()),
            "median_cosine_path": float(df["cosine_path_mean"].median()),
            "q25_gene": float(df["cosine_gene_mean"].quantile(0.25)),
            "q75_gene": float(df["cosine_gene_mean"].quantile(0.75)),
            "q25_path": float(df["cosine_path_mean"].quantile(0.25)),
            "q75_path": float(df["cosine_path_mean"].quantile(0.75)),
        }

    summary = pd.DataFrame(
        [
            _summary_row("same_cell_target", by_cell_target),
            _summary_row("same_cell_target_dose_time", by_cell_target_dt),
        ]
    )
    summary.to_csv(out_dir / "lincs_internal_consistency_summary.csv", index=False)
    return by_cell_target, by_cell_target_dt, summary


def build_cell_bias_sensitivity(
    per_pair: pd.DataFrame,
    retrieval_per_query: pd.DataFrame,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Pairwise: cell-stratified
    cell_pair = (
        per_pair.groupby(["cell_std", "modality"], as_index=False)
        .agg(
            n_pairs=("cosine_gene", "size"),
            mean_cosine_gene=("cosine_gene", "mean"),
            mean_cosine_path=("cosine_path", "mean"),
        )
        .sort_values(["modality", "n_pairs"], ascending=[True, False])
    )
    cell_pair.to_csv(out_dir / "cell_stratified_pairwise.csv", index=False)

    # Pairwise: leave-one-cell-out
    cells = sorted(per_pair["cell_std"].dropna().astype(str).unique().tolist())
    loo_pair_rows: list[dict] = []
    for excluded in [None] + cells:
        sub = per_pair if excluded is None else per_pair[per_pair["cell_std"] != excluded]
        for modality in ["All", "Chemical", "Genetic"]:
            ss = sub if modality == "All" else sub[sub["modality"] == modality]
            if ss.empty:
                continue
            loo_pair_rows.append(
                {
                    "excluded_cell": "ALL" if excluded is None else excluded,
                    "modality": modality,
                    "n_pairs": int(len(ss)),
                    "mean_cosine_gene": float(ss["cosine_gene"].mean()),
                    "mean_cosine_path": float(ss["cosine_path"].mean()),
                }
            )
    loo_pair = pd.DataFrame(loo_pair_rows)
    loo_pair.to_csv(out_dir / "leave_one_cell_out_pairwise.csv", index=False)

    # Retrieval: cell-stratified (query cell)
    rq = retrieval_per_query.copy()
    rq["balanced_valid"] = rq["balanced_valid"].fillna(False).astype(bool)
    cell_ret = (
        rq.groupby(["query_cell_std", "track", "direction", "modality"], as_index=False)
        .agg(
            n_query=("mrr", "size"),
            mrr=("mrr", "mean"),
            top1=("success_top1", "mean"),
            mrr_balanced=("mrr_balanced", "mean"),
            top1_balanced=("success_top1_balanced", "mean"),
        )
        .sort_values(["n_query"], ascending=False)
    )
    cell_ret.to_csv(out_dir / "cell_stratified_retrieval.csv", index=False)

    # Retrieval: leave-one-cell-out (query-side sensitivity)
    cells_q = sorted(rq["query_cell_std"].dropna().astype(str).unique().tolist())
    loo_ret_rows: list[dict] = []
    for excluded in [None] + cells_q:
        sub = rq if excluded is None else rq[rq["query_cell_std"] != excluded]
        grouped = sub.groupby(["track", "direction", "modality"], as_index=False)
        agg = grouped.agg(
            n_query=("mrr", "size"),
            mrr=("mrr", "mean"),
            top1=("success_top1", "mean"),
            mrr_balanced=("mrr_balanced", "mean"),
            top1_balanced=("success_top1_balanced", "mean"),
        )
        agg["excluded_cell"] = "ALL" if excluded is None else excluded
        loo_ret_rows.extend(agg.to_dict("records"))
    loo_ret = pd.DataFrame(loo_ret_rows)
    loo_ret.to_csv(out_dir / "leave_one_cell_out_retrieval.csv", index=False)

    return cell_pair, loo_pair, cell_ret, loo_ret


def build_strict_subset_composition(per_pair: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    df = per_pair.copy()
    df["is_strict_protocol_pair"] = df["is_strict_protocol_pair"].fillna(False).astype(bool)
    summary = (
        df.groupby(["modality", "is_strict_protocol_pair"], as_index=False)
        .agg(n_pairs=("pair_id_task1", "size"))
        .sort_values(["is_strict_protocol_pair", "modality"])
    )
    totals = summary.groupby("is_strict_protocol_pair")["n_pairs"].transform("sum")
    summary["fraction_within_strict_flag"] = summary["n_pairs"] / totals
    summary.to_csv(out_dir / "strict_subset_composition.csv", index=False)
    return summary


def build_benchmarkability_map(
    context_overlap: pd.DataFrame,
    per_pair: pd.DataFrame,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    key_cols = ["cell_std", "modality", "target_std"]
    ctx = context_overlap[key_cols + ["LINCS", "scPerturb", "is_overlap_context"]].copy()
    ctx["is_overlap_context"] = ctx["is_overlap_context"].astype(bool)

    matched_ctx = (
        per_pair.groupby(key_cols, as_index=False)
        .agg(
            n_matched_pairs=("pair_id_task1", "size"),
            strict_pair_rate=("is_strict_protocol_pair", "mean"),
            has_strict_pair=("is_strict_protocol_pair", "max"),
        )
    )
    merged = ctx.merge(matched_ctx, on=key_cols, how="left")
    merged["n_matched_pairs"] = merged["n_matched_pairs"].fillna(0).astype(int)
    merged["strict_pair_rate"] = merged["strict_pair_rate"].fillna(0.0)
    merged["has_strict_pair"] = merged["has_strict_pair"].fillna(False).astype(bool)

    def _zone(row: pd.Series) -> str:
        if not bool(row["is_overlap_context"]):
            return "Z0_NoCrossSourceOverlap"
        if str(row["modality"]) == "Genetic":
            return "Z3_GeneticExactComparable"
        if str(row["modality"]) == "Chemical" and bool(row["has_strict_pair"]):
            return "Z2_ChemicalStrictComparable"
        if str(row["modality"]) == "Chemical":
            return "Z1_ChemicalNearestOnly"
        return "Zx_Other"

    merged["benchmark_zone"] = merged.apply(_zone, axis=1)
    merged.to_csv(out_dir / "benchmarkability_context_map.csv", index=False)

    summary = (
        merged.groupby(["benchmark_zone", "modality"], as_index=False)
        .agg(
            n_contexts=("target_std", "size"),
            n_lincs_rows=("LINCS", "sum"),
            n_sc_rows=("scPerturb", "sum"),
            n_matched_pairs=("n_matched_pairs", "sum"),
        )
        .sort_values(["benchmark_zone", "modality"])
    )
    summary.to_csv(out_dir / "benchmarkability_zone_summary.csv", index=False)
    return merged, summary


# ---------------------------------------------------------------------
# Section 3. Main
# ---------------------------------------------------------------------


def run_task1_reviewer_addons(cfg: Task1ReviewerAddonsConfig) -> None:
    t0 = time.time()
    _ensure_dirs([cfg.output_dir, cfg.analysis_dir, cfg.qc_dir])

    # Load inputs
    candidates = _read_csv_or_parquet(cfg.m1_candidates_path)
    unified_meta = _read_unified_meta_minimal(cfg.task0_unified_meta_path)
    per_pair = _read_csv_or_parquet(cfg.per_pair_path)
    retrieval_summary = _read_csv_or_parquet(cfg.retrieval_summary_path)
    retrieval_null = _read_csv_or_parquet(cfg.retrieval_null_path)
    retrieval_per_query = _read_csv_or_parquet(cfg.retrieval_per_query_path)
    context_overlap = _read_csv_or_parquet(cfg.context_overlap_path)

    # Harmonize minimal columns
    if "balanced_valid" not in retrieval_per_query.columns:
        retrieval_per_query["balanced_valid"] = False

    # Retrieval dual report
    retrieval_dual, retrieval_compact = build_retrieval_dual_report(
        retrieval_summary=retrieval_summary,
        retrieval_null=retrieval_null,
        retrieval_per_query=retrieval_per_query,
        out_dir=Path(cfg.analysis_dir),
    )

    # Effect-size calibration (requires candidate embeddings)
    y_gene, y_path = _fetch_embeddings_for_candidates(candidates.reset_index(drop=True), cfg=cfg)
    ctx_calib, summary_calib = build_effect_size_calibration(
        candidates=candidates,
        y_gene=y_gene,
        y_path=y_path,
        per_pair=per_pair,
        cfg=cfg,
        out_dir=Path(cfg.analysis_dir),
    )

    # Protocol continuous sensitivity
    tests, partial_tests, curves, heat = build_protocol_continuous_sensitivity(
        per_pair=per_pair,
        cfg=cfg,
        out_dir=Path(cfg.analysis_dir),
    )

    # LINCS internal consistency (within-bulk repeatability sanity)
    lincs_ctx, lincs_ctx_dt, lincs_summary = build_lincs_internal_consistency(
        unified_meta=unified_meta,
        cfg=cfg,
        out_dir=Path(cfg.analysis_dir),
    )

    # Cell-bias sensitivity
    cell_pair, loo_pair, cell_ret, loo_ret = build_cell_bias_sensitivity(
        per_pair=per_pair,
        retrieval_per_query=retrieval_per_query,
        out_dir=Path(cfg.analysis_dir),
    )

    # Strict-subset composition (for explicit denominator/caveat reporting)
    strict_comp = build_strict_subset_composition(per_pair=per_pair, out_dir=Path(cfg.analysis_dir))

    # Benchmarkability map
    map_df, map_summary = build_benchmarkability_map(
        context_overlap=context_overlap,
        per_pair=per_pair,
        out_dir=Path(cfg.analysis_dir),
    )

    # Quick machine-readable key summary
    key_summary = pd.DataFrame(
        [
            {"metric": "n_task1_candidates", "value": int(len(candidates))},
            {"metric": "n_task1_matched_pairs", "value": int(len(per_pair))},
            {"metric": "n_retrieval_groups", "value": int(len(retrieval_dual))},
            {"metric": "n_calibration_context_rows", "value": int(len(ctx_calib))},
            {"metric": "n_protocol_spearman_tests", "value": int(len(tests))},
            {"metric": "n_protocol_partial_tests", "value": int(len(partial_tests))},
            {"metric": "n_lincs_internal_cell_target_rows", "value": int(len(lincs_ctx))},
            {"metric": "n_lincs_internal_cell_target_dose_time_rows", "value": int(len(lincs_ctx_dt))},
            {"metric": "n_benchmarkability_contexts", "value": int(len(map_df))},
            {"metric": "n_strict_subset_comp_rows", "value": int(len(strict_comp))},
            {"metric": "n_leave_one_cell_out_pair_rows", "value": int(len(loo_pair))},
            {"metric": "n_leave_one_cell_out_retrieval_rows", "value": int(len(loo_ret))},
        ]
    )
    key_summary.to_csv(Path(cfg.analysis_dir) / "task1_reviewer_addons_key_summary.csv", index=False)

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "retrieval_groups": int(len(retrieval_dual)),
            "calibration_rows": int(len(ctx_calib)),
            "protocol_tests": int(len(tests)),
            "protocol_partial_tests": int(len(partial_tests)),
            "lincs_internal_summary_rows": int(len(lincs_summary)),
            "cell_pair_rows": int(len(cell_pair)),
            "map_context_rows": int(len(map_df)),
        },
        "outputs": {
            "analysis_dir": cfg.analysis_dir,
            "qc_dir": cfg.qc_dir,
        },
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")
    print(f"[Task1 Reviewer Addons] done. Outputs: {cfg.analysis_dir}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Task1 reviewer add-on analyses")
    parser.add_argument("--task0-unified-meta-path", type=str, default=Task1ReviewerAddonsConfig.task0_unified_meta_path)
    parser.add_argument("--m1-candidates-path", type=str, default=Task1ReviewerAddonsConfig.m1_candidates_path)
    parser.add_argument("--m1-matched-pairs-path", type=str, default=Task1ReviewerAddonsConfig.m1_matched_pairs_path)
    parser.add_argument("--per-pair-path", type=str, default=Task1ReviewerAddonsConfig.per_pair_path)
    parser.add_argument("--retrieval-summary-path", type=str, default=Task1ReviewerAddonsConfig.retrieval_summary_path)
    parser.add_argument("--retrieval-null-path", type=str, default=Task1ReviewerAddonsConfig.retrieval_null_path)
    parser.add_argument("--retrieval-per-query-path", type=str, default=Task1ReviewerAddonsConfig.retrieval_per_query_path)
    parser.add_argument("--context-overlap-path", type=str, default=Task1ReviewerAddonsConfig.context_overlap_path)
    parser.add_argument("--processed-dir", type=str, default=Task1ReviewerAddonsConfig.processed_dir)
    parser.add_argument("--lincs-tensor-name", type=str, default=Task1ReviewerAddonsConfig.lincs_tensor_name)
    parser.add_argument("--output-dir", type=str, default=Task1ReviewerAddonsConfig.output_dir)
    parser.add_argument("--max-within-pair-samples-per-context", type=int, default=Task1ReviewerAddonsConfig.max_within_pair_samples_per_context)
    parser.add_argument("--protocol-n-bins", type=int, default=Task1ReviewerAddonsConfig.protocol_n_bins)
    parser.add_argument("--random-seed", type=int, default=Task1ReviewerAddonsConfig.random_seed)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task1ReviewerAddonsConfig(
        task0_unified_meta_path=args.task0_unified_meta_path,
        m1_candidates_path=args.m1_candidates_path,
        m1_matched_pairs_path=args.m1_matched_pairs_path,
        per_pair_path=args.per_pair_path,
        retrieval_summary_path=args.retrieval_summary_path,
        retrieval_null_path=args.retrieval_null_path,
        retrieval_per_query_path=args.retrieval_per_query_path,
        context_overlap_path=args.context_overlap_path,
        processed_dir=args.processed_dir,
        lincs_tensor_name=args.lincs_tensor_name,
        output_dir=args.output_dir,
        max_within_pair_samples_per_context=args.max_within_pair_samples_per_context,
        protocol_n_bins=args.protocol_n_bins,
        random_seed=args.random_seed,
    )
    run_task1_reviewer_addons(cfg)


if __name__ == "__main__":
    main()
