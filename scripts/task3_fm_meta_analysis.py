#!/usr/bin/env python3
"""
Task3 foundation-model meta-analysis for M2M-Bench.

Goal
----
Relate model performance to two attributes requested for reviewer-facing analyses:
1) model size class,
2) perturbation-trained vs non-perturbation-trained.

This module is not a pure leaderboard. The aggregate score
(`mean_scaled_best`) is explicitly defined and exported.

Input
-----
Primary input is viz_scoreboard_long.csv from Task3 unified outputs.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


# ---------------------------------------------------------------------
# Section 0. Configuration
# ---------------------------------------------------------------------


@dataclass
class Task3MetaConfig:
    scoreboard_long_path: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified/viz_scoreboard_long.csv"
    output_dir: str = "./outputs/task3_meta"
    model_meta_map_path: str | None = None

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task3_meta.json")


# ---------------------------------------------------------------------
# Section 1. Model metadata
# ---------------------------------------------------------------------


DEFAULT_MODEL_META = {
    # Feature/PCA baselines
    "Gene": {
        "model_family": "feature_space",
        "size_class": "baseline",
        "size_code": 0,
        "perturbation_trained": False,
    },
    "Pathway_HALLMARK": {
        "model_family": "feature_space",
        "size_class": "baseline",
        "size_code": 0,
        "perturbation_trained": False,
    },
    "PCA50": {
        "model_family": "linear_baseline",
        "size_class": "baseline",
        "size_code": 0,
        "perturbation_trained": False,
    },
    "PCA100": {
        "model_family": "linear_baseline",
        "size_class": "baseline",
        "size_code": 0,
        "perturbation_trained": False,
    },
    "PCA200": {
        "model_family": "linear_baseline",
        "size_class": "baseline",
        "size_code": 0,
        "perturbation_trained": False,
    },
    # General foundation models
    "scGPT": {
        "model_family": "foundation_model",
        "size_class": "fm_standard",
        "size_code": 1,
        "perturbation_trained": False,
    },
    "scBERT": {
        "model_family": "foundation_model",
        "size_class": "fm_standard",
        "size_code": 1,
        "perturbation_trained": False,
    },
    "Geneformer": {
        "model_family": "foundation_model",
        "size_class": "fm_standard",
        "size_code": 1,
        "perturbation_trained": False,
    },
    "scFoundation": {
        "model_family": "foundation_model",
        "size_class": "fm_standard",
        "size_code": 1,
        "perturbation_trained": False,
    },
    "UCE": {
        "model_family": "foundation_model",
        "size_class": "fm_standard",
        "size_code": 1,
        "perturbation_trained": False,
    },
    # Perturbation-oriented models
    "STATE": {
        "model_family": "foundation_model",
        "size_class": "fm_standard",
        "size_code": 1,
        "perturbation_trained": True,
    },
    "TahoeX1_3b": {
        "model_family": "foundation_model",
        "size_class": "fm_large",
        "size_code": 2,
        "perturbation_trained": True,
    },
}


def _load_model_meta(optional_csv: str | None) -> pd.DataFrame:
    base = pd.DataFrame([{"Track": k, **v} for k, v in DEFAULT_MODEL_META.items()])
    base["meta_source"] = "default_internal"

    if optional_csv is None:
        return base

    ext = pd.read_csv(optional_csv)
    required = {"Track", "model_family", "size_class", "size_code", "perturbation_trained"}
    missing = required - set(ext.columns)
    if missing:
        raise ValueError(f"External model meta missing columns: {sorted(missing)}")
    ext = ext[list(required)].copy()
    ext["meta_source"] = "external_csv"
    ext = ext.drop_duplicates(subset=["Track"], keep="first")

    merged = base.set_index("Track")
    ext = ext.set_index("Track")
    merged.update(ext)
    merged = merged.reset_index()
    return merged


# ---------------------------------------------------------------------
# Section 2. Analyses
# ---------------------------------------------------------------------


def build_model_scoreboard(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate model scores from per-metric slices.
    Use Scaled_Best as the primary comparable metric across tasks.
    """
    out = (
        df.groupby(
            ["Track", "model_family", "size_class", "size_code", "perturbation_trained"],
            as_index=False,
        )
        .agg(
            n_metric_slices=("Metric", "size"),
            n_metric_groups=("MetricGroup", "nunique"),
            n_views=("View", "nunique"),
            mean_scaled_best=("Scaled_Best", "mean"),
            median_scaled_best=("Scaled_Best", "median"),
            mean_rank_best=("Rank_Best", "mean"),
            median_rank_best=("Rank_Best", "median"),
            mean_value=("Value", "mean"),
        )
        .sort_values(["mean_scaled_best", "mean_rank_best"], ascending=[False, True])
    )
    out["score_definition_id"] = "equal_weight_metric_slice_mean_scaled_best"
    return out


def build_perturbation_training_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    group_cols = ["MetricGroup", "View", "Direction", "Metric"]
    for keys, sub in df.groupby(group_cols, dropna=False, sort=False):
        mgroup, view, direction, metric = keys
        x = sub[sub["perturbation_trained"] == True]["Scaled_Best"].to_numpy(dtype=float)  # noqa: E712
        y = sub[sub["perturbation_trained"] == False]["Scaled_Best"].to_numpy(dtype=float)  # noqa: E712
        if len(x) < 2 or len(y) < 2:
            continue
        stat, pval = mannwhitneyu(x, y, alternative="two-sided")
        rows.append(
            {
                "MetricGroup": mgroup,
                "View": view,
                "Direction": direction if pd.notna(direction) else "NA",
                "Metric": metric,
                "n_perturb_trained": int(len(x)),
                "n_not_perturb_trained": int(len(y)),
                "mean_scaled_perturb_trained": float(np.mean(x)),
                "mean_scaled_not_perturb_trained": float(np.mean(y)),
                "delta_mean_scaled": float(np.mean(x) - np.mean(y)),
                "p_value_mannwhitney": float(pval),
                "u_statistic": float(stat),
            }
        )
    return pd.DataFrame(rows)


def build_size_trend_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    group_cols = ["MetricGroup", "View", "Direction", "Metric"]
    for keys, sub in df.groupby(group_cols, dropna=False, sort=False):
        mgroup, view, direction, metric = keys
        x = pd.to_numeric(sub["size_code"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub["Scaled_Best"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
        if len(x) < 3 or len(np.unique(x)) < 2:
            continue
        rho, pval = spearmanr(x, y)
        rows.append(
            {
                "MetricGroup": mgroup,
                "View": view,
                "Direction": direction if pd.notna(direction) else "NA",
                "Metric": metric,
                "n_models": int(len(x)),
                "spearman_rho_size_vs_scaled": float(rho),
                "p_value_spearman": float(pval),
            }
        )
    return pd.DataFrame(rows)


def build_size_class_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = (
        df.groupby(["size_class", "perturbation_trained"], as_index=False)
        .agg(
            n_rows=("Metric", "size"),
            n_models=("Track", "nunique"),
            mean_scaled_best=("Scaled_Best", "mean"),
            median_scaled_best=("Scaled_Best", "median"),
            mean_rank_best=("Rank_Best", "mean"),
        )
        .sort_values(["size_class", "perturbation_trained"])
    )
    return out


# ---------------------------------------------------------------------
# Section 3. Main
# ---------------------------------------------------------------------


def run_task3_meta(cfg: Task3MetaConfig) -> None:
    t0 = time.time()
    Path(cfg.analysis_dir).mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(cfg.scoreboard_long_path)
    required = {"Track", "MetricGroup", "View", "Direction", "Metric", "Value", "Scaled_Best", "Rank_Best"}
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"scoreboard_long missing columns: {sorted(missing)}")

    meta = _load_model_meta(cfg.model_meta_map_path)
    merged = raw.merge(meta, on="Track", how="left")
    for col, default in [
        ("model_family", "unknown"),
        ("size_class", "unknown"),
        ("size_code", 0),
        ("perturbation_trained", False),
        ("meta_source", "unmapped"),
    ]:
        merged[col] = merged[col].fillna(default)

    model_scoreboard = build_model_scoreboard(merged)
    perturb_tests = build_perturbation_training_tests(merged)
    size_tests = build_size_trend_tests(merged)
    size_summary = build_size_class_summary(merged)
    score_def = pd.DataFrame(
        [
            {
                "score_definition_id": "equal_weight_metric_slice_mean_scaled_best",
                "score_name": "mean_scaled_best",
                "formula": "mean(Scaled_Best) across all metric slices for each model",
                "weighting": "equal weight per metric-slice row in scoreboard_long",
                "interpretation_note": "Cross-slice aggregate for trend analysis; not a single-task primary endpoint.",
            }
        ]
    )

    out = Path(cfg.analysis_dir)
    meta.to_csv(out / "task3_model_metadata_used.csv", index=False)
    merged.to_csv(out / "task3_model_metric_joined.csv", index=False)
    model_scoreboard.to_csv(out / "task3_model_scoreboard_meta.csv", index=False)
    score_def.to_csv(out / "task3_mean_scaled_best_definition.csv", index=False)
    perturb_tests.to_csv(out / "task3_meta_perturbation_training_tests.csv", index=False)
    size_tests.to_csv(out / "task3_meta_size_trend_tests.csv", index=False)
    size_summary.to_csv(out / "task3_meta_size_class_summary.csv", index=False)

    key_summary = pd.DataFrame(
        [
            {"metric": "n_score_rows", "value": int(len(merged))},
            {"metric": "n_models", "value": int(merged["Track"].nunique())},
            {"metric": "n_metric_groups", "value": int(merged["MetricGroup"].nunique())},
            {"metric": "n_perturbation_training_tests", "value": int(len(perturb_tests))},
            {"metric": "n_size_trend_tests", "value": int(len(size_tests))},
        ]
    )
    key_summary.to_csv(out / "task3_meta_key_summary.csv", index=False)

    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_models": int(merged["Track"].nunique()),
            "n_rows_joined": int(len(merged)),
            "n_perturb_tests": int(len(perturb_tests)),
            "n_size_tests": int(len(size_tests)),
            "score_definition_id": "equal_weight_metric_slice_mean_scaled_best",
        },
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[Task3 Meta] done. Outputs: {cfg.analysis_dir}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Task3 FM meta-analysis")
    parser.add_argument("--scoreboard-long-path", type=str, default=Task3MetaConfig.scoreboard_long_path)
    parser.add_argument("--output-dir", type=str, default=Task3MetaConfig.output_dir)
    parser.add_argument("--model-meta-map-path", type=str, default=Task3MetaConfig.model_meta_map_path)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task3MetaConfig(
        scoreboard_long_path=args.scoreboard_long_path,
        output_dir=args.output_dir,
        model_meta_map_path=args.model_meta_map_path,
    )
    run_task3_meta(cfg)


if __name__ == "__main__":
    main()
