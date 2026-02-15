#!/usr/bin/env python3
"""
Task1 retrieval sensitivity analysis.

Purpose
-------
Evaluate robustness of balanced retrieval conclusions across:
1) balanced gallery size,
2) random seeds (to quantify stochastic variability).

This script reuses the existing Task1 retrieval module and writes
compact, reviewer-facing sensitivity tables with bootstrap CIs.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from _script_bootstrap import setup_project_imports

setup_project_imports()

from m2m_bench.task1.retrieval_instance import Task1RetrievalConfig, run_task1_retrieval


@dataclass
class Task1RetrievalSensitivityConfig:
    m1_candidates_path: str = "./outputs/task1/data/m1_candidates.csv"
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    output_dir: str = "./outputs/task1_reviewer_fixes/retrieval_sensitivity"
    gallery_sizes: str = "64,128,256,512"
    seeds: str = "42,43"
    balanced_true_per_query: int = 1
    balanced_n_repeats: int = 30
    n_perm: int = 30
    tracks: str = "gene,path"
    directions: str = "LINCS->scPerturb,scPerturb->LINCS"
    topk: str = "1,5,10"
    min_queries_per_group: int = 10
    n_bootstrap: int = 1000

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def runs_dir(self) -> str:
        return str(Path(self.output_dir) / "runs")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task1_retrieval_sensitivity.json")


def _parse_int_list(x: str) -> list[int]:
    return [int(v.strip()) for v in str(x).split(",") if v.strip()]


def _bootstrap_ci_mean(
    values: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
    alpha: float = 0.05,
) -> tuple[float, float]:
    v = values[np.isfinite(values)]
    if len(v) == 0:
        return np.nan, np.nan
    if len(v) == 1:
        return float(v[0]), float(v[0])
    idx = rng.integers(0, len(v), size=(n_bootstrap, len(v)))
    means = v[idx].mean(axis=1)
    lo = float(np.quantile(means, alpha / 2))
    hi = float(np.quantile(means, 1 - alpha / 2))
    return lo, hi


def run_task1_retrieval_sensitivity(cfg: Task1RetrievalSensitivityConfig) -> None:
    t0 = time.time()
    Path(cfg.analysis_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.runs_dir).mkdir(parents=True, exist_ok=True)

    gallery_sizes = _parse_int_list(cfg.gallery_sizes)
    seeds = _parse_int_list(cfg.seeds)
    rows: list[dict] = []

    # Include raw baseline from existing Task1 retrieval summary
    base_summary_path = Path("./outputs/task1/retrieval/analysis/retrieval_summary.csv")
    if base_summary_path.exists():
        base = pd.read_csv(base_summary_path)
        base["setting"] = "raw_unbalanced_all_true"
        base["gallery_size_balanced"] = np.nan
        base["seed"] = np.nan
        base["mrr_metric"] = base["mrr"]
        base["top1_metric"] = base["top1"]
        base[
            [
                "track",
                "direction",
                "modality",
                "setting",
                "gallery_size_balanced",
                "seed",
                "n_query",
                "mrr_metric",
                "top1_metric",
            ]
        ].to_csv(Path(cfg.analysis_dir) / "retrieval_sensitivity_raw_baseline.csv", index=False)

    for g in gallery_sizes:
        for seed in seeds:
            run_dir = Path(cfg.runs_dir) / f"g{g}_seed{seed}"
            run_cfg = Task1RetrievalConfig(
                m1_candidates_path=cfg.m1_candidates_path,
                processed_dir=cfg.processed_dir,
                output_dir=str(run_dir),
                tracks=cfg.tracks,
                directions=cfg.directions,
                topk=cfg.topk,
                n_perm=cfg.n_perm,
                min_queries_per_group=cfg.min_queries_per_group,
                random_seed=seed,
                run_balanced_eval=True,
                balanced_gallery_size=g,
                balanced_true_per_query=cfg.balanced_true_per_query,
                balanced_n_repeats=cfg.balanced_n_repeats,
            )
            run_task1_retrieval(run_cfg)

            summary = pd.read_csv(run_dir / "analysis" / "retrieval_summary.csv")
            for _, r in summary.iterrows():
                rows.append(
                    {
                        "track": str(r["track"]),
                        "direction": str(r["direction"]),
                        "modality": str(r["modality"]),
                        "setting": "balanced_true1",
                        "gallery_size_balanced": int(g),
                        "seed": int(seed),
                        "n_query": int(r["n_query"]),
                        "n_query_balanced": int(r.get("n_query_balanced", 0)),
                        "mrr_metric": float(r["mrr_balanced"]),
                        "top1_metric": float(r["top1_balanced"]),
                    }
                )

    long_df = pd.DataFrame(rows)
    long_df.to_csv(Path(cfg.analysis_dir) / "retrieval_sensitivity_long.csv", index=False)

    rng = np.random.default_rng(20260214)
    agg_rows: list[dict] = []
    for keys, grp in long_df.groupby(["track", "direction", "modality", "setting", "gallery_size_balanced"], dropna=False):
        track, direction, modality, setting, g = keys
        mrr_vals = grp["mrr_metric"].to_numpy(dtype=float)
        top1_vals = grp["top1_metric"].to_numpy(dtype=float)
        mrr_lo, mrr_hi = _bootstrap_ci_mean(mrr_vals, n_bootstrap=cfg.n_bootstrap, rng=rng)
        t1_lo, t1_hi = _bootstrap_ci_mean(top1_vals, n_bootstrap=cfg.n_bootstrap, rng=rng)
        agg_rows.append(
            {
                "track": track,
                "direction": direction,
                "modality": modality,
                "setting": setting,
                "gallery_size_balanced": g,
                "n_seed_runs": int(len(grp)),
                "mrr_mean": float(np.nanmean(mrr_vals)),
                "mrr_sd": float(np.nanstd(mrr_vals, ddof=1)) if len(grp) > 1 else 0.0,
                "mrr_ci95_low": mrr_lo,
                "mrr_ci95_high": mrr_hi,
                "top1_mean": float(np.nanmean(top1_vals)),
                "top1_sd": float(np.nanstd(top1_vals, ddof=1)) if len(grp) > 1 else 0.0,
                "top1_ci95_low": t1_lo,
                "top1_ci95_high": t1_hi,
                "n_query_mean": float(np.nanmean(grp["n_query"].to_numpy(dtype=float))),
                "n_query_balanced_mean": float(np.nanmean(grp["n_query_balanced"].to_numpy(dtype=float))),
            }
        )
    agg = pd.DataFrame(agg_rows).sort_values(["modality", "track", "direction", "gallery_size_balanced"])
    agg.to_csv(Path(cfg.analysis_dir) / "retrieval_sensitivity_summary.csv", index=False)

    # Compact modality-level trajectory for manuscript
    mod = (
        agg[agg["setting"] == "balanced_true1"]
        .groupby(["modality", "gallery_size_balanced"], as_index=False)
        .agg(
            mean_mrr=("mrr_mean", "mean"),
            mean_top1=("top1_mean", "mean"),
            mean_n_query_balanced=("n_query_balanced_mean", "mean"),
        )
        .sort_values(["modality", "gallery_size_balanced"])
    )
    mod.to_csv(Path(cfg.analysis_dir) / "retrieval_sensitivity_modality_compact.csv", index=False)

    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_gallery_sizes": int(len(gallery_sizes)),
            "n_seeds": int(len(seeds)),
            "n_long_rows": int(len(long_df)),
            "n_summary_rows": int(len(agg)),
        },
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[Task1 Retrieval Sensitivity] done. Outputs: {cfg.analysis_dir}")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser("Task1 retrieval sensitivity")
    p.add_argument("--m1-candidates-path", type=str, default=Task1RetrievalSensitivityConfig.m1_candidates_path)
    p.add_argument("--processed-dir", type=str, default=Task1RetrievalSensitivityConfig.processed_dir)
    p.add_argument("--output-dir", type=str, default=Task1RetrievalSensitivityConfig.output_dir)
    p.add_argument("--gallery-sizes", type=str, default=Task1RetrievalSensitivityConfig.gallery_sizes)
    p.add_argument("--seeds", type=str, default=Task1RetrievalSensitivityConfig.seeds)
    p.add_argument("--balanced-true-per-query", type=int, default=Task1RetrievalSensitivityConfig.balanced_true_per_query)
    p.add_argument("--balanced-n-repeats", type=int, default=Task1RetrievalSensitivityConfig.balanced_n_repeats)
    p.add_argument("--n-perm", type=int, default=Task1RetrievalSensitivityConfig.n_perm)
    p.add_argument("--tracks", type=str, default=Task1RetrievalSensitivityConfig.tracks)
    p.add_argument("--directions", type=str, default=Task1RetrievalSensitivityConfig.directions)
    p.add_argument("--topk", type=str, default=Task1RetrievalSensitivityConfig.topk)
    p.add_argument("--min-queries-per-group", type=int, default=Task1RetrievalSensitivityConfig.min_queries_per_group)
    p.add_argument("--n-bootstrap", type=int, default=Task1RetrievalSensitivityConfig.n_bootstrap)
    return p


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task1RetrievalSensitivityConfig(
        m1_candidates_path=args.m1_candidates_path,
        processed_dir=args.processed_dir,
        output_dir=args.output_dir,
        gallery_sizes=args.gallery_sizes,
        seeds=args.seeds,
        balanced_true_per_query=args.balanced_true_per_query,
        balanced_n_repeats=args.balanced_n_repeats,
        n_perm=args.n_perm,
        tracks=args.tracks,
        directions=args.directions,
        topk=args.topk,
        min_queries_per_group=args.min_queries_per_group,
        n_bootstrap=args.n_bootstrap,
    )
    run_task1_retrieval_sensitivity(cfg)


if __name__ == "__main__":
    main()
