#!/usr/bin/env python3
"""
Task3 result audit for M2M-Bench.

Purpose
-------
Perform objective, reproducible sanity checks on finished Task3 outputs:
1) file completeness,
2) schema and value-range checks,
3) aggregation consistency (per-query vs summary; per-target vs summary),
4) cross-file key alignment (analysis vs R-ready tables),
5) concise markdown report for quick review.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Section 0. Configuration
# ---------------------------------------------------------------------


@dataclass
class Task3AuditConfig:
    deltas_dir: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Deltas_K562_CellLevel"
    analysis_dir: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Analysis"
    unified_dir: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified"
    output_dir: str = "./outputs/task3_audit"
    tolerance: float = 1e-9

    @property
    def out_analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task3_audit.json")


# ---------------------------------------------------------------------
# Section 1. IO helpers
# ---------------------------------------------------------------------


def _ensure_dir(path_like: str) -> None:
    Path(path_like).mkdir(parents=True, exist_ok=True)


def _read_required_csv(path: Path, name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {name} ({path})")
    return pd.read_csv(path)


def _write_md(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _check_required_columns(df: pd.DataFrame, required: list[str]) -> tuple[bool, str]:
    missing = [c for c in required if c not in df.columns]
    if missing:
        return False, f"missing={missing}"
    return True, "ok"


def _record_check(rows: list[dict], check_name: str, passed: bool, detail: str) -> None:
    rows.append(
        {
            "check_name": check_name,
            "passed": bool(passed),
            "detail": str(detail),
        }
    )


# ---------------------------------------------------------------------
# Section 2. Core audit
# ---------------------------------------------------------------------


def run_task3_audit(cfg: Task3AuditConfig) -> None:
    t0 = time.time()
    _ensure_dir(cfg.out_analysis_dir)
    out = Path(cfg.out_analysis_dir)

    # ------------------------------------------------------------
    # Step A. File inventory
    # ------------------------------------------------------------
    files = [
        Path(cfg.deltas_dir) / "manifest.json",
        Path(cfg.analysis_dir) / "Task3_Retrieval_PerQuery.csv",
        Path(cfg.analysis_dir) / "Task3_Retrieval_Summary.csv",
        Path(cfg.analysis_dir) / "Task3_Pairwise_PerTarget.csv",
        Path(cfg.analysis_dir) / "Task3_Pairwise_Summary.csv",
        Path(cfg.analysis_dir) / "Task3_QC_NaNFiltering.csv",
        Path(cfg.unified_dir) / "annotated_pairwise_per_target.csv",
        Path(cfg.unified_dir) / "annotated_retrieval_per_query.csv",
        Path(cfg.unified_dir) / "viz_fig2_lollipop.csv",
        Path(cfg.unified_dir) / "viz_scoreboard_long.csv",
        Path(cfg.unified_dir) / "viz_scoreboard_wide.csv",
        Path(cfg.unified_dir) / "viz_target_gain.csv",
        Path(cfg.unified_dir) / "viz_stats_deltas_query.csv",
    ]
    file_rows: list[dict] = []
    for path in files:
        exists = path.exists()
        file_rows.append(
            {
                "path": str(path),
                "exists": bool(exists),
                "size_bytes": int(path.stat().st_size) if exists else np.nan,
            }
        )
    file_manifest = pd.DataFrame(file_rows)
    file_manifest.to_csv(out / "task3_file_manifest.csv", index=False)

    # ------------------------------------------------------------
    # Step B. Load core tables
    # ------------------------------------------------------------
    retrieval_per_query = _read_required_csv(Path(cfg.analysis_dir) / "Task3_Retrieval_PerQuery.csv", "Task3_Retrieval_PerQuery.csv")
    retrieval_summary = _read_required_csv(Path(cfg.analysis_dir) / "Task3_Retrieval_Summary.csv", "Task3_Retrieval_Summary.csv")
    pairwise_per_target = _read_required_csv(Path(cfg.analysis_dir) / "Task3_Pairwise_PerTarget.csv", "Task3_Pairwise_PerTarget.csv")
    pairwise_summary = _read_required_csv(Path(cfg.analysis_dir) / "Task3_Pairwise_Summary.csv", "Task3_Pairwise_Summary.csv")
    qc_nan = _read_required_csv(Path(cfg.analysis_dir) / "Task3_QC_NaNFiltering.csv", "Task3_QC_NaNFiltering.csv")
    annotated_pairwise = _read_required_csv(Path(cfg.unified_dir) / "annotated_pairwise_per_target.csv", "annotated_pairwise_per_target.csv")
    annotated_retrieval = _read_required_csv(Path(cfg.unified_dir) / "annotated_retrieval_per_query.csv", "annotated_retrieval_per_query.csv")

    check_rows: list[dict] = []

    # ------------------------------------------------------------
    # Step C. Schema checks
    # ------------------------------------------------------------
    ok, detail = _check_required_columns(
        retrieval_per_query,
        ["Track", "View", "Direction", "Query_ID", "True_Target", "True_Rank", "MRR", "Success_Score"],
    )
    _record_check(check_rows, "retrieval_per_query_required_columns", ok, detail)

    ok, detail = _check_required_columns(
        retrieval_summary,
        ["Track", "View", "Direction", "n_queries", "mean_mrr", "mean_success", "median_rank"],
    )
    _record_check(check_rows, "retrieval_summary_required_columns", ok, detail)

    ok, detail = _check_required_columns(
        pairwise_per_target,
        ["Track", "View", "Target", "centroid_cosine", "edist_mean"],
    )
    _record_check(check_rows, "pairwise_per_target_required_columns", ok, detail)

    ok, detail = _check_required_columns(
        pairwise_summary,
        ["Track", "View", "n_targets", "mean_centroid_cos", "mean_edist"],
    )
    _record_check(check_rows, "pairwise_summary_required_columns", ok, detail)

    # ------------------------------------------------------------
    # Step D. Value-range checks
    # ------------------------------------------------------------
    mrr = pd.to_numeric(retrieval_per_query["MRR"], errors="coerce")
    success = pd.to_numeric(retrieval_per_query["Success_Score"], errors="coerce")
    true_rank = pd.to_numeric(retrieval_per_query["True_Rank"], errors="coerce")
    _record_check(
        check_rows,
        "retrieval_MRR_in_0_1",
        bool(np.all((mrr >= 0) & (mrr <= 1))),
        f"min={float(np.nanmin(mrr))}, max={float(np.nanmax(mrr))}",
    )
    _record_check(
        check_rows,
        "retrieval_SuccessScore_in_0_1",
        bool(np.all((success >= 0) & (success <= 1))),
        f"min={float(np.nanmin(success))}, max={float(np.nanmax(success))}",
    )
    _record_check(
        check_rows,
        "retrieval_TrueRank_positive",
        bool(np.nanmin(true_rank) >= 1),
        f"min={float(np.nanmin(true_rank))}",
    )

    # ------------------------------------------------------------
    # Step E. Aggregation consistency checks
    # ------------------------------------------------------------
    agg_retrieval = (
        retrieval_per_query.groupby(["Track", "View", "Direction"], as_index=False)
        .agg(
            n_queries=("Query_ID", "count"),
            mean_mrr=("MRR", "mean"),
            mean_success=("Success_Score", "mean"),
            median_rank=("True_Rank", "median"),
        )
    )
    merged_retrieval = agg_retrieval.merge(
        retrieval_summary,
        on=["Track", "View", "Direction"],
        how="outer",
        suffixes=("_calc", "_csv"),
        indicator=True,
    )
    _record_check(
        check_rows,
        "retrieval_summary_group_key_match",
        bool((merged_retrieval["_merge"] == "both").all()),
        merged_retrieval["_merge"].value_counts(dropna=False).to_dict(),
    )
    if (merged_retrieval["_merge"] == "both").any():
        both = merged_retrieval[merged_retrieval["_merge"] == "both"].copy()
        diff_n = (both["n_queries_calc"] - both["n_queries_csv"]).abs().max()
        diff_mrr = (both["mean_mrr_calc"] - both["mean_mrr_csv"]).abs().max()
        diff_success = (both["mean_success_calc"] - both["mean_success_csv"]).abs().max()
        _record_check(check_rows, "retrieval_summary_n_queries_consistent", bool(diff_n <= cfg.tolerance), f"max_abs_diff={float(diff_n)}")
        _record_check(check_rows, "retrieval_summary_mean_mrr_consistent", bool(diff_mrr <= cfg.tolerance), f"max_abs_diff={float(diff_mrr)}")
        _record_check(
            check_rows,
            "retrieval_summary_mean_success_consistent",
            bool(diff_success <= cfg.tolerance),
            f"max_abs_diff={float(diff_success)}",
        )

    agg_pairwise = (
        pairwise_per_target.groupby(["Track", "View"], as_index=False)
        .agg(
            n_targets=("Target", "nunique"),
            mean_centroid_cos=("centroid_cosine", "mean"),
            mean_edist=("edist_mean", "mean"),
        )
    )
    merged_pairwise = agg_pairwise.merge(
        pairwise_summary,
        on=["Track", "View"],
        how="outer",
        suffixes=("_calc", "_csv"),
        indicator=True,
    )
    _record_check(
        check_rows,
        "pairwise_summary_group_key_match",
        bool((merged_pairwise["_merge"] == "both").all()),
        merged_pairwise["_merge"].value_counts(dropna=False).to_dict(),
    )
    if (merged_pairwise["_merge"] == "both").any():
        both = merged_pairwise[merged_pairwise["_merge"] == "both"].copy()
        diff_n = (both["n_targets_calc"] - both["n_targets_csv"]).abs().max()
        diff_cos = (both["mean_centroid_cos_calc"] - both["mean_centroid_cos_csv"]).abs().max()
        diff_ed = (both["mean_edist_calc"] - both["mean_edist_csv"]).abs().max()
        _record_check(check_rows, "pairwise_summary_n_targets_consistent", bool(diff_n <= cfg.tolerance), f"max_abs_diff={float(diff_n)}")
        _record_check(
            check_rows,
            "pairwise_summary_mean_centroid_cos_consistent",
            bool(diff_cos <= cfg.tolerance),
            f"max_abs_diff={float(diff_cos)}",
        )
        _record_check(
            check_rows,
            "pairwise_summary_mean_edist_consistent",
            bool(diff_ed <= cfg.tolerance),
            f"max_abs_diff={float(diff_ed)}",
        )

    # ------------------------------------------------------------
    # Step F. Key uniqueness and cross-file alignment
    # ------------------------------------------------------------
    dup_count = int(retrieval_per_query.duplicated(["Track", "View", "Direction", "Query_ID"]).sum())
    _record_check(check_rows, "retrieval_query_key_unique", bool(dup_count == 0), f"duplicate_rows={dup_count}")

    key_cols = ["Track", "View", "Direction", "Query_ID"]
    a_keys = set(tuple(v) for v in retrieval_per_query[key_cols].itertuples(index=False, name=None))
    b_keys = set(tuple(v) for v in annotated_retrieval[key_cols].itertuples(index=False, name=None))
    _record_check(
        check_rows,
        "analysis_retrieval_key_subset_of_annotated",
        bool(a_keys.issubset(b_keys)),
        f"missing_in_annotated={len(a_keys - b_keys)}, extra_in_annotated={len(b_keys - a_keys)}",
    )
    _record_check(
        check_rows,
        "annotated_retrieval_row_count_match",
        bool(len(annotated_retrieval) == len(retrieval_per_query)),
        f"analysis_rows={len(retrieval_per_query)}, annotated_rows={len(annotated_retrieval)}",
    )

    # ------------------------------------------------------------
    # Step G. Structured summaries for interpretation
    # ------------------------------------------------------------
    key_summary_rows = [
        {"metric": "n_retrieval_per_query_rows", "value": int(len(retrieval_per_query))},
        {"metric": "n_retrieval_summary_rows", "value": int(len(retrieval_summary))},
        {"metric": "n_pairwise_per_target_rows", "value": int(len(pairwise_per_target))},
        {"metric": "n_pairwise_summary_rows", "value": int(len(pairwise_summary))},
        {"metric": "n_tracks_retrieval", "value": int(retrieval_per_query["Track"].nunique())},
        {"metric": "n_views_retrieval", "value": int(retrieval_per_query["View"].nunique())},
        {"metric": "n_directions_retrieval", "value": int(retrieval_per_query["Direction"].nunique())},
        {"metric": "n_tracks_pairwise", "value": int(pairwise_per_target["Track"].nunique())},
        {"metric": "n_views_pairwise", "value": int(pairwise_per_target["View"].nunique())},
        {"metric": "retrieval_mrr_min", "value": float(mrr.min())},
        {"metric": "retrieval_mrr_max", "value": float(mrr.max())},
        {"metric": "retrieval_success_min", "value": float(success.min())},
        {"metric": "retrieval_success_max", "value": float(success.max())},
    ]
    key_summary = pd.DataFrame(key_summary_rows)
    key_summary.to_csv(out / "task3_key_summary.csv", index=False)

    track_view_dir = (
        retrieval_per_query.groupby(["Track", "View", "Direction"], as_index=False)
        .agg(
            n_queries=("Query_ID", "count"),
            mean_mrr=("MRR", "mean"),
            mean_success=("Success_Score", "mean"),
            median_rank=("True_Rank", "median"),
        )
        .sort_values(["Track", "View", "Direction"])
    )
    track_view_dir.to_csv(out / "task3_track_view_direction_summary.csv", index=False)

    nan_summary = (
        qc_nan.groupby(["phase", "track", "view", "role"], dropna=False, as_index=False)
        .agg(
            n_rows=("n_total", "sum"),
            n_kept=("n_kept", "sum"),
            n_nonfinite_rows=("n_nonfinite_rows", "sum"),
            n_empty_target=("n_empty_target", "sum"),
        )
        .sort_values(["phase", "track", "view", "role"])
    )
    nan_summary.to_csv(out / "task3_nan_filter_summary.csv", index=False)

    checks = pd.DataFrame(check_rows)
    checks.to_csv(out / "task3_consistency_checks.csv", index=False)

    # ------------------------------------------------------------
    # Step H. Markdown report + run manifest
    # ------------------------------------------------------------
    n_pass = int(checks["passed"].sum())
    n_total = int(len(checks))
    n_fail = int(n_total - n_pass)

    report_lines = [
        "# Task3 Audit Report",
        "",
        "## Headline",
        f"- Checks passed: {n_pass}/{n_total}",
        f"- Checks failed: {n_fail}",
        f"- Retrieval rows: {len(retrieval_per_query):,}",
        f"- Pairwise rows: {len(pairwise_per_target):,}",
        "",
        "## Files",
        f"- File manifest: `{out / 'task3_file_manifest.csv'}`",
        "",
        "## Core checks",
    ]
    for _, row in checks.iterrows():
        icon = "PASS" if bool(row["passed"]) else "FAIL"
        report_lines.append(f"- [{icon}] {row['check_name']}: {row['detail']}")
    report_lines.extend(
        [
            "",
            "## Output tables",
            f"- `{out / 'task3_consistency_checks.csv'}`",
            f"- `{out / 'task3_key_summary.csv'}`",
            f"- `{out / 'task3_track_view_direction_summary.csv'}`",
            f"- `{out / 'task3_nan_filter_summary.csv'}`",
        ]
    )
    _write_md(out / "task3_audit_report.md", report_lines)

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_checks": n_total,
            "n_pass": n_pass,
            "n_fail": n_fail,
            "n_retrieval_rows": int(len(retrieval_per_query)),
            "n_pairwise_rows": int(len(pairwise_per_target)),
        },
        "outputs": {
            "analysis_dir": str(out),
        },
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------
# Section 3. CLI
# ---------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task3 results audit")
    parser.add_argument("--deltas-dir", type=str, default=Task3AuditConfig.deltas_dir)
    parser.add_argument("--analysis-dir", type=str, default=Task3AuditConfig.analysis_dir)
    parser.add_argument("--unified-dir", type=str, default=Task3AuditConfig.unified_dir)
    parser.add_argument("--output-dir", type=str, default=Task3AuditConfig.output_dir)
    parser.add_argument("--tolerance", type=float, default=Task3AuditConfig.tolerance)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task3AuditConfig(
        deltas_dir=args.deltas_dir,
        analysis_dir=args.analysis_dir,
        unified_dir=args.unified_dir,
        output_dir=args.output_dir,
        tolerance=args.tolerance,
    )
    run_task3_audit(cfg)


if __name__ == "__main__":
    main()
