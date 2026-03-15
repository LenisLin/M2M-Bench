#!/usr/bin/env python3
"""
Phase 1 manuscript-support skeleton for A1: Task1 internal vs Task1 cross.

Analysis purpose:
- Compare Task1 internal versus Task1 cross performance on the lawful
  common-scope slice.

Why this script exists in the manuscript workflow:
- The manuscript needs a downstream comparison scaffold that distinguishes
  between preferred exact pairing and explicit summary fallback without
  overstating what the frozen local outputs support.

Lawful input scope:
- reviewed local S1 `task1_leaderboard_long.csv`
- reviewed local S1 `task1_retrieval_summary.csv`
- reviewed local S1 `task1_retrieval_per_query.parquet` only when present and
  schema-compatible
- reviewed local S2 `task1_cross_leaderboard_long.csv`
- reviewed local S2 `task1_cross_retrieval_summary.csv`
- reviewed local S2 `task1_cross_retrieval_per_query.parquet`
- reviewed local S2 `task1_group_cross.parquet`

Prohibited behavior:
- do not invent pseudo-pairs
- do not label summary fallback as if it were exact pairing
- do not re-export or mutate upstream benchmark artifacts
- do not treat absent local exact-pair support as permission to guess

Intended output semantics:
- one row per lawful comparison unit after implementation review
- rows must carry `comparison_mode` that distinguishes exact pairing from
  summary fallback
- rows must carry explicit status and notes when exact pairing is unavailable

What remains TODO in later implementation:
- review whether exact lawful pairing can be supported from current frozen
  per-query outputs alone
- if not, implement the highest-validity summary fallback with explicit labels
- freeze the final row granularity for fallback comparisons
- confirm whether any last-resort upstream internal-group export is needed;
  this scaffold does not assume it is allowed
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
NAS_RUNS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs")
MANUSCRIPT_PHASE1_ROOT = NAS_RUNS_ROOT / "manuscript_phase1" / "a1_task1_internal_cross"

DEFAULT_INTERNAL_LEADERBOARD_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv"
)
DEFAULT_INTERNAL_SUMMARY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_summary.csv"
)
DEFAULT_INTERNAL_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_per_query.parquet"
)
DEFAULT_CROSS_LEADERBOARD_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_leaderboard_long.csv"
)
DEFAULT_CROSS_SUMMARY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_summary.csv"
)
DEFAULT_CROSS_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_per_query.parquet"
)
DEFAULT_GROUP_CROSS_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_group_cross.parquet"
)
DEFAULT_OUTPUT_PATH = MANUSCRIPT_PHASE1_ROOT / "task1_internal_vs_cross_common_scope_comparison.csv"

COMMON_SCOPE_REPRESENTATIONS = {"Gene", "Pathway"}
COMMON_SCOPE_PERTURBATION_TYPE = "Genetic"

LEADERBOARD_REQUIRED_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "metric_name",
    "metric_value",
    "n_total",
    "n_valid",
    "n_excluded",
    "N_gallery_max",
    "cross_alignment_contract",
]
SUMMARY_REQUIRED_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "n_total",
    "n_valid",
    "n_excluded_missing_metric_or_mpos0",
    "N_gallery_max",
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
]
PER_QUERY_REQUIRED_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "query_uid",
    "cell_line",
    "target_token",
    "N_gallery",
    "N_gallery_max",
    "m_pos",
    "mrr_corrected",
    "hit1_corrected",
    "hit5_corrected",
    "hit10_corrected",
    "loo_policy",
    "cross_alignment_contract",
]
GROUP_CROSS_REQUIRED_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "group_id",
    "cell_line",
    "target_token",
    "n_L",
    "n_S",
    "n_L_sub",
    "n_S_sub",
    "cosine",
    "pcc",
    "edist",
    "cross_alignment_contract",
]
OUTPUT_COLUMNS = [
    "representation",
    "metric_name",
    "internal_dataset_or_direction",
    "cross_dataset_or_direction",
    "perturbation_type",
    "comparison_mode",
    "internal_value",
    "cross_value",
    "delta_value",
    "internal_n_total",
    "internal_n_valid",
    "cross_n_total",
    "cross_n_valid",
    "status",
    "pairing_notes",
]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the A1 scaffold."""
    parser = argparse.ArgumentParser(description="Phase 1 skeleton for A1 Task1 internal vs cross.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--internal-leaderboard-path", type=Path, default=DEFAULT_INTERNAL_LEADERBOARD_PATH)
    parser.add_argument("--internal-summary-path", type=Path, default=DEFAULT_INTERNAL_SUMMARY_PATH)
    parser.add_argument("--internal-per-query-path", type=Path, default=DEFAULT_INTERNAL_PER_QUERY_PATH)
    parser.add_argument("--cross-leaderboard-path", type=Path, default=DEFAULT_CROSS_LEADERBOARD_PATH)
    parser.add_argument("--cross-summary-path", type=Path, default=DEFAULT_CROSS_SUMMARY_PATH)
    parser.add_argument("--cross-per-query-path", type=Path, default=DEFAULT_CROSS_PER_QUERY_PATH)
    parser.add_argument("--group-cross-path", type=Path, default=DEFAULT_GROUP_CROSS_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    return parser.parse_args()


def resolve_path(project_root: Path, raw_path: Path) -> Path:
    """Resolve a CLI path against the requested project root."""
    return raw_path if raw_path.is_absolute() else (project_root / raw_path)


def load_table(path: Path, label: str) -> pd.DataFrame:
    """Load a CSV or parquet table without hidden fallback behavior."""
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if path.suffix == ".csv":
        return pd.read_csv(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"{label} must be CSV or parquet, got: {path}")


def maybe_load_table(path: Path, label: str) -> pd.DataFrame | None:
    """Load an optional table when present; return None when absent."""
    if not path.exists():
        return None
    return load_table(path, label)


def ensure_required_columns(frame: pd.DataFrame, required: Iterable[str], label: str) -> None:
    """Fail fast when a required frozen column is absent."""
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def ensure_non_null_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    """Reject nulls in core scope columns rather than filling them silently."""
    if frame[columns].isnull().any().any():
        raise ValueError(f"{label} contains null values in required columns: {columns}")


def filter_to_common_scope(frame: pd.DataFrame) -> pd.DataFrame:
    """
    Apply only the coarse common-scope labels.

    This is not the final lawful common-scope anchor. The later implementation
    must retain only the reviewed cross-supported slice rather than relying on
    broad labels alone.
    """
    return frame.loc[
        frame["representation"].isin(COMMON_SCOPE_REPRESENTATIONS)
        & frame["perturbation_type"].eq(COMMON_SCOPE_PERTURBATION_TYPE)
    ].copy()


def load_reviewed_inputs(args: argparse.Namespace, project_root: Path) -> dict[str, pd.DataFrame | None]:
    """Load the reviewed local A1 inputs."""
    return {
        "internal_leaderboard": load_table(
            resolve_path(project_root, args.internal_leaderboard_path),
            "S1 internal leaderboard",
        ),
        "internal_summary": load_table(
            resolve_path(project_root, args.internal_summary_path),
            "S1 internal retrieval summary",
        ),
        "internal_per_query": maybe_load_table(
            resolve_path(project_root, args.internal_per_query_path),
            "S1 internal retrieval per-query",
        ),
        "cross_leaderboard": load_table(
            resolve_path(project_root, args.cross_leaderboard_path),
            "S2 cross leaderboard",
        ),
        "cross_summary": load_table(
            resolve_path(project_root, args.cross_summary_path),
            "S2 cross retrieval summary",
        ),
        "cross_per_query": load_table(
            resolve_path(project_root, args.cross_per_query_path),
            "S2 cross retrieval per-query",
        ),
        "group_cross": load_table(resolve_path(project_root, args.group_cross_path), "S2 group cross table"),
    }


def validate_inputs(inputs: dict[str, pd.DataFrame | None]) -> None:
    """Run the reviewed schema checks that are safe to freeze now."""
    ensure_required_columns(inputs["internal_leaderboard"], LEADERBOARD_REQUIRED_COLUMNS, "S1 internal leaderboard")
    ensure_required_columns(inputs["internal_summary"], SUMMARY_REQUIRED_COLUMNS, "S1 internal retrieval summary")
    ensure_required_columns(inputs["cross_leaderboard"], LEADERBOARD_REQUIRED_COLUMNS, "S2 cross leaderboard")
    ensure_required_columns(inputs["cross_summary"], SUMMARY_REQUIRED_COLUMNS, "S2 cross retrieval summary")
    ensure_required_columns(inputs["cross_per_query"], PER_QUERY_REQUIRED_COLUMNS, "S2 cross retrieval per-query")
    ensure_required_columns(inputs["group_cross"], GROUP_CROSS_REQUIRED_COLUMNS, "S2 group cross table")

    if inputs["internal_per_query"] is not None:
        ensure_required_columns(inputs["internal_per_query"], PER_QUERY_REQUIRED_COLUMNS, "S1 internal retrieval per-query")

    ensure_non_null_columns(
        inputs["internal_summary"],
        ["scope", "dataset_or_direction", "perturbation_type", "representation"],
        "S1 internal retrieval summary",
    )
    ensure_non_null_columns(
        inputs["cross_summary"],
        ["scope", "dataset_or_direction", "perturbation_type", "representation"],
        "S2 cross retrieval summary",
    )


def build_internal_cross_comparison(inputs: dict[str, pd.DataFrame | None]) -> pd.DataFrame:
    """
    Build the future A1 comparison table.

    The scaffold applies only a coarse common-scope prefilter and does not yet
    decide whether exact pairing is supported from the frozen local outputs.
    """
    internal_summary = filter_to_common_scope(inputs["internal_summary"])
    cross_summary = filter_to_common_scope(inputs["cross_summary"])
    internal_per_query = None
    if inputs["internal_per_query"] is not None:
        internal_per_query = filter_to_common_scope(inputs["internal_per_query"])
    cross_per_query = filter_to_common_scope(inputs["cross_per_query"])
    group_cross = filter_to_common_scope(inputs["group_cross"])

    if internal_summary.empty:
        raise ValueError("S1 internal summary has no common-scope rows after lawful filtering")
    if cross_summary.empty:
        raise ValueError("S2 cross summary has no common-scope rows after lawful filtering")
    if group_cross.empty:
        raise ValueError("S2 group cross table has no common-scope rows after lawful filtering")

    _ = internal_per_query
    _ = cross_per_query

    # TODO: Anchor the final lawful common-scope slice to the reviewed
    # cross-supported slice. The broad `Gene`/`Pathway` plus `Genetic` filter is
    # only a coarse prefilter and is not sufficient by itself.
    #
    # TODO: Review whether current frozen per-query outputs support lawful exact
    # pairing on a shared query identity. The planning docs prefer exact pairing,
    # but this scaffold must not assume exact pairing merely because an optional
    # internal per-query table happens to be present.
    #
    # TODO: If exact pairing is lawful, build rows labeled
    # `comparison_mode=exact_pairing`.
    #
    # TODO: If exact pairing is unavailable because the optional internal
    # per-query input is absent, schema-incompatible, or otherwise insufficient,
    # degrade only to an explicit summary fallback and label rows as
    # `comparison_mode=summary_fallback`.
    #
    # TODO: Do not infer exactness from summary coincidences, and do not invent
    # pseudo-pairs from unmatched denominators.
    return pd.DataFrame(columns=OUTPUT_COLUMNS)


def write_output(frame: pd.DataFrame, output_path: Path) -> None:
    """Write the manuscript-support comparison table."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)


def main() -> int:
    """Run the A1 skeleton through loading, validation, and placeholder output."""
    args = parse_args()
    project_root = resolve_path(Path.cwd(), args.project_root)
    inputs = load_reviewed_inputs(args, project_root)
    validate_inputs(inputs)
    output = build_internal_cross_comparison(inputs)
    write_output(output, resolve_path(project_root, args.output_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
