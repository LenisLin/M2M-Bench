#!/usr/bin/env python3
"""
Phase 1 manuscript-support skeleton for A3: Task2 direction robustness.

Analysis purpose:
- Audit direction-specific Task2 retrieval support without collapsing `C2G`
  and `G2C`.

Why this script exists in the manuscript workflow:
- The manuscript needs a downstream robustness table that keeps corrected
  retrieval metrics primary while making gallery and `m_pos` context explicit.
- This script is manuscript support only. It must consume frozen benchmark
  outputs rather than recompute or mutate upstream S5/S6 artifacts.

Lawful input scope:
- reviewed local outputs from S5 corrected multisource retrieval:
  `task2_retrieval_summary.csv`
- reviewed local outputs from S5 corrected multisource retrieval:
  `task2_retrieval_per_query.parquet`
- reviewed local outputs from S5 corrected multisource retrieval:
  `task2_chance_identity_check.csv`
- reviewed local outputs from S6 corrected multisource synthesis:
  `task2_retrieval_leaderboard.csv`

Prohibited behavior:
- do not collapse `C2G` and `G2C`
- do not replace corrected metrics with raw metrics
- do not silently pool across datasets or cell lines
- do not silently redefine S5/S6 contracts
- do not write back into benchmark-stage run directories

Intended output semantics:
- one row per lawful robustness slice after implementation review
- corrected retrieval fields remain primary
- gallery-definition and `m_pos` context remain explicit
- output rows carry status and provenance notes when a slice is partial,
  unavailable, or blocked

What remains TODO in later implementation:
- freeze the final reviewed robustness-slice labels
- build denominator-conserving robustness slices from S5 per-query evidence
- carry chance-check coverage onto the exact audited slice keys
- decide which partial-output rows are lawful for still-unreviewed slice logic
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
NAS_RUNS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs")
MANUSCRIPT_PHASE1_ROOT = NAS_RUNS_ROOT / "manuscript_phase1" / "a3_direction_robustness"

DEFAULT_SUMMARY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv"
)
DEFAULT_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet"
)
DEFAULT_CHANCE_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_chance_identity_check.csv"
)
DEFAULT_LEADERBOARD_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv"
)
DEFAULT_OUTPUT_PATH = MANUSCRIPT_PHASE1_ROOT / "task2_direction_robustness_audit.csv"

FLOAT_COMPARE_TOL = 1e-12
REVIEWED_SLICE_LABEL = "reviewed_primary_direction_slice"

SLICE_KEY_COLUMNS = ["dataset", "cell_line", "direction", "representation"]
SUMMARY_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "gallery_definition_id",
    "pos_definition_id",
    "n_total",
    "n_valid",
    "n_excluded_missing_metric_or_mpos0",
    "N_gallery_mean",
    "N_gallery_max",
    "m_pos_mean",
    "m_pos_p50",
    "m_pos_p90",
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
]
PER_QUERY_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "query_row_id",
    "query_uid",
    "gallery_definition_id",
    "pos_definition_id",
    "gallery_ids_hash",
    "N_gallery",
    "m_pos",
    "mrr_corrected",
    "hit1_corrected",
    "hit5_corrected",
    "hit10_corrected",
]
CHANCE_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "abs_delta_mrr",
    "abs_delta_hit1",
    "abs_delta_hit5",
    "abs_delta_hit10",
]
LEADERBOARD_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "gallery_definition_id",
    "pos_definition_id",
    "n_total",
    "n_valid",
    "N_gallery_mean",
    "N_gallery_max",
    "m_pos_mean",
    "m_pos_p50",
    "m_pos_p90",
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
]
OUTPUT_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "robustness_slice_label",
    "gallery_definition_id",
    "pos_definition_id",
    "rank_by_mean_mrr_corrected",
    "n_total",
    "n_valid",
    "n_valid_observed",
    "n_excluded_missing_metric_or_mpos0",
    "N_gallery_mean",
    "N_gallery_min_observed",
    "N_gallery_max",
    "m_pos_mean",
    "m_pos_p50",
    "m_pos_p90",
    "m_pos_min_observed",
    "m_pos_max_observed",
    "gallery_ids_hash_unique_n",
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
    "chance_check_available_bool",
    "chance_check_max_abs_delta",
    "status",
    "status_reason",
    "provenance_notes",
]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the manuscript-support scaffold."""
    parser = argparse.ArgumentParser(description="Phase 1 skeleton for A3 Task2 direction robustness.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--summary-path", type=Path, default=DEFAULT_SUMMARY_PATH)
    parser.add_argument("--per-query-path", type=Path, default=DEFAULT_PER_QUERY_PATH)
    parser.add_argument("--chance-path", type=Path, default=DEFAULT_CHANCE_PATH)
    parser.add_argument("--leaderboard-path", type=Path, default=DEFAULT_LEADERBOARD_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    return parser.parse_args()


def resolve_path(project_root: Path, raw_path: Path) -> Path:
    """Resolve a CLI path against the requested project root."""
    return raw_path if raw_path.is_absolute() else (project_root / raw_path)


def load_table(path: Path, label: str) -> pd.DataFrame:
    """Load a CSV or parquet table without coercing schema drift."""
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if path.suffix == ".csv":
        return pd.read_csv(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"{label} must be CSV or parquet, got: {path}")


def ensure_required_columns(frame: pd.DataFrame, required: Iterable[str], label: str) -> None:
    """Fail fast when a required frozen column is absent."""
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def ensure_non_null_keys(frame: pd.DataFrame, key_columns: list[str], label: str) -> None:
    """Reject null key material instead of letting joins drop rows silently."""
    if frame[key_columns].isnull().any().any():
        raise ValueError(f"{label} contains null values in key columns: {key_columns}")


def ensure_unique_keys(frame: pd.DataFrame, key_columns: list[str], label: str) -> None:
    """Check frozen key uniqueness on summary-like tables."""
    duplicated = frame.duplicated(subset=key_columns, keep=False)
    if duplicated.any():
        raise ValueError(f"{label} has duplicate key rows on {key_columns}")


def compare_key_sets(left: pd.DataFrame, right: pd.DataFrame, key_columns: list[str], label: str) -> None:
    """Check key-set consistency across frozen upstream outputs."""
    left_keys = set(map(tuple, left[key_columns].itertuples(index=False, name=None)))
    right_keys = set(map(tuple, right[key_columns].itertuples(index=False, name=None)))
    if left_keys != right_keys:
        missing_from_right = sorted(left_keys - right_keys)[:5]
        missing_from_left = sorted(right_keys - left_keys)[:5]
        raise ValueError(
            f"{label} key mismatch on {key_columns}; "
            f"missing_from_right={missing_from_right}; missing_from_left={missing_from_left}"
        )


def validate_summary_leaderboard_consistency(summary: pd.DataFrame, leaderboard: pd.DataFrame) -> None:
    """Validate the shared S5/S6 retrieval slice contract on exact keys."""
    compare_key_sets(summary, leaderboard, SLICE_KEY_COLUMNS, "S5 summary vs S6 leaderboard")

    merged = summary.merge(
        leaderboard,
        on=SLICE_KEY_COLUMNS,
        how="inner",
        validate="one_to_one",
        suffixes=("_summary", "_leaderboard"),
    )
    exact_columns = [
        "gallery_definition_id",
        "pos_definition_id",
        "n_total",
        "n_valid",
        "N_gallery_max",
    ]
    float_columns = [
        "N_gallery_mean",
        "m_pos_mean",
        "m_pos_p50",
        "m_pos_p90",
        "mean_mrr_corrected",
        "mean_hit1_corrected",
        "mean_hit5_corrected",
        "mean_hit10_corrected",
    ]

    for column in exact_columns:
        summary_col = f"{column}_summary"
        leaderboard_col = f"{column}_leaderboard"
        if not merged[summary_col].equals(merged[leaderboard_col]):
            raise ValueError(f"S5 summary vs S6 leaderboard mismatch for column {column}")
    for column in float_columns:
        assert_close_series(
            merged[f"{column}_summary"],
            merged[f"{column}_leaderboard"],
            f"S5 summary vs S6 leaderboard {column}",
        )


def validate_summary_internal_consistency(summary: pd.DataFrame) -> None:
    """Validate denominator conservation inside the frozen S5 summary table."""
    expected_total = summary["n_valid"] + summary["n_excluded_missing_metric_or_mpos0"]
    if not expected_total.equals(summary["n_total"]):
        raise ValueError("S5 retrieval summary violates n_total == n_valid + n_excluded_missing_metric_or_mpos0")


def assert_close_series(
    left: pd.Series,
    right: pd.Series,
    label: str,
    tol: float = FLOAT_COMPARE_TOL,
) -> None:
    """Compare aligned numeric series using a strict absolute tolerance."""
    delta = (left - right).abs()
    if not (delta <= tol).all():
        raise ValueError(f"{label} mismatch exceeds tolerance {tol}: max_abs_delta={float(delta.max())}")


def build_per_query_observed_context(per_query: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-query evidence back to the reviewed slice keys."""
    grouped = per_query.groupby(SLICE_KEY_COLUMNS, dropna=False)
    observed = grouped.agg(
        gallery_definition_id_observed=("gallery_definition_id", "first"),
        pos_definition_id_observed=("pos_definition_id", "first"),
        gallery_definition_id_nunique=("gallery_definition_id", "nunique"),
        pos_definition_id_nunique=("pos_definition_id", "nunique"),
        n_valid_observed=("query_uid", "size"),
        N_gallery_mean_observed=("N_gallery", "mean"),
        N_gallery_min_observed=("N_gallery", "min"),
        N_gallery_max_observed=("N_gallery", "max"),
        m_pos_mean_observed=("m_pos", "mean"),
        m_pos_min_observed=("m_pos", "min"),
        m_pos_max_observed=("m_pos", "max"),
        gallery_ids_hash_unique_n=("gallery_ids_hash", "nunique"),
        mean_mrr_corrected_observed=("mrr_corrected", "mean"),
        mean_hit1_corrected_observed=("hit1_corrected", "mean"),
        mean_hit5_corrected_observed=("hit5_corrected", "mean"),
        mean_hit10_corrected_observed=("hit10_corrected", "mean"),
    ).reset_index()

    observed["m_pos_p50_observed"] = grouped["m_pos"].quantile(0.5).to_numpy()
    observed["m_pos_p90_observed"] = grouped["m_pos"].quantile(0.9).to_numpy()

    if (observed["gallery_definition_id_nunique"] != 1).any():
        raise ValueError("Per-query gallery_definition_id is not single-valued within a reviewed slice")
    if (observed["pos_definition_id_nunique"] != 1).any():
        raise ValueError("Per-query pos_definition_id is not single-valued within a reviewed slice")

    return observed


def validate_per_query_against_summary(summary: pd.DataFrame, observed: pd.DataFrame) -> pd.DataFrame:
    """Validate that per-query evidence reproduces the frozen S5 slice values."""
    compare_key_sets(summary, observed, SLICE_KEY_COLUMNS, "S5 summary vs S5 per-query observed context")

    merged = summary.merge(
        observed,
        on=SLICE_KEY_COLUMNS,
        how="inner",
        validate="one_to_one",
    )

    if not merged["gallery_definition_id"].equals(merged["gallery_definition_id_observed"]):
        raise ValueError("S5 summary vs per-query mismatch for gallery_definition_id")
    if not merged["pos_definition_id"].equals(merged["pos_definition_id_observed"]):
        raise ValueError("S5 summary vs per-query mismatch for pos_definition_id")
    if not merged["n_valid"].equals(merged["n_valid_observed"]):
        raise ValueError("S5 summary vs per-query mismatch for n_valid")
    if not merged["N_gallery_max"].equals(merged["N_gallery_max_observed"]):
        raise ValueError("S5 summary vs per-query mismatch for N_gallery_max")

    assert_close_series(merged["N_gallery_mean"], merged["N_gallery_mean_observed"], "N_gallery_mean")
    assert_close_series(merged["m_pos_mean"], merged["m_pos_mean_observed"], "m_pos_mean")
    assert_close_series(merged["m_pos_p50"], merged["m_pos_p50_observed"], "m_pos_p50")
    assert_close_series(merged["m_pos_p90"], merged["m_pos_p90_observed"], "m_pos_p90")
    assert_close_series(merged["mean_mrr_corrected"], merged["mean_mrr_corrected_observed"], "mean_mrr_corrected")
    assert_close_series(merged["mean_hit1_corrected"], merged["mean_hit1_corrected_observed"], "mean_hit1_corrected")
    assert_close_series(merged["mean_hit5_corrected"], merged["mean_hit5_corrected_observed"], "mean_hit5_corrected")
    assert_close_series(merged["mean_hit10_corrected"], merged["mean_hit10_corrected_observed"], "mean_hit10_corrected")

    return merged


def load_reviewed_inputs(args: argparse.Namespace, project_root: Path) -> dict[str, pd.DataFrame]:
    """Load the reviewed frozen inputs used by the A3 scaffold."""
    return {
        "summary": load_table(resolve_path(project_root, args.summary_path), "S5 retrieval summary"),
        "per_query": load_table(resolve_path(project_root, args.per_query_path), "S5 retrieval per-query"),
        "chance": load_table(resolve_path(project_root, args.chance_path), "S5 chance identity check"),
        "leaderboard": load_table(resolve_path(project_root, args.leaderboard_path), "S6 retrieval leaderboard"),
    }


def validate_inputs(inputs: dict[str, pd.DataFrame]) -> None:
    """Run the frozen schema and key checks that are already reviewed."""
    ensure_required_columns(inputs["summary"], SUMMARY_REQUIRED_COLUMNS, "S5 retrieval summary")
    ensure_required_columns(inputs["per_query"], PER_QUERY_REQUIRED_COLUMNS, "S5 retrieval per-query")
    ensure_required_columns(inputs["chance"], CHANCE_REQUIRED_COLUMNS, "S5 chance identity check")
    ensure_required_columns(inputs["leaderboard"], LEADERBOARD_REQUIRED_COLUMNS, "S6 retrieval leaderboard")

    ensure_non_null_keys(inputs["summary"], SLICE_KEY_COLUMNS, "S5 retrieval summary")
    ensure_non_null_keys(inputs["chance"], SLICE_KEY_COLUMNS, "S5 chance identity check")
    ensure_non_null_keys(inputs["leaderboard"], SLICE_KEY_COLUMNS, "S6 retrieval leaderboard")

    ensure_unique_keys(inputs["summary"], SLICE_KEY_COLUMNS, "S5 retrieval summary")
    ensure_unique_keys(inputs["chance"], SLICE_KEY_COLUMNS, "S5 chance identity check")
    ensure_unique_keys(inputs["leaderboard"], SLICE_KEY_COLUMNS, "S6 retrieval leaderboard")

    validate_summary_internal_consistency(inputs["summary"])
    compare_key_sets(inputs["summary"], inputs["chance"], SLICE_KEY_COLUMNS, "S5 summary vs chance check")
    validate_summary_leaderboard_consistency(inputs["summary"], inputs["leaderboard"])


def build_robustness_audit_table(inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Build the conservative A3 robustness table.

    This implementation emits one row per frozen
    `dataset, cell_line, direction, representation` slice. It does not invent
    additional robustness sub-slices whose logic is still under review.
    """
    observed = build_per_query_observed_context(inputs["per_query"])
    validated = validate_per_query_against_summary(inputs["summary"], observed)

    chance = inputs["chance"].copy()
    chance["chance_check_max_abs_delta"] = chance[
        ["abs_delta_mrr", "abs_delta_hit1", "abs_delta_hit5", "abs_delta_hit10"]
    ].max(axis=1)
    chance["chance_check_available_bool"] = True

    merged = validated.merge(
        chance[SLICE_KEY_COLUMNS + ["chance_check_available_bool", "chance_check_max_abs_delta"]],
        on=SLICE_KEY_COLUMNS,
        how="inner",
        validate="one_to_one",
    ).merge(
        inputs["leaderboard"][SLICE_KEY_COLUMNS + ["rank_by_mean_mrr_corrected"]],
        on=SLICE_KEY_COLUMNS,
        how="inner",
        validate="one_to_one",
    )

    output = pd.DataFrame(
        {
            "dataset": merged["dataset"],
            "cell_line": merged["cell_line"],
            "direction": merged["direction"],
            "representation": merged["representation"],
            "robustness_slice_label": REVIEWED_SLICE_LABEL,
            "gallery_definition_id": merged["gallery_definition_id"],
            "pos_definition_id": merged["pos_definition_id"],
            "rank_by_mean_mrr_corrected": merged["rank_by_mean_mrr_corrected"],
            "n_total": merged["n_total"],
            "n_valid": merged["n_valid"],
            "n_valid_observed": merged["n_valid_observed"],
            "n_excluded_missing_metric_or_mpos0": merged["n_excluded_missing_metric_or_mpos0"],
            "N_gallery_mean": merged["N_gallery_mean"],
            "N_gallery_min_observed": merged["N_gallery_min_observed"],
            "N_gallery_max": merged["N_gallery_max"],
            "m_pos_mean": merged["m_pos_mean"],
            "m_pos_p50": merged["m_pos_p50"],
            "m_pos_p90": merged["m_pos_p90"],
            "m_pos_min_observed": merged["m_pos_min_observed"],
            "m_pos_max_observed": merged["m_pos_max_observed"],
            "gallery_ids_hash_unique_n": merged["gallery_ids_hash_unique_n"],
            "mean_mrr_corrected": merged["mean_mrr_corrected"],
            "mean_hit1_corrected": merged["mean_hit1_corrected"],
            "mean_hit5_corrected": merged["mean_hit5_corrected"],
            "mean_hit10_corrected": merged["mean_hit10_corrected"],
            "chance_check_available_bool": merged["chance_check_available_bool"],
            "chance_check_max_abs_delta": merged["chance_check_max_abs_delta"],
            "status": "successful",
            "status_reason": "validated_reviewed_primary_direction_slice",
            "provenance_notes": (
                "Conservative A3 implementation: one row per frozen "
                "dataset-cell_line-direction-representation slice from S5/S6. "
                "Additional robustness sub-slice logic remains TODO pending "
                "implementation review."
            ),
        }
    )
    return output.sort_values(SLICE_KEY_COLUMNS).reset_index(drop=True)[OUTPUT_COLUMNS]


def write_output(frame: pd.DataFrame, output_path: Path) -> None:
    """Write the manuscript-support table without touching upstream runs."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)


def main() -> int:
    """Run the A3 skeleton through input loading, validation, and placeholder output."""
    args = parse_args()
    project_root = resolve_path(Path.cwd(), args.project_root)
    inputs = load_reviewed_inputs(args, project_root)
    validate_inputs(inputs)
    output = build_robustness_audit_table(inputs)
    write_output(output, resolve_path(project_root, args.output_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
