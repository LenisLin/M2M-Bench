#!/usr/bin/env python3
"""
Phase 1 manuscript-support skeleton for B6: exact dose/time overlap audit.

Analysis purpose:
- Determine whether exact dose/time covariate overlap exists in evaluable Task2
  query pools before any controlled comparison is attempted.

Why this script exists in the manuscript workflow:
- The manuscript needs an explicit gate before any dose/time-controlled claim.
- This script is an overlap audit only. It exists to show whether a lawful
  comparison is supportable or whether the manuscript must stop at a limitation.

Lawful input scope:
- reviewed local S5 corrected multisource retrieval per-query output
- reviewed local LINCS Task2 `delta_meta.csv`
- reviewed local scPerturb Task2 `delta_meta.csv`

Prohibited behavior:
- do not implement the controlled comparison in this pass
- do not bin, round, infer, or harmonize dose/time values
- do not rescue missing covariates through hidden fallback rules
- do not mutate any upstream benchmark-stage artifacts

Intended output semantics:
- one audit row per lawful slice after explicit query-to-metadata alignment
- missing covariates remain explicit
- overlap status must use the reviewed vocabulary:
  `no_overlap`, `partial_overlap`, `usable_overlap`,
  `blocked_missing_covariates`
- limitation-only output is lawful if exact overlap is insufficient

What remains TODO in later implementation:
- implement the explicit `query_row_id -> row_id` metadata join
- implement missing-covariate accounting without silent rescue logic
- implement exact overlap counting on `dataset, cell_line, time, dose_value`
- freeze the support threshold that separates `partial_overlap` from
  `usable_overlap`
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
NAS_RUNS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs")
MANUSCRIPT_PHASE1_ROOT = NAS_RUNS_ROOT / "manuscript_phase1" / "b6_dose_time_overlap"

DEFAULT_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet"
)
DEFAULT_LINCS_META_PATH = Path("data/task2_snapshot_v2/lincs/derived/delta_meta.csv")
DEFAULT_SCPERTURB_META_PATH = Path("data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv")
DEFAULT_OUTPUT_PATH = MANUSCRIPT_PHASE1_ROOT / "task2_dose_time_overlap_audit.csv"

USABLE_OVERLAP_MIN_STRATA = 3
USABLE_OVERLAP_MIN_QUERY_UNITS = 10

PER_QUERY_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "query_row_id",
    "query_uid",
]
METADATA_REQUIRED_COLUMNS = [
    "row_id",
    "cell_line",
    "time",
    "dose_value",
]
ALIGNMENT_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "query_row_id",
    "query_uid",
    "metadata_row_found_bool",
    "metadata_cell_line",
    "metadata_cell_line_matches_query_bool",
    "time",
    "dose_value",
    "time_missing_bool",
    "dose_value_missing_bool",
]
OUTPUT_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "n_query_rows_total",
    "n_query_rows_joined",
    "n_query_rows_missing_metadata",
    "n_query_rows_missing_time",
    "n_query_rows_missing_dose_value",
    "n_query_rows_complete_covariates",
    "n_exact_overlap_strata",
    "n_query_units_in_exact_overlap_strata",
    "overlap_status",
    "blocked_reason",
    "support_count_basis",
    "notes",
]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the B6 scaffold."""
    parser = argparse.ArgumentParser(description="Phase 1 skeleton for B6 dose/time overlap audit.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--per-query-path", type=Path, default=DEFAULT_PER_QUERY_PATH)
    parser.add_argument("--lincs-meta-path", type=Path, default=DEFAULT_LINCS_META_PATH)
    parser.add_argument("--scperturb-meta-path", type=Path, default=DEFAULT_SCPERTURB_META_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    return parser.parse_args()


def resolve_path(project_root: Path, raw_path: Path) -> Path:
    """Resolve a CLI path against the requested project root."""
    return raw_path if raw_path.is_absolute() else (project_root / raw_path)


def load_table(path: Path, label: str) -> pd.DataFrame:
    """Load a CSV or parquet table without hidden coercion."""
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


def ensure_non_null_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    """Reject nulls in key columns before any join logic is attempted."""
    if frame[columns].isnull().any().any():
        raise ValueError(f"{label} contains null values in required key columns: {columns}")


def ensure_unique_key(frame: pd.DataFrame, key_columns: list[str], label: str) -> None:
    """Validate uniqueness where the frozen contract requires it."""
    duplicated = frame.duplicated(subset=key_columns, keep=False)
    if duplicated.any():
        raise ValueError(f"{label} has duplicate key rows on {key_columns}")


def load_reviewed_inputs(args: argparse.Namespace, project_root: Path) -> dict[str, pd.DataFrame]:
    """Load the reviewed local B6 inputs."""
    return {
        "per_query": load_table(resolve_path(project_root, args.per_query_path), "S5 retrieval per-query"),
        "lincs_meta": load_table(resolve_path(project_root, args.lincs_meta_path), "LINCS Task2 delta_meta"),
        "scperturb_meta": load_table(resolve_path(project_root, args.scperturb_meta_path), "scPerturb Task2 delta_meta"),
    }


def validate_inputs(inputs: dict[str, pd.DataFrame]) -> None:
    """Run the reviewed schema checks that are safe to freeze now."""
    ensure_required_columns(inputs["per_query"], PER_QUERY_REQUIRED_COLUMNS, "S5 retrieval per-query")
    ensure_required_columns(inputs["lincs_meta"], METADATA_REQUIRED_COLUMNS, "LINCS Task2 delta_meta")
    ensure_required_columns(inputs["scperturb_meta"], METADATA_REQUIRED_COLUMNS, "scPerturb Task2 delta_meta")

    ensure_non_null_columns(inputs["per_query"], ["dataset", "query_row_id"], "S5 retrieval per-query")
    ensure_non_null_columns(inputs["lincs_meta"], ["row_id"], "LINCS Task2 delta_meta")
    ensure_non_null_columns(inputs["scperturb_meta"], ["row_id"], "scPerturb Task2 delta_meta")

    ensure_unique_key(inputs["lincs_meta"], ["row_id"], "LINCS Task2 delta_meta")
    ensure_unique_key(inputs["scperturb_meta"], ["row_id"], "scPerturb Task2 delta_meta")


def build_query_metadata_alignment(inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Build the explicit dataset-routed query-to-metadata alignment table.

    This uses a left join only and keeps all query rows, including unmatched
    rows, so missing metadata remains auditable instead of being dropped.
    """
    per_query = inputs["per_query"]
    metadata_by_dataset = {
        "LINCS": inputs["lincs_meta"][["row_id", "cell_line", "time", "dose_value"]].copy(),
        "scPerturb": inputs["scperturb_meta"][["row_id", "cell_line", "time", "dose_value"]].copy(),
    }

    unknown_datasets = sorted(set(per_query["dataset"].unique()) - set(metadata_by_dataset))
    if unknown_datasets:
        raise ValueError(f"Unsupported datasets in S5 retrieval per-query: {unknown_datasets}")

    aligned_parts: list[pd.DataFrame] = []
    for dataset, metadata in metadata_by_dataset.items():
        query_subset = per_query.loc[per_query["dataset"].eq(dataset)].copy()
        if query_subset.empty:
            continue

        metadata_subset = metadata.rename(columns={"row_id": "query_row_id", "cell_line": "metadata_cell_line"})
        aligned = query_subset.merge(
            metadata_subset,
            on="query_row_id",
            how="left",
            validate="many_to_one",
            indicator=True,
        )
        if len(aligned) != len(query_subset):
            raise ValueError(f"{dataset} left join changed the query-row count unexpectedly")

        aligned["metadata_row_found_bool"] = aligned["_merge"].eq("both")
        aligned["metadata_cell_line_matches_query_bool"] = (
            aligned["metadata_row_found_bool"] & aligned["metadata_cell_line"].eq(aligned["cell_line"])
        )
        if (~aligned.loc[aligned["metadata_row_found_bool"], "metadata_cell_line_matches_query_bool"]).any():
            raise ValueError(f"{dataset} metadata join produced cell-line mismatches for matched rows")

        aligned["time_missing_bool"] = aligned["metadata_row_found_bool"] & aligned["time"].isna()
        aligned["dose_value_missing_bool"] = aligned["metadata_row_found_bool"] & aligned["dose_value"].isna()
        aligned_parts.append(aligned[ALIGNMENT_COLUMNS].copy())

    if not aligned_parts:
        return pd.DataFrame(columns=ALIGNMENT_COLUMNS)

    return pd.concat(aligned_parts, ignore_index=True).sort_values(
        ["dataset", "cell_line", "direction", "representation", "query_uid"]
    ).reset_index(drop=True)


def classify_overlap_status(
    n_complete_covariates: int,
    n_exact_overlap_strata: int,
    n_query_units_in_exact_overlap_strata: int,
) -> tuple[str, str]:
    """Apply the explicit conservative overlap-status rule for this pass."""
    if n_complete_covariates == 0:
        return "blocked_missing_covariates", "no_complete_time_and_dose_after_explicit_metadata_join"
    if n_exact_overlap_strata == 0:
        return "no_overlap", "not_blocked"
    if (
        n_exact_overlap_strata >= USABLE_OVERLAP_MIN_STRATA
        and n_query_units_in_exact_overlap_strata >= USABLE_OVERLAP_MIN_QUERY_UNITS
    ):
        return "usable_overlap", "not_blocked"
    return "partial_overlap", "not_blocked"


def build_overlap_audit_table(alignment: pd.DataFrame) -> pd.DataFrame:
    """
    Build the B6 overlap audit table.

    Support accounting is deduplicated at the query level across
    representations to avoid inflating overlap support. Representation is then
    retained in the output only as a labeling axis.
    """
    if alignment.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    query_unit_keys = ["dataset", "cell_line", "direction", "query_row_id", "query_uid"]
    base_slice_keys = ["dataset", "cell_line", "direction"]

    # Deduplicate at the query-unit level before overlap counting so the same
    # query repeated across representations does not inflate support.
    query_level = alignment.groupby(query_unit_keys, dropna=False).agg(
        metadata_row_found_min=("metadata_row_found_bool", "min"),
        metadata_row_found_max=("metadata_row_found_bool", "max"),
        metadata_cell_line_matches_query_min=("metadata_cell_line_matches_query_bool", "min"),
        metadata_cell_line_matches_query_max=("metadata_cell_line_matches_query_bool", "max"),
        time_first=("time", "first"),
        dose_value_first=("dose_value", "first"),
        time_nunique_dropna_false=("time", lambda s: s.nunique(dropna=False)),
        dose_value_nunique_dropna_false=("dose_value", lambda s: s.nunique(dropna=False)),
        representation_nunique=("representation", "nunique"),
    ).reset_index()

    if not query_level["metadata_row_found_min"].equals(query_level["metadata_row_found_max"]):
        raise ValueError("metadata_row_found_bool is not stable across representation duplicates for a query unit")
    if not query_level["metadata_cell_line_matches_query_min"].equals(
        query_level["metadata_cell_line_matches_query_max"]
    ):
        raise ValueError("metadata cell-line match flag is not stable across representation duplicates for a query unit")
    if (query_level["time_nunique_dropna_false"] > 1).any():
        raise ValueError("Joined `time` is not stable across representation duplicates for a query unit")
    if (query_level["dose_value_nunique_dropna_false"] > 1).any():
        raise ValueError("Joined `dose_value` is not stable across representation duplicates for a query unit")

    query_level["metadata_row_found_bool"] = query_level["metadata_row_found_min"].astype(bool)
    query_level["time"] = query_level["time_first"]
    query_level["dose_value"] = query_level["dose_value_first"]
    query_level["time_missing_bool"] = query_level["metadata_row_found_bool"] & query_level["time"].isna()
    query_level["dose_value_missing_bool"] = (
        query_level["metadata_row_found_bool"] & query_level["dose_value"].isna()
    )
    query_level["complete_covariates_bool"] = (
        query_level["metadata_row_found_bool"]
        & ~query_level["time_missing_bool"]
        & ~query_level["dose_value_missing_bool"]
    )

    base_rows: list[dict[str, object]] = []
    for slice_key, slice_frame in query_level.groupby(base_slice_keys, dropna=False):
        dataset, cell_line, direction = slice_key
        n_total = int(len(slice_frame))
        n_joined = int(slice_frame["metadata_row_found_bool"].sum())
        n_missing_metadata = int((~slice_frame["metadata_row_found_bool"]).sum())
        n_missing_time = int(
            (slice_frame["metadata_row_found_bool"] & slice_frame["time_missing_bool"]).sum()
        )
        n_missing_dose_value = int(
            (slice_frame["metadata_row_found_bool"] & slice_frame["dose_value_missing_bool"]).sum()
        )
        n_complete = int(slice_frame["complete_covariates_bool"].sum())

        complete_frame = slice_frame.loc[slice_frame["complete_covariates_bool"]].copy()
        if complete_frame.empty:
            n_exact_overlap_strata = 0
            n_query_units_in_exact_overlap_strata = 0
        else:
            overlap_counts = complete_frame.groupby(
                ["dataset", "cell_line", "time", "dose_value"],
                dropna=False,
            ).size()
            overlapping_counts = overlap_counts.loc[overlap_counts >= 2]
            n_exact_overlap_strata = int(len(overlapping_counts))
            n_query_units_in_exact_overlap_strata = int(overlapping_counts.sum())

        overlap_status, blocked_reason = classify_overlap_status(
            n_complete_covariates=n_complete,
            n_exact_overlap_strata=n_exact_overlap_strata,
            n_query_units_in_exact_overlap_strata=n_query_units_in_exact_overlap_strata,
        )
        base_rows.append(
            {
                "dataset": dataset,
                "cell_line": cell_line,
                "direction": direction,
                "n_query_rows_total": n_total,
                "n_query_rows_joined": n_joined,
                "n_query_rows_missing_metadata": n_missing_metadata,
                "n_query_rows_missing_time": n_missing_time,
                "n_query_rows_missing_dose_value": n_missing_dose_value,
                "n_query_rows_complete_covariates": n_complete,
                "n_exact_overlap_strata": n_exact_overlap_strata,
                "n_query_units_in_exact_overlap_strata": n_query_units_in_exact_overlap_strata,
                "overlap_status": overlap_status,
                "blocked_reason": blocked_reason,
                "support_count_basis": "deduplicated_query_pool_within_dataset_cell_line_direction",
                "notes": (
                    "Support counts are computed on a query-level deduplicated "
                    "pool across representations to avoid representation-driven "
                    "inflation. Exact overlap is checked only on "
                    "`dataset, cell_line, time, dose_value` with exact equality."
                ),
            }
        )

    base_output = pd.DataFrame(base_rows)
    representation_grid = alignment[["dataset", "cell_line", "direction", "representation"]].drop_duplicates()
    output = representation_grid.merge(
        base_output,
        on=base_slice_keys,
        how="left",
        validate="many_to_one",
    )
    return output.sort_values(["dataset", "cell_line", "direction", "representation"]).reset_index(drop=True)[
        OUTPUT_COLUMNS
    ]


def write_output(frame: pd.DataFrame, output_path: Path) -> None:
    """Write the manuscript-support audit table."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)


def main() -> int:
    """Run the B6 skeleton through loading, validation, and placeholder output."""
    args = parse_args()
    project_root = resolve_path(Path.cwd(), args.project_root)
    inputs = load_reviewed_inputs(args, project_root)
    validate_inputs(inputs)
    alignment = build_query_metadata_alignment(inputs)
    output = build_overlap_audit_table(alignment)
    write_output(output, resolve_path(project_root, args.output_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
