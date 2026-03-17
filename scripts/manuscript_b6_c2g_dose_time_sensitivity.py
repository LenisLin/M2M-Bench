#!/usr/bin/env python3
"""
Build the active B6 C2G-only dose/time sensitivity table.

This is a downstream manuscript-support table only. It keeps one long-form row
per C2G query and corrected retrieval metric after an explicit
`query_row_id -> row_id` metadata join.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet"
)
DEFAULT_LINCS_META_PATH = Path("data/task2_snapshot_v2/lincs/derived/delta_meta.csv")
DEFAULT_SCPERTURB_META_PATH = Path("data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv")
DEFAULT_OUTPUT_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/sensitivity/task2_c2g_dose_time_sensitivity.csv"
)

RETRIEVAL_METRICS = [
    "mrr_corrected",
    "hit1_corrected",
    "hit5_corrected",
    "hit10_corrected",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the B6 C2G dose/time sensitivity table.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--per-query-path", type=Path, default=DEFAULT_PER_QUERY_PATH)
    parser.add_argument("--lincs-meta-path", type=Path, default=DEFAULT_LINCS_META_PATH)
    parser.add_argument("--scperturb-meta-path", type=Path, default=DEFAULT_SCPERTURB_META_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--summary-output-path", type=Path, default=None)
    return parser.parse_args()


def resolve_path(project_root: Path, raw_path: Path) -> Path:
    return raw_path if raw_path.is_absolute() else (project_root / raw_path)


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")


def derive_summary_output_path(output_path: Path) -> Path:
    return output_path.with_name(f"{output_path.stem}_summary.csv")


def normalize_target_label(target_token: object, target_tokens: object) -> str:
    target_token_str = "" if pd.isna(target_token) else str(target_token).strip()
    if target_token_str and target_token_str != "NA":
        return target_token_str
    target_tokens_str = "" if pd.isna(target_tokens) else str(target_tokens).strip()
    return target_tokens_str


def same_numeric_value(left: object, right: object) -> bool:
    if pd.isna(left) and pd.isna(right):
        return True
    if pd.isna(left) or pd.isna(right):
        return False
    return math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=1e-12)


def load_metadata(path: Path, label: str) -> pd.DataFrame:
    require_path(path, label)
    frame = pd.read_csv(path)
    required_columns = {"row_id", "cell_line", "time", "dose_value"}
    missing = sorted(required_columns - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")
    if frame.duplicated(subset=["row_id"]).any():
        raise ValueError(f"{label} has duplicate `row_id` values")
    return frame[["row_id", "cell_line", "time", "dose_value"]].rename(
        columns={
            "row_id": "query_row_id",
            "cell_line": "metadata_cell_line",
            "time": "metadata_time",
            "dose_value": "metadata_dose_value",
        }
    )


def build_batch_rows(frame: pd.DataFrame, metadata_lookup: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frame = frame.loc[frame["direction"].eq("C2G")].copy()
    if frame.empty:
        return pd.DataFrame()

    parts: list[pd.DataFrame] = []
    for dataset, metadata in metadata_lookup.items():
        subset = frame.loc[frame["dataset"].eq(dataset)].copy()
        if subset.empty:
            continue

        joined = subset.merge(metadata, on="query_row_id", how="left", validate="many_to_one", indicator=True)
        joined["metadata_row_found_bool"] = joined["_merge"].eq("both")
        joined["cell_line_match_bool"] = (
            joined["metadata_row_found_bool"] & joined["cell_line"].eq(joined["metadata_cell_line"])
        )
        if (~joined.loc[joined["metadata_row_found_bool"], "cell_line_match_bool"]).any():
            raise ValueError(f"{dataset} metadata join produced cell-line mismatches")

        query_time_present = joined["query_time"].notna()
        query_dose_present = joined["query_dose_value"].notna()
        time_conflict = joined["metadata_row_found_bool"] & query_time_present & ~joined.apply(
            lambda row: same_numeric_value(row["query_time"], row["metadata_time"]),
            axis=1,
        )
        dose_conflict = joined["metadata_row_found_bool"] & query_dose_present & ~joined.apply(
            lambda row: same_numeric_value(row["query_dose_value"], row["metadata_dose_value"]),
            axis=1,
        )
        if time_conflict.any():
            raise ValueError(f"{dataset} metadata join produced query/meta time conflicts")
        if dose_conflict.any():
            raise ValueError(f"{dataset} metadata join produced query/meta dose conflicts")

        joined["target"] = [
            normalize_target_label(token, tokens)
            for token, tokens in zip(joined["query_target_token"], joined["query_target_tokens"])
        ]
        joined["time"] = joined["metadata_time"]
        joined["dose_value"] = joined["metadata_dose_value"]
        joined["covariate_status"] = "complete"
        joined.loc[~joined["metadata_row_found_bool"], "covariate_status"] = "metadata_missing"
        joined.loc[
            joined["metadata_row_found_bool"] & joined["time"].isna() & joined["dose_value"].isna(),
            "covariate_status",
        ] = "missing_time_and_dose"
        joined.loc[
            joined["metadata_row_found_bool"] & joined["time"].isna() & joined["dose_value"].notna(),
            "covariate_status",
        ] = "missing_time"
        joined.loc[
            joined["metadata_row_found_bool"] & joined["time"].notna() & joined["dose_value"].isna(),
            "covariate_status",
        ] = "missing_dose"

        joined["covariate_slice_id"] = joined.apply(
            lambda row: (
                f"{row['dataset']}||{row['cell_line']}||{row['representation']}||"
                f"time={row['time']}||dose={row['dose_value']}"
            )
            if row["covariate_status"] == "complete"
            else pd.NA,
            axis=1,
        )

        parts.append(
            joined[
                [
                    "dataset",
                    "cell_line",
                    "direction",
                    "representation",
                    "query_row_id",
                    "query_uid",
                    "target",
                    "time",
                    "dose_value",
                    "covariate_status",
                    "covariate_slice_id",
                    *RETRIEVAL_METRICS,
                ]
            ].copy()
        )

    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def long_format_with_slice_counts(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame()

    complete_slice_counts = (
        frame.loc[frame["covariate_status"].eq("complete")]
        .groupby(["dataset", "cell_line", "representation", "covariate_slice_id"], dropna=False)
        .size()
        .rename("n_queries_in_covariate_slice")
        .reset_index()
    )
    merged = frame.merge(
        complete_slice_counts,
        on=["dataset", "cell_line", "representation", "covariate_slice_id"],
        how="left",
        validate="many_to_one",
    )
    long_frame = merged.melt(
        id_vars=[
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "query_row_id",
            "query_uid",
            "target",
            "time",
            "dose_value",
            "covariate_status",
            "covariate_slice_id",
            "n_queries_in_covariate_slice",
        ],
        value_vars=RETRIEVAL_METRICS,
        var_name="retrieval_metric_name",
        value_name="retrieval_value",
    )
    long_frame["sensitivity_scope"] = "C2G_only_downstream_dose_time_sensitivity"
    long_frame["sensitivity_only_bool"] = True
    long_frame["not_benchmark_wide_bool"] = True
    long_frame["provenance_note"] = (
        "Rows come from S5 C2G retrieval per-query output after a dataset-routed "
        "`query_row_id -> row_id` join into frozen Task2 chemical metadata. "
        "The `target` label uses `query_target_token` when present and otherwise "
        "falls back to the exact `query_target_tokens` string."
    )
    return long_frame.sort_values(
        ["dataset", "cell_line", "representation", "query_row_id", "retrieval_metric_name"]
    ).reset_index(drop=True)


def build_summary_table(sensitivity_table: pd.DataFrame) -> pd.DataFrame:
    summary = (
        sensitivity_table.groupby(
            [
                "dataset",
                "cell_line",
                "target",
                "retrieval_metric_name",
                "dose_value",
                "time",
                "representation",
            ],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_queries=("query_row_id", "size"),
            mean_value=("retrieval_value", "mean"),
            median_value=("retrieval_value", "median"),
        )
        .reset_index()
    )
    return summary[
        [
            "dataset",
            "cell_line",
            "target",
            "retrieval_metric_name",
            "dose_value",
            "time",
            "n_queries",
            "mean_value",
            "median_value",
            "representation",
        ]
    ]


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    per_query_path = resolve_path(project_root, args.per_query_path)
    lincs_meta_path = resolve_path(project_root, args.lincs_meta_path)
    scperturb_meta_path = resolve_path(project_root, args.scperturb_meta_path)
    output_path = resolve_path(project_root, args.output_path)
    summary_output_path = (
        resolve_path(project_root, args.summary_output_path)
        if args.summary_output_path is not None
        else derive_summary_output_path(output_path)
    )

    require_path(per_query_path, "Task2 retrieval per-query parquet")
    metadata_lookup = {
        "LINCS": load_metadata(lincs_meta_path, "LINCS Task2 delta_meta"),
        "scPerturb": load_metadata(scperturb_meta_path, "scPerturb Task2 delta_meta"),
    }

    parquet_file = pq.ParquetFile(per_query_path)
    columns = [
        "dataset",
        "cell_line",
        "direction",
        "representation",
        "query_row_id",
        "query_uid",
        "query_target_token",
        "query_target_tokens",
        "query_time",
        "query_dose_value",
        *RETRIEVAL_METRICS,
    ]
    batch_frames: list[pd.DataFrame] = []
    for batch in parquet_file.iter_batches(batch_size=200_000, columns=columns):
        batch_frame = batch.to_pandas()
        batch_rows = build_batch_rows(batch_frame, metadata_lookup)
        if not batch_rows.empty:
            batch_frames.append(batch_rows)

    if not batch_frames:
        raise ValueError("B6 produced no C2G query rows from the reviewed Task2 retrieval output.")

    c2g_queries = pd.concat(batch_frames, ignore_index=True)
    sensitivity_table = long_format_with_slice_counts(c2g_queries)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sensitivity_table.to_csv(output_path, index=False)
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_table = build_summary_table(sensitivity_table)
    summary_table.to_csv(summary_output_path, index=False)

    print(f"Wrote {len(sensitivity_table):,} B6 sensitivity rows to {output_path}")
    print(f"Wrote {len(summary_table):,} B6 summary rows to {summary_output_path}")
    print("C2G queries:", c2g_queries["query_row_id"].nunique())
    print(
        "Complete covariate slices:",
        c2g_queries.loc[c2g_queries["covariate_status"].eq("complete"), "covariate_slice_id"].nunique(),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
