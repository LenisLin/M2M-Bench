#!/usr/bin/env python3
"""
Build the active A1 Task1 internal-vs-cross bridge table.

Status:
- support-only manuscript builder

Consumed by:
- Figure 2 internal-to-cross support and retained A1 backfills

Architecture:
- not a canonical Figure 2 object builder; see `scripts/ARCHITECTURE.md`

The frozen comparison unit is the shared `(dataset, cell_line, target)` group.
`representation` is retained only as source detail for explicit row-level
provenance and never defines the shared-group universe.

This retained manuscript-support builder is live for Figure 2 only.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from manuscript_task1_group_support import (
    COMMON_REPRESENTATIONS,
    DETAIL_COLUMNS,
    GROUP_METRICS,
    TRIPLET_COLUMNS,
    compute_task1_internal_group_metrics,
    load_task1_cross_group_metrics,
    source_dataset_from_cross_direction,
)


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_INTERNAL_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_per_query.parquet"
)
DEFAULT_CROSS_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_per_query.parquet"
)
DEFAULT_TASK1_GROUP_CROSS_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_group_cross.parquet"
)
DEFAULT_OUTPUT_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/group_bridge/task1_internal_vs_cross_group_bridge.csv"
)

RETRIEVAL_METRICS = [
    "mrr_corrected",
    "hit1_corrected",
    "hit5_corrected",
    "hit10_corrected",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the Figure 2 A1 Task1 internal-vs-cross group bridge table.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--internal-per-query-path", type=Path, default=DEFAULT_INTERNAL_PER_QUERY_PATH)
    parser.add_argument("--cross-per-query-path", type=Path, default=DEFAULT_CROSS_PER_QUERY_PATH)
    parser.add_argument("--task1-group-cross-path", type=Path, default=DEFAULT_TASK1_GROUP_CROSS_PATH)
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


def init_metric_state() -> dict[str, float | int]:
    state: dict[str, float | int] = {"n_queries": 0}
    for metric_name in RETRIEVAL_METRICS:
        state[f"{metric_name}_sum"] = 0.0
        state[f"{metric_name}_count"] = 0
    return state


def update_state(
    aggregate: dict[tuple[str, str, str, str], dict[str, float | int]],
    batch_summary: pd.DataFrame,
) -> None:
    for row in batch_summary.itertuples(index=False):
        key = (row.dataset, row.cell_line, row.target, row.representation_detail)
        state = aggregate.setdefault(key, init_metric_state())
        state["n_queries"] += int(row.n_queries)
        for metric_name in RETRIEVAL_METRICS:
            metric_sum = getattr(row, f"{metric_name}_sum")
            metric_count = getattr(row, f"{metric_name}_count")
            if not math.isnan(float(metric_sum)):
                state[f"{metric_name}_sum"] += float(metric_sum)
            state[f"{metric_name}_count"] += int(metric_count)


def finalize_aggregate(
    aggregate: dict[tuple[str, str, str, str], dict[str, float | int]],
    n_queries_column: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for key, state in aggregate.items():
        dataset, cell_line, target, representation_detail = key
        row: dict[str, object] = {
            "dataset": dataset,
            "cell_line": cell_line,
            "target": target,
            "representation_detail": representation_detail,
            n_queries_column: int(state["n_queries"]),
        }
        for metric_name in RETRIEVAL_METRICS:
            metric_sum = float(state[f"{metric_name}_sum"])
            metric_count = int(state[f"{metric_name}_count"])
            row[metric_name] = metric_sum / metric_count if metric_count else pd.NA
            row[f"{metric_name}_n_valid"] = metric_count
        rows.append(row)

    return pd.DataFrame(rows)


def summarize_batch(frame: pd.DataFrame, dataset_mapper) -> pd.DataFrame:
    frame = frame.loc[
        frame["perturbation_type"].eq("Genetic")
        & frame["representation"].isin(COMMON_REPRESENTATIONS)
    ].copy()
    if frame.empty:
        return pd.DataFrame()

    frame["dataset"] = frame["dataset_or_direction"].map(dataset_mapper)
    frame = frame.loc[frame["dataset"].notna()].copy()
    if frame.empty:
        return pd.DataFrame()

    frame["target"] = frame["target_token"].astype(str)
    frame["representation_detail"] = frame["representation"].astype(str)
    grouped = frame.groupby(DETAIL_COLUMNS, dropna=False)
    summary = grouped.size().rename("n_queries").to_frame()
    for metric_name in RETRIEVAL_METRICS:
        summary[f"{metric_name}_sum"] = grouped[metric_name].sum(min_count=1)
        summary[f"{metric_name}_count"] = grouped[metric_name].count()
    return summary.reset_index()


def aggregate_retrieval_from_parquet(path: Path, dataset_mapper) -> pd.DataFrame:
    require_path(path, f"retrieval parquet: {path.name}")
    parquet_file = pq.ParquetFile(path)
    columns = [
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "cell_line",
        "target_token",
        *RETRIEVAL_METRICS,
    ]
    aggregate: dict[tuple[str, str, str, str], dict[str, float | int]] = {}
    for batch in parquet_file.iter_batches(batch_size=200_000, columns=columns):
        batch_frame = batch.to_pandas()
        batch_summary = summarize_batch(batch_frame, dataset_mapper)
        if not batch_summary.empty:
            update_state(aggregate, batch_summary)
    frame = finalize_aggregate(aggregate, "n_queries")
    if frame.empty:
        return frame
    return frame.sort_values(DETAIL_COLUMNS, kind="mergesort").reset_index(drop=True)


def restrict_to_shared_triplets(left: pd.DataFrame, right: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    shared_triplets = left[TRIPLET_COLUMNS].drop_duplicates().merge(
        right[TRIPLET_COLUMNS].drop_duplicates(),
        on=TRIPLET_COLUMNS,
        how="inner",
        validate="one_to_one",
    )
    if shared_triplets.empty:
        return left.iloc[0:0].copy(), right.iloc[0:0].copy()
    left_shared = left.merge(shared_triplets, on=TRIPLET_COLUMNS, how="inner")
    right_shared = right.merge(shared_triplets, on=TRIPLET_COLUMNS, how="inner")
    return left_shared, right_shared


def conservative_group_support(*values: object) -> object:
    finite = [int(value) for value in values if pd.notna(value)]
    if not finite:
        return pd.NA
    return min(finite)


def build_retrieval_bridge_rows(internal_summary: pd.DataFrame, cross_summary: pd.DataFrame) -> list[dict[str, object]]:
    internal_shared, cross_shared = restrict_to_shared_triplets(internal_summary, cross_summary)
    merged = internal_shared.merge(
        cross_shared,
        on=DETAIL_COLUMNS,
        how="inner",
        suffixes=("_internal", "_cross"),
        validate="one_to_one",
    )

    rows: list[dict[str, object]] = []
    for row in merged.itertuples(index=False):
        for metric_name in RETRIEVAL_METRICS:
            internal_value = getattr(row, f"{metric_name}_internal")
            cross_value = getattr(row, f"{metric_name}_cross")
            rows.append(
                {
                    "dataset": row.dataset,
                    "cell_line": row.cell_line,
                    "target": row.target,
                    "representation_detail": row.representation_detail,
                    "metric_family": "retrieval_corrected",
                    "metric_name": metric_name,
                    "internal_value": internal_value,
                    "cross_value": cross_value,
                    "cross_minus_internal": (
                        cross_value - internal_value
                        if pd.notna(internal_value) and pd.notna(cross_value)
                        else pd.NA
                    ),
                    "comparison_scope": "dataset_cell_line_target_group_bridge",
                    "support_n_internal": row.n_queries_internal,
                    "support_n_cross": row.n_queries_cross,
                    "internal_support_detail": f"n_queries={row.n_queries_internal}",
                    "cross_support_detail": f"n_queries={row.n_queries_cross}",
                    "internal_n_valid": getattr(row, f"{metric_name}_n_valid_internal"),
                    "cross_n_valid": getattr(row, f"{metric_name}_n_valid_cross"),
                    "provenance_note": (
                        "Task1 internal retrieval values are grouped corrected genetic-query "
                        "means over shared triplets; Task1 cross retrieval values use the same "
                        "triplet universe after source-dataset mapping. "
                        f"representation_detail={row.representation_detail} is retained only as source detail."
                    ),
                }
            )
    return rows


def build_group_bridge_rows(internal_group: pd.DataFrame, cross_group: pd.DataFrame) -> list[dict[str, object]]:
    internal_shared, cross_shared = restrict_to_shared_triplets(internal_group, cross_group)
    merged = internal_shared.merge(
        cross_shared,
        on=DETAIL_COLUMNS,
        how="inner",
        suffixes=("_internal", "_cross"),
        validate="one_to_one",
    )

    rows: list[dict[str, object]] = []
    for row in merged.itertuples(index=False):
        for metric_name in GROUP_METRICS:
            internal_value = getattr(row, f"{metric_name}_internal")
            cross_value = getattr(row, f"{metric_name}_cross")
            rows.append(
                {
                    "dataset": row.dataset,
                    "cell_line": row.cell_line,
                    "target": row.target,
                    "representation_detail": row.representation_detail,
                    "metric_family": "group_concordance",
                    "metric_name": metric_name,
                    "internal_value": internal_value,
                    "cross_value": cross_value,
                    "cross_minus_internal": (
                        cross_value - internal_value
                        if pd.notna(internal_value) and pd.notna(cross_value)
                        else pd.NA
                    ),
                    "comparison_scope": "dataset_cell_line_target_group_bridge",
                    "support_n_internal": conservative_group_support(row.n_A, row.n_B),
                    "support_n_cross": conservative_group_support(row.n_L, row.n_S),
                    "internal_support_detail": (
                        f"n_total={row.n_total};n_A={row.n_A};n_B={row.n_B};"
                        f"n_A_sub={row.n_A_sub};n_B_sub={row.n_B_sub}"
                    ),
                    "cross_support_detail": (
                        f"pair_direction={row.task1_pair_direction};n_L={row.n_L};n_S={row.n_S};"
                        f"n_L_sub={row.n_L_sub};n_S_sub={row.n_S_sub}"
                    ),
                    "internal_n_valid": int(pd.notna(internal_value)),
                    "cross_n_valid": int(pd.notna(cross_value)),
                    "provenance_note": (
                        "Task1 internal group metrics use the existing S1 split-half group logic "
                        "recomputed from the frozen Task1 snapshot; Task1 cross group metrics come "
                        "from S2 cross-group rows duplicated onto each paired dataset only after "
                        "the shared triplet universe is fixed. "
                        f"representation_detail={row.representation_detail} is retained only as source detail."
                    ),
                }
            )
    return rows


def build_bridge_table(
    internal_retrieval: pd.DataFrame,
    cross_retrieval: pd.DataFrame,
    internal_group: pd.DataFrame,
    cross_group: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rows.extend(build_retrieval_bridge_rows(internal_retrieval, cross_retrieval))
    rows.extend(build_group_bridge_rows(internal_group, cross_group))
    if not rows:
        raise ValueError("A1 bridge produced no lawful shared `(dataset, cell_line, target)` rows.")

    return (
        pd.DataFrame(rows)
        .sort_values(
            ["dataset", "cell_line", "target", "metric_family", "metric_name", "representation_detail"],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )


def build_summary_table(bridge_table: pd.DataFrame) -> pd.DataFrame:
    summary_source = bridge_table.loc[
        bridge_table[["internal_value", "cross_value", "cross_minus_internal"]].notna().all(axis=1)
    ].copy()
    summary = (
        summary_source.groupby(
            ["metric_family", "metric_name", "dataset", "cell_line", "representation_detail"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_groups=("target", "size"),
            internal_mean=("internal_value", "mean"),
            cross_mean=("cross_value", "mean"),
            delta_mean=("cross_minus_internal", "mean"),
            delta_median=("cross_minus_internal", "median"),
        )
        .reset_index()
    )
    return summary[
        [
            "metric_name",
            "dataset",
            "cell_line",
            "n_groups",
            "internal_mean",
            "cross_mean",
            "delta_mean",
            "delta_median",
            "metric_family",
            "representation_detail",
        ]
    ]


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    internal_path = resolve_path(project_root, args.internal_per_query_path)
    cross_path = resolve_path(project_root, args.cross_per_query_path)
    task1_group_cross_path = resolve_path(project_root, args.task1_group_cross_path)
    output_path = resolve_path(project_root, args.output_path)
    summary_output_path = (
        resolve_path(project_root, args.summary_output_path)
        if args.summary_output_path is not None
        else derive_summary_output_path(output_path)
    )

    internal_retrieval = aggregate_retrieval_from_parquet(internal_path, lambda value: value)
    cross_retrieval = aggregate_retrieval_from_parquet(cross_path, source_dataset_from_cross_direction)
    require_path(task1_group_cross_path, "Task1 cross group parquet")
    cross_group = load_task1_cross_group_metrics(task1_group_cross_path)
    internal_group = compute_task1_internal_group_metrics(
        project_root,
        allowed_triplets=cross_group[TRIPLET_COLUMNS].drop_duplicates(),
    )

    bridge_table = build_bridge_table(
        internal_retrieval=internal_retrieval,
        cross_retrieval=cross_retrieval,
        internal_group=internal_group,
        cross_group=cross_group,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    bridge_table.to_csv(output_path, index=False)
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_table = build_summary_table(bridge_table)
    summary_table.to_csv(summary_output_path, index=False)

    print(f"Wrote {len(bridge_table):,} A1 bridge rows to {output_path}")
    print(f"Wrote {len(summary_table):,} A1 summary rows to {summary_output_path}")
    print("Shared triplets:", bridge_table[TRIPLET_COLUMNS].drop_duplicates().shape[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
