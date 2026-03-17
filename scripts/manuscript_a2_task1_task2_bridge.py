#!/usr/bin/env python3
"""
Build the active A2 Task1-vs-Task2 bridge table.

The frozen comparison unit is the shared `(dataset, cell_line, target)` group.
`representation` is retained only as source detail for explicit row-level
provenance and never defines the shared-group universe.
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
    TRIPLET_COLUMNS,
    load_task1_cross_group_metrics,
    source_dataset_from_cross_direction,
)


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_TASK1_INTERNAL_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_per_query.parquet"
)
DEFAULT_TASK1_CROSS_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_per_query.parquet"
)
DEFAULT_TASK1_GROUP_CROSS_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_group_cross.parquet"
)
DEFAULT_TASK2_GROUP_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv"
)
DEFAULT_TASK2_RETRIEVAL_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet"
)
DEFAULT_OUTPUT_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/group_bridge/task1_task2_group_bridge.csv"
)

RETRIEVAL_METRICS = [
    "mrr_corrected",
    "hit1_corrected",
    "hit5_corrected",
    "hit10_corrected",
]
GROUP_METRIC_PAIRS = [
    ("cosine", "cosine_centroid", "cosine_valid_bool"),
    ("pcc", "pcc_centroid", "pcc_valid_bool"),
    ("edist", "edist_biascorr", "edist_valid_bool"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the A2 Task1-vs-Task2 group bridge table.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--task1-internal-per-query-path", type=Path, default=DEFAULT_TASK1_INTERNAL_PER_QUERY_PATH)
    parser.add_argument("--task1-cross-per-query-path", type=Path, default=DEFAULT_TASK1_CROSS_PER_QUERY_PATH)
    parser.add_argument("--task1-group-cross-path", type=Path, default=DEFAULT_TASK1_GROUP_CROSS_PATH)
    parser.add_argument("--task2-group-path", type=Path, default=DEFAULT_TASK2_GROUP_PATH)
    parser.add_argument("--task2-retrieval-per-query-path", type=Path, default=DEFAULT_TASK2_RETRIEVAL_PER_QUERY_PATH)
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


def summarize_task1_batch(frame: pd.DataFrame, dataset_mapper) -> pd.DataFrame:
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


def aggregate_task1_retrieval(path: Path, dataset_mapper, n_queries_column: str) -> pd.DataFrame:
    require_path(path, f"Task1 retrieval parquet: {path.name}")
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
        batch_summary = summarize_task1_batch(batch_frame, dataset_mapper)
        if not batch_summary.empty:
            update_state(aggregate, batch_summary)
    frame = finalize_aggregate(aggregate, n_queries_column)
    if frame.empty:
        return frame
    return frame.sort_values(DETAIL_COLUMNS, kind="mergesort").reset_index(drop=True)


def summarize_task2_retrieval_batch(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.loc[
        frame["direction"].eq("G2C")
        & frame["representation"].isin(COMMON_REPRESENTATIONS)
    ].copy()
    if frame.empty:
        return pd.DataFrame()

    frame["target"] = frame["query_target_token"].astype(str)
    frame = frame.loc[
        frame["target"].notna()
        & frame["target"].ne("NA")
        & frame["target"].ne("")
    ].copy()
    if frame.empty:
        return pd.DataFrame()

    frame["representation_detail"] = frame["representation"].astype(str)
    grouped = frame.groupby(DETAIL_COLUMNS, dropna=False)
    summary = grouped.size().rename("n_queries").to_frame()
    for metric_name in RETRIEVAL_METRICS:
        summary[f"{metric_name}_sum"] = grouped[metric_name].sum(min_count=1)
        summary[f"{metric_name}_count"] = grouped[metric_name].count()
    return summary.reset_index()


def aggregate_task2_retrieval(path: Path) -> pd.DataFrame:
    require_path(path, f"Task2 retrieval parquet: {path.name}")
    parquet_file = pq.ParquetFile(path)
    columns = [
        "dataset",
        "cell_line",
        "direction",
        "representation",
        "query_target_token",
        *RETRIEVAL_METRICS,
    ]
    aggregate: dict[tuple[str, str, str, str], dict[str, float | int]] = {}
    for batch in parquet_file.iter_batches(batch_size=200_000, columns=columns):
        batch_frame = batch.to_pandas()
        batch_summary = summarize_task2_retrieval_batch(batch_frame)
        if not batch_summary.empty:
            update_state(aggregate, batch_summary)
    frame = finalize_aggregate(aggregate, "task2_n_queries")
    if frame.empty:
        return frame
    return frame.sort_values(DETAIL_COLUMNS, kind="mergesort").reset_index(drop=True)


def load_task2_group(path: Path) -> pd.DataFrame:
    require_path(path, "Task2 group concordance CSV")
    frame = pd.read_csv(path)
    required_columns = {
        "dataset",
        "cell_line",
        "target_token",
        "representation",
        "n_chem_sub",
        "n_gen_sub",
        "cosine_centroid",
        "pcc_centroid",
        "edist_biascorr",
        "cosine_valid_bool",
        "pcc_valid_bool",
        "edist_valid_bool",
    }
    missing = sorted(required_columns - set(frame.columns))
    if missing:
        raise ValueError(f"Task2 group concordance is missing required columns: {missing}")
    frame = frame.loc[frame["representation"].isin(COMMON_REPRESENTATIONS)].copy()
    frame = frame.rename(
        columns={
            "target_token": "target",
            "representation": "representation_detail",
        }
    )
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


def build_retrieval_bridge_rows(
    task1_frame: pd.DataFrame,
    task1_scope: str,
    task2_frame: pd.DataFrame,
) -> list[dict[str, object]]:
    task1_shared, task2_shared = restrict_to_shared_triplets(task1_frame, task2_frame)
    merged = task1_shared.merge(
        task2_shared,
        on=DETAIL_COLUMNS,
        how="inner",
        suffixes=("_task1", "_task2"),
        validate="one_to_one",
    )

    rows: list[dict[str, object]] = []
    for row in merged.itertuples(index=False):
        for metric_name in RETRIEVAL_METRICS:
            task1_value = getattr(row, f"{metric_name}_task1")
            task2_value = getattr(row, f"{metric_name}_task2")
            rows.append(
                {
                    "dataset": row.dataset,
                    "cell_line": row.cell_line,
                    "target": row.target,
                    "representation_detail": row.representation_detail,
                    "metric_family": "retrieval_corrected",
                    "task1_scope": task1_scope,
                    "task2_direction": "G2C",
                    "task1_metric_name": metric_name,
                    "task1_value": task1_value,
                    "task2_metric_name": metric_name,
                    "task2_value": task2_value,
                    "task2_minus_task1": (
                        task2_value - task1_value
                        if pd.notna(task1_value) and pd.notna(task2_value)
                        else pd.NA
                    ),
                    "support_n_task1": row.task1_n_queries,
                    "support_n_task2": row.task2_n_queries,
                    "task1_support_detail": f"n_queries={row.task1_n_queries}",
                    "task2_support_detail": f"n_queries={row.task2_n_queries}",
                    "bridge_scope": "dataset_cell_line_target_group_bridge",
                    "provenance_note": (
                        f"Task1 {task1_scope} retrieval values and Task2 G2C retrieval values are "
                        "paired only after fixing the shared triplet universe; "
                        f"representation_detail={row.representation_detail} is retained only as source detail."
                    ),
                }
            )
    return rows


def build_group_bridge_rows(
    task1_group: pd.DataFrame,
    task1_scope: str,
    task2_group: pd.DataFrame,
) -> list[dict[str, object]]:
    task1_shared, task2_shared = restrict_to_shared_triplets(task1_group, task2_group)
    merged = task1_shared.merge(
        task2_shared,
        on=DETAIL_COLUMNS,
        how="inner",
        suffixes=("_task1", "_task2"),
        validate="one_to_one",
    )

    rows: list[dict[str, object]] = []
    for row in merged.itertuples(index=False):
        for task1_metric_name, task2_metric_name, task2_valid_column in GROUP_METRIC_PAIRS:
            if not bool(getattr(row, task2_valid_column)):
                continue
            task1_value = getattr(row, task1_metric_name)
            task2_value = getattr(row, task2_metric_name)
            rows.append(
                {
                    "dataset": row.dataset,
                    "cell_line": row.cell_line,
                    "target": row.target,
                    "representation_detail": row.representation_detail,
                    "metric_family": "group_concordance",
                    "task1_scope": task1_scope,
                    "task2_direction": pd.NA,
                    "task1_metric_name": task1_metric_name,
                    "task1_value": task1_value,
                    "task2_metric_name": task2_metric_name,
                    "task2_value": task2_value,
                    "task2_minus_task1": (
                        task2_value - task1_value
                        if pd.notna(task1_value) and pd.notna(task2_value)
                        else pd.NA
                    ),
                    "support_n_task1": (
                        conservative_group_support(row.n_A, row.n_B)
                        if task1_scope == "internal"
                        else conservative_group_support(row.n_L, row.n_S)
                    ),
                    "support_n_task2": conservative_group_support(row.n_chem_sub, row.n_gen_sub),
                    "task1_support_detail": (
                        f"n_total={row.n_total};n_A={row.n_A};n_B={row.n_B};"
                        f"n_A_sub={row.n_A_sub};n_B_sub={row.n_B_sub}"
                        if task1_scope == "internal"
                        else (
                            f"pair_direction={row.task1_pair_direction};n_L={row.n_L};n_S={row.n_S};"
                            f"n_L_sub={row.n_L_sub};n_S_sub={row.n_S_sub}"
                        )
                    ),
                    "task2_support_detail": f"n_chem_sub={row.n_chem_sub};n_gen_sub={row.n_gen_sub}",
                    "bridge_scope": "dataset_cell_line_target_group_bridge",
                    "provenance_note": (
                        f"Task1 {task1_scope} group metrics are paired to Task2 group concordance "
                        "only after fixing the shared triplet universe; "
                        f"representation_detail={row.representation_detail} is retained only as source detail."
                    ),
                }
            )
    return rows


def build_summary_table(bridge_table: pd.DataFrame) -> pd.DataFrame:
    summary_source = bridge_table.loc[
        bridge_table[["task1_value", "task2_value", "task2_minus_task1"]].notna().all(axis=1)
    ].copy()
    summary = (
        summary_source.groupby(
            [
                "metric_family",
                "task1_scope",
                "task1_metric_name",
                "task2_metric_name",
                "dataset",
                "cell_line",
                "representation_detail",
            ],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_groups=("target", "size"),
            task1_mean=("task1_value", "mean"),
            task2_mean=("task2_value", "mean"),
            delta_mean=("task2_minus_task1", "mean"),
            delta_median=("task2_minus_task1", "median"),
        )
        .reset_index()
    )
    return summary[
        [
            "task1_metric_name",
            "task2_metric_name",
            "dataset",
            "cell_line",
            "n_groups",
            "task1_mean",
            "task2_mean",
            "delta_mean",
            "delta_median",
            "metric_family",
            "task1_scope",
            "representation_detail",
        ]
    ]


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    task1_internal_path = resolve_path(project_root, args.task1_internal_per_query_path)
    task1_cross_path = resolve_path(project_root, args.task1_cross_per_query_path)
    task1_group_cross_path = resolve_path(project_root, args.task1_group_cross_path)
    task2_group_path = resolve_path(project_root, args.task2_group_path)
    task2_retrieval_path = resolve_path(project_root, args.task2_retrieval_per_query_path)
    output_path = resolve_path(project_root, args.output_path)
    summary_output_path = (
        resolve_path(project_root, args.summary_output_path)
        if args.summary_output_path is not None
        else derive_summary_output_path(output_path)
    )

    task1_internal_retrieval = aggregate_task1_retrieval(task1_internal_path, lambda value: value, "task1_n_queries")
    task1_cross_retrieval = aggregate_task1_retrieval(
        task1_cross_path,
        source_dataset_from_cross_direction,
        "task1_n_queries",
    )
    require_path(task1_group_cross_path, "Task1 cross group parquet")
    task1_cross_group = load_task1_cross_group_metrics(task1_group_cross_path)
    task2_group = load_task2_group(task2_group_path)
    task2_retrieval = aggregate_task2_retrieval(task2_retrieval_path)

    rows: list[dict[str, object]] = []
    rows.extend(build_retrieval_bridge_rows(task1_internal_retrieval, "internal", task2_retrieval))
    rows.extend(build_retrieval_bridge_rows(task1_cross_retrieval, "cross", task2_retrieval))
    rows.extend(build_group_bridge_rows(task1_cross_group, "cross", task2_group))

    if not rows:
        raise ValueError("A2 bridge produced no lawful shared `(dataset, cell_line, target)` rows.")

    bridge_table = (
        pd.DataFrame(rows)
        .sort_values(
            [
                "dataset",
                "cell_line",
                "target",
                "metric_family",
                "task1_scope",
                "task1_metric_name",
                "representation_detail",
            ],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    bridge_table.to_csv(output_path, index=False)
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_table = build_summary_table(bridge_table)
    summary_table.to_csv(summary_output_path, index=False)

    print(f"Wrote {len(bridge_table):,} A2 bridge rows to {output_path}")
    print(f"Wrote {len(summary_table):,} A2 summary rows to {summary_output_path}")
    print("Shared triplets:", bridge_table[TRIPLET_COLUMNS].drop_duplicates().shape[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
