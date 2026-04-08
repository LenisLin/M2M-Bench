#!/usr/bin/env python3
"""
Materialize manuscript-facing comparison statistics under the frozen contract.

Status:
- active manuscript canonical builder

Manuscript role:
- builds downstream comparison statistics over the frozen manuscript object set

Architecture:
- see `scripts/ARCHITECTURE.md` for script-family classification

This layer is downstream only. It must not modify S1-S7 benchmark logic or
reinterpret metric orientation. It consumes the frozen manuscript-facing
objects plus the explicitly allowed lawful upstream files where the canonical
objects do not expose the correct test unit.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

try:
    import pyarrow.parquet as pq
except ModuleNotFoundError:
    pq = None

ANALYSIS_OUTPUT_FILENAMES = {"framework_analysis_manifest": "framework_analysis_manifest.json"}


def annotate_canonical_objects(entries: list[dict[str, object]]) -> list[dict[str, object]]:
    try:
        from manuscript_framework_analysis_objects import annotate_canonical_objects as _annotate
    except ModuleNotFoundError:
        from scripts.manuscript_framework_analysis_objects import annotate_canonical_objects as _annotate
    return _annotate(entries)


def build_task2_target_support_long(task2_group_root: Path, task2_retrieval_root: Path) -> pd.DataFrame:
    try:
        from manuscript_framework_analysis_objects import build_task2_target_support_long as _build
    except ModuleNotFoundError:
        from scripts.manuscript_framework_analysis_objects import build_task2_target_support_long as _build
    return _build(task2_group_root, task2_retrieval_root)


NAS_RUNS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs")
DEFAULT_MANUSCRIPT_ACTIVE_ROOT = NAS_RUNS_ROOT / "manuscript_active"
DEFAULT_ANALYSIS_ROOT = DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis"
DEFAULT_TASK2_RETRIEVAL_PER_QUERY = (
    NAS_RUNS_ROOT
    / "s5_multisource_impl_verify_20260311_a"
    / "s5_task2_retrieval_multisource"
    / "task2_retrieval_per_query.parquet"
)
DEFAULT_TASK1_BRIDGE_SUMMARY = (
    DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "group_bridge" / "task1_internal_vs_cross_group_bridge_summary.csv"
)
DEFAULT_TASK1_BRIDGE_DETAIL = (
    DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "group_bridge" / "task1_internal_vs_cross_group_bridge.csv"
)
DEFAULT_TASK2_GROUP_ROOT = (
    NAS_RUNS_ROOT
    / "s4_multisource_impl_verify_20260310_c"
    / "s4_task2_group_concordance_multisource"
)

LOW_N_THRESHOLD = 5
REQUIRED_OUTPUT_COLUMNS = [
    "comparison_spec_id",
    "comparison_family",
    "question_id",
    "figure_id",
    "analysis_family",
    "metric_name",
    "group_a",
    "group_b",
    "pairing_type",
    "pairing_key_definition",
    "n_test_units",
    "raw_p",
    "bh_q",
    "significant_bool",
    "effect_direction",
    "median_delta",
    "test_name",
    "test_status",
    "notes",
    "n_units",
    "p",
    "effect_size",
]
OUTPUT_COLUMNS = [
    *REQUIRED_OUTPUT_COLUMNS,
    "primary_source_object",
    "dataset_scope",
    "representation_scope",
    "direction_scope",
]
TASK2_RETRIEVAL_METRICS = ["mrr_corrected", "hit1_corrected", "hit5_corrected", "hit10_corrected"]


@dataclass(frozen=True)
class ComparisonSpec:
    comparison_spec_id: str
    comparison_family: str
    figure_id: str
    primary_source_object: str
    pairing_type: str
    pairing_key_definition: str
    question_id_resolver: Callable[[dict[str, object]], tuple[str, ...]]
    default_analysis_family: str | None = None
    notes_prefix: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build manuscript-facing comparison statistics.")
    parser.add_argument("--manuscript-active-root", type=Path, default=DEFAULT_MANUSCRIPT_ACTIVE_ROOT)
    parser.add_argument("--analysis-root", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    parser.add_argument("--task2-retrieval-per-query-path", type=Path, default=DEFAULT_TASK2_RETRIEVAL_PER_QUERY)
    parser.add_argument("--task2-group-root", type=Path, default=DEFAULT_TASK2_GROUP_ROOT)
    parser.add_argument("--task1-bridge-summary-path", type=Path, default=DEFAULT_TASK1_BRIDGE_SUMMARY)
    parser.add_argument("--task1-bridge-detail-path", type=Path, default=DEFAULT_TASK1_BRIDGE_DETAIL)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_ANALYSIS_ROOT / "manuscript_comparison_statistics.csv")
    parser.add_argument(
        "--manifest-path",
        type=Path,
        default=DEFAULT_ANALYSIS_ROOT / ANALYSIS_OUTPUT_FILENAMES["framework_analysis_manifest"],
    )
    return parser.parse_args()


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")


def read_csv_required(path: Path, required: set[str], label: str) -> pd.DataFrame:
    require_path(path, label)
    frame = pd.read_csv(path)
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")
    return frame


def dataset_scope_label(value: object) -> str:
    if value is None or pd.isna(value):
        return "ALL"
    text = str(value)
    return text if text else "ALL"


def direction_scope_label(value: object) -> str:
    if value is None or pd.isna(value):
        return "ALL"
    text = str(value)
    return text if text else "ALL"


def task1_metric_family_to_analysis_family(metric_family: object) -> str:
    text = str(metric_family)
    if text == "group_concordance":
        return "group_concordance"
    if text == "retrieval_corrected":
        return "retrieval"
    raise ValueError(f"Unexpected Task1 metric_family: {metric_family}")


def effect_direction_from_delta(median_delta: float) -> str:
    if median_delta > 0:
        return "group_a_gt_group_b"
    if median_delta < 0:
        return "group_b_gt_group_a"
    return "tie_or_zero_delta"


def bh_adjust(series: pd.Series) -> pd.Series:
    if series.empty:
        return pd.Series(dtype="float64")
    order = np.argsort(series.to_numpy(), kind="mergesort")
    ranked = series.to_numpy()[order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1, dtype="float64")
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype="float64")
    out[order] = adjusted
    return pd.Series(out, index=series.index, dtype="float64")


SIGNED_LOG_METRICS = {"edist", "edist_biascorr", "mean_edist_biascorr"}


def signed_log10_1p(values: pd.Series | np.ndarray | Sequence[float]) -> pd.Series:
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").astype("float64")
    return np.sign(numeric) * np.log10(1.0 + np.abs(numeric))


def transform_metric_series(metric_name: object, values: pd.Series) -> pd.Series:
    if str(metric_name) in SIGNED_LOG_METRICS:
        return signed_log10_1p(values)
    return pd.to_numeric(values, errors="coerce")


def apply_metric_transform(long_frame: pd.DataFrame, *, metric_column: str, value_column: str) -> pd.DataFrame:
    out = long_frame.copy()
    transformed = pd.Series(np.nan, index=out.index, dtype="float64")
    for metric_name, idx in out.groupby(metric_column, dropna=False, sort=False).groups.items():
        transformed.loc[idx] = transform_metric_series(metric_name, out.loc[idx, value_column]).astype("float64")
    out[value_column] = transformed
    return out


def assert_unique_rows(frame: pd.DataFrame, key_columns: Sequence[str], label: str) -> None:
    if frame.empty:
        return
    dup_mask = frame.duplicated(subset=list(key_columns), keep=False)
    if not dup_mask.any():
        return
    examples = frame.loc[dup_mask, list(key_columns)].head(10).to_dict("records")
    raise ValueError(f"{label} has duplicate rows on key={list(key_columns)}; examples={examples}")


def validate_output_schema(frame: pd.DataFrame) -> None:
    missing = [column for column in OUTPUT_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"Comparison statistics output is missing required columns: {missing}")


def validate_no_group_vs_retrieval_rows(frame: pd.DataFrame) -> None:
    illegal = frame.loc[
        frame["comparison_family"].astype(str).str.contains("group_vs_retrieval", na=False)
        | frame["comparison_spec_id"].astype(str).str.contains("group_vs_retrieval", na=False)
    ]
    if not illegal.empty:
        raise ValueError("Group-vs-retrieval rows are not allowed by the frozen contract")


def apply_family_bh(unique_results: pd.DataFrame) -> pd.DataFrame:
    out = unique_results.copy()
    out["bh_q"] = np.nan
    out["significant_bool"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
    tested_mask = out["test_status"].eq("tested")
    for family, family_idx in out.loc[tested_mask].groupby("comparison_family", sort=False).groups.items():
        _ = family  # keep readable groupby pattern without relying on the key later
        q_values = bh_adjust(out.loc[family_idx, "raw_p"].astype("float64"))
        out.loc[family_idx, "bh_q"] = q_values
        out.loc[family_idx, "significant_bool"] = q_values.le(0.05).astype("boolean")
    return out


def finalize_unique_result(
    spec: ComparisonSpec,
    fixed_context: dict[str, object],
    group_a: str,
    group_b: str,
    values_a: pd.Series,
    values_b: pd.Series,
) -> dict[str, object]:
    n_units = int(len(values_a))
    notes_parts = [spec.notes_prefix] if spec.notes_prefix else []
    notes_parts.append("effect_direction reflects the sign of median(group_a - group_b); metric interpretation is deferred.")
    analysis_family = (
        str(fixed_context["analysis_family"])
        if "analysis_family" in fixed_context and pd.notna(fixed_context["analysis_family"])
        else str(spec.default_analysis_family)
    )
    if analysis_family == "None":
        raise ValueError(f"Missing analysis_family for comparison spec {spec.comparison_spec_id}")

    result = {
        "comparison_spec_id": spec.comparison_spec_id,
        "comparison_family": spec.comparison_family,
        "figure_id": spec.figure_id,
        "analysis_family": analysis_family,
        "metric_name": str(fixed_context["metric_name"]),
        "group_a": group_a,
        "group_b": group_b,
        "pairing_type": spec.pairing_type,
        "pairing_key_definition": spec.pairing_key_definition,
        "n_test_units": n_units,
        "raw_p": np.nan,
        "bh_q": np.nan,
        "significant_bool": pd.NA,
        "effect_direction": pd.NA,
        "median_delta": np.nan,
        "test_name": "wilcoxon_signed_rank" if spec.pairing_type == "paired" else "mann_whitney_u",
        "test_status": "not_tested_low_n" if n_units < LOW_N_THRESHOLD else "tested",
        "notes": " ".join(part for part in notes_parts if part),
        "n_units": n_units,
        "p": np.nan,
        "effect_size": np.nan,
        "primary_source_object": spec.primary_source_object,
        "dataset_scope": dataset_scope_label(fixed_context.get("dataset")),
        "representation_scope": (
            str(
                fixed_context.get(
                    "representation_scope",
                    fixed_context.get("representation_detail", fixed_context.get("representation", "ALL")),
                )
            )
            if pd.notna(
                fixed_context.get(
                    "representation_scope",
                    fixed_context.get("representation_detail", fixed_context.get("representation", "ALL")),
                )
            )
            else "ALL"
        ),
        "direction_scope": direction_scope_label(fixed_context.get("direction")),
        "question_ids": spec.question_id_resolver(fixed_context),
    }

    if n_units < LOW_N_THRESHOLD:
        return result

    deltas = values_a.to_numpy(dtype="float64") - values_b.to_numpy(dtype="float64")
    result["median_delta"] = float(np.median(deltas))
    result["effect_direction"] = effect_direction_from_delta(float(result["median_delta"]))
    result["effect_size"] = float(abs(result["median_delta"]))

    if spec.pairing_type == "paired":
        if np.allclose(deltas, 0.0, rtol=0.0, atol=0.0):
            result["raw_p"] = 1.0
            result["notes"] = (
                f"{result['notes']} all paired deltas were zero; raw_p was set to 1.0."
            ).strip()
        else:
            test_result = wilcoxon(
                values_a.to_numpy(dtype="float64"),
                values_b.to_numpy(dtype="float64"),
                alternative="two-sided",
                zero_method="wilcox",
                method="auto",
            )
            result["raw_p"] = float(test_result.pvalue)
    else:
        test_result = mannwhitneyu(
            values_a.to_numpy(dtype="float64"),
            values_b.to_numpy(dtype="float64"),
            alternative="two-sided",
        )
        result["raw_p"] = float(test_result.pvalue)
    result["p"] = result["raw_p"]

    return result


def run_paired_tests_from_long(
    *,
    long_frame: pd.DataFrame,
    spec: ComparisonSpec,
    fixed_columns: Sequence[str],
    unit_columns: Sequence[str],
    group_column: str,
    value_column: str,
    group_pairs: Sequence[tuple[str, str]],
    duplicate_label: str,
) -> pd.DataFrame:
    assert_unique_rows(
        long_frame,
        [*fixed_columns, *unit_columns, group_column],
        duplicate_label,
    )
    unique_rows: list[dict[str, object]] = []

    for group_a, group_b in group_pairs:
        subset = long_frame.loc[long_frame[group_column].isin([group_a, group_b])].copy()
        if subset.empty:
            continue
        wide = subset.pivot(
            index=[*fixed_columns, *unit_columns],
            columns=group_column,
            values=value_column,
        ).reset_index()
        if group_a not in wide.columns or group_b not in wide.columns:
            continue

        groupby_columns = list(fixed_columns)
        grouped = wide.groupby(groupby_columns, dropna=False, sort=True) if groupby_columns else [((), wide)]
        for group_key, strata in grouped:
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            fixed_context = {
                column: value
                for column, value in zip(groupby_columns, group_key)
            }
            paired = strata.loc[strata[[group_a, group_b]].notna().all(axis=1), [group_a, group_b]]
            unique_rows.append(
                finalize_unique_result(
                    spec=spec,
                    fixed_context=fixed_context,
                    group_a=group_a,
                    group_b=group_b,
                    values_a=paired[group_a] if not paired.empty else pd.Series(dtype="float64"),
                    values_b=paired[group_b] if not paired.empty else pd.Series(dtype="float64"),
                )
            )

    if not unique_rows:
        return pd.DataFrame(
            columns=[
                "comparison_spec_id",
                "comparison_family",
                "figure_id",
                "analysis_family",
                "metric_name",
                "group_a",
                "group_b",
                "pairing_type",
                "pairing_key_definition",
                "n_test_units",
                "raw_p",
                "bh_q",
                "significant_bool",
                "effect_direction",
                "median_delta",
                "test_name",
                "test_status",
                "notes",
                "primary_source_object",
                "dataset_scope",
                "direction_scope",
                "question_ids",
            ]
        )
    return pd.DataFrame(unique_rows)


def expand_question_rows(unique_results: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in unique_results.to_dict("records"):
        question_ids = tuple(row.pop("question_ids"))
        for question_id in question_ids:
            expanded = dict(row)
            expanded["question_id"] = question_id
            rows.append(expanded)
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    for column in OUTPUT_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    out = out[OUTPUT_COLUMNS].copy()
    out["significant_bool"] = out["significant_bool"].astype("boolean")
    return out


def build_task1_bridge_summary_long(path: Path, value_column: str) -> pd.DataFrame:
    frame = read_csv_required(
        path,
        {
            "dataset",
            "cell_line",
            "metric_family",
            "metric_name",
            "representation_detail",
            value_column,
        },
        "Task1 bridge summary",
    )
    frame = frame.rename(columns={value_column: "metric_value"}).copy()
    frame["analysis_family"] = frame["metric_family"].map(task1_metric_family_to_analysis_family)
    frame["direction"] = pd.NA
    out = frame[
        [
            "dataset",
            "cell_line",
            "analysis_family",
            "direction",
            "metric_name",
            "representation_detail",
            "metric_value",
        ]
    ].copy()
    out = out.rename(columns={"representation_detail": "group_label"})
    out = out.loc[out["metric_value"].notna()].copy()
    return out.sort_values(
        ["dataset", "cell_line", "analysis_family", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task1_bridge_detail_long(path: Path) -> pd.DataFrame:
    frame = read_csv_required(
        path,
        {
            "dataset",
            "cell_line",
            "target",
            "representation_detail",
            "metric_family",
            "metric_name",
            "internal_value",
            "cross_value",
        },
        "Task1 bridge detail",
    )
    rows: list[dict[str, object]] = []
    for group_label, source_column in [("internal", "internal_value"), ("cross", "cross_value")]:
        subset = frame[
            [
                "dataset",
                "cell_line",
                "target",
                "representation_detail",
                "metric_family",
                "metric_name",
                source_column,
            ]
        ].copy()
        subset = subset.rename(columns={source_column: "metric_value"})
        subset["group_label"] = group_label
        rows.append(subset)
    out = pd.concat(rows, ignore_index=True)
    out["analysis_family"] = out["metric_family"].map(task1_metric_family_to_analysis_family)
    out["direction"] = pd.NA
    out = out.loc[out["metric_value"].notna()].copy()
    return out[
        [
            "dataset",
            "cell_line",
            "target",
            "representation_detail",
            "analysis_family",
            "direction",
            "metric_name",
            "group_label",
            "metric_value",
        ]
    ].sort_values(
        ["dataset", "cell_line", "target", "representation_detail", "analysis_family", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task2_direction_long(path: Path) -> pd.DataFrame:
    frame = read_csv_required(
        path,
        {
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "mean_mrr_corrected",
            "mean_hit1_corrected",
            "mean_hit5_corrected",
            "mean_hit10_corrected",
        },
        "Figure 3 direction support summary",
    )
    long_frame = frame.melt(
        id_vars=["dataset", "cell_line", "direction", "representation"],
        value_vars=[f"mean_{metric}" for metric in TASK2_RETRIEVAL_METRICS],
        var_name="metric_name",
        value_name="metric_value",
    )
    long_frame["metric_name"] = long_frame["metric_name"].str.removeprefix("mean_")
    long_frame["analysis_family"] = "retrieval"
    long_frame = long_frame.loc[long_frame["metric_value"].notna()].copy()
    return long_frame.rename(columns={"representation": "representation_label"}).sort_values(
        ["dataset", "cell_line", "direction", "representation_label", "metric_name"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task2_performance_long(path: Path) -> pd.DataFrame:
    frame = read_csv_required(
        path,
        {
            "analysis_family",
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "metric_name",
            "metric_value",
        },
        "Figure 3 Task2 performance structure",
    )
    frame["direction"] = frame["direction"].fillna("")
    frame = frame.loc[frame["metric_value"].notna()].copy()
    return frame.rename(columns={"representation": "group_label"}).sort_values(
        ["analysis_family", "dataset", "cell_line", "direction", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task2_target_level_performance_long(
    *,
    task2_group_root: Path,
    task2_retrieval_per_query_path: Path,
) -> pd.DataFrame:
    long_frame = build_task2_target_support_long(task2_group_root, task2_retrieval_per_query_path.parent).copy()
    if long_frame.empty:
        return pd.DataFrame(
            columns=[
                "dataset",
                "cell_line",
                "target",
                "analysis_family",
                "direction",
                "metric_name",
                "group_label",
                "metric_value",
            ]
        )
    long_frame["direction"] = long_frame["direction"].fillna("")
    out = long_frame.rename(columns={"representation": "group_label"}).copy()
    return out[
        [
            "dataset",
            "cell_line",
            "target",
            "analysis_family",
            "direction",
            "metric_name",
            "group_label",
            "metric_value",
        ]
    ].sort_values(
        ["dataset", "cell_line", "target", "analysis_family", "direction", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task2_cell_line_pattern_long(path: Path) -> pd.DataFrame:
    frame = read_csv_required(
        path,
        {
            "dataset",
            "cell_line",
            "analysis_family",
            "direction",
            "representation",
            "metric_name",
            "mean_metric_value",
        },
        "Figure 3 Task2 cell-line pattern summary",
    )
    frame["direction"] = frame["direction"].fillna("")
    frame = frame.rename(columns={"mean_metric_value": "metric_value", "representation": "group_label"})
    frame = frame.loc[frame["metric_value"].notna()].copy()
    return frame.sort_values(
        ["dataset", "cell_line", "analysis_family", "direction", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task2_target_pattern_long(path: Path) -> pd.DataFrame:
    frame = read_csv_required(
        path,
        {
            "dataset",
            "target",
            "analysis_family",
            "direction",
            "representation",
            "metric_name",
            "mean_metric_value",
        },
        "Figure 3 Task2 target pattern summary",
    )
    frame["direction"] = frame["direction"].fillna("")
    frame = frame.rename(columns={"mean_metric_value": "metric_value", "representation": "group_label"})
    frame = frame.loc[frame["metric_value"].notna()].copy()
    return frame.sort_values(
        ["dataset", "target", "analysis_family", "direction", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task1_contextual_long(path: Path) -> pd.DataFrame:
    frame = read_csv_required(
        path,
        {
            "dataset",
            "cell_line",
            "target",
            "perturbation_type",
            "representation",
            "analysis_family",
            "metric_name",
            "metric_value",
        },
        "Figure 3 Task1 internal contextual support summary",
    )
    frame = frame.loc[frame["metric_value"].notna()].copy()
    frame = (
        frame.groupby(
            ["dataset", "cell_line", "target", "representation", "analysis_family", "metric_name"],
            dropna=False,
            sort=True,
            as_index=False,
        )["metric_value"]
        .mean()
    )
    frame["direction"] = pd.NA
    return frame[
        [
            "dataset",
            "cell_line",
            "target",
            "analysis_family",
            "direction",
            "metric_name",
            "representation",
            "metric_value",
        ]
    ].rename(columns={"representation": "group_label"}).sort_values(
        ["dataset", "cell_line", "target", "analysis_family", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_task2_k562_c2g_query_long(path: Path) -> pd.DataFrame:
    require_path(path, "Task2 retrieval per-query parquet")
    if pq is None:
        raise ModuleNotFoundError(
            "pyarrow is required to read Task2 retrieval parquet inputs for K562 C2G query-long support."
        )
    parquet_file = pq.ParquetFile(path)
    rows: list[pd.DataFrame] = []
    columns = ["dataset", "cell_line", "direction", "representation", "query_row_id", *TASK2_RETRIEVAL_METRICS]
    for batch in parquet_file.iter_batches(batch_size=200_000, columns=columns):
        frame = batch.to_pandas()
        frame = frame.loc[
            frame["dataset"].eq("scPerturb")
            & frame["cell_line"].eq("K562")
            & frame["direction"].eq("C2G")
        ].copy()
        if frame.empty:
            continue
        rows.append(frame)
    if not rows:
        return pd.DataFrame(columns=["analysis_family", "direction", "metric_name", "group_label", "metric_value", "query_row_id"])
    frame = pd.concat(rows, ignore_index=True)
    long_frame = frame.melt(
        id_vars=["dataset", "cell_line", "direction", "representation", "query_row_id"],
        value_vars=TASK2_RETRIEVAL_METRICS,
        var_name="metric_name",
        value_name="metric_value",
    )
    long_frame["analysis_family"] = "retrieval"
    long_frame = long_frame.loc[long_frame["metric_value"].notna()].copy()
    return long_frame.rename(columns={"representation": "group_label"}).sort_values(
        ["dataset", "cell_line", "direction", "query_row_id", "metric_name", "group_label"],
        kind="mergesort",
    ).reset_index(drop=True)


def spec1_questions(context: dict[str, object]) -> tuple[str, ...]:
    dataset = str(context["dataset"])
    if context["analysis_family"] == "group_concordance":
        if dataset == "LINCS":
            return ("C2_F2a_internal_group_lincs",)
        if dataset == "scPerturb":
            return ("C2_F2a_internal_group_scperturb_common_scope",)
    else:
        if dataset == "LINCS":
            return ("C2_F2b_internal_retrieval_lincs",)
        if dataset == "scPerturb":
            return ("C2_F2a_internal_retrieval_scperturb_common_scope",)
    raise ValueError(f"Unexpected dataset for Figure 2 internal common-scope questions: {dataset}")


def spec2_questions(context: dict[str, object]) -> tuple[str, ...]:
    if context["analysis_family"] == "group_concordance":
        return ("C2_F2e_cross_group_genetic",)
    return ("C2_F2f_cross_retrieval_genetic",)


def spec5_questions(_: dict[str, object]) -> tuple[str, ...]:
    return ("C2_F2i_internal_to_cross_degradation_summary",)


def spec6_questions(_: dict[str, object]) -> tuple[str, ...]:
    return ("C3_F3b_retrieval_c2g_common_scope", "C3_F3c_retrieval_g2c_common_scope")


def spec7_questions(context: dict[str, object]) -> tuple[str, ...]:
    if context["analysis_family"] == "group_concordance":
        return ("C3_F3a_group_common_scope",)
    if context["direction"] == "C2G":
        return ("C3_F3b_retrieval_c2g_common_scope",)
    return ("C3_F3c_retrieval_g2c_common_scope",)


def spec8_questions(_: dict[str, object]) -> tuple[str, ...]:
    return ("C3_F3h_cell_line_pattern_summary",)


def spec9_questions(_: dict[str, object]) -> tuple[str, ...]:
    return ("C3_F3i_target_pattern_summary",)


def spec10_questions(_: dict[str, object]) -> tuple[str, ...]:
    return ("C2_F3g_task1_internal_contextual_support",)


def spec12_questions(context: dict[str, object]) -> tuple[str, ...]:
    if context["analysis_family"] == "group_concordance":
        return ("C4_F3e_fm_group_k562_local",)
    return ("C4_F3f_fm_retrieval_k562_local",)


def spec13_questions(_: dict[str, object]) -> tuple[str, ...]:
    return ("C4_F3f_fm_retrieval_k562_local",)


def expand_fm_vs_base_tests(
    *,
    long_frame: pd.DataFrame,
    spec: ComparisonSpec,
    fixed_columns: Sequence[str],
    unit_columns: Sequence[str],
    duplicate_label: str,
) -> pd.DataFrame:
    fm_groups = sorted(value for value in long_frame["group_label"].dropna().astype(str).unique() if value.startswith("FM:") or value in {"geneformer", "scbert", "scfoundation", "scgpt", "state", "tahoe-x1", "uce"})
    rows: list[pd.DataFrame] = []
    for fm_group in fm_groups:
        for base_group in ("Gene", "Pathway"):
            subset = long_frame.loc[long_frame["group_label"].isin([fm_group, base_group])].copy()
            if subset.empty:
                continue
            rows.append(
                run_paired_tests_from_long(
                    long_frame=subset,
                    spec=spec,
                    fixed_columns=fixed_columns,
                    unit_columns=unit_columns,
                    group_column="group_label",
                    value_column="metric_value",
                    group_pairs=[(fm_group, base_group)],
                    duplicate_label=duplicate_label,
                )
            )
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def build_unique_results(
    *,
    task1_bridge_summary_path: Path,
    task1_bridge_detail_path: Path,
    task2_retrieval_per_query_path: Path,
    task2_group_root: Path,
    analysis_root: Path,
) -> pd.DataFrame:
    results: list[pd.DataFrame] = []

    spec1 = ComparisonSpec(
        comparison_spec_id="f2_task1_internal_common_representation",
        comparison_family="figure2_task1_internal_common_representation",
        figure_id="F2",
        primary_source_object=str(task1_bridge_detail_path),
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line + target",
        question_id_resolver=spec1_questions,
        notes_prefix="Task1 internal common-scope Gene vs Pathway from bridge detail internal triplets.",
    )
    bridge_internal = build_task1_bridge_detail_long(task1_bridge_detail_path)
    bridge_internal = bridge_internal.loc[bridge_internal["group_label"].eq("internal")].copy()
    bridge_internal["group_label"] = bridge_internal["representation_detail"]
    results.append(
        run_paired_tests_from_long(
            long_frame=bridge_internal,
            spec=spec1,
            fixed_columns=["dataset", "analysis_family", "metric_name"],
            unit_columns=["cell_line", "target"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="Task1 internal bridge detail comparison units",
        )
    )

    spec2 = ComparisonSpec(
        comparison_spec_id="f2_task1_cross_common_representation",
        comparison_family="figure2_task1_cross_common_representation",
        figure_id="F2",
        primary_source_object=str(task1_bridge_detail_path),
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line + target",
        question_id_resolver=spec2_questions,
        notes_prefix="Task1 cross common-scope Gene vs Pathway from bridge detail shared genetic triplets.",
    )
    bridge_cross = build_task1_bridge_detail_long(task1_bridge_detail_path)
    bridge_cross = bridge_cross.loc[bridge_cross["group_label"].eq("cross")].copy()
    bridge_cross["group_label"] = bridge_cross["representation_detail"]
    results.append(
        run_paired_tests_from_long(
            long_frame=bridge_cross,
            spec=spec2,
            fixed_columns=["dataset", "analysis_family", "metric_name"],
            unit_columns=["cell_line", "target"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="Task1 cross bridge detail comparison units",
        )
    )

    spec5 = ComparisonSpec(
        comparison_spec_id="f2_task1_internal_to_cross_degradation",
        comparison_family="figure2_task1_internal_to_cross_degradation",
        figure_id="F2",
        primary_source_object=str(task1_bridge_detail_path),
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line + target",
        question_id_resolver=spec5_questions,
        notes_prefix="Task1 internal vs cross degradation from bridge detail.",
    )
    bridge_detail_long = build_task1_bridge_detail_long(task1_bridge_detail_path)
    bridge_detail_long = apply_metric_transform(
        bridge_detail_long,
        metric_column="metric_name",
        value_column="metric_value",
    )
    results.append(
        run_paired_tests_from_long(
            long_frame=bridge_detail_long,
            spec=spec5,
            fixed_columns=["dataset", "representation_detail", "analysis_family", "metric_name"],
            unit_columns=["cell_line", "target"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("internal", "cross")],
            duplicate_label="Task1 internal-to-cross bridge detail units",
        )
    )

    spec6 = ComparisonSpec(
        comparison_spec_id="f3_task2_direction_retrieval",
        comparison_family="figure3_task2_direction_retrieval",
        figure_id="F3",
        primary_source_object=str(analysis_root / "figure3_task2_direction_support_summary.csv"),
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line",
        question_id_resolver=spec6_questions,
        default_analysis_family="retrieval",
        notes_prefix=(
            "Task2 C2G vs G2C retrieval from Figure 3 direction support summary. "
            "Supporting-only in current C2G-only manuscript phase."
        ),
    )
    task2_direction_long = build_task2_direction_long(analysis_root / "figure3_task2_direction_support_summary.csv")
    results.append(
        run_paired_tests_from_long(
            long_frame=task2_direction_long.rename(columns={"representation_label": "representation"}),
            spec=spec6,
            fixed_columns=["analysis_family", "representation", "metric_name"],
            unit_columns=["dataset", "cell_line"],
            group_column="direction",
            value_column="metric_value",
            group_pairs=[("C2G", "G2C")],
            duplicate_label="Task2 direction support units",
        )
    )

    spec7 = ComparisonSpec(
        comparison_spec_id="f3_task2_common_representation",
        comparison_family="figure3_task2_common_representation",
        figure_id="F3",
        primary_source_object=str(analysis_root / "figure3_task2_performance_structure.csv"),
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line",
        question_id_resolver=spec7_questions,
        notes_prefix="Task2 common-scope Gene vs Pathway from Figure 3 performance structure.",
    )
    task2_performance_long = build_task2_performance_long(analysis_root / "figure3_task2_performance_structure.csv")
    results.append(
        run_paired_tests_from_long(
            long_frame=task2_performance_long.loc[
                task2_performance_long["group_label"].isin(["Gene", "Pathway"])
            ].copy(),
            spec=spec7,
            fixed_columns=["dataset", "analysis_family", "direction", "metric_name"],
            unit_columns=["cell_line"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="Task2 performance-structure common representation units",
        )
    )

    spec7b = ComparisonSpec(
        comparison_spec_id="f3_task2_performance_target_level_representation",
        comparison_family="figure3_task2_performance_target_level_representation",
        figure_id="F3",
        primary_source_object=(
            f"{task2_group_root / 'task2_group_concordance.csv'}|{task2_retrieval_per_query_path}"
        ),
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line + target",
        question_id_resolver=spec7_questions,
        notes_prefix=(
            "Task2 target-level Gene vs Pathway from canonical Task2 target-anchored support rows. "
            "Both LINCS and scPerturb use matched (cell_line, target) units; retrieval remains C2G single-target only."
        ),
    )
    task2_target_level_performance = build_task2_target_level_performance_long(
        task2_group_root=task2_group_root,
        task2_retrieval_per_query_path=task2_retrieval_per_query_path,
    )
    task2_target_level_performance = apply_metric_transform(
        task2_target_level_performance,
        metric_column="metric_name",
        value_column="metric_value",
    )
    results.append(
        run_paired_tests_from_long(
            long_frame=task2_target_level_performance.loc[
                task2_target_level_performance["group_label"].isin(["Gene", "Pathway"])
            ].copy(),
            spec=spec7b,
            fixed_columns=["dataset", "analysis_family", "direction", "metric_name"],
            unit_columns=["cell_line", "target"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="Task2 target-level performance common representation units",
        )
    )

    spec8 = ComparisonSpec(
        comparison_spec_id="f3_task2_cell_line_pattern_representation",
        comparison_family="figure3_task2_cell_line_pattern_representation",
        figure_id="F3",
        primary_source_object=str(analysis_root / "figure3_task2_cell_line_pattern_summary.csv"),
        pairing_type="paired",
        pairing_key_definition="shared cell_line",
        question_id_resolver=spec8_questions,
        notes_prefix="Task2 cell-line pattern Gene vs Pathway from Figure 3 cell-line pattern summary.",
    )
    task2_cell_line_pattern = build_task2_cell_line_pattern_long(analysis_root / "figure3_task2_cell_line_pattern_summary.csv")
    results.append(
        run_paired_tests_from_long(
            long_frame=task2_cell_line_pattern.loc[
                task2_cell_line_pattern["group_label"].isin(["Gene", "Pathway"])
            ].copy(),
            spec=spec8,
            fixed_columns=["dataset", "analysis_family", "direction", "metric_name"],
            unit_columns=["cell_line"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="Task2 cell-line pattern units",
        )
    )

    spec9 = ComparisonSpec(
        comparison_spec_id="f3_task2_target_pattern_representation",
        comparison_family="figure3_task2_target_pattern_representation",
        figure_id="F3",
        primary_source_object=str(analysis_root / "figure3_task2_target_pattern_summary.csv"),
        pairing_type="paired",
        pairing_key_definition="shared target",
        question_id_resolver=spec9_questions,
        notes_prefix="Task2 target pattern Gene vs Pathway from Figure 3 target pattern summary.",
    )
    task2_target_pattern = build_task2_target_pattern_long(analysis_root / "figure3_task2_target_pattern_summary.csv")
    results.append(
        run_paired_tests_from_long(
            long_frame=task2_target_pattern.loc[
                task2_target_pattern["group_label"].isin(["Gene", "Pathway"])
            ].copy(),
            spec=spec9,
            fixed_columns=["dataset", "analysis_family", "direction", "metric_name"],
            unit_columns=["target"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="Task2 target pattern units",
        )
    )

    spec10 = ComparisonSpec(
        comparison_spec_id="f3_task1_contextual_common_representation",
        comparison_family="figure3_task1_contextual_common_representation",
        figure_id="F3",
        primary_source_object=str(analysis_root / "figure3_task1_internal_contextual_support_summary.csv"),
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line + target",
        question_id_resolver=spec10_questions,
        notes_prefix="Task1 internal contextual Gene vs Pathway from Figure 3 contextual support summary; repeated perturbation_type rows are collapsed by mean within the frozen manuscript unit.",
    )
    task1_contextual_long = build_task1_contextual_long(analysis_root / "figure3_task1_internal_contextual_support_summary.csv")
    results.append(
        run_paired_tests_from_long(
            long_frame=task1_contextual_long.loc[
                task1_contextual_long["group_label"].isin(["Gene", "Pathway"])
            ].copy(),
            spec=spec10,
            fixed_columns=["analysis_family", "metric_name"],
            unit_columns=["dataset", "cell_line", "target"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="Task1 contextual common representation units",
        )
    )

    spec12 = ComparisonSpec(
        comparison_spec_id="f3_task2_k562_fm_target_anchored",
        comparison_family="figure3_task2_k562_fm_target_anchored",
        figure_id="F3",
        primary_source_object=str(analysis_root / "figure3_task2_target_pattern_summary.csv"),
        pairing_type="paired",
        pairing_key_definition="shared target",
        question_id_resolver=spec12_questions,
        notes_prefix="Task2 K562-local FM vs Gene/Pathway target-anchored surface from Figure 3 target pattern summary.",
    )
    task2_k562_target_pattern = task2_target_pattern.loc[task2_target_pattern["dataset"].eq("scPerturb")].copy()
    results.append(
        expand_fm_vs_base_tests(
            long_frame=task2_k562_target_pattern,
            spec=spec12,
            fixed_columns=["dataset", "analysis_family", "direction", "metric_name"],
            unit_columns=["target"],
            duplicate_label="Task2 K562-local target-anchored FM units",
        )
    )

    spec13 = ComparisonSpec(
        comparison_spec_id="f3_task2_k562_fm_c2g_query_anchored",
        comparison_family="figure3_task2_k562_fm_c2g_query_anchored",
        figure_id="F3",
        primary_source_object=str(task2_retrieval_per_query_path),
        pairing_type="paired",
        pairing_key_definition="shared query_row_id",
        question_id_resolver=spec13_questions,
        default_analysis_family="retrieval",
        notes_prefix="Task2 K562-local FM vs Gene/Pathway C2G retrieval from per-query task2 retrieval parquet.",
    )
    task2_k562_query_long = build_task2_k562_c2g_query_long(task2_retrieval_per_query_path)
    results.append(
        expand_fm_vs_base_tests(
            long_frame=task2_k562_query_long,
            spec=spec13,
            fixed_columns=["analysis_family", "direction", "metric_name"],
            unit_columns=["dataset", "cell_line", "query_row_id"],
            duplicate_label="Task2 K562-local C2G FM query-anchored units",
        )
    )

    combined = pd.concat([frame for frame in results if not frame.empty], ignore_index=True)
    if combined.empty:
        raise ValueError("No manuscript comparison statistics rows were generated")
    return combined


def sort_output(frame: pd.DataFrame) -> pd.DataFrame:
    sort_columns = [
        "comparison_family",
        "comparison_spec_id",
        "question_id",
        "analysis_family",
        "metric_name",
        "group_a",
        "group_b",
        "dataset_scope",
        "direction_scope",
        "representation_scope",
    ]
    return frame.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)


def build_manifest_entry(output_path: Path) -> dict[str, object]:
    return annotate_canonical_objects(
        [
            {
                "object_name": "Manuscript comparison statistics",
                "manuscript_questions": [
                    "C2_F2a_internal_group_lincs",
                    "C2_F2b_internal_retrieval_lincs",
                    "C2_F2a_internal_group_scperturb_common_scope",
                    "C2_F2a_internal_retrieval_scperturb_common_scope",
                    "C2_F2e_cross_group_genetic",
                    "C2_F2f_cross_retrieval_genetic",
                    "C2_F2i_internal_to_cross_degradation_summary",
                    "C3_F3a_group_common_scope",
                    "C3_F3b_retrieval_c2g_common_scope",
                    "C3_F3c_retrieval_g2c_common_scope",
                    "C3_F3h_cell_line_pattern_summary",
                    "C3_F3i_target_pattern_summary",
                    "C2_F3g_task1_internal_contextual_support",
                    "C4_F3e_fm_group_k562_local",
                    "C4_F3f_fm_retrieval_k562_local",
                ],
                "source_inputs": [
                    str(DEFAULT_ANALYSIS_ROOT / "figure3_task2_direction_support_summary.csv"),
                    str(DEFAULT_ANALYSIS_ROOT / "figure3_task2_performance_structure.csv"),
                    str(DEFAULT_ANALYSIS_ROOT / "figure3_task2_cell_line_pattern_summary.csv"),
                    str(DEFAULT_ANALYSIS_ROOT / "figure3_task2_target_pattern_summary.csv"),
                    str(DEFAULT_ANALYSIS_ROOT / "figure3_task1_internal_contextual_support_summary.csv"),
                    str(DEFAULT_TASK2_GROUP_ROOT / "task2_group_concordance.csv"),
                    str(DEFAULT_TASK1_BRIDGE_SUMMARY),
                    str(DEFAULT_TASK1_BRIDGE_DETAIL),
                    str(DEFAULT_TASK2_RETRIEVAL_PER_QUERY),
                ],
                "output_files": [str(output_path)],
                "key_columns": [
                    "comparison_spec_id",
                    "question_id",
                    "analysis_family",
                    "metric_name",
                    "group_a",
                    "group_b",
                ],
                "representation_scope": "Frozen manuscript-facing Figure 2 and Figure 3 comparison surfaces only.",
                "included_logic": "Formal paired comparison tests with family-local BH correction across the frozen manuscript comparison specs.",
                "excluded_logic": "No group-vs-retrieval formal comparisons and no unsupported target-anchored C2G family.",
                "fm_scope": "Figure 3 Task2 scPerturb-local FM only; Task1-FM comparison surfaces are excluded in this phase.",
                "justification": "Comparative manuscript claims require formal statistical backing under the frozen current-phase contract.",
            }
        ]
    )[0]


def update_manifest_registration(manifest_path: Path, output_path: Path) -> None:
    require_path(manifest_path, "Framework analysis manifest")
    require_path(output_path, "Manuscript comparison statistics output")
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    canonical_objects = payload.get("canonical_objects", [])
    output_name = output_path.name
    canonical_objects = [
        entry
        for entry in canonical_objects
        if (
            str(output_path) not in {str(path) for path in entry.get("output_files", [])}
            and str(entry.get("canonical_filename", "")) != output_name
            and output_name not in {Path(str(path)).name for path in entry.get("output_files", [])}
        )
    ]
    canonical_objects.append(build_manifest_entry(output_path))
    payload["canonical_objects"] = canonical_objects
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    manuscript_active_root = args.manuscript_active_root.resolve()
    analysis_root = args.analysis_root.resolve()
    task2_retrieval_per_query_path = args.task2_retrieval_per_query_path.resolve()
    task2_group_root = args.task2_group_root.resolve()
    task1_bridge_summary_path = args.task1_bridge_summary_path.resolve()
    task1_bridge_detail_path = args.task1_bridge_detail_path.resolve()
    output_path = args.output_path.resolve()
    manifest_path = args.manifest_path.resolve()

    require_path(manuscript_active_root, "Manuscript active root")
    require_path(analysis_root, "Manuscript analysis root")
    require_path(task2_retrieval_per_query_path, "Task2 retrieval per-query parquet")
    require_path(task2_group_root, "Task2 group-concordance root")
    require_path(task1_bridge_summary_path, "Task1 bridge summary")
    require_path(task1_bridge_detail_path, "Task1 bridge detail")
    require_path(manifest_path, "Framework analysis manifest")

    unique_results = build_unique_results(
        task1_bridge_summary_path=task1_bridge_summary_path,
        task1_bridge_detail_path=task1_bridge_detail_path,
        task2_retrieval_per_query_path=task2_retrieval_per_query_path,
        task2_group_root=task2_group_root,
        analysis_root=analysis_root,
    )
    unique_results = apply_family_bh(unique_results)
    output = expand_question_rows(unique_results)
    validate_output_schema(output)
    validate_no_group_vs_retrieval_rows(output)
    output = sort_output(output)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    update_manifest_registration(manifest_path, output_path)

    print(output_path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
