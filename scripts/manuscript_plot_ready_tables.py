#!/usr/bin/env python3
"""
Build support-only plot-ready tables for the frozen manuscript plotting phase.

Status:
- support-only manuscript builder

Manuscript role:
- materializes plot-ready registries for Figure 1
- materializes panel-specific reductions for large downstream support tables
- keeps Python-side scope filtering and pre-aggregation separate from R plotting

Architecture:
- see `scripts/ARCHITECTURE.md` for script-family classification

This layer is downstream only. It does not create new canonical benchmark
evidence and it does not modify the frozen Figure 2/Figure 3 analysis objects.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, wilcoxon
try:
    import pyarrow.parquet as pq
except ModuleNotFoundError:
    pq = None

try:
    from path_policy import DEFAULT_ARCHIVE_ROOT, DEFAULT_RUNS_ROOT, DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT
except ModuleNotFoundError:
    from scripts.path_policy import DEFAULT_ARCHIVE_ROOT, DEFAULT_RUNS_ROOT, DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT

try:
    from manuscript_comparison_statistics import (
        build_task1_bridge_detail_long,
        build_task1_bridge_summary_long,
        build_task1_contextual_long,
        build_task2_cell_line_pattern_long,
        build_task2_performance_long,
        build_task2_target_pattern_long,
    )
except ModuleNotFoundError:
    from scripts.manuscript_comparison_statistics import (
        build_task1_bridge_detail_long,
        build_task1_bridge_summary_long,
        build_task1_contextual_long,
        build_task2_cell_line_pattern_long,
        build_task2_performance_long,
        build_task2_target_pattern_long,
    )


ROOT = Path(__file__).resolve().parents[1]
NAS_RUNS_ROOT = DEFAULT_RUNS_ROOT
NAS_ARCHIVE_ROOT = DEFAULT_ARCHIVE_ROOT
DEFAULT_MANUSCRIPT_ACTIVE_ROOT = NAS_RUNS_ROOT / "manuscript_active"
DEFAULT_MANUSCRIPT_SUPPORT_ROOT = NAS_RUNS_ROOT / "manuscript_support"
DEFAULT_MANUSCRIPT_HISTORY_ROOT = NAS_ARCHIVE_ROOT / "manuscript_history"
DEFAULT_ANALYSIS_ROOT = DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis"
DEFAULT_SUPPORT_ANALYSIS_ROOT = DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "analysis"
DEFAULT_PLOT_READY_ROOT = DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "plot_ready"
DEFAULT_CURRENT_S0_ROOT = DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "upstream_s0_current" / "s0_build_data_inventory"
DEFAULT_CURRENT_S4_ROOT = DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "upstream_s4_current" / "s4_task2_group_concordance_multisource"
DEFAULT_CURRENT_S5_ROOT = DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "upstream_s5_current" / "s5_task2_retrieval_multisource"

DEFAULT_TASK1_INVENTORY_PATH = DEFAULT_CURRENT_S0_ROOT / "task1_data_inventory_long.csv"
DEFAULT_TASK2_PAIRS_COVERAGE_PATH = DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT / "task2_pairs_coverage.csv"
DEFAULT_REPRESENTATION_REGISTRY_PATH = (
    DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT / "representation_availability_registry.csv"
)
DEFAULT_FRAMEWORK_MANIFEST_PATH = DEFAULT_ANALYSIS_ROOT / "framework_analysis_manifest.json"

DEFAULT_FIGURE2_CELL_HIGH_CONCORDANCE_PATH = (
    DEFAULT_SUPPORT_ANALYSIS_ROOT / "figure2_task1_cell_line_high_concordance_summary.csv"
)
DEFAULT_FIGURE2_TARGET_HIGH_CONCORDANCE_PATH = (
    DEFAULT_SUPPORT_ANALYSIS_ROOT / "figure2_task1_target_high_concordance_summary.csv"
)
DEFAULT_TASK1_BRIDGE_SUMMARY_PATH = (
    DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "group_bridge" / "task1_internal_vs_cross_group_bridge_summary.csv"
)
DEFAULT_TASK1_BRIDGE_DETAIL_PATH = (
    DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "group_bridge" / "task1_internal_vs_cross_group_bridge.csv"
)
DEFAULT_FIGURE3_SCOPE_PATH = DEFAULT_ANALYSIS_ROOT / "figure3_task2_scope_summary.csv"
DEFAULT_FIGURE3_TARGET_PATTERN_PATH = DEFAULT_ANALYSIS_ROOT / "figure3_task2_target_pattern_summary.csv"
DEFAULT_FIGURE3_PERFORMANCE_PATH = DEFAULT_ANALYSIS_ROOT / "figure3_task2_performance_structure.csv"
DEFAULT_FIGURE3_DIRECTION_SUPPORT_PATH = DEFAULT_ANALYSIS_ROOT / "figure3_task2_direction_support_summary.csv"
DEFAULT_FIGURE3_CONTEXTUAL_SUPPORT_PATH = (
    DEFAULT_ANALYSIS_ROOT / "figure3_task1_internal_contextual_support_summary.csv"
)
DEFAULT_COMPARISON_STATS_PATH = DEFAULT_ANALYSIS_ROOT / "manuscript_comparison_statistics.csv"
DEFAULT_TASK2_RETRIEVAL_PER_QUERY_PATH = DEFAULT_CURRENT_S5_ROOT / "task2_retrieval_per_query.parquet"
DEFAULT_TASK2_DRUG_META_PATH = (
    DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT / "scperturb_k562" / "Drug_meta.csv"
)
DEFAULT_TASK2_GROUP_ROOT = DEFAULT_CURRENT_S4_ROOT

MAIN_FIGURE1_RELATIVE_OUTPUTS = {
    "figure1_panel_1b_dataset_context_coverage.csv": Path("figure1") / "figure1_panel_1b_dataset_context_coverage.csv",
    "figure1_panel_1c_lawful_scope_matrix.csv": Path("figure1") / "figure1_panel_1c_lawful_scope_matrix.csv",
    "figure1_panel_1d_representation_modifier_availability.csv": (
        Path("figure1") / "figure1_panel_1d_representation_modifier_availability.csv"
    ),
    "figure1_panel_1e_result_object_map.csv": Path("figure1") / "figure1_panel_1e_result_object_map.csv",
    "extended_figure1_support_registry.csv": Path("extended") / "extended_figure1_support_registry.csv",
}

PLOT_READY_RELATIVE_OUTPUTS = {
    **MAIN_FIGURE1_RELATIVE_OUTPUTS,
    "figure2_panel_2a_task1_scope.csv": Path("figure2") / "figure2_panel_2a_task1_scope.csv",
    "figure2_panel_2b_internal_performance_overview.csv": (
        Path("figure2") / "figure2_panel_2b_internal_performance_overview.csv"
    ),
    "figure2_panel_2c_gene_vs_pathway_matched_units.csv": (
        Path("figure2") / "figure2_panel_2c_gene_vs_pathway_matched_units.csv"
    ),
    "figure2_panel_2d_internal_to_cross_degradation.csv": (
        Path("figure2") / "figure2_panel_2d_internal_to_cross_degradation.csv"
    ),
    "figure2_panel_2e_cell_line_high_concordance_summary.csv": (
        Path("figure2") / "figure2_panel_2e_cell_line_high_concordance_summary.csv"
    ),
    "figure2_panel_2f_target_high_concordance_summary.csv": (
        Path("figure2") / "figure2_panel_2f_target_high_concordance_summary.csv"
    ),
    "extended_figure4_cell_line_high_concordance_full.csv": (
        Path("extended") / "extended_figure4_cell_line_high_concordance_full.csv"
    ),
    "extended_figure4_target_high_concordance_full.csv": (
        Path("extended") / "extended_figure4_target_high_concordance_full.csv"
    ),
    "figure3_panel_3b_c2g_performance_overview.csv": (
        Path("figure3") / "figure3_panel_3b_c2g_performance_overview.csv"
    ),
    "figure3_panel_3c_cell_line_pattern.csv": Path("figure3") / "figure3_panel_3c_cell_line_pattern.csv",
    "figure3_panel_3d_target_pattern_summary.csv": (
        Path("figure3") / "figure3_panel_3d_target_pattern_summary.csv"
    ),
    "figure3_panel_3e_gene_vs_pathway_paired.csv": Path("figure3") / "figure3_panel_3e_gene_vs_pathway_paired.csv",
    "extended_figure6_target_pattern_full.csv": (
        Path("extended") / "extended_figure6_target_pattern_full.csv"
    ),
    "figure3_panel_3f_fm_local_tradeoff.csv": (
        Path("figure3") / "figure3_panel_3f_fm_local_tradeoff.csv"
    ),
    "figure3_panel_3j_task1_contextual_support_reference.csv": (
        Path("figure3") / "figure3_panel_3j_task1_contextual_support_reference.csv"
    ),
    "extended_figure8_contextual_support_full.csv": (
        Path("extended") / "extended_figure8_contextual_support_full.csv"
    ),
    "extended_figure9_support_vs_suitability.csv": (
        Path("extended") / "extended_figure9_support_vs_suitability.csv"
    ),
    "extended_figure9_target_exemplars.csv": (
        Path("extended") / "extended_figure9_target_exemplars.csv"
    ),
}

P0_PLOT_READY_OUTPUTS = [
    "figure1_panel_1b_dataset_context_coverage.csv",
    "figure1_panel_1c_lawful_scope_matrix.csv",
    "figure1_panel_1d_representation_modifier_availability.csv",
    "figure1_panel_1e_result_object_map.csv",
    "figure2_panel_2e_cell_line_high_concordance_summary.csv",
    "figure2_panel_2f_target_high_concordance_summary.csv",
    "figure3_panel_3d_target_pattern_summary.csv",
    "figure3_panel_3f_fm_local_tradeoff.csv",
    "figure3_panel_3j_task1_contextual_support_reference.csv",
]

SPATIAL_REVISION_PLOT_READY_OUTPUTS = [
    "figure2_panel_2a_task1_scope.csv",
    "figure2_panel_2b_internal_performance_overview.csv",
    "figure2_panel_2c_gene_vs_pathway_matched_units.csv",
    "figure2_panel_2d_internal_to_cross_degradation.csv",
    "figure2_panel_2e_cell_line_high_concordance_summary.csv",
    "figure2_panel_2f_target_high_concordance_summary.csv",
    "figure3_panel_3b_c2g_performance_overview.csv",
    "figure3_panel_3c_cell_line_pattern.csv",
    "figure3_panel_3d_target_pattern_summary.csv",
    "figure3_panel_3e_gene_vs_pathway_paired.csv",
    "figure3_panel_3f_fm_local_tradeoff.csv",
    "figure1_panel_1b_dataset_context_coverage.csv",
]

FIGURE2_FAMILY_COLUMNS = ["scope", "dataset_or_direction", "perturbation_type", "representation"]
HIGH_CONCORDANCE_COMMON_SCOPE_REPRESENTATIONS = {"Gene", "Pathway"}
HIGH_CONCORDANCE_UNDERPOWERED_THRESHOLD = 20
TASK1_HIGH_CONCORDANCE_SUCCESS_DEFINITION = "(mrr_corrected > 0) AND (hit10_corrected > 0)"
TASK2_COMMON_BASE_REPRESENTATIONS = {"Gene", "Pathway"}
TASK2_FROZEN_REPRESENTATION_CLASSES = ("Gene", "Pathway", "FM")
FM_LOCAL_DATASET = "scPerturb"
FM_LOCAL_CELL_LINE = "K562"
EXCLUDED_METRIC_NAMES = {"edist", "edist_biascorr", "mean_edist_biascorr"}
TASK2_RETRIEVAL_VALUE_COLUMNS = ("mrr_corrected", "hit1_corrected", "hit5_corrected", "hit10_corrected")
TASK2_TARGET_SUPPORT_COLUMNS = [
    "dataset",
    "cell_line",
    "target",
    "representation",
    "analysis_family",
    "direction",
    "metric_name",
    "metric_value",
    "n_queries",
    "N_gallery_mean",
    "m_pos_mean",
    "n_chem_instances_used",
    "n_gen_instances_used",
    "n_chem_sub",
    "n_gen_sub",
]

FM_REPRESENTATIONS = {"geneformer", "scbert", "scfoundation", "scgpt", "state", "tahoe-x1", "uce"}
FIGURE1_OBJECT_MAP_SUPPORT_OBJECTS = [
    {
        "object_name": "Task1 cell-line high-concordance summary",
        "object_status": "support_only",
        "figure_target": "F2",
        "panel_target": "2E|EF4A",
        "canonical_filename": "figure2_task1_cell_line_high_concordance_summary.csv",
        "source_class": "support_only_analysis",
        "disposition": "supporting",
        "notes": "Support-only enrichment summary; underpowered strata must be filtered before plotting.",
        "path": DEFAULT_FIGURE2_CELL_HIGH_CONCORDANCE_PATH,
    },
    {
        "object_name": "Task1 target high-concordance summary",
        "object_status": "support_only",
        "figure_target": "F2",
        "panel_target": "2F|EF4B",
        "canonical_filename": "figure2_task1_target_high_concordance_summary.csv",
        "source_class": "support_only_analysis",
        "disposition": "supporting",
        "notes": "Support-only enrichment summary; main panel must use pooled summary form rather than the raw target surface.",
        "path": DEFAULT_FIGURE2_TARGET_HIGH_CONCORDANCE_PATH,
    },
    {
        "object_name": "Task2 time modifier summary",
        "object_status": "support_only",
        "figure_target": "EF7",
        "panel_target": "EF7A",
        "canonical_filename": "figure3_task2_c2g_time_effect_summary.csv",
        "source_class": "support_only_analysis",
        "disposition": "supporting",
        "notes": "Supplementary C2G-only modifier result family.",
        "path": DEFAULT_ANALYSIS_ROOT / "figure3_task2_c2g_time_effect_summary.csv",
    },
    {
        "object_name": "Task2 dose modifier summary",
        "object_status": "support_only",
        "figure_target": "EF7",
        "panel_target": "EF7B",
        "canonical_filename": "figure3_task2_c2g_dose_effect_summary.csv",
        "source_class": "support_only_analysis",
        "disposition": "supporting",
        "notes": "Supplementary C2G-only modifier result family.",
        "path": DEFAULT_ANALYSIS_ROOT / "figure3_task2_c2g_dose_effect_summary.csv",
    },
    {
        "object_name": "Task2 target multiplicity modifier summary",
        "object_status": "support_only",
        "figure_target": "EF7",
        "panel_target": "EF7C",
        "canonical_filename": "figure3_task2_c2g_target_multiplicity_summary.csv",
        "source_class": "support_only_analysis",
        "disposition": "supporting",
        "notes": "Supplementary C2G-only modifier result family.",
        "path": DEFAULT_ANALYSIS_ROOT / "figure3_task2_c2g_target_multiplicity_summary.csv",
    },
    {
        "object_name": "Task2 modifier comparison statistics",
        "object_status": "support_only",
        "figure_target": "EF7",
        "panel_target": "EF7D",
        "canonical_filename": "figure3_task2_c2g_modifier_comparison_statistics.csv",
        "source_class": "support_only_analysis",
        "disposition": "supporting",
        "notes": "Supplementary BH-FDR comparison layer across modifier families.",
        "path": DEFAULT_ANALYSIS_ROOT / "figure3_task2_c2g_modifier_comparison_statistics.csv",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build support-only plot-ready manuscript tables.")
    parser.add_argument("--plot-ready-root", type=Path, default=DEFAULT_PLOT_READY_ROOT)
    parser.add_argument("--analysis-root", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    parser.add_argument("--task1-inventory-path", type=Path, default=DEFAULT_TASK1_INVENTORY_PATH)
    parser.add_argument("--task2-pairs-coverage-path", type=Path, default=DEFAULT_TASK2_PAIRS_COVERAGE_PATH)
    parser.add_argument("--representation-registry-path", type=Path, default=DEFAULT_REPRESENTATION_REGISTRY_PATH)
    parser.add_argument("--framework-manifest-path", type=Path, default=DEFAULT_FRAMEWORK_MANIFEST_PATH)
    parser.add_argument(
        "--figure2-cell-high-concordance-path",
        type=Path,
        default=DEFAULT_FIGURE2_CELL_HIGH_CONCORDANCE_PATH,
    )
    parser.add_argument(
        "--figure2-target-high-concordance-path",
        type=Path,
        default=DEFAULT_FIGURE2_TARGET_HIGH_CONCORDANCE_PATH,
    )
    parser.add_argument("--task1-bridge-summary-path", type=Path, default=DEFAULT_TASK1_BRIDGE_SUMMARY_PATH)
    parser.add_argument("--task1-bridge-detail-path", type=Path, default=DEFAULT_TASK1_BRIDGE_DETAIL_PATH)
    parser.add_argument("--figure3-scope-path", type=Path, default=DEFAULT_FIGURE3_SCOPE_PATH)
    parser.add_argument("--figure3-target-pattern-path", type=Path, default=DEFAULT_FIGURE3_TARGET_PATTERN_PATH)
    parser.add_argument("--figure3-performance-path", type=Path, default=DEFAULT_FIGURE3_PERFORMANCE_PATH)
    parser.add_argument(
        "--figure3-direction-support-path",
        type=Path,
        default=DEFAULT_FIGURE3_DIRECTION_SUPPORT_PATH,
    )
    parser.add_argument("--task2-drug-meta-path", type=Path, default=DEFAULT_TASK2_DRUG_META_PATH)
    parser.add_argument(
        "--figure3-contextual-support-path",
        type=Path,
        default=DEFAULT_FIGURE3_CONTEXTUAL_SUPPORT_PATH,
    )
    parser.add_argument("--comparison-stats-path", type=Path, default=DEFAULT_COMPARISON_STATS_PATH)
    parser.add_argument("--task2-retrieval-per-query-path", type=Path, default=DEFAULT_TASK2_RETRIEVAL_PER_QUERY_PATH)
    parser.add_argument("--task2-group-root", type=Path, default=DEFAULT_TASK2_GROUP_ROOT)
    parser.add_argument(
        "--build-target",
        default="all",
        choices=["all", "p0", "spatial_revision", *sorted(PLOT_READY_RELATIVE_OUTPUTS.keys())],
        help="Build one named plot-ready output, the frozen P0 subset, the staged spatial subset, or every declared output.",
    )
    return parser.parse_args()


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")


def read_csv_required(path: Path, required_columns: set[str], label: str) -> pd.DataFrame:
    require_path(path, label)
    frame = pd.read_csv(path)
    missing = sorted(required_columns - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")
    return frame


def write_csv(frame: pd.DataFrame, output_path: Path, expected_columns: Sequence[str]) -> Path:
    missing = [column for column in expected_columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{output_path.name} is missing required columns: {missing}")
    out = frame.reindex(columns=list(expected_columns)).copy()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    return output_path


def aggregate_task2_c2g_target_means(retrieval_per_query_path: Path) -> pd.DataFrame:
    require_path(retrieval_per_query_path, "Task2 retrieval per-query parquet")
    if pq is None:
        raise ModuleNotFoundError("pyarrow is required to read Task2 retrieval parquet inputs.")

    parquet_file = pq.ParquetFile(retrieval_per_query_path)
    key_columns = ["dataset", "cell_line", "target", "representation"]
    tracked_columns = [*TASK2_RETRIEVAL_VALUE_COLUMNS, "N_gallery", "m_pos"]
    aggregate: dict[tuple[object, ...], dict[str, float | int]] = {}
    columns = [
        "dataset",
        "cell_line",
        "direction",
        "representation",
        "query_target_tokens",
        "query_n_targets",
        "N_gallery",
        "m_pos",
        *TASK2_RETRIEVAL_VALUE_COLUMNS,
    ]

    for batch in parquet_file.iter_batches(batch_size=200_000, columns=columns):
        frame = batch.to_pandas()
        frame = frame.loc[frame["direction"].eq("C2G") & frame["query_n_targets"].eq(1)].copy()
        if frame.empty:
            continue

        frame["target"] = frame["query_target_tokens"].astype(str).str.split(";").str[0]
        frame = frame.loc[frame["target"].ne("") & frame["target"].ne("NA")].copy()
        if frame.empty:
            continue

        for column in tracked_columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

        grouped = frame.groupby(key_columns, dropna=False, observed=True, sort=False)
        summary = grouped.size().rename("n_queries").to_frame()
        for column in tracked_columns:
            summary[f"{column}_sum"] = grouped[column].sum(min_count=1)
            summary[f"{column}_count"] = grouped[column].count()
        summary = summary.reset_index()

        for row in summary.itertuples(index=False):
            key = tuple(getattr(row, column) for column in key_columns)
            state = aggregate.setdefault(key, {"n_queries": 0})
            state["n_queries"] += int(row.n_queries)
            for column in tracked_columns:
                state.setdefault(f"{column}_sum", 0.0)
                state.setdefault(f"{column}_count", 0)
                sum_value = getattr(row, f"{column}_sum")
                count_value = getattr(row, f"{column}_count")
                if pd.notna(sum_value):
                    state[f"{column}_sum"] += float(sum_value)
                state[f"{column}_count"] += int(count_value)

    if not aggregate:
        return pd.DataFrame(
            columns=[
                *key_columns,
                "n_queries",
                *TASK2_RETRIEVAL_VALUE_COLUMNS,
                *(f"{column}_n_valid" for column in tracked_columns),
                "N_gallery",
                "m_pos",
            ]
        )

    rows: list[dict[str, object]] = []
    for key, state in aggregate.items():
        row = dict(zip(key_columns, key, strict=True))
        row["n_queries"] = int(state["n_queries"])
        for column in tracked_columns:
            count = int(state[f"{column}_count"])
            row[column] = float(state[f"{column}_sum"]) / count if count else pd.NA
            row[f"{column}_n_valid"] = count
        rows.append(row)

    return pd.DataFrame(rows).sort_values(key_columns, kind="mergesort").reset_index(drop=True)


def build_task2_target_support_long(task2_group_root: Path, task2_retrieval_root: Path) -> pd.DataFrame:
    group = read_csv_required(
        task2_group_root / "task2_group_concordance.csv",
        {
            "dataset",
            "cell_line",
            "target_token",
            "representation",
            "n_chem_instances_used",
            "n_gen_instances_used",
            "n_chem_sub",
            "n_gen_sub",
            "cosine_centroid",
            "pcc_centroid",
            "edist_biascorr",
            "cosine_valid_bool",
            "pcc_valid_bool",
            "edist_valid_bool",
        },
        "Task2 group concordance",
    ).rename(columns={"target_token": "target"})

    group_rows: list[dict[str, object]] = []
    metric_validity = {
        "cosine_centroid": "cosine_valid_bool",
        "pcc_centroid": "pcc_valid_bool",
        "edist_biascorr": "edist_valid_bool",
    }
    for row in group.itertuples(index=False):
        for metric_name, valid_column in metric_validity.items():
            if not bool(getattr(row, valid_column)):
                continue
            group_rows.append(
                {
                    "dataset": row.dataset,
                    "cell_line": row.cell_line,
                    "target": row.target,
                    "representation": row.representation,
                    "analysis_family": "group_concordance",
                    "direction": pd.NA,
                    "metric_name": metric_name,
                    "metric_value": getattr(row, metric_name),
                    "n_queries": pd.NA,
                    "N_gallery_mean": pd.NA,
                    "m_pos_mean": pd.NA,
                    "n_chem_instances_used": row.n_chem_instances_used,
                    "n_gen_instances_used": row.n_gen_instances_used,
                    "n_chem_sub": row.n_chem_sub,
                    "n_gen_sub": row.n_gen_sub,
                }
            )
    group_frame = pd.DataFrame(group_rows, columns=TASK2_TARGET_SUPPORT_COLUMNS)

    retrieval_target = aggregate_task2_c2g_target_means(task2_retrieval_root / "task2_retrieval_per_query.parquet")
    if retrieval_target.empty:
        return group_frame.sort_values(
            ["dataset", "cell_line", "target", "analysis_family", "direction", "representation", "metric_name"],
            kind="mergesort",
            na_position="last",
        ).reset_index(drop=True)

    retrieval_frame = retrieval_target.melt(
        id_vars=["dataset", "cell_line", "target", "representation", "n_queries", "N_gallery", "m_pos"],
        value_vars=list(TASK2_RETRIEVAL_VALUE_COLUMNS),
        var_name="metric_name",
        value_name="metric_value",
    )
    retrieval_frame["analysis_family"] = "retrieval"
    retrieval_frame["direction"] = "C2G"
    retrieval_frame["N_gallery_mean"] = retrieval_frame["N_gallery"]
    retrieval_frame["m_pos_mean"] = retrieval_frame["m_pos"]
    retrieval_frame["n_chem_instances_used"] = pd.NA
    retrieval_frame["n_gen_instances_used"] = pd.NA
    retrieval_frame["n_chem_sub"] = pd.NA
    retrieval_frame["n_gen_sub"] = pd.NA
    retrieval_frame = retrieval_frame[TASK2_TARGET_SUPPORT_COLUMNS]

    combined = pd.concat([group_frame, retrieval_frame], ignore_index=True, sort=False)
    return combined.sort_values(
        ["dataset", "cell_line", "target", "analysis_family", "direction", "representation", "metric_name"],
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)


def stringify_series(values: Iterable[object]) -> str:
    cleaned = [str(value) for value in values if pd.notna(value) and str(value)]
    return "|".join(sorted(dict.fromkeys(cleaned)))


def representation_class(value: object) -> str:
    text = str(value)
    if text in TASK2_COMMON_BASE_REPRESENTATIONS:
        return text
    return "FM"


def collapsed_status_label(values: Iterable[object]) -> str:
    series = pd.Series(list(values)).dropna().astype(str)
    if series.empty:
        return "not_applicable_scope"
    if series.eq("materialized").any():
        return "materialized"
    if series.eq("available").any():
        return "available"
    return str(series.iloc[0])


def status_available_numerator(value: object) -> int:
    return int(str(value) in {"available", "materialized"})


def task1_scope_coverage_denominator(analysis_family: object, dataset_or_direction: object) -> int:
    if str(analysis_family) == "cross_eligibility":
        return 1
    if str(dataset_or_direction) == "scPerturb":
        return 3
    return 2


def metric_is_excluded(metric_name: object) -> bool:
    return str(metric_name).strip() in EXCLUDED_METRIC_NAMES


def filter_excluded_metrics(frame: pd.DataFrame, metric_column: str = "metric_name") -> pd.DataFrame:
    if metric_column not in frame.columns:
        return frame.copy()
    return frame.loc[~frame[metric_column].map(metric_is_excluded)].copy()


SIGNED_LOG_METRICS = {"edist", "edist_biascorr", "mean_edist_biascorr"}


def signed_log10_1p(values: pd.Series | Sequence[float]) -> pd.Series:
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").astype("float64")
    return np.sign(numeric) * np.log10(1.0 + np.abs(numeric))


def transform_metric_values(metric_name: object, values: pd.Series) -> pd.Series:
    if str(metric_name) in SIGNED_LOG_METRICS:
        return signed_log10_1p(values)
    return pd.to_numeric(values, errors="coerce")


def apply_metric_transform(frame: pd.DataFrame, *, metric_column: str, value_column: str) -> pd.DataFrame:
    out = frame.copy()
    transformed = pd.Series(np.nan, index=out.index, dtype="float64")
    transform_label = pd.Series("identity", index=out.index, dtype="object")
    for metric_name, idx in out.groupby(metric_column, dropna=False, sort=False).groups.items():
        transformed.loc[idx] = transform_metric_values(metric_name, out.loc[idx, value_column]).astype("float64")
        if str(metric_name) in SIGNED_LOG_METRICS:
            transform_label.loc[idx] = "signed_log10_1p_abs"
    out[value_column] = transformed
    out["value_transform"] = transform_label
    return out


def split_target_tokens(value: object) -> list[str]:
    text = str(value).strip()
    if text in {"", "NA"}:
        return []
    parts: list[str] = []
    for chunk in text.split("|"):
        parts.extend(piece.strip() for piece in chunk.split("_"))
    return sorted({part for part in parts if part and part != "NA"})


def target_family_label(value: object) -> str:
    tokens = split_target_tokens(value)
    if not tokens:
        return str(value).strip()
    if all(token.startswith("HDAC") for token in tokens):
        return "HDAC"
    if all(token.startswith("PRKC") for token in tokens):
        return "PRKC"
    if all(token.startswith("KCNJ") for token in tokens):
        return "KCNJ"
    if all(token == "ABCC8" or token.startswith("KCNJ") for token in tokens):
        return "KATP-Channel"
    if len(tokens) == 1:
        match = re.match(r"^([A-Z]+)", tokens[0])
        if match and len(match.group(1)) >= 3:
            return match.group(1)
        return tokens[0]
    return target_display_label("|".join(tokens))


def target_display_label(value: object) -> str:
    text = str(value).strip()
    if not text or text == "NA":
        return text
    tokens = [token.strip() for token in text.split("|") if token.strip()]
    if len(tokens) <= 1:
        return text
    if all(token.startswith("HDAC") for token in tokens):
        return "HDAC-Family"
    if all(token.startswith("PRKC") for token in tokens):
        return "PKC-Family"
    if all(token.startswith("KCNJ") for token in tokens):
        return "KCNJ-Family"
    if all(token == "ABCC8" or token.startswith("KCNJ") for token in tokens):
        return "KATP-Channel"
    compact = "|".join(tokens)
    if len(compact) <= 32 and len(tokens) <= 3:
        return compact
    return f"{tokens[0]} +{len(tokens) - 1} targets"


def canonicalize_target_token_label(value: object) -> str:
    text = str(value).strip()
    if text in {"", "NA"}:
        return text
    parts: list[str] = []
    for chunk in text.split("|"):
        parts.extend(piece.strip() for piece in chunk.split("_"))
    cleaned = sorted({part for part in parts if part and part != "NA"})
    if not cleaned:
        return text
    if len(cleaned) == 1:
        return cleaned[0]
    return "|".join(cleaned)


def entity_display_label(entity_column: str, value: object) -> str:
    if entity_column == "target_token":
        return target_display_label(value)
    return str(value)


def entity_display_exclusion(entity_column: str, value: object) -> bool:
    if entity_column != "cell_line":
        return False
    return str(value).strip().upper() == "INTRAHEPATIC CHOLANGIOCYTE ORGANOIDS"


def fm_local_bool(value: object) -> bool:
    return representation_class(value) == "FM"


def max_or_na(values: Iterable[object]) -> float | pd._libs.missing.NAType:
    numeric = pd.to_numeric(pd.Series(list(values)), errors="coerce")
    if numeric.notna().sum() == 0:
        return pd.NA
    return float(numeric.max())


def any_true_bool(values: Iterable[object]) -> bool:
    series = pd.Series(list(values), dtype="boolean")
    return bool(series.fillna(False).any())


def assert_unique_rows(frame: pd.DataFrame, key_columns: Sequence[str], label: str) -> None:
    if frame.empty:
        return
    dup_mask = frame.duplicated(subset=list(key_columns), keep=False)
    if not dup_mask.any():
        return
    examples = frame.loc[dup_mask, list(key_columns)].head(10).to_dict("records")
    raise ValueError(f"{label} has duplicate rows on key={list(key_columns)}; examples={examples}")


def dense_rank_desc(frame: pd.DataFrame, group_columns: Sequence[str], value_column: str, output_column: str) -> pd.DataFrame:
    out = frame.copy()
    out[output_column] = (
        out.groupby(list(group_columns), dropna=False)[value_column]
        .rank(method="first", ascending=False)
        .astype("Int64")
    )
    return out


def percentile_from_rank(frame: pd.DataFrame, group_columns: Sequence[str], rank_column: str, output_column: str) -> pd.DataFrame:
    out = frame.copy()
    group_sizes = out.groupby(list(group_columns), dropna=False)[rank_column].transform("max").astype("float64")
    out[output_column] = np.where(
        group_sizes.gt(0.0),
        1.0 - ((out[rank_column].astype("float64") - 1.0) / group_sizes),
        np.nan,
    )
    return out


def percent_rank_within(frame: pd.DataFrame, group_columns: Sequence[str], value_column: str, output_column: str) -> pd.DataFrame:
    out = frame.copy()
    ranked = out.groupby(list(group_columns), dropna=False)[value_column].rank(
        method="average",
        ascending=True,
        na_option="keep",
    )
    valid_n = out.groupby(list(group_columns), dropna=False)[value_column].transform(lambda values: values.notna().sum())
    out[output_column] = np.where(
        valid_n.gt(0),
        ranked.astype("float64") / valid_n.astype("float64"),
        np.nan,
    )
    return out


def display_band_from_percentile(percentile_value: float) -> str:
    if pd.isna(percentile_value):
        return "unranked"
    if percentile_value >= 0.90:
        return "top_10pct"
    if percentile_value >= 0.75:
        return "top_25pct"
    if percentile_value >= 0.25:
        return "middle_50pct"
    return "bottom_25pct"


def pooled_odds_ratio(
    n_positive: int,
    n_negative: int,
    n_background_positive: int,
    n_background_negative: int,
) -> float | pd._libs.missing.NAType:
    if (n_background_positive + n_background_negative) == 0:
        return pd.NA
    odds_ratio, _ = fisher_exact(
        [[n_positive, n_negative], [n_background_positive, n_background_negative]]
    )
    return float(odds_ratio)


def build_figure1_coverage(
    *,
    task1_inventory_path: Path,
    task2_pairs_coverage_path: Path,
    output_path: Path,
) -> Path:
    task1_inventory = read_csv_required(
        task1_inventory_path,
        {"dataset", "perturbation_type", "cell_line", "target_token", "n_instances"},
        "Task1 inventory",
    ).copy()
    task1_inventory["target_token"] = task1_inventory["target_token"].map(canonicalize_target_token_label)
    task1 = (
        task1_inventory.groupby(
            ["dataset", "perturbation_type", "cell_line", "target_token"],
            dropna=False,
            sort=True,
        )
        .agg(
            leaf_weight=("n_instances", "sum"),
            n_instances=("n_instances", "sum"),
        )
        .reset_index()
    )
    task1["task"] = "Task1"
    task1["task_order"] = 1
    task1["dataset_order"] = task1["dataset"].map({"LINCS": 1, "scPerturb": 2}).fillna(99).astype(int)
    task1["eligible_target_bool"] = True
    task1["n_chem_instances"] = pd.NA
    task1["n_gen_instances"] = pd.NA
    task1["source_table"] = str(task1_inventory_path)
    task1["scope_note"] = "Task1 S0 ecosystem/context inventory."

    task2_pairs = read_csv_required(
        task2_pairs_coverage_path,
        {"dataset", "cell_line", "target_token", "n_chem_instances", "n_gen_instances", "is_eligible_bool"},
        "Task2 pairs coverage",
    ).copy()
    task2_pairs["target_token"] = task2_pairs["target_token"].map(canonicalize_target_token_label)
    task2_pairs["eligible_target_bool"] = task2_pairs["is_eligible_bool"].astype(str).str.lower().eq("true")
    task2_pairs = task2_pairs.loc[task2_pairs["eligible_target_bool"]].copy()
    task2_rows: list[pd.DataFrame] = []
    for perturbation_type, source_column in [("Chemical", "n_chem_instances"), ("Genetic", "n_gen_instances")]:
        subset = task2_pairs[["dataset", "cell_line", "target_token", "n_chem_instances", "n_gen_instances"]].copy()
        subset["leaf_weight"] = pd.to_numeric(task2_pairs[source_column], errors="coerce").fillna(0)
        subset = subset.loc[subset["leaf_weight"].gt(0)].copy()
        if subset.empty:
            continue
        grouped = (
            subset.groupby(["dataset", "cell_line", "target_token"], dropna=False, sort=True)
            .agg(
                leaf_weight=("leaf_weight", "sum"),
                n_chem_instances=("n_chem_instances", "sum"),
                n_gen_instances=("n_gen_instances", "sum"),
            )
            .reset_index()
        )
        grouped["n_instances"] = grouped["leaf_weight"]
        grouped["task"] = "Task2"
        grouped["task_order"] = 2
        grouped["dataset_order"] = grouped["dataset"].map({"LINCS": 1, "scPerturb": 2}).fillna(99).astype(int)
        grouped["perturbation_type"] = perturbation_type
        grouped["eligible_target_bool"] = True
        grouped["source_table"] = str(task2_pairs_coverage_path)
        grouped["scope_note"] = (
            "Task2 corrected snapshot target coverage; paired rows are split into Chemical and Genetic ring contributions."
        )
        task2_rows.append(grouped)
    task2 = pd.concat(task2_rows, ignore_index=True, sort=False) if task2_rows else pd.DataFrame(columns=task1.columns)

    combined = pd.concat([task1, task2], ignore_index=True, sort=False)
    combined["target_display_label"] = combined["target_token"].map(target_display_label)
    slice_totals = (
        combined.groupby(["dataset", "task", "perturbation_type"], dropna=False, sort=True)
        .agg(
            cell_line_count=("cell_line", "nunique"),
            target_count=("target_token", "nunique"),
            instance_count=("leaf_weight", "sum"),
            eligible_target_count=("target_token", "nunique"),
        )
        .reset_index()
    )
    dataset_totals = (
        combined.groupby(["dataset"], dropna=False, sort=True)
        .agg(
            dataset_cell_line_count=("cell_line", "nunique"),
            dataset_target_count=("target_token", "nunique"),
            dataset_leaf_weight=("leaf_weight", "sum"),
        )
        .reset_index()
    )
    task_totals = (
        combined.groupby(["dataset", "task"], dropna=False, sort=True)
        .agg(task_weight=("leaf_weight", "sum"))
        .reset_index()
    )
    perturbation_totals = (
        combined.groupby(["dataset", "task", "perturbation_type"], dropna=False, sort=True)
        .agg(perturbation_weight=("leaf_weight", "sum"))
        .reset_index()
        .sort_values(
            ["dataset", "task", "perturbation_type"],
            ascending=[True, True, True],
            kind="mergesort",
        )
    )
    perturbation_totals["perturbation_rank_dataset"] = (
        perturbation_totals.groupby(["dataset", "task"], dropna=False).cumcount() + 1
    ).astype(int)
    cell_line_totals = (
        combined.groupby(["dataset", "task", "perturbation_type", "cell_line"], dropna=False, sort=True)
        .agg(
            cell_line_weight=("leaf_weight", "sum"),
            cell_line_target_count=("target_token", "nunique"),
            cell_line_instance_count=("n_instances", "sum"),
            cell_line_n_chem_instances=("n_chem_instances", "sum"),
            cell_line_n_gen_instances=("n_gen_instances", "sum"),
            source_table=("source_table", stringify_series),
            scope_note=("scope_note", stringify_series),
        )
        .reset_index()
        .sort_values(
            ["dataset", "task", "perturbation_type", "cell_line_weight", "cell_line"],
            ascending=[True, True, True, False, True],
            kind="mergesort",
        )
    )
    parent_totals = cell_line_totals.groupby(
        ["dataset", "task", "perturbation_type"],
        dropna=False,
        sort=False,
    )["cell_line_weight"].transform("sum")
    cell_line_totals["share_within_parent"] = np.where(
        parent_totals.gt(0),
        pd.to_numeric(cell_line_totals["cell_line_weight"], errors="coerce") / parent_totals,
        np.nan,
    )
    cell_line_totals["cell_line_rank_within_parent"] = (
        cell_line_totals.groupby(["dataset", "task", "perturbation_type"], dropna=False).cumcount() + 1
    ).astype(int)
    cell_line_totals["context_keep_bool"] = (
        cell_line_totals["cell_line_rank_within_parent"].le(6)
        | cell_line_totals["share_within_parent"].ge(0.04)
    )
    cell_line_totals["context_base_label"] = cell_line_totals["cell_line"].astype(str)
    cell_line_totals["context_label"] = np.where(
        cell_line_totals["context_keep_bool"],
        cell_line_totals["context_base_label"],
        "Others",
    )
    context_totals = (
        cell_line_totals.groupby(
            ["dataset", "task", "perturbation_type", "context_label"],
            dropna=False,
            sort=True,
        )
        .agg(
            context_weight=("cell_line_weight", "sum"),
            context_cell_line_count=("cell_line", "nunique"),
            context_target_count=("cell_line_target_count", "sum"),
            n_instances=("cell_line_instance_count", "sum"),
            n_chem_instances=("cell_line_n_chem_instances", "sum"),
            n_gen_instances=("cell_line_n_gen_instances", "sum"),
            source_table=("source_table", stringify_series),
            scope_note=("scope_note", stringify_series),
        )
        .reset_index()
        .sort_values(
            ["dataset", "task", "perturbation_type", "context_weight", "context_label"],
            ascending=[True, True, True, False, True],
            kind="mergesort",
        )
    )
    context_totals["context_rank_within_parent"] = (
        context_totals.groupby(["dataset", "task", "perturbation_type"], dropna=False).cumcount() + 1
    ).astype(int)
    context_totals["context_folded_bool"] = context_totals["context_label"].eq("Others")
    context_totals["context_label_keep_bool"] = True
    context_totals["context_note"] = np.where(
        context_totals["context_folded_bool"],
        context_totals["context_cell_line_count"].map(lambda value: f"{int(value)} cell lines"),
        context_totals["context_target_count"].map(lambda value: f"{int(value)} targets"),
    )
    context_totals["context_summary_label"] = np.where(
        context_totals["context_folded_bool"],
        context_totals["context_label"].astype(str),
        context_totals["context_label"].astype(str) + "\n" + context_totals["context_note"].astype(str),
    )
    target_totals = (
        combined.groupby(["dataset", "target_token"], dropna=False, sort=True)
        .agg(target_weight=("leaf_weight", "sum"))
        .reset_index()
        .sort_values(["dataset", "target_weight", "target_token"], ascending=[True, False, True], kind="mergesort")
    )
    target_totals["target_rank_dataset"] = (
        target_totals.groupby("dataset", dropna=False).cumcount() + 1
    ).astype(int)
    target_totals["target_display_label"] = target_totals["target_token"].map(target_display_label)
    target_label_totals = (
        target_totals.groupby(["dataset", "target_display_label"], dropna=False, sort=True)
        .agg(target_display_weight=("target_weight", "sum"))
        .reset_index()
        .sort_values(
            ["dataset", "target_display_weight", "target_display_label"],
            ascending=[True, False, True],
            kind="mergesort",
        )
    )
    target_label_totals["target_display_rank_dataset"] = (
        target_label_totals.groupby("dataset", dropna=False).cumcount() + 1
    ).astype(int)
    combined = combined.merge(
        slice_totals,
        on=["dataset", "task", "perturbation_type"],
        how="left",
        validate="many_to_one",
    ).merge(
        dataset_totals,
        on=["dataset"],
        how="left",
        validate="many_to_one",
    ).merge(
        task_totals,
        on=["dataset", "task"],
        how="left",
        validate="many_to_one",
    ).merge(
        perturbation_totals,
        on=["dataset", "task", "perturbation_type"],
        how="left",
        validate="many_to_one",
    ).merge(
        target_totals[["dataset", "target_token", "target_weight", "target_rank_dataset"]],
        on=["dataset", "target_token"],
        how="left",
        validate="many_to_one",
    ).merge(
        target_label_totals,
        on=["dataset", "target_display_label"],
        how="left",
        validate="many_to_one",
    )

    figure1_context = context_totals.merge(
        slice_totals,
        on=["dataset", "task", "perturbation_type"],
        how="left",
        validate="many_to_one",
    ).merge(
        dataset_totals,
        on=["dataset"],
        how="left",
        validate="many_to_one",
    ).merge(
        task_totals,
        on=["dataset", "task"],
        how="left",
        validate="many_to_one",
    ).merge(
        perturbation_totals,
        on=["dataset", "task", "perturbation_type"],
        how="left",
        validate="many_to_one",
    )
    figure1_context["wedge_value"] = figure1_context["context_weight"]
    figure1_context["task_order"] = figure1_context["task"].map({"Task1": 1, "Task2": 2}).fillna(99).astype(int)
    figure1_context["dataset_order"] = figure1_context["dataset"].map({"LINCS": 1, "scPerturb": 2}).fillna(99).astype(int)
    figure1_context["cell_line"] = figure1_context["context_label"]
    figure1_context["target_token"] = np.where(
        figure1_context["context_folded_bool"],
        "Others",
        figure1_context["context_label"].astype(str),
    )
    figure1_context["target_display_label"] = figure1_context["context_note"]
    figure1_context["cell_line_weight"] = figure1_context["context_weight"]
    figure1_context["leaf_weight"] = figure1_context["context_weight"]
    figure1_context["target_weight"] = figure1_context["context_weight"]
    figure1_context["target_display_weight"] = figure1_context["context_weight"]
    figure1_context["target_rank_dataset"] = figure1_context["context_rank_within_parent"]
    figure1_context["target_display_rank_dataset"] = figure1_context["context_rank_within_parent"]
    figure1_context["task_rank_dataset"] = (
        figure1_context["task"].map({"Task1": 1, "Task2": 2}).fillna(99).astype(int)
    )
    figure1_context["hierarchy_level1"] = figure1_context["task"]
    figure1_context["hierarchy_level2"] = figure1_context["perturbation_type"]
    figure1_context["hierarchy_level3"] = figure1_context["context_label"]
    figure1_context["hierarchy_level4"] = pd.NA
    figure1_context["hierarchy_path"] = (
        figure1_context["dataset"].astype(str)
        + "|"
        + figure1_context["hierarchy_level1"].astype(str)
        + "|"
        + figure1_context["hierarchy_level2"].astype(str)
        + "|"
        + figure1_context["hierarchy_level3"].astype(str)
    )
    figure1_context["outer_label"] = figure1_context["context_summary_label"]
    figure1_context["cell_line_label"] = figure1_context["context_summary_label"]
    figure1_context["target_label_keep_bool"] = figure1_context["context_label_keep_bool"]
    figure1_context["cell_line_label_keep_bool"] = figure1_context["context_label_keep_bool"]
    figure1_context["legend_ring_label"] = figure1_context["perturbation_type"].astype(str)
    figure1_context["target_label_rule"] = (
        "Outer ring shows ranked cell-line contexts; low-share contexts are folded into Others within dataset × task × perturbation."
    )
    figure1_context["cell_line_label_rule"] = "Prefer direct on-arc labels; Others remains an explicit folded slice."
    figure1_context = figure1_context[
        [
            "task",
            "task_order",
            "dataset",
            "dataset_order",
            "perturbation_type",
            "task_weight",
            "cell_line",
            "target_token",
            "target_display_label",
            "cell_line_count",
            "target_count",
            "instance_count",
            "eligible_target_count",
            "dataset_cell_line_count",
            "dataset_target_count",
            "dataset_leaf_weight",
            "perturbation_weight",
            "leaf_weight",
            "wedge_value",
            "cell_line_weight",
            "target_weight",
            "target_display_weight",
            "perturbation_rank_dataset",
            "context_rank_within_parent",
            "target_rank_dataset",
            "target_display_rank_dataset",
            "target_label_keep_bool",
            "cell_line_label_keep_bool",
            "n_instances",
            "n_chem_instances",
            "n_gen_instances",
            "context_label",
            "context_note",
            "context_summary_label",
            "context_folded_bool",
            "context_cell_line_count",
            "context_target_count",
            "hierarchy_level1",
            "hierarchy_level2",
            "hierarchy_level3",
            "hierarchy_level4",
            "hierarchy_path",
            "outer_label",
            "cell_line_label",
            "legend_ring_label",
            "target_label_rule",
            "cell_line_label_rule",
            "source_table",
            "scope_note",
        ]
    ].sort_values(
        [
            "dataset_order",
            "task_order",
            "perturbation_rank_dataset",
            "context_rank_within_parent",
            "context_label",
        ],
        kind="mergesort",
    )
    assert_unique_rows(
        figure1_context,
        ["task", "dataset", "perturbation_type", "context_label"],
        "Figure 1 coverage",
    )
    return write_csv(
        figure1_context,
        output_path,
        [
            "task",
            "task_order",
            "dataset",
            "dataset_order",
            "perturbation_type",
            "task_weight",
            "cell_line",
            "target_token",
            "target_display_label",
            "cell_line_count",
            "target_count",
            "instance_count",
            "eligible_target_count",
            "dataset_cell_line_count",
            "dataset_target_count",
            "dataset_leaf_weight",
            "perturbation_weight",
            "leaf_weight",
            "wedge_value",
            "cell_line_weight",
            "target_weight",
            "target_display_weight",
            "perturbation_rank_dataset",
            "context_rank_within_parent",
            "target_rank_dataset",
            "target_display_rank_dataset",
            "target_label_keep_bool",
            "cell_line_label_keep_bool",
            "n_instances",
            "n_chem_instances",
            "n_gen_instances",
            "context_label",
            "context_note",
            "context_summary_label",
            "context_folded_bool",
            "context_cell_line_count",
            "context_target_count",
            "hierarchy_level1",
            "hierarchy_level2",
            "hierarchy_level3",
            "hierarchy_level4",
            "hierarchy_path",
            "outer_label",
            "cell_line_label",
            "legend_ring_label",
            "target_label_rule",
            "cell_line_label_rule",
            "source_table",
            "scope_note",
        ],
    )


def build_figure1_scope_matrix(output_path: Path) -> Path:
    source_ref = (
        "docs/manuscript_master.md|docs/contracts/project-positioning.md|"
        "docs/contracts/task1_spec.md|docs/contracts/task2_spec.md|"
        "docs/contracts/exclusions-and-policies.md"
    )
    rows = [
        ("Task1", "internal", "LINCS", "Chemical", "Gene", "materialized", "lawful_internal_task1", "Task1 internal LINCS base representation."),
        ("Task1", "internal", "LINCS", "Chemical", "Pathway", "materialized", "lawful_internal_task1", "Task1 internal LINCS base representation."),
        ("Task1", "internal", "LINCS", "Chemical", "FM", "not_applicable_scope", "fm_not_supported_for_lincs", "FM excluded because the slice involves LINCS."),
        ("Task1", "internal", "LINCS", "Genetic", "Gene", "materialized", "lawful_internal_task1", "Task1 internal LINCS base representation."),
        ("Task1", "internal", "LINCS", "Genetic", "Pathway", "materialized", "lawful_internal_task1", "Task1 internal LINCS base representation."),
        ("Task1", "internal", "LINCS", "Genetic", "FM", "not_applicable_scope", "fm_not_supported_for_lincs", "FM excluded because the slice involves LINCS."),
        ("Task1", "internal", "scPerturb", "Chemical", "Gene", "materialized", "lawful_internal_task1", "Task1 internal scPerturb base representation."),
        ("Task1", "internal", "scPerturb", "Chemical", "Pathway", "materialized", "lawful_internal_task1", "Task1 internal scPerturb base representation."),
        ("Task1", "internal", "scPerturb", "Chemical", "FM", "materialized", "fm_local_scope", "FM included because the slice is scPerturb-only."),
        ("Task1", "internal", "scPerturb", "Genetic", "Gene", "materialized", "lawful_internal_task1", "Task1 internal scPerturb base representation."),
        ("Task1", "internal", "scPerturb", "Genetic", "Pathway", "materialized", "lawful_internal_task1", "Task1 internal scPerturb base representation."),
        ("Task1", "internal", "scPerturb", "Genetic", "FM", "materialized", "fm_local_scope", "FM included because the slice is scPerturb-only."),
        ("Task1", "cross", "LINCS_to_scPerturb", "Chemical", "Gene", "excluded_by_support_gate", "cross_chemical_policy_exclusion", "Cross chemical is policy-excluded and must not be framed as attrition."),
        ("Task1", "cross", "LINCS_to_scPerturb", "Chemical", "Pathway", "excluded_by_support_gate", "cross_chemical_policy_exclusion", "Cross chemical is policy-excluded and must not be framed as attrition."),
        ("Task1", "cross", "LINCS_to_scPerturb", "Chemical", "FM", "not_applicable_scope", "task1_cross_gene_pathway_only", "FM excluded because Task1 cross is frozen to Gene/Pathway only."),
        ("Task1", "cross", "LINCS_to_scPerturb", "Genetic", "Gene", "materialized", "cross_genetic_common_scope", "Task1 cross genetic common-scope row."),
        ("Task1", "cross", "LINCS_to_scPerturb", "Genetic", "Pathway", "materialized", "cross_genetic_common_scope", "Task1 cross genetic common-scope row."),
        ("Task1", "cross", "LINCS_to_scPerturb", "Genetic", "FM", "not_applicable_scope", "task1_cross_gene_pathway_only", "FM excluded because Task1 cross is frozen to Gene/Pathway only."),
        ("Task1", "cross", "scPerturb_to_LINCS", "Genetic", "Gene", "materialized", "cross_genetic_common_scope", "Task1 cross genetic common-scope row."),
        ("Task1", "cross", "scPerturb_to_LINCS", "Genetic", "Pathway", "materialized", "cross_genetic_common_scope", "Task1 cross genetic common-scope row."),
        ("Task1", "cross", "scPerturb_to_LINCS", "Genetic", "FM", "not_applicable_scope", "task1_cross_gene_pathway_only", "FM excluded because Task1 cross is frozen to Gene/Pathway only."),
        ("Task2", "dataset_cell_line", "LINCS", "Paired", "Gene", "materialized", "corrected_task2_scope", "Corrected Task2 LINCS base representation."),
        ("Task2", "dataset_cell_line", "LINCS", "Paired", "Pathway", "materialized", "corrected_task2_scope", "Corrected Task2 LINCS base representation."),
        ("Task2", "dataset_cell_line", "LINCS", "Paired", "FM", "not_applicable_scope", "fm_not_supported_for_lincs", "FM excluded because the slice involves LINCS."),
        ("Task2", "dataset_cell_line", "scPerturb_K562", "Paired", "Gene", "materialized", "corrected_task2_scope", "Corrected Task2 scPerturb K562 base representation."),
        ("Task2", "dataset_cell_line", "scPerturb_K562", "Paired", "Pathway", "materialized", "corrected_task2_scope", "Corrected Task2 scPerturb K562 base representation."),
        ("Task2", "dataset_cell_line", "scPerturb_K562", "Paired", "FM", "materialized", "fm_local_scope", "FM included because the slice is scPerturb K562 only."),
    ]
    frame = pd.DataFrame(
        rows,
        columns=[
            "task",
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation_class",
            "status",
            "status_reason",
            "scope_note",
        ],
    )
    frame["source_ref"] = source_ref
    assert_unique_rows(
        frame,
        ["task", "scope", "dataset_or_direction", "perturbation_type", "representation_class"],
        "Figure 1 lawful scope matrix",
    )
    return write_csv(
        frame.sort_values(
            ["task", "scope", "dataset_or_direction", "perturbation_type", "representation_class"],
            kind="mergesort",
        ),
        output_path,
        [
            "task",
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation_class",
            "status",
            "status_reason",
            "source_ref",
            "scope_note",
        ],
    )


def build_figure1_rep_modifier_matrix(
    *,
    representation_registry_path: Path,
    figure3_scope_path: Path,
    output_path: Path,
) -> Path:
    registry = read_csv_required(
        representation_registry_path,
        {"dataset", "cell_line", "representation", "availability_status", "availability_reason"},
        "Representation availability registry",
    )
    scope = read_csv_required(
        figure3_scope_path,
        {
            "dataset",
            "cell_line",
            "representation",
            "availability_status",
            "availability_reason",
            "group_slice_materialized_bool",
            "retrieval_c2g_materialized_bool",
            "retrieval_g2c_materialized_bool",
        },
        "Figure 3 scope summary",
    )
    registry = registry.copy()
    registry["representation"] = registry["representation"].map(representation_class)
    registry["cell_line_scope"] = np.where(
        registry["dataset"].eq("scPerturb"),
        "K562_only",
        "multi_cell_line",
    )
    task2_rep = (
        registry.groupby(["dataset", "cell_line_scope", "representation"], dropna=False, sort=True)
        .agg(
            representation_status=("availability_status", lambda values: "available" if (pd.Series(values) == "available").any() else "not_applicable_scope"),
            availability_reason=("availability_reason", stringify_series),
        )
        .reset_index()
    )
    task2_rep["task_surface"] = "Task2_common"
    task2_rep["dataset_scope"] = task2_rep["dataset"]

    scope = scope.copy()
    scope["representation"] = scope["representation"].map(representation_class)
    scope["cell_line_scope"] = np.where(scope["dataset"].eq("scPerturb"), "K562_only", "multi_cell_line")
    modifier = (
        scope.groupby(["dataset", "cell_line_scope", "representation"], dropna=False, sort=True)
        .agg(
            any_c2g=("retrieval_c2g_materialized_bool", "any"),
        )
        .reset_index()
    )
    modifier["task_surface"] = "Task2_c2g_modifiers"
    modifier["dataset_scope"] = modifier["dataset"]
    modifier["representation_status"] = np.where(modifier["any_c2g"], "materialized", "not_applicable_scope")
    modifier["availability_reason"] = np.where(
        modifier["any_c2g"],
        "supplementary_c2g_modifier_surface_materialized",
        "modifier_not_applicable_for_scope",
    )

    task1_rows = pd.DataFrame(
        [
            ("Task1_internal", "LINCS", "multi_cell_line", "Gene", "available", "task1_internal_base_representation"),
            ("Task1_internal", "LINCS", "multi_cell_line", "Pathway", "available", "task1_internal_base_representation"),
            ("Task1_internal", "LINCS", "multi_cell_line", "FM", "not_applicable_scope", "fm_not_supported_for_lincs"),
            ("Task1_internal", "scPerturb", "multi_cell_line", "Gene", "available", "task1_internal_base_representation"),
            ("Task1_internal", "scPerturb", "multi_cell_line", "Pathway", "available", "task1_internal_base_representation"),
            ("Task1_internal", "scPerturb", "multi_cell_line", "FM", "available", "fm_local_scope"),
            ("Task1_cross", "cross_genetic", "paired_direction", "Gene", "materialized", "task1_cross_gene_pathway_only"),
            ("Task1_cross", "cross_genetic", "paired_direction", "Pathway", "materialized", "task1_cross_gene_pathway_only"),
            ("Task1_cross", "cross_genetic", "paired_direction", "FM", "not_applicable_scope", "task1_cross_gene_pathway_only"),
        ],
        columns=[
            "task_surface",
            "dataset_scope",
            "cell_line_scope",
            "representation",
            "representation_status",
            "availability_reason",
        ],
    )
    task1_rows["time_modifier_status"] = "not_applicable_scope"
    task1_rows["dose_modifier_status"] = "not_applicable_scope"
    task1_rows["target_multiplicity_status"] = "not_applicable_scope"
    task1_rows["scope_note"] = "Task1 rows carry representation availability only; modifier families are Task2 C2G supplementary only."

    task2_rep["time_modifier_status"] = "not_applicable_scope"
    task2_rep["dose_modifier_status"] = "not_applicable_scope"
    task2_rep["target_multiplicity_status"] = "not_applicable_scope"
    task2_rep["scope_note"] = "Task2 common representation availability under the corrected dataset/cell_line scope."

    modifier["time_modifier_status"] = modifier["representation_status"]
    modifier["dose_modifier_status"] = modifier["representation_status"]
    modifier["target_multiplicity_status"] = modifier["representation_status"]
    modifier["scope_note"] = "Modifier families are supplementary, C2G-only, and materialized only where the local retrieval surface exists."

    combined = pd.concat(
        [
            task1_rows,
            task2_rep[
                [
                    "task_surface",
                    "dataset_scope",
                    "cell_line_scope",
                    "representation",
                    "representation_status",
                    "availability_reason",
                    "time_modifier_status",
                    "dose_modifier_status",
                    "target_multiplicity_status",
                    "scope_note",
                ]
            ],
            modifier[
                [
                    "task_surface",
                    "dataset_scope",
                    "cell_line_scope",
                    "representation",
                    "representation_status",
                    "availability_reason",
                    "time_modifier_status",
                    "dose_modifier_status",
                    "target_multiplicity_status",
                    "scope_note",
                ]
            ],
        ],
        ignore_index=True,
        sort=False,
    ).sort_values(
        ["task_surface", "dataset_scope", "cell_line_scope", "representation"],
        kind="mergesort",
    )
    assert_unique_rows(
        combined,
        ["task_surface", "dataset_scope", "cell_line_scope", "representation"],
        "Figure 1 representation/modifier availability",
    )
    return write_csv(
        combined,
        output_path,
        [
            "task_surface",
            "dataset_scope",
            "cell_line_scope",
            "representation",
            "representation_status",
            "time_modifier_status",
            "dose_modifier_status",
            "target_multiplicity_status",
            "availability_reason",
            "scope_note",
        ],
    )


def build_figure1_object_map(*, manifest_path: Path, output_path: Path) -> Path:
    require_path(manifest_path, "Framework manifest")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows: list[dict[str, object]] = []

    canonical_to_panel = {
        "figure2_task1_scope_summary.csv": ("F2", "2A", "main"),
        "figure2_task1_performance_structure.csv": ("F2", "2B|2C", "main"),
        "figure2_task1_internal_to_cross_degradation_summary.csv": ("F2", "2D", "main"),
        "figure3_task2_scope_summary.csv": ("F3", "1D", "supporting"),
        "figure3_task2_direction_support_summary.csv": ("F3", "3B|EF5B|EF5C", "main"),
        "figure3_task2_performance_structure.csv": ("F3", "3B|EF5A", "main"),
        "figure3_task2_cell_line_pattern_summary.csv": ("F3", "3C|EF6A", "main"),
        "figure3_task2_target_pattern_summary.csv": ("F3", "3D|3F|EF6B|EF8B", "main"),
        "figure3_task1_internal_contextual_support_summary.csv": ("F3", "3J|EF8C", "supporting"),
        "manuscript_comparison_statistics.csv": ("F2|F3|EF8", "2C|3E|3F|EF8A|EF8B", "supporting"),
    }

    for entry in manifest.get("canonical_objects", []):
        filename = str(entry.get("canonical_filename", ""))
        figure_target, panel_target, disposition = canonical_to_panel.get(
            filename, ("unmapped", "unmapped", "supporting")
        )
        rows.append(
            {
                "object_name": entry.get("object_name", filename),
                "object_status": "canonical",
                "figure_target": figure_target,
                "panel_target": panel_target,
                "canonical_filename": filename,
                "source_class": "canonical_analysis",
                "disposition": disposition,
                "notes": entry.get("justification", ""),
            }
        )

    for support_entry in FIGURE1_OBJECT_MAP_SUPPORT_OBJECTS:
        rows.append(
            {
                "object_name": support_entry["object_name"],
                "object_status": support_entry["object_status"],
                "figure_target": support_entry["figure_target"],
                "panel_target": support_entry["panel_target"],
                "canonical_filename": support_entry["canonical_filename"],
                "source_class": support_entry["source_class"],
                "disposition": support_entry["disposition"],
                "notes": support_entry["notes"],
            }
        )

    for entry in manifest.get("historical_noncanonical_objects", []):
        rows.append(
            {
                "object_name": entry.get("object_name", Path(str(entry.get("path", ""))).name),
                "object_status": "historical_only",
                "figure_target": "historical",
                "panel_target": "historical",
                "canonical_filename": Path(str(entry.get("path", ""))).name,
                "source_class": "historical_noncanonical",
                "disposition": entry.get("disposition", "historical_only"),
                "notes": entry.get("reason_no_longer_canonical", ""),
            }
        )

    frame = pd.DataFrame(rows).sort_values(["object_status", "figure_target", "panel_target", "canonical_filename"], kind="mergesort")
    assert_unique_rows(frame, ["object_status", "canonical_filename", "panel_target"], "Figure 1 object map")
    return write_csv(
        frame,
        output_path,
        [
            "object_name",
            "object_status",
            "figure_target",
            "panel_target",
            "canonical_filename",
            "source_class",
            "disposition",
            "notes",
        ],
    )


def build_extended_figure1_registry(
    *,
    coverage_path: Path,
    scope_path: Path,
    rep_modifier_path: Path,
    object_map_path: Path,
    output_path: Path,
) -> Path:
    coverage = pd.read_csv(coverage_path)
    scope = pd.read_csv(scope_path)
    rep_modifier = pd.read_csv(rep_modifier_path)
    object_map = pd.read_csv(object_map_path)

    rows: list[dict[str, object]] = []
    coverage_registry = coverage[
        ["task", "dataset", "perturbation_type", "source_table", "scope_note"]
    ].drop_duplicates()
    for row in coverage_registry.itertuples(index=False):
        item_id = f"{row.task}|{row.dataset}|{row.perturbation_type}"
        rows.append(
            {
                "figure_id": "EF1",
                "panel_id": "EF1B",
                "registry_type": "coverage",
                "item_id": item_id,
                "item_label": item_id,
                "source_file": row.source_table,
                "status": "materialized",
                "scope_note": row.scope_note,
                "governing_doc": "docs/manuscript_master.md",
            }
        )
    for row in scope.itertuples(index=False):
        item_id = "|".join(
            [row.task, row.scope, row.dataset_or_direction, row.perturbation_type, row.representation_class]
        )
        rows.append(
            {
                "figure_id": "EF1",
                "panel_id": "EF1C",
                "registry_type": "lawful_scope",
                "item_id": item_id,
                "item_label": item_id,
                "source_file": row.source_ref,
                "status": row.status,
                "scope_note": row.scope_note,
                "governing_doc": "docs/contracts/exclusions-and-policies.md",
            }
        )
    for row in rep_modifier.itertuples(index=False):
        item_id = "|".join([row.task_surface, row.dataset_scope, row.cell_line_scope, row.representation])
        rows.append(
            {
                "figure_id": "EF1",
                "panel_id": "EF1D",
                "registry_type": "representation_modifier",
                "item_id": item_id,
                "item_label": item_id,
                "source_file": str(DEFAULT_REPRESENTATION_REGISTRY_PATH),
                "status": row.representation_status,
                "scope_note": row.scope_note,
                "governing_doc": "docs/manuscript_master.md",
            }
        )
    for row in object_map.itertuples(index=False):
        rows.append(
            {
                "figure_id": "EF1",
                "panel_id": "EF1E",
                "registry_type": "result_object",
                "item_id": row.canonical_filename,
                "item_label": row.object_name,
                "source_file": str(DEFAULT_FRAMEWORK_MANIFEST_PATH),
                "status": row.object_status,
                "scope_note": row.notes,
                "governing_doc": "scripts/ARCHITECTURE.md",
            }
        )
    frame = pd.DataFrame(rows).sort_values(["panel_id", "registry_type", "item_id"], kind="mergesort")
    assert_unique_rows(frame, ["panel_id", "registry_type", "item_id"], "Extended Figure 1 support registry")
    return write_csv(
        frame,
        output_path,
        [
            "figure_id",
            "panel_id",
            "registry_type",
            "item_id",
            "item_label",
            "source_file",
            "status",
            "scope_note",
            "governing_doc",
        ],
    )


def load_filtered_high_concordance(path: Path, entity_column: str) -> pd.DataFrame:
    required = {
        *FIGURE2_FAMILY_COLUMNS,
        entity_column,
        "n_queries_total",
        "n_positive",
        "n_negative",
        "n_background_positive",
        "n_background_negative",
        "success_rate",
        "background_success_rate",
        "odds_ratio",
        "raw_p",
        "bh_q",
        "significant_bool",
        "underpowered_strata_bool",
    }
    frame = read_csv_required(path, required, f"Task1 high-concordance {entity_column}")
    frame = frame.copy()
    if entity_column == "target_token":
        invalid = frame["target_token"].fillna("").astype(str).str.contains(r"[|_]", regex=True)
        if invalid.any():
            raise ValueError(
                "Task1 target high-concordance summary must contain atomic target_token rows before plot-ready staging."
            )
    frame["underpowered_strata_bool"] = frame["underpowered_strata_bool"].astype(bool)
    filtered = frame.loc[
        (~frame["underpowered_strata_bool"])
        & frame["representation"].isin(HIGH_CONCORDANCE_COMMON_SCOPE_REPRESENTATIONS)
        & (
            frame["scope"].eq("internal")
            | (frame["scope"].eq("cross") & frame["perturbation_type"].eq("Genetic"))
        )
    ].copy()
    return filtered


def build_repaired_target_high_concordance_frame(source_path: Path) -> pd.DataFrame:
    required = {
        *FIGURE2_FAMILY_COLUMNS,
        "target_token",
        "n_queries_total",
        "n_positive",
        "n_negative",
        "n_background_positive",
        "n_background_negative",
        "success_rate",
        "background_success_rate",
        "odds_ratio",
        "raw_p",
        "bh_q",
        "significant_bool",
        "success_definition",
        "source_path",
        "underpowered_strata_bool",
    }
    out = read_csv_required(source_path, required, "Task1 target high-concordance summary").copy()
    invalid = out["target_token"].fillna("").astype(str).str.contains(r"[|_]", regex=True)
    if invalid.any():
        raise ValueError(
            "Task1 target high-concordance summary still contains multi-target target_token rows; repair the analysis-layer builder first."
        )
    return out.sort_values([*FIGURE2_FAMILY_COLUMNS, "target_token"], kind="mergesort")


def repair_target_high_concordance_summary(*, source_path: Path, output_path: Path) -> Path:
    out = build_repaired_target_high_concordance_frame(source_path)
    return write_csv(
        out,
        output_path,
        [
            *FIGURE2_FAMILY_COLUMNS,
            "target_token",
            "n_queries_total",
            "n_positive",
            "n_negative",
            "n_background_positive",
            "n_background_negative",
            "success_rate",
            "background_success_rate",
            "odds_ratio",
            "raw_p",
            "bh_q",
            "significant_bool",
            "success_definition",
            "source_path",
            "underpowered_strata_bool",
        ],
    )


def pooled_high_concordance_summary(frame: pd.DataFrame, entity_column: str) -> pd.DataFrame:
    grouped = (
        frame.groupby([entity_column, "representation"], dropna=False, sort=True)
        .agg(
            n_queries_total=("n_queries_total", "sum"),
            n_positive=("n_positive", "sum"),
            n_negative=("n_negative", "sum"),
            n_background_positive=("n_background_positive", "sum"),
            n_background_negative=("n_background_negative", "sum"),
            min_bh_q=("bh_q", "min"),
            any_significant=("significant_bool", lambda values: bool(pd.Series(values).fillna(False).astype(bool).any())),
            n_strata_included=(entity_column, "size"),
            scope_set=("scope", stringify_series),
            dataset_set=("dataset_or_direction", stringify_series),
        )
        .reset_index()
    )
    grouped["success_rate_pooled"] = np.where(
        grouped["n_queries_total"].gt(0),
        grouped["n_positive"] / grouped["n_queries_total"],
        np.nan,
    )
    background_totals = grouped["n_background_positive"] + grouped["n_background_negative"]
    grouped["background_success_rate_pooled"] = np.where(
        background_totals.gt(0),
        grouped["n_background_positive"] / background_totals,
        np.nan,
    )
    grouped["odds_ratio_pooled"] = [
        pooled_odds_ratio(a, b, c, d)
        for a, b, c, d in zip(
            grouped["n_positive"],
            grouped["n_negative"],
            grouped["n_background_positive"],
            grouped["n_background_negative"],
        )
    ]
    grouped = grouped.sort_values(
        ["representation", "success_rate_pooled", "odds_ratio_pooled", entity_column],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    grouped["rank_within_representation"] = (
        grouped.groupby("representation", dropna=False).cumcount() + 1
    ).astype("Int64")
    return grouped


def add_label_flag(frame: pd.DataFrame, group_column: str, top_n: int = 12) -> pd.DataFrame:
    out = frame.copy()
    out["label_bool"] = out.groupby(group_column, dropna=False)["rank_within_representation"].transform(
        lambda values: values.le(top_n)
    ).astype(bool)
    return out


def load_task1_scope_summary(path: Path) -> pd.DataFrame:
    required = {
        "scope",
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "analysis_family",
        "scope_status",
        "metric_names",
        "metric_count",
        "n_total",
        "n_valid",
        "n_excluded",
        "N_gallery_max",
        "cross_alignment_contract",
        "n_matched_keys",
        "eligible_bool",
        "exclusion_reason",
        "scope_note",
        "fm_scope_note",
    }
    return read_csv_required(path, required, "Task1 scope summary")


def task1_scope_panel_axis(frame: pd.DataFrame) -> pd.Series:
    representation = frame["representation"].astype(str)
    out = np.where(
        frame["analysis_family"].eq("cross_eligibility") | representation.eq("ALL"),
        "ALL",
        np.where(representation.isin(["Gene", "Pathway"]), representation, "FM"),
    )
    return pd.Series(out, index=frame.index, dtype="object")


def task1_scope_bar_value(frame: pd.DataFrame) -> pd.Series:
    return pd.Series(
        np.where(
            frame["analysis_family"].eq("cross_eligibility"),
            pd.to_numeric(frame["n_matched_keys"], errors="coerce"),
            pd.to_numeric(frame["n_valid"], errors="coerce"),
        ),
        index=frame.index,
        dtype="float64",
    )


def task1_cross_direction_label(dataset_value: object) -> str:
    text = str(dataset_value)
    if text in {"LINCS", "LINCS_to_scPerturb"}:
        return "LINCS -> scPerturb"
    if text in {"scPerturb", "scPerturb_to_LINCS"}:
        return "scPerturb -> LINCS"
    return text


def simple_ci_bounds(series: pd.Series) -> tuple[float | pd._libs.missing.NAType, float | pd._libs.missing.NAType]:
    numeric = pd.to_numeric(series, errors="coerce").dropna()
    if numeric.empty:
        return (pd.NA, pd.NA)
    if len(numeric) == 1:
        value = float(numeric.iloc[0])
        return (value, value)
    mean_value = float(numeric.mean())
    se_value = float(numeric.std(ddof=1) / math.sqrt(len(numeric)))
    delta = 1.96 * se_value
    return (mean_value - delta, mean_value + delta)


def summary_quantiles(series: pd.Series) -> dict[str, float | pd._libs.missing.NAType]:
    numeric = pd.to_numeric(series, errors="coerce").dropna()
    if numeric.empty:
        return {
            "summary_mean": pd.NA,
            "summary_median": pd.NA,
            "q10_value": pd.NA,
            "q25_value": pd.NA,
            "q50_value": pd.NA,
            "q75_value": pd.NA,
            "q90_value": pd.NA,
        }
    quantiles = numeric.quantile([0.10, 0.25, 0.50, 0.75, 0.90])
    return {
        "summary_mean": float(numeric.mean()),
        "summary_median": float(numeric.median()),
        "q10_value": float(quantiles.loc[0.10]),
        "q25_value": float(quantiles.loc[0.25]),
        "q50_value": float(quantiles.loc[0.50]),
        "q75_value": float(quantiles.loc[0.75]),
        "q90_value": float(quantiles.loc[0.90]),
    }


def summarize_quantiles_frame(
    frame: pd.DataFrame,
    *,
    group_columns: Sequence[str],
    value_column: str,
) -> pd.DataFrame:
    grouped = frame.groupby(list(group_columns), dropna=False, sort=True)[value_column]
    return grouped.agg(
        summary_mean=lambda values: summary_quantiles(values)["summary_mean"],
        summary_median=lambda values: summary_quantiles(values)["summary_median"],
        q10_value=lambda values: summary_quantiles(values)["q10_value"],
        q25_value=lambda values: summary_quantiles(values)["q25_value"],
        q50_value=lambda values: summary_quantiles(values)["q50_value"],
        q75_value=lambda values: summary_quantiles(values)["q75_value"],
        q90_value=lambda values: summary_quantiles(values)["q90_value"],
    ).reset_index()


def task1_family_order(scope: object, dataset_or_direction: object, perturbation_type: object) -> int:
    key = (str(scope), str(dataset_or_direction), str(perturbation_type))
    order = {
        ("internal", "LINCS", "Chemical"): 1,
        ("internal", "LINCS", "Genetic"): 2,
        ("internal", "scPerturb", "Chemical"): 3,
        ("internal", "scPerturb", "Genetic"): 4,
        ("cross", "LINCS_to_scPerturb", "Genetic"): 5,
        ("cross", "scPerturb_to_LINCS", "Genetic"): 6,
    }
    return order.get(key, 999)


def task1_family_label(scope: object, dataset_or_direction: object, perturbation_type: object) -> str:
    scope_text = str(scope)
    dataset_text = str(dataset_or_direction)
    perturbation_text = str(perturbation_type)
    if scope_text == "internal":
        return f"Internal | {dataset_text} | {perturbation_text}"
    return f"Cross | {task1_cross_direction_label(dataset_text)} | {perturbation_text}"


def rank_bin_series(values: pd.Series, *, n_bins: int = 5) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series(pd.array([pd.NA] * len(values), dtype="Int64"), index=values.index)
    valid = numeric.notna()
    if valid.sum() == 0:
        return out
    ranks = numeric.loc[valid].rank(method="first", ascending=True)
    bins = np.ceil(ranks / len(ranks) * n_bins).astype(int)
    bins = bins.clip(1, n_bins)
    out.loc[valid] = bins.astype("Int64")
    return out


def log2_enrichment_score(frame: pd.DataFrame) -> pd.Series:
    odds = (
        (pd.to_numeric(frame["n_positive"], errors="coerce") + 0.5)
        / (pd.to_numeric(frame["n_negative"], errors="coerce") + 0.5)
    ) / (
        (pd.to_numeric(frame["n_background_positive"], errors="coerce") + 0.5)
        / (pd.to_numeric(frame["n_background_negative"], errors="coerce") + 0.5)
    )
    return np.log2(odds.astype("float64"))


def fisher_exact_raw_p(
    n_positive: object,
    n_negative: object,
    n_background_positive: object,
    n_background_negative: object,
) -> float | pd._libs.missing.NAType:
    counts = [
        pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        for value in (n_positive, n_negative, n_background_positive, n_background_negative)
    ]
    if any(pd.isna(value) for value in counts):
        return pd.NA
    n_pos, n_neg, bg_pos, bg_neg = [int(value) for value in counts]
    if min(n_pos, n_neg, bg_pos, bg_neg) < 0:
        return pd.NA
    if (n_pos + n_neg) == 0 or (bg_pos + bg_neg) == 0:
        return pd.NA
    _, raw_p = fisher_exact([[n_pos, n_neg], [bg_pos, bg_neg]])
    return float(raw_p)


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


def add_enrichment_statistics(
    frame: pd.DataFrame,
    *,
    adjust_group_columns: Sequence[str],
) -> pd.DataFrame:
    out = frame.copy()
    out["success_rate"] = np.where(
        (pd.to_numeric(out["n_positive"], errors="coerce") + pd.to_numeric(out["n_negative"], errors="coerce")).gt(0),
        pd.to_numeric(out["n_positive"], errors="coerce")
        / (pd.to_numeric(out["n_positive"], errors="coerce") + pd.to_numeric(out["n_negative"], errors="coerce")),
        np.nan,
    )
    background_total = pd.to_numeric(out["n_background_positive"], errors="coerce") + pd.to_numeric(
        out["n_background_negative"], errors="coerce"
    )
    out["background_success_rate"] = np.where(
        background_total.gt(0),
        pd.to_numeric(out["n_background_positive"], errors="coerce") / background_total,
        np.nan,
    )
    out["odds_ratio"] = [
        pooled_odds_ratio(a, b, c, d)
        for a, b, c, d in zip(
            out["n_positive"],
            out["n_negative"],
            out["n_background_positive"],
            out["n_background_negative"],
        )
    ]
    out["raw_p"] = [
        fisher_exact_raw_p(a, b, c, d)
        for a, b, c, d in zip(
            out["n_positive"],
            out["n_negative"],
            out["n_background_positive"],
            out["n_background_negative"],
        )
    ]
    out["bh_q"] = pd.Series(pd.array([pd.NA] * len(out), dtype="Float64"), index=out.index)
    for _, index in out.groupby(list(adjust_group_columns), dropna=False, sort=False).groups.items():
        raw = pd.to_numeric(out.loc[index, "raw_p"], errors="coerce")
        valid = raw.notna()
        if valid.any():
            q_values = bh_adjust(raw.loc[valid].astype("float64"))
            out.loc[raw.loc[valid].index, "bh_q"] = q_values.astype("Float64")
    out["significant_bool"] = pd.to_numeric(out["bh_q"], errors="coerce").le(0.05).fillna(False)
    out["significance_label"] = np.where(
        pd.to_numeric(out["bh_q"], errors="coerce").le(0.05),
        "*",
        "",
    )
    out["enrichment_score"] = log2_enrichment_score(out)
    return out


def aggregation_support_label(value: object) -> str:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return ""
    return f"n={format(int(round(float(numeric))), ',')}"


def ordered_exemplar_candidates(
    frame: pd.DataFrame,
    *,
    score_column: str,
    direction: str,
) -> pd.DataFrame:
    ascending = direction == "low"
    return frame.sort_values(
        [score_column, "max_abs_enrichment", "support_value", "entity_display_label"],
        ascending=[ascending, False, False, True],
        kind="mergesort",
    )


def pick_direction_indexes(
    frame: pd.DataFrame,
    *,
    score_column: str,
    direction: str,
    n: int,
    excluded_indexes: set[int] | None = None,
) -> list[int]:
    excluded_indexes = excluded_indexes or set()
    ordered = ordered_exemplar_candidates(frame, score_column=score_column, direction=direction)
    picked: list[int] = []
    for row in ordered.itertuples():
        if int(row.Index) in excluded_indexes:
            continue
        picked.append(int(row.Index))
        if len(picked) >= n:
            break
    return picked


def select_dual_tail_exemplars(
    frame: pd.DataFrame,
    *,
    entity_column: str,
    dataset_column: str,
    score_column: str,
    high_n: int,
    low_n: int,
    shared_column: str,
    min_shared: int,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    selection_direction_order = {"high": 1, "low": 2}
    for _, subset in frame.groupby(dataset_column, dropna=False, sort=True):
        subset = subset.copy()
        if subset.empty:
            continue
        high_indexes = pick_direction_indexes(
            subset,
            score_column=score_column,
            direction="high",
            n=high_n,
        )
        low_indexes = pick_direction_indexes(
            subset,
            score_column=score_column,
            direction="low",
            n=low_n,
            excluded_indexes=set(high_indexes),
        )
        selected_map = {
            "high": high_indexes,
            "low": low_indexes,
        }
        selection_notes: dict[int, str] = {}
        available_shared = int(subset[shared_column].fillna(False).astype(bool).sum())
        target_shared = min(min_shared, available_shared)
        selected_indexes = set(high_indexes + low_indexes)
        current_shared = int(subset.loc[list(selected_indexes), shared_column].fillna(False).astype(bool).sum()) if selected_indexes else 0
        if target_shared > current_shared:
            for direction in ("high", "low", "high", "low"):
                if current_shared >= target_shared:
                    break
                current_direction = selected_map[direction]
                if not current_direction:
                    continue
                selected_frame = subset.loc[current_direction].copy()
                replaceable = selected_frame.loc[~selected_frame[shared_column].fillna(False).astype(bool)].copy()
                if replaceable.empty:
                    continue
                candidate_indexes = pick_direction_indexes(
                    subset.loc[subset[shared_column].fillna(False).astype(bool)].copy(),
                    score_column=score_column,
                    direction=direction,
                    n=len(subset),
                    excluded_indexes=selected_indexes,
                )
                if not candidate_indexes:
                    continue
                candidate_index = candidate_indexes[0]
                replace_order = ordered_exemplar_candidates(replaceable, score_column=score_column, direction=direction)
                remove_index = int(replace_order.tail(1).index[0])
                selected_map[direction] = [candidate_index if idx == remove_index else idx for idx in current_direction]
                selected_indexes.discard(remove_index)
                selected_indexes.add(candidate_index)
                selection_notes[candidate_index] = "shared_quota_backfill"
                current_shared = int(subset.loc[list(selected_indexes), shared_column].fillna(False).astype(bool).sum())

        selected_frames: list[pd.DataFrame] = []
        for direction in ("high", "low"):
            direction_indexes = selected_map[direction]
            if not direction_indexes:
                continue
            direction_frame = subset.loc[direction_indexes].copy()
            direction_frame = ordered_exemplar_candidates(
                direction_frame,
                score_column=score_column,
                direction=direction,
            ).reset_index(drop=False)
            direction_frame["selection_direction"] = direction
            direction_frame["selection_rank_within_direction"] = np.arange(1, len(direction_frame) + 1, dtype=int)
            direction_frame["selection_note"] = direction_frame["index"].map(selection_notes).fillna("")
            direction_frame = direction_frame.drop(columns=["index"])
            selected_frames.append(direction_frame)
        selected = pd.concat(selected_frames, ignore_index=True, sort=False) if selected_frames else subset.head(0).copy()
        selected["selection_direction_order"] = selected["selection_direction"].map(selection_direction_order).astype(int)
        selected = selected.sort_values(
            ["selection_direction_order", "selection_rank_within_direction", "entity_display_label"],
            kind="mergesort",
        ).reset_index(drop=True)
        selected["row_order"] = np.arange(1, len(selected) + 1, dtype=int)
        rows.append(selected.drop(columns=["selection_direction_order"]))
    out = pd.concat(rows, ignore_index=True, sort=False) if rows else frame.head(0).copy()
    if not out.empty:
        out["row_id"] = (
            out[dataset_column].astype(str)
            + "::"
            + out["selection_direction"].astype(str)
            + "::"
            + out[entity_column].astype(str)
        )
    else:
        out["row_id"] = pd.Series(dtype="object")
    return out


def combine_pair_annotations(values: pd.Series) -> str:
    cleaned = [value for value in values.astype(str).tolist() if value and value != "nan"]
    return "|".join(sorted(dict.fromkeys(cleaned)))


def ranked_entity_order(frame: pd.DataFrame, *, score_column: str) -> pd.DataFrame:
    return frame.sort_values(
        [score_column, "max_abs_enrichment", "support_value", "entity_display_label"],
        ascending=[False, False, False, True],
        kind="mergesort",
    )


def select_ranked_entities(
    frame: pd.DataFrame,
    *,
    dataset_column: str,
    score_column: str,
    top_n: int | None,
    shared_column: str | None = None,
    min_shared: int = 0,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for _, subset in frame.groupby(dataset_column, dropna=False, sort=True):
        subset = ranked_entity_order(subset.copy(), score_column=score_column)
        if subset.empty:
            continue
        selected = subset.copy() if top_n is None or top_n <= 0 else subset.head(top_n).copy()
        if shared_column is not None and min_shared > 0 and shared_column in subset.columns and not subset.empty:
            available_shared = subset.loc[subset[shared_column].fillna(False).astype(bool)].copy()
            target_shared = min(int(len(available_shared)), int(min_shared))
            current_shared = int(selected[shared_column].fillna(False).astype(bool).sum())
            if target_shared > current_shared and not available_shared.empty:
                selected_indices = list(selected.index)
                for candidate in available_shared.itertuples():
                    if candidate.Index in selected_indices:
                        continue
                    replaceable = selected.loc[~selected[shared_column].fillna(False).astype(bool)].copy()
                    if replaceable.empty:
                        break
                    remove_index = int(replaceable.tail(1).index[0])
                    selected_indices = [idx for idx in selected_indices if idx != remove_index] + [int(candidate.Index)]
                    selected = subset.loc[selected_indices].copy()
                    selected = ranked_entity_order(selected, score_column=score_column)
                    current_shared = int(selected[shared_column].fillna(False).astype(bool).sum())
                    if current_shared >= target_shared:
                        break
        selected = ranked_entity_order(selected, score_column=score_column).reset_index(drop=True)
        selected["row_rank"] = np.arange(1, len(selected) + 1, dtype=int)
        rows.append(selected)
    return pd.concat(rows, ignore_index=True, sort=False) if rows else frame.head(0).copy()


def build_paired_task1_exemplar_panel(
    *,
    source_path: Path,
    entity_column: str,
    output_path: Path,
) -> Path:
    filtered = load_filtered_high_concordance(source_path, entity_column)
    filtered = filtered.loc[filtered["scope"].eq("internal")].copy()
    filtered["dataset"] = filtered["dataset_or_direction"].astype(str)
    filtered["entity_display_label"] = filtered[entity_column].map(lambda value: entity_display_label(entity_column, value))
    filtered["display_excluded_bool"] = filtered[entity_column].map(lambda value: entity_display_exclusion(entity_column, value)).astype(bool)
    grouped = (
        filtered.groupby(
            ["dataset", entity_column, "representation", "entity_display_label", "display_excluded_bool"],
            dropna=False,
            sort=True,
        )
        .agg(
            n_queries_total=("n_queries_total", "sum"),
            n_positive=("n_positive", "sum"),
            n_negative=("n_negative", "sum"),
            n_background_positive=("n_background_positive", "sum"),
            n_background_negative=("n_background_negative", "sum"),
            perturbation_type_label=("perturbation_type", stringify_series),
            n_source_strata=("representation", "size"),
        )
        .reset_index()
    )
    grouped = add_enrichment_statistics(grouped, adjust_group_columns=["dataset", "representation"])
    wide_index = ["dataset", entity_column, "entity_display_label", "display_excluded_bool"]
    value_columns = [
        "n_queries_total",
        "n_positive",
        "n_negative",
        "n_background_positive",
        "n_background_negative",
        "success_rate",
        "background_success_rate",
        "odds_ratio",
        "raw_p",
        "bh_q",
        "significant_bool",
        "significance_label",
        "enrichment_score",
        "perturbation_type_label",
        "n_source_strata",
    ]
    gene = grouped.loc[grouped["representation"].eq("Gene"), wide_index + value_columns].copy()
    pathway = grouped.loc[grouped["representation"].eq("Pathway"), wide_index + value_columns].copy()
    gene = gene.rename(columns={column: f"gene_{column}" for column in value_columns})
    pathway = pathway.rename(columns={column: f"pathway_{column}" for column in value_columns})
    wide = gene.merge(pathway, on=wide_index, how="inner", validate="one_to_one")
    wide = wide.loc[~wide["display_excluded_bool"]].copy()
    wide["scope"] = "internal"
    wide["dataset_or_direction"] = wide["dataset"]
    wide["perturbation_type_label"] = (
        wide[["gene_perturbation_type_label", "pathway_perturbation_type_label"]]
        .apply(combine_pair_annotations, axis=1)
    )
    wide["support_value"] = pd.concat(
        [
            pd.to_numeric(wide["gene_n_queries_total"], errors="coerce"),
            pd.to_numeric(wide["pathway_n_queries_total"], errors="coerce"),
        ],
        axis=1,
    ).max(axis=1, skipna=True)
    wide["support_label"] = wide["support_value"].map(aggregation_support_label)
    wide["pair_mean_enrichment"] = pd.concat(
        [
            pd.to_numeric(wide["gene_enrichment_score"], errors="coerce"),
            pd.to_numeric(wide["pathway_enrichment_score"], errors="coerce"),
        ],
        axis=1,
    ).mean(axis=1, skipna=True)
    wide["max_abs_enrichment"] = pd.concat(
        [
            pd.to_numeric(wide["gene_enrichment_score"], errors="coerce").abs(),
            pd.to_numeric(wide["pathway_enrichment_score"], errors="coerce").abs(),
        ],
        axis=1,
    ).max(axis=1, skipna=True)
    shared_entities = (
        wide.groupby(entity_column, dropna=False)["dataset"]
        .nunique()
        .ge(2)
        .rename("shared_across_modalities_bool")
        .reset_index()
    )
    wide = wide.merge(shared_entities, on=entity_column, how="left", validate="many_to_one")
    selected = select_dual_tail_exemplars(
        wide,
        entity_column=entity_column,
        dataset_column="dataset",
        score_column="pair_mean_enrichment",
        high_n=3,
        low_n=3,
        shared_column="shared_across_modalities_bool",
        min_shared=2,
    )
    selected = selected.sort_values(["dataset", "row_order"], kind="mergesort").reset_index(drop=True)
    selected["global_row_order"] = np.arange(1, len(selected) + 1, dtype=int)
    return write_csv(
        selected,
        output_path,
        [
            "dataset",
            entity_column,
            "row_id",
            "global_row_order",
            "row_order",
            "entity_display_label",
            "scope",
            "dataset_or_direction",
            "perturbation_type_label",
            "selection_direction",
            "selection_rank_within_direction",
            "selection_note",
            "shared_across_modalities_bool",
            "gene_n_queries_total",
            "pathway_n_queries_total",
            "gene_n_positive",
            "pathway_n_positive",
            "gene_n_negative",
            "pathway_n_negative",
            "gene_n_background_positive",
            "pathway_n_background_positive",
            "gene_n_background_negative",
            "pathway_n_background_negative",
            "gene_success_rate",
            "pathway_success_rate",
            "gene_background_success_rate",
            "pathway_background_success_rate",
            "gene_odds_ratio",
            "pathway_odds_ratio",
            "gene_raw_p",
            "pathway_raw_p",
            "gene_bh_q",
            "pathway_bh_q",
            "gene_significant_bool",
            "pathway_significant_bool",
            "gene_significance_label",
            "pathway_significance_label",
            "gene_enrichment_score",
            "pathway_enrichment_score",
            "pair_mean_enrichment",
            "max_abs_enrichment",
            "support_value",
            "support_label",
            "gene_n_source_strata",
            "pathway_n_source_strata",
        ],
    )


def build_task2_pattern_exemplar_panel(
    *,
    source_path: Path,
    entity_column: str,
    output_path: Path,
    keep_dataset: str | None,
    support_group_column: str,
    min_shared: int,
) -> Path:
    wide = build_task2_pattern_surface_wide(
        source_path=source_path,
        entity_column=entity_column,
        support_group_column=support_group_column,
    )
    if keep_dataset is not None:
        wide_for_selection = wide.loc[wide["dataset"].eq(keep_dataset)].copy()
    else:
        wide_for_selection = wide.copy()
    selected = select_ranked_entities(
        wide_for_selection,
        dataset_column="dataset",
        score_column="suitability_score",
        top_n=(25 if keep_dataset is not None else 10),
        shared_column="shared_flag",
        min_shared=min_shared,
    )
    selected = selected.sort_values(["dataset", "row_rank"], kind="mergesort").reset_index(drop=True)
    selected["global_row_order"] = np.arange(1, len(selected) + 1, dtype=int)
    selected["row_id"] = (
        selected["dataset"].astype(str)
        + "::"
        + selected[entity_column].astype(str)
    )
    return write_csv(
        selected,
        output_path,
        [
            "dataset",
            entity_column,
            "row_id",
            "global_row_order",
            "row_rank",
            "entity_display_label",
            "suitability_score",
            "shared_across_modalities_bool",
            "shared_flag",
            "gene_n_positive",
            "pathway_n_positive",
            "gene_n_negative",
            "pathway_n_negative",
            "gene_n_background_positive",
            "pathway_n_background_positive",
            "gene_n_background_negative",
            "pathway_n_background_negative",
            "gene_success_rate",
            "pathway_success_rate",
            "gene_background_success_rate",
            "pathway_background_success_rate",
            "gene_odds_ratio",
            "pathway_odds_ratio",
            "gene_raw_p",
            "pathway_raw_p",
            "gene_bh_q",
            "pathway_bh_q",
            "gene_significant_bool",
            "pathway_significant_bool",
            "gene_significance_label",
            "pathway_significance_label",
            "gene_enrichment_score",
            "pathway_enrichment_score",
            "pair_mean_enrichment",
            "max_abs_enrichment",
            "support_value",
            "support_n",
            "support_label",
            "gene_n_strata_included",
            "pathway_n_strata_included",
            "gene_analysis_family_set",
            "pathway_analysis_family_set",
            "gene_direction_set",
            "pathway_direction_set",
            "gene_metric_name_set",
            "pathway_metric_name_set",
        ],
    )


def build_task2_pattern_surface_wide(
    *,
    source_path: Path,
    entity_column: str,
    support_group_column: str,
) -> pd.DataFrame:
    required = {
        "dataset",
        entity_column,
        "analysis_family",
        "direction",
        "representation",
        "metric_name",
        support_group_column,
        "n_queries_total",
        "mean_metric_value",
        "pattern_summary_scope",
        "pattern_note",
        "fm_scope_note",
    }
    frame = read_csv_required(source_path, required, f"Task2 {entity_column} pattern summary").copy()
    frame["direction"] = frame["direction"].fillna("ALL")
    frame = frame.loc[
        frame["representation"].isin(TASK2_COMMON_BASE_REPRESENTATIONS)
        & frame["analysis_family"].isin(["group_concordance", "retrieval"])
        & ((frame["analysis_family"].ne("retrieval")) | frame["direction"].eq("C2G"))
    ].copy()
    frame = filter_excluded_metrics(frame)
    frame["metric_value"] = pd.to_numeric(frame["mean_metric_value"], errors="coerce")
    stratum_columns = ["dataset", "analysis_family", "direction", "metric_name"]
    valid_count = frame.groupby(stratum_columns, dropna=False)["metric_value"].transform(lambda values: values.notna().sum())
    success_cutoff = np.ceil(valid_count.astype("float64") * 0.25)
    rank_desc = frame.groupby(stratum_columns, dropna=False)["metric_value"].rank(
        method="first",
        ascending=False,
        na_option="keep",
    )
    frame["success_bool"] = rank_desc.le(success_cutoff).fillna(False)
    frame["n_positive"] = frame["success_bool"].astype(int)
    frame["n_negative"] = frame["metric_value"].notna().astype(int) - frame["n_positive"]
    frame["stratum_positive_total"] = frame.groupby(stratum_columns, dropna=False)["n_positive"].transform("sum")
    frame["stratum_valid_total"] = frame.groupby(stratum_columns, dropna=False)["metric_value"].transform(
        lambda values: values.notna().sum()
    )
    frame["n_background_positive"] = frame["stratum_positive_total"] - frame["n_positive"]
    frame["n_background_negative"] = (
        frame["stratum_valid_total"] - frame["stratum_positive_total"] - frame["n_negative"]
    )
    frame["support_component"] = np.where(
        frame["analysis_family"].eq("group_concordance"),
        pd.to_numeric(frame[support_group_column], errors="coerce"),
        pd.to_numeric(frame["n_queries_total"], errors="coerce"),
    )
    frame["entity_display_label"] = frame[entity_column].map(lambda value: entity_display_label(entity_column, value))
    grouped = (
        frame.groupby(["dataset", entity_column, "representation", "entity_display_label"], dropna=False, sort=True)
        .agg(
            n_positive=("n_positive", "sum"),
            n_negative=("n_negative", "sum"),
            n_background_positive=("n_background_positive", "sum"),
            n_background_negative=("n_background_negative", "sum"),
            support_value=("support_component", "max"),
            n_strata_included=("metric_name", "size"),
            analysis_family_set=("analysis_family", stringify_series),
            direction_set=("direction", stringify_series),
            metric_name_set=("metric_name", stringify_series),
        )
        .reset_index()
    )
    grouped = add_enrichment_statistics(grouped, adjust_group_columns=["dataset", "representation"])
    wide_index = ["dataset", entity_column, "entity_display_label"]
    value_columns = [
        "n_positive",
        "n_negative",
        "n_background_positive",
        "n_background_negative",
        "success_rate",
        "background_success_rate",
        "odds_ratio",
        "raw_p",
        "bh_q",
        "significant_bool",
        "significance_label",
        "enrichment_score",
        "support_value",
        "n_strata_included",
        "analysis_family_set",
        "direction_set",
        "metric_name_set",
    ]
    gene = grouped.loc[grouped["representation"].eq("Gene"), wide_index + value_columns].copy()
    pathway = grouped.loc[grouped["representation"].eq("Pathway"), wide_index + value_columns].copy()
    gene = gene.rename(columns={column: f"gene_{column}" for column in value_columns})
    pathway = pathway.rename(columns={column: f"pathway_{column}" for column in value_columns})
    wide = gene.merge(pathway, on=wide_index, how="inner", validate="one_to_one")
    wide["support_value"] = pd.concat(
        [
            pd.to_numeric(wide["gene_support_value"], errors="coerce"),
            pd.to_numeric(wide["pathway_support_value"], errors="coerce"),
        ],
        axis=1,
    ).max(axis=1, skipna=True)
    wide["support_label"] = wide["support_value"].map(aggregation_support_label)
    wide["pair_mean_enrichment"] = pd.concat(
        [
            pd.to_numeric(wide["gene_enrichment_score"], errors="coerce"),
            pd.to_numeric(wide["pathway_enrichment_score"], errors="coerce"),
        ],
        axis=1,
    ).mean(axis=1, skipna=True)
    wide["max_abs_enrichment"] = pd.concat(
        [
            pd.to_numeric(wide["gene_enrichment_score"], errors="coerce").abs(),
            pd.to_numeric(wide["pathway_enrichment_score"], errors="coerce").abs(),
        ],
        axis=1,
    ).max(axis=1, skipna=True)
    shared_entities = (
        wide.groupby(entity_column, dropna=False)["dataset"]
        .nunique()
        .ge(2)
        .rename("shared_across_modalities_bool")
        .reset_index()
    )
    wide = wide.merge(shared_entities, on=entity_column, how="left", validate="many_to_one")
    wide["shared_flag"] = wide["shared_across_modalities_bool"].fillna(False).astype(bool)
    wide["support_n"] = pd.to_numeric(wide["support_value"], errors="coerce").round().astype("Int64")
    wide["suitability_score"] = wide["pair_mean_enrichment"]
    return wide.sort_values(["dataset", entity_column], kind="mergesort").reset_index(drop=True)


def select_heatmap_entities(frame: pd.DataFrame, entity_column: str, top_n: int) -> pd.DataFrame:
    ranked = (
        frame.groupby([entity_column, "representation"], dropna=False, sort=True)["rank_within_family"]
        .min()
        .reset_index(name="best_rank")
    )
    entity_rank = (
        ranked.groupby(entity_column, dropna=False)
        .agg(best_rank=("best_rank", "min"))
        .sort_values(["best_rank", entity_column], kind="mergesort")
        .reset_index()
    )
    selected = entity_rank.head(top_n)[entity_column]
    return frame.loc[frame[entity_column].isin(selected)].copy()


def select_diverse_top_rows(
    frame: pd.DataFrame,
    *,
    family_column: str,
    score_column: str,
    top_n: int,
) -> pd.DataFrame:
    if frame.empty:
        out = frame.copy()
        out["family_rank"] = pd.Series(dtype="Int64")
        return out
    ordered = frame.sort_values(
        [family_column, score_column, "min_bh_q", "entity_display_label"],
        ascending=[True, False, True, True],
        kind="mergesort",
    ).copy()
    ordered["family_rank"] = ordered.groupby(family_column, dropna=False).cumcount() + 1
    return ordered.loc[ordered["family_rank"].le(top_n)].copy()


def paired_wilcoxon_from_long(
    frame: pd.DataFrame,
    *,
    group_columns: Sequence[str],
    unit_column: str | None = None,
    unit_columns: Sequence[str] | None = None,
    value_column: str,
    representation_column: str = "representation",
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if unit_columns is None:
        if unit_column is None:
            raise ValueError("paired_wilcoxon_from_long requires unit_column or unit_columns")
        normalized_unit_columns = [str(unit_column)]
    else:
        normalized_unit_columns = [str(column) for column in unit_columns]
    if frame.empty:
        return pd.DataFrame(
            columns=[
                *group_columns,
                "bh_q",
                "test_status",
                "effect_direction",
                "median_delta",
                "n_test_units",
            ]
        )
    for group_key, subset in frame.groupby(list(group_columns), dropna=False, sort=True):
        key_map = dict(zip(group_columns, group_key if isinstance(group_key, tuple) else (group_key,)))
        gene = subset.loc[subset[representation_column].eq("Gene"), [*normalized_unit_columns, value_column]].rename(
            columns={value_column: "gene_value"}
        )
        pathway = subset.loc[
            subset[representation_column].eq("Pathway"), [*normalized_unit_columns, value_column]
        ].rename(columns={value_column: "pathway_value"})
        paired = gene.merge(pathway, on=normalized_unit_columns, how="inner", validate="one_to_one")
        paired["delta"] = pd.to_numeric(paired["pathway_value"], errors="coerce") - pd.to_numeric(
            paired["gene_value"], errors="coerce"
        )
        paired = paired.loc[paired["delta"].notna()].copy()
        n_units = int(len(paired))
        if n_units < 2:
            rows.append(
                {
                    **key_map,
                    "bh_q": pd.NA,
                    "test_status": "not_tested_low_n",
                    "effect_direction": "undetermined",
                    "median_delta": pd.NA,
                    "n_test_units": n_units,
                }
            )
            continue
        delta = paired["delta"].astype("float64")
        nonzero = delta.loc[delta.ne(0)]
        if nonzero.empty:
            raw_p = 1.0
        else:
            _, raw_p = wilcoxon(delta, zero_method="wilcox", alternative="two-sided", mode="auto")
        median_delta = float(delta.median())
        rows.append(
            {
                **key_map,
                "raw_p": float(raw_p),
                "test_status": "tested",
                "effect_direction": (
                    "group_a_gt_group_b"
                    if median_delta > 0
                    else ("group_b_gt_group_a" if median_delta < 0 else "no_directional_shift")
                ),
                "median_delta": median_delta,
                "n_test_units": n_units,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["bh_q"] = pd.Series(pd.NA, index=out.index, dtype="Float64")
    valid_mask = out.get("raw_p", pd.Series(index=out.index, dtype="float64")).notna()
    if valid_mask.any():
        for _, idx in out.loc[valid_mask].groupby(
            [column for column in group_columns if column != "representation"], sort=False
        ).groups.items():
            q_values = bh_adjust(out.loc[idx, "raw_p"].astype("float64"))
            out.loc[idx, "bh_q"] = q_values.astype("Float64")
    if "raw_p" not in out.columns:
        out["raw_p"] = pd.NA
    return out


def build_task1_enrichment_pairs(
    *,
    source_path: Path,
    entity_column: str,
    exclusive_top_n: int,
    shared_top_n: int,
    selection_mode: str = "main",
    output_path: Path,
) -> Path:
    filtered = load_filtered_high_concordance(source_path, entity_column)
    filtered = filtered.loc[
        filtered["representation"].isin(["Gene", "Pathway"])
        & filtered["scope"].eq("internal")
    ].copy()
    filtered["dataset"] = filtered["dataset_or_direction"].astype(str)
    filtered["entity_display_label"] = filtered[entity_column].map(
        lambda value: entity_display_label(entity_column, value)
    )
    filtered["display_excluded_bool"] = filtered[entity_column].map(
        lambda value: entity_display_exclusion(entity_column, value)
    ).astype(bool)
    filtered["enrichment_score"] = log2_enrichment_score(filtered)
    q_values = pd.to_numeric(filtered["bh_q"], errors="coerce")
    filtered["significance_label"] = np.where(
        q_values.notna() & q_values.le(0.05),
        "*",
        "",
    )

    index_columns = [
        "dataset",
        entity_column,
        "entity_display_label",
        "display_excluded_bool",
        "scope",
        "perturbation_type",
    ]
    value_columns = [
        "n_queries_total",
        "n_positive",
        "n_negative",
        "success_rate",
        "background_success_rate",
        "odds_ratio",
        "raw_p",
        "bh_q",
        "significant_bool",
        "enrichment_score",
        "significance_label",
    ]
    gene = filtered.loc[filtered["representation"].eq("Gene"), index_columns + value_columns].copy()
    pathway = filtered.loc[filtered["representation"].eq("Pathway"), index_columns + value_columns].copy()
    gene = gene.rename(columns={column: f"gene_{column}" for column in value_columns})
    pathway = pathway.rename(columns={column: f"pathway_{column}" for column in value_columns})
    wide = gene.merge(pathway, on=index_columns, how="outer", validate="one_to_one")
    wide["pair_available_bool"] = wide["gene_enrichment_score"].notna() & wide["pathway_enrichment_score"].notna()
    wide = wide.loc[wide["pair_available_bool"]].copy()
    wide["support_value"] = pd.concat(
        [
            pd.to_numeric(wide["gene_n_queries_total"], errors="coerce"),
            pd.to_numeric(wide["pathway_n_queries_total"], errors="coerce"),
        ],
        axis=1,
    ).max(axis=1, skipna=True)
    wide["support_n"] = pd.to_numeric(wide["support_value"], errors="coerce").round().astype("Int64")
    wide["pair_mean_enrichment"] = pd.concat(
        [
            pd.to_numeric(wide["gene_enrichment_score"], errors="coerce"),
            pd.to_numeric(wide["pathway_enrichment_score"], errors="coerce"),
        ],
        axis=1,
    ).mean(axis=1, skipna=True)
    wide["max_abs_enrichment"] = pd.concat(
        [
            pd.to_numeric(wide["gene_enrichment_score"], errors="coerce").abs(),
            pd.to_numeric(wide["pathway_enrichment_score"], errors="coerce").abs(),
        ],
        axis=1,
    ).max(axis=1, skipna=True)
    wide["min_bh_q"] = pd.concat(
        [
            pd.to_numeric(wide["gene_bh_q"], errors="coerce"),
            pd.to_numeric(wide["pathway_bh_q"], errors="coerce"),
        ],
        axis=1,
    ).min(axis=1, skipna=True)
    wide = wide.loc[~wide["display_excluded_bool"]].copy()
    shared_entities = (
        wide.groupby([entity_column, "perturbation_type"], dropna=False)["dataset"]
        .nunique()
        .ge(2)
        .rename("shared_flag")
        .reset_index()
    )
    wide = wide.merge(shared_entities, on=[entity_column, "perturbation_type"], how="left", validate="many_to_one")
    wide["shared_flag"] = wide["shared_flag"].fillna(False).astype(bool)
    wide["shared_across_modalities_bool"] = wide["shared_flag"]
    wide["row_block"] = wide["perturbation_type"].astype(str)
    wide["col_block"] = wide["dataset"].astype(str)
    wide["block_id"] = wide["col_block"].astype(str) + "::" + wide["row_block"].astype(str)
    wide["scope_label"] = "Internal"
    wide["dataset_or_direction"] = wide["dataset"]
    wide["dataset_or_direction_label"] = wide["dataset"]
    wide["left_annotation_scope"] = wide["scope_label"]
    wide["left_annotation_dataset"] = wide["dataset_or_direction_label"]
    wide["left_annotation_perturbation"] = wide["perturbation_type"].astype(str)
    wide["perturbation_type_label"] = wide["perturbation_type"].astype(str)

    block_order = {"Chemical": 1, "Genetic": 2}
    col_order = {"LINCS": 1, "scPerturb": 2}
    if selection_mode == "main":
        selected = select_dual_tail_exemplars(
            wide,
            entity_column=entity_column,
            dataset_column="block_id",
            score_column="pair_mean_enrichment",
            high_n=3,
            low_n=3,
            shared_column="shared_flag",
            min_shared=1 if shared_top_n >= 0 else 0,
        )
        wide = selected.rename(columns={"selection_direction": "selection_tail"}).copy()
    elif selection_mode == "full":
        wide = wide.sort_values(
            ["row_block", "col_block", "pair_mean_enrichment", "max_abs_enrichment", "support_value", "entity_display_label"],
            ascending=[True, True, False, False, False, True],
            kind="mergesort",
        ).copy()
        wide["selection_tail"] = ""
        wide["selection_direction"] = ""
        wide["selection_rank_within_direction"] = pd.NA
        wide["selection_note"] = ""
        wide["row_order"] = (
            wide.groupby(["row_block", "col_block"], dropna=False).cumcount() + 1
        ).astype(int)
        wide["row_id"] = (
            wide["col_block"].astype(str)
            + "::"
            + wide["row_block"].astype(str)
            + "::"
            + wide[entity_column].astype(str)
        )
    else:
        raise ValueError(f"Unsupported Task1 enrichment selection_mode: {selection_mode}")

    wide["display_block_order"] = wide["row_block"].map(block_order).fillna(99).astype(int)
    wide["display_col_order"] = wide["col_block"].map(col_order).fillna(99).astype(int)
    wide = wide.sort_values(
        ["display_block_order", "display_col_order", "row_order", "entity_display_label"],
        ascending=[True, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    wide["display_rank_in_block"] = (
        wide.groupby(["row_block", "col_block"], dropna=False).cumcount() + 1
    ).astype(int)
    wide["selection_direction"] = wide["selection_tail"]
    wide["support_label"] = wide["support_value"].map(
        lambda value: f"n={format(int(round(value)), ',')}" if pd.notna(value) else ""
    )
    if "global_row_order" not in wide.columns:
        wide["global_row_order"] = np.arange(1, len(wide) + 1, dtype=int)
    else:
        wide["global_row_order"] = np.arange(1, len(wide) + 1, dtype=int)
    out_columns = [
        "dataset",
        entity_column,
        "row_id",
        "global_row_order",
        "row_order",
        "entity_display_label",
        "scope",
        "scope_label",
        "dataset_or_direction",
        "dataset_or_direction_label",
        "perturbation_type",
        "perturbation_type_label",
        "row_block",
        "col_block",
        "left_annotation_scope",
        "left_annotation_dataset",
        "left_annotation_perturbation",
        "gene_n_queries_total",
        "pathway_n_queries_total",
        "gene_n_positive",
        "pathway_n_positive",
        "gene_n_negative",
        "pathway_n_negative",
        "gene_success_rate",
        "pathway_success_rate",
        "gene_background_success_rate",
        "pathway_background_success_rate",
        "gene_odds_ratio",
        "pathway_odds_ratio",
        "gene_raw_p",
        "pathway_raw_p",
        "gene_bh_q",
        "pathway_bh_q",
        "gene_significant_bool",
        "pathway_significant_bool",
        "gene_significance_label",
        "pathway_significance_label",
        "gene_enrichment_score",
        "pathway_enrichment_score",
        "pair_mean_enrichment",
        "support_value",
        "support_n",
        "support_label",
        "max_abs_enrichment",
        "min_bh_q",
        "pair_available_bool",
        "shared_flag",
        "shared_across_modalities_bool",
        "selection_direction",
        "selection_tail",
        "display_rank_in_block",
        "selection_rank_within_direction",
        "selection_note",
    ]
    return write_csv(wide.drop(columns=["display_block_order"]), output_path, out_columns)


def task2_metric_block(analysis_family: object, direction: object, metric_name: object) -> str:
    family = str(analysis_family)
    metric_map = {
        "cosine_centroid": "Cosine",
        "pcc_centroid": "Pearson",
        "mrr_corrected": "MRR",
        "hit1_corrected": "Hit@1",
        "hit5_corrected": "Hit@5",
        "hit10_corrected": "Hit@10",
    }
    metric = metric_map.get(str(metric_name), str(metric_name))
    direction_text = str(direction) if pd.notna(direction) and str(direction) else "ALL"
    if family == "group_concordance":
        return f"Group | {metric}"
    return f"{direction_text} retrieval | {metric}"


def build_figure2_panel_2a(*, source_path: Path, output_path: Path) -> Path:
    frame = load_task1_scope_summary(source_path).copy()
    frame["representation_class"] = task1_scope_panel_axis(frame)
    frame["slice_label"] = np.where(
        frame["analysis_family"].eq("cross_eligibility"),
        "Cross eligibility | " + frame["perturbation_type"].astype(str),
        frame["dataset_or_direction"].astype(str) + " | " + frame["perturbation_type"].astype(str),
    )
    grouped = (
        frame.groupby(
            [
                "scope",
                "analysis_family",
                "slice_label",
                "dataset_or_direction",
                "perturbation_type",
                "representation_class",
            ],
            dropna=False,
            sort=True,
        )
        .agg(
            scope_status=("scope_status", collapsed_status_label),
            n_total=("n_total", max_or_na),
            n_valid=("n_valid", max_or_na),
            n_excluded=("n_excluded", max_or_na),
            n_matched_keys=("n_matched_keys", max_or_na),
            eligible_bool=("eligible_bool", any_true_bool),
            exclusion_reason=("exclusion_reason", stringify_series),
            scope_note=("scope_note", stringify_series),
            fm_scope_note=("fm_scope_note", stringify_series),
        )
        .reset_index()
    )
    grouped["coverage_denominator"] = [
        task1_scope_coverage_denominator(row.analysis_family, row.dataset_or_direction)
        for row in grouped.itertuples(index=False)
    ]
    grouped["coverage_numerator"] = grouped["scope_status"].map(status_available_numerator).astype(int)
    grouped["coverage_fraction"] = grouped["coverage_numerator"] / grouped["coverage_denominator"]
    grouped["bar_value"] = grouped["coverage_fraction"]
    grouped["count_annotation_type"] = np.where(
        grouped["analysis_family"].eq("cross_eligibility"),
        "n_matched_keys",
        "n_valid",
    )
    grouped["count_annotation"] = pd.to_numeric(
        np.where(
            grouped["analysis_family"].eq("cross_eligibility"),
            grouped["n_matched_keys"],
            grouped["n_valid"],
        ),
        errors="coerce",
    )
    out = grouped.sort_values(
        ["scope", "analysis_family", "dataset_or_direction", "perturbation_type", "representation_class"],
        kind="mergesort",
    )
    return write_csv(
        out,
        output_path,
        [
            "scope",
            "analysis_family",
            "slice_label",
            "dataset_or_direction",
            "perturbation_type",
            "representation_class",
            "scope_status",
            "bar_value",
            "coverage_denominator",
            "coverage_numerator",
            "coverage_fraction",
            "count_annotation",
            "count_annotation_type",
            "n_total",
            "n_valid",
            "n_excluded",
            "n_matched_keys",
            "eligible_bool",
            "exclusion_reason",
            "scope_note",
            "fm_scope_note",
        ],
    )


def load_task1_performance_structure(path: Path) -> pd.DataFrame:
    required = {
        "scope",
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "analysis_family",
        "metric_name",
        "value_variant",
        "value",
        "n_total",
        "n_valid",
        "n_excluded",
        "N_gallery_max",
        "cross_alignment_contract",
        "chance_check_available_bool",
        "chance_abs_delta",
        "fm_scope_note",
    }
    return read_csv_required(path, required, "Task1 performance structure")


def build_task1_performance_overview(*, source_path: Path, scope_value: str, output_path: Path) -> Path:
    frame = load_task1_performance_structure(source_path).copy()
    keep = (
        (frame["analysis_family"].eq("group_concordance") & frame["value_variant"].eq("group"))
        | (frame["analysis_family"].eq("retrieval") & frame["value_variant"].eq("retrieval_corrected"))
    )
    frame = frame.loc[frame["scope"].eq(scope_value) & keep].copy()
    if scope_value == "cross":
        frame = frame.loc[frame["perturbation_type"].eq("Genetic")].copy()
    frame = filter_excluded_metrics(frame)
    frame = percent_rank_within(
        frame,
        ["dataset_or_direction", "perturbation_type", "analysis_family", "metric_name"],
        "value",
        "fill_value",
    )
    out = frame.sort_values(
        ["dataset_or_direction", "perturbation_type", "analysis_family", "metric_name", "representation"],
        kind="mergesort",
    )
    return write_csv(
        out,
        output_path,
        [
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation",
            "analysis_family",
            "metric_name",
            "value_variant",
            "value",
            "fill_value",
            "n_total",
            "n_valid",
            "n_excluded",
            "N_gallery_max",
            "cross_alignment_contract",
            "chance_check_available_bool",
            "chance_abs_delta",
            "fm_scope_note",
        ],
    )


def build_task1_pathway_vs_gene_panel(
    *,
    bridge_detail_path: Path,
    comparison_stats_path: Path,
    output_path: Path,
) -> Path:
    raw = build_task1_bridge_detail_long(bridge_detail_path)
    raw = raw.loc[raw["representation_detail"].isin(["Gene", "Pathway"])].copy()
    raw = filter_excluded_metrics(raw)
    raw["comparison_context"] = np.where(
        raw["group_label"].eq("internal"),
        raw["dataset"].map({"LINCS": "Internal LINCS", "scPerturb": "Internal scPerturb"}),
        raw["dataset"].map(task1_cross_direction_label),
    )
    raw["representation"] = raw["representation_detail"]
    raw["unit_id"] = (
        raw["dataset"].astype(str)
        + "::"
        + raw["cell_line"].astype(str)
        + "::"
        + raw["target"].astype(str)
    )

    stats_frame = read_csv_required(
        comparison_stats_path,
        {
            "comparison_family",
            "analysis_family",
            "metric_name",
            "dataset_scope",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
        },
        "Manuscript comparison statistics",
    )
    internal_stats = stats_frame.loc[
        stats_frame["comparison_family"].eq("figure2_task1_internal_common_representation")
    ].copy()
    internal_stats["comparison_context"] = internal_stats["dataset_scope"].map(
        {"LINCS": "Internal LINCS", "scPerturb": "Internal scPerturb"}
    )
    cross_stats = stats_frame.loc[
        stats_frame["comparison_family"].eq("figure2_task1_cross_common_representation")
    ].copy()
    cross_stats["comparison_context"] = cross_stats["dataset_scope"].map(task1_cross_direction_label)
    stats_subset = pd.concat([internal_stats, cross_stats], ignore_index=True, sort=False)[
        [
            "comparison_context",
            "analysis_family",
            "metric_name",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
        ]
    ].drop_duplicates()
    stats_subset = filter_excluded_metrics(stats_subset)

    raw = raw.merge(
        stats_subset,
        on=["comparison_context", "analysis_family", "metric_name"],
        how="left",
        validate="many_to_one",
    )
    raw["row_kind"] = "raw"
    raw["mean_value"] = pd.NA
    raw["summary_mean"] = pd.NA
    raw["summary_median"] = pd.NA
    raw["q10_value"] = pd.NA
    raw["q25_value"] = pd.NA
    raw["q50_value"] = pd.NA
    raw["q75_value"] = pd.NA
    raw["q90_value"] = pd.NA
    raw["ci_low"] = pd.NA
    raw["ci_high"] = pd.NA
    raw["n_units"] = pd.NA
    raw["cross_direction"] = np.where(raw["group_label"].eq("cross"), raw["comparison_context"], pd.NA)

    summary = (
        raw.groupby(["comparison_context", "analysis_family", "metric_name", "representation"], dropna=False, sort=True)
        .agg(n_units=("unit_id", "nunique"))
        .reset_index()
    )
    summary_quantile_frame = summarize_quantiles_frame(
        raw,
        group_columns=["comparison_context", "analysis_family", "metric_name", "representation"],
        value_column="metric_value",
    )
    summary = summary.merge(
        summary_quantile_frame,
        on=["comparison_context", "analysis_family", "metric_name", "representation"],
        how="left",
        validate="one_to_one",
    )
    summary["mean_value"] = summary["summary_mean"]
    ci_bounds = raw.groupby(
        ["comparison_context", "analysis_family", "metric_name", "representation"],
        dropna=False,
        sort=True,
    )["metric_value"].apply(simple_ci_bounds)
    ci_frame = ci_bounds.reset_index(name="ci_bounds")
    ci_frame["ci_low"] = ci_frame["ci_bounds"].map(lambda value: value[0])
    ci_frame["ci_high"] = ci_frame["ci_bounds"].map(lambda value: value[1])
    ci_frame = ci_frame.drop(columns=["ci_bounds"])
    summary = summary.merge(
        ci_frame,
        on=["comparison_context", "analysis_family", "metric_name", "representation"],
        how="left",
        validate="one_to_one",
    ).merge(
        stats_subset,
        on=["comparison_context", "analysis_family", "metric_name"],
        how="left",
        validate="many_to_one",
    )
    summary["row_kind"] = "summary"
    summary["dataset"] = pd.NA
    summary["cell_line"] = pd.NA
    summary["target"] = pd.NA
    summary["direction"] = pd.NA
    summary["cross_direction"] = np.where(summary["comparison_context"].str.contains("->"), summary["comparison_context"], pd.NA)
    summary["unit_id"] = "summary"
    summary["metric_value"] = summary["q50_value"]

    combined = pd.concat(
        [
            raw[
                [
                    "row_kind",
                    "comparison_context",
                    "cross_direction",
                    "dataset",
                    "cell_line",
                    "target",
                    "direction",
                    "unit_id",
                    "analysis_family",
                    "metric_name",
                    "representation",
                    "metric_value",
                    "mean_value",
                    "summary_mean",
                    "summary_median",
                    "q10_value",
                    "q25_value",
                    "q50_value",
                    "q75_value",
                    "q90_value",
                    "ci_low",
                    "ci_high",
                    "n_units",
                    "bh_q",
                    "test_status",
                    "effect_direction",
                    "median_delta",
                    "n_test_units",
                ]
            ],
            summary[
                [
                    "row_kind",
                    "comparison_context",
                    "cross_direction",
                    "dataset",
                    "cell_line",
                    "target",
                    "direction",
                    "unit_id",
                    "analysis_family",
                    "metric_name",
                    "representation",
                    "metric_value",
                    "mean_value",
                    "summary_mean",
                    "summary_median",
                    "q10_value",
                    "q25_value",
                    "q50_value",
                    "q75_value",
                    "q90_value",
                    "ci_low",
                    "ci_high",
                    "n_units",
                    "bh_q",
                    "test_status",
                    "effect_direction",
                    "median_delta",
                    "n_test_units",
                ]
            ],
        ],
        ignore_index=True,
        sort=False,
    ).sort_values(
        ["analysis_family", "metric_name", "comparison_context", "representation", "row_kind", "unit_id"],
        kind="mergesort",
    )
    return write_csv(
        combined,
        output_path,
        [
            "row_kind",
            "comparison_context",
            "cross_direction",
            "dataset",
            "cell_line",
            "target",
            "direction",
            "unit_id",
            "analysis_family",
            "metric_name",
            "representation",
            "metric_value",
            "mean_value",
            "summary_mean",
            "summary_median",
            "q10_value",
            "q25_value",
            "q50_value",
            "q75_value",
            "q90_value",
            "ci_low",
            "ci_high",
            "n_units",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
        ],
    )


def build_figure2_internal_to_cross_degradation_panel(
    *,
    source_path: Path,
    comparison_stats_path: Path,
    output_path: Path,
) -> Path:
    required = {
        "dataset",
        "cell_line",
        "target",
        "representation_detail",
        "metric_family",
        "metric_name",
        "internal_value",
        "cross_value",
    }
    frame = read_csv_required(source_path, required, "Task1 internal-to-cross bridge detail").copy()
    frame = frame.loc[
        frame["representation_detail"].isin(["Gene", "Pathway"])
        & frame["internal_value"].notna()
        & frame["cross_value"].notna()
    ].copy()
    frame["analysis_family"] = frame["metric_family"].map(
        {
            "group_concordance": "group_concordance",
            "retrieval_corrected": "retrieval",
        }
    )
    frame["unit_id"] = (
        frame["dataset"].astype(str)
        + "::"
        + frame["cell_line"].astype(str)
        + "::"
        + frame["target"].astype(str)
    )
    summary_long = pd.concat(
        [
            frame.assign(scope="Internal", metric_value=frame["internal_value"]),
            frame.assign(scope="Cross", metric_value=frame["cross_value"]),
        ],
        ignore_index=True,
        sort=False,
    )
    summary_long = apply_metric_transform(
        summary_long.loc[summary_long["metric_value"].notna()].copy(),
        metric_column="metric_name",
        value_column="metric_value",
    )
    summary_long = summary_long.rename(columns={"representation_detail": "representation"})
    summary_long["family_block"] = summary_long["analysis_family"].map(
        {
            "group_concordance": "Group",
            "retrieval": "Retrieval",
        }
    )
    summary = (
        summary_long.groupby(
            ["dataset", "family_block", "analysis_family", "metric_name", "representation", "scope", "value_transform"],
            dropna=False,
            sort=True,
        )
        .agg(n_units=("unit_id", "nunique"))
        .reset_index()
    )
    quantiles = summarize_quantiles_frame(
        summary_long,
        group_columns=["dataset", "family_block", "analysis_family", "metric_name", "representation", "scope", "value_transform"],
        value_column="metric_value",
    )
    summary = summary.merge(
        quantiles,
        on=["dataset", "family_block", "analysis_family", "metric_name", "representation", "scope", "value_transform"],
        how="left",
        validate="one_to_one",
    )
    comparison = read_csv_required(
        comparison_stats_path,
        {
            "comparison_family",
            "dataset_scope",
            "representation_scope",
            "analysis_family",
            "metric_name",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
        },
        "Manuscript comparison statistics",
    )
    comparison = comparison.loc[
        comparison["comparison_family"].eq("figure2_task1_internal_to_cross_degradation")
    ].copy()
    comparison = comparison.rename(
        columns={
            "dataset_scope": "dataset",
            "representation_scope": "representation",
        }
    )[
        [
            "dataset",
            "representation",
            "analysis_family",
            "metric_name",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
        ]
    ].drop_duplicates()
    summary = summary.merge(
        comparison,
        on=["dataset", "representation", "analysis_family", "metric_name"],
        how="left",
        validate="many_to_one",
    )
    summary["degradation_scope_note"] = "Direct internal-versus-cross paired comparison over shared Task1 bridge triplets."
    summary["fm_scope_note"] = "FM excluded because Task1 cross excludes FM."
    out = summary.sort_values(
        ["family_block", "dataset", "metric_name", "representation", "scope"],
        kind="mergesort",
    )
    return write_csv(
        out,
        output_path,
        [
            "dataset",
            "family_block",
            "analysis_family",
            "metric_name",
            "representation",
            "scope",
            "value_transform",
            "summary_mean",
            "summary_median",
            "q10_value",
            "q25_value",
            "q50_value",
            "q75_value",
            "q90_value",
            "n_units",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
            "degradation_scope_note",
            "fm_scope_note",
        ],
    )


def build_figure2_panel_2c(
    *,
    bridge_detail_path: Path,
    comparison_stats_path: Path,
    output_path: Path,
) -> Path:
    return build_task1_pathway_vs_gene_panel(
        bridge_detail_path=bridge_detail_path,
        comparison_stats_path=comparison_stats_path,
        output_path=output_path,
    )


def build_figure2_panel_2d(
    *,
    source_path: Path,
    comparison_stats_path: Path,
    output_path: Path,
) -> Path:
    return build_figure2_internal_to_cross_degradation_panel(
        source_path=source_path,
        comparison_stats_path=comparison_stats_path,
        output_path=output_path,
    )


def build_figure2_panel_2e(*, source_path: Path, output_path: Path) -> Path:
    return build_task1_enrichment_pairs(
        source_path=source_path,
        entity_column="cell_line",
        exclusive_top_n=0,
        shared_top_n=0,
        selection_mode="main",
        output_path=output_path,
    )


def build_figure2_panel_2f(*, source_path: Path, output_path: Path) -> Path:
    return build_task1_enrichment_pairs(
        source_path=source_path,
        entity_column="target_token",
        exclusive_top_n=0,
        shared_top_n=0,
        selection_mode="main",
        output_path=output_path,
    )


def build_figure2_panel_2g(
    *,
    bridge_detail_path: Path,
    comparison_stats_path: Path,
    output_path: Path,
) -> Path:
    return build_task1_pathway_vs_gene_panel(
        bridge_detail_path=bridge_detail_path,
        comparison_stats_path=comparison_stats_path,
        output_path=output_path,
    )


def extended_high_concordance_full(frame: pd.DataFrame, entity_column: str) -> pd.DataFrame:
    out = frame.copy()
    out["family_id"] = out[FIGURE2_FAMILY_COLUMNS].astype(str).agg("|".join, axis=1)
    out = out.sort_values(
        ["family_id", "odds_ratio", "success_rate", entity_column],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    out["rank_within_family"] = out.groupby("family_id", dropna=False).cumcount().add(1).astype("Int64")
    out["label_bool"] = out["rank_within_family"].le(10)
    return out


def build_extended_figure4_cell_line(*, source_path: Path, output_path: Path) -> Path:
    return build_task1_enrichment_pairs(
        source_path=source_path,
        entity_column="cell_line",
        exclusive_top_n=0,
        shared_top_n=0,
        selection_mode="full",
        output_path=output_path,
    )


def build_extended_figure4_target(*, source_path: Path, output_path: Path) -> Path:
    return build_task1_enrichment_pairs(
        source_path=source_path,
        entity_column="target_token",
        exclusive_top_n=0,
        shared_top_n=0,
        selection_mode="full",
        output_path=output_path,
    )


def build_panel_comparison_extract(
    *,
    source_path: Path,
    comparison_families: Sequence[str],
    output_path: Path,
) -> Path:
    required = {
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
        "primary_source_object",
        "dataset_scope",
        "direction_scope",
    }
    frame = read_csv_required(source_path, required, "Manuscript comparison statistics").copy()
    out = frame.loc[frame["comparison_family"].isin(comparison_families)].copy()
    out = out.sort_values(
        ["comparison_family", "analysis_family", "metric_name", "dataset_scope", "direction_scope"],
        kind="mergesort",
    )
    return write_csv(out, output_path, list(frame.columns))


def collapse_task2_scope_by_class(frame: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        frame.assign(representation_class=frame["representation"].map(representation_class))
        .groupby(["dataset", "cell_line", "representation_class"], dropna=False, sort=True)
        .agg(
            availability_status=(
                "availability_status",
                lambda values: "available" if pd.Series(values).isin(["available", "materialized"]).any() else str(pd.Series(values).iloc[0]),
            ),
            availability_reason=("availability_reason", stringify_series),
            n_targets_eligible=("n_targets_eligible", max_or_na),
            n_targets_ineligible=("n_targets_ineligible", max_or_na),
            group_slice_materialized_bool=("group_slice_materialized_bool", "any"),
            retrieval_c2g_materialized_bool=("retrieval_c2g_materialized_bool", "any"),
            retrieval_g2c_materialized_bool=("retrieval_g2c_materialized_bool", "any"),
            retrieval_directions_present=("retrieval_directions_present", stringify_series),
            n_targets_total=("n_targets_total", max_or_na),
            n_attrition_target_rows=("n_attrition_target_rows", max_or_na),
            group_attrition_rows=("group_attrition_rows", max_or_na),
            scope_note=("scope_note", stringify_series),
            fm_scope_note=("fm_scope_note", stringify_series),
        )
        .reset_index()
    )
    return grouped


def build_figure3_panel_3a(*, source_path: Path, output_path: Path) -> Path:
    required = {
        "dataset",
        "cell_line",
        "representation",
        "availability_status",
        "availability_reason",
        "n_targets_eligible",
        "n_targets_ineligible",
        "group_slice_materialized_bool",
        "retrieval_c2g_materialized_bool",
        "retrieval_g2c_materialized_bool",
        "retrieval_directions_present",
        "n_targets_total",
        "n_attrition_target_rows",
        "group_attrition_rows",
        "scope_note",
        "fm_scope_note",
    }
    frame = read_csv_required(source_path, required, "Task2 scope summary")
    collapsed = collapse_task2_scope_by_class(frame)
    eligible_counts = (
        collapsed.loc[collapsed["availability_status"].isin(["available", "materialized"])]
        .groupby(["dataset", "cell_line"], dropna=False, sort=True)["n_targets_eligible"]
        .max()
        .reset_index(name="eligible_target_count")
    )
    cohort_template = collapsed[["dataset", "cell_line"]].drop_duplicates().copy()
    cohort_template["_key"] = 1
    class_template = pd.DataFrame(
        {
            "representation_class": list(TASK2_FROZEN_REPRESENTATION_CLASSES),
            "_key": 1,
        }
    )
    out = cohort_template.merge(class_template, on="_key", how="inner").drop(columns="_key")
    out = out.merge(
        collapsed,
        on=["dataset", "cell_line", "representation_class"],
        how="left",
        validate="one_to_one",
    )
    out["availability_status"] = out["availability_status"].fillna("not_applicable_scope")
    out["status_group"] = np.where(
        out["availability_status"].isin(["available", "materialized"]),
        "available",
        "not_applicable_scope",
    )
    out["coverage_denominator"] = 3
    out["coverage_numerator"] = out["status_group"].eq("available").astype(int)
    out["coverage_fraction"] = out["coverage_numerator"] / out["coverage_denominator"]
    out["representation_count"] = out["coverage_numerator"]
    out = out.merge(eligible_counts, on=["dataset", "cell_line"], how="left", validate="many_to_one")
    out["count_annotation"] = out["eligible_target_count"]
    out["count_annotation_type"] = "eligible_target_count"
    return write_csv(
        out,
        output_path,
        [
            "dataset",
            "cell_line",
            "representation_class",
            "status_group",
            "representation_count",
            "eligible_target_count",
            "coverage_denominator",
            "coverage_numerator",
            "coverage_fraction",
            "count_annotation",
            "count_annotation_type",
            "availability_reason",
            "scope_note",
            "fm_scope_note",
        ],
    )


def direction_support_long(frame: pd.DataFrame) -> pd.DataFrame:
    metric_columns = {
        "mean_mrr_corrected": "mean_mrr_corrected",
        "mean_hit1_corrected": "mean_hit1_corrected",
        "mean_hit5_corrected": "mean_hit5_corrected",
        "mean_hit10_corrected": "mean_hit10_corrected",
    }
    rows: list[pd.DataFrame] = []
    for metric_name, column in metric_columns.items():
        rows.append(
            pd.DataFrame(
                {
                    "analysis_family": "retrieval",
                    "dataset": frame["dataset"],
                    "cell_line": frame["cell_line"],
                    "direction": frame["direction"],
                    "representation": frame["representation"],
                    "metric_name": metric_name,
                    "metric_value": frame[column],
                    "n_total": frame["n_total"],
                    "n_valid": frame["n_valid"],
                    "n_excluded_missing_metric_or_mpos0": frame["n_excluded_missing_metric_or_mpos0"],
                    "N_gallery_mean": frame["N_gallery_mean"],
                    "N_gallery_max": frame["N_gallery_max"],
                    "m_pos_mean": frame["m_pos_mean"],
                    "m_pos_p50": frame["m_pos_p50"],
                    "m_pos_p90": frame["m_pos_p90"],
                    "source_note": "C2G retrieval row from figure3_task2_direction_support_summary.csv.",
                    "fm_scope_note": frame["fm_scope_note"],
                }
            )
        )
    return pd.concat(rows, ignore_index=True, sort=False)


def scperturb_k562_tier_retrieval_long(
    *,
    task2_retrieval_per_query_path: Path,
    drug_meta_path: Path,
) -> pd.DataFrame:
    require_path(task2_retrieval_per_query_path, "Task2 retrieval per-query parquet")
    require_path(drug_meta_path, "Drug_meta.csv")
    if pq is None:
        raise ModuleNotFoundError("pyarrow is required to read Task2 retrieval parquet inputs.")

    available_columns = set(pq.ParquetFile(task2_retrieval_per_query_path).schema.names)
    has_query_specificity_tier = "query_specificity_tier" in available_columns
    parquet_columns = [
        "dataset",
        "cell_line",
        "direction",
        "representation",
        "query_row_id",
        "treated_cell_id",
        "N_gallery",
        "m_pos",
        *TASK2_RETRIEVAL_VALUE_COLUMNS,
    ]
    if has_query_specificity_tier:
        parquet_columns.insert(6, "query_specificity_tier")

    per_query = pd.read_parquet(
        task2_retrieval_per_query_path,
        columns=parquet_columns,
    )
    per_query = per_query.loc[
        per_query["dataset"].eq(FM_LOCAL_DATASET)
        & per_query["cell_line"].eq(FM_LOCAL_CELL_LINE)
        & per_query["direction"].eq("C2G")
    ].copy()
    if per_query.empty:
        return pd.DataFrame(
            columns=[
                "analysis_family",
                "dataset",
                "cell_line",
                "direction",
                "specificity_tier",
                "representation",
                "metric_name",
                "metric_value",
                "n_targets_total",
                "n_targets_metric_valid",
                "n_total",
                "n_valid",
                "n_excluded_missing_metric_or_mpos0",
                "N_gallery_mean",
                "N_gallery_max",
                "m_pos_mean",
                "m_pos_p50",
                "m_pos_p90",
                "source_note",
                "fm_scope_note",
            ]
        )

    drug_meta = pd.read_csv(drug_meta_path)
    if "treated_cell_id" not in drug_meta.columns and "Unnamed: 0" in drug_meta.columns:
        drug_meta = drug_meta.rename(columns={"Unnamed: 0": "treated_cell_id"})
    required_columns = {"treated_cell_id", "benchmark_group", "specificity_tier"}
    missing_columns = sorted(required_columns - set(drug_meta.columns))
    if missing_columns:
        raise ValueError(f"Drug_meta.csv is missing required columns: {missing_columns}")
    drug_meta = drug_meta.copy()
    drug_meta["treated_cell_id"] = drug_meta["treated_cell_id"].astype(str)
    drug_meta["benchmark_group"] = drug_meta["benchmark_group"].map(lambda value: "" if pd.isna(value) else str(value))
    drug_meta["specificity_tier"] = drug_meta["specificity_tier"].map(lambda value: "" if pd.isna(value) else str(value))
    if drug_meta["treated_cell_id"].duplicated().any():
        dupes = drug_meta.loc[drug_meta["treated_cell_id"].duplicated(), "treated_cell_id"].head(5).tolist()
        raise ValueError(f"Drug_meta.csv treated_cell_id must be unique, examples={dupes}")

    per_query["treated_cell_id"] = per_query["treated_cell_id"].astype(str)
    if has_query_specificity_tier:
        per_query["query_specificity_tier"] = per_query["query_specificity_tier"].map(
            lambda value: "" if pd.isna(value) else str(value)
        )
    else:
        per_query["query_specificity_tier"] = ""
    joined = per_query.merge(
        drug_meta[["treated_cell_id", "benchmark_group", "specificity_tier"]],
        on="treated_cell_id",
        how="left",
        validate="many_to_one",
    )
    if joined["benchmark_group"].isna().any():
        examples = joined.loc[joined["benchmark_group"].isna(), "treated_cell_id"].head(5).tolist()
        raise ValueError(f"Drug_meta.csv join missing benchmark_group for treated_cell_id examples={examples}")
    if joined["specificity_tier"].isna().any():
        examples = joined.loc[joined["specificity_tier"].isna(), "treated_cell_id"].head(5).tolist()
        raise ValueError(f"Drug_meta.csv join missing specificity_tier for treated_cell_id examples={examples}")

    joined = joined.loc[~joined["benchmark_group"].eq("Control")].copy()
    if joined.empty:
        return pd.DataFrame(
            columns=[
                "analysis_family",
                "dataset",
                "cell_line",
                "direction",
                "specificity_tier",
                "representation",
                "metric_name",
                "metric_value",
                "n_targets_total",
                "n_targets_metric_valid",
                "n_total",
                "n_valid",
                "n_excluded_missing_metric_or_mpos0",
                "N_gallery_mean",
                "N_gallery_max",
                "m_pos_mean",
                "m_pos_p50",
                "m_pos_p90",
                "source_note",
                "fm_scope_note",
            ]
        )

    tier_match = joined["query_specificity_tier"].eq(joined["specificity_tier"])
    if has_query_specificity_tier and not bool(tier_match.all()):
        examples = joined.loc[~tier_match, ["treated_cell_id", "query_specificity_tier", "specificity_tier"]]
        raise ValueError(
            "Drug_meta.csv specificity_tier mismatch against per-query query_specificity_tier, "
            f"examples={examples.head(5).to_dict(orient='records')}"
        )

    metric_columns = {
        "mean_mrr_corrected": "mrr_corrected",
        "mean_hit1_corrected": "hit1_corrected",
        "mean_hit5_corrected": "hit5_corrected",
        "mean_hit10_corrected": "hit10_corrected",
    }
    long_rows: list[pd.DataFrame] = []
    for metric_name, metric_column in metric_columns.items():
        long_rows.append(
            pd.DataFrame(
                {
                    "dataset": joined["dataset"],
                    "cell_line": joined["cell_line"],
                    "direction": joined["direction"],
                    "representation": joined["representation"],
                    "query_row_id": joined["query_row_id"],
                    "treated_cell_id": joined["treated_cell_id"],
                    "specificity_tier": joined["specificity_tier"],
                    "N_gallery": joined["N_gallery"],
                    "m_pos": joined["m_pos"],
                    "metric_name": metric_name,
                    "metric_value": joined[metric_column],
                }
            )
        )
    long_frame = pd.concat(long_rows, ignore_index=True, sort=False)
    summary = (
        long_frame.groupby(
            ["representation", "specificity_tier", "metric_name"],
            dropna=False,
            sort=True,
        )
        .agg(
            metric_value=("metric_value", "mean"),
            n_total=("query_row_id", "size"),
            n_valid=("metric_value", "count"),
            n_excluded_missing_metric_or_mpos0=("metric_value", lambda values: int(values.isna().sum())),
            N_gallery_mean=("N_gallery", "mean"),
            N_gallery_max=("N_gallery", "max"),
            m_pos_mean=("m_pos", "mean"),
            m_pos_p50=("m_pos", "median"),
            m_pos_p90=("m_pos", lambda values: values.quantile(0.9)),
        )
        .reset_index()
    )
    summary["analysis_family"] = "retrieval"
    summary["dataset"] = FM_LOCAL_DATASET
    summary["cell_line"] = FM_LOCAL_CELL_LINE
    summary["direction"] = "C2G"
    summary["n_targets_total"] = pd.NA
    summary["n_targets_metric_valid"] = pd.NA
    summary["source_note"] = (
        "C2G retrieval row from task2_retrieval_per_query.parquet joined to Drug_meta.csv by treated_cell_id."
    )
    summary["fm_scope_note"] = "scPerturb/K562 local tier-aware retrieval rows; Control benchmark_group excluded."
    return summary


def build_figure3_panel_3b(
    *,
    performance_path: Path,
    direction_support_path: Path,
    task2_retrieval_per_query_path: Path,
    drug_meta_path: Path,
    output_path: Path,
) -> Path:
    ordered_columns = [
        "analysis_family",
        "dataset",
        "cell_line",
        "direction",
        "specificity_tier",
        "representation",
        "metric_name",
        "metric_value",
        "n_units",
        "n_queries",
        "n_targets_total",
        "n_targets_metric_valid",
        "n_total",
        "n_valid",
        "n_excluded_missing_metric_or_mpos0",
        "N_gallery_mean",
        "N_gallery_max",
        "m_pos_mean",
        "m_pos_p50",
        "m_pos_p90",
        "source_note",
        "fm_scope_note",
    ]
    performance = read_csv_required(
        performance_path,
        {
            "analysis_family",
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "metric_name",
            "metric_value",
            "n_targets_total",
            "n_targets_metric_valid",
            "performance_scope_note",
            "fm_scope_note",
        },
        "Task2 performance structure",
    ).copy()
    group_rows = performance.loc[
        performance["dataset"].eq("LINCS")
        & performance["analysis_family"].eq("group_concordance")
        & performance["representation"].isin(TASK2_COMMON_BASE_REPRESENTATIONS)
    ].copy()
    group_rows = filter_excluded_metrics(group_rows)
    group_rows["direction"] = group_rows["direction"].fillna("ALL")
    group_rows["specificity_tier"] = pd.NA
    group_rows["n_units"] = group_rows["n_targets_metric_valid"]
    group_rows["n_queries"] = pd.Series(pd.array([pd.NA] * len(group_rows), dtype="Float64"), index=group_rows.index)
    group_rows["source_note"] = group_rows["performance_scope_note"]
    for column in [
        "n_total",
        "n_valid",
        "n_excluded_missing_metric_or_mpos0",
        "N_gallery_mean",
        "N_gallery_max",
        "m_pos_mean",
        "m_pos_p50",
        "m_pos_p90",
    ]:
        group_rows[column] = pd.Series(pd.array([pd.NA] * len(group_rows), dtype="Float64"), index=group_rows.index)
    group_rows = group_rows[ordered_columns]

    direction_support = read_csv_required(
        direction_support_path,
        {
            "dataset",
            "cell_line",
            "direction",
            "representation",
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
            "fm_scope_note",
        },
        "Task2 direction support summary",
    )
    direction_support = direction_support.loc[
        direction_support["dataset"].eq("LINCS")
        & direction_support["direction"].eq("C2G")
        & direction_support["representation"].isin(TASK2_COMMON_BASE_REPRESENTATIONS)
    ].copy()
    retrieval_rows = direction_support_long(direction_support)
    retrieval_rows["specificity_tier"] = pd.NA
    retrieval_rows["n_units"] = pd.Series(
        pd.array([pd.NA] * len(retrieval_rows), dtype="Float64"),
        index=retrieval_rows.index,
    )
    retrieval_rows["n_queries"] = retrieval_rows["n_total"]
    for column in ["n_targets_total", "n_targets_metric_valid"]:
        retrieval_rows[column] = pd.Series(
            pd.array([pd.NA] * len(retrieval_rows), dtype="Float64"),
            index=retrieval_rows.index,
        )
    retrieval_rows = retrieval_rows[ordered_columns]

    scperturb_retrieval_rows = scperturb_k562_tier_retrieval_long(
        task2_retrieval_per_query_path=task2_retrieval_per_query_path,
        drug_meta_path=drug_meta_path,
    )
    scperturb_retrieval_rows = scperturb_retrieval_rows.loc[
        scperturb_retrieval_rows["representation"].isin(TASK2_COMMON_BASE_REPRESENTATIONS)
    ].copy()
    scperturb_retrieval_rows["n_units"] = pd.Series(
        pd.array([pd.NA] * len(scperturb_retrieval_rows), dtype="Float64"),
        index=scperturb_retrieval_rows.index,
    )
    scperturb_retrieval_rows["n_queries"] = scperturb_retrieval_rows["n_total"]
    scperturb_retrieval_rows["n_targets_total"] = pd.Series(
        pd.array([pd.NA] * len(scperturb_retrieval_rows), dtype="Float64"),
        index=scperturb_retrieval_rows.index,
    )
    scperturb_retrieval_rows["n_targets_metric_valid"] = pd.Series(
        pd.array([pd.NA] * len(scperturb_retrieval_rows), dtype="Float64"),
        index=scperturb_retrieval_rows.index,
    )
    scperturb_retrieval_rows = scperturb_retrieval_rows[ordered_columns]

    combined = pd.DataFrame.from_records(
        group_rows.to_dict(orient="records")
        + retrieval_rows.to_dict(orient="records")
        + scperturb_retrieval_rows.to_dict(orient="records")
    )
    combined = combined.reindex(columns=ordered_columns)
    combined = filter_excluded_metrics(combined)
    combined = percent_rank_within(
        combined,
        ["dataset", "analysis_family", "representation", "metric_name"],
        "metric_value",
        "fill_value",
    )
    out = combined.sort_values(
        ["dataset", "analysis_family", "representation", "metric_name", "cell_line", "specificity_tier"],
        kind="mergesort",
    )
    return write_csv(
        out,
        output_path,
        [
            *ordered_columns[:7],
            "fill_value",
            *ordered_columns[7:],
        ],
    )


def load_target_pattern(path: Path) -> pd.DataFrame:
    required = {
        "dataset",
        "target",
        "analysis_family",
        "direction",
        "representation",
        "metric_name",
        "n_cell_lines",
        "n_source_rows",
        "mean_metric_value",
        "n_queries_total",
        "fm_scope_note",
    }
    frame = read_csv_required(path, required, "Figure 3 target pattern summary")
    frame = frame.copy()
    frame["direction"] = frame["direction"].fillna("ALL")
    frame["metric_value"] = frame["mean_metric_value"]
    return filter_excluded_metrics(frame)


def build_figure3_panel_3d(*, source_path: Path, output_path: Path) -> Path:
    return build_task2_pattern_exemplar_panel(
        source_path=source_path,
        entity_column="target",
        output_path=output_path,
        keep_dataset=None,
        support_group_column="n_cell_lines",
        min_shared=2,
    )


def build_extended_figure6_target(*, source_path: Path, output_path: Path) -> Path:
    frame = load_target_pattern(source_path)
    frame = frame.loc[
        frame["analysis_family"].isin(["group_concordance", "retrieval"])
        & (
            frame["analysis_family"].ne("retrieval")
            | frame["direction"].eq("C2G")
        )
    ].copy()
    frame = dense_rank_desc(
        frame,
        ["dataset", "analysis_family", "direction", "representation", "metric_name"],
        "metric_value",
        "rank_within_representation",
    )
    frame = percentile_from_rank(
        frame,
        ["dataset", "analysis_family", "direction", "representation", "metric_name"],
        "rank_within_representation",
        "rank_percentile",
    )
    frame["fm_local_bool"] = frame["representation"].map(fm_local_bool)
    out = frame.sort_values(
        ["dataset", "analysis_family", "direction", "representation", "metric_name", "rank_within_representation", "target"],
        kind="mergesort",
    )
    expected = list(frame.columns) + ["fm_local_bool"]
    expected = [column for column in expected if column in out.columns]
    return write_csv(
        out,
        output_path,
        [
            "dataset",
            "target",
            "analysis_family",
            "direction",
            "representation",
            "metric_name",
            "n_cell_lines",
            "n_source_rows",
            "mean_metric_value",
            "median_metric_value",
            "min_metric_value",
            "max_metric_value",
            "n_queries_total",
            "n_queries_mean",
            "n_chem_instances_used_total",
            "n_gen_instances_used_total",
            "n_chem_sub_total",
            "n_gen_sub_total",
            "pattern_summary_scope",
            "pattern_note",
            "fm_scope_note",
            "metric_value",
            "rank_within_representation",
            "rank_percentile",
            "fm_local_bool",
        ],
    )


def select_family_diverse_targets(
    frame: pd.DataFrame,
    *,
    n: int,
    score_column: str,
    family_column: str,
    excluded_targets: set[str] | None = None,
) -> pd.DataFrame:
    excluded_targets = excluded_targets or set()
    if frame.empty or n <= 0:
        return frame.head(0).copy()
    available = frame.loc[~frame["target"].astype(str).isin(excluded_targets)].copy()
    available = available.sort_values([score_column, "target"], ascending=[False, True], kind="mergesort")
    selected_indexes: list[int] = []
    used_families: set[str] = set()
    used_targets: set[str] = set()
    for row in available.itertuples():
        if row.target in used_targets or getattr(row, family_column) in used_families:
            continue
        selected_indexes.append(int(row.Index))
        used_targets.add(str(row.target))
        used_families.add(str(getattr(row, family_column)))
        if len(selected_indexes) >= n:
            break
    if len(selected_indexes) < n:
        for row in available.itertuples():
            if row.target in used_targets:
                continue
            selected_indexes.append(int(row.Index))
            used_targets.add(str(row.target))
            if len(selected_indexes) >= n:
                break
    out = available.loc[selected_indexes].copy()
    out = out.sort_values([score_column, "target"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    out["family_diversity_rank"] = np.arange(1, len(out) + 1, dtype=int)
    return out


def build_figure3_panel_3c(*, source_path: Path, output_path: Path) -> Path:
    return build_task2_pattern_exemplar_panel(
        source_path=source_path,
        entity_column="cell_line",
        output_path=output_path,
        keep_dataset="LINCS",
        support_group_column="n_targets",
        min_shared=0,
    )


def build_figure3_panel_3e(
    *,
    task2_group_root: Path,
    task2_retrieval_per_query_path: Path,
    comparison_stats_path: Path,
    output_path: Path,
) -> Path:
    performance = build_task2_target_support_long(
        task2_group_root,
        task2_retrieval_per_query_path.parent,
    ).copy()
    performance["direction"] = performance["direction"].fillna("ALL")
    performance = performance.loc[
        performance["representation"].isin(["Gene", "Pathway"])
        & performance["metric_value"].notna()
        & (
            (
                performance["analysis_family"].eq("group_concordance")
                & performance["direction"].eq("ALL")
            )
            | (
                performance["analysis_family"].eq("retrieval")
                & performance["direction"].eq("C2G")
            )
        )
    ].copy()
    performance = filter_excluded_metrics(performance)
    performance = apply_metric_transform(
        performance,
        metric_column="metric_name",
        value_column="metric_value",
    )
    performance["unit_source"] = "target_level_performance"
    performance["unit_id"] = (
        performance["dataset"].astype(str)
        + "::"
        + performance["cell_line"].astype(str)
        + "::"
        + performance["target"].astype(str)
    )

    raw = performance.copy()
    raw["row_kind"] = "raw"
    raw["mean_value"] = pd.NA
    raw["ci_low"] = pd.NA
    raw["ci_high"] = pd.NA
    raw["n_units"] = pd.NA
    raw["summary_mean"] = pd.NA
    raw["summary_median"] = pd.NA
    raw["q10_value"] = pd.NA
    raw["q25_value"] = pd.NA
    raw["q50_value"] = pd.NA
    raw["q75_value"] = pd.NA
    raw["q90_value"] = pd.NA

    stats_frame = read_csv_required(
        comparison_stats_path,
        {
            "comparison_family",
            "dataset_scope",
            "analysis_family",
            "direction_scope",
            "metric_name",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
        },
        "Manuscript comparison statistics",
    )
    stats_subset = (
        stats_frame.loc[
            stats_frame["comparison_family"].eq("figure3_task2_performance_target_level_representation")
            & (
                (
                    stats_frame["analysis_family"].eq("group_concordance")
                    & stats_frame["direction_scope"].eq("ALL")
                )
                | (
                    stats_frame["analysis_family"].eq("retrieval")
                    & stats_frame["direction_scope"].eq("C2G")
                )
            )
        ][
            [
                "dataset_scope",
                "analysis_family",
                "direction_scope",
                "metric_name",
                "bh_q",
                "test_status",
                "effect_direction",
                "median_delta",
                "n_test_units",
            ]
        ]
        .drop_duplicates()
        .rename(columns={"dataset_scope": "dataset", "direction_scope": "direction"})
    )
    raw = raw.merge(
        stats_subset,
        on=["dataset", "analysis_family", "direction", "metric_name"],
        how="left",
        validate="many_to_one",
    )

    summary = (
        raw.groupby(
            ["dataset", "analysis_family", "direction", "metric_name", "representation"],
            dropna=False,
            sort=True,
        )
        .agg(n_units=("unit_id", "nunique"))
        .reset_index()
    )
    quantiles = summarize_quantiles_frame(
        raw,
        group_columns=["dataset", "analysis_family", "direction", "metric_name", "representation"],
        value_column="metric_value",
    )
    summary = summary.merge(
        quantiles,
        on=["dataset", "analysis_family", "direction", "metric_name", "representation"],
        how="left",
        validate="one_to_one",
    )
    summary["mean_value"] = summary["summary_mean"]
    summary["ci_low"] = summary["q25_value"]
    summary["ci_high"] = summary["q75_value"]
    summary["metric_value"] = summary["q50_value"]
    summary = summary.merge(
        stats_subset,
        on=["dataset", "analysis_family", "direction", "metric_name"],
        how="left",
        validate="many_to_one",
    )
    summary["row_kind"] = "summary"
    summary["cell_line"] = pd.NA
    summary["target"] = pd.NA
    summary["unit_source"] = "target_level_performance"
    summary["unit_id"] = "summary"

    combined = pd.concat(
        [
            raw[
                [
                    "row_kind",
                    "dataset",
                    "cell_line",
                    "target",
                    "direction",
                    "unit_source",
                    "unit_id",
                    "analysis_family",
                    "metric_name",
                    "representation",
                    "metric_value",
                    "mean_value",
                    "summary_mean",
                    "summary_median",
                    "q10_value",
                    "q25_value",
                    "q50_value",
                    "q75_value",
                    "q90_value",
                    "ci_low",
                    "ci_high",
                    "n_units",
                    "bh_q",
                    "test_status",
                    "effect_direction",
                    "median_delta",
                    "n_test_units",
                ]
            ],
            summary[
                [
                    "row_kind",
                    "dataset",
                    "cell_line",
                    "target",
                    "direction",
                    "unit_source",
                    "unit_id",
                    "analysis_family",
                    "metric_name",
                    "representation",
                    "metric_value",
                    "mean_value",
                    "summary_mean",
                    "summary_median",
                    "q10_value",
                    "q25_value",
                    "q50_value",
                    "q75_value",
                    "q90_value",
                    "ci_low",
                    "ci_high",
                    "n_units",
                    "bh_q",
                    "test_status",
                    "effect_direction",
                    "median_delta",
                    "n_test_units",
                ]
            ],
        ],
        ignore_index=True,
        sort=False,
    ).sort_values(
        ["analysis_family", "metric_name", "dataset", "direction", "representation", "row_kind", "unit_source", "unit_id"],
        kind="mergesort",
    )
    return write_csv(
        combined,
        output_path,
        [
            "row_kind",
            "dataset",
            "cell_line",
            "target",
            "direction",
            "unit_source",
            "unit_id",
            "analysis_family",
            "metric_name",
            "representation",
            "metric_value",
            "mean_value",
            "summary_mean",
            "summary_median",
            "q10_value",
            "q25_value",
            "q50_value",
            "q75_value",
            "q90_value",
            "ci_low",
            "ci_high",
            "n_units",
            "bh_q",
            "test_status",
            "effect_direction",
            "median_delta",
            "n_test_units",
        ],
    )


def build_figure3_panel_3f_tradeoff(
    *,
    target_pattern_path: Path,
    performance_path: Path,
    direction_support_path: Path,
    comparison_stats_path: Path,
    figure3_scope_path: Path,
    task2_retrieval_per_query_path: Path,
    output_path: Path,
) -> Path:
    _ = performance_path
    _ = direction_support_path
    _ = task2_retrieval_per_query_path
    scope_frame = read_csv_required(
        figure3_scope_path,
        {"dataset", "cell_line", "representation", "availability_status"},
        "Figure 3 scope summary",
    )
    valid_representations = set(
        scope_frame.loc[
            scope_frame["dataset"].eq(FM_LOCAL_DATASET)
            & scope_frame["cell_line"].eq(FM_LOCAL_CELL_LINE)
            & scope_frame["availability_status"].eq("available"),
            "representation",
        ].astype(str)
    )
    target_pattern = load_target_pattern(target_pattern_path)
    target_pattern = target_pattern.loc[
        target_pattern["dataset"].eq(FM_LOCAL_DATASET)
        & target_pattern["representation"].isin(valid_representations)
        & target_pattern["analysis_family"].isin(["group_concordance", "retrieval"])
        & (
            target_pattern["analysis_family"].ne("retrieval")
            | target_pattern["direction"].eq("C2G")
        )
        & target_pattern["metric_name"].isin(
            [
                "cosine_centroid",
                "pcc_centroid",
                "mrr_corrected",
                "hit1_corrected",
                "hit5_corrected",
                "hit10_corrected",
            ]
        )
    ].copy()
    summary = (
        target_pattern.groupby(["representation", "analysis_family", "metric_name"], dropna=False, sort=True)
        .agg(
            metric_value=("metric_value", "mean"),
            n_units=("target", "nunique"),
        )
        .reset_index()
    )
    summary_ci = target_pattern.groupby(
        ["representation", "analysis_family", "metric_name"],
        dropna=False,
        sort=True,
    )["metric_value"].apply(simple_ci_bounds)
    summary_ci = summary_ci.reset_index(name="ci_bounds")
    summary_ci["ci_low"] = summary_ci["ci_bounds"].map(lambda value: value[0])
    summary_ci["ci_high"] = summary_ci["ci_bounds"].map(lambda value: value[1])
    summary_ci = summary_ci.drop(columns=["ci_bounds"])
    summary = summary.merge(
        summary_ci,
        on=["representation", "analysis_family", "metric_name"],
        how="left",
        validate="one_to_one",
    )

    comparison = read_csv_required(
        comparison_stats_path,
        {
            "comparison_family",
            "analysis_family",
            "group_a",
            "group_b",
            "bh_q",
            "metric_name",
            "test_status",
            "dataset_scope",
            "direction_scope",
        },
        "Manuscript comparison statistics",
    )
    fm_vs_gene = comparison.loc[
        comparison["comparison_family"].isin(
            [
                "figure3_task2_k562_fm_target_anchored",
                "figure3_task2_k562_fm_c2g_query_anchored",
            ]
        )
        & comparison["group_b"].eq("Gene")
    ].copy()
    fm_vs_gene["representation"] = fm_vs_gene["group_a"].astype(str)
    fm_vs_gene["comparison_source"] = fm_vs_gene["comparison_family"].astype(str)

    pathway_vs_gene = comparison.loc[
        comparison["comparison_family"].eq("figure3_task2_target_pattern_representation")
        & comparison["dataset_scope"].eq(FM_LOCAL_DATASET)
        & (
            (
                comparison["analysis_family"].eq("group_concordance")
                & comparison["direction_scope"].eq("ALL")
            )
            | (
                comparison["analysis_family"].eq("retrieval")
                & comparison["direction_scope"].eq("C2G")
            )
        )
        & (
            (
                comparison["group_a"].eq("Gene")
                & comparison["group_b"].eq("Pathway")
            )
            | (
                comparison["group_a"].eq("Pathway")
                & comparison["group_b"].eq("Gene")
            )
        )
    ].copy()
    pathway_vs_gene["representation"] = np.where(
        pathway_vs_gene["group_a"].eq("Gene"),
        pathway_vs_gene["group_b"],
        pathway_vs_gene["group_a"],
    )
    pathway_vs_gene["comparison_source"] = pathway_vs_gene["comparison_family"].astype(str)

    q_values = pd.concat([fm_vs_gene, pathway_vs_gene], ignore_index=True, sort=False)
    q_values = q_values[
        [
            "analysis_family",
            "metric_name",
            "representation",
            "bh_q",
            "test_status",
            "comparison_source",
        ]
    ].rename(
        columns={
            "bh_q": "bh_q_vs_gene",
            "test_status": "test_status_vs_gene",
        }
    )
    q_values = q_values.drop_duplicates(
        subset=["analysis_family", "metric_name", "representation"],
        keep="first",
    )

    out = summary.merge(
        q_values,
        on=["analysis_family", "metric_name", "representation"],
        how="left",
        validate="one_to_one",
    )
    out["dataset"] = FM_LOCAL_DATASET
    out["cell_line"] = FM_LOCAL_CELL_LINE
    out["metric_panel"] = [
        task2_metric_block(row.analysis_family, "C2G" if row.analysis_family == "retrieval" else "ALL", row.metric_name)
        for row in out.itertuples(index=False)
    ]
    out["is_gene_baseline"] = out["representation"].eq("Gene")
    gene_reference = (
        out.loc[out["is_gene_baseline"], ["analysis_family", "metric_name", "metric_value"]]
        .drop_duplicates()
        .rename(columns={"metric_value": "gene_reference_metric_value"})
    )
    out = out.merge(
        gene_reference,
        on=["analysis_family", "metric_name"],
        how="left",
        validate="many_to_one",
    )
    out["gene_reference_plot_value"] = out["gene_reference_metric_value"]
    out["comparison_source"] = out["comparison_source"].fillna("baseline_or_not_available")
    out["fm_scope_note"] = "scPerturb/K562-local only; not benchmark-wide."
    if not out["dataset"].eq(FM_LOCAL_DATASET).all():
        raise ValueError("Figure 3F staged output must remain scPerturb-only.")
    if not out["cell_line"].eq(FM_LOCAL_CELL_LINE).all():
        raise ValueError("Figure 3F staged output must remain K562-only.")
    out = out.sort_values(
        ["analysis_family", "metric_name", "representation"],
        kind="mergesort",
    )
    return write_csv(
        out,
        output_path,
        [
            "dataset",
            "cell_line",
            "representation",
            "analysis_family",
            "metric_name",
            "metric_panel",
            "metric_value",
            "gene_reference_metric_value",
            "gene_reference_plot_value",
            "ci_low",
            "ci_high",
            "n_units",
            "bh_q_vs_gene",
            "test_status_vs_gene",
            "comparison_source",
            "is_gene_baseline",
            "fm_scope_note",
        ],
    )


def build_figure3_panel_3j(*, source_path: Path, output_path: Path) -> Path:
    long_frame = build_task1_contextual_long(source_path)
    long_frame = filter_excluded_metrics(long_frame)
    common = long_frame.loc[long_frame["group_label"].isin(TASK2_COMMON_BASE_REPRESENTATIONS)].copy()
    grouped = (
        common.groupby(["dataset", "analysis_family", "metric_name", "group_label"], dropna=False, sort=True)
        .agg(
            n_triplets=("target", "size"),
            n_cell_lines=("cell_line", "nunique"),
            n_targets=("target", "nunique"),
            mean_metric_value=("metric_value", "mean"),
            median_metric_value=("metric_value", "median"),
            iqr_low=("metric_value", lambda values: float(pd.Series(values).quantile(0.25))),
            iqr_high=("metric_value", lambda values: float(pd.Series(values).quantile(0.75))),
            min_metric_value=("metric_value", "min"),
            max_metric_value=("metric_value", "max"),
        )
        .reset_index()
        .rename(columns={"group_label": "representation"})
    )
    grouped["reference_only_note"] = "background/reference only; not bridge/ceiling."
    return write_csv(
        grouped,
        output_path,
        [
            "dataset",
            "analysis_family",
            "metric_name",
            "representation",
            "n_triplets",
            "n_cell_lines",
            "n_targets",
            "mean_metric_value",
            "median_metric_value",
            "iqr_low",
            "iqr_high",
            "min_metric_value",
            "max_metric_value",
            "reference_only_note",
        ],
    )


def build_extended_figure8_context(*, source_path: Path, output_path: Path) -> Path:
    long_frame = build_task1_contextual_long(source_path)
    long_frame = filter_excluded_metrics(long_frame)
    long_frame = long_frame.rename(columns={"group_label": "representation"}).copy()
    long_frame = dense_rank_desc(
        long_frame,
        ["dataset", "analysis_family", "metric_name", "representation"],
        "metric_value",
        "rank_within_dataset_metric_representation",
    )
    long_frame["fm_local_bool"] = long_frame["representation"].map(fm_local_bool)
    long_frame["context_support_scope"] = "task1_internal_contextual_support"
    long_frame["context_note"] = (
        "Repeated perturbation_type rows collapsed by mean within the frozen manuscript contextual unit."
    )
    long_frame["reference_only_note"] = "background/reference only; not bridge/ceiling."
    return write_csv(
        long_frame.sort_values(
            [
                "dataset",
                "analysis_family",
                "metric_name",
                "representation",
                "rank_within_dataset_metric_representation",
                "cell_line",
                "target",
            ],
            kind="mergesort",
        ),
        output_path,
        [
            "dataset",
            "cell_line",
            "target",
            "analysis_family",
            "direction",
            "metric_name",
            "representation",
            "metric_value",
            "rank_within_dataset_metric_representation",
            "fm_local_bool",
            "context_support_scope",
            "context_note",
            "reference_only_note",
        ],
    )


def select_best_available_row(
    frame: pd.DataFrame,
    *,
    score_column: str,
    ascending: bool,
    excluded_targets: set[str],
) -> pd.Series | None:
    if frame.empty:
        return None
    available = frame.loc[~frame["target"].astype(str).isin(excluded_targets)].copy()
    if available.empty:
        return None
    ordered = available.sort_values(
        [score_column, "target"],
        ascending=[ascending, True],
        kind="mergesort",
    )
    if ordered.empty:
        return None
    return ordered.iloc[0]


def build_extended_figure9_support_diagnostic(
    *,
    cell_line_source_path: Path,
    target_source_path: Path,
    output_path: Path,
) -> Path:
    cell_surface = build_task2_pattern_surface_wide(
        source_path=cell_line_source_path,
        entity_column="cell_line",
        support_group_column="n_targets",
    ).copy()
    cell_surface["surface_type"] = "cell_line"
    cell_surface["entity_key"] = cell_surface["cell_line"].astype(str)

    target_surface = build_task2_pattern_surface_wide(
        source_path=target_source_path,
        entity_column="target",
        support_group_column="n_cell_lines",
    ).copy()
    target_surface["surface_type"] = "target"
    target_surface["entity_key"] = target_surface["target"].astype(str)

    combined = pd.concat([cell_surface, target_surface], ignore_index=True, sort=False)
    combined["support_log10"] = np.log10(pd.to_numeric(combined["support_value"], errors="coerce").clip(lower=1.0))
    combined["representation_gap"] = pd.to_numeric(
        combined["pathway_enrichment_score"], errors="coerce"
    ) - pd.to_numeric(combined["gene_enrichment_score"], errors="coerce")
    combined["shared_label"] = np.where(combined["shared_flag"], "shared", "dataset_specific")
    combined = dense_rank_desc(
        combined,
        ["surface_type", "dataset"],
        "suitability_score",
        "rank_within_dataset_surface",
    )
    combined["surface_note"] = (
        "Diagnostic-only support: relates suitability summaries to support size; not a causal biology claim."
    )
    return write_csv(
        combined.sort_values(
            ["surface_type", "dataset", "rank_within_dataset_surface", "entity_key"],
            kind="mergesort",
        ),
        output_path,
        [
            "surface_type",
            "dataset",
            "entity_key",
            "entity_display_label",
            "suitability_score",
            "support_value",
            "support_n",
            "support_log10",
            "shared_across_modalities_bool",
            "shared_flag",
            "shared_label",
            "representation_gap",
            "pair_mean_enrichment",
            "max_abs_enrichment",
            "gene_enrichment_score",
            "pathway_enrichment_score",
            "gene_significance_label",
            "pathway_significance_label",
            "rank_within_dataset_surface",
            "surface_note",
        ],
    )


def build_extended_figure9_target_exemplars(
    *,
    source_path: Path,
    output_path: Path,
) -> Path:
    target_surface = build_task2_pattern_surface_wide(
        source_path=source_path,
        entity_column="target",
        support_group_column="n_cell_lines",
    ).copy()
    target_surface["target"] = target_surface["target"].astype(str)
    target_surface["representation_gap"] = pd.to_numeric(
        target_surface["pathway_enrichment_score"], errors="coerce"
    ) - pd.to_numeric(target_surface["gene_enrichment_score"], errors="coerce")

    shared = target_surface.loc[target_surface["shared_flag"]].copy()
    shared_summary = (
        shared.groupby("target", dropna=False, sort=True)
        .agg(
            entity_display_label=("entity_display_label", "first"),
            shared_mean_suitability=("suitability_score", "mean"),
            shared_abs_gap=("suitability_score", lambda values: float(values.max() - values.min())),
            gene_advantage_mean=("representation_gap", lambda values: float((-pd.to_numeric(values, errors="coerce")).mean())),
            pathway_advantage_mean=("representation_gap", lambda values: float(pd.to_numeric(values, errors="coerce").mean())),
            dataset_count=("dataset", "nunique"),
        )
        .reset_index()
    )
    shared_summary = shared_summary.loc[shared_summary["dataset_count"].ge(2)].copy()
    lincs_only = target_surface.loc[target_surface["dataset"].eq("LINCS") & ~target_surface["shared_flag"]].copy()
    scperturb_only = target_surface.loc[
        target_surface["dataset"].eq("scPerturb") & ~target_surface["shared_flag"]
    ].copy()

    category_specs = [
        (
            "shared_consensus_high",
            "Shared high-suitability anchor",
            "Highest mean suitability across shared targets.",
            shared_summary,
            "shared_mean_suitability",
            False,
        ),
        (
            "shared_cross_dataset_gap",
            "Shared cross-dataset gap",
            "Largest suitability difference between LINCS and scPerturb among shared targets.",
            shared_summary,
            "shared_abs_gap",
            False,
        ),
        (
            "shared_gene_favored",
            "Shared Gene-favored target",
            "Largest Gene-over-Pathway enrichment advantage among shared targets.",
            shared_summary,
            "gene_advantage_mean",
            False,
        ),
        (
            "shared_pathway_favored",
            "Shared Pathway-favored target",
            "Largest Pathway-over-Gene enrichment advantage among shared targets.",
            shared_summary,
            "pathway_advantage_mean",
            False,
        ),
        (
            "lincs_dataset_specific_high",
            "LINCS-specific high target",
            "Top dataset-specific target in LINCS by suitability.",
            lincs_only,
            "suitability_score",
            False,
        ),
        (
            "scperturb_dataset_specific_high",
            "scPerturb-specific high target",
            "Top dataset-specific target in scPerturb by suitability.",
            scperturb_only,
            "suitability_score",
            False,
        ),
    ]

    selected_meta: list[dict[str, object]] = []
    excluded_targets: set[str] = set()
    for category_id, category_label, selection_reason, frame, score_column, ascending in category_specs:
        selected = select_best_available_row(
            frame,
            score_column=score_column,
            ascending=ascending,
            excluded_targets=excluded_targets,
        )
        if selected is None:
            continue
        target = str(selected["target"])
        excluded_targets.add(target)
        selected_meta.append(
            {
                "selection_category": category_id,
                "selection_label": category_label,
                "selection_reason": selection_reason,
                "target": target,
            }
        )

    if not selected_meta:
        raise ValueError("Extended Figure 9 exemplar selection produced no targets.")

    meta = pd.DataFrame(selected_meta)
    meta["selection_rank"] = np.arange(1, len(meta) + 1, dtype=int)
    selected_rows = target_surface.merge(meta, on="target", how="inner", validate="many_to_one")
    selected_rows["row_display_label"] = (
        selected_rows["selection_rank"].astype(str)
        + ". "
        + selected_rows["selection_label"].astype(str)
        + " | "
        + selected_rows["entity_display_label"].astype(str)
    )

    long_rows: list[dict[str, object]] = []
    for row in selected_rows.to_dict("records"):
        for representation in ("Gene", "Pathway"):
            enrichment_score = row[f"{representation.lower()}_enrichment_score"]
            significance_label = row[f"{representation.lower()}_significance_label"]
            long_rows.append(
                {
                    "selection_rank": row["selection_rank"],
                    "selection_category": row["selection_category"],
                    "selection_label": row["selection_label"],
                    "selection_reason": row["selection_reason"],
                    "target": row["target"],
                    "dataset": row["dataset"],
                    "entity_display_label": row["entity_display_label"],
                    "row_display_label": row["row_display_label"],
                    "representation": representation,
                    "enrichment_score": enrichment_score,
                    "significance_label": significance_label,
                    "support_value": row["support_value"],
                    "support_n": row["support_n"],
                    "suitability_score": row["suitability_score"],
                    "pair_mean_enrichment": row["pair_mean_enrichment"],
                    "max_abs_enrichment": row["max_abs_enrichment"],
                    "representation_gap": row["representation_gap"],
                    "shared_flag": row["shared_flag"],
                }
            )
    out = pd.DataFrame(long_rows)
    out["shared_label"] = np.where(out["shared_flag"], "shared", "dataset_specific")
    out["exemplar_note"] = (
        "Illustrative target-level archetypes chosen from the frozen pattern surface; interpret as benchmark-supporting exemplars, not causal mechanism proof."
    )
    return write_csv(
        out.sort_values(["selection_rank", "dataset", "representation"], kind="mergesort"),
        output_path,
        [
            "selection_rank",
            "selection_category",
            "selection_label",
            "selection_reason",
            "target",
            "dataset",
            "entity_display_label",
            "row_display_label",
            "representation",
            "enrichment_score",
            "significance_label",
            "support_value",
            "support_n",
            "suitability_score",
            "pair_mean_enrichment",
            "max_abs_enrichment",
            "representation_gap",
            "shared_flag",
            "shared_label",
            "exemplar_note",
        ],
    )


def build_one(name: str, args: argparse.Namespace) -> Path:
    output_path = args.plot_ready_root / PLOT_READY_RELATIVE_OUTPUTS[name]
    if name == "figure1_panel_1b_dataset_context_coverage.csv":
        return build_figure1_coverage(
            task1_inventory_path=args.task1_inventory_path,
            task2_pairs_coverage_path=args.task2_pairs_coverage_path,
            output_path=output_path,
        )
    if name == "figure1_panel_1c_lawful_scope_matrix.csv":
        return build_figure1_scope_matrix(output_path=output_path)
    if name == "figure1_panel_1d_representation_modifier_availability.csv":
        return build_figure1_rep_modifier_matrix(
            representation_registry_path=args.representation_registry_path,
            figure3_scope_path=args.figure3_scope_path,
            output_path=output_path,
        )
    if name == "figure1_panel_1e_result_object_map.csv":
        return build_figure1_object_map(manifest_path=args.framework_manifest_path, output_path=output_path)
    if name == "extended_figure1_support_registry.csv":
        figure1_dir = args.plot_ready_root / "figure1"
        coverage_path = figure1_dir / "figure1_panel_1b_dataset_context_coverage.csv"
        scope_path = figure1_dir / "figure1_panel_1c_lawful_scope_matrix.csv"
        rep_modifier_path = figure1_dir / "figure1_panel_1d_representation_modifier_availability.csv"
        object_map_path = figure1_dir / "figure1_panel_1e_result_object_map.csv"
        for dep_name, dep_path in [
            ("coverage", coverage_path),
            ("scope matrix", scope_path),
            ("representation/modifier matrix", rep_modifier_path),
            ("object map", object_map_path),
        ]:
            require_path(dep_path, f"Extended Figure 1 dependency ({dep_name})")
        return build_extended_figure1_registry(
            coverage_path=coverage_path,
            scope_path=scope_path,
            rep_modifier_path=rep_modifier_path,
            object_map_path=object_map_path,
            output_path=output_path,
        )
    if name == "figure2_panel_2a_task1_scope.csv":
        return build_figure2_panel_2a(
            source_path=args.analysis_root / "figure2_task1_scope_summary.csv",
            output_path=output_path,
        )
    if name == "figure2_panel_2b_internal_performance_overview.csv":
        return build_task1_performance_overview(
            source_path=args.analysis_root / "figure2_task1_performance_structure.csv",
            scope_value="internal",
            output_path=output_path,
        )
    if name == "figure2_panel_2c_gene_vs_pathway_matched_units.csv":
        return build_figure2_panel_2c(
            bridge_detail_path=args.task1_bridge_detail_path,
            comparison_stats_path=args.comparison_stats_path,
            output_path=output_path,
        )
    if name == "figure2_panel_2d_internal_to_cross_degradation.csv":
        return build_figure2_panel_2d(
            source_path=args.task1_bridge_detail_path,
            comparison_stats_path=args.comparison_stats_path,
            output_path=output_path,
        )
    if name == "figure2_panel_2e_cell_line_high_concordance_summary.csv":
        return build_figure2_panel_2e(
            source_path=args.figure2_cell_high_concordance_path,
            output_path=output_path,
        )
    if name == "figure2_panel_2f_target_high_concordance_summary.csv":
        return build_figure2_panel_2f(
            source_path=args.figure2_target_high_concordance_path,
            output_path=output_path,
        )
    if name == "extended_figure4_cell_line_high_concordance_full.csv":
        return build_extended_figure4_cell_line(
            source_path=args.figure2_cell_high_concordance_path,
            output_path=output_path,
        )
    if name == "extended_figure4_target_high_concordance_full.csv":
        return build_extended_figure4_target(
            source_path=args.figure2_target_high_concordance_path,
            output_path=output_path,
        )
    if name == "figure3_panel_3b_c2g_performance_overview.csv":
        return build_figure3_panel_3b(
            performance_path=args.figure3_performance_path,
            direction_support_path=args.figure3_direction_support_path,
            task2_retrieval_per_query_path=args.task2_retrieval_per_query_path,
            drug_meta_path=args.task2_drug_meta_path,
            output_path=output_path,
        )
    if name == "figure3_panel_3c_cell_line_pattern.csv":
        return build_figure3_panel_3c(
            source_path=args.analysis_root / "figure3_task2_cell_line_pattern_summary.csv",
            output_path=output_path,
        )
    if name == "figure3_panel_3d_target_pattern_summary.csv":
        return build_figure3_panel_3d(
            source_path=args.figure3_target_pattern_path,
            output_path=output_path,
        )
    if name == "figure3_panel_3e_gene_vs_pathway_paired.csv":
        return build_figure3_panel_3e(
            task2_group_root=args.task2_group_root,
            task2_retrieval_per_query_path=args.task2_retrieval_per_query_path,
            comparison_stats_path=args.comparison_stats_path,
            output_path=output_path,
        )
    if name == "extended_figure6_target_pattern_full.csv":
        return build_extended_figure6_target(
            source_path=args.figure3_target_pattern_path,
            output_path=output_path,
        )
    if name == "figure3_panel_3f_fm_local_tradeoff.csv":
        return build_figure3_panel_3f_tradeoff(
            target_pattern_path=args.figure3_target_pattern_path,
            performance_path=args.figure3_performance_path,
            direction_support_path=args.figure3_direction_support_path,
            comparison_stats_path=args.comparison_stats_path,
            figure3_scope_path=args.figure3_scope_path,
            task2_retrieval_per_query_path=args.task2_retrieval_per_query_path,
            output_path=output_path,
        )
    if name == "figure3_panel_3j_task1_contextual_support_reference.csv":
        return build_figure3_panel_3j(
            source_path=args.figure3_contextual_support_path,
            output_path=output_path,
        )
    if name == "extended_figure8_contextual_support_full.csv":
        return build_extended_figure8_context(
            source_path=args.figure3_contextual_support_path,
            output_path=output_path,
        )
    if name == "extended_figure9_support_vs_suitability.csv":
        return build_extended_figure9_support_diagnostic(
            cell_line_source_path=args.analysis_root / "figure3_task2_cell_line_pattern_summary.csv",
            target_source_path=args.figure3_target_pattern_path,
            output_path=output_path,
        )
    if name == "extended_figure9_target_exemplars.csv":
        return build_extended_figure9_target_exemplars(
            source_path=args.figure3_target_pattern_path,
            output_path=output_path,
        )
    raise ValueError(f"Unhandled build target: {name}")


def main() -> None:
    args = parse_args()
    targets = (
        list(PLOT_READY_RELATIVE_OUTPUTS.keys())
        if args.build_target == "all"
        else (
            SPATIAL_REVISION_PLOT_READY_OUTPUTS
            if args.build_target == "spatial_revision"
            else (P0_PLOT_READY_OUTPUTS if args.build_target == "p0" else [args.build_target])
        )
    )

    # Figure 1 registry rows are upstream dependencies for the EF1 concatenation.
    if args.build_target == "all":
        targets = [
            "figure1_panel_1b_dataset_context_coverage.csv",
            "figure1_panel_1c_lawful_scope_matrix.csv",
            "figure1_panel_1d_representation_modifier_availability.csv",
            "figure1_panel_1e_result_object_map.csv",
            "extended_figure1_support_registry.csv",
            "figure2_panel_2a_task1_scope.csv",
            "figure2_panel_2b_internal_performance_overview.csv",
            "figure2_panel_2c_gene_vs_pathway_matched_units.csv",
            "figure2_panel_2d_internal_to_cross_degradation.csv",
            "figure2_panel_2e_cell_line_high_concordance_summary.csv",
            "figure2_panel_2f_target_high_concordance_summary.csv",
            "figure3_panel_3b_c2g_performance_overview.csv",
            "figure3_panel_3c_cell_line_pattern.csv",
            "extended_figure4_cell_line_high_concordance_full.csv",
            "extended_figure4_target_high_concordance_full.csv",
            "figure3_panel_3d_target_pattern_summary.csv",
            "figure3_panel_3e_gene_vs_pathway_paired.csv",
            "extended_figure6_target_pattern_full.csv",
            "figure3_panel_3f_fm_local_tradeoff.csv",
            "figure3_panel_3j_task1_contextual_support_reference.csv",
            "extended_figure8_contextual_support_full.csv",
            "extended_figure9_support_vs_suitability.csv",
            "extended_figure9_target_exemplars.csv",
        ]
    elif args.build_target == "p0":
        targets = list(P0_PLOT_READY_OUTPUTS)
    elif args.build_target == "spatial_revision":
        targets = list(SPATIAL_REVISION_PLOT_READY_OUTPUTS)

    for target in targets:
        output_path = build_one(target, args)
        frame = pd.read_csv(output_path)
        print(f"{output_path}: {len(frame)} rows")


if __name__ == "__main__":
    main()
