#!/usr/bin/env python3
"""
Assemble the frozen downstream manuscript analysis objects for manuscript support.

Status:
- active manuscript canonical builder

Manuscript role:
- materializes the frozen 10-object downstream backbone and registry
- does not create a dedicated Figure 1 object

Architecture:
- see `scripts/ARCHITECTURE.md` for script-family classification

Canonical downstream contract for the current frozen phase:
- Figure 2
  - `figure2_task1_scope_summary.csv`
  - `figure2_task1_performance_structure.csv`
  - `figure2_task1_internal_to_cross_degradation_summary.csv`
- Figure 3
  - `figure3_task2_scope_summary.csv`
  - `figure3_task2_direction_support_summary.csv`
  - `figure3_task2_performance_structure.csv`
  - `figure3_task2_cell_line_pattern_summary.csv`
  - `figure3_task2_target_pattern_summary.csv`
  - `figure3_task1_internal_contextual_support_summary.csv`

Explicit exclusions from the canonical contract:
- No Task1 <-> Task2 pairwise bridge object is canonical in this phase.
- Figure 3 uses Task1 internal contextual support only.
- Figure 3 Task1 support must therefore not use Task1 cross.

FM scope boundaries:
- scPerturb-only slices keep FM alongside Gene and Pathway.
- Any LINCS-involving slice excludes FM.
- Task1 cross excludes FM because the frozen cross scope is Gene/Pathway only.
- Figure 3 target-anchored summaries use canonical Task2 target metadata only;
  retrieval target summaries use C2G single-target queries in current manuscript
  mode (queries where query_target_token is non-null and non-compound).

Execution surface:
- Inputs resolve from explicit upstream roots.
- Canonical manuscript analysis outputs materialize under
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`.
- Active support surfaces live under
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support`.
- Historical/non-canonical artifacts live under
  `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history`.
- The authoritative canonical registry is `analysis/framework_analysis_manifest.json`.
- No canonical source or output path depends on a repo-root `runs/` directory.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
try:
    import pyarrow.parquet as pq
except ModuleNotFoundError:
    pq = None

try:
    from path_policy import (
        DEFAULT_ARCHIVE_ROOT,
        DEFAULT_MANUSCRIPT_ACTIVE_ROOT,
        DEFAULT_MANUSCRIPT_HISTORY_ROOT,
        DEFAULT_MANUSCRIPT_PHASE1_ROOT,
        DEFAULT_MANUSCRIPT_SUPPORT_ROOT,
        DEFAULT_RUNS_ROOT,
        DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT,
        resolve_path,
    )
except ModuleNotFoundError:
    from scripts.path_policy import (
        DEFAULT_ARCHIVE_ROOT,
        DEFAULT_MANUSCRIPT_ACTIVE_ROOT,
        DEFAULT_MANUSCRIPT_HISTORY_ROOT,
        DEFAULT_MANUSCRIPT_PHASE1_ROOT,
        DEFAULT_MANUSCRIPT_SUPPORT_ROOT,
        DEFAULT_RUNS_ROOT,
        DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT,
        resolve_path,
    )

try:
    from s1_task1_internal_metrics import (
        EXPECTED_TASK1_SNAPSHOT,
        GLOBAL_SEED,
        build_cohorts,
        cosine_similarity,
        deterministic_subsample_by_uid,
        energy_distance_biascorr,
        init_global_rng,
        make_group_id,
        pearson_corr,
    )
except ModuleNotFoundError:
    from scripts.s1_task1_internal_metrics import (
        EXPECTED_TASK1_SNAPSHOT,
        GLOBAL_SEED,
        build_cohorts,
        cosine_similarity,
        deterministic_subsample_by_uid,
        energy_distance_biascorr,
        init_global_rng,
        make_group_id,
        pearson_corr,
    )


ROOT = Path(__file__).resolve().parents[1]
NAS_RUNS_ROOT = DEFAULT_RUNS_ROOT
NAS_ARCHIVE_ROOT = DEFAULT_ARCHIVE_ROOT

DEFAULT_TASK1_INTERNAL_ROOT = NAS_RUNS_ROOT / "s1_task1_internal_metrics_0303" / "s1_task1_internal_metrics"
DEFAULT_TASK1_CROSS_ROOT = NAS_RUNS_ROOT / "s2_task1_cross_metrics_0303" / "s2_task1_cross_metrics"
DEFAULT_TASK2_GROUP_ROOT = (
    NAS_RUNS_ROOT / "s4_multisource_impl_verify_20260310_c" / "s4_task2_group_concordance_multisource"
)
DEFAULT_TASK2_RETRIEVAL_ROOT = (
    NAS_RUNS_ROOT / "s5_multisource_impl_verify_20260311_a" / "s5_task2_retrieval_multisource"
)
DEFAULT_TASK2_SYNTHESIS_ROOT = (
    NAS_RUNS_ROOT / "s6_multisource_impl_verify_20260311_a" / "s6_task2_result_synthesis_multisource"
)
DEFAULT_MANUSCRIPT_ACTIVE_ROOT = DEFAULT_MANUSCRIPT_ACTIVE_ROOT
DEFAULT_MANUSCRIPT_SUPPORT_ROOT = DEFAULT_MANUSCRIPT_SUPPORT_ROOT
DEFAULT_MANUSCRIPT_HISTORY_ROOT = DEFAULT_MANUSCRIPT_HISTORY_ROOT
DEFAULT_MANUSCRIPT_PHASE1_ROOT = DEFAULT_MANUSCRIPT_PHASE1_ROOT
DEFAULT_TASK1_SNAPSHOT_ROOT = resolve_path(ROOT, EXPECTED_TASK1_SNAPSHOT)
DEFAULT_TASK2_SNAPSHOT_ROOT = DEFAULT_TASK2_CORRECTED_SNAPSHOT_ROOT

ANALYSIS_OUTPUT_FILENAMES = {
    "figure2_task1_scope_summary": "figure2_task1_scope_summary.csv",
    "figure2_task1_performance_structure": "figure2_task1_performance_structure.csv",
    "figure2_task1_internal_to_cross_degradation_summary": "figure2_task1_internal_to_cross_degradation_summary.csv",
    "figure3_task2_scope_summary": "figure3_task2_scope_summary.csv",
    "figure3_task2_direction_support_summary": "figure3_task2_direction_support_summary.csv",
    "figure3_task2_performance_structure": "figure3_task2_performance_structure.csv",
    "figure3_task2_cell_line_pattern_summary": "figure3_task2_cell_line_pattern_summary.csv",
    "figure3_task2_target_pattern_summary": "figure3_task2_target_pattern_summary.csv",
    "figure3_task1_internal_contextual_support_summary": "figure3_task1_internal_contextual_support_summary.csv",
    "manuscript_comparison_statistics": "manuscript_comparison_statistics.csv",
    "framework_analysis_manifest": "framework_analysis_manifest.json",
}

HISTORICAL_NONCANONICAL_CANDIDATES = [
    {
        "object_name": "Figure 2 A1 metric overview",
        "relative_path": Path("analysis") / "figure2_a1_metric_overview.csv",
        "reason_no_longer_canonical": "Figure 2 now uses the explicit scope, performance, and degradation summaries instead of legacy A1 overview tables.",
    },
    {
        "object_name": "Figure 2 A1 dataset-cellline overview",
        "relative_path": Path("analysis") / "figure2_a1_dataset_cellline_overview.csv",
        "reason_no_longer_canonical": "Figure 2 now uses the explicit scope, performance, and degradation summaries instead of legacy A1 overview tables.",
    },
    {
        "object_name": "Task1-Task2 pairwise bridge table",
        "relative_path": Path("group_bridge") / "task1_task2_group_bridge.csv",
        "reason_no_longer_canonical": "Frozen downstream contract cancels canonical Task1 <-> Task2 pairwise bridge analysis.",
    },
    {
        "object_name": "Task1-Task2 pairwise bridge summary",
        "relative_path": Path("group_bridge") / "task1_task2_group_bridge_summary.csv",
        "reason_no_longer_canonical": "Frozen downstream contract cancels canonical Task1 <-> Task2 pairwise bridge analysis.",
    },
    {
        "object_name": "Figure 3 A2 bridge support matrix",
        "relative_path": Path("analysis") / "figure3_a2_bridge_support_matrix.csv",
        "reason_no_longer_canonical": "Figure 3 now uses Task1 internal contextual support only, not A2 bridge logic.",
    },
    {
        "object_name": "Figure 3 A2 metricpair overview",
        "relative_path": Path("analysis") / "figure3_a2_metricpair_overview.csv",
        "reason_no_longer_canonical": "Pairwise Task1/Task2 bridge overviews are outside the frozen current phase.",
    },
    {
        "object_name": "Figure 3 A2 dataset-cellline overview",
        "relative_path": Path("analysis") / "figure3_a2_dataset_cellline_overview.csv",
        "reason_no_longer_canonical": "Pairwise Task1/Task2 bridge overviews are outside the frozen current phase.",
    },
    {
        "object_name": "Figure 3 Task2 C2G covariate summary",
        "relative_path": Path("analysis") / "figure3_task2_c2g_covariate_summary.csv",
        "reason_no_longer_canonical": "All C2G-side manuscript experiment objects were removed from the canonical current-phase contract.",
    },
    {
        "object_name": "Figure 3 B6 dataset-cellline overview",
        "relative_path": Path("analysis") / "figure3_b6_dataset_cellline_overview.csv",
        "reason_no_longer_canonical": "C2G-oriented Figure 3 B6 overview outputs are historical only under the frozen scope.",
    },
    {
        "object_name": "Figure 3 B6 metric-dose-time overview",
        "relative_path": Path("analysis") / "figure3_b6_metric_dose_time_overview.csv",
        "reason_no_longer_canonical": "C2G-oriented Figure 3 B6 overview outputs are historical only under the frozen scope.",
    },
    {
        "object_name": "Figure 3 B6 supported slices",
        "relative_path": Path("analysis") / "figure3_b6_supported_slices.csv",
        "reason_no_longer_canonical": "C2G-oriented Figure 3 B6 overview outputs are historical only under the frozen scope.",
    },
    {
        "object_name": "Task2 C2G dose/time sensitivity detail",
        "relative_path": Path("sensitivity") / "task2_c2g_dose_time_sensitivity.csv",
        "reason_no_longer_canonical": "C2G-side sensitivity artifacts may remain on disk as historical support but are non-canonical in this phase.",
    },
    {
        "object_name": "Task2 C2G dose/time sensitivity summary",
        "relative_path": Path("sensitivity") / "task2_c2g_dose_time_sensitivity_summary.csv",
        "reason_no_longer_canonical": "C2G-side sensitivity artifacts may remain on disk as historical support but are non-canonical in this phase.",
    },
    {
        "object_name": "Figure 3 Task2 merged descriptive support summary",
        "relative_path": Path("analysis") / "figure3_task2_preference_support_summary.csv",
        "reason_no_longer_canonical": "Replaced by explicit Task2 cell_line and target pattern summaries.",
    },
    {
        "object_name": "Legacy manuscript-active manifest",
        "relative_path": Path("manifest.json"),
        "reason_no_longer_canonical": "Superseded by the framework analysis manifest, which declares the frozen canonical and historical object sets explicitly.",
    },
    {
        "object_name": "Legacy figure analysis update note",
        "relative_path": Path("analysis") / "figure_analysis_update.md",
        "reason_no_longer_canonical": "Summarizes legacy figure-input backfills that are no longer canonical current-phase outputs.",
    },
    {
        "object_name": "Legacy manuscript-active experiment update note",
        "relative_path": Path("notes") / "manuscript_active_experiment_update.md",
        "reason_no_longer_canonical": "References historical A2 bridge and C2G sensitivity surfaces that are non-canonical in the frozen current phase.",
    },
]

MANIFEST_CONTRACT_NOTES = [
    "Figure 3 Task1 support is internal-only contextual support.",
    "Figure 3 Task1 support includes both Task1 internal group and Task1 internal retrieval surfaces.",
    "Figure 3 Task2 direction support and Task2 performance structure are separate canonical objects.",
    "Task1/Task2 pairwise bridge outputs are historical and non-canonical for this phase.",
    "C2G-side B6 sensitivity artifacts are historical; C2G single-target retrieval is now canonical for pattern tables.",
    "Legacy manuscript-active manifest and figure-input notes are historical only; the framework analysis manifest is the authoritative current-phase registry.",
    "scPerturb-only slices include FM; any LINCS-involving slice excludes FM.",
    "Task2 target-anchored summaries use canonical target metadata only and therefore use C2G single-target retrieval in current manuscript mode.",
    "No canonical output depends on a repo-root runs/ directory.",
]

TASK1_GROUP_METRICS = {"mean_cosine_centroid", "mean_pcc_centroid", "mean_edist_biascorr"}
TASK1_RETRIEVAL_METRICS = ("mrr_corrected", "hit1_corrected", "hit5_corrected", "hit10_corrected")
TASK1_RETRIEVAL_BASES = ("mrr", "hit1", "hit5", "hit10")
TASK1_CONTEXT_GROUP_METRICS = ("cosine", "pcc", "edist")
TASK1_CONTEXT_GROUP_METRIC_NAME_MAP = {
    "cosine": "cosine_centroid",
    "pcc": "pcc_centroid",
    "edist": "edist_biascorr",
}
TASK1_CONTEXT_GROUP_OUTPUT_COLUMNS = [
    "dataset",
    "cell_line",
    "target",
    "perturbation_type",
    "representation",
    "analysis_family",
    "metric_name",
    "metric_value",
    "group_id",
    "n_queries",
    "n_total",
    "n_A",
    "n_B",
    "n_A_sub",
    "n_B_sub",
    "N_gallery_mean",
    "m_pos_mean",
    "context_support_scope",
    "context_note",
    "fm_scope_note",
]

TASK2_GROUP_VALUE_COLUMNS = ("cosine_centroid", "pcc_centroid", "edist_biascorr")
TASK2_RETRIEVAL_VALUE_COLUMNS = ("mrr_corrected", "hit1_corrected", "hit5_corrected", "hit10_corrected")
TASK2_BENCHMARK_ALLOWED_ANALYSIS_FAMILIES = {"group_concordance", "retrieval"}

TASK2_TARGET_SUPPORT_GROUP_COLUMNS = [
    "dataset",
    "cell_line",
    "target",
    "representation",
    "analysis_family",
    "direction",
    "metric_name",
]
TASK1_CONTEXT_GROUP_COLUMNS = [
    "dataset",
    "cell_line",
    "target",
    "perturbation_type",
    "representation",
]
TASK1_CONTEXT_RETRIEVAL_CONTEXT_COLUMNS = ("N_gallery", "m_pos")

GENE_PATHWAY_REPRESENTATIONS = {"Gene", "Pathway"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build active manuscript framework analysis objects.")
    parser.add_argument("--task1-snapshot-root", type=Path, default=DEFAULT_TASK1_SNAPSHOT_ROOT)
    parser.add_argument("--task1-internal-root", type=Path, default=DEFAULT_TASK1_INTERNAL_ROOT)
    parser.add_argument("--task1-cross-root", type=Path, default=DEFAULT_TASK1_CROSS_ROOT)
    parser.add_argument("--task2-group-root", type=Path, default=DEFAULT_TASK2_GROUP_ROOT)
    parser.add_argument("--task2-retrieval-root", type=Path, default=DEFAULT_TASK2_RETRIEVAL_ROOT)
    parser.add_argument("--task2-synthesis-root", type=Path, default=DEFAULT_TASK2_SYNTHESIS_ROOT)
    parser.add_argument("--manuscript-active-root", type=Path, default=DEFAULT_MANUSCRIPT_ACTIVE_ROOT)
    parser.add_argument("--manuscript-support-root", type=Path, default=DEFAULT_MANUSCRIPT_SUPPORT_ROOT)
    parser.add_argument("--manuscript-history-root", type=Path, default=DEFAULT_MANUSCRIPT_HISTORY_ROOT)
    parser.add_argument("--manuscript-phase1-root", type=Path, default=DEFAULT_MANUSCRIPT_PHASE1_ROOT)
    parser.add_argument("--task2-snapshot-root", type=Path, default=DEFAULT_TASK2_SNAPSHOT_ROOT)
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


def join_unique(values: pd.Series) -> str:
    cleaned = sorted({str(value) for value in values.dropna().tolist() if str(value) != ""})
    return "|".join(cleaned)


def nullable_str(value: object) -> object:
    if pd.isna(value):
        return pd.NA
    text = str(value)
    return text if text else pd.NA


def sum_min_count(series: pd.Series) -> object:
    return series.sum(min_count=1)


def is_fm_representation(representation: object) -> bool:
    return str(representation) not in GENE_PATHWAY_REPRESENTATIONS


def task2_fm_scope_note(dataset: object, representation: object) -> str:
    if str(dataset) == "LINCS":
        return "FM excluded because the Task2 slice involves LINCS."
    if is_fm_representation(representation):
        return "FM included because the Task2 slice is scPerturb-only."
    return "Gene/Pathway base representation retained under frozen Task2 scope."


def task1_fm_scope_note(dataset: object, representation: object) -> str:
    if str(dataset) == "LINCS":
        return "FM excluded because LINCS Task1 internal outputs do not carry FM."
    if is_fm_representation(representation):
        return "FM included because the Task1 slice is scPerturb internal."
    return "Gene/Pathway base representation retained under frozen Task1 scope."


def make_grouped_mean_summary(
    frame: pd.DataFrame,
    key_columns: list[str],
    value_columns: tuple[str, ...] | list[str],
    mean_context_columns: tuple[str, ...] | list[str] = (),
    count_column: str = "n_queries",
) -> pd.DataFrame:
    if frame.empty:
        empty_cols = list(key_columns) + [count_column]
        for column in [*value_columns, *mean_context_columns]:
            empty_cols.extend([f"{column}_sum", f"{column}_count"])
        return pd.DataFrame(columns=empty_cols)

    grouped = frame.groupby(key_columns, dropna=False, observed=True, sort=True)
    summary = grouped.size().rename(count_column).to_frame()
    for column in [*value_columns, *mean_context_columns]:
        summary[f"{column}_sum"] = grouped[column].sum(min_count=1)
        summary[f"{column}_count"] = grouped[column].count()
    return summary.reset_index()


def update_mean_state(
    aggregate: dict[tuple[object, ...], dict[str, float | int]],
    summary_frame: pd.DataFrame,
    key_columns: list[str],
    value_columns: tuple[str, ...] | list[str],
    mean_context_columns: tuple[str, ...] | list[str] = (),
    count_column: str = "n_queries",
) -> None:
    tracked_columns = [*value_columns, *mean_context_columns]
    for row in summary_frame.itertuples(index=False):
        key = tuple(getattr(row, column) for column in key_columns)
        state = aggregate.setdefault(key, {count_column: 0})
        state[count_column] += int(getattr(row, count_column))
        for column in tracked_columns:
            state.setdefault(f"{column}_sum", 0.0)
            state.setdefault(f"{column}_count", 0)
            value_sum = getattr(row, f"{column}_sum")
            value_count = getattr(row, f"{column}_count")
            if pd.notna(value_sum):
                state[f"{column}_sum"] += float(value_sum)
            state[f"{column}_count"] += int(value_count)


def finalize_mean_state(
    aggregate: dict[tuple[object, ...], dict[str, float | int]],
    key_columns: list[str],
    value_columns: tuple[str, ...] | list[str],
    mean_context_columns: tuple[str, ...] | list[str] = (),
    count_column: str = "n_queries",
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    tracked_columns = [*value_columns, *mean_context_columns]
    for key, state in aggregate.items():
        row = dict(zip(key_columns, key, strict=True))
        row[count_column] = int(state[count_column])
        for column in tracked_columns:
            total = float(state[f"{column}_sum"])
            observed = int(state[f"{column}_count"])
            row[column] = total / observed if observed else pd.NA
            row[f"{column}_n_valid"] = observed
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_task2_c2g_target_means(retrieval_per_query_path: Path) -> pd.DataFrame:
    require_path(retrieval_per_query_path, "Task2 retrieval per-query parquet")
    if pq is None:
        raise ModuleNotFoundError("pyarrow is required to read Task2 retrieval parquet inputs.")
    parquet_file = pq.ParquetFile(retrieval_per_query_path)
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
    aggregate: dict[tuple[object, ...], dict[str, float | int]] = {}
    key_columns = ["dataset", "cell_line", "target", "representation"]

    for batch in parquet_file.iter_batches(batch_size=200_000, columns=columns):
        frame = batch.to_pandas()
        frame = frame.loc[frame["direction"].eq("C2G")].copy()
        if frame.empty:
            continue

        frame = frame.loc[frame["query_n_targets"].eq(1)].copy()
        if frame.empty:
            continue
        frame["target"] = frame["query_target_tokens"].astype(str).str.split(";").str[0]
        frame = frame.loc[frame["target"].ne("") & frame["target"].ne("NA")].copy()
        if frame.empty:
            continue

        summary = make_grouped_mean_summary(
            frame=frame,
            key_columns=key_columns,
            value_columns=TASK2_RETRIEVAL_VALUE_COLUMNS,
            mean_context_columns=("N_gallery", "m_pos"),
            count_column="n_queries",
        )
        update_mean_state(
            aggregate=aggregate,
            summary_frame=summary,
            key_columns=key_columns,
            value_columns=TASK2_RETRIEVAL_VALUE_COLUMNS,
            mean_context_columns=("N_gallery", "m_pos"),
            count_column="n_queries",
        )

    out = finalize_mean_state(
        aggregate=aggregate,
        key_columns=key_columns,
        value_columns=TASK2_RETRIEVAL_VALUE_COLUMNS,
        mean_context_columns=("N_gallery", "m_pos"),
        count_column="n_queries",
    )
    if out.empty:
        return out
    return out.sort_values(key_columns, kind="mergesort").reset_index(drop=True)


def load_eligible_task2_triplets(task2_snapshot_root: Path) -> pd.DataFrame:
    coverage = read_csv_required(
        task2_snapshot_root / "task2_pairs_coverage.csv",
        {"dataset", "cell_line", "target_token", "is_eligible_bool"},
        "Task2 pairs coverage",
    )
    eligible = coverage.loc[coverage["is_eligible_bool"]].copy()
    eligible["target"] = eligible["target_token"].astype(str)
    eligible = eligible.loc[eligible["target"].ne("") & eligible["target"].ne("NA")].copy()
    if eligible.empty:
        return pd.DataFrame(columns=["dataset", "cell_line", "target"])
    return (
        eligible[["dataset", "cell_line", "target"]]
        .drop_duplicates()
        .sort_values(["dataset", "cell_line", "target"], kind="mergesort")
        .reset_index(drop=True)
    )


def aggregate_task1_internal_context_retrieval_means(
    task1_retrieval_per_query_path: Path,
    eligible_task2_triplets: pd.DataFrame,
) -> pd.DataFrame:
    require_path(task1_retrieval_per_query_path, "Task1 internal retrieval per-query parquet")
    if pq is None:
        raise ModuleNotFoundError("pyarrow is required to read Task1 retrieval parquet inputs.")
    parquet_file = pq.ParquetFile(task1_retrieval_per_query_path)
    columns = [
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "cell_line",
        "target_token",
        "N_gallery",
        "m_pos",
        *TASK1_RETRIEVAL_METRICS,
    ]
    aggregate: dict[tuple[object, ...], dict[str, float | int]] = {}

    for batch in parquet_file.iter_batches(batch_size=200_000, columns=columns):
        frame = batch.to_pandas()
        frame["dataset"] = frame["dataset_or_direction"].astype(str)
        frame["target"] = frame["target_token"].astype(str)
        frame = frame.loc[frame["target"].ne("") & frame["target"].ne("NA")].copy()
        if frame.empty:
            continue

        frame = frame.merge(
            eligible_task2_triplets,
            on=["dataset", "cell_line", "target"],
            how="inner",
            validate="many_to_one",
        )
        if frame.empty:
            continue

        summary = make_grouped_mean_summary(
            frame=frame,
            key_columns=TASK1_CONTEXT_GROUP_COLUMNS,
            value_columns=TASK1_RETRIEVAL_METRICS,
            mean_context_columns=TASK1_CONTEXT_RETRIEVAL_CONTEXT_COLUMNS,
            count_column="n_queries",
        )
        update_mean_state(
            aggregate=aggregate,
            summary_frame=summary,
            key_columns=TASK1_CONTEXT_GROUP_COLUMNS,
            value_columns=TASK1_RETRIEVAL_METRICS,
            mean_context_columns=TASK1_CONTEXT_RETRIEVAL_CONTEXT_COLUMNS,
            count_column="n_queries",
        )

    out = finalize_mean_state(
        aggregate=aggregate,
        key_columns=TASK1_CONTEXT_GROUP_COLUMNS,
        value_columns=TASK1_RETRIEVAL_METRICS,
        mean_context_columns=TASK1_CONTEXT_RETRIEVAL_CONTEXT_COLUMNS,
        count_column="n_queries",
    )
    if out.empty:
        return pd.DataFrame(
            columns=[
                *TASK1_CONTEXT_GROUP_COLUMNS,
                "n_queries",
                *TASK1_RETRIEVAL_METRICS,
                *(f"{metric}_n_valid" for metric in TASK1_RETRIEVAL_METRICS),
                *TASK1_CONTEXT_RETRIEVAL_CONTEXT_COLUMNS,
                *(f"{column}_n_valid" for column in TASK1_CONTEXT_RETRIEVAL_CONTEXT_COLUMNS),
            ]
        )
    return out.sort_values(TASK1_CONTEXT_GROUP_COLUMNS, kind="mergesort").reset_index(drop=True)


def compute_task1_internal_context_group_metrics(
    task1_snapshot_root: Path,
    eligible_task2_triplets: pd.DataFrame,
) -> pd.DataFrame:
    require_path(task1_snapshot_root, "Task1 snapshot root")
    if eligible_task2_triplets.empty:
        return pd.DataFrame(
            columns=[
                "dataset",
                "cell_line",
                "target",
                "perturbation_type",
                "representation",
                "group_id",
                "n_total",
                "n_A",
                "n_B",
                "n_A_sub",
                "n_B_sub",
                *TASK1_CONTEXT_GROUP_METRICS,
            ]
        )

    eligible_by_dataset: dict[str, pd.DataFrame] = {}
    for dataset, subset in eligible_task2_triplets.groupby("dataset", sort=False):
        eligible_by_dataset[str(dataset)] = subset[["cell_line", "target"]].drop_duplicates().copy()

    rng = init_global_rng(GLOBAL_SEED)
    cohorts, _, _ = build_cohorts(task1_snapshot_root)
    rows: list[dict[str, object]] = []

    # Figure 3 context uses Task1 internal only. We therefore compute split-half
    # Task1 internal group surfaces directly from S1 inputs and never touch Task1 cross.
    for cohort in cohorts:
        dataset = cohort.key.dataset_or_direction
        allowed_subset = eligible_by_dataset.get(str(dataset))
        if allowed_subset is None or allowed_subset.empty:
            continue

        meta = cohort.meta.copy()
        meta["cell_line"] = meta["cell_line"].astype(str)
        meta["target"] = meta["target"].astype(str)
        meta = meta.loc[meta["target"].ne("") & meta["target"].ne("NA")].copy()
        meta = meta.merge(allowed_subset, on=["cell_line", "target"], how="inner", validate="many_to_one")
        if meta.empty:
            continue

        meta["delta_valid_bool"] = meta["delta_valid_bool"].fillna(False).astype(bool)
        meta = meta.loc[meta["delta_valid_bool"]].copy().reset_index(drop=True)
        if meta.empty:
            continue

        meta["group_id"] = [
            make_group_id(cell_line=row.cell_line, target=row.target)
            for row in meta.itertuples(index=False)
        ]

        for group_id, idx in sorted(meta.groupby("group_id", sort=False).groups.items()):
            pos_all = pd.Index(sorted(idx)).to_numpy(dtype="int64")
            row_idx = meta.loc[pos_all, "delta_row_index"].to_numpy(dtype="int64")
            vectors = cohort.vector_store.load_rows(row_idx)
            finite_mask = np.isfinite(vectors).all(axis=1)
            if not finite_mask.any():
                continue

            pos_valid = pos_all[finite_mask]
            vectors_valid = vectors[finite_mask]
            uids_valid = meta.loc[pos_valid, "canonical_query_uid"].astype(str).to_numpy()
            first_row = meta.loc[pos_valid[0]]
            n_total = int(vectors_valid.shape[0])

            row: dict[str, object] = {
                "dataset": dataset,
                "cell_line": str(first_row["cell_line"]),
                "target": str(first_row["target"]),
                "perturbation_type": cohort.key.perturbation_type,
                "representation": cohort.key.representation,
                "group_id": group_id,
                "n_total": n_total,
                "n_A": pd.NA,
                "n_B": pd.NA,
                "n_A_sub": pd.NA,
                "n_B_sub": pd.NA,
                "cosine": pd.NA,
                "pcc": pd.NA,
                "edist": pd.NA,
            }

            if n_total >= 4:
                sort_idx = uids_valid.argsort(kind="mergesort")
                vectors_sorted = vectors_valid[sort_idx]
                uids_sorted = uids_valid[sort_idx]
                perm = rng.permutation(n_total)
                half = n_total // 2
                a_idx = perm[:half]
                b_idx = perm[half:]
                a_vec = vectors_sorted[a_idx]
                b_vec = vectors_sorted[b_idx]
                a_uid = uids_sorted[a_idx]
                b_uid = uids_sorted[b_idx]
                row["n_A"] = int(len(a_vec))
                row["n_B"] = int(len(b_vec))
                row["cosine"] = float(cosine_similarity(a_vec.mean(axis=0), b_vec.mean(axis=0)))
                row["pcc"] = float(pearson_corr(a_vec.mean(axis=0), b_vec.mean(axis=0)))

                if len(a_vec) >= 2 and len(b_vec) >= 2:
                    a_sub, n_a_sub = deterministic_subsample_by_uid(a_vec, a_uid, 256, rng)
                    b_sub, n_b_sub = deterministic_subsample_by_uid(b_vec, b_uid, 256, rng)
                    row["n_A_sub"] = int(n_a_sub)
                    row["n_B_sub"] = int(n_b_sub)
                    if len(a_sub) >= 2 and len(b_sub) >= 2:
                        row["edist"] = float(energy_distance_biascorr(a_sub, b_sub))

            rows.append(row)

    if not rows:
        return pd.DataFrame(
            columns=[
                "dataset",
                "cell_line",
                "target",
                "perturbation_type",
                "representation",
                "group_id",
                "n_total",
                "n_A",
                "n_B",
                "n_A_sub",
                "n_B_sub",
                *TASK1_CONTEXT_GROUP_METRICS,
            ]
        )

    return pd.DataFrame(rows).sort_values(
        ["dataset", "cell_line", "target", "perturbation_type", "representation", "group_id"],
        kind="mergesort",
    ).reset_index(drop=True)


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
    )
    group = group.rename(columns={"target_token": "target"})
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
    group_frame = pd.DataFrame(group_rows)

    retrieval_target = aggregate_task2_c2g_target_means(task2_retrieval_root / "task2_retrieval_per_query.parquet")
    # C2G single-target retrieval rows (query_n_targets==1) are recovered via query_target_tokens.
    # Multi-target C2G queries remain excluded by design (query_n_targets > 1).
    # The empty guard below handles edge cases but is not expected to trigger in normal operation.
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
    retrieval_frame = retrieval_frame[
        [
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
    ]

    combined = pd.concat([group_frame, retrieval_frame], ignore_index=True, sort=False)
    return combined.sort_values(
        ["dataset", "cell_line", "target", "analysis_family", "direction", "representation", "metric_name"],
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)


def build_task1_scope_summary(internal_root: Path, cross_root: Path, output_path: Path) -> Path:
    internal_group = read_csv_required(
        internal_root / "task1_leaderboard_long.csv",
        {
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation",
            "metric_name",
            "n_total",
            "n_valid",
            "n_excluded",
        },
        "Task1 internal leaderboard",
    )
    internal_group = internal_group.loc[internal_group["metric_name"].isin(TASK1_GROUP_METRICS)].copy()
    internal_group_rows = (
        internal_group.groupby(
            ["scope", "dataset_or_direction", "perturbation_type", "representation"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            metric_names=("metric_name", join_unique),
            metric_count=("metric_name", "nunique"),
            n_total=("n_total", "max"),
            n_valid=("n_valid", "max"),
            n_excluded=("n_excluded", "max"),
        )
        .reset_index()
    )
    internal_group_rows["analysis_family"] = "group_concordance"
    internal_group_rows["scope_status"] = "materialized"
    internal_group_rows["cross_alignment_contract"] = pd.NA
    internal_group_rows["n_matched_keys"] = pd.NA
    internal_group_rows["eligible_bool"] = pd.NA
    internal_group_rows["exclusion_reason"] = pd.NA
    internal_group_rows["scope_note"] = "internal lawful slice from S1 group leaderboard"
    internal_group_rows["fm_scope_note"] = internal_group_rows.apply(
        lambda row: task1_fm_scope_note(row["dataset_or_direction"], row["representation"]),
        axis=1,
    )

    internal_retrieval = read_csv_required(
        internal_root / "task1_retrieval_summary.csv",
        {
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation",
            "n_total",
            "n_valid",
            "n_excluded_missing_metric_or_mpos0",
            "N_gallery_max",
        },
        "Task1 internal retrieval summary",
    )
    internal_retrieval_rows = internal_retrieval.copy()
    internal_retrieval_rows["analysis_family"] = "retrieval"
    internal_retrieval_rows["metric_names"] = "mrr|hit1|hit5|hit10"
    internal_retrieval_rows["metric_count"] = len(TASK1_RETRIEVAL_BASES)
    internal_retrieval_rows["scope_status"] = "materialized"
    internal_retrieval_rows["cross_alignment_contract"] = pd.NA
    internal_retrieval_rows["n_matched_keys"] = pd.NA
    internal_retrieval_rows["eligible_bool"] = pd.NA
    internal_retrieval_rows["exclusion_reason"] = pd.NA
    internal_retrieval_rows["scope_note"] = "internal lawful slice from S1 retrieval summary"
    internal_retrieval_rows["fm_scope_note"] = internal_retrieval_rows.apply(
        lambda row: task1_fm_scope_note(row["dataset_or_direction"], row["representation"]),
        axis=1,
    )
    internal_retrieval_rows = internal_retrieval_rows.rename(
        columns={"n_excluded_missing_metric_or_mpos0": "n_excluded"}
    )

    cross_group = read_csv_required(
        cross_root / "task1_cross_leaderboard_long.csv",
        {
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation",
            "metric_name",
            "n_total",
            "n_valid",
            "n_excluded",
            "cross_alignment_contract",
        },
        "Task1 cross leaderboard",
    )
    cross_group = cross_group.loc[cross_group["metric_name"].isin(TASK1_GROUP_METRICS)].copy()
    cross_group_rows = (
        cross_group.groupby(
            ["scope", "dataset_or_direction", "perturbation_type", "representation", "cross_alignment_contract"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            metric_names=("metric_name", join_unique),
            metric_count=("metric_name", "nunique"),
            n_total=("n_total", "max"),
            n_valid=("n_valid", "max"),
            n_excluded=("n_excluded", "max"),
        )
        .reset_index()
    )
    cross_group_rows["analysis_family"] = "group_concordance"
    cross_group_rows["scope_status"] = "materialized"
    cross_group_rows["n_matched_keys"] = pd.NA
    cross_group_rows["eligible_bool"] = pd.NA
    cross_group_rows["exclusion_reason"] = pd.NA
    cross_group_rows["scope_note"] = "cross lawful slice from S2 group leaderboard"
    cross_group_rows["fm_scope_note"] = "FM excluded because Task1 cross is frozen to Gene/Pathway only."

    cross_retrieval = read_csv_required(
        cross_root / "task1_cross_retrieval_summary.csv",
        {
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation",
            "n_total",
            "n_valid",
            "n_excluded_missing_metric_or_mpos0",
            "N_gallery_max",
        },
        "Task1 cross retrieval summary",
    )
    cross_retrieval_rows = cross_retrieval.copy()
    cross_retrieval_rows["analysis_family"] = "retrieval"
    cross_retrieval_rows["metric_names"] = "mrr|hit1|hit5|hit10"
    cross_retrieval_rows["metric_count"] = len(TASK1_RETRIEVAL_BASES)
    cross_retrieval_rows["scope_status"] = "materialized"
    cross_retrieval_rows["cross_alignment_contract"] = "global_idx_lincs + sc_delta_row_idx"
    cross_retrieval_rows["n_matched_keys"] = pd.NA
    cross_retrieval_rows["eligible_bool"] = pd.NA
    cross_retrieval_rows["exclusion_reason"] = pd.NA
    cross_retrieval_rows["scope_note"] = "cross lawful slice from S2 retrieval summary"
    cross_retrieval_rows["fm_scope_note"] = "FM excluded because Task1 cross is frozen to Gene/Pathway only."
    cross_retrieval_rows = cross_retrieval_rows.rename(
        columns={"n_excluded_missing_metric_or_mpos0": "n_excluded"}
    )

    alignment = read_csv_required(
        cross_root / "task1_cross_alignment_proof.csv",
        {
            "cross_alignment_contract",
            "perturbation_type",
            "n_matched_keys",
            "eligible_bool",
            "excluded_reason",
        },
        "Task1 cross alignment proof",
    )
    attrition = read_csv_required(
        cross_root / "task1_cross_attrition.csv",
        {"perturbation_type", "representation", "reason", "n_dropped", "notes"},
        "Task1 cross attrition",
    )
    cross_gate_rows: list[dict[str, object]] = []
    for row in alignment.itertuples(index=False):
        reason_match = attrition.loc[
            attrition["perturbation_type"].eq(row.perturbation_type)
            & attrition["reason"].eq(row.excluded_reason if pd.notna(row.excluded_reason) else ""),
            "notes",
        ]
        gate_note = reason_match.iloc[0] if not reason_match.empty else pd.NA
        cross_gate_rows.append(
            {
                "scope": "cross",
                "dataset_or_direction": "cross",
                "perturbation_type": row.perturbation_type,
                "representation": "ALL",
                "analysis_family": "cross_eligibility",
                "metric_names": "eligibility",
                "metric_count": 1,
                "n_total": pd.NA,
                "n_valid": pd.NA,
                "n_excluded": pd.NA,
                "N_gallery_max": pd.NA,
                "cross_alignment_contract": row.cross_alignment_contract,
                "n_matched_keys": int(row.n_matched_keys),
                "eligible_bool": bool(row.eligible_bool),
                "exclusion_reason": nullable_str(row.excluded_reason),
                "scope_status": "materialized" if bool(row.eligible_bool) else "excluded_by_support_gate",
                "scope_note": gate_note,
                "fm_scope_note": "FM excluded because Task1 cross is frozen to Gene/Pathway only.",
            }
        )
    cross_gate_frame = pd.DataFrame(cross_gate_rows)

    combined = pd.concat(
        [
            internal_group_rows,
            internal_retrieval_rows[
                [
                    "scope",
                    "dataset_or_direction",
                    "perturbation_type",
                    "representation",
                    "analysis_family",
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
                    "scope_status",
                    "scope_note",
                    "fm_scope_note",
                ]
            ],
            cross_group_rows,
            cross_retrieval_rows[
                [
                    "scope",
                    "dataset_or_direction",
                    "perturbation_type",
                    "representation",
                    "analysis_family",
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
                    "scope_status",
                    "scope_note",
                    "fm_scope_note",
                ]
            ],
            cross_gate_frame,
        ],
        ignore_index=True,
        sort=False,
    )
    combined = combined[
        [
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
        ]
    ].sort_values(
        ["scope", "dataset_or_direction", "perturbation_type", "representation", "analysis_family"],
        kind="mergesort",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, index=False)
    return output_path


def build_task1_performance_structure(internal_root: Path, cross_root: Path, output_path: Path) -> Path:
    leaderboard_required = {
        "scope",
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "metric_name",
        "metric_value",
        "n_total",
        "n_valid",
        "n_excluded",
        "cross_alignment_contract",
    }
    summary_required = {
        "scope",
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "n_total",
        "n_valid",
        "n_excluded_missing_metric_or_mpos0",
        "N_gallery_max",
        "mean_mrr_raw",
        "mean_expected_mrr_chance",
        "mean_mrr_corrected",
        "mean_hit1_raw",
        "mean_expected_hit1_chance",
        "mean_hit1_corrected",
        "mean_hit5_raw",
        "mean_expected_hit5_chance",
        "mean_hit5_corrected",
        "mean_hit10_raw",
        "mean_expected_hit10_chance",
        "mean_hit10_corrected",
    }
    chance_required = {
        "scope",
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "abs_delta_mrr",
        "abs_delta_hit1",
        "abs_delta_hit5",
        "abs_delta_hit10",
    }

    internal_group = read_csv_required(
        internal_root / "task1_leaderboard_long.csv",
        leaderboard_required,
        "Task1 internal leaderboard",
    )
    cross_group = read_csv_required(
        cross_root / "task1_cross_leaderboard_long.csv",
        leaderboard_required,
        "Task1 cross leaderboard",
    )
    group_frame = pd.concat([internal_group, cross_group], ignore_index=True)
    group_frame = group_frame.loc[group_frame["metric_name"].isin(TASK1_GROUP_METRICS)].copy()
    group_rows = group_frame.rename(columns={"metric_value": "value"}).copy()
    group_rows["analysis_family"] = "group_concordance"
    group_rows["value_variant"] = "group"
    group_rows["chance_check_available_bool"] = False
    group_rows["chance_abs_delta"] = pd.NA

    internal_summary = read_csv_required(
        internal_root / "task1_retrieval_summary.csv",
        summary_required,
        "Task1 internal retrieval summary",
    )
    cross_summary = read_csv_required(
        cross_root / "task1_cross_retrieval_summary.csv",
        summary_required,
        "Task1 cross retrieval summary",
    )
    retrieval_summary = pd.concat([internal_summary, cross_summary], ignore_index=True)

    internal_chance = read_csv_required(
        internal_root / "task1_chance_identity_check.csv",
        chance_required,
        "Task1 internal chance identity check",
    )
    cross_chance = read_csv_required(
        cross_root / "task1_cross_chance_identity_check.csv",
        chance_required,
        "Task1 cross chance identity check",
    )
    chance = pd.concat([internal_chance, cross_chance], ignore_index=True)
    chance_key = ["scope", "dataset_or_direction", "perturbation_type", "representation"]
    chance = chance.rename(
        columns={
            "abs_delta_mrr": "mrr_chance_abs_delta",
            "abs_delta_hit1": "hit1_chance_abs_delta",
            "abs_delta_hit5": "hit5_chance_abs_delta",
            "abs_delta_hit10": "hit10_chance_abs_delta",
        }
    )
    retrieval_summary = retrieval_summary.merge(chance, on=chance_key, how="left", validate="one_to_one")

    retrieval_rows: list[dict[str, object]] = []
    retrieval_variants = {
        "mrr": {
            "retrieval_raw": "mean_mrr_raw",
            "retrieval_expected": "mean_expected_mrr_chance",
            "retrieval_corrected": "mean_mrr_corrected",
            "chance_abs_delta": "mrr_chance_abs_delta",
        },
        "hit1": {
            "retrieval_raw": "mean_hit1_raw",
            "retrieval_expected": "mean_expected_hit1_chance",
            "retrieval_corrected": "mean_hit1_corrected",
            "chance_abs_delta": "hit1_chance_abs_delta",
        },
        "hit5": {
            "retrieval_raw": "mean_hit5_raw",
            "retrieval_expected": "mean_expected_hit5_chance",
            "retrieval_corrected": "mean_hit5_corrected",
            "chance_abs_delta": "hit5_chance_abs_delta",
        },
        "hit10": {
            "retrieval_raw": "mean_hit10_raw",
            "retrieval_expected": "mean_expected_hit10_chance",
            "retrieval_corrected": "mean_hit10_corrected",
            "chance_abs_delta": "hit10_chance_abs_delta",
        },
    }
    for row in retrieval_summary.itertuples(index=False):
        for metric_name, variant_map in retrieval_variants.items():
            for value_variant, value_column in variant_map.items():
                if value_variant == "chance_abs_delta":
                    continue
                retrieval_rows.append(
                    {
                        "scope": row.scope,
                        "dataset_or_direction": row.dataset_or_direction,
                        "perturbation_type": row.perturbation_type,
                        "representation": row.representation,
                        "analysis_family": "retrieval",
                        "metric_name": metric_name,
                        "value_variant": value_variant,
                        "value": getattr(row, value_column),
                        "n_total": row.n_total,
                        "n_valid": row.n_valid,
                        "n_excluded": row.n_excluded_missing_metric_or_mpos0,
                        "N_gallery_max": row.N_gallery_max,
                        "cross_alignment_contract": pd.NA,
                        "chance_check_available_bool": pd.notna(getattr(row, variant_map["chance_abs_delta"])),
                        "chance_abs_delta": getattr(row, variant_map["chance_abs_delta"]),
                    }
                )
    retrieval_frame = pd.DataFrame(retrieval_rows)
    combined = pd.concat(
        [
            group_rows[
                [
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
                ]
            ],
            retrieval_frame,
        ],
        ignore_index=True,
        sort=False,
    ).sort_values(
        [
            "scope",
            "dataset_or_direction",
            "perturbation_type",
            "representation",
            "analysis_family",
            "metric_name",
            "value_variant",
        ],
        kind="mergesort",
    )
    combined["fm_scope_note"] = combined.apply(
        lambda row: (
            "FM excluded because Task1 cross is frozen to Gene/Pathway only."
            if row["scope"] == "cross"
            else task1_fm_scope_note(row["dataset_or_direction"], row["representation"])
        ),
        axis=1,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, index=False)
    return output_path


def build_task1_internal_to_cross_degradation_summary(manuscript_support_root: Path, output_path: Path) -> Path:
    summary = read_csv_required(
        manuscript_support_root / "group_bridge" / "task1_internal_vs_cross_group_bridge_summary.csv",
        {
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
        },
        "Task1 internal-vs-cross bridge summary",
    )
    summary = summary.rename(columns={"representation_detail": "representation"})
    summary["degradation_scope_note"] = (
        "Downstream degradation summary from Task1 internal-vs-cross manuscript support only; not a Figure 3 bridge."
    )
    summary["fm_scope_note"] = "FM excluded because Task1 cross is frozen to Gene/Pathway only."
    summary = summary.sort_values(
        ["metric_family", "metric_name", "dataset", "cell_line", "representation"],
        kind="mergesort",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, index=False)
    return output_path


def build_task2_scope_summary(
    task2_group_root: Path,
    task2_retrieval_root: Path,
    task2_synthesis_root: Path,
    task2_snapshot_root: Path,
    output_path: Path,
) -> Path:
    coverage = read_csv_required(
        task2_snapshot_root / "task2_pairs_coverage.csv",
        {"dataset", "cell_line", "target_token", "is_eligible_bool"},
        "Task2 pairs coverage",
    )
    registry = read_csv_required(
        task2_snapshot_root / "representation_availability_registry.csv",
        {"dataset", "cell_line", "representation", "availability_status", "availability_reason"},
        "Task2 representation availability registry",
    )
    group_attrition = read_csv_required(
        task2_group_root / "task2_group_attrition.csv",
        {"dataset", "cell_line", "target_token", "representation", "metric_name", "reason"},
        "Task2 group attrition",
    )
    retrieval_summary = read_csv_required(
        task2_retrieval_root / "task2_retrieval_summary.csv",
        {"dataset", "cell_line", "direction", "representation"},
        "Task2 retrieval summary",
    )
    group_leaderboard = read_csv_required(
        task2_synthesis_root / "task2_group_leaderboard.csv",
        {"dataset", "cell_line", "representation", "n_targets_total", "n_attrition_target_rows"},
        "Task2 group leaderboard",
    )
    retrieval_leaderboard = read_csv_required(
        task2_synthesis_root / "task2_retrieval_leaderboard.csv",
        {"dataset", "cell_line", "direction", "representation"},
        "Task2 retrieval leaderboard",
    )

    eligible_counts = (
        coverage.groupby(["dataset", "cell_line", "is_eligible_bool"], dropna=False)
        .size()
        .rename("n_targets")
        .reset_index()
        .pivot_table(
            index=["dataset", "cell_line"],
            columns="is_eligible_bool",
            values="n_targets",
            fill_value=0,
            aggfunc="sum",
        )
        .reset_index()
        .rename(columns={False: "n_targets_ineligible", True: "n_targets_eligible"})
    )
    for required_col in ("n_targets_eligible", "n_targets_ineligible"):
        if required_col not in eligible_counts.columns:
            eligible_counts[required_col] = 0

    attrition_counts = (
        group_attrition.groupby(["dataset", "cell_line", "representation"], dropna=False)
        .size()
        .rename("group_attrition_rows")
        .reset_index()
    )
    retrieval_direction_presence = (
        retrieval_summary.groupby(["dataset", "cell_line", "representation"], dropna=False)["direction"]
        .agg(join_unique)
        .rename("retrieval_directions_present")
        .reset_index()
    )
    c2g_presence = (
        retrieval_leaderboard.loc[retrieval_leaderboard["direction"].eq("C2G")]
        .assign(retrieval_c2g_materialized_bool=True)
        [["dataset", "cell_line", "representation", "retrieval_c2g_materialized_bool"]]
        .drop_duplicates()
    )
    g2c_presence = (
        retrieval_leaderboard.loc[retrieval_leaderboard["direction"].eq("G2C")]
        .assign(retrieval_g2c_materialized_bool=True)
        [["dataset", "cell_line", "representation", "retrieval_g2c_materialized_bool"]]
        .drop_duplicates()
    )
    group_presence = (
        group_leaderboard.assign(group_slice_materialized_bool=True)[
            ["dataset", "cell_line", "representation", "group_slice_materialized_bool", "n_targets_total", "n_attrition_target_rows"]
        ]
        .drop_duplicates()
    )

    out = registry.merge(eligible_counts, on=["dataset", "cell_line"], how="left")
    out = out.merge(attrition_counts, on=["dataset", "cell_line", "representation"], how="left")
    out = out.merge(retrieval_direction_presence, on=["dataset", "cell_line", "representation"], how="left")
    out = out.merge(group_presence, on=["dataset", "cell_line", "representation"], how="left")
    out = out.merge(c2g_presence, on=["dataset", "cell_line", "representation"], how="left")
    out = out.merge(g2c_presence, on=["dataset", "cell_line", "representation"], how="left")
    out["n_targets_eligible"] = out["n_targets_eligible"].fillna(0).astype(int)
    out["n_targets_ineligible"] = out["n_targets_ineligible"].fillna(0).astype(int)
    out["group_attrition_rows"] = out["group_attrition_rows"].fillna(0).astype(int)
    out["group_slice_materialized_bool"] = out["group_slice_materialized_bool"].fillna(False).astype(bool)
    out["retrieval_c2g_materialized_bool"] = out["retrieval_c2g_materialized_bool"].fillna(False).astype(bool)
    out["retrieval_g2c_materialized_bool"] = out["retrieval_g2c_materialized_bool"].fillna(False).astype(bool)
    out["scope_note"] = out.apply(
        lambda row: (
            "representation absent by scope policy"
            if row["availability_status"] == "not_applicable_scope"
            else "representation available under corrected Task2 scope"
        ),
        axis=1,
    )
    out["fm_scope_note"] = out.apply(
        lambda row: task2_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    out = out[
        [
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
        ]
    ].sort_values(["dataset", "cell_line", "representation"], kind="mergesort")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    return output_path


def build_task2_direction_summary_file(
    task2_retrieval_root: Path,
    manuscript_phase1_root: Path,
    output_path: Path,
) -> Path:
    summary = read_csv_required(
        task2_retrieval_root / "task2_retrieval_summary.csv",
        {
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
        },
        "Task2 retrieval summary",
    )
    chance = read_csv_required(
        task2_retrieval_root / "task2_chance_identity_check.csv",
        {
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "abs_delta_mrr",
            "abs_delta_hit1",
            "abs_delta_hit5",
            "abs_delta_hit10",
        },
        "Task2 retrieval chance identity check",
    )
    a3 = read_csv_required(
        manuscript_phase1_root / "a3_direction_robustness" / "task2_direction_robustness_audit.csv",
        {
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "status",
            "status_reason",
            "chance_check_available_bool",
            "chance_check_max_abs_delta",
        },
        "A3 direction robustness audit",
    )
    chance["chance_check_max_abs_delta"] = chance[
        ["abs_delta_mrr", "abs_delta_hit1", "abs_delta_hit5", "abs_delta_hit10"]
    ].max(axis=1)
    chance["chance_check_available_bool"] = True
    chance = chance[
        [
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "chance_check_available_bool",
            "chance_check_max_abs_delta",
        ]
    ]
    out = summary.merge(
        chance,
        on=["dataset", "cell_line", "direction", "representation"],
        how="left",
        validate="one_to_one",
    )
    out = out.merge(
        a3.rename(columns={"status": "a3_status", "status_reason": "a3_status_reason"})[
            [
                "dataset",
                "cell_line",
                "direction",
                "representation",
                "a3_status",
                "a3_status_reason",
            ]
        ],
        on=["dataset", "cell_line", "direction", "representation"],
        how="left",
        validate="one_to_one",
    )
    out["chance_check_available_bool"] = out["chance_check_available_bool"].fillna(False).astype(bool)
    out["fm_scope_note"] = out.apply(
        lambda row: task2_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    out = out.sort_values(["dataset", "cell_line", "direction", "representation"], kind="mergesort")
    # Flag current-phase primary direction: C2G is the main axis; G2C is supporting-only.
    out["current_phase_primary_direction_bool"] = out["direction"].eq("C2G")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    return output_path


def build_task2_performance_structure_file(task2_synthesis_root: Path, output_path: Path) -> Path:
    group_leaderboard = read_csv_required(
        task2_synthesis_root / "task2_group_leaderboard.csv",
        {"dataset", "cell_line", "representation"},
        "Task2 group leaderboard",
    )
    retrieval_leaderboard = read_csv_required(
        task2_synthesis_root / "task2_retrieval_leaderboard.csv",
        {"dataset", "cell_line", "direction", "representation"},
        "Task2 retrieval leaderboard",
    )
    benchmark = read_csv_required(
        task2_synthesis_root / "task2_benchmark_summary_long.csv",
        {
            "analysis_family",
            "dataset",
            "cell_line",
            "direction",
            "representation",
            "metric_name",
            "metric_value",
            "rank_value",
            "rank_basis_metric_name",
            "leaderboard_eligible_bool",
            "cross_representation_comparable_bool",
            "source_table",
            "source_stage",
            "source_run_id",
        },
        "Task2 benchmark summary long",
    )
    benchmark = benchmark.loc[
        benchmark["analysis_family"].isin(TASK2_BENCHMARK_ALLOWED_ANALYSIS_FAMILIES)
    ].copy()
    group_keys = set(
        tuple(row)
        for row in group_leaderboard[["dataset", "cell_line", "representation"]].drop_duplicates().itertuples(
            index=False, name=None
        )
    )
    retrieval_keys = set(
        tuple(row)
        for row in retrieval_leaderboard[
            ["dataset", "cell_line", "direction", "representation"]
        ].drop_duplicates().itertuples(index=False, name=None)
    )

    def has_source_row(row: pd.Series) -> bool:
        if row["analysis_family"] == "group_concordance":
            return (row["dataset"], row["cell_line"], row["representation"]) in group_keys
        return (row["dataset"], row["cell_line"], row["direction"], row["representation"]) in retrieval_keys

    benchmark["source_key_present_bool"] = benchmark.apply(has_source_row, axis=1)
    benchmark["performance_scope_note"] = benchmark["analysis_family"].map(
        {
            "group_concordance": "group leaderboard row from S6 benchmark summary long",
            "retrieval": "direction-specific retrieval row from S6 benchmark summary long",
        }
    )
    benchmark["fm_scope_note"] = benchmark.apply(
        lambda row: task2_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    # C2G-only current manuscript mode: keep group_concordance (direction=NA) and C2G
    # retrieval rows only; G2C retrieval rows are excluded from this canonical object.
    benchmark = benchmark.loc[
        benchmark["analysis_family"].ne("retrieval") | benchmark["direction"].eq("C2G")
    ].copy()
    benchmark = benchmark.sort_values(
        ["analysis_family", "dataset", "cell_line", "direction", "representation", "metric_name"],
        kind="mergesort",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    benchmark.to_csv(output_path, index=False)
    return output_path


def build_task2_direction_performance_summary(
    task2_retrieval_root: Path,
    task2_synthesis_root: Path,
    manuscript_phase1_root: Path,
    analysis_root: Path,
) -> tuple[Path, Path]:
    direction_path = build_task2_direction_summary_file(
        task2_retrieval_root,
        manuscript_phase1_root,
        analysis_root / "figure3_task2_direction_support_summary.csv",
    )
    performance_path = build_task2_performance_structure_file(
        task2_synthesis_root,
        analysis_root / "figure3_task2_performance_structure.csv",
    )
    return direction_path, performance_path


# Historical helper retained for explicit non-canonical inventory only.
# It must not be wired into the canonical registry for the frozen current phase.
def build_task2_c2g_covariate_summary(manuscript_history_root: Path, output_path: Path) -> Path:
    summary = read_csv_required(
        manuscript_history_root / "sensitivity" / "task2_c2g_dose_time_sensitivity_summary.csv",
        {
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
        },
        "Task2 C2G dose/time sensitivity summary",
    )
    summary = summary.rename(columns={"retrieval_metric_name": "metric_name"})
    summary["covariate_scope_note"] = "C2G-only downstream covariate summary from existing B6 manuscript support."
    summary["fm_scope_note"] = summary.apply(
        lambda row: task2_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    summary = summary.sort_values(
        ["dataset", "cell_line", "target", "metric_name", "dose_value", "time", "representation"],
        kind="mergesort",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, index=False)
    return output_path


def build_task2_cell_line_pattern_summary(
    task2_group_root: Path,
    task2_retrieval_root: Path,
    output_path: Path,
) -> Path:
    support_long = build_task2_target_support_long(task2_group_root, task2_retrieval_root)
    out = (
        support_long.groupby(
            ["dataset", "cell_line", "analysis_family", "direction", "representation", "metric_name"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_targets=("target", "nunique"),
            n_source_rows=("target", "size"),
            mean_metric_value=("metric_value", "mean"),
            median_metric_value=("metric_value", "median"),
            min_metric_value=("metric_value", "min"),
            max_metric_value=("metric_value", "max"),
            n_queries_total=("n_queries", sum_min_count),
            n_queries_mean=("n_queries", "mean"),
            n_chem_instances_used_total=("n_chem_instances_used", sum_min_count),
            n_gen_instances_used_total=("n_gen_instances_used", sum_min_count),
            n_chem_sub_total=("n_chem_sub", sum_min_count),
            n_gen_sub_total=("n_gen_sub", sum_min_count),
        )
        .reset_index()
    )
    out["pattern_summary_scope"] = "cell_line_over_canonical_target_aggregation"
    out["pattern_note"] = out["analysis_family"].map(
        {
            "group_concordance": "Aggregated across canonical Task2 target_token rows within each cell_line.",
            "retrieval": "C2G single-target retrieval rows present (query_n_targets==1); multi-target C2G queries excluded by design.",
        }
    )
    out["fm_scope_note"] = out.apply(
        lambda row: task2_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    out = out.sort_values(
        ["dataset", "cell_line", "analysis_family", "direction", "representation", "metric_name"],
        kind="mergesort",
        na_position="last",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    return output_path


def build_task2_target_pattern_summary(
    task2_group_root: Path,
    task2_retrieval_root: Path,
    output_path: Path,
) -> Path:
    support_long = build_task2_target_support_long(task2_group_root, task2_retrieval_root)
    out = (
        support_long.groupby(
            ["dataset", "target", "analysis_family", "direction", "representation", "metric_name"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_cell_lines=("cell_line", "nunique"),
            n_source_rows=("cell_line", "size"),
            mean_metric_value=("metric_value", "mean"),
            median_metric_value=("metric_value", "median"),
            min_metric_value=("metric_value", "min"),
            max_metric_value=("metric_value", "max"),
            n_queries_total=("n_queries", sum_min_count),
            n_queries_mean=("n_queries", "mean"),
            n_chem_instances_used_total=("n_chem_instances_used", sum_min_count),
            n_gen_instances_used_total=("n_gen_instances_used", sum_min_count),
            n_chem_sub_total=("n_chem_sub", sum_min_count),
            n_gen_sub_total=("n_gen_sub", sum_min_count),
        )
        .reset_index()
    )
    out["pattern_summary_scope"] = "target_over_cell_line_aggregation"
    out["pattern_note"] = out["analysis_family"].map(
        {
            "group_concordance": "Aggregated across Task2 cell_line rows for each canonical target_token.",
            "retrieval": "C2G single-target retrieval rows present (query_n_targets==1); multi-target C2G queries excluded by design.",
        }
    )
    out["fm_scope_note"] = out.apply(
        lambda row: task2_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    out = out.sort_values(
        ["dataset", "target", "analysis_family", "direction", "representation", "metric_name"],
        kind="mergesort",
        na_position="last",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    return output_path


def build_task1_internal_contextual_support_summary(
    task1_snapshot_root: Path,
    task1_internal_root: Path,
    task2_snapshot_root: Path,
    output_path: Path,
) -> Path:
    eligible_task2_triplets = load_eligible_task2_triplets(task2_snapshot_root)

    retrieval_summary = aggregate_task1_internal_context_retrieval_means(
        task1_internal_root / "task1_retrieval_per_query.parquet",
        eligible_task2_triplets,
    )
    group_summary = compute_task1_internal_context_group_metrics(
        task1_snapshot_root,
        eligible_task2_triplets,
    )

    if retrieval_summary.empty and group_summary.empty:
        raise ValueError("Task1 internal contextual support summary produced no metadata-aligned rows.")

    retrieval_long = retrieval_summary.melt(
        id_vars=[
            "dataset",
            "cell_line",
            "target",
            "perturbation_type",
            "representation",
            "n_queries",
            "N_gallery",
            "m_pos",
        ],
        value_vars=list(TASK1_RETRIEVAL_METRICS),
        var_name="metric_name",
        value_name="metric_value",
    )
    retrieval_long["analysis_family"] = "retrieval"
    retrieval_long["group_id"] = pd.NA
    retrieval_long["n_total"] = pd.NA
    retrieval_long["n_A"] = pd.NA
    retrieval_long["n_B"] = pd.NA
    retrieval_long["n_A_sub"] = pd.NA
    retrieval_long["n_B_sub"] = pd.NA
    retrieval_long["N_gallery_mean"] = retrieval_long["N_gallery"]
    retrieval_long["m_pos_mean"] = retrieval_long["m_pos"]
    retrieval_long["context_support_scope"] = (
        "task1_internal_only_metadata_aligned_to_task2_evaluable_triplets"
    )
    retrieval_long["context_note"] = (
        "Task1 internal retrieval context only. Rows are filtered to eligible Task2 "
        "(dataset, cell_line, target) metadata keys without Task1 cross and without "
        "any Task1/Task2 pairwise bridge construction."
    )
    retrieval_long["fm_scope_note"] = retrieval_long.apply(
        lambda row: task1_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    retrieval_long = retrieval_long[
        TASK1_CONTEXT_GROUP_OUTPUT_COLUMNS
    ]

    group_long = group_summary.melt(
        id_vars=[
            "dataset",
            "cell_line",
            "target",
            "perturbation_type",
            "representation",
            "group_id",
            "n_total",
            "n_A",
            "n_B",
            "n_A_sub",
            "n_B_sub",
        ],
        value_vars=list(TASK1_CONTEXT_GROUP_METRICS),
        var_name="metric_name",
        value_name="metric_value",
    )
    group_long["metric_name"] = group_long["metric_name"].map(TASK1_CONTEXT_GROUP_METRIC_NAME_MAP)
    group_long["analysis_family"] = "group_concordance"
    group_long["n_queries"] = pd.NA
    group_long["N_gallery_mean"] = pd.NA
    group_long["m_pos_mean"] = pd.NA
    group_long["context_support_scope"] = (
        "task1_internal_only_metadata_aligned_to_task2_evaluable_triplets"
    )
    group_long["context_note"] = (
        "Task1 internal group context only. Split-half group metrics are filtered to eligible Task2 "
        "(dataset, cell_line, target) metadata keys without Task1 cross and without "
        "any Task1/Task2 pairwise bridge construction."
    )
    group_long["fm_scope_note"] = group_long.apply(
        lambda row: task1_fm_scope_note(row["dataset"], row["representation"]),
        axis=1,
    )
    group_long = group_long[
        TASK1_CONTEXT_GROUP_OUTPUT_COLUMNS
    ]

    long_frame = pd.concat([group_long, retrieval_long], ignore_index=True, sort=False)
    long_frame = long_frame.sort_values(
        [
            "dataset",
            "cell_line",
            "target",
            "perturbation_type",
            "representation",
            "analysis_family",
            "metric_name",
        ],
        kind="mergesort",
    ).reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    long_frame.to_csv(output_path, index=False)
    return output_path


def materialization_status_for_path(path_str: str) -> str:
    return "present" if Path(path_str).exists() else "absent"


def annotate_canonical_objects(canonical_objects: list[dict[str, object]]) -> list[dict[str, object]]:
    annotated: list[dict[str, object]] = []
    for entry in canonical_objects:
        record = dict(entry)
        output_files = [str(path) for path in record.get("output_files", [])]
        if len(output_files) == 1:
            record["canonical_filename"] = Path(output_files[0]).name
        record["output_materialization"] = [
            {
                "path": path,
                "status": materialization_status_for_path(path),
            }
            for path in output_files
        ]
        annotated.append(record)
    return annotated


def build_historical_noncanonical_objects(manuscript_history_root: Path) -> list[dict[str, object]]:
    objects: list[dict[str, object]] = []
    for candidate in HISTORICAL_NONCANONICAL_CANDIDATES:
        path = manuscript_history_root / candidate["relative_path"]
        objects.append(
            {
                "object_name": candidate["object_name"],
                "path": str(path),
                "reason_no_longer_canonical": candidate["reason_no_longer_canonical"],
                "disposition": "historical_only",
                "materialization_status": "present" if path.exists() else "absent",
            }
        )
    return objects


def build_framework_manifest(
    output_path: Path,
    output_root: Path,
    source_roots: dict[str, str],
    canonical_objects: list[dict[str, object]],
    historical_noncanonical_objects: list[dict[str, object]],
) -> Path:
    canonical_payload = annotate_canonical_objects(canonical_objects)
    payload = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "output_root": str(output_root),
        "source_roots": source_roots,
        "contract_notes": MANIFEST_CONTRACT_NOTES,
        "canonical_objects": canonical_payload,
        "historical_noncanonical_objects": historical_noncanonical_objects,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return output_path


def build_manuscript_comparison_statistics_manifest_entry(output_path: Path) -> dict[str, object]:
    return {
        "object_name": "Manuscript comparison statistics",
        "manuscript_questions": [
            "C2_F2a_internal_group_lincs",
            "C2_F2b_internal_retrieval_lincs",
            "C2_F2a_internal_group_scperturb_common_scope",
            "C2_F2a_internal_retrieval_scperturb_common_scope",
            "C2_F2c_internal_group_scperturb_fm",
            "C2_F2d_internal_retrieval_scperturb_fm",
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
            str(DEFAULT_TASK1_SNAPSHOT_ROOT),
            str(DEFAULT_TASK1_INTERNAL_ROOT / "task1_retrieval_per_query.parquet"),
            str(DEFAULT_TASK2_RETRIEVAL_ROOT / "task2_retrieval_per_query.parquet"),
            str(DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "group_bridge" / "task1_internal_vs_cross_group_bridge_summary.csv"),
            str(DEFAULT_MANUSCRIPT_SUPPORT_ROOT / "group_bridge" / "task1_internal_vs_cross_group_bridge.csv"),
            str(DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis" / "figure3_task2_direction_support_summary.csv"),
            str(DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis" / "figure3_task2_performance_structure.csv"),
            str(DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis" / "figure3_task2_cell_line_pattern_summary.csv"),
            str(DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis" / "figure3_task2_target_pattern_summary.csv"),
            str(DEFAULT_MANUSCRIPT_ACTIVE_ROOT / "analysis" / "figure3_task1_internal_contextual_support_summary.csv"),
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
        "included_logic": "Formal paired manuscript comparison tests with family-local BH correction.",
        "excluded_logic": "No group-vs-retrieval formal comparisons and no unsupported target-anchored C2G comparison family.",
        "fm_scope": "Figure 2 scPerturb-internal FM plus Figure 3 scPerturb-local FM only.",
        "justification": "Comparative manuscript claims require formal statistical backing under the frozen current-phase contract.",
    }


def main() -> int:
    args = parse_args()
    task1_snapshot_root = args.task1_snapshot_root.resolve()
    task1_internal_root = args.task1_internal_root.resolve()
    task1_cross_root = args.task1_cross_root.resolve()
    task2_group_root = args.task2_group_root.resolve()
    task2_retrieval_root = args.task2_retrieval_root.resolve()
    task2_synthesis_root = args.task2_synthesis_root.resolve()
    manuscript_active_root = args.manuscript_active_root.resolve()
    manuscript_support_root = args.manuscript_support_root.resolve()
    manuscript_history_root = args.manuscript_history_root.resolve()
    manuscript_phase1_root = args.manuscript_phase1_root.resolve()
    task2_snapshot_root = args.task2_snapshot_root.resolve()

    analysis_root = manuscript_active_root / "analysis"

    require_path(
        manuscript_support_root / "group_bridge" / "task1_internal_vs_cross_group_bridge_summary.csv",
        "A1 Task1 internal-vs-cross bridge summary",
    )
    require_path(
        manuscript_phase1_root / "a3_direction_robustness" / "task2_direction_robustness_audit.csv",
        "A3 direction robustness audit",
    )

    canonical_objects: list[dict[str, object]] = []

    figure2_scope_path = build_task1_scope_summary(
        task1_internal_root,
        task1_cross_root,
        analysis_root / "figure2_task1_scope_summary.csv",
    )
    canonical_objects.append(
        {
            "object_name": "Task1 scope summary",
            "manuscript_questions": [
                "Figure 2 Task1 evaluable scope",
            ],
            "source_inputs": [
                str(task1_internal_root / "task1_leaderboard_long.csv"),
                str(task1_internal_root / "task1_retrieval_summary.csv"),
                str(task1_cross_root / "task1_cross_leaderboard_long.csv"),
                str(task1_cross_root / "task1_cross_retrieval_summary.csv"),
                str(task1_cross_root / "task1_cross_alignment_proof.csv"),
                str(task1_cross_root / "task1_cross_attrition.csv"),
            ],
            "output_files": [str(figure2_scope_path)],
            "key_columns": ["scope", "dataset_or_direction", "perturbation_type", "representation", "analysis_family"],
            "representation_scope": "Task1 internal keeps Gene/Pathway/FM where supported; Task1 cross keeps Gene/Pathway only.",
            "included_logic": "Internal and cross evaluable scope rows, plus explicit cross eligibility gate rows.",
            "excluded_logic": "No Task1/Task2 bridge logic.",
            "fm_scope": "Internal scPerturb keeps FM; Task1 cross excludes FM by frozen scope.",
            "justification": "Figure 2 needs an explicit Task1 evaluable scope object with the frozen FM and cross-alignment boundaries.",
        }
    )

    figure2_perf_path = build_task1_performance_structure(
        task1_internal_root,
        task1_cross_root,
        analysis_root / "figure2_task1_performance_structure.csv",
    )
    canonical_objects.append(
        {
            "object_name": "Task1 performance-structure summary",
            "manuscript_questions": [
                "Figure 2 Task1 performance structure",
            ],
            "source_inputs": [
                str(task1_internal_root / "task1_leaderboard_long.csv"),
                str(task1_internal_root / "task1_retrieval_summary.csv"),
                str(task1_internal_root / "task1_chance_identity_check.csv"),
                str(task1_cross_root / "task1_cross_leaderboard_long.csv"),
                str(task1_cross_root / "task1_cross_retrieval_summary.csv"),
                str(task1_cross_root / "task1_cross_chance_identity_check.csv"),
            ],
            "output_files": [str(figure2_perf_path)],
            "key_columns": [
                "scope",
                "dataset_or_direction",
                "perturbation_type",
                "representation",
                "analysis_family",
                "metric_name",
                "value_variant",
            ],
            "representation_scope": "Task1 internal and cross performance rows within the frozen Task1 scope.",
            "included_logic": "Group and retrieval performance rows with chance-check context.",
            "excluded_logic": "No Figure 3 contextual support and no Task1/Task2 bridge logic.",
            "fm_scope": "Internal scPerturb keeps FM; Task1 cross excludes FM by frozen scope.",
            "justification": "Figure 2 needs the within-Task1 performance structure without introducing any Task1/Task2 bridge surface.",
        }
    )

    figure2_deg_path = build_task1_internal_to_cross_degradation_summary(
        manuscript_support_root,
        analysis_root / "figure2_task1_internal_to_cross_degradation_summary.csv",
    )
    canonical_objects.append(
        {
            "object_name": "Task1 internal-to-cross degradation summary",
            "manuscript_questions": [
                "Figure 2 Task1 internal-to-cross degradation",
            ],
            "source_inputs": [
                str(manuscript_support_root / "group_bridge" / "task1_internal_vs_cross_group_bridge_summary.csv"),
            ],
            "output_files": [str(figure2_deg_path)],
            "key_columns": ["metric_family", "metric_name", "dataset", "cell_line", "representation"],
            "representation_scope": "Shared Task1 internal/cross Gene/Pathway representation surface only.",
            "included_logic": "A1 internal-vs-cross degradation summary only.",
            "excluded_logic": "No Figure 3 bridge usage and no FM rows.",
            "fm_scope": "FM excluded because Task1 cross excludes FM.",
            "justification": "Figure 2 needs a downstream degradation object anchored to the existing A1 internal-vs-cross summary only.",
        }
    )

    figure3_scope_path = build_task2_scope_summary(
        task2_group_root,
        task2_retrieval_root,
        task2_synthesis_root,
        task2_snapshot_root,
        analysis_root / "figure3_task2_scope_summary.csv",
    )
    canonical_objects.append(
        {
            "object_name": "Task2 scope summary",
            "manuscript_questions": [
                "Figure 3 Task2 evaluable scope",
            ],
            "source_inputs": [
                str(task2_snapshot_root / "task2_pairs_coverage.csv"),
                str(task2_snapshot_root / "representation_availability_registry.csv"),
                str(task2_group_root / "task2_group_attrition.csv"),
                str(task2_retrieval_root / "task2_retrieval_summary.csv"),
                str(task2_synthesis_root / "task2_group_leaderboard.csv"),
                str(task2_synthesis_root / "task2_retrieval_leaderboard.csv"),
            ],
            "output_files": [str(figure3_scope_path)],
            "key_columns": ["dataset", "cell_line", "representation"],
            "representation_scope": "Task2 lawful scope by dataset, cell_line, representation.",
            "included_logic": "Coverage, availability, attrition, and materialized direction presence.",
            "excluded_logic": "No bridge or target-pattern logic.",
            "fm_scope": "scPerturb-only scope keeps FM; LINCS scope excludes FM.",
            "justification": "Figure 3 needs an explicit Task2 evaluable-scope object with frozen dataset and FM boundaries.",
        }
    )

    figure3_direction_path, figure3_performance_path = build_task2_direction_performance_summary(
        task2_retrieval_root,
        task2_synthesis_root,
        manuscript_phase1_root,
        analysis_root,
    )
    canonical_objects.append(
        {
            "object_name": "Task2 direction support summary",
            "manuscript_questions": [
                "Figure 3 Task2 directionality",
            ],
            "source_inputs": [
                str(task2_retrieval_root / "task2_retrieval_summary.csv"),
                str(task2_retrieval_root / "task2_chance_identity_check.csv"),
                str(manuscript_phase1_root / "a3_direction_robustness" / "task2_direction_robustness_audit.csv"),
            ],
            "output_files": [str(figure3_direction_path)],
            "key_columns": [
                "dataset",
                "cell_line",
                "direction",
                "representation",
                "gallery_definition_id",
                "pos_definition_id",
            ],
            "representation_scope": "Direction-specific Task2 retrieval support within the lawful Task2 representation surface.",
            "included_logic": "Direction-specific retrieval support with corrected metrics, denominator context, and reviewed A3 status.",
            "excluded_logic": "No pairwise Task1/Task2 bridge logic.",
            "fm_scope": "scPerturb-only rows keep FM; LINCS rows exclude FM.",
            "justification": "Figure 3 needs an explicit direction-support object separate from the broader Task2 performance structure.",
        }
    )
    canonical_objects.append(
        {
            "object_name": "Task2 performance structure",
            "manuscript_questions": [
                "Figure 3 Task2 performance structure",
            ],
            "source_inputs": [
                str(task2_synthesis_root / "task2_group_leaderboard.csv"),
                str(task2_synthesis_root / "task2_retrieval_leaderboard.csv"),
                str(task2_synthesis_root / "task2_benchmark_summary_long.csv"),
            ],
            "output_files": [str(figure3_performance_path)],
            "key_columns": [
                "analysis_family",
                "dataset",
                "cell_line",
                "direction",
                "representation",
                "metric_name",
            ],
            "representation_scope": "Task2 group and retrieval performance within the lawful Task2 representation surface.",
            "included_logic": "Benchmark performance structure rows only.",
            "excluded_logic": "No pairwise Task1/Task2 bridge logic and no standalone C2G experiment object.",
            "fm_scope": "scPerturb-only rows keep FM; LINCS rows exclude FM.",
            "justification": "Figure 3 needs a benchmark performance-structure object separate from the direction-support table.",
        }
    )

    figure3_cell_line_path = build_task2_cell_line_pattern_summary(
        task2_group_root,
        task2_retrieval_root,
        analysis_root / "figure3_task2_cell_line_pattern_summary.csv",
    )
    canonical_objects.append(
        {
            "object_name": "Task2 cell_line pattern summary",
            "manuscript_questions": [
                "Figure 3 Task2 cell_line pattern summary",
            ],
            "source_inputs": [
                str(task2_group_root / "task2_group_concordance.csv"),
                str(task2_retrieval_root / "task2_retrieval_per_query.parquet"),
            ],
            "output_files": [str(figure3_cell_line_path)],
            "key_columns": ["dataset", "cell_line", "analysis_family", "direction", "representation", "metric_name"],
            "representation_scope": "Target-anchored Task2 group rows plus C2G single-target retrieval rows.",
            "included_logic": "Aggregation over canonical target rows within each cell_line.",
            "excluded_logic": "G2C excluded in current C2G-only manuscript phase; no Task1/Task2 bridge logic.",
            "fm_scope": "scPerturb-only rows keep FM; LINCS rows exclude FM.",
            "justification": "Figure 3 needs a cell_line summary built only from canonical Task2 target-anchored surfaces.",
        }
    )

    figure3_target_path = build_task2_target_pattern_summary(
        task2_group_root,
        task2_retrieval_root,
        analysis_root / "figure3_task2_target_pattern_summary.csv",
    )
    canonical_objects.append(
        {
            "object_name": "Task2 target pattern summary",
            "manuscript_questions": [
                "Figure 3 Task2 target pattern summary",
            ],
            "source_inputs": [
                str(task2_group_root / "task2_group_concordance.csv"),
                str(task2_retrieval_root / "task2_retrieval_per_query.parquet"),
            ],
            "output_files": [str(figure3_target_path)],
            "key_columns": ["dataset", "target", "analysis_family", "direction", "representation", "metric_name"],
            "representation_scope": "Target-anchored Task2 group rows plus C2G single-target retrieval rows.",
            "included_logic": "Aggregation over cell_line rows for each canonical target.",
            "excluded_logic": "G2C excluded in current C2G-only manuscript phase; no Task1/Task2 bridge logic.",
            "fm_scope": "scPerturb-only rows keep FM; LINCS rows exclude FM.",
            "justification": "Figure 3 needs a target summary built only from canonical Task2 target-anchored surfaces.",
        }
    )

    figure3_task1_context_path = build_task1_internal_contextual_support_summary(
        task1_snapshot_root,
        task1_internal_root,
        task2_snapshot_root,
        analysis_root / "figure3_task1_internal_contextual_support_summary.csv",
    )
    canonical_objects.append(
        {
            "object_name": "Task1 internal contextual support summary",
            "manuscript_questions": [
                "Figure 3 Task1 internal contextual support",
            ],
            "source_inputs": [
                str(task1_snapshot_root),
                str(task1_internal_root / "task1_retrieval_per_query.parquet"),
                str(task2_snapshot_root / "task2_pairs_coverage.csv"),
            ],
            "output_files": [str(figure3_task1_context_path)],
            "key_columns": [
                "dataset",
                "cell_line",
                "target",
                "perturbation_type",
                "representation",
                "analysis_family",
                "metric_name",
            ],
            "representation_scope": "Task1 internal-only group and retrieval surfaces aligned to eligible Task2 biological metadata keys.",
            "included_logic": "Metadata-aligned Task1 internal group surface plus Task1 internal retrieval surface only.",
            "excluded_logic": "No Task1 cross rows and no Task1/Task2 pairwise bridge logic.",
            "fm_scope": "scPerturb internal keeps FM; LINCS internal excludes FM.",
            "justification": "Figure 3 discussion now requires Task1 internal context from both group and retrieval surfaces without recreating any bridge object.",
        }
    )

    stats_path = analysis_root / ANALYSIS_OUTPUT_FILENAMES["manuscript_comparison_statistics"]
    if stats_path.exists():
        canonical_objects.append(build_manuscript_comparison_statistics_manifest_entry(stats_path))

    historical_noncanonical_objects = build_historical_noncanonical_objects(manuscript_history_root)

    manifest_path = build_framework_manifest(
        output_path=analysis_root / ANALYSIS_OUTPUT_FILENAMES["framework_analysis_manifest"],
        output_root=analysis_root,
        source_roots={
            "task1_snapshot_root": str(task1_snapshot_root),
            "task1_internal_root": str(task1_internal_root),
            "task1_cross_root": str(task1_cross_root),
            "task2_group_root": str(task2_group_root),
            "task2_retrieval_root": str(task2_retrieval_root),
            "task2_synthesis_root": str(task2_synthesis_root),
            "task2_snapshot_root": str(task2_snapshot_root),
            "manuscript_active_root": str(manuscript_active_root),
            "manuscript_support_root": str(manuscript_support_root),
            "manuscript_history_root": str(manuscript_history_root),
            "manuscript_phase1_root": str(manuscript_phase1_root),
        },
        canonical_objects=canonical_objects,
        historical_noncanonical_objects=historical_noncanonical_objects,
    )

    print(f"Wrote {len(canonical_objects)} canonical manuscript analysis objects under {analysis_root}")
    for entry in canonical_objects:
        for path in entry["output_files"]:
            print(path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
