#!/usr/bin/env python3
"""
Build focused Task1 high-concordance enrichment summaries for Figure 2 support.

Status:
- support-only manuscript builder

Consumed by:
- Figure 2 cell-line and target high-concordance result families

Architecture:
- one focused local support script; does not modify S1/S2 benchmark logic

This script reads the existing Task1 internal and cross retrieval per-query
outputs, applies a frozen success definition, and materializes enrichment-style
summary tables by cell line and target.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_INTERNAL_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_per_query.parquet"
)
DEFAULT_CROSS_PER_QUERY_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_per_query.parquet"
)
DEFAULT_ANALYSIS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/analysis")
DEFAULT_CELL_LINE_OUTPUT_PATH = DEFAULT_ANALYSIS_ROOT / "figure2_task1_cell_line_high_concordance_summary.csv"
DEFAULT_TARGET_OUTPUT_PATH = DEFAULT_ANALYSIS_ROOT / "figure2_task1_target_high_concordance_summary.csv"

FAMILY_COLUMNS = ["scope", "dataset_or_direction", "perturbation_type", "representation"]
SUCCESS_DEFINITION = "(mrr_corrected > 0) AND (hit10_corrected > 0)"
REQUIRED_COLUMNS = [
    *FAMILY_COLUMNS,
    "cell_line",
    "target_token",
    "mrr_corrected",
    "hit10_corrected",
]
OUTPUT_COLUMNS_CELL_LINE = [
    *FAMILY_COLUMNS,
    "cell_line",
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
]
OUTPUT_COLUMNS_TARGET = [
    *FAMILY_COLUMNS,
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
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Task1 high-concordance summaries by cell line and target.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--internal-per-query-path", type=Path, default=DEFAULT_INTERNAL_PER_QUERY_PATH)
    parser.add_argument("--cross-per-query-path", type=Path, default=DEFAULT_CROSS_PER_QUERY_PATH)
    parser.add_argument("--cell-line-output-path", type=Path, default=DEFAULT_CELL_LINE_OUTPUT_PATH)
    parser.add_argument("--target-output-path", type=Path, default=DEFAULT_TARGET_OUTPUT_PATH)
    return parser.parse_args()


def resolve_path(project_root: Path, raw_path: Path) -> Path:
    return raw_path if raw_path.is_absolute() else (project_root / raw_path)


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")


def split_target_tokens(value: object) -> list[str]:
    text = str(value).strip()
    if text in {"", "NA"}:
        return []
    parts: list[str] = []
    for chunk in text.split("|"):
        parts.extend(piece.strip() for piece in chunk.split("_"))
    return sorted({part for part in parts if part and part != "NA"})


def canonicalize_target_token(value: object) -> str:
    cleaned = split_target_tokens(value)
    if not cleaned:
        return str(value).strip()
    if len(cleaned) == 1:
        return cleaned[0]
    return "|".join(cleaned)


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


def update_entity_counts(
    entity_counts: dict[tuple[str, ...], list[int]],
    frame: pd.DataFrame,
    entity_column: str,
) -> None:
    if frame.empty:
        return
    group_columns = [*FAMILY_COLUMNS, entity_column]
    entity_summary = (
        frame.groupby(group_columns, dropna=False, observed=True)
        .agg(
            n_queries_total=("success_bool", "size"),
            n_positive=("success_bool", "sum"),
        )
        .reset_index()
    )
    for row in entity_summary.itertuples(index=False):
        entity_key = tuple(str(getattr(row, column)) for column in group_columns)
        if entity_key in entity_counts:
            entity_counts[entity_key][0] += int(row.n_queries_total)
            entity_counts[entity_key][1] += int(row.n_positive)
        else:
            entity_counts[entity_key] = [int(row.n_queries_total), int(row.n_positive)]


def update_family_counts(
    family_counts: dict[tuple[str, ...], list[int]],
    frame: pd.DataFrame,
) -> None:
    if frame.empty:
        return
    family_summary = (
        frame.groupby(FAMILY_COLUMNS, dropna=False, observed=True)
        .agg(
            n_queries_total=("success_bool", "size"),
            n_positive=("success_bool", "sum"),
        )
        .reset_index()
    )
    for row in family_summary.itertuples(index=False):
        family_key = tuple(str(getattr(row, column)) for column in FAMILY_COLUMNS)
        if family_key in family_counts:
            family_counts[family_key][0] += int(row.n_queries_total)
            family_counts[family_key][1] += int(row.n_positive)
        else:
            family_counts[family_key] = [int(row.n_queries_total), int(row.n_positive)]


def explode_target_rows(frame: pd.DataFrame) -> pd.DataFrame:
    expanded = frame.copy()
    expanded["target_token"] = expanded["target_token"].map(split_target_tokens)
    expanded = expanded.loc[expanded["target_token"].map(bool)].copy()
    if expanded.empty:
        return expanded
    expanded = expanded.explode("target_token", ignore_index=True)
    expanded["target_token"] = expanded["target_token"].astype(str)
    return expanded


def aggregate_counts(path: Path, entity_column: str) -> tuple[dict[tuple[str, ...], list[int]], dict[tuple[str, ...], list[int]]]:
    require_path(path, f"{entity_column} source parquet")
    parquet_file = pq.ParquetFile(path)
    entity_counts: dict[tuple[str, ...], list[int]] = {}
    family_counts: dict[tuple[str, ...], list[int]] = {}
    for batch in parquet_file.iter_batches(batch_size=250_000, columns=REQUIRED_COLUMNS):
        batch_frame = batch.to_pandas()
        batch_frame["success_bool"] = (
            batch_frame["mrr_corrected"].astype("float64").gt(0.0)
            & batch_frame["hit10_corrected"].astype("float64").gt(0.0)
        )
        for column in FAMILY_COLUMNS:
            batch_frame[column] = batch_frame[column].astype(str)
        update_family_counts(family_counts, batch_frame)
        if entity_column == "target_token":
            entity_frame = explode_target_rows(batch_frame)
        else:
            entity_frame = batch_frame.copy()
            entity_frame[entity_column] = entity_frame[entity_column].astype(str)
        update_entity_counts(entity_counts, entity_frame, entity_column)
    return entity_counts, family_counts


def build_summary_rows(
    entity_counts: dict[tuple[str, ...], list[int]],
    family_counts: dict[tuple[str, ...], list[int]],
    entity_column: str,
    source_path: Path,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for entity_key, counts in entity_counts.items():
        family_key = entity_key[: len(FAMILY_COLUMNS)]
        family_total, family_positive = family_counts[family_key]
        n_queries_total, n_positive = counts
        n_negative = n_queries_total - n_positive
        family_negative = family_total - family_positive
        n_background_positive = family_positive - n_positive
        n_background_negative = family_negative - n_negative

        raw_p = pd.NA
        odds_ratio = pd.NA
        if (n_background_positive + n_background_negative) > 0:
            odds_ratio_value, raw_p_value = fisher_exact(
                [
                    [n_positive, n_negative],
                    [n_background_positive, n_background_negative],
                ]
            )
            odds_ratio = float(odds_ratio_value)
            raw_p = float(raw_p_value)

        row: dict[str, object] = {column: entity_key[idx] for idx, column in enumerate([*FAMILY_COLUMNS, entity_column])}
        row["n_queries_total"] = n_queries_total
        row["n_positive"] = n_positive
        row["n_negative"] = n_negative
        row["n_background_positive"] = n_background_positive
        row["n_background_negative"] = n_background_negative
        row["success_rate"] = (n_positive / n_queries_total) if n_queries_total else pd.NA
        background_total = n_background_positive + n_background_negative
        row["background_success_rate"] = (n_background_positive / background_total) if background_total else pd.NA
        row["odds_ratio"] = odds_ratio
        row["raw_p"] = raw_p
        row["bh_q"] = pd.NA
        row["significant_bool"] = pd.NA
        row["success_definition"] = SUCCESS_DEFINITION
        row["source_path"] = str(source_path)
        rows.append(row)

    return pd.DataFrame(rows)


def apply_family_bh(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    out = frame.copy()
    out["bh_q"] = pd.Series(pd.NA, index=out.index, dtype="Float64")
    out["significant_bool"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
    valid_mask = out["raw_p"].notna()
    for _, family_idx in out.loc[valid_mask].groupby(FAMILY_COLUMNS, sort=False).groups.items():
        q_values = bh_adjust(out.loc[family_idx, "raw_p"].astype("float64"))
        out.loc[family_idx, "bh_q"] = q_values.astype("Float64")
        out.loc[family_idx, "significant_bool"] = q_values.le(0.05).astype("boolean")
    return out


def build_output(entity_column: str, input_paths: list[Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in input_paths:
        entity_counts, family_counts = aggregate_counts(path, entity_column)
        frame = build_summary_rows(entity_counts, family_counts, entity_column, path)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=[*FAMILY_COLUMNS, entity_column])
    out = pd.concat(frames, ignore_index=True)
    out = apply_family_bh(out)
    sort_columns = [*FAMILY_COLUMNS, entity_column]
    out = out.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)
    if entity_column == "cell_line":
        return out.reindex(columns=OUTPUT_COLUMNS_CELL_LINE)
    return out.reindex(columns=OUTPUT_COLUMNS_TARGET)


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    input_paths = [
        resolve_path(project_root, args.internal_per_query_path),
        resolve_path(project_root, args.cross_per_query_path),
    ]
    cell_line_output_path = resolve_path(project_root, args.cell_line_output_path)
    target_output_path = resolve_path(project_root, args.target_output_path)

    cell_line_output_path.parent.mkdir(parents=True, exist_ok=True)
    target_output_path.parent.mkdir(parents=True, exist_ok=True)

    cell_line_frame = build_output("cell_line", input_paths)
    target_frame = build_output("target_token", input_paths)

    # Mark strata where n_queries_total is too small for reliable interpretation.
    # This is a reporting/interpretation constraint flag; rows are never deleted.
    _underpowered_threshold = 20
    cell_line_frame["underpowered_strata_bool"] = cell_line_frame["n_queries_total"] < _underpowered_threshold
    target_frame["underpowered_strata_bool"] = target_frame["n_queries_total"] < _underpowered_threshold

    cell_line_frame.to_csv(cell_line_output_path, index=False)
    target_frame.to_csv(target_output_path, index=False)

    print(f"{cell_line_output_path}: {len(cell_line_frame)} rows")
    print(f"{target_output_path}: {len(target_frame)} rows")


if __name__ == "__main__":
    main()
