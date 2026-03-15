#!/usr/bin/env python3
# SCRIPT_HEADER_CONTRACT
# Script: scripts/s6_task2_result_synthesis_multisource.py
# Purpose: Synthesize corrected multisource Task2 S4/S5 outputs into
#   benchmark-ready summary tables without recomputing any Task2 metric.
# Inputs:
#   - audited S4 stage dir with:
#     - task2_group_concordance.csv
#     - task2_group_attrition.csv
#     - run_manifest.json
#     - audit_assertions.json
#     - manifest.json
#   - audited S5 stage dir with:
#     - task2_retrieval_summary.csv
#     - task2_retrieval_summary_long.csv
#     - task2_retrieval_attrition.csv
#     - task2_chance_identity_check.csv
#     - run_manifest.json
#     - audit_assertions.json
#     - manifest.json
# Outputs:
#   - task2_group_concordance_long.csv: runs/<run_id>/s6_task2_result_synthesis_multisource/
#   - task2_group_leaderboard.csv: runs/<run_id>/s6_task2_result_synthesis_multisource/
#   - task2_retrieval_leaderboard.csv: runs/<run_id>/s6_task2_result_synthesis_multisource/
#   - task2_benchmark_summary_long.csv: runs/<run_id>/s6_task2_result_synthesis_multisource/
#   - AVCP artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s6_task2_result_synthesis_multisource/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s6_task2_result_synthesis_multisource.py --run-id <run_id> --s4-stage-dir <path> --s5-stage-dir <path> --seed 619
# Failure Modes:
#   - Missing or unaudited upstream files -> exit non-zero
#   - Snapshot / schema / scope / rank-order / chance-tolerance drift -> exit non-zero
#   - Attempted edist ranking or benchmark summary contract violation -> exit non-zero
# Last Updated: 2026-03-11

"""
Inputs:
- audited corrected multisource S4/S5 stage outputs only

Outputs:
- task2_group_concordance_long.csv
- task2_group_leaderboard.csv
- task2_retrieval_leaderboard.csv
- task2_benchmark_summary_long.csv
- run_manifest.json
- audit_assertions.json
- manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- Dataset order = LINCS, scPerturb
- Direction order = C2G, G2C
- Representation order = Gene, Pathway, scgpt, geneformer, scbert, scfoundation, uce, state, tahoe-x1
- target_membership_source = delta_meta.target_tokens

Interpretation rules:
- S6 is synthesis only and does not recompute S4/S5 metrics
- edist_biascorr remains informational only and never receives leaderboard ranks
- dataset, cell_line, direction, and representation stratification must remain intact
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml

STAGE = "s6_task2_result_synthesis_multisource"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v2")
GLOBAL_SEED = 619
CHANCE_IDENTITY_TOL = 1e-12
MAX_COUNTEREXAMPLES = 5

EXPECTED_S4_STAGE = "s4_task2_group_concordance_multisource"
EXPECTED_S5_STAGE = "s5_task2_retrieval_multisource"
TARGET_MEMBERSHIP_SOURCE = "delta_meta.target_tokens"

DATASET_ORDER: Tuple[str, ...] = ("LINCS", "scPerturb")
DIRECTION_ORDER: Tuple[str, ...] = ("C2G", "G2C")
BENCHMARK_DIRECTION_ORDER: Tuple[str, ...] = ("NA",) + DIRECTION_ORDER
REPRESENTATION_ORDER: Tuple[str, ...] = (
    "Gene",
    "Pathway",
    "scgpt",
    "geneformer",
    "scbert",
    "scfoundation",
    "uce",
    "state",
    "tahoe-x1",
)
ALLOWED_SCOPE: Mapping[str, Tuple[str, ...]] = {
    "LINCS": ("Gene", "Pathway"),
    "scPerturb": REPRESENTATION_ORDER,
}
GROUP_METRIC_SPECS: Tuple[Tuple[str, str, str, str], ...] = (
    ("cosine_centroid", "mean_cosine_centroid", "cosine_valid_bool", "cosine_na_reason"),
    ("pcc_centroid", "mean_pcc_centroid", "pcc_valid_bool", "pcc_na_reason"),
    ("edist_biascorr", "mean_edist_biascorr", "edist_valid_bool", "edist_na_reason"),
)
GROUP_METRIC_ORDER: Tuple[str, ...] = tuple(spec[0] for spec in GROUP_METRIC_SPECS)
GROUP_SUMMARY_METRIC_ORDER: Tuple[str, ...] = tuple(spec[1] for spec in GROUP_METRIC_SPECS)
RETRIEVAL_METRICS: Tuple[str, ...] = (
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
)
ANALYSIS_FAMILY_ORDER: Tuple[str, ...] = ("group_concordance", "retrieval")

C2G_GALLERY_DEFINITION = "C2G_genetic_target_centroid_gallery_v1"
G2C_GALLERY_DEFINITION = "G2C_chemical_target_centroid_gallery_v1"
C2G_POS_DEFINITION = "C2G_tokens_intersect_genetic_target_gallery_v1"
G2C_POS_DEFINITION = "G2C_atomic_token_matches_chemical_target_gallery_v1"

S4_REQUIRED_FILES: Tuple[str, ...] = (
    "task2_group_concordance.csv",
    "task2_group_attrition.csv",
    "run_manifest.json",
    "audit_assertions.json",
    "manifest.json",
)
S5_REQUIRED_FILES: Tuple[str, ...] = (
    "task2_retrieval_summary.csv",
    "task2_retrieval_summary_long.csv",
    "task2_retrieval_attrition.csv",
    "task2_chance_identity_check.csv",
    "run_manifest.json",
    "audit_assertions.json",
    "manifest.json",
)

S4_CONCORDANCE_COLUMNS = [
    "dataset",
    "cell_line",
    "target_token",
    "group_id",
    "representation",
    "n_chem_instances_total",
    "n_gen_instances_total",
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
    "cosine_na_reason",
    "pcc_na_reason",
    "edist_na_reason",
]
S4_ATTRITION_COLUMNS = [
    "dataset",
    "cell_line",
    "target_token",
    "representation",
    "metric_name",
    "reason",
    "n_chem_before",
    "n_chem_after",
    "n_gen_before",
    "n_gen_after",
    "n_chem_removed",
    "n_gen_removed",
    "notes",
]
S5_SUMMARY_COLUMNS = [
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
]
S5_SUMMARY_LONG_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "metric_name",
    "metric_value",
    "n_valid",
    "N_gallery_mean",
    "N_gallery_max",
    "m_pos_mean",
    "m_pos_p50",
    "m_pos_p90",
]
S5_ATTRITION_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "reason",
    "n_query_rows_before",
    "n_query_rows_after",
    "n_query_rows_removed",
    "n_gallery_member_rows_before",
    "n_gallery_member_rows_after",
    "n_gallery_member_rows_removed",
    "n_gallery_items_before",
    "n_gallery_items_after",
    "n_gallery_items_removed",
    "notes",
]
S5_CHANCE_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "delta_mrr",
    "abs_delta_mrr",
    "delta_hit1",
    "abs_delta_hit1",
    "delta_hit5",
    "abs_delta_hit5",
    "delta_hit10",
    "abs_delta_hit10",
]

GROUP_LONG_COLUMNS = [
    "dataset",
    "cell_line",
    "target_token",
    "group_id",
    "representation",
    "metric_name",
    "metric_value",
    "metric_valid_bool",
    "metric_na_reason",
    "n_chem_instances_total",
    "n_gen_instances_total",
    "n_chem_instances_used",
    "n_gen_instances_used",
    "n_chem_sub",
    "n_gen_sub",
]
GROUP_LEADERBOARD_COLUMNS = [
    "dataset",
    "cell_line",
    "representation",
    "n_targets_total",
    "n_targets_cosine_valid",
    "n_targets_pcc_valid",
    "n_targets_edist_valid",
    "mean_cosine_centroid",
    "mean_pcc_centroid",
    "mean_edist_biascorr",
    "rank_by_mean_cosine_centroid",
    "rank_by_mean_pcc_centroid",
    "mean_edist_biascorr_leaderboard_eligible_bool",
    "mean_edist_biascorr_cross_representation_comparable_bool",
    "n_invalid_rows_unique",
    "n_invalid_drug_rows_unique",
    "n_attrition_target_rows",
    "n_attrition_chem_rows_removed_membership",
    "n_attrition_gen_rows_removed_membership",
]
RETRIEVAL_LEADERBOARD_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "rank_by_mean_mrr_corrected",
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
BENCHMARK_SUMMARY_LONG_COLUMNS = [
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
    "caution_codes",
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
]


@dataclass(frozen=True)
class StageBundle:
    stage_dir: Path
    run_manifest: Dict[str, Any]
    audit_assertions: Dict[str, Any]
    manifest: Dict[str, Any]
    tables: Dict[str, pd.DataFrame]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S6 corrected Task2 multisource result synthesis")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--s4-stage-dir", type=Path, required=True)
    parser.add_argument("--s5-stage-dir", type=Path, required=True)
    return parser.parse_args()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def resolve_config_path(project_root: Path, raw_path: str) -> Path:
    path = Path(str(raw_path))
    if path.is_absolute():
        return path.resolve()
    return (project_root / path).resolve()


def safe_git_head(project_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=project_root,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except Exception:
        return "unavailable"


def compute_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=True)
        handle.write("\n")


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def ensure_required_columns(frame: pd.DataFrame, required: Sequence[str], name: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns for {name}: {missing}")


def assert_unique_key(frame: pd.DataFrame, columns: Sequence[str], name: str) -> None:
    if frame.empty:
        return
    duplicates = frame.loc[frame.duplicated(list(columns), keep=False), list(columns)]
    if not duplicates.empty:
        examples = duplicates.drop_duplicates().head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"{name} is not unique on key={list(columns)}: {examples}")


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def snapshot_matches_expected(raw_path: object) -> bool:
    value = Path(str(raw_path)).as_posix().rstrip("/")
    expected = EXPECTED_TASK2_SNAPSHOT.as_posix()
    return value == expected or value.endswith(f"/{expected}")


def numeric_columns_match(left: pd.Series, right: pd.Series, tol: float = 1e-12) -> pd.Series:
    left_num = pd.to_numeric(left, errors="coerce")
    right_num = pd.to_numeric(right, errors="coerce")
    both_na = left_num.isna() & right_num.isna()
    close = (left_num - right_num).abs().le(tol)
    close = close.fillna(False)
    return both_na | close


def rank_desc(series: pd.Series) -> pd.Series:
    out = pd.Series(pd.NA, index=series.index, dtype="Int64")
    valid = series.notna()
    if valid.any():
        out.loc[valid] = series.loc[valid].rank(method="min", ascending=False).astype("Int64")
    return out


def require_stage_files(stage_dir: Path, required_files: Sequence[str]) -> Dict[str, Path]:
    if not stage_dir.exists():
        raise FileNotFoundError(f"Stage directory does not exist: {stage_dir}")
    if not stage_dir.is_dir():
        raise NotADirectoryError(f"Stage path is not a directory: {stage_dir}")

    out: Dict[str, Path] = {}
    for filename in required_files:
        path = (stage_dir / filename).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Missing required upstream file: {path}")
        out[filename] = path
    return out


def load_stage_bundle(stage_dir: Path, required_files: Sequence[str], csv_files: Sequence[str]) -> StageBundle:
    paths = require_stage_files(stage_dir, required_files)
    tables = {name: pd.read_csv(paths[name]) for name in csv_files}
    return StageBundle(
        stage_dir=stage_dir.resolve(),
        run_manifest=load_json(paths["run_manifest.json"]),
        audit_assertions=load_json(paths["audit_assertions.json"]),
        manifest=load_json(paths["manifest.json"]),
        tables=tables,
    )


def validate_manifest_file_list(bundle: StageBundle, required_files: Sequence[str], stage_name: str) -> None:
    manifest_stage = bundle.manifest.get("stage")
    if manifest_stage != stage_name:
        raise ValueError(f"{stage_name} manifest stage mismatch: {manifest_stage}")
    files = bundle.manifest.get("files")
    if not isinstance(files, list):
        raise ValueError(f"{stage_name} manifest missing files list")
    available = {str(item.get("relative_path", "")) for item in files if isinstance(item, Mapping)}
    missing = [filename for filename in required_files if filename != "manifest.json" and filename not in available]
    if missing:
        raise ValueError(f"{stage_name} manifest.json missing required files: {missing}")


def assert_all_upstream_assertions_pass(bundle: StageBundle, stage_name: str) -> None:
    raw_assertions = bundle.audit_assertions.get("assertions")
    if not isinstance(raw_assertions, list):
        raise ValueError(f"{stage_name} audit_assertions.json missing 'assertions' list")
    failing = [item for item in raw_assertions if not bool(item.get("pass", False))]
    if failing:
        names = [str(item.get("name", "<unnamed>")) for item in failing[:MAX_COUNTEREXAMPLES]]
        raise ValueError(f"{stage_name} has failing upstream assertions: {names}")


def build_representation_details_map(run_manifest: Mapping[str, Any]) -> Dict[Tuple[str, str], Dict[str, Any]]:
    summary = run_manifest.get("summary")
    if not isinstance(summary, Mapping):
        raise ValueError("run_manifest missing summary block")
    details = summary.get("representation_details")
    if not isinstance(details, list):
        raise ValueError("run_manifest missing representation_details")
    out: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for detail in details:
        if not isinstance(detail, Mapping):
            continue
        dataset = str(detail.get("dataset", ""))
        representation = str(detail.get("representation", ""))
        if not dataset or not representation:
            continue
        key = (dataset, representation)
        if key in out:
            raise ValueError(f"Duplicate representation_details key: {key}")
        out[key] = dict(detail)
    return out


def extract_non_blocking_notes(run_manifest: Mapping[str, Any]) -> List[str]:
    notes: List[str] = []
    config = run_manifest.get("config")
    if isinstance(config, Mapping):
        for key, value in config.items():
            if "note" in str(key).lower() and isinstance(value, str) and value.strip():
                notes.append(f"{key}: {value.strip()}")
    return notes


def first_seen_values(series: Iterable[str]) -> List[str]:
    seen: List[str] = []
    for value in series:
        if value not in seen:
            seen.append(value)
    return seen


def extract_cell_line_order_by_dataset(frame: pd.DataFrame) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for dataset in DATASET_ORDER:
        subset = frame.loc[frame["dataset"].eq(dataset), "cell_line"]
        out[dataset] = first_seen_values(subset.tolist())
    return out


def build_cell_line_order_map(s4_concordance: pd.DataFrame, s5_summary: pd.DataFrame) -> Dict[str, Dict[str, int]]:
    s4_order = extract_cell_line_order_by_dataset(s4_concordance[["dataset", "cell_line"]].drop_duplicates())
    s5_order = extract_cell_line_order_by_dataset(s5_summary[["dataset", "cell_line"]].drop_duplicates())
    if s4_order != s5_order:
        raise ValueError(f"S4/S5 cell_line order drift detected: s4={s4_order} s5={s5_order}")
    return {
        dataset: {cell_line: idx for idx, cell_line in enumerate(cell_lines)}
        for dataset, cell_lines in s4_order.items()
    }


def extract_target_order_by_slice(concordance: pd.DataFrame) -> Dict[Tuple[str, str], List[str]]:
    target_frame = concordance[["dataset", "cell_line", "target_token"]].drop_duplicates().reset_index(drop=True)
    out: Dict[Tuple[str, str], List[str]] = {}
    for row in target_frame.itertuples(index=False):
        key = (str(row.dataset), str(row.cell_line))
        out.setdefault(key, []).append(str(row.target_token))
    return out


def validate_dataset_and_representation_scope(frame: pd.DataFrame, name: str, include_direction: bool = False) -> None:
    actual_datasets = first_seen_values(frame["dataset"].tolist())
    if actual_datasets != list(DATASET_ORDER):
        raise ValueError(f"{name} dataset order drift: {actual_datasets}")
    if include_direction:
        actual_directions = first_seen_values(frame["direction"].tolist())
        if actual_directions != list(DIRECTION_ORDER):
            raise ValueError(f"{name} direction order drift: {actual_directions}")
    for dataset in DATASET_ORDER:
        reps = first_seen_values(frame.loc[frame["dataset"].eq(dataset), "representation"].tolist())
        expected_reps = [rep for rep in REPRESENTATION_ORDER if rep in ALLOWED_SCOPE[dataset]]
        if reps != expected_reps:
            raise ValueError(f"{name} representation scope drift for dataset={dataset}: {reps}")
    scperturb_cells = first_seen_values(frame.loc[frame["dataset"].eq("scPerturb"), "cell_line"].tolist())
    if scperturb_cells != ["K562"]:
        raise ValueError(f"{name} scPerturb cell_line drift: {scperturb_cells}")


def validate_group_id_contract(concordance: pd.DataFrame) -> None:
    expected = (
        concordance["dataset"].astype(str)
        + "||"
        + concordance["cell_line"].astype(str)
        + "||"
        + concordance["target_token"].astype(str)
    )
    bad = concordance.loc[expected.ne(concordance["group_id"].astype(str))]
    if not bad.empty:
        raise ValueError(
            "S4 group_id contract drift: "
            f"{bad[['dataset', 'cell_line', 'target_token', 'group_id']].head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )


def validate_representation_detail_scope(
    details: Mapping[Tuple[str, str], Mapping[str, Any]],
    stage_name: str,
) -> None:
    expected_keys = {
        (dataset, representation)
        for dataset in DATASET_ORDER
        for representation in ALLOWED_SCOPE[dataset]
    }
    actual_keys = set(details.keys())
    if actual_keys != expected_keys:
        raise ValueError(
            f"{stage_name} representation_details keys drift: missing={sorted(expected_keys - actual_keys)} "
            f"extra={sorted(actual_keys - expected_keys)}"
        )


def validate_s5_long_matches_summary(summary: pd.DataFrame, summary_long: pd.DataFrame) -> None:
    key_cols = ["dataset", "cell_line", "direction", "representation"]
    metric_names = summary_long["metric_name"].drop_duplicates().tolist()
    if metric_names != list(RETRIEVAL_METRICS):
        raise ValueError(f"S5 summary_long metric drift: {metric_names}")

    metric_frame = (
        summary_long[key_cols + ["metric_name", "metric_value"]]
        .pivot(index=key_cols, columns="metric_name", values="metric_value")
        .reset_index()
    )
    metric_frame.columns.name = None
    denom_cols = ["n_valid", "N_gallery_mean", "N_gallery_max", "m_pos_mean", "m_pos_p50", "m_pos_p90"]
    denom_frame = summary_long[key_cols + denom_cols].drop_duplicates().reset_index(drop=True)
    assert_unique_key(metric_frame, key_cols, "S5 task2_retrieval_summary_long metric pivot")
    assert_unique_key(denom_frame, key_cols, "S5 task2_retrieval_summary_long denominators")

    merged = summary[key_cols + denom_cols + list(RETRIEVAL_METRICS)].merge(
        metric_frame,
        on=key_cols,
        how="outer",
        sort=False,
        suffixes=("_summary", "_long"),
        indicator=True,
    )
    if not (merged["_merge"] == "both").all():
        raise ValueError(
            "S5 summary_long key grid drift: "
            f"{merged.loc[merged['_merge'].ne('both'), key_cols + ['_merge']].head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )
    merged = merged.drop(columns=["_merge"]).merge(denom_frame, on=key_cols, how="inner", sort=False, suffixes=("", "_long"))

    for column in denom_cols:
        check = numeric_columns_match(merged[column], merged[f"{column}_long"])
        if not bool(check.all()):
            bad = merged.loc[~check, key_cols + [column, f"{column}_long"]]
            raise ValueError(
                f"S5 summary_long denominator drift for {column}: "
                f"{bad.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
            )

    for metric_name in RETRIEVAL_METRICS:
        check = numeric_columns_match(merged[f"{metric_name}_summary"], merged[f"{metric_name}_long"])
        if not bool(check.all()):
            bad = merged.loc[~check, key_cols + [f"{metric_name}_summary", f"{metric_name}_long"]]
            raise ValueError(
                f"S5 summary_long metric drift for {metric_name}: "
                f"{bad.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
            )


def validate_s4_bundle(bundle: StageBundle) -> Dict[Tuple[str, str], Dict[str, Any]]:
    if bundle.run_manifest.get("stage") != EXPECTED_S4_STAGE:
        raise ValueError(f"S4 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S4_REQUIRED_FILES, EXPECTED_S4_STAGE)
    assert_all_upstream_assertions_pass(bundle, EXPECTED_S4_STAGE)

    config = bundle.run_manifest.get("config")
    summary = bundle.run_manifest.get("summary")
    if not isinstance(config, Mapping):
        raise ValueError("S4 run_manifest missing config block")
    if not isinstance(summary, Mapping):
        raise ValueError("S4 run_manifest missing summary block")
    if int(config.get("seed", -1)) != GLOBAL_SEED:
        raise ValueError(f"S4 seed drift detected: {config.get('seed')}")
    if list(config.get("dataset_order", [])) != list(DATASET_ORDER):
        raise ValueError("S4 dataset_order drift detected")
    if list(config.get("representation_order", [])) != list(REPRESENTATION_ORDER):
        raise ValueError("S4 representation_order drift detected")
    if config.get("target_membership_source") != TARGET_MEMBERSHIP_SOURCE:
        raise ValueError("S4 target_membership_source drift detected")
    if not snapshot_matches_expected(config.get("task2_snapshot", "")):
        raise ValueError(f"S4 task2_snapshot drift detected: {config.get('task2_snapshot')}")

    concordance = bundle.tables["task2_group_concordance.csv"]
    attrition = bundle.tables["task2_group_attrition.csv"]
    ensure_required_columns(concordance, S4_CONCORDANCE_COLUMNS, "S4 task2_group_concordance.csv")
    ensure_required_columns(attrition, S4_ATTRITION_COLUMNS, "S4 task2_group_attrition.csv")
    assert_unique_key(concordance, ["dataset", "cell_line", "target_token", "representation"], "S4 task2_group_concordance.csv")
    assert_unique_key(
        attrition,
        ["dataset", "cell_line", "target_token", "representation", "metric_name", "reason"],
        "S4 task2_group_attrition.csv",
    )
    validate_dataset_and_representation_scope(concordance, "S4 task2_group_concordance.csv")
    validate_group_id_contract(concordance)

    rep_details = build_representation_details_map(bundle.run_manifest)
    validate_representation_detail_scope(rep_details, EXPECTED_S4_STAGE)

    if int(summary.get("n_concordance_rows", -1)) != int(len(concordance)):
        raise ValueError("S4 summary n_concordance_rows drift detected")
    if int(summary.get("n_attrition_rows", -1)) != int(len(attrition)):
        raise ValueError("S4 summary n_attrition_rows drift detected")

    actual_grid = (
        concordance.groupby(["dataset", "representation"], sort=False).size().reset_index(name="n_rows")
    )
    expected_grid = pd.DataFrame(summary.get("supported_grid_breakdown", []))
    if expected_grid.empty or not actual_grid.equals(expected_grid[["dataset", "representation", "n_rows"]]):
        raise ValueError("S4 supported_grid_breakdown drift detected")

    actual_mech = (
        concordance[["dataset", "cell_line", "target_token"]]
        .drop_duplicates()
        .groupby("dataset", sort=False)
        .size()
        .to_dict()
    )
    expected_mech = {str(key): int(value) for key, value in summary.get("eligible_mech_keys_by_dataset", {}).items()}
    if actual_mech != expected_mech:
        raise ValueError(f"S4 eligible_mech_keys_by_dataset drift: actual={actual_mech} expected={expected_mech}")

    return rep_details


def validate_s5_bundle(bundle: StageBundle) -> Dict[Tuple[str, str], Dict[str, Any]]:
    if bundle.run_manifest.get("stage") != EXPECTED_S5_STAGE:
        raise ValueError(f"S5 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S5_REQUIRED_FILES, EXPECTED_S5_STAGE)
    assert_all_upstream_assertions_pass(bundle, EXPECTED_S5_STAGE)

    config = bundle.run_manifest.get("config")
    summary_meta = bundle.run_manifest.get("summary")
    if not isinstance(config, Mapping):
        raise ValueError("S5 run_manifest missing config block")
    if not isinstance(summary_meta, Mapping):
        raise ValueError("S5 run_manifest missing summary block")
    if int(config.get("seed", -1)) != GLOBAL_SEED:
        raise ValueError(f"S5 seed drift detected: {config.get('seed')}")
    if list(config.get("dataset_order", [])) != list(DATASET_ORDER):
        raise ValueError("S5 dataset_order drift detected")
    if list(config.get("direction_order", [])) != list(DIRECTION_ORDER):
        raise ValueError("S5 direction_order drift detected")
    if list(config.get("representation_order", [])) != list(REPRESENTATION_ORDER):
        raise ValueError("S5 representation_order drift detected")
    if config.get("target_membership_source") != TARGET_MEMBERSHIP_SOURCE:
        raise ValueError("S5 target_membership_source drift detected")
    if not snapshot_matches_expected(config.get("task2_snapshot", "")):
        raise ValueError(f"S5 task2_snapshot drift detected: {config.get('task2_snapshot')}")

    expected_gallery_ids = {"C2G": C2G_GALLERY_DEFINITION, "G2C": G2C_GALLERY_DEFINITION}
    expected_pos_ids = {"C2G": C2G_POS_DEFINITION, "G2C": G2C_POS_DEFINITION}
    if dict(config.get("gallery_definition_ids", {})) != expected_gallery_ids:
        raise ValueError("S5 gallery_definition_ids drift detected")
    if dict(config.get("pos_definition_ids", {})) != expected_pos_ids:
        raise ValueError("S5 pos_definition_ids drift detected")
    if float(config.get("chance_identity_tol", -1.0)) != CHANCE_IDENTITY_TOL:
        raise ValueError("S5 chance_identity_tol drift detected")

    summary = bundle.tables["task2_retrieval_summary.csv"]
    summary_long = bundle.tables["task2_retrieval_summary_long.csv"]
    attrition = bundle.tables["task2_retrieval_attrition.csv"]
    chance = bundle.tables["task2_chance_identity_check.csv"]
    ensure_required_columns(summary, S5_SUMMARY_COLUMNS, "S5 task2_retrieval_summary.csv")
    ensure_required_columns(summary_long, S5_SUMMARY_LONG_COLUMNS, "S5 task2_retrieval_summary_long.csv")
    ensure_required_columns(attrition, S5_ATTRITION_COLUMNS, "S5 task2_retrieval_attrition.csv")
    ensure_required_columns(chance, S5_CHANCE_COLUMNS, "S5 task2_chance_identity_check.csv")
    assert_unique_key(summary, ["dataset", "cell_line", "direction", "representation"], "S5 task2_retrieval_summary.csv")
    assert_unique_key(
        summary_long,
        ["dataset", "cell_line", "direction", "representation", "metric_name"],
        "S5 task2_retrieval_summary_long.csv",
    )
    assert_unique_key(
        attrition,
        ["dataset", "cell_line", "direction", "representation", "reason"],
        "S5 task2_retrieval_attrition.csv",
    )
    assert_unique_key(
        chance,
        ["dataset", "cell_line", "direction", "representation"],
        "S5 task2_chance_identity_check.csv",
    )
    validate_dataset_and_representation_scope(summary, "S5 task2_retrieval_summary.csv", include_direction=True)
    validate_s5_long_matches_summary(summary, summary_long)

    if int(summary_meta.get("n_summary_rows", -1)) != int(len(summary)):
        raise ValueError("S5 summary n_summary_rows drift detected")
    if int(summary_meta.get("n_summary_long_rows", -1)) != int(len(summary_long)):
        raise ValueError("S5 summary n_summary_long_rows drift detected")
    if int(summary_meta.get("n_attrition_rows", -1)) != int(len(attrition)):
        raise ValueError("S5 summary n_attrition_rows drift detected")
    if int(summary_meta.get("n_chance_rows", -1)) != int(len(chance)):
        raise ValueError("S5 summary n_chance_rows drift detected")
    if int(summary_meta.get("expected_summary_rows", -1)) != int(len(summary)):
        raise ValueError("S5 expected_summary_rows drift detected")
    if int(summary_meta.get("expected_summary_long_rows", -1)) != int(len(summary_long)):
        raise ValueError("S5 expected_summary_long_rows drift detected")

    actual_by_dataset = {str(key): int(value) for key, value in summary.groupby("dataset", sort=False).size().to_dict().items()}
    expected_by_dataset = {str(key): int(value) for key, value in summary_meta.get("actual_by_dataset", {}).items()}
    if actual_by_dataset != expected_by_dataset:
        raise ValueError(f"S5 actual_by_dataset drift: actual={actual_by_dataset} expected={expected_by_dataset}")

    tol_cols = ["abs_delta_mrr", "abs_delta_hit1", "abs_delta_hit5", "abs_delta_hit10"]
    for column in tol_cols:
        if (pd.to_numeric(chance[column], errors="coerce").abs() > CHANCE_IDENTITY_TOL).any():
            raise ValueError(f"S5 chance identity tolerance exceeded for {column}")

    rep_details = build_representation_details_map(bundle.run_manifest)
    validate_representation_detail_scope(rep_details, EXPECTED_S5_STAGE)
    return rep_details


def validate_cross_stage_contracts(
    s4_bundle: StageBundle,
    s5_bundle: StageBundle,
    s4_rep_details: Mapping[Tuple[str, str], Mapping[str, Any]],
    s5_rep_details: Mapping[Tuple[str, str], Mapping[str, Any]],
) -> Dict[str, Dict[str, int]]:
    s4_concordance = s4_bundle.tables["task2_group_concordance.csv"]
    s5_summary = s5_bundle.tables["task2_retrieval_summary.csv"]

    s4_slice_grid = s4_concordance[["dataset", "cell_line", "representation"]].drop_duplicates().reset_index(drop=True)
    s5_slice_grid = s5_summary[["dataset", "cell_line", "representation"]].drop_duplicates().reset_index(drop=True)
    assert_unique_key(s4_slice_grid, ["dataset", "cell_line", "representation"], "S4 slice grid")
    assert_unique_key(s5_slice_grid, ["dataset", "cell_line", "representation"], "S5 slice grid")
    if not s4_slice_grid.equals(s5_slice_grid):
        raise ValueError(
            "S4/S5 slice grid drift detected: "
            f"s4_head={s4_slice_grid.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')} "
            f"s5_head={s5_slice_grid.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )

    if set(s4_rep_details.keys()) != set(s5_rep_details.keys()):
        raise ValueError("S4/S5 representation_details key drift detected")
    for key in sorted(s4_rep_details.keys()):
        s4_detail = s4_rep_details[key]
        s5_detail = s5_rep_details[key]
        for field_name in ("family", "n_rows", "dim", "n_valid", "n_invalid"):
            if s4_detail.get(field_name) != s5_detail.get(field_name):
                raise ValueError(f"S4/S5 representation_details drift for {key} field={field_name}")
        if dict(s4_detail.get("invalid_by_side", {})) != dict(s5_detail.get("invalid_by_side", {})):
            raise ValueError(f"S4/S5 invalid_by_side drift for {key}")

    return build_cell_line_order_map(s4_concordance, s5_summary)


def build_group_concordance_long(concordance: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for row in concordance.to_dict(orient="records"):
        for metric_name, _, valid_col, na_reason_col in GROUP_METRIC_SPECS:
            rows.append(
                {
                    "dataset": row["dataset"],
                    "cell_line": row["cell_line"],
                    "target_token": row["target_token"],
                    "group_id": row["group_id"],
                    "representation": row["representation"],
                    "metric_name": metric_name,
                    "metric_value": row[metric_name],
                    "metric_valid_bool": bool(row[valid_col]),
                    "metric_na_reason": "" if pd.isna(row[na_reason_col]) else str(row[na_reason_col]),
                    "n_chem_instances_total": row["n_chem_instances_total"],
                    "n_gen_instances_total": row["n_gen_instances_total"],
                    "n_chem_instances_used": row["n_chem_instances_used"],
                    "n_gen_instances_used": row["n_gen_instances_used"],
                    "n_chem_sub": row["n_chem_sub"],
                    "n_gen_sub": row["n_gen_sub"],
                }
            )
    out = pd.DataFrame(rows)
    return out[GROUP_LONG_COLUMNS].reset_index(drop=True)


def build_group_leaderboard(
    concordance: pd.DataFrame,
    attrition: pd.DataFrame,
    rep_details: Mapping[Tuple[str, str], Mapping[str, Any]],
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    slice_grid = concordance[["dataset", "cell_line", "representation"]].drop_duplicates().reset_index(drop=True)
    for row in slice_grid.itertuples(index=False):
        dataset = str(row.dataset)
        cell_line = str(row.cell_line)
        representation = str(row.representation)
        subset = concordance.loc[
            concordance["dataset"].eq(dataset)
            & concordance["cell_line"].eq(cell_line)
            & concordance["representation"].eq(representation)
        ].reset_index(drop=True)
        if subset.empty:
            raise ValueError(f"Missing S4 rows for slice={(dataset, cell_line, representation)}")
        attr = attrition.loc[
            attrition["dataset"].eq(dataset)
            & attrition["cell_line"].eq(cell_line)
            & attrition["representation"].eq(representation)
        ].reset_index(drop=True)
        detail = dict(rep_details[(dataset, representation)])
        invalid_by_side = dict(detail.get("invalid_by_side", {}))
        rows.append(
            {
                "dataset": dataset,
                "cell_line": cell_line,
                "representation": representation,
                "n_targets_total": int(len(subset)),
                "n_targets_cosine_valid": int(subset["cosine_valid_bool"].astype(bool).sum()),
                "n_targets_pcc_valid": int(subset["pcc_valid_bool"].astype(bool).sum()),
                "n_targets_edist_valid": int(subset["edist_valid_bool"].astype(bool).sum()),
                "mean_cosine_centroid": subset.loc[subset["cosine_valid_bool"].astype(bool), "cosine_centroid"].mean(),
                "mean_pcc_centroid": subset.loc[subset["pcc_valid_bool"].astype(bool), "pcc_centroid"].mean(),
                "mean_edist_biascorr": subset.loc[subset["edist_valid_bool"].astype(bool), "edist_biascorr"].mean(),
                "mean_edist_biascorr_leaderboard_eligible_bool": False,
                "mean_edist_biascorr_cross_representation_comparable_bool": False,
                "n_invalid_rows_unique": int(detail.get("n_invalid", 0)),
                "n_invalid_drug_rows_unique": int(invalid_by_side.get("DRUG", 0)),
                "n_attrition_target_rows": int(len(attr)),
                "n_attrition_chem_rows_removed_membership": int(pd.to_numeric(attr["n_chem_removed"], errors="coerce").fillna(0).sum()) if not attr.empty else 0,
                "n_attrition_gen_rows_removed_membership": int(pd.to_numeric(attr["n_gen_removed"], errors="coerce").fillna(0).sum()) if not attr.empty else 0,
            }
        )

    out = pd.DataFrame(rows)
    out["rank_by_mean_cosine_centroid"] = (
        out.groupby(["dataset", "cell_line"], sort=False)["mean_cosine_centroid"].transform(rank_desc).astype("Int64")
    )
    out["rank_by_mean_pcc_centroid"] = (
        out.groupby(["dataset", "cell_line"], sort=False)["mean_pcc_centroid"].transform(rank_desc).astype("Int64")
    )

    dataset_order = {dataset: idx for idx, dataset in enumerate(DATASET_ORDER)}
    rep_order = {representation: idx for idx, representation in enumerate(REPRESENTATION_ORDER)}
    out["_dataset_order"] = out["dataset"].map(dataset_order)
    out["_cell_line_order"] = [
        cell_line_order_map[dataset][cell_line]
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_rank_order"] = out["rank_by_mean_cosine_centroid"]
    out["_representation_order"] = out["representation"].map(rep_order)
    out = out.sort_values(
        ["_dataset_order", "_cell_line_order", "_rank_order", "_representation_order"],
        kind="stable",
        na_position="last",
    )
    return out.drop(columns=["_dataset_order", "_cell_line_order", "_rank_order", "_representation_order"]).reset_index(
        drop=True
    )[GROUP_LEADERBOARD_COLUMNS]


def build_retrieval_leaderboard(
    summary: pd.DataFrame,
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    out = summary[
        [
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
    ].copy()
    out["rank_by_mean_mrr_corrected"] = (
        out.groupby(["dataset", "cell_line", "direction"], sort=False)["mean_mrr_corrected"].transform(rank_desc).astype("Int64")
    )

    dataset_order = {dataset: idx for idx, dataset in enumerate(DATASET_ORDER)}
    direction_order = {direction: idx for idx, direction in enumerate(DIRECTION_ORDER)}
    rep_order = {representation: idx for idx, representation in enumerate(REPRESENTATION_ORDER)}
    out["_dataset_order"] = out["dataset"].map(dataset_order)
    out["_cell_line_order"] = [
        cell_line_order_map[dataset][cell_line]
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_direction_order"] = out["direction"].map(direction_order)
    out["_rank_order"] = out["rank_by_mean_mrr_corrected"]
    out["_representation_order"] = out["representation"].map(rep_order)
    out = out.sort_values(
        ["_dataset_order", "_cell_line_order", "_direction_order", "_rank_order", "_representation_order"],
        kind="stable",
        na_position="last",
    )
    return out.drop(
        columns=["_dataset_order", "_cell_line_order", "_direction_order", "_rank_order", "_representation_order"]
    ).reset_index(drop=True)[RETRIEVAL_LEADERBOARD_COLUMNS]


def join_caution_codes(*codes: Optional[str]) -> str:
    ordered = [str(code) for code in codes if code]
    return ";".join(ordered)


def build_benchmark_summary_long(
    group_leaderboard: pd.DataFrame,
    retrieval_leaderboard: pd.DataFrame,
    s4_rep_details: Mapping[Tuple[str, str], Mapping[str, Any]],
    s5_rep_details: Mapping[Tuple[str, str], Mapping[str, Any]],
    s4_run_id: str,
    s5_run_id: str,
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    for row in group_leaderboard.to_dict(orient="records"):
        dataset = str(row["dataset"])
        cell_line = str(row["cell_line"])
        representation = str(row["representation"])
        base_codes: List[str] = []
        if representation == "uce" and int(row["n_invalid_rows_unique"]) > 0:
            base_codes.append("UCE_ATTRITION_PRESENT")

        group_metric_valid_counts = {
            "mean_cosine_centroid": int(row["n_targets_cosine_valid"]),
            "mean_pcc_centroid": int(row["n_targets_pcc_valid"]),
            "mean_edist_biascorr": int(row["n_targets_edist_valid"]),
        }
        metric_rows = [
            ("mean_cosine_centroid", row["mean_cosine_centroid"], row["rank_by_mean_cosine_centroid"], "mean_cosine_centroid", True, True),
            ("mean_pcc_centroid", row["mean_pcc_centroid"], row["rank_by_mean_pcc_centroid"], "mean_pcc_centroid", True, True),
            ("mean_edist_biascorr", row["mean_edist_biascorr"], pd.NA, "", False, False),
        ]
        for metric_name, metric_value, rank_value, rank_basis, eligible, comparable in metric_rows:
            metric_codes = list(base_codes)
            if metric_name == "mean_edist_biascorr":
                metric_codes.append("GROUP_EDIST_NOT_CROSS_REP")
                if int(group_metric_valid_counts[metric_name]) < int(row["n_targets_total"]):
                    metric_codes.append("GROUP_EDIST_PARTIAL_TARGET_VALIDITY")
            rows.append(
                {
                    "analysis_family": "group_concordance",
                    "dataset": dataset,
                    "cell_line": cell_line,
                    "direction": "NA",
                    "representation": representation,
                    "metric_name": metric_name,
                    "metric_value": metric_value,
                    "rank_value": rank_value,
                    "rank_basis_metric_name": rank_basis,
                    "leaderboard_eligible_bool": eligible,
                    "cross_representation_comparable_bool": comparable,
                    "source_table": "task2_group_leaderboard.csv",
                    "source_stage": EXPECTED_S4_STAGE,
                    "source_run_id": s4_run_id,
                    "caution_codes": join_caution_codes(*metric_codes),
                    "n_targets_total": row["n_targets_total"],
                    "n_targets_metric_valid": group_metric_valid_counts[metric_name],
                    "n_total": pd.NA,
                    "n_valid": pd.NA,
                    "n_excluded_missing_metric_or_mpos0": pd.NA,
                    "N_gallery_mean": pd.NA,
                    "N_gallery_max": pd.NA,
                    "m_pos_mean": pd.NA,
                    "m_pos_p50": pd.NA,
                    "m_pos_p90": pd.NA,
                }
            )

    invalid_retrieval_keys = {
        key
        for key, detail in s5_rep_details.items()
        if int(detail.get("n_invalid", 0)) > 0
    }
    for row in retrieval_leaderboard.to_dict(orient="records"):
        dataset = str(row["dataset"])
        representation = str(row["representation"])
        base_codes = ["RETRIEVAL_DIRECTION_SPECIFIC"]
        if (dataset, representation) in invalid_retrieval_keys and representation == "uce":
            base_codes.append("UCE_ATTRITION_PRESENT")
        for metric_name in RETRIEVAL_METRICS:
            rows.append(
                {
                    "analysis_family": "retrieval",
                    "dataset": row["dataset"],
                    "cell_line": row["cell_line"],
                    "direction": row["direction"],
                    "representation": representation,
                    "metric_name": metric_name,
                    "metric_value": row[metric_name],
                    "rank_value": row["rank_by_mean_mrr_corrected"],
                    "rank_basis_metric_name": "mean_mrr_corrected",
                    "leaderboard_eligible_bool": True,
                    "cross_representation_comparable_bool": True,
                    "source_table": "task2_retrieval_leaderboard.csv",
                    "source_stage": EXPECTED_S5_STAGE,
                    "source_run_id": s5_run_id,
                    "caution_codes": join_caution_codes(*base_codes),
                    "n_targets_total": pd.NA,
                    "n_targets_metric_valid": pd.NA,
                    "n_total": row["n_total"],
                    "n_valid": row["n_valid"],
                    "n_excluded_missing_metric_or_mpos0": row["n_excluded_missing_metric_or_mpos0"],
                    "N_gallery_mean": row["N_gallery_mean"],
                    "N_gallery_max": row["N_gallery_max"],
                    "m_pos_mean": row["m_pos_mean"],
                    "m_pos_p50": row["m_pos_p50"],
                    "m_pos_p90": row["m_pos_p90"],
                }
            )

    out = pd.DataFrame(rows)
    dataset_order = {dataset: idx for idx, dataset in enumerate(DATASET_ORDER)}
    family_order = {family: idx for idx, family in enumerate(ANALYSIS_FAMILY_ORDER)}
    direction_order = {direction: idx for idx, direction in enumerate(BENCHMARK_DIRECTION_ORDER)}
    rep_order = {representation: idx for idx, representation in enumerate(REPRESENTATION_ORDER)}
    metric_order = {name: idx for idx, name in enumerate(GROUP_SUMMARY_METRIC_ORDER + RETRIEVAL_METRICS)}

    out["_dataset_order"] = out["dataset"].map(dataset_order)
    out["_cell_line_order"] = [
        cell_line_order_map[dataset][cell_line]
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_family_order"] = out["analysis_family"].map(family_order)
    out["_direction_order"] = out["direction"].map(direction_order)
    out["_metric_order"] = out["metric_name"].map(metric_order)
    out["_representation_order"] = out["representation"].map(rep_order)
    out = out.sort_values(
        [
            "_dataset_order",
            "_cell_line_order",
            "_family_order",
            "_direction_order",
            "_metric_order",
            "rank_value",
            "_representation_order",
        ],
        kind="stable",
        na_position="last",
    )
    return out.drop(
        columns=[
            "_dataset_order",
            "_cell_line_order",
            "_family_order",
            "_direction_order",
            "_metric_order",
            "_representation_order",
        ]
    ).reset_index(drop=True)[BENCHMARK_SUMMARY_LONG_COLUMNS]


def build_stage_manifest(stage_dir: Path) -> List[Dict[str, object]]:
    entries: List[Dict[str, object]] = []
    for file_path in sorted(stage_dir.iterdir()):
        if file_path.is_file() and file_path.name != "manifest.json":
            entries.append(
                {
                    "relative_path": file_path.name,
                    "size_bytes": int(file_path.stat().st_size),
                    "sha256": compute_sha256(file_path),
                }
            )
    return entries


def build_output_assertions(
    *,
    s4_bundle: StageBundle,
    s5_bundle: StageBundle,
    group_long: pd.DataFrame,
    group_leaderboard: pd.DataFrame,
    retrieval_leaderboard: pd.DataFrame,
    benchmark_summary_long: pd.DataFrame,
    non_blocking_notes: Sequence[str],
    stage_dir: Path,
) -> List[Dict[str, Any]]:
    assertions: List[Dict[str, Any]] = []

    expected_group_long_rows = int(len(s4_bundle.tables["task2_group_concordance.csv"]) * len(GROUP_METRIC_ORDER))
    expected_group_leaderboard_rows = int(
        len(s4_bundle.tables["task2_group_concordance.csv"][["dataset", "cell_line", "representation"]].drop_duplicates())
    )
    expected_retrieval_leaderboard_rows = int(len(s5_bundle.tables["task2_retrieval_summary.csv"]))
    expected_group_benchmark_rows = int(expected_group_leaderboard_rows * len(GROUP_SUMMARY_METRIC_ORDER))
    expected_retrieval_benchmark_rows = int(expected_retrieval_leaderboard_rows * len(RETRIEVAL_METRICS))
    expected_benchmark_rows = int(expected_group_benchmark_rows + expected_retrieval_benchmark_rows)

    input_paths = [
        *(str((s4_bundle.stage_dir / name).resolve()) for name in S4_REQUIRED_FILES),
        *(str((s5_bundle.stage_dir / name).resolve()) for name in S5_REQUIRED_FILES),
    ]
    inputs_under_stage_dirs = all(
        is_relative_to(Path(path), s4_bundle.stage_dir) or is_relative_to(Path(path), s5_bundle.stage_dir)
        for path in input_paths
    )
    assertions.append(
        {
            "name": "upstream_stage_contract_locked",
            "pass": True,
            "details": {
                "s4_stage": s4_bundle.run_manifest.get("stage"),
                "s5_stage": s5_bundle.run_manifest.get("stage"),
                "dataset_order": list(DATASET_ORDER),
                "direction_order": list(DIRECTION_ORDER),
                "representation_order": list(REPRESENTATION_ORDER),
                "target_membership_source": TARGET_MEMBERSHIP_SOURCE,
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": bool(inputs_under_stage_dirs),
            "details": {
                "s4_stage_dir": str(s4_bundle.stage_dir),
                "s5_stage_dir": str(s5_bundle.stage_dir),
                "n_inputs": len(input_paths),
            },
            "counterexamples": [] if inputs_under_stage_dirs else [{"inputs": input_paths[:MAX_COUNTEREXAMPLES]}],
        }
    )

    s4_slice_grid = {
        tuple(row)
        for row in s4_bundle.tables["task2_group_concordance.csv"][
            ["dataset", "cell_line", "representation"]
        ].drop_duplicates().itertuples(index=False, name=None)
    }
    s5_slice_grid = {
        tuple(row)
        for row in s5_bundle.tables["task2_retrieval_summary.csv"][
            ["dataset", "cell_line", "representation"]
        ].drop_duplicates().itertuples(index=False, name=None)
    }
    only_in_s4 = sorted(s4_slice_grid - s5_slice_grid)
    only_in_s5 = sorted(s5_slice_grid - s4_slice_grid)
    cross_stage_slice_grid_pass = not only_in_s4 and not only_in_s5
    assertions.append(
        {
            "name": "cross_stage_slice_grid_match",
            "pass": bool(cross_stage_slice_grid_pass),
            "details": {
                "key": ["dataset", "cell_line", "representation"],
                "n_s4_slices": int(len(s4_slice_grid)),
                "n_s5_slices": int(len(s5_slice_grid)),
                "n_intersection": int(len(s4_slice_grid & s5_slice_grid)),
                "n_only_in_s4": int(len(only_in_s4)),
                "n_only_in_s5": int(len(only_in_s5)),
                "sample_only_in_s4": [
                    {"dataset": key[0], "cell_line": key[1], "representation": key[2]}
                    for key in only_in_s4[:MAX_COUNTEREXAMPLES]
                ],
                "sample_only_in_s5": [
                    {"dataset": key[0], "cell_line": key[1], "representation": key[2]}
                    for key in only_in_s5[:MAX_COUNTEREXAMPLES]
                ],
            },
            "counterexamples": (
                [
                    {"source": "s4_only", "dataset": key[0], "cell_line": key[1], "representation": key[2]}
                    for key in only_in_s4[:MAX_COUNTEREXAMPLES]
                ]
                + [
                    {"source": "s5_only", "dataset": key[0], "cell_line": key[1], "representation": key[2]}
                    for key in only_in_s5[:MAX_COUNTEREXAMPLES]
                ]
            ),
        }
    )

    group_long_pass = len(group_long) == expected_group_long_rows
    assertions.append(
        {
            "name": "group_concordance_long_contract",
            "pass": bool(group_long_pass),
            "details": {
                "expected_rows": expected_group_long_rows,
                "actual_rows": int(len(group_long)),
                "sort_order": "dataset_order -> cell_line_order -> target_order -> representation_order -> metric_order",
            },
            "counterexamples": [] if group_long_pass else group_long.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    edist_rank_column_absent = "rank_by_mean_edist_biascorr" not in group_leaderboard.columns
    edist_flags = (
        (group_leaderboard["mean_edist_biascorr_leaderboard_eligible_bool"] == False).all()
        and (group_leaderboard["mean_edist_biascorr_cross_representation_comparable_bool"] == False).all()
    )
    group_leaderboard_pass = len(group_leaderboard) == expected_group_leaderboard_rows and edist_rank_column_absent and edist_flags
    assertions.append(
        {
            "name": "group_leaderboard_contract",
            "pass": bool(group_leaderboard_pass),
            "details": {
                "expected_rows": expected_group_leaderboard_rows,
                "actual_rows": int(len(group_leaderboard)),
                "edist_rank_column_present": not edist_rank_column_absent,
            },
            "counterexamples": [] if group_leaderboard_pass else group_leaderboard.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    retrieval_leaderboard_pass = len(retrieval_leaderboard) == expected_retrieval_leaderboard_rows
    assertions.append(
        {
            "name": "retrieval_leaderboard_contract",
            "pass": bool(retrieval_leaderboard_pass),
            "details": {
                "expected_rows": expected_retrieval_leaderboard_rows,
                "actual_rows": int(len(retrieval_leaderboard)),
                "rank_basis_metric": "mean_mrr_corrected",
            },
            "counterexamples": [] if retrieval_leaderboard_pass else retrieval_leaderboard.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    group_rows = int((benchmark_summary_long["analysis_family"] == "group_concordance").sum())
    retrieval_rows = int((benchmark_summary_long["analysis_family"] == "retrieval").sum())
    edist_rows = benchmark_summary_long.loc[benchmark_summary_long["metric_name"].eq("mean_edist_biascorr")].reset_index(drop=True)
    edist_ok = (
        edist_rows["rank_value"].isna().all()
        and (edist_rows["leaderboard_eligible_bool"] == False).all()
        and (edist_rows["cross_representation_comparable_bool"] == False).all()
        and edist_rows["caution_codes"].astype(str).str.contains("GROUP_EDIST_NOT_CROSS_REP", regex=False).all()
    )
    retrieval_code_ok = benchmark_summary_long.loc[
        benchmark_summary_long["analysis_family"].eq("retrieval"), "caution_codes"
    ].astype(str).str.contains("RETRIEVAL_DIRECTION_SPECIFIC", regex=False).all()
    benchmark_pass = (
        len(benchmark_summary_long) == expected_benchmark_rows
        and group_rows == expected_group_benchmark_rows
        and retrieval_rows == expected_retrieval_benchmark_rows
        and edist_ok
        and retrieval_code_ok
    )
    assertions.append(
        {
            "name": "benchmark_summary_long_contract",
            "pass": bool(benchmark_pass),
            "details": {
                "expected_rows": expected_benchmark_rows,
                "actual_rows": int(len(benchmark_summary_long)),
                "expected_group_rows": expected_group_benchmark_rows,
                "actual_group_rows": group_rows,
                "expected_retrieval_rows": expected_retrieval_benchmark_rows,
                "actual_retrieval_rows": retrieval_rows,
            },
            "counterexamples": [] if benchmark_pass else benchmark_summary_long.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    edist_metric_valid = pd.to_numeric(benchmark_summary_long["n_targets_metric_valid"], errors="coerce")
    edist_targets_total = pd.to_numeric(benchmark_summary_long["n_targets_total"], errors="coerce")
    partial_edist_mask = (
        benchmark_summary_long["metric_name"].eq("mean_edist_biascorr")
        & edist_metric_valid.notna()
        & edist_targets_total.notna()
        & edist_metric_valid.lt(edist_targets_total)
    )
    partial_edist_rows = benchmark_summary_long.loc[partial_edist_mask].reset_index(drop=True)
    partial_edist_codes_ok = (
        partial_edist_rows.empty
        or partial_edist_rows["caution_codes"].astype(str).str.contains("GROUP_EDIST_PARTIAL_TARGET_VALIDITY", regex=False).all()
    )
    assertions.append(
        {
            "name": "edist_partial_validity_note_preserved",
            "pass": bool(partial_edist_codes_ok),
            "details": {"n_partial_edist_rows": int(len(partial_edist_rows))},
            "counterexamples": [] if partial_edist_codes_ok else partial_edist_rows.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    uce_rows = benchmark_summary_long.loc[benchmark_summary_long["representation"].eq("uce")].reset_index(drop=True)
    uce_code_ok = uce_rows.empty or uce_rows["caution_codes"].astype(str).str.contains("UCE_ATTRITION_PRESENT", regex=False).all()
    assertions.append(
        {
            "name": "uce_attrition_note_preserved",
            "pass": bool(uce_code_ok),
            "details": {"n_uce_rows": int(len(uce_rows))},
            "counterexamples": [] if uce_code_ok else uce_rows.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    outputs_root_pass = all(
        is_relative_to(path, stage_dir)
        for path in [
            stage_dir / "task2_group_concordance_long.csv",
            stage_dir / "task2_group_leaderboard.csv",
            stage_dir / "task2_retrieval_leaderboard.csv",
            stage_dir / "task2_benchmark_summary_long.csv",
        ]
    )
    assertions.append(
        {
            "name": "outputs_under_stage_dir",
            "pass": bool(outputs_root_pass),
            "details": {"stage_dir": str(stage_dir)},
            "counterexamples": [] if outputs_root_pass else [{"stage_dir": str(stage_dir)}],
        }
    )

    assertions.append(
        {
            "name": "non_blocking_provenance_notes_surfaced",
            "pass": True,
            "details": {"notes": list(non_blocking_notes)},
            "counterexamples": [],
        }
    )
    return assertions


def load_config(project_root: Path) -> Dict[str, Any]:
    with (project_root / CONFIG_PATH).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    started_at = utc_now_iso()

    config = load_config(project_root)
    seed = int(args.seed if args.seed is not None else config["project"]["seed"])
    if seed != GLOBAL_SEED:
        raise ValueError(f"GLOBAL_SEED must remain {GLOBAL_SEED}, got {seed}")

    runs_dir = resolve_config_path(project_root, str(config["paths"]["runs_dir"]))
    s4_stage_dir = (project_root / args.s4_stage_dir).resolve() if not args.s4_stage_dir.is_absolute() else args.s4_stage_dir.resolve()
    s5_stage_dir = (project_root / args.s5_stage_dir).resolve() if not args.s5_stage_dir.is_absolute() else args.s5_stage_dir.resolve()
    s4_bundle = load_stage_bundle(
        s4_stage_dir,
        S4_REQUIRED_FILES,
        ("task2_group_concordance.csv", "task2_group_attrition.csv"),
    )
    s5_bundle = load_stage_bundle(
        s5_stage_dir,
        S5_REQUIRED_FILES,
        (
            "task2_retrieval_summary.csv",
            "task2_retrieval_summary_long.csv",
            "task2_retrieval_attrition.csv",
            "task2_chance_identity_check.csv",
        ),
    )

    s4_rep_details = validate_s4_bundle(s4_bundle)
    s5_rep_details = validate_s5_bundle(s5_bundle)
    cell_line_order_map = validate_cross_stage_contracts(s4_bundle, s5_bundle, s4_rep_details, s5_rep_details)

    s4_concordance = s4_bundle.tables["task2_group_concordance.csv"]
    s4_attrition = s4_bundle.tables["task2_group_attrition.csv"]
    s5_summary = s5_bundle.tables["task2_retrieval_summary.csv"]
    s5_summary_long = s5_bundle.tables["task2_retrieval_summary_long.csv"]
    s5_attrition = s5_bundle.tables["task2_retrieval_attrition.csv"]
    _ = s5_summary_long, s5_attrition

    group_long = build_group_concordance_long(s4_concordance)
    assert_unique_key(group_long, ["dataset", "cell_line", "target_token", "representation", "metric_name"], "task2_group_concordance_long.csv")
    group_leaderboard = build_group_leaderboard(s4_concordance, s4_attrition, s4_rep_details, cell_line_order_map)
    assert_unique_key(group_leaderboard, ["dataset", "cell_line", "representation"], "task2_group_leaderboard.csv")
    retrieval_leaderboard = build_retrieval_leaderboard(s5_summary, cell_line_order_map)
    assert_unique_key(retrieval_leaderboard, ["dataset", "cell_line", "direction", "representation"], "task2_retrieval_leaderboard.csv")
    benchmark_summary_long = build_benchmark_summary_long(
        group_leaderboard=group_leaderboard,
        retrieval_leaderboard=retrieval_leaderboard,
        s4_rep_details=s4_rep_details,
        s5_rep_details=s5_rep_details,
        s4_run_id=str(s4_bundle.run_manifest["run_id"]),
        s5_run_id=str(s5_bundle.run_manifest["run_id"]),
        cell_line_order_map=cell_line_order_map,
    )
    assert_unique_key(
        benchmark_summary_long,
        ["analysis_family", "dataset", "cell_line", "direction", "representation", "metric_name"],
        "task2_benchmark_summary_long.csv",
    )

    expected_group_long_rows = int(len(s4_concordance) * len(GROUP_METRIC_ORDER))
    expected_group_leaderboard_rows = int(len(s4_concordance[["dataset", "cell_line", "representation"]].drop_duplicates()))
    expected_retrieval_leaderboard_rows = int(len(s5_summary))
    expected_benchmark_rows = int(expected_group_leaderboard_rows * len(GROUP_SUMMARY_METRIC_ORDER) + expected_retrieval_leaderboard_rows * len(RETRIEVAL_METRICS))
    if len(group_long) != expected_group_long_rows:
        raise ValueError(f"task2_group_concordance_long.csv row count mismatch: {len(group_long)}")
    if len(group_leaderboard) != expected_group_leaderboard_rows:
        raise ValueError(f"task2_group_leaderboard.csv row count mismatch: {len(group_leaderboard)}")
    if len(retrieval_leaderboard) != expected_retrieval_leaderboard_rows:
        raise ValueError(f"task2_retrieval_leaderboard.csv row count mismatch: {len(retrieval_leaderboard)}")
    if len(benchmark_summary_long) != expected_benchmark_rows:
        raise ValueError(f"task2_benchmark_summary_long.csv row count mismatch: {len(benchmark_summary_long)}")

    stage_dir = (runs_dir / args.run_id / STAGE).resolve()
    stage_dir.mkdir(parents=True, exist_ok=True)

    group_long_path = stage_dir / "task2_group_concordance_long.csv"
    group_leaderboard_path = stage_dir / "task2_group_leaderboard.csv"
    retrieval_leaderboard_path = stage_dir / "task2_retrieval_leaderboard.csv"
    benchmark_summary_long_path = stage_dir / "task2_benchmark_summary_long.csv"
    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    manifest_path = stage_dir / "manifest.json"

    write_csv(group_long, group_long_path)
    write_csv(group_leaderboard, group_leaderboard_path)
    write_csv(retrieval_leaderboard, retrieval_leaderboard_path)
    write_csv(benchmark_summary_long, benchmark_summary_long_path)

    non_blocking_notes = extract_non_blocking_notes(s4_bundle.run_manifest) + extract_non_blocking_notes(s5_bundle.run_manifest)
    if s4_bundle.run_manifest.get("git_head") != s5_bundle.run_manifest.get("git_head"):
        non_blocking_notes.append("upstream_git_head_mismatch: S4 and S5 audited runs were generated from different git heads")

    assertions = build_output_assertions(
        s4_bundle=s4_bundle,
        s5_bundle=s5_bundle,
        group_long=group_long,
        group_leaderboard=group_leaderboard,
        retrieval_leaderboard=retrieval_leaderboard,
        benchmark_summary_long=benchmark_summary_long,
        non_blocking_notes=non_blocking_notes,
        stage_dir=stage_dir,
    )
    failed_assertions = [item["name"] for item in assertions if not bool(item.get("pass", False))]
    if failed_assertions:
        raise ValueError(f"Internal Script6 multisource assertion failure: {failed_assertions[:MAX_COUNTEREXAMPLES]}")

    cell_line_order_by_dataset = {
        dataset: [
            cell_line
            for cell_line, _ in sorted(order_map.items(), key=lambda item: item[1])
        ]
        for dataset, order_map in cell_line_order_map.items()
    }
    target_order_by_slice = {
        f"{dataset}||{cell_line}": targets
        for (dataset, cell_line), targets in extract_target_order_by_slice(s4_concordance).items()
    }

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts" / "s6_task2_result_synthesis_multisource.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": utc_now_iso(),
        "config": {
            "seed": seed,
            "runs_dir": str(runs_dir),
            "task2_snapshot_contract": EXPECTED_TASK2_SNAPSHOT.as_posix(),
            "dataset_order": list(DATASET_ORDER),
            "direction_order": list(DIRECTION_ORDER),
            "representation_order": list(REPRESENTATION_ORDER),
            "group_metric_order": list(GROUP_METRIC_ORDER),
            "group_summary_metric_order": list(GROUP_SUMMARY_METRIC_ORDER),
            "retrieval_metric_order": list(RETRIEVAL_METRICS),
            "target_membership_source": TARGET_MEMBERSHIP_SOURCE,
            "s4_stage_dir": str(s4_bundle.stage_dir),
            "s5_stage_dir": str(s5_bundle.stage_dir),
            "s4_source_stage": EXPECTED_S4_STAGE,
            "s5_source_stage": EXPECTED_S5_STAGE,
            "s4_source_run_id": str(s4_bundle.run_manifest["run_id"]),
            "s5_source_run_id": str(s5_bundle.run_manifest["run_id"]),
            "expected_group_concordance_long_rows": expected_group_long_rows,
            "expected_group_leaderboard_rows": expected_group_leaderboard_rows,
            "expected_retrieval_leaderboard_rows": expected_retrieval_leaderboard_rows,
            "expected_benchmark_summary_long_rows": expected_benchmark_rows,
            "non_blocking_provenance_notes": non_blocking_notes,
        },
        "inputs": [
            str((s4_bundle.stage_dir / name).resolve()) for name in S4_REQUIRED_FILES
        ]
        + [
            str((s5_bundle.stage_dir / name).resolve()) for name in S5_REQUIRED_FILES
        ],
        "outputs": [
            str(group_long_path),
            str(group_leaderboard_path),
            str(retrieval_leaderboard_path),
            str(benchmark_summary_long_path),
            str(run_manifest_path),
            str(audit_assertions_path),
            str(manifest_path),
        ],
        "summary": {
            "n_group_concordance_rows_source": int(len(s4_concordance)),
            "n_group_attrition_rows_source": int(len(s4_attrition)),
            "n_group_concordance_long_rows": int(len(group_long)),
            "n_group_leaderboard_rows": int(len(group_leaderboard)),
            "n_retrieval_summary_rows_source": int(len(s5_summary)),
            "n_retrieval_summary_long_rows_source": int(len(s5_summary_long)),
            "n_retrieval_attrition_rows_source": int(len(s5_attrition)),
            "n_retrieval_leaderboard_rows": int(len(retrieval_leaderboard)),
            "n_benchmark_summary_long_rows": int(len(benchmark_summary_long)),
            "benchmark_group_rows": int((benchmark_summary_long["analysis_family"] == "group_concordance").sum()),
            "benchmark_retrieval_rows": int((benchmark_summary_long["analysis_family"] == "retrieval").sum()),
            "cell_line_order_by_dataset": cell_line_order_by_dataset,
            "target_order_by_slice": target_order_by_slice,
            "non_blocking_provenance_notes": non_blocking_notes,
        },
    }
    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})
    write_json(manifest_path, {"stage": STAGE, "files": build_stage_manifest(stage_dir)})


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)
