#!/usr/bin/env python3
# SCRIPT_HEADER_CONTRACT
# Script: scripts/s7_project_benchmark_synthesis.py
# Purpose: Build project-level Task1+Task2 benchmark synthesis tables from
#   audited Task1 S1/S2 reporting outputs and corrected Task2 S6 reporting
#   outputs without recomputing any upstream metric.
# Inputs:
#   - audited S1 stage dir with:
#     - task1_leaderboard_long.csv
#     - task1_retrieval_summary.csv
#     - task1_chance_identity_check.csv
#     - task1_attrition.csv
#     - run_manifest.json
#     - audit_assertions.json
#     - manifest.json
#   - audited S2 stage dir with:
#     - task1_cross_leaderboard_long.csv
#     - task1_cross_retrieval_summary.csv
#     - task1_cross_chance_identity_check.csv
#     - task1_cross_alignment_proof.csv
#     - task1_cross_attrition.csv
#     - run_manifest.json
#     - audit_assertions.json
#     - manifest.json
#   - audited S6 stage dir with:
#     - task2_benchmark_summary_long.csv
#     - task2_group_leaderboard.csv
#     - task2_retrieval_leaderboard.csv
#     - task2_group_concordance_long.csv
#     - run_manifest.json
#     - audit_assertions.json
#     - manifest.json
# Outputs:
#   - project_input_registry.csv: runs/<run_id>/s7_project_benchmark_synthesis/
#   - project_benchmark_summary_long.csv: runs/<run_id>/s7_project_benchmark_synthesis/
#   - project_axis_score_inputs_long.csv: runs/<run_id>/s7_project_benchmark_synthesis/
#   - project_representation_scorecard.csv: runs/<run_id>/s7_project_benchmark_synthesis/
#   - AVCP artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s7_project_benchmark_synthesis/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s7_project_benchmark_synthesis.py --run-id <run_id> --task1-s1-stage-dir <path> --task1-s2-stage-dir <path> --task2-s6-stage-dir <path> --seed 619
# Failure Modes:
#   - Missing or unaudited upstream files -> exit non-zero
#   - Reporting-schema harmonization drift -> exit non-zero
#   - Observed-vs-contract-exclusion leakage -> exit non-zero
#   - Unique-key / sort-order / percentile / caution aggregation drift -> exit non-zero
# Last Updated: 2026-03-11

"""
Inputs:
- audited Task1 S1/S2 reporting-stage outputs only
- audited corrected Task2 S6 reporting-stage outputs only

Outputs:
- project_input_registry.csv
- project_benchmark_summary_long.csv
- project_axis_score_inputs_long.csv
- project_representation_scorecard.csv
- run_manifest.json
- audit_assertions.json
- manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- Task order = Task1, Task2
- Task scope order = internal, cross, mechanism
- Analysis family order = group_concordance, retrieval
- Dataset order = LINCS, scPerturb, NA
- Cross direction order = LINCS_to_scPerturb, scPerturb_to_LINCS, NA
- Direction order = NA, C2G, G2C
- Perturbation type order = Chemical, Genetic, NA
- Metric order = mean_cosine_centroid, mean_pcc_centroid, mean_edist_biascorr, mean_mrr_corrected, mean_hit1_corrected, mean_hit5_corrected, mean_hit10_corrected
- representation_canonical = prefix-strip-only on FM:* tokens
- scorecard mean = equal-family macro average over observed family-level mean_mrr_corrected

Interpretation rules:
- S7 is synthesis only and does not recompute any Task1/Task2 metric
- project_benchmark_summary_long contains observed upstream reporting rows only
- project_axis_score_inputs_long is the only output allowed to materialize contract-exclusion rows
- project_representation_scorecard is derived only from project_axis_score_inputs_long
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml

STAGE = "s7_project_benchmark_synthesis"
CONFIG_PATH = Path("config/config.yaml")
GLOBAL_SEED = 619
MAX_COUNTEREXAMPLES = 5
S1_SOFT_ASSERTION = "fm_policy_presence_soft"
S1_SOFT_CAUTION = "INHERITED_STAGE_SOFT_AUDIT_S1"
TASK1_CROSS_EXCLUSION_CAUTION = "TASK1_CROSS_CHEMICAL_CONTRACT_EXCLUDED"
PERCENTILE_TIE_RULE = "average_rank_desc_high_is_better"

EXPECTED_S1_STAGE = "s1_task1_internal_metrics"
EXPECTED_S2_STAGE = "s2_task1_cross_metrics"
EXPECTED_S4_STAGE = "s4_task2_group_concordance_multisource"
EXPECTED_S5_STAGE = "s5_task2_retrieval_multisource"
EXPECTED_S6_STAGE = "s6_task2_result_synthesis_multisource"

TASK_ORDER: Tuple[str, ...] = ("Task1", "Task2")
TASK_SCOPE_ORDER: Tuple[str, ...] = ("internal", "cross", "mechanism")
ANALYSIS_FAMILY_ORDER: Tuple[str, ...] = ("group_concordance", "retrieval")
DATASET_ORDER: Tuple[str, ...] = ("LINCS", "scPerturb", "NA")
CROSS_DIRECTION_ORDER: Tuple[str, ...] = ("LINCS_to_scPerturb", "scPerturb_to_LINCS", "NA")
DIRECTION_ORDER: Tuple[str, ...] = ("NA", "C2G", "G2C")
PERTURBATION_TYPE_ORDER: Tuple[str, ...] = ("Chemical", "Genetic", "NA")
METRIC_ORDER: Tuple[str, ...] = (
    "mean_cosine_centroid",
    "mean_pcc_centroid",
    "mean_edist_biascorr",
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
)
PROVENANCE_ROLE_ORDER: Tuple[str, ...] = ("direct_ingest_stage", "transitive_provenance_stage")
ROW_ORIGIN_ORDER: Tuple[str, ...] = ("observed_upstream", "contract_exclusion")
STAGE_ORDER: Tuple[str, ...] = (
    EXPECTED_S1_STAGE,
    EXPECTED_S2_STAGE,
    EXPECTED_S6_STAGE,
    EXPECTED_S4_STAGE,
    EXPECTED_S5_STAGE,
)

TASK1_GROUP_METRICS = {"mean_cosine_centroid", "mean_pcc_centroid", "mean_edist_biascorr"}
TASK1_RETRIEVAL_METRICS = {"mean_mrr_corrected", "mean_hit1_corrected", "mean_hit5_corrected", "mean_hit10_corrected"}
S2_ALLOWED_REPRESENTATIONS: Tuple[str, ...] = ("Gene", "Pathway")
TASK1_CROSS_ALLOWED_DIRECTIONS: Tuple[str, ...] = ("LINCS_to_scPerturb", "scPerturb_to_LINCS")
S6_ALLOWED_SOURCE_STAGES = {EXPECTED_S4_STAGE, EXPECTED_S5_STAGE}
S6_ALLOWED_SOURCE_TABLES = {"task2_group_leaderboard.csv", "task2_retrieval_leaderboard.csv"}
AXIS_ALLOWED_COVERAGE_STATUS = {"observed", "contract_excluded"}
SCORECARD_AXIS_BY_TASK = {"Task1": "Task1", "Task2": "Task2"}

S1_REQUIRED_FILES: Tuple[str, ...] = (
    "task1_leaderboard_long.csv",
    "task1_retrieval_summary.csv",
    "task1_chance_identity_check.csv",
    "task1_attrition.csv",
    "run_manifest.json",
    "audit_assertions.json",
    "manifest.json",
)
S2_REQUIRED_FILES: Tuple[str, ...] = (
    "task1_cross_leaderboard_long.csv",
    "task1_cross_retrieval_summary.csv",
    "task1_cross_chance_identity_check.csv",
    "task1_cross_alignment_proof.csv",
    "task1_cross_attrition.csv",
    "run_manifest.json",
    "audit_assertions.json",
    "manifest.json",
)
S4_REQUIRED_FILES: Tuple[str, ...] = ("run_manifest.json", "audit_assertions.json", "manifest.json")
S5_REQUIRED_FILES: Tuple[str, ...] = ("run_manifest.json", "audit_assertions.json", "manifest.json")
S6_REQUIRED_FILES: Tuple[str, ...] = (
    "task2_benchmark_summary_long.csv",
    "task2_group_leaderboard.csv",
    "task2_retrieval_leaderboard.csv",
    "task2_group_concordance_long.csv",
    "run_manifest.json",
    "audit_assertions.json",
    "manifest.json",
)

TASK1_LEADERBOARD_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "metric_name",
    "metric_value",
    "n_total",
    "n_valid",
    "n_excluded",
    "N_gallery_max",
    "cross_alignment_contract",
]
TASK1_RETRIEVAL_SUMMARY_COLUMNS = [
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
]
TASK1_CHANCE_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
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
TASK1_ATTRITION_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "reason",
    "n_dropped",
    "n_total_before",
    "notes",
]
TASK1_CROSS_ALIGNMENT_COLUMNS = [
    "cross_alignment_contract",
    "perturbation_type",
    "n_matched_keys",
    "eligible_bool",
    "excluded_reason",
]

S6_GROUP_LEADERBOARD_COLUMNS = [
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
S6_RETRIEVAL_LEADERBOARD_COLUMNS = [
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
S6_GROUP_LONG_COLUMNS = [
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
S6_BENCHMARK_COLUMNS = [
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

PROJECT_INPUT_REGISTRY_COLUMNS = [
    "provenance_role",
    "task",
    "stage_name",
    "run_id",
    "stage_dir",
    "reporting_tables_json",
    "validation_tables_json",
    "manifest_path",
    "audit_assertions_path",
    "stage_manifest_path",
    "audit_status",
    "non_pass_assertions_json",
    "notes",
]
PROJECT_BENCHMARK_SUMMARY_COLUMNS = [
    "task",
    "task_scope",
    "analysis_family",
    "dataset",
    "cross_direction",
    "cell_line",
    "perturbation_type",
    "direction",
    "representation_raw",
    "representation_canonical",
    "metric_name",
    "metric_value",
    "row_origin",
    "rank_value",
    "rank_basis_metric_name",
    "leaderboard_eligible_bool",
    "cross_representation_comparable_bool",
    "n_total",
    "n_valid",
    "n_excluded",
    "n_excluded_missing_metric_or_mpos0",
    "n_targets_total",
    "n_targets_metric_valid",
    "N_gallery_mean",
    "N_gallery_max",
    "m_pos_mean",
    "m_pos_p50",
    "m_pos_p90",
    "cross_alignment_contract",
    "ingest_table",
    "ingest_stage",
    "ingest_run_id",
    "source_table",
    "source_stage",
    "source_run_id",
    "caution_codes",
]
PROJECT_AXIS_COLUMNS = [
    "task_axis",
    "family_id",
    "row_origin",
    "task",
    "task_scope",
    "dataset",
    "cross_direction",
    "cell_line",
    "perturbation_type",
    "direction",
    "representation_raw",
    "representation_canonical",
    "metric_name",
    "metric_value",
    "coverage_status",
    "aggregation_weight",
    "cross_alignment_contract",
    "ingest_table",
    "ingest_stage",
    "ingest_run_id",
    "source_table",
    "source_stage",
    "source_run_id",
    "caution_codes",
]
PROJECT_SCORECARD_COLUMNS = [
    "representation_canonical",
    "task1_axis_value_raw",
    "task1_axis_observed_families",
    "task1_axis_contract_exclusion_families",
    "task1_axis_expected_families",
    "task2_axis_value_raw",
    "task2_axis_observed_families",
    "task2_axis_contract_exclusion_families",
    "task2_axis_expected_families",
    "x_percentile",
    "y_percentile",
    "scorecard_eligible_bool",
    "caution_codes",
]


@dataclass(frozen=True)
class StageBundle:
    stage_dir: Path
    paths: Dict[str, Path]
    run_manifest: Dict[str, Any]
    audit_assertions: Dict[str, Any]
    manifest: Dict[str, Any]
    tables: Dict[str, pd.DataFrame]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S7 project-level benchmark synthesis")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--task1-s1-stage-dir", type=Path, required=True)
    parser.add_argument("--task1-s2-stage-dir", type=Path, required=True)
    parser.add_argument("--task2-s6-stage-dir", type=Path, required=True)
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


def fill_null_sentinel(frame: pd.DataFrame, columns: Sequence[str], sentinel: str = "__NULL__") -> pd.DataFrame:
    out = frame.loc[:, list(columns)].copy()
    for column in columns:
        out[column] = out[column].astype(object).where(~pd.isna(out[column]), sentinel)
    return out


def assert_unique_key(frame: pd.DataFrame, columns: Sequence[str], name: str) -> None:
    if frame.empty:
        return
    norm = fill_null_sentinel(frame, columns)
    duplicated = norm.duplicated(list(columns), keep=False)
    if bool(duplicated.any()):
        examples = frame.loc[duplicated, list(columns)].drop_duplicates().head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"{name} is not unique on key={list(columns)}: {examples}")


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
    paths = require_stage_files(stage_dir.resolve(), required_files)
    tables = {name: pd.read_csv(paths[name]) for name in csv_files}
    return StageBundle(
        stage_dir=stage_dir.resolve(),
        paths=paths,
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


def extract_assertions(bundle: StageBundle, stage_name: str) -> List[Dict[str, Any]]:
    raw = bundle.audit_assertions.get("assertions")
    if not isinstance(raw, list):
        raise ValueError(f"{stage_name} audit_assertions.json missing 'assertions' list")
    return [dict(item) for item in raw if isinstance(item, Mapping)]


def validate_all_assertions_pass(bundle: StageBundle, stage_name: str) -> None:
    failing = [item for item in extract_assertions(bundle, stage_name) if not bool(item.get("pass", False))]
    if failing:
        names = [str(item.get("name", "<unnamed>")) for item in failing[:MAX_COUNTEREXAMPLES]]
        raise ValueError(f"{stage_name} has failing upstream assertions: {names}")


def validate_s1_assertions(bundle: StageBundle) -> Tuple[List[str], bool]:
    failing = [item for item in extract_assertions(bundle, EXPECTED_S1_STAGE) if not bool(item.get("pass", False))]
    names = [str(item.get("name", "<unnamed>")) for item in failing]
    unexpected = [name for name in names if name != S1_SOFT_ASSERTION]
    if unexpected:
        raise ValueError(f"{EXPECTED_S1_STAGE} has unexpected failing assertions: {unexpected}")
    return names, S1_SOFT_ASSERTION in names


def first_seen_values(series: Iterable[str]) -> List[str]:
    seen: List[str] = []
    for value in series:
        if value not in seen:
            seen.append(value)
    return seen


def parse_caution_codes(value: object) -> List[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    raw = str(value).strip()
    if not raw:
        return []
    if raw.startswith("["):
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError:
            parsed = None
        if isinstance(parsed, list):
            return [str(item) for item in parsed if str(item).strip()]
    return [part.strip() for part in raw.split(";") if part.strip()]


def json_array_string(values: Iterable[str]) -> str:
    items = sorted({str(value) for value in values if str(value).strip()})
    return json.dumps(items, ensure_ascii=True)


def canonicalize_representation(raw: str) -> str:
    value = str(raw)
    if value.startswith("FM:"):
        return value[3:]
    return value


def read_config(project_root: Path) -> Tuple[int, Path]:
    with (project_root / CONFIG_PATH).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, Mapping):
        raise ValueError("config/config.yaml must decode to a mapping")
    project = config.get("project")
    paths = config.get("paths")
    if not isinstance(project, Mapping) or not isinstance(paths, Mapping):
        raise ValueError("config/config.yaml missing project/paths mapping blocks")
    seed = int(project.get("seed", -1))
    if seed != GLOBAL_SEED:
        raise ValueError(f"config seed drift detected: {seed}")
    runs_dir = resolve_config_path(project_root, str(paths.get("runs_dir", "")))
    return seed, runs_dir


def assert_chance_tolerance(frame: pd.DataFrame, name: str, tol: float = 1e-12) -> None:
    for column in ("abs_delta_mrr", "abs_delta_hit1", "abs_delta_hit5", "abs_delta_hit10"):
        if column not in frame.columns:
            raise ValueError(f"{name} missing chance tolerance column: {column}")
        numeric = pd.to_numeric(frame[column], errors="coerce")
        if bool((numeric > tol).fillna(False).any()):
            bad = frame.loc[(numeric > tol).fillna(False), :].head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
            raise ValueError(f"{name} chance tolerance exceeded for {column}: {bad}")


def validate_s1_bundle(bundle: StageBundle) -> bool:
    if bundle.run_manifest.get("stage") != EXPECTED_S1_STAGE:
        raise ValueError(f"S1 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S1_REQUIRED_FILES, EXPECTED_S1_STAGE)
    non_pass, has_soft = validate_s1_assertions(bundle)
    if non_pass and non_pass != [S1_SOFT_ASSERTION]:
        raise ValueError(f"S1 non-pass assertions drift detected: {non_pass}")
    config = bundle.run_manifest.get("config")
    if not isinstance(config, Mapping):
        raise ValueError("S1 run_manifest missing config block")
    if int(config.get("seed", -1)) != GLOBAL_SEED:
        raise ValueError(f"S1 seed drift detected: {config.get('seed')}")

    leaderboard = bundle.tables["task1_leaderboard_long.csv"]
    summary = bundle.tables["task1_retrieval_summary.csv"]
    chance = bundle.tables["task1_chance_identity_check.csv"]
    attrition = bundle.tables["task1_attrition.csv"]
    ensure_required_columns(leaderboard, TASK1_LEADERBOARD_COLUMNS, "S1 task1_leaderboard_long.csv")
    ensure_required_columns(summary, TASK1_RETRIEVAL_SUMMARY_COLUMNS, "S1 task1_retrieval_summary.csv")
    ensure_required_columns(chance, TASK1_CHANCE_COLUMNS, "S1 task1_chance_identity_check.csv")
    ensure_required_columns(attrition, TASK1_ATTRITION_COLUMNS, "S1 task1_attrition.csv")
    assert_unique_key(leaderboard, ["scope", "dataset_or_direction", "perturbation_type", "representation", "metric_name"], "S1 task1_leaderboard_long.csv")
    assert_unique_key(summary, ["scope", "dataset_or_direction", "perturbation_type", "representation"], "S1 task1_retrieval_summary.csv")
    assert_unique_key(chance, ["scope", "dataset_or_direction", "perturbation_type", "representation"], "S1 task1_chance_identity_check.csv")
    assert_unique_key(attrition, ["scope", "dataset_or_direction", "perturbation_type", "representation", "reason"], "S1 task1_attrition.csv")
    if set(leaderboard["scope"].astype(str).unique()) != {"internal"}:
        raise ValueError("S1 leaderboard scope drift detected")
    if set(summary["scope"].astype(str).unique()) != {"internal"}:
        raise ValueError("S1 retrieval_summary scope drift detected")
    metric_names = set(leaderboard["metric_name"].astype(str).unique())
    if metric_names != TASK1_GROUP_METRICS | TASK1_RETRIEVAL_METRICS:
        raise ValueError(f"S1 leaderboard metric drift detected: {sorted(metric_names)}")
    assert_chance_tolerance(chance, "S1 task1_chance_identity_check.csv")
    return has_soft


def validate_s2_bundle(bundle: StageBundle) -> None:
    if bundle.run_manifest.get("stage") != EXPECTED_S2_STAGE:
        raise ValueError(f"S2 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S2_REQUIRED_FILES, EXPECTED_S2_STAGE)
    validate_all_assertions_pass(bundle, EXPECTED_S2_STAGE)
    config = bundle.run_manifest.get("config")
    if not isinstance(config, Mapping):
        raise ValueError("S2 run_manifest missing config block")
    if int(config.get("seed", -1)) != GLOBAL_SEED:
        raise ValueError(f"S2 seed drift detected: {config.get('seed')}")

    leaderboard = bundle.tables["task1_cross_leaderboard_long.csv"]
    summary = bundle.tables["task1_cross_retrieval_summary.csv"]
    chance = bundle.tables["task1_cross_chance_identity_check.csv"]
    alignment = bundle.tables["task1_cross_alignment_proof.csv"]
    attrition = bundle.tables["task1_cross_attrition.csv"]
    ensure_required_columns(leaderboard, TASK1_LEADERBOARD_COLUMNS, "S2 task1_cross_leaderboard_long.csv")
    ensure_required_columns(summary, TASK1_RETRIEVAL_SUMMARY_COLUMNS, "S2 task1_cross_retrieval_summary.csv")
    ensure_required_columns(chance, TASK1_CHANCE_COLUMNS, "S2 task1_cross_chance_identity_check.csv")
    ensure_required_columns(alignment, TASK1_CROSS_ALIGNMENT_COLUMNS, "S2 task1_cross_alignment_proof.csv")
    ensure_required_columns(attrition, TASK1_ATTRITION_COLUMNS, "S2 task1_cross_attrition.csv")
    assert_unique_key(leaderboard, ["scope", "dataset_or_direction", "perturbation_type", "representation", "metric_name"], "S2 task1_cross_leaderboard_long.csv")
    assert_unique_key(summary, ["scope", "dataset_or_direction", "perturbation_type", "representation"], "S2 task1_cross_retrieval_summary.csv")
    assert_unique_key(chance, ["scope", "dataset_or_direction", "perturbation_type", "representation"], "S2 task1_cross_chance_identity_check.csv")
    assert_unique_key(alignment, ["cross_alignment_contract", "perturbation_type"], "S2 task1_cross_alignment_proof.csv")
    assert_unique_key(attrition, ["scope", "dataset_or_direction", "perturbation_type", "representation", "reason"], "S2 task1_cross_attrition.csv")
    if set(leaderboard["scope"].astype(str).unique()) != {"cross"}:
        raise ValueError("S2 leaderboard scope drift detected")
    if set(summary["scope"].astype(str).unique()) != {"cross"}:
        raise ValueError("S2 retrieval_summary scope drift detected")
    reps = first_seen_values(leaderboard["representation"].astype(str).tolist())
    if reps != list(S2_ALLOWED_REPRESENTATIONS):
        raise ValueError(f"S2 representation drift detected: {reps}")
    directions = first_seen_values(summary["dataset_or_direction"].astype(str).tolist())
    if directions != list(TASK1_CROSS_ALLOWED_DIRECTIONS):
        raise ValueError(f"S2 cross direction drift detected: {directions}")
    metric_names = set(leaderboard["metric_name"].astype(str).unique())
    if metric_names != TASK1_GROUP_METRICS | TASK1_RETRIEVAL_METRICS:
        raise ValueError(f"S2 leaderboard metric drift detected: {sorted(metric_names)}")
    assert_chance_tolerance(chance, "S2 task1_cross_chance_identity_check.csv")
    chemical_rows = alignment.loc[alignment["perturbation_type"].astype(str).eq("Chemical")].reset_index(drop=True)
    genetic_rows = alignment.loc[alignment["perturbation_type"].astype(str).eq("Genetic")].reset_index(drop=True)
    if len(chemical_rows) != 1 or len(genetic_rows) != 1:
        raise ValueError("S2 alignment proof perturbation_type grid drift detected")
    chemical_row = chemical_rows.iloc[0].to_dict()
    genetic_row = genetic_rows.iloc[0].to_dict()
    if bool(chemical_row["eligible_bool"]) or str(chemical_row["excluded_reason"]) != "matched_keys_lt5":
        raise ValueError(f"S2 chemical exclusion contract drift detected: {chemical_row}")
    if not bool(genetic_row["eligible_bool"]):
        raise ValueError(f"S2 genetic eligibility drift detected: {genetic_row}")


def validate_s4_transitive_bundle(bundle: StageBundle) -> None:
    if bundle.run_manifest.get("stage") != EXPECTED_S4_STAGE:
        raise ValueError(f"S4 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S4_REQUIRED_FILES, EXPECTED_S4_STAGE)
    validate_all_assertions_pass(bundle, EXPECTED_S4_STAGE)


def validate_s5_transitive_bundle(bundle: StageBundle) -> None:
    if bundle.run_manifest.get("stage") != EXPECTED_S5_STAGE:
        raise ValueError(f"S5 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S5_REQUIRED_FILES, EXPECTED_S5_STAGE)
    validate_all_assertions_pass(bundle, EXPECTED_S5_STAGE)


def validate_s6_bundle(bundle: StageBundle) -> Tuple[Path, Path, Dict[str, List[str]]]:
    if bundle.run_manifest.get("stage") != EXPECTED_S6_STAGE:
        raise ValueError(f"S6 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S6_REQUIRED_FILES, EXPECTED_S6_STAGE)
    validate_all_assertions_pass(bundle, EXPECTED_S6_STAGE)
    config = bundle.run_manifest.get("config")
    summary = bundle.run_manifest.get("summary")
    if not isinstance(config, Mapping):
        raise ValueError("S6 run_manifest missing config block")
    if not isinstance(summary, Mapping):
        raise ValueError("S6 run_manifest missing summary block")
    if int(config.get("seed", -1)) != GLOBAL_SEED:
        raise ValueError(f"S6 seed drift detected: {config.get('seed')}")

    benchmark = bundle.tables["task2_benchmark_summary_long.csv"]
    group_leaderboard = bundle.tables["task2_group_leaderboard.csv"]
    retrieval_leaderboard = bundle.tables["task2_retrieval_leaderboard.csv"]
    group_long = bundle.tables["task2_group_concordance_long.csv"]
    ensure_required_columns(benchmark, S6_BENCHMARK_COLUMNS, "S6 task2_benchmark_summary_long.csv")
    ensure_required_columns(group_leaderboard, S6_GROUP_LEADERBOARD_COLUMNS, "S6 task2_group_leaderboard.csv")
    ensure_required_columns(retrieval_leaderboard, S6_RETRIEVAL_LEADERBOARD_COLUMNS, "S6 task2_retrieval_leaderboard.csv")
    ensure_required_columns(group_long, S6_GROUP_LONG_COLUMNS, "S6 task2_group_concordance_long.csv")
    assert_unique_key(benchmark, ["analysis_family", "dataset", "cell_line", "direction", "representation", "metric_name"], "S6 task2_benchmark_summary_long.csv")
    assert_unique_key(group_leaderboard, ["dataset", "cell_line", "representation"], "S6 task2_group_leaderboard.csv")
    assert_unique_key(retrieval_leaderboard, ["dataset", "cell_line", "direction", "representation"], "S6 task2_retrieval_leaderboard.csv")
    assert_unique_key(group_long, ["dataset", "cell_line", "target_token", "representation", "metric_name"], "S6 task2_group_concordance_long.csv")

    if int(summary.get("n_benchmark_summary_long_rows", -1)) != int(len(benchmark)):
        raise ValueError("S6 benchmark row-count drift detected")
    if int(summary.get("n_group_leaderboard_rows", -1)) != int(len(group_leaderboard)):
        raise ValueError("S6 group_leaderboard row-count drift detected")
    if int(summary.get("n_retrieval_leaderboard_rows", -1)) != int(len(retrieval_leaderboard)):
        raise ValueError("S6 retrieval_leaderboard row-count drift detected")
    if int(summary.get("n_group_concordance_long_rows", -1)) != int(len(group_long)):
        raise ValueError("S6 group_concordance_long row-count drift detected")
    expected_benchmark_rows = int(len(group_leaderboard) * 3 + len(retrieval_leaderboard) * 4)
    if len(benchmark) != expected_benchmark_rows:
        raise ValueError(f"S6 benchmark expected rows drift: actual={len(benchmark)} expected={expected_benchmark_rows}")
    source_stage_values = set(benchmark["source_stage"].astype(str).unique())
    if not source_stage_values.issubset(S6_ALLOWED_SOURCE_STAGES):
        raise ValueError(f"S6 source_stage drift detected: {sorted(source_stage_values)}")
    source_table_values = set(benchmark["source_table"].astype(str).unique())
    if not source_table_values.issubset(S6_ALLOWED_SOURCE_TABLES):
        raise ValueError(f"S6 source_table drift detected: {sorted(source_table_values)}")
    s4_stage_dir = Path(str(config.get("s4_stage_dir", ""))).resolve()
    s5_stage_dir = Path(str(config.get("s5_stage_dir", ""))).resolve()
    if config.get("s4_source_stage") != EXPECTED_S4_STAGE:
        raise ValueError("S6 s4_source_stage drift detected")
    if config.get("s5_source_stage") != EXPECTED_S5_STAGE:
        raise ValueError("S6 s5_source_stage drift detected")
    cell_line_order_by_dataset = summary.get("cell_line_order_by_dataset")
    if not isinstance(cell_line_order_by_dataset, Mapping):
        raise ValueError("S6 run_manifest missing cell_line_order_by_dataset")
    out: Dict[str, List[str]] = {}
    for dataset, values in cell_line_order_by_dataset.items():
        if not isinstance(values, list):
            raise ValueError(f"S6 cell_line_order_by_dataset[{dataset}] must be a list")
        out[str(dataset)] = [str(value) for value in values]
    return s4_stage_dir, s5_stage_dir, out


def build_representation_canonical_order(
    s1_leaderboard: pd.DataFrame,
    s2_leaderboard: pd.DataFrame,
    s6_benchmark: pd.DataFrame,
) -> List[str]:
    seen: List[str] = []
    for frame in (s1_leaderboard, s2_leaderboard):
        for raw in frame["representation"].astype(str).tolist():
            canonical = canonicalize_representation(raw)
            if canonical not in seen:
                seen.append(canonical)
    for raw in s6_benchmark["representation"].astype(str).tolist():
        canonical = canonicalize_representation(raw)
        if canonical not in seen:
            seen.append(canonical)
    return seen


def build_cell_line_order_map(cell_line_order_by_dataset: Mapping[str, List[str]]) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    for dataset, values in cell_line_order_by_dataset.items():
        out[str(dataset)] = {str(value): idx for idx, value in enumerate(values)}
    return out


def metric_to_family(metric_name: str) -> str:
    if metric_name in TASK1_GROUP_METRICS:
        return "group_concordance"
    if metric_name in TASK1_RETRIEVAL_METRICS:
        return "retrieval"
    raise ValueError(f"Unexpected Task1 metric_name for analysis_family mapping: {metric_name}")


def build_s1_benchmark_rows(bundle: StageBundle, has_soft_caution: bool) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    caution = S1_SOFT_CAUTION if has_soft_caution else ""
    for row in bundle.tables["task1_leaderboard_long.csv"].to_dict(orient="records"):
        task_scope = str(row["scope"])
        if task_scope != "internal":
            raise ValueError(f"S1 unexpected scope row: {row}")
        raw_representation = str(row["representation"])
        rows.append(
            {
                "task": "Task1",
                "task_scope": "internal",
                "analysis_family": metric_to_family(str(row["metric_name"])),
                "dataset": str(row["dataset_or_direction"]),
                "cross_direction": None,
                "cell_line": None,
                "perturbation_type": str(row["perturbation_type"]),
                "direction": None,
                "representation_raw": raw_representation,
                "representation_canonical": canonicalize_representation(raw_representation),
                "metric_name": str(row["metric_name"]),
                "metric_value": row["metric_value"],
                "row_origin": "observed_upstream",
                "rank_value": None,
                "rank_basis_metric_name": None,
                "leaderboard_eligible_bool": None,
                "cross_representation_comparable_bool": None,
                "n_total": row["n_total"],
                "n_valid": row["n_valid"],
                "n_excluded": row["n_excluded"],
                "n_excluded_missing_metric_or_mpos0": None,
                "n_targets_total": None,
                "n_targets_metric_valid": None,
                "N_gallery_mean": None,
                "N_gallery_max": row["N_gallery_max"],
                "m_pos_mean": None,
                "m_pos_p50": None,
                "m_pos_p90": None,
                "cross_alignment_contract": row["cross_alignment_contract"],
                "ingest_table": "task1_leaderboard_long.csv",
                "ingest_stage": EXPECTED_S1_STAGE,
                "ingest_run_id": str(bundle.run_manifest["run_id"]),
                "source_table": "task1_leaderboard_long.csv",
                "source_stage": EXPECTED_S1_STAGE,
                "source_run_id": str(bundle.run_manifest["run_id"]),
                "caution_codes": caution,
            }
        )
    return rows


def build_s2_benchmark_rows(bundle: StageBundle) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for row in bundle.tables["task1_cross_leaderboard_long.csv"].to_dict(orient="records"):
        task_scope = str(row["scope"])
        if task_scope != "cross":
            raise ValueError(f"S2 unexpected scope row: {row}")
        raw_representation = str(row["representation"])
        rows.append(
            {
                "task": "Task1",
                "task_scope": "cross",
                "analysis_family": metric_to_family(str(row["metric_name"])),
                "dataset": None,
                "cross_direction": str(row["dataset_or_direction"]),
                "cell_line": None,
                "perturbation_type": str(row["perturbation_type"]),
                "direction": None,
                "representation_raw": raw_representation,
                "representation_canonical": canonicalize_representation(raw_representation),
                "metric_name": str(row["metric_name"]),
                "metric_value": row["metric_value"],
                "row_origin": "observed_upstream",
                "rank_value": None,
                "rank_basis_metric_name": None,
                "leaderboard_eligible_bool": None,
                "cross_representation_comparable_bool": None,
                "n_total": row["n_total"],
                "n_valid": row["n_valid"],
                "n_excluded": row["n_excluded"],
                "n_excluded_missing_metric_or_mpos0": None,
                "n_targets_total": None,
                "n_targets_metric_valid": None,
                "N_gallery_mean": None,
                "N_gallery_max": row["N_gallery_max"],
                "m_pos_mean": None,
                "m_pos_p50": None,
                "m_pos_p90": None,
                "cross_alignment_contract": row["cross_alignment_contract"],
                "ingest_table": "task1_cross_leaderboard_long.csv",
                "ingest_stage": EXPECTED_S2_STAGE,
                "ingest_run_id": str(bundle.run_manifest["run_id"]),
                "source_table": "task1_cross_leaderboard_long.csv",
                "source_stage": EXPECTED_S2_STAGE,
                "source_run_id": str(bundle.run_manifest["run_id"]),
                "caution_codes": "",
            }
        )
    return rows


def build_s6_benchmark_rows(bundle: StageBundle) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for row in bundle.tables["task2_benchmark_summary_long.csv"].to_dict(orient="records"):
        raw_representation = str(row["representation"])
        rows.append(
            {
                "task": "Task2",
                "task_scope": "mechanism",
                "analysis_family": str(row["analysis_family"]),
                "dataset": str(row["dataset"]),
                "cross_direction": None,
                "cell_line": str(row["cell_line"]),
                "perturbation_type": None,
                "direction": str(row["direction"]) if str(row["direction"]).strip() else None,
                "representation_raw": raw_representation,
                "representation_canonical": canonicalize_representation(raw_representation),
                "metric_name": str(row["metric_name"]),
                "metric_value": row["metric_value"],
                "row_origin": "observed_upstream",
                "rank_value": row["rank_value"],
                "rank_basis_metric_name": row["rank_basis_metric_name"] if str(row["rank_basis_metric_name"]).strip() else None,
                "leaderboard_eligible_bool": row["leaderboard_eligible_bool"],
                "cross_representation_comparable_bool": row["cross_representation_comparable_bool"],
                "n_total": row["n_total"],
                "n_valid": row["n_valid"],
                "n_excluded": None,
                "n_excluded_missing_metric_or_mpos0": row["n_excluded_missing_metric_or_mpos0"],
                "n_targets_total": row["n_targets_total"],
                "n_targets_metric_valid": row["n_targets_metric_valid"],
                "N_gallery_mean": row["N_gallery_mean"],
                "N_gallery_max": row["N_gallery_max"],
                "m_pos_mean": row["m_pos_mean"],
                "m_pos_p50": row["m_pos_p50"],
                "m_pos_p90": row["m_pos_p90"],
                "cross_alignment_contract": None,
                "ingest_table": "task2_benchmark_summary_long.csv",
                "ingest_stage": EXPECTED_S6_STAGE,
                "ingest_run_id": str(bundle.run_manifest["run_id"]),
                "source_table": str(row["source_table"]),
                "source_stage": str(row["source_stage"]),
                "source_run_id": str(row["source_run_id"]),
                "caution_codes": str(row["caution_codes"]) if str(row["caution_codes"]).strip() else "",
            }
        )
    return rows


def dataset_sort_value(value: object, order_map: Mapping[str, int]) -> int:
    if value is None or pd.isna(value):
        return order_map["NA"]
    token = str(value)
    return order_map.get(token, len(order_map) + 1000)


def cell_line_sort_value(dataset: object, cell_line: object, cell_line_order_map: Mapping[str, Mapping[str, int]]) -> int:
    if cell_line is None or pd.isna(cell_line):
        return -1
    dataset_token = "NA" if dataset is None or pd.isna(dataset) else str(dataset)
    mapping = cell_line_order_map.get(dataset_token, {})
    return mapping.get(str(cell_line), len(mapping) + 1000)


def sort_project_benchmark_summary(
    frame: pd.DataFrame,
    representation_order: Mapping[str, int],
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    out = frame.copy()
    task_order = {value: idx for idx, value in enumerate(TASK_ORDER)}
    task_scope_order = {value: idx for idx, value in enumerate(TASK_SCOPE_ORDER)}
    analysis_family_order = {value: idx for idx, value in enumerate(ANALYSIS_FAMILY_ORDER)}
    dataset_order = {value: idx for idx, value in enumerate(DATASET_ORDER)}
    cross_direction_order = {value: idx for idx, value in enumerate(CROSS_DIRECTION_ORDER)}
    direction_order = {value: idx for idx, value in enumerate(DIRECTION_ORDER)}
    perturbation_order = {value: idx for idx, value in enumerate(PERTURBATION_TYPE_ORDER)}
    metric_order = {value: idx for idx, value in enumerate(METRIC_ORDER)}
    out["_task_order"] = out["task"].map(task_order)
    out["_task_scope_order"] = out["task_scope"].map(task_scope_order)
    out["_analysis_family_order"] = out["analysis_family"].map(analysis_family_order)
    out["_dataset_order"] = [dataset_sort_value(value, dataset_order) for value in out["dataset"].tolist()]
    out["_cross_direction_order"] = [dataset_sort_value(value, cross_direction_order) for value in out["cross_direction"].tolist()]
    out["_cell_line_order"] = [
        cell_line_sort_value(dataset, cell_line, cell_line_order_map)
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_perturbation_type_order"] = [dataset_sort_value(value, perturbation_order) for value in out["perturbation_type"].tolist()]
    out["_direction_order"] = [dataset_sort_value(value, direction_order) for value in out["direction"].tolist()]
    out["_metric_order"] = out["metric_name"].map(metric_order)
    out["_representation_order"] = out["representation_canonical"].map(representation_order)
    out = out.sort_values(
        [
            "_task_order",
            "_task_scope_order",
            "_analysis_family_order",
            "_dataset_order",
            "_cross_direction_order",
            "_cell_line_order",
            "_perturbation_type_order",
            "_direction_order",
            "_metric_order",
            "_representation_order",
            "representation_raw",
        ],
        kind="stable",
        na_position="last",
    )
    return out.drop(
        columns=[
            "_task_order",
            "_task_scope_order",
            "_analysis_family_order",
            "_dataset_order",
            "_cross_direction_order",
            "_cell_line_order",
            "_perturbation_type_order",
            "_direction_order",
            "_metric_order",
            "_representation_order",
        ]
    ).reset_index(drop=True)[PROJECT_BENCHMARK_SUMMARY_COLUMNS]


def build_project_benchmark_summary(
    s1_bundle: StageBundle,
    s2_bundle: StageBundle,
    s6_bundle: StageBundle,
    s1_has_soft: bool,
    representation_order: Mapping[str, int],
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    rows = build_s1_benchmark_rows(s1_bundle, s1_has_soft)
    rows.extend(build_s2_benchmark_rows(s2_bundle))
    rows.extend(build_s6_benchmark_rows(s6_bundle))
    out = pd.DataFrame(rows)
    assert_unique_key(
        out,
        ["task", "task_scope", "analysis_family", "dataset", "cross_direction", "cell_line", "perturbation_type", "direction", "representation_raw", "metric_name"],
        "project_benchmark_summary_long.csv",
    )
    out = sort_project_benchmark_summary(out, representation_order, cell_line_order_map)
    return out


def task1_family_id(row: Mapping[str, Any]) -> str:
    dataset_token = "NA" if row.get("dataset") is None or pd.isna(row.get("dataset")) else str(row.get("dataset"))
    cross_direction_token = "NA" if row.get("cross_direction") is None or pd.isna(row.get("cross_direction")) else str(row.get("cross_direction"))
    return (
        f"Task1|scope={row['task_scope']}|cross_direction={cross_direction_token}|dataset={dataset_token}"
        f"|perturbation_type={row['perturbation_type']}|representation_raw={row['representation_raw']}"
    )


def task2_family_id(row: Mapping[str, Any]) -> str:
    return (
        f"Task2|dataset={row['dataset']}|cell_line={row['cell_line']}"
        f"|direction={row['direction']}|representation_raw={row['representation_raw']}"
    )


def build_contract_exclusion_rows(s2_bundle: StageBundle, benchmark_summary: pd.DataFrame) -> List[Dict[str, Any]]:
    alignment = s2_bundle.tables["task1_cross_alignment_proof.csv"]
    chemical = alignment.loc[alignment["perturbation_type"].astype(str).eq("Chemical")].reset_index(drop=True)
    if len(chemical) != 1:
        raise ValueError("Expected exactly one Chemical row in task1_cross_alignment_proof.csv")
    chemical_row = chemical.iloc[0]
    if bool(chemical_row["eligible_bool"]) or str(chemical_row["excluded_reason"]) != "matched_keys_lt5":
        raise ValueError(f"Chemical cross exclusion contract drift detected: {chemical_row.to_dict()}")

    observed_cross = benchmark_summary.loc[
        benchmark_summary["task"].eq("Task1")
        & benchmark_summary["task_scope"].eq("cross")
        & benchmark_summary["analysis_family"].eq("retrieval")
        & benchmark_summary["metric_name"].eq("mean_mrr_corrected")
        & benchmark_summary["row_origin"].eq("observed_upstream")
    ].reset_index(drop=True)
    cross_directions = first_seen_values(observed_cross["cross_direction"].dropna().astype(str).tolist())
    representations = first_seen_values(observed_cross["representation_raw"].astype(str).tolist())
    if cross_directions != list(TASK1_CROSS_ALLOWED_DIRECTIONS):
        raise ValueError(f"Cross direction universe for exclusion materialization drift detected: {cross_directions}")
    if representations != list(S2_ALLOWED_REPRESENTATIONS):
        raise ValueError(f"Cross representation universe for exclusion materialization drift detected: {representations}")

    rows: List[Dict[str, Any]] = []
    for cross_direction in cross_directions:
        for raw_representation in representations:
            rows.append(
                {
                    "task_axis": SCORECARD_AXIS_BY_TASK["Task1"],
                    "family_id": (
                        f"Task1|scope=cross|cross_direction={cross_direction}|dataset=NA"
                        f"|perturbation_type=Chemical|representation_raw={raw_representation}"
                    ),
                    "row_origin": "contract_exclusion",
                    "task": "Task1",
                    "task_scope": "cross",
                    "dataset": None,
                    "cross_direction": cross_direction,
                    "cell_line": None,
                    "perturbation_type": "Chemical",
                    "direction": None,
                    "representation_raw": raw_representation,
                    "representation_canonical": canonicalize_representation(raw_representation),
                    "metric_name": "mean_mrr_corrected",
                    "metric_value": None,
                    "coverage_status": "contract_excluded",
                    "aggregation_weight": 0.0,
                    "cross_alignment_contract": str(chemical_row["cross_alignment_contract"]),
                    "ingest_table": "task1_cross_alignment_proof.csv",
                    "ingest_stage": EXPECTED_S2_STAGE,
                    "ingest_run_id": str(s2_bundle.run_manifest["run_id"]),
                    "source_table": "task1_cross_alignment_proof.csv",
                    "source_stage": EXPECTED_S2_STAGE,
                    "source_run_id": str(s2_bundle.run_manifest["run_id"]),
                    "caution_codes": TASK1_CROSS_EXCLUSION_CAUTION,
                }
            )
    return rows


def sort_project_axis_inputs(
    frame: pd.DataFrame,
    representation_order: Mapping[str, int],
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    out = frame.copy()
    task_order = {value: idx for idx, value in enumerate(TASK_ORDER)}
    task_scope_order = {value: idx for idx, value in enumerate(TASK_SCOPE_ORDER)}
    dataset_order = {value: idx for idx, value in enumerate(DATASET_ORDER)}
    cross_direction_order = {value: idx for idx, value in enumerate(CROSS_DIRECTION_ORDER)}
    direction_order = {value: idx for idx, value in enumerate(DIRECTION_ORDER)}
    perturbation_order = {value: idx for idx, value in enumerate(PERTURBATION_TYPE_ORDER)}
    row_origin_order = {value: idx for idx, value in enumerate(ROW_ORIGIN_ORDER)}
    out["_task_order"] = out["task"].map(task_order)
    out["_row_origin_order"] = out["row_origin"].map(row_origin_order)
    out["_task_scope_order"] = out["task_scope"].map(task_scope_order)
    out["_dataset_order"] = [dataset_sort_value(value, dataset_order) for value in out["dataset"].tolist()]
    out["_cross_direction_order"] = [dataset_sort_value(value, cross_direction_order) for value in out["cross_direction"].tolist()]
    out["_cell_line_order"] = [
        cell_line_sort_value(dataset, cell_line, cell_line_order_map)
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_perturbation_type_order"] = [dataset_sort_value(value, perturbation_order) for value in out["perturbation_type"].tolist()]
    out["_direction_order"] = [dataset_sort_value(value, direction_order) for value in out["direction"].tolist()]
    out["_representation_order"] = out["representation_canonical"].map(representation_order)
    out = out.sort_values(
        [
            "_task_order",
            "_row_origin_order",
            "_task_scope_order",
            "_dataset_order",
            "_cross_direction_order",
            "_cell_line_order",
            "_perturbation_type_order",
            "_direction_order",
            "_representation_order",
            "representation_raw",
        ],
        kind="stable",
        na_position="last",
    )
    return out.drop(
        columns=[
            "_task_order",
            "_row_origin_order",
            "_task_scope_order",
            "_dataset_order",
            "_cross_direction_order",
            "_cell_line_order",
            "_perturbation_type_order",
            "_direction_order",
            "_representation_order",
        ]
    ).reset_index(drop=True)[PROJECT_AXIS_COLUMNS]


def build_project_axis_inputs(
    benchmark_summary: pd.DataFrame,
    s2_bundle: StageBundle,
    representation_order: Mapping[str, int],
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    observed_rows: List[Dict[str, Any]] = []
    retrieval_mrr = benchmark_summary.loc[
        benchmark_summary["analysis_family"].eq("retrieval")
        & benchmark_summary["metric_name"].eq("mean_mrr_corrected")
        & benchmark_summary["row_origin"].eq("observed_upstream")
    ].reset_index(drop=True)
    for row in retrieval_mrr.to_dict(orient="records"):
        family_id = task1_family_id(row) if row["task"] == "Task1" else task2_family_id(row)
        observed_rows.append(
            {
                "task_axis": SCORECARD_AXIS_BY_TASK[str(row["task"])],
                "family_id": family_id,
                "row_origin": "observed_upstream",
                "task": row["task"],
                "task_scope": row["task_scope"],
                "dataset": row["dataset"],
                "cross_direction": row["cross_direction"],
                "cell_line": row["cell_line"],
                "perturbation_type": row["perturbation_type"],
                "direction": row["direction"],
                "representation_raw": row["representation_raw"],
                "representation_canonical": row["representation_canonical"],
                "metric_name": row["metric_name"],
                "metric_value": row["metric_value"],
                "coverage_status": "observed",
                "aggregation_weight": 1.0,
                "cross_alignment_contract": row["cross_alignment_contract"],
                "ingest_table": row["ingest_table"],
                "ingest_stage": row["ingest_stage"],
                "ingest_run_id": row["ingest_run_id"],
                "source_table": row["source_table"],
                "source_stage": row["source_stage"],
                "source_run_id": row["source_run_id"],
                "caution_codes": row["caution_codes"],
            }
        )
    rows = observed_rows + build_contract_exclusion_rows(s2_bundle, benchmark_summary)
    out = pd.DataFrame(rows)
    assert_unique_key(out, ["task_axis", "row_origin", "family_id"], "project_axis_score_inputs_long.csv")
    out = sort_project_axis_inputs(out, representation_order, cell_line_order_map)
    return out


def percentile_from_scores(series: pd.Series) -> pd.Series:
    out = pd.Series(pd.NA, index=series.index, dtype="Float64")
    valid = series.notna()
    if not bool(valid.any()):
        return out
    valid_series = pd.to_numeric(series.loc[valid], errors="coerce")
    if not bool(valid_series.notna().any()):
        return out
    valid_series = valid_series.dropna()
    if len(valid_series) == 1:
        out.loc[valid_series.index] = 1.0
        return out
    ranks = valid_series.rank(method="average", ascending=False)
    percentiles = 1.0 - (ranks - 1.0) / float(len(valid_series) - 1)
    out.loc[percentiles.index] = percentiles.astype(float)
    return out


def sort_project_scorecard(frame: pd.DataFrame, representation_order: Mapping[str, int]) -> pd.DataFrame:
    out = frame.copy()
    out["_representation_order"] = out["representation_canonical"].map(representation_order)
    out = out.sort_values(["_representation_order", "representation_canonical"], kind="stable", na_position="last")
    return out.drop(columns=["_representation_order"]).reset_index(drop=True)[PROJECT_SCORECARD_COLUMNS]


def build_project_scorecard(axis_inputs: pd.DataFrame, representation_order: Mapping[str, int]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for representation in first_seen_values(axis_inputs["representation_canonical"].astype(str).tolist()):
        subset = axis_inputs.loc[axis_inputs["representation_canonical"].astype(str).eq(representation)].reset_index(drop=True)
        task1_rows = subset.loc[subset["task_axis"].eq("Task1")].reset_index(drop=True)
        task2_rows = subset.loc[subset["task_axis"].eq("Task2")].reset_index(drop=True)
        task1_observed = task1_rows.loc[task1_rows["row_origin"].eq("observed_upstream")].reset_index(drop=True)
        task1_excluded = task1_rows.loc[task1_rows["row_origin"].eq("contract_exclusion")].reset_index(drop=True)
        task2_observed = task2_rows.loc[task2_rows["row_origin"].eq("observed_upstream")].reset_index(drop=True)
        task2_excluded = task2_rows.loc[task2_rows["row_origin"].eq("contract_exclusion")].reset_index(drop=True)

        task1_axis_value_raw = (
            float(pd.to_numeric(task1_observed["metric_value"], errors="coerce").mean())
            if not task1_observed.empty
            else None
        )
        task2_axis_value_raw = (
            float(pd.to_numeric(task2_observed["metric_value"], errors="coerce").mean())
            if not task2_observed.empty
            else None
        )
        task1_observed_families = int(len(task1_observed))
        task1_excluded_families = int(len(task1_excluded))
        task2_observed_families = int(len(task2_observed))
        task2_excluded_families = int(len(task2_excluded))
        eligible = task1_observed_families > 0 and task2_observed_families > 0
        caution_union: List[str] = []
        for raw in subset["caution_codes"].tolist():
            caution_union.extend(parse_caution_codes(raw))

        rows.append(
            {
                "representation_canonical": representation,
                "task1_axis_value_raw": task1_axis_value_raw,
                "task1_axis_observed_families": task1_observed_families,
                "task1_axis_contract_exclusion_families": task1_excluded_families,
                "task1_axis_expected_families": task1_observed_families + task1_excluded_families,
                "task2_axis_value_raw": task2_axis_value_raw,
                "task2_axis_observed_families": task2_observed_families,
                "task2_axis_contract_exclusion_families": task2_excluded_families,
                "task2_axis_expected_families": task2_observed_families + task2_excluded_families,
                "x_percentile": None,
                "y_percentile": None,
                "scorecard_eligible_bool": eligible,
                "caution_codes": json_array_string(caution_union),
            }
        )

    out = pd.DataFrame(rows)
    eligible_mask = out["scorecard_eligible_bool"].astype(bool)
    eligible_with_x = eligible_mask & out["task1_axis_value_raw"].notna()
    eligible_with_y = eligible_mask & out["task2_axis_value_raw"].notna()
    out.loc[:, "x_percentile"] = percentile_from_scores(out.loc[eligible_with_x, "task1_axis_value_raw"]).reindex(out.index)
    out.loc[:, "y_percentile"] = percentile_from_scores(out.loc[eligible_with_y, "task2_axis_value_raw"]).reindex(out.index)
    out.loc[~eligible_mask, "x_percentile"] = pd.NA
    out.loc[~eligible_mask, "y_percentile"] = pd.NA
    out.loc[out["task1_axis_value_raw"].isna(), "x_percentile"] = pd.NA
    out.loc[out["task2_axis_value_raw"].isna(), "y_percentile"] = pd.NA
    assert_unique_key(out, ["representation_canonical"], "project_representation_scorecard.csv")
    out = sort_project_scorecard(out, representation_order)
    return out


def build_input_registry(
    s1_bundle: StageBundle,
    s2_bundle: StageBundle,
    s4_bundle: StageBundle,
    s5_bundle: StageBundle,
    s6_bundle: StageBundle,
    s1_non_pass_names: Sequence[str],
) -> pd.DataFrame:
    rows = [
        {
            "provenance_role": "direct_ingest_stage",
            "task": "Task1",
            "stage_name": EXPECTED_S1_STAGE,
            "run_id": str(s1_bundle.run_manifest["run_id"]),
            "stage_dir": str(s1_bundle.stage_dir),
            "reporting_tables_json": json.dumps(["task1_leaderboard_long.csv"], ensure_ascii=True),
            "validation_tables_json": json.dumps(
                ["task1_retrieval_summary.csv", "task1_chance_identity_check.csv", "task1_attrition.csv"], ensure_ascii=True
            ),
            "manifest_path": str(s1_bundle.paths["run_manifest.json"]),
            "audit_assertions_path": str(s1_bundle.paths["audit_assertions.json"]),
            "stage_manifest_path": str(s1_bundle.paths["manifest.json"]),
            "audit_status": "pass_with_whitelisted_nonpass" if s1_non_pass_names else "pass",
            "non_pass_assertions_json": json.dumps(list(s1_non_pass_names), ensure_ascii=True),
            "notes": "whitelisted_soft_assertion" if s1_non_pass_names else "",
        },
        {
            "provenance_role": "direct_ingest_stage",
            "task": "Task1",
            "stage_name": EXPECTED_S2_STAGE,
            "run_id": str(s2_bundle.run_manifest["run_id"]),
            "stage_dir": str(s2_bundle.stage_dir),
            "reporting_tables_json": json.dumps(["task1_cross_leaderboard_long.csv"], ensure_ascii=True),
            "validation_tables_json": json.dumps(
                [
                    "task1_cross_retrieval_summary.csv",
                    "task1_cross_chance_identity_check.csv",
                    "task1_cross_alignment_proof.csv",
                    "task1_cross_attrition.csv",
                ],
                ensure_ascii=True,
            ),
            "manifest_path": str(s2_bundle.paths["run_manifest.json"]),
            "audit_assertions_path": str(s2_bundle.paths["audit_assertions.json"]),
            "stage_manifest_path": str(s2_bundle.paths["manifest.json"]),
            "audit_status": "pass",
            "non_pass_assertions_json": json.dumps([], ensure_ascii=True),
            "notes": "",
        },
        {
            "provenance_role": "direct_ingest_stage",
            "task": "Task2",
            "stage_name": EXPECTED_S6_STAGE,
            "run_id": str(s6_bundle.run_manifest["run_id"]),
            "stage_dir": str(s6_bundle.stage_dir),
            "reporting_tables_json": json.dumps(["task2_benchmark_summary_long.csv"], ensure_ascii=True),
            "validation_tables_json": json.dumps(
                ["task2_group_leaderboard.csv", "task2_retrieval_leaderboard.csv", "task2_group_concordance_long.csv"],
                ensure_ascii=True,
            ),
            "manifest_path": str(s6_bundle.paths["run_manifest.json"]),
            "audit_assertions_path": str(s6_bundle.paths["audit_assertions.json"]),
            "stage_manifest_path": str(s6_bundle.paths["manifest.json"]),
            "audit_status": "pass",
            "non_pass_assertions_json": json.dumps([], ensure_ascii=True),
            "notes": "",
        },
        {
            "provenance_role": "transitive_provenance_stage",
            "task": "Task2",
            "stage_name": EXPECTED_S4_STAGE,
            "run_id": str(s4_bundle.run_manifest["run_id"]),
            "stage_dir": str(s4_bundle.stage_dir),
            "reporting_tables_json": json.dumps([], ensure_ascii=True),
            "validation_tables_json": json.dumps([], ensure_ascii=True),
            "manifest_path": str(s4_bundle.paths["run_manifest.json"]),
            "audit_assertions_path": str(s4_bundle.paths["audit_assertions.json"]),
            "stage_manifest_path": str(s4_bundle.paths["manifest.json"]),
            "audit_status": "pass",
            "non_pass_assertions_json": json.dumps([], ensure_ascii=True),
            "notes": "transitive_provenance_only_via_s6_source_fields",
        },
        {
            "provenance_role": "transitive_provenance_stage",
            "task": "Task2",
            "stage_name": EXPECTED_S5_STAGE,
            "run_id": str(s5_bundle.run_manifest["run_id"]),
            "stage_dir": str(s5_bundle.stage_dir),
            "reporting_tables_json": json.dumps([], ensure_ascii=True),
            "validation_tables_json": json.dumps([], ensure_ascii=True),
            "manifest_path": str(s5_bundle.paths["run_manifest.json"]),
            "audit_assertions_path": str(s5_bundle.paths["audit_assertions.json"]),
            "stage_manifest_path": str(s5_bundle.paths["manifest.json"]),
            "audit_status": "pass",
            "non_pass_assertions_json": json.dumps([], ensure_ascii=True),
            "notes": "transitive_provenance_only_via_s6_source_fields",
        },
    ]
    out = pd.DataFrame(rows)
    assert_unique_key(out, ["provenance_role", "task", "stage_name", "run_id"], "project_input_registry.csv")
    provenance_role_order = {value: idx for idx, value in enumerate(PROVENANCE_ROLE_ORDER)}
    task_order = {value: idx for idx, value in enumerate(TASK_ORDER)}
    stage_order = {value: idx for idx, value in enumerate(STAGE_ORDER)}
    out["_provenance_role_order"] = out["provenance_role"].map(provenance_role_order)
    out["_task_order"] = out["task"].map(task_order)
    out["_stage_order"] = out["stage_name"].map(stage_order)
    out = out.sort_values(
        ["_provenance_role_order", "_task_order", "_stage_order", "run_id"],
        kind="stable",
        na_position="last",
    )
    return out.drop(columns=["_provenance_role_order", "_task_order", "_stage_order"]).reset_index(drop=True)[
        PROJECT_INPUT_REGISTRY_COLUMNS
    ]


def frame_is_sorted(frame: pd.DataFrame, sorted_frame: pd.DataFrame) -> bool:
    if list(frame.columns) != list(sorted_frame.columns):
        return False
    if len(frame) != len(sorted_frame):
        return False
    return frame.reset_index(drop=True).equals(sorted_frame.reset_index(drop=True))


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


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    seed, runs_dir = read_config(project_root)
    if args.seed is not None and int(args.seed) != seed:
        raise ValueError(f"Seed mismatch: cli={args.seed} config={seed}")

    stage_dir = (runs_dir / args.run_id / STAGE).resolve()
    if stage_dir.exists():
        raise FileExistsError(f"Output stage directory already exists: {stage_dir}")
    stage_dir.mkdir(parents=True, exist_ok=False)

    s1_bundle = load_stage_bundle(
        args.task1_s1_stage_dir,
        S1_REQUIRED_FILES,
        ["task1_leaderboard_long.csv", "task1_retrieval_summary.csv", "task1_chance_identity_check.csv", "task1_attrition.csv"],
    )
    s2_bundle = load_stage_bundle(
        args.task1_s2_stage_dir,
        S2_REQUIRED_FILES,
        [
            "task1_cross_leaderboard_long.csv",
            "task1_cross_retrieval_summary.csv",
            "task1_cross_chance_identity_check.csv",
            "task1_cross_alignment_proof.csv",
            "task1_cross_attrition.csv",
        ],
    )
    s6_bundle = load_stage_bundle(
        args.task2_s6_stage_dir,
        S6_REQUIRED_FILES,
        [
            "task2_benchmark_summary_long.csv",
            "task2_group_leaderboard.csv",
            "task2_retrieval_leaderboard.csv",
            "task2_group_concordance_long.csv",
        ],
    )

    s1_has_soft = validate_s1_bundle(s1_bundle)
    s2_bundle_validation = validate_s2_bundle(s2_bundle)
    _ = s2_bundle_validation
    s4_stage_dir, s5_stage_dir, cell_line_order_by_dataset = validate_s6_bundle(s6_bundle)

    s4_bundle = load_stage_bundle(s4_stage_dir, S4_REQUIRED_FILES, [])
    s5_bundle = load_stage_bundle(s5_stage_dir, S5_REQUIRED_FILES, [])
    validate_s4_transitive_bundle(s4_bundle)
    validate_s5_transitive_bundle(s5_bundle)

    s1_non_pass_names, _ = validate_s1_assertions(s1_bundle)
    input_registry = build_input_registry(s1_bundle, s2_bundle, s4_bundle, s5_bundle, s6_bundle, s1_non_pass_names)
    representation_canonical_order = build_representation_canonical_order(
        s1_bundle.tables["task1_leaderboard_long.csv"],
        s2_bundle.tables["task1_cross_leaderboard_long.csv"],
        s6_bundle.tables["task2_benchmark_summary_long.csv"],
    )
    representation_order_map = {value: idx for idx, value in enumerate(representation_canonical_order)}
    cell_line_order_map = build_cell_line_order_map(cell_line_order_by_dataset)

    benchmark_summary = build_project_benchmark_summary(
        s1_bundle,
        s2_bundle,
        s6_bundle,
        s1_has_soft,
        representation_order_map,
        cell_line_order_map,
    )
    axis_inputs = build_project_axis_inputs(benchmark_summary, s2_bundle, representation_order_map, cell_line_order_map)
    scorecard = build_project_scorecard(axis_inputs, representation_order_map)

    expected_registry_rows = 5
    expected_benchmark_rows = int(
        len(s1_bundle.tables["task1_leaderboard_long.csv"])
        + len(s2_bundle.tables["task1_cross_leaderboard_long.csv"])
        + len(s6_bundle.tables["task2_benchmark_summary_long.csv"])
    )
    expected_axis_rows = int(
        len(
            benchmark_summary.loc[
                benchmark_summary["analysis_family"].eq("retrieval")
                & benchmark_summary["metric_name"].eq("mean_mrr_corrected")
                & benchmark_summary["row_origin"].eq("observed_upstream")
            ]
        )
        + 4
    )
    expected_scorecard_rows = int(axis_inputs["representation_canonical"].astype(str).nunique())

    if len(input_registry) != expected_registry_rows:
        raise ValueError(f"project_input_registry.csv row count mismatch: {len(input_registry)}")
    if len(benchmark_summary) != expected_benchmark_rows:
        raise ValueError(f"project_benchmark_summary_long.csv row count mismatch: {len(benchmark_summary)}")
    if len(axis_inputs) != expected_axis_rows:
        raise ValueError(f"project_axis_score_inputs_long.csv row count mismatch: {len(axis_inputs)}")
    if len(scorecard) != expected_scorecard_rows:
        raise ValueError(f"project_representation_scorecard.csv row count mismatch: {len(scorecard)}")

    benchmark_path = stage_dir / "project_benchmark_summary_long.csv"
    axis_path = stage_dir / "project_axis_score_inputs_long.csv"
    scorecard_path = stage_dir / "project_representation_scorecard.csv"
    registry_path = stage_dir / "project_input_registry.csv"
    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    manifest_path = stage_dir / "manifest.json"

    write_csv(input_registry, registry_path)
    write_csv(benchmark_summary, benchmark_path)
    write_csv(axis_inputs, axis_path)
    write_csv(scorecard, scorecard_path)

    input_paths = [
        str(path)
        for bundle in (s1_bundle, s2_bundle, s6_bundle, s4_bundle, s5_bundle)
        for path in bundle.paths.values()
    ]
    output_paths = [
        str(registry_path),
        str(benchmark_path),
        str(axis_path),
        str(scorecard_path),
        str(run_manifest_path),
        str(audit_assertions_path),
        str(manifest_path),
    ]

    benchmark_row_origin_ok = bool(benchmark_summary["row_origin"].eq("observed_upstream").all())
    axis_row_origin_ok = set(axis_inputs["row_origin"].astype(str).unique()) == {"observed_upstream", "contract_exclusion"}
    axis_coverage_ok = bool(axis_inputs["coverage_status"].isin(AXIS_ALLOWED_COVERAGE_STATUS).all()) and bool(
        (
            (
                axis_inputs["row_origin"].eq("observed_upstream")
                & axis_inputs["coverage_status"].eq("observed")
            )
            | (
                axis_inputs["row_origin"].eq("contract_exclusion")
                & axis_inputs["coverage_status"].eq("contract_excluded")
            )
        ).all()
    )
    scorecard_eligibility_expected = (
        (pd.to_numeric(scorecard["task1_axis_observed_families"], errors="coerce") > 0)
        & (pd.to_numeric(scorecard["task2_axis_observed_families"], errors="coerce") > 0)
    )
    scorecard_eligibility_ok = bool(
        scorecard["scorecard_eligible_bool"].astype(bool).reset_index(drop=True).equals(scorecard_eligibility_expected.reset_index(drop=True))
    )
    eligible_mask = scorecard["scorecard_eligible_bool"].astype(bool)
    percentile_universe_ok = bool(scorecard.loc[~eligible_mask, "x_percentile"].isna().all()) and bool(
        scorecard.loc[~eligible_mask, "y_percentile"].isna().all()
    )
    percentile_universe_ok = percentile_universe_ok and bool(scorecard.loc[scorecard["task1_axis_value_raw"].isna(), "x_percentile"].isna().all())
    percentile_universe_ok = percentile_universe_ok and bool(scorecard.loc[scorecard["task2_axis_value_raw"].isna(), "y_percentile"].isna().all())
    null_sentinel_leak_ok = True
    for frame in (benchmark_summary, axis_inputs, scorecard):
        stringified = frame.astype(str)
        if bool(stringified.eq("__NULL__").any().any()):
            null_sentinel_leak_ok = False
            break
    s1_rows_have_caution = bool(
        benchmark_summary.loc[benchmark_summary["ingest_stage"].eq(EXPECTED_S1_STAGE), "caution_codes"].eq(S1_SOFT_CAUTION).all()
    )
    s1_axis_rows_have_caution = bool(
        axis_inputs.loc[axis_inputs["ingest_stage"].eq(EXPECTED_S1_STAGE), "caution_codes"].eq(S1_SOFT_CAUTION).all()
    )
    scorecard_caution_ok = True
    for row in scorecard.to_dict(orient="records"):
        representation = str(row["representation_canonical"])
        subset = axis_inputs.loc[axis_inputs["representation_canonical"].astype(str).eq(representation)].reset_index(drop=True)
        expected_json = json_array_string(code for raw in subset["caution_codes"].tolist() for code in parse_caution_codes(raw))
        if str(row["caution_codes"]) != expected_json:
            scorecard_caution_ok = False
            break

    benchmark_sorted = sort_project_benchmark_summary(benchmark_summary, representation_order_map, cell_line_order_map)
    axis_sorted = sort_project_axis_inputs(axis_inputs, representation_order_map, cell_line_order_map)
    scorecard_sorted = sort_project_scorecard(scorecard, representation_order_map)
    input_registry_sorted = build_input_registry(s1_bundle, s2_bundle, s4_bundle, s5_bundle, s6_bundle, s1_non_pass_names)
    sort_ok = all(
        [
            frame_is_sorted(input_registry, input_registry_sorted),
            frame_is_sorted(benchmark_summary, benchmark_sorted),
            frame_is_sorted(axis_inputs, axis_sorted),
            frame_is_sorted(scorecard, scorecard_sorted),
        ]
    )

    task1_harmonization_ok = True
    task1_rows = benchmark_summary.loc[benchmark_summary["task"].eq("Task1")].reset_index(drop=True)
    if not task1_rows.empty:
        task1_harmonization_ok = task1_harmonization_ok and bool(
            (
                task1_rows["analysis_family"].eq("group_concordance")
                == task1_rows["metric_name"].isin(TASK1_GROUP_METRICS)
            ).all()
        )
        task1_harmonization_ok = task1_harmonization_ok and bool(
            (
                task1_rows["analysis_family"].eq("retrieval")
                == task1_rows["metric_name"].isin(TASK1_RETRIEVAL_METRICS)
            ).all()
        )
        task1_harmonization_ok = task1_harmonization_ok and bool(
            task1_rows.loc[task1_rows["task_scope"].eq("internal"), "dataset"].notna().all()
        )
        task1_harmonization_ok = task1_harmonization_ok and bool(
            task1_rows.loc[task1_rows["task_scope"].eq("cross"), "cross_direction"].notna().all()
        )
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["direction"].isna().all())
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["cell_line"].isna().all())
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["rank_value"].isna().all())
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["leaderboard_eligible_bool"].isna().all())
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["cross_representation_comparable_bool"].isna().all())
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["n_excluded_missing_metric_or_mpos0"].isna().all())
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["n_targets_total"].isna().all())
        task1_harmonization_ok = task1_harmonization_ok and bool(task1_rows["n_targets_metric_valid"].isna().all())

    task2_harmonization_ok = True
    task2_rows = benchmark_summary.loc[benchmark_summary["task"].eq("Task2")].reset_index(drop=True)
    if not task2_rows.empty:
        task2_harmonization_ok = task2_harmonization_ok and bool(task2_rows["task_scope"].eq("mechanism").all())
        task2_harmonization_ok = task2_harmonization_ok and bool(task2_rows["dataset"].notna().all())
        task2_harmonization_ok = task2_harmonization_ok and bool(task2_rows["cross_direction"].isna().all())
        task2_harmonization_ok = task2_harmonization_ok and bool(task2_rows["perturbation_type"].isna().all())
        task2_harmonization_ok = task2_harmonization_ok and bool(task2_rows["n_excluded"].isna().all())
    harmonization_ok = task1_harmonization_ok and task2_harmonization_ok

    legacy_contamination_ok = bool(benchmark_summary["source_stage"].astype(str).isin({EXPECTED_S1_STAGE, EXPECTED_S2_STAGE, EXPECTED_S4_STAGE, EXPECTED_S5_STAGE}).all())
    legacy_contamination_ok = legacy_contamination_ok and bool(input_registry["stage_name"].astype(str).isin(STAGE_ORDER).all())
    edist_not_in_axis_ok = not bool(axis_inputs["metric_name"].eq("mean_edist_biascorr").any())
    scorecard_source_ok = set(scorecard["representation_canonical"].astype(str).tolist()) == set(
        axis_inputs["representation_canonical"].astype(str).unique().tolist()
    )
    no_hidden_recompute_ok = all("per_query" not in path for path in input_paths) and all("/data/task" not in path for path in input_paths)
    source_provenance_ok = bool(
        benchmark_summary[["ingest_table", "ingest_stage", "ingest_run_id", "source_table", "source_stage", "source_run_id"]].notna().all().all()
    )
    not_applicable_vs_exclusion_ok = bool(
        axis_inputs.loc[axis_inputs["row_origin"].eq("contract_exclusion"), "coverage_status"].eq("contract_excluded").all()
    )
    representation_canonicalization_ok = bool(
        benchmark_summary["representation_canonical"].astype(str).reset_index(drop=True).equals(
            benchmark_summary["representation_raw"].astype(str).map(canonicalize_representation).reset_index(drop=True)
        )
    )
    representation_canonicalization_ok = representation_canonicalization_ok and {"tahoex1_3b", "tahoe-x1"}.issubset(
        set(scorecard["representation_canonical"].astype(str).tolist())
    )

    assertions = [
        {
            "name": "input_stage_contract_locked",
            "pass": bool(
                input_registry["stage_name"].astype(str).tolist()
                == [EXPECTED_S1_STAGE, EXPECTED_S2_STAGE, EXPECTED_S6_STAGE, EXPECTED_S4_STAGE, EXPECTED_S5_STAGE]
            ),
            "details": {
                "expected_stage_order": list(STAGE_ORDER),
                "actual_stage_order": input_registry["stage_name"].astype(str).tolist(),
            },
            "counterexamples": [],
        },
        {
            "name": "whitelisted_soft_assertions_only",
            "pass": s1_non_pass_names in ([], [S1_SOFT_ASSERTION]),
            "details": {"s1_non_pass_assertions": list(s1_non_pass_names)},
            "counterexamples": [] if s1_non_pass_names in ([], [S1_SOFT_ASSERTION]) else [{"s1_non_pass_assertions": s1_non_pass_names}],
        },
        {
            "name": "task2_successor_path_locked",
            "pass": bool(
                s6_bundle.run_manifest.get("config", {}).get("s4_source_stage") == EXPECTED_S4_STAGE
                and s6_bundle.run_manifest.get("config", {}).get("s5_source_stage") == EXPECTED_S5_STAGE
            ),
            "details": {
                "s4_stage_dir": str(s4_stage_dir),
                "s5_stage_dir": str(s5_stage_dir),
            },
            "counterexamples": [],
        },
        {
            "name": "task1_task2_reporting_schema_harmonization_contract",
            "pass": bool(harmonization_ok),
            "details": {
                "task1_rows": int(len(task1_rows)),
                "task2_rows": int(len(task2_rows)),
            },
            "counterexamples": [] if harmonization_ok else benchmark_summary.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "observed_vs_contract_exclusion_row_separation",
            "pass": bool(benchmark_row_origin_ok and axis_row_origin_ok),
            "details": {
                "benchmark_row_origin_values": sorted(benchmark_summary["row_origin"].astype(str).unique().tolist()),
                "axis_row_origin_values": sorted(axis_inputs["row_origin"].astype(str).unique().tolist()),
            },
            "counterexamples": []
            if benchmark_row_origin_ok and axis_row_origin_ok
            else {
                "benchmark_bad_rows": benchmark_summary.loc[benchmark_summary["row_origin"].ne("observed_upstream")].head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
                "axis_bad_rows": axis_inputs.loc[~axis_inputs["row_origin"].isin(["observed_upstream", "contract_exclusion"])].head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
            },
        },
        {
            "name": "project_summary_union_contract",
            "pass": bool(len(benchmark_summary) == expected_benchmark_rows),
            "details": {
                "expected_rows": expected_benchmark_rows,
                "actual_rows": int(len(benchmark_summary)),
                "s1_rows": int(len(s1_bundle.tables["task1_leaderboard_long.csv"])),
                "s2_rows": int(len(s2_bundle.tables["task1_cross_leaderboard_long.csv"])),
                "s6_rows": int(len(s6_bundle.tables["task2_benchmark_summary_long.csv"])),
            },
            "counterexamples": [] if len(benchmark_summary) == expected_benchmark_rows else benchmark_summary.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "scorecard_axis_contract",
            "pass": bool(
                len(axis_inputs) == expected_axis_rows
                and set(axis_inputs["metric_name"].astype(str).unique()) == {"mean_mrr_corrected"}
                and edist_not_in_axis_ok
                and scorecard_source_ok
            ),
            "details": {
                "expected_rows": expected_axis_rows,
                "actual_rows": int(len(axis_inputs)),
                "row_origin_values": sorted(axis_inputs["row_origin"].astype(str).unique().tolist()),
                "coverage_status_values": sorted(axis_inputs["coverage_status"].astype(str).unique().tolist()),
            },
            "counterexamples": []
            if len(axis_inputs) == expected_axis_rows and set(axis_inputs["metric_name"].astype(str).unique()) == {"mean_mrr_corrected"} and edist_not_in_axis_ok and scorecard_source_ok
            else axis_inputs.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "representation_canonicalization_contract",
            "pass": bool(representation_canonicalization_ok),
            "details": {
                "representation_canonical_order": representation_canonical_order,
            },
            "counterexamples": [] if representation_canonicalization_ok else benchmark_summary.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "output_unique_key_and_sort_contract",
            "pass": bool(sort_ok),
            "details": {
                "registry_rows": int(len(input_registry)),
                "benchmark_rows": int(len(benchmark_summary)),
                "axis_rows": int(len(axis_inputs)),
                "scorecard_rows": int(len(scorecard)),
            },
            "counterexamples": []
            if sort_ok
            else {
                "registry_sorted": frame_is_sorted(input_registry, input_registry_sorted),
                "benchmark_sorted": frame_is_sorted(benchmark_summary, benchmark_sorted),
                "axis_sorted": frame_is_sorted(axis_inputs, axis_sorted),
                "scorecard_sorted": frame_is_sorted(scorecard, scorecard_sorted),
            },
        },
        {
            "name": "scorecard_eligibility_contract",
            "pass": bool(scorecard_eligibility_ok),
            "details": {
                "eligible_rows": int(eligible_mask.sum()),
                "total_rows": int(len(scorecard)),
            },
            "counterexamples": []
            if scorecard_eligibility_ok
            else scorecard.loc[scorecard["scorecard_eligible_bool"].astype(bool) != scorecard_eligibility_expected].head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "percentile_universe_contract",
            "pass": bool(percentile_universe_ok),
            "details": {
                "tie_rule": PERCENTILE_TIE_RULE,
                "eligible_rows": int(eligible_mask.sum()),
            },
            "counterexamples": []
            if percentile_universe_ok
            else scorecard.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "axis_coverage_status_contract",
            "pass": bool(axis_coverage_ok),
            "details": {
                "allowed_values": sorted(AXIS_ALLOWED_COVERAGE_STATUS),
                "actual_values": sorted(axis_inputs["coverage_status"].astype(str).unique().tolist()),
            },
            "counterexamples": []
            if axis_coverage_ok
            else axis_inputs.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "null_vs_sort_sentinel_contract",
            "pass": bool(null_sentinel_leak_ok),
            "details": {
                "sentinel": "__NULL__",
            },
            "counterexamples": [] if null_sentinel_leak_ok else [{"sentinel_found": "__NULL__"}],
        },
        {
            "name": "s1_soft_audit_caution_propagation_contract",
            "pass": bool((not s1_has_soft) or (s1_rows_have_caution and s1_axis_rows_have_caution)),
            "details": {
                "s1_has_soft_assertion": bool(s1_has_soft),
                "benchmark_rows_with_s1_caution": bool(s1_rows_have_caution),
                "axis_rows_with_s1_caution": bool(s1_axis_rows_have_caution),
            },
            "counterexamples": []
            if (not s1_has_soft) or (s1_rows_have_caution and s1_axis_rows_have_caution)
            else benchmark_summary.loc[benchmark_summary["ingest_stage"].eq(EXPECTED_S1_STAGE)].head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "scorecard_caution_aggregation_contract",
            "pass": bool(scorecard_caution_ok),
            "details": {
                "serialization": "json_array_string",
            },
            "counterexamples": [] if scorecard_caution_ok else scorecard.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "no_hidden_recomputation",
            "pass": bool(no_hidden_recompute_ok),
            "details": {
                "n_inputs": len(input_paths),
                "reads_only_from_runs": True,
            },
            "counterexamples": [] if no_hidden_recompute_ok else [{"input_paths": input_paths[:MAX_COUNTEREXAMPLES]}],
        },
        {
            "name": "source_provenance_fields_complete",
            "pass": bool(source_provenance_ok),
            "details": {
                "required_fields": [
                    "ingest_table",
                    "ingest_stage",
                    "ingest_run_id",
                    "source_table",
                    "source_stage",
                    "source_run_id",
                ],
            },
            "counterexamples": [] if source_provenance_ok else benchmark_summary.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "legacy_contamination_absent",
            "pass": bool(legacy_contamination_ok),
            "details": {
                "source_stages": sorted(benchmark_summary["source_stage"].astype(str).unique().tolist()),
            },
            "counterexamples": [] if legacy_contamination_ok else benchmark_summary.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "task2_edist_informational_only_preserved",
            "pass": bool(edist_not_in_axis_ok),
            "details": {
                "axis_metric_names": sorted(axis_inputs["metric_name"].astype(str).unique().tolist()),
            },
            "counterexamples": [] if edist_not_in_axis_ok else axis_inputs.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
        {
            "name": "not_applicable_vs_exclusion_semantics_preserved",
            "pass": bool(not_applicable_vs_exclusion_ok),
            "details": {
                "contract_exclusion_rows": int(axis_inputs["row_origin"].eq("contract_exclusion").sum()),
            },
            "counterexamples": [] if not_applicable_vs_exclusion_ok else axis_inputs.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        },
    ]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts" / "s7_project_benchmark_synthesis.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": utc_now_iso(),
        "completed_at": utc_now_iso(),
        "config": {
            "seed": seed,
            "runs_dir": str(runs_dir),
            "task1_s1_stage_dir": str(s1_bundle.stage_dir),
            "task1_s2_stage_dir": str(s2_bundle.stage_dir),
            "task2_s6_stage_dir": str(s6_bundle.stage_dir),
            "task2_s4_stage_dir_transitive": str(s4_bundle.stage_dir),
            "task2_s5_stage_dir_transitive": str(s5_bundle.stage_dir),
            "task_order": list(TASK_ORDER),
            "task_scope_order": list(TASK_SCOPE_ORDER),
            "analysis_family_order": list(ANALYSIS_FAMILY_ORDER),
            "dataset_order": list(DATASET_ORDER),
            "cross_direction_order": list(CROSS_DIRECTION_ORDER),
            "direction_order": list(DIRECTION_ORDER),
            "perturbation_type_order": list(PERTURBATION_TYPE_ORDER),
            "metric_order": list(METRIC_ORDER),
            "representation_canonical_order": list(representation_canonical_order),
            "percentile_tie_rule": PERCENTILE_TIE_RULE,
            "frozen_policy_defaults_requiring_explicit_human_acceptance": [
                "equal_family_macro_average_over_observed_mean_mrr_corrected",
                "prefix_strip_only_representation_canonicalization",
            ],
            "whitelisted_soft_assertion": S1_SOFT_ASSERTION,
            "whitelisted_soft_assertion_present": bool(s1_has_soft),
        },
        "inputs": input_paths,
        "outputs": output_paths,
        "summary": {
            "expected_project_input_registry_rows": expected_registry_rows,
            "n_project_input_registry_rows": int(len(input_registry)),
            "expected_project_benchmark_summary_rows": expected_benchmark_rows,
            "n_project_benchmark_summary_rows": int(len(benchmark_summary)),
            "expected_project_axis_score_inputs_rows": expected_axis_rows,
            "n_project_axis_score_inputs_rows": int(len(axis_inputs)),
            "expected_project_representation_scorecard_rows": expected_scorecard_rows,
            "n_project_representation_scorecard_rows": int(len(scorecard)),
            "representation_canonical_order": list(representation_canonical_order),
            "task1_axis_observed_rows": int(axis_inputs["task_axis"].eq("Task1").sum()),
            "task2_axis_observed_rows": int(axis_inputs["task_axis"].eq("Task2").sum()),
            "contract_exclusion_rows": int(axis_inputs["row_origin"].eq("contract_exclusion").sum()),
            "scorecard_eligible_rows": int(scorecard["scorecard_eligible_bool"].astype(bool).sum()),
        },
    }

    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})
    write_json(manifest_path, {"stage": STAGE, "files": build_stage_manifest(stage_dir)})


if __name__ == "__main__":
    main()
