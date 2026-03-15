# SCRIPT_HEADER_CONTRACT
# Script: scripts/s6_task2_result_synthesis.py
# Legacy note:
#   - Preserved historical scPerturb-K562 S6 artifact.
#   - Corrected multisource Task2 S6 synthesis now lives in scripts/s6_task2_result_synthesis_multisource.py.
# Purpose: Synthesize legacy/interim audited Task2 S4/S5 outputs into benchmark-ready summary tables without recomputing metrics.
# Inputs:
#   - S4 audited stage dir with:
#     - task2_group_concordance.csv
#     - task2_group_attrition.csv
#     - run_manifest.json
#     - audit_assertions.json
#     - manifest.json
#   - S5 audited stage dir with:
#     - task2_retrieval_summary.csv
#     - task2_retrieval_summary_long.csv
#     - task2_retrieval_attrition.csv
#     - task2_chance_identity_check.csv
#     - run_manifest.json
#     - audit_assertions.json
#     - manifest.json
# Outputs:
#   - task2_group_concordance_long.csv: runs/<run_id>/s6_task2_result_synthesis/
#   - task2_group_leaderboard.csv: runs/<run_id>/s6_task2_result_synthesis/
#   - task2_retrieval_leaderboard.csv: runs/<run_id>/s6_task2_result_synthesis/
#   - task2_benchmark_summary_long.csv: runs/<run_id>/s6_task2_result_synthesis/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s6_task2_result_synthesis/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s6_task2_result_synthesis.py --run-id <run_id> --s4-stage-dir <path> --s5-stage-dir <path> --seed 619
# Failure Modes:
#   - Missing required upstream files -> exit non-zero
#   - Schema / key / representation-order / target-order / assertion drift -> exit non-zero
#   - Attempted edist ranking or incorrect benchmark summary row count -> exit non-zero
# Last Updated: 2026-03-13

"""
Inputs:
- audited S4/S5 stage outputs only (no task2_retrieval_per_query.parquet required)

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
- Representation order = Gene, Pathway, scgpt, geneformer, scbert, scfoundation, uce, state, tahoe-x1
- Canonical target order = file order from Common_Targets_K562.csv
- target_membership_source = delta_meta.target_tokens
- Benchmark summary row count = (9 x 3) + (9 x 2 x 4) = 99

Interpretation rules:
- mean_edist_biascorr is informational only
- mean_edist_biascorr is not leaderboard-eligible, not rankable, not cross-representation comparable
- C2G and G2C remain separate; no collapsed retrieval score is emitted
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

import pandas as pd
import yaml

STAGE = "s6_task2_result_synthesis"
CONFIG_PATH = Path("config/config.yaml")
GLOBAL_SEED = 619
CHANCE_IDENTITY_TOL = 1e-12
MAX_COUNTEREXAMPLES = 5

EXPECTED_S4_STAGE = "s4_task2_group_concordance"
EXPECTED_S5_STAGE = "s5_task2_retrieval"
EXPECTED_DATASET = "scPerturb"
EXPECTED_CELL_LINE = "K562"
TARGET_MEMBERSHIP_SOURCE = "delta_meta.target_tokens"

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
TARGET_ORDER: Tuple[str, ...] = (
    "BCL2L1",
    "BCR",
    "BRD4",
    "CA5A",
    "CDK1",
    "CDK2",
    "CDK9",
    "CHEK1",
    "DNMT1",
    "EEF2",
    "EGLN2",
    "HDAC7",
    "KDM2A",
    "MTOR",
    "PTK2",
    "STAT5A",
    "TUBB",
)
DIRECTION_ORDER: Tuple[str, ...] = ("C2G", "G2C")
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
EXPECTED_BENCHMARK_ROWS = (len(REPRESENTATION_ORDER) * len(GROUP_SUMMARY_METRIC_ORDER)) + (
    len(REPRESENTATION_ORDER) * len(DIRECTION_ORDER) * len(RETRIEVAL_METRICS)
)

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
    "direction",
    "dataset",
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
    parser = argparse.ArgumentParser(description="S6 Task2 result synthesis")
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


def canonical_rep_order(values: Iterable[str]) -> List[str]:
    seen: List[str] = []
    for value in values:
        if value not in seen:
            seen.append(value)
    return seen


def require_stage_files(stage_dir: Path, required_files: Sequence[str]) -> Dict[str, Path]:
    if not stage_dir.exists():
        raise FileNotFoundError(f"Stage directory does not exist: {stage_dir}")
    if not stage_dir.is_dir():
        raise NotADirectoryError(f"Stage path is not a directory: {stage_dir}")

    paths: Dict[str, Path] = {}
    for filename in required_files:
        path = (stage_dir / filename).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Missing required upstream file: {path}")
        paths[filename] = path
    return paths


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


def extract_non_blocking_notes(run_manifest: Mapping[str, Any]) -> List[str]:
    notes: List[str] = []
    config = run_manifest.get("config")
    if isinstance(config, Mapping):
        for key, value in config.items():
            if "note" in str(key).lower() and isinstance(value, str) and value.strip():
                notes.append(f"{key}: {value.strip()}")
    return notes


def assert_all_upstream_assertions_pass(bundle: StageBundle, stage_name: str) -> None:
    raw_assertions = bundle.audit_assertions.get("assertions")
    if not isinstance(raw_assertions, list):
        raise ValueError(f"{stage_name} audit_assertions.json missing 'assertions' list")

    failing = [item for item in raw_assertions if not bool(item.get("pass", False))]
    if failing:
        names = [str(item.get("name", "<unnamed>")) for item in failing[:MAX_COUNTEREXAMPLES]]
        raise ValueError(f"{stage_name} has failing required upstream assertions: {names}")


def validate_manifest_file_list(bundle: StageBundle, required_files: Sequence[str], stage_name: str) -> None:
    manifest_stage = bundle.manifest.get("stage")
    if manifest_stage != stage_name:
        raise ValueError(f"{stage_name} manifest stage mismatch: {manifest_stage}")

    files = bundle.manifest.get("files")
    if not isinstance(files, list):
        raise ValueError(f"{stage_name} manifest missing file list")
    available = {str(item.get("relative_path", "")) for item in files if isinstance(item, Mapping)}
    missing = [filename for filename in required_files if filename != "manifest.json" and filename not in available]
    if missing:
        raise ValueError(f"{stage_name} manifest.json missing required files: {missing}")


def build_representation_details_map(run_manifest: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    summary = run_manifest.get("summary")
    if not isinstance(summary, Mapping):
        raise ValueError("run_manifest missing summary block")
    details = summary.get("representation_details")
    if not isinstance(details, list):
        raise ValueError("run_manifest missing representation_details")
    out: Dict[str, Dict[str, Any]] = {}
    for detail in details:
        if isinstance(detail, Mapping) and "representation" in detail:
            out[str(detail["representation"])] = dict(detail)
    return out


def validate_s4_bundle(bundle: StageBundle) -> Dict[str, Dict[str, Any]]:
    if bundle.run_manifest.get("stage") != EXPECTED_S4_STAGE:
        raise ValueError(f"S4 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S4_REQUIRED_FILES, EXPECTED_S4_STAGE)
    assert_all_upstream_assertions_pass(bundle, EXPECTED_S4_STAGE)

    config = bundle.run_manifest.get("config")
    if not isinstance(config, Mapping):
        raise ValueError("S4 run_manifest missing config block")
    if list(config.get("representation_order", [])) != list(REPRESENTATION_ORDER):
        raise ValueError("S4 representation_order drift detected")
    if list(config.get("canonical_target_order", [])) != list(TARGET_ORDER):
        raise ValueError("S4 canonical_target_order drift detected")
    if config.get("target_membership_source") != TARGET_MEMBERSHIP_SOURCE:
        raise ValueError("S4 target_membership_source drift detected")

    concordance = bundle.tables["task2_group_concordance.csv"]
    attrition = bundle.tables["task2_group_attrition.csv"]
    ensure_required_columns(concordance, S4_CONCORDANCE_COLUMNS, "S4 task2_group_concordance.csv")
    ensure_required_columns(attrition, S4_ATTRITION_COLUMNS, "S4 task2_group_attrition.csv")

    if len(concordance) != len(REPRESENTATION_ORDER) * len(TARGET_ORDER):
        raise ValueError(f"S4 concordance row count mismatch: {len(concordance)}")
    if set(concordance["dataset"].dropna().unique()) != {EXPECTED_DATASET}:
        raise ValueError("S4 dataset drift detected")
    if set(concordance["cell_line"].dropna().unique()) != {EXPECTED_CELL_LINE}:
        raise ValueError("S4 cell_line drift detected")

    rep_order = canonical_rep_order(concordance["representation"].tolist())
    if rep_order != list(REPRESENTATION_ORDER):
        raise ValueError(f"S4 output representation order drift: {rep_order}")
    for representation in REPRESENTATION_ORDER:
        target_seq = concordance.loc[
            concordance["representation"] == representation, "target_token"
        ].tolist()
        if target_seq != list(TARGET_ORDER):
            raise ValueError(f"S4 target order drift for representation={representation}: {target_seq}")

    if not attrition.empty:
        if set(attrition["representation"].dropna().unique()) != {"uce"}:
            raise ValueError("S4 attrition is no longer UCE-only")
        if set(attrition["dataset"].dropna().unique()) != {EXPECTED_DATASET}:
            raise ValueError("S4 attrition dataset drift detected")
        if set(attrition["cell_line"].dropna().unique()) != {EXPECTED_CELL_LINE}:
            raise ValueError("S4 attrition cell_line drift detected")

    return build_representation_details_map(bundle.run_manifest)


def validate_s5_bundle(bundle: StageBundle) -> Dict[str, Dict[str, Any]]:
    if bundle.run_manifest.get("stage") != EXPECTED_S5_STAGE:
        raise ValueError(f"S5 stage mismatch: {bundle.run_manifest.get('stage')}")
    validate_manifest_file_list(bundle, S5_REQUIRED_FILES, EXPECTED_S5_STAGE)
    assert_all_upstream_assertions_pass(bundle, EXPECTED_S5_STAGE)

    config = bundle.run_manifest.get("config")
    if not isinstance(config, Mapping):
        raise ValueError("S5 run_manifest missing config block")
    if list(config.get("representation_order", [])) != list(REPRESENTATION_ORDER):
        raise ValueError("S5 representation_order drift detected")
    if list(config.get("direction_order", [])) != list(DIRECTION_ORDER):
        raise ValueError("S5 direction_order drift detected")
    if list(config.get("canonical_target_order", [])) != list(TARGET_ORDER):
        raise ValueError("S5 canonical_target_order drift detected")
    if config.get("target_membership_source") != TARGET_MEMBERSHIP_SOURCE:
        raise ValueError("S5 target_membership_source drift detected")

    summary = bundle.tables["task2_retrieval_summary.csv"]
    summary_long = bundle.tables["task2_retrieval_summary_long.csv"]
    attrition = bundle.tables["task2_retrieval_attrition.csv"]
    chance = bundle.tables["task2_chance_identity_check.csv"]
    ensure_required_columns(summary, S5_SUMMARY_COLUMNS, "S5 task2_retrieval_summary.csv")
    ensure_required_columns(summary_long, S5_SUMMARY_LONG_COLUMNS, "S5 task2_retrieval_summary_long.csv")
    ensure_required_columns(attrition, S5_ATTRITION_COLUMNS, "S5 task2_retrieval_attrition.csv")
    ensure_required_columns(chance, S5_CHANCE_COLUMNS, "S5 task2_chance_identity_check.csv")

    if len(summary) != len(REPRESENTATION_ORDER) * len(DIRECTION_ORDER):
        raise ValueError(f"S5 summary row count mismatch: {len(summary)}")
    if len(summary_long) != len(REPRESENTATION_ORDER) * len(DIRECTION_ORDER) * len(RETRIEVAL_METRICS):
        raise ValueError(f"S5 summary_long row count mismatch: {len(summary_long)}")
    if len(chance) != len(REPRESENTATION_ORDER) * len(DIRECTION_ORDER):
        raise ValueError(f"S5 chance row count mismatch: {len(chance)}")
    if set(summary["dataset"].dropna().unique()) != {EXPECTED_DATASET}:
        raise ValueError("S5 dataset drift detected")
    if set(summary["cell_line"].dropna().unique()) != {EXPECTED_CELL_LINE}:
        raise ValueError("S5 cell_line drift detected")

    actual_directions = canonical_rep_order(summary["direction"].tolist())
    if actual_directions != list(DIRECTION_ORDER):
        raise ValueError(f"S5 direction order drift: {actual_directions}")
    for direction in DIRECTION_ORDER:
        rep_seq = summary.loc[summary["direction"] == direction, "representation"].tolist()
        if rep_seq != list(REPRESENTATION_ORDER):
            raise ValueError(f"S5 representation order drift for direction={direction}: {rep_seq}")

    long_metric_names = summary_long["metric_name"].drop_duplicates().tolist()
    if long_metric_names != list(RETRIEVAL_METRICS):
        raise ValueError(f"S5 summary_long metric drift: {long_metric_names}")
    if set(attrition["representation"].dropna().unique()) - {"uce"}:
        raise ValueError("S5 attrition is no longer UCE-only")

    for column in ["abs_delta_mrr", "abs_delta_hit1", "abs_delta_hit5", "abs_delta_hit10"]:
        if (chance[column].astype(float).abs() > CHANCE_IDENTITY_TOL).any():
            raise ValueError(f"S5 chance identity tolerance exceeded for {column}")

    return build_representation_details_map(bundle.run_manifest)


def validate_cross_stage_contracts(
    s4_rep_details: Mapping[str, Mapping[str, Any]],
    s5_rep_details: Mapping[str, Mapping[str, Any]],
) -> None:
    for representation in REPRESENTATION_ORDER:
        if representation not in s4_rep_details:
            raise ValueError(f"S4 representation_details missing {representation}")
        if representation not in s5_rep_details:
            raise ValueError(f"S5 representation_details missing {representation}")

    for stage_name, details in ((EXPECTED_S4_STAGE, s4_rep_details), (EXPECTED_S5_STAGE, s5_rep_details)):
        uce = details["uce"]
        n_invalid = int(uce.get("n_invalid", -1))
        invalid_by_side = dict(uce.get("invalid_by_side", {}))
        if n_invalid != 384:
            raise ValueError(f"{stage_name} UCE invalid row drift: {n_invalid}")
        if int(invalid_by_side.get("DRUG", -1)) != 384:
            raise ValueError(f"{stage_name} UCE invalid DRUG row drift: {invalid_by_side}")


def rank_desc(series: pd.Series) -> pd.Series:
    out = pd.Series(pd.NA, index=series.index, dtype="Int64")
    valid = series.notna()
    if valid.any():
        out.loc[valid] = series.loc[valid].rank(method="min", ascending=False).astype("Int64")
    return out


def build_group_concordance_long(concordance: pd.DataFrame) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    for metric_name, _, valid_col, na_reason_col in GROUP_METRIC_SPECS:
        frame = concordance[
            [
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
                metric_name,
                valid_col,
                na_reason_col,
            ]
        ].copy()
        frame["metric_name"] = metric_name
        frame["metric_value"] = frame[metric_name]
        frame["metric_valid_bool"] = frame[valid_col].astype(bool)
        frame["metric_na_reason"] = frame[na_reason_col].fillna("").astype(str)
        frame = frame.drop(columns=[metric_name, valid_col, na_reason_col])
        rows.append(frame)

    out = pd.concat(rows, ignore_index=True)
    rep_index = {name: idx for idx, name in enumerate(REPRESENTATION_ORDER)}
    target_index = {name: idx for idx, name in enumerate(TARGET_ORDER)}
    metric_index = {name: idx for idx, name in enumerate(GROUP_METRIC_ORDER)}
    out["_rep_idx"] = out["representation"].map(rep_index)
    out["_target_idx"] = out["target_token"].map(target_index)
    out["_metric_idx"] = out["metric_name"].map(metric_index)
    out = out.sort_values(["_rep_idx", "_target_idx", "_metric_idx"], kind="stable")
    return out.drop(columns=["_rep_idx", "_target_idx", "_metric_idx"]).reset_index(drop=True)[GROUP_LONG_COLUMNS]


def build_group_leaderboard(
    concordance: pd.DataFrame,
    attrition: pd.DataFrame,
    rep_details: Mapping[str, Mapping[str, Any]],
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for representation in REPRESENTATION_ORDER:
        subset = concordance.loc[concordance["representation"] == representation].reset_index(drop=True)
        if subset.empty:
            raise ValueError(f"Missing S4 rows for representation={representation}")
        attr = attrition.loc[attrition["representation"] == representation].reset_index(drop=True)
        detail = rep_details[representation]
        invalid_by_side = dict(detail.get("invalid_by_side", {}))
        rows.append(
            {
                "dataset": subset.at[0, "dataset"],
                "cell_line": subset.at[0, "cell_line"],
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
                "n_attrition_chem_rows_removed_membership": int(attr["n_chem_removed"].sum()) if not attr.empty else 0,
                "n_attrition_gen_rows_removed_membership": int(attr["n_gen_removed"].sum()) if not attr.empty else 0,
            }
        )

    out = pd.DataFrame(rows)
    out["rank_by_mean_cosine_centroid"] = rank_desc(out["mean_cosine_centroid"])
    out["rank_by_mean_pcc_centroid"] = rank_desc(out["mean_pcc_centroid"])
    rep_index = {name: idx for idx, name in enumerate(REPRESENTATION_ORDER)}
    out["_rep_idx"] = out["representation"].map(rep_index)
    out = out.sort_values(["rank_by_mean_cosine_centroid", "_rep_idx"], kind="stable", na_position="last")
    return out.drop(columns=["_rep_idx"]).reset_index(drop=True)[GROUP_LEADERBOARD_COLUMNS]


def build_retrieval_leaderboard(summary: pd.DataFrame) -> pd.DataFrame:
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
        out.groupby("direction", sort=False)["mean_mrr_corrected"].transform(rank_desc).astype("Int64")
    )
    rep_index = {name: idx for idx, name in enumerate(REPRESENTATION_ORDER)}
    dir_index = {name: idx for idx, name in enumerate(DIRECTION_ORDER)}
    out["_dir_idx"] = out["direction"].map(dir_index)
    out["_rep_idx"] = out["representation"].map(rep_index)
    out = out.sort_values(
        ["_dir_idx", "rank_by_mean_mrr_corrected", "_rep_idx"],
        kind="stable",
        na_position="last",
    )
    return out.drop(columns=["_dir_idx", "_rep_idx"]).reset_index(drop=True)[RETRIEVAL_LEADERBOARD_COLUMNS]


def join_caution_codes(*codes: Optional[str]) -> str:
    ordered = [code for code in codes if code]
    return ";".join(ordered)


def build_benchmark_summary_long(
    group_leaderboard: pd.DataFrame,
    retrieval_leaderboard: pd.DataFrame,
    s4_run_id: str,
    s5_run_id: str,
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    for row in group_leaderboard.to_dict(orient="records"):
        base_codes = ["UCE_ATTRITION_PRESENT"] if int(row["n_invalid_rows_unique"]) > 0 else []

        rows.append(
            {
                "analysis_family": "group_concordance",
                "dataset": row["dataset"],
                "cell_line": row["cell_line"],
                "direction": "NA",
                "representation": row["representation"],
                "metric_name": "mean_cosine_centroid",
                "metric_value": row["mean_cosine_centroid"],
                "rank_value": row["rank_by_mean_cosine_centroid"],
                "rank_basis_metric_name": "mean_cosine_centroid",
                "leaderboard_eligible_bool": True,
                "cross_representation_comparable_bool": True,
                "source_table": "task2_group_leaderboard.csv",
                "source_stage": EXPECTED_S4_STAGE,
                "source_run_id": s4_run_id,
                "caution_codes": join_caution_codes(*base_codes),
                "n_targets_total": row["n_targets_total"],
                "n_targets_metric_valid": row["n_targets_cosine_valid"],
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
        rows.append(
            {
                "analysis_family": "group_concordance",
                "dataset": row["dataset"],
                "cell_line": row["cell_line"],
                "direction": "NA",
                "representation": row["representation"],
                "metric_name": "mean_pcc_centroid",
                "metric_value": row["mean_pcc_centroid"],
                "rank_value": row["rank_by_mean_pcc_centroid"],
                "rank_basis_metric_name": "mean_pcc_centroid",
                "leaderboard_eligible_bool": True,
                "cross_representation_comparable_bool": True,
                "source_table": "task2_group_leaderboard.csv",
                "source_stage": EXPECTED_S4_STAGE,
                "source_run_id": s4_run_id,
                "caution_codes": join_caution_codes(*base_codes),
                "n_targets_total": row["n_targets_total"],
                "n_targets_metric_valid": row["n_targets_pcc_valid"],
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
        rows.append(
            {
                "analysis_family": "group_concordance",
                "dataset": row["dataset"],
                "cell_line": row["cell_line"],
                "direction": "NA",
                "representation": row["representation"],
                "metric_name": "mean_edist_biascorr",
                "metric_value": row["mean_edist_biascorr"],
                "rank_value": pd.NA,
                "rank_basis_metric_name": "",
                "leaderboard_eligible_bool": False,
                "cross_representation_comparable_bool": False,
                "source_table": "task2_group_leaderboard.csv",
                "source_stage": EXPECTED_S4_STAGE,
                "source_run_id": s4_run_id,
                "caution_codes": join_caution_codes("GROUP_EDIST_NOT_CROSS_REP", *base_codes),
                "n_targets_total": row["n_targets_total"],
                "n_targets_metric_valid": row["n_targets_edist_valid"],
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

    for row in retrieval_leaderboard.to_dict(orient="records"):
        base_codes = ["RETRIEVAL_DIRECTION_SPECIFIC"]
        if row["representation"] == "uce":
            base_codes.append("UCE_ATTRITION_PRESENT")
        for metric_name in RETRIEVAL_METRICS:
            rows.append(
                {
                    "analysis_family": "retrieval",
                    "dataset": row["dataset"],
                    "cell_line": row["cell_line"],
                    "direction": row["direction"],
                    "representation": row["representation"],
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
    family_index = {"group_concordance": 0, "retrieval": 1}
    direction_index = {"NA": 0, **{name: idx + 1 for idx, name in enumerate(DIRECTION_ORDER)}}
    metric_index = {name: idx for idx, name in enumerate(GROUP_SUMMARY_METRIC_ORDER)}
    metric_index.update({name: idx for idx, name in enumerate(RETRIEVAL_METRICS, start=len(metric_index))})
    rep_index = {name: idx for idx, name in enumerate(REPRESENTATION_ORDER)}
    out["_family_idx"] = out["analysis_family"].map(family_index)
    out["_direction_idx"] = out["direction"].map(direction_index)
    out["_metric_idx"] = out["metric_name"].map(metric_index)
    out["_rep_idx"] = out["representation"].map(rep_index)
    out = out.sort_values(
        ["_family_idx", "_direction_idx", "_metric_idx", "rank_value", "_rep_idx"],
        kind="stable",
        na_position="last",
    )
    return out.drop(columns=["_family_idx", "_direction_idx", "_metric_idx", "_rep_idx"]).reset_index(drop=True)[
        BENCHMARK_SUMMARY_LONG_COLUMNS
    ]


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

    s4_required_present = all((s4_bundle.stage_dir / name).exists() for name in S4_REQUIRED_FILES)
    s5_required_present = all((s5_bundle.stage_dir / name).exists() for name in S5_REQUIRED_FILES)
    assertions.append(
        {
            "name": "required_upstream_files_present",
            "pass": bool(s4_required_present and s5_required_present),
            "details": {
                "s4_stage_dir": str(s4_bundle.stage_dir),
                "s5_stage_dir": str(s5_bundle.stage_dir),
                "s4_required_files": list(S4_REQUIRED_FILES),
                "s5_required_files": list(S5_REQUIRED_FILES),
                "task2_retrieval_per_query_required": False,
            },
            "counterexamples": [] if (s4_required_present and s5_required_present) else [{"s4": s4_required_present, "s5": s5_required_present}],
        }
    )

    assertions.append(
        {
            "name": "representation_registry_frozen",
            "pass": True,
            "details": {
                "representation_order": list(REPRESENTATION_ORDER),
                "s4_representation_order": list(s4_bundle.run_manifest["config"]["representation_order"]),
                "s5_representation_order": list(s5_bundle.run_manifest["config"]["representation_order"]),
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "canonical_target_order_frozen",
            "pass": True,
            "details": {
                "canonical_target_order": list(TARGET_ORDER),
                "s4_target_order": list(s4_bundle.run_manifest["config"]["canonical_target_order"]),
                "s5_target_order": list(s5_bundle.run_manifest["config"]["canonical_target_order"]),
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "target_membership_source_locked",
            "pass": True,
            "details": {
                "target_membership_source": TARGET_MEMBERSHIP_SOURCE,
                "s4_target_membership_source": s4_bundle.run_manifest["config"]["target_membership_source"],
                "s5_target_membership_source": s5_bundle.run_manifest["config"]["target_membership_source"],
            },
            "counterexamples": [],
        }
    )

    group_long_pass = len(group_long) == len(REPRESENTATION_ORDER) * len(TARGET_ORDER) * len(GROUP_METRIC_ORDER)
    assertions.append(
        {
            "name": "group_concordance_long_contract",
            "pass": bool(group_long_pass),
            "details": {
                "source_table": "task2_group_concordance.csv",
                "output_table": "task2_group_concordance_long.csv",
                "expected_rows": len(REPRESENTATION_ORDER) * len(TARGET_ORDER) * len(GROUP_METRIC_ORDER),
                "actual_rows": int(len(group_long)),
                "sort_order": "representation_order -> canonical_target_order -> metric_order",
            },
            "counterexamples": [] if group_long_pass else [{"actual_rows": int(len(group_long))}],
        }
    )

    edist_cols_present = "rank_by_mean_edist_biascorr" not in group_leaderboard.columns
    edist_flags = (
        (group_leaderboard["mean_edist_biascorr_leaderboard_eligible_bool"] == False).all()  # noqa: E712
        and (group_leaderboard["mean_edist_biascorr_cross_representation_comparable_bool"] == False).all()  # noqa: E712
    )
    group_leaderboard_pass = len(group_leaderboard) == len(REPRESENTATION_ORDER) and edist_cols_present and edist_flags
    assertions.append(
        {
            "name": "group_leaderboard_edist_contract",
            "pass": bool(group_leaderboard_pass),
            "details": {
                "expected_rows": len(REPRESENTATION_ORDER),
                "actual_rows": int(len(group_leaderboard)),
                "edist_leaderboard_eligible": False,
                "edist_cross_representation_comparable": False,
                "edist_rank_column_present": not edist_cols_present,
            },
            "counterexamples": [] if group_leaderboard_pass else group_leaderboard.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    retrieval_leaderboard_pass = len(retrieval_leaderboard) == len(REPRESENTATION_ORDER) * len(DIRECTION_ORDER)
    assertions.append(
        {
            "name": "retrieval_leaderboard_contract",
            "pass": bool(retrieval_leaderboard_pass),
            "details": {
                "expected_rows": len(REPRESENTATION_ORDER) * len(DIRECTION_ORDER),
                "actual_rows": int(len(retrieval_leaderboard)),
                "directions": list(DIRECTION_ORDER),
                "rank_basis": "mean_mrr_corrected",
            },
            "counterexamples": [] if retrieval_leaderboard_pass else retrieval_leaderboard.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    benchmark_pass = len(benchmark_summary_long) == EXPECTED_BENCHMARK_ROWS
    group_rows = int((benchmark_summary_long["analysis_family"] == "group_concordance").sum())
    retrieval_rows = int((benchmark_summary_long["analysis_family"] == "retrieval").sum())
    edist_rows = benchmark_summary_long.loc[
        benchmark_summary_long["metric_name"] == "mean_edist_biascorr"
    ].reset_index(drop=True)
    edist_unranked = edist_rows["rank_value"].isna().all()
    edist_ineligible = (edist_rows["leaderboard_eligible_bool"] == False).all()  # noqa: E712
    edist_non_comparable = (edist_rows["cross_representation_comparable_bool"] == False).all()  # noqa: E712
    edist_source_ok = set(edist_rows["source_table"].tolist()) == {"task2_group_leaderboard.csv"}
    benchmark_pass = benchmark_pass and group_rows == 27 and retrieval_rows == 72 and edist_unranked and edist_ineligible and edist_non_comparable and edist_source_ok
    assertions.append(
        {
            "name": "benchmark_summary_long_contract",
            "pass": bool(benchmark_pass),
            "details": {
                "expected_rows": EXPECTED_BENCHMARK_ROWS,
                "actual_rows": int(len(benchmark_summary_long)),
                "expected_group_rows": 27,
                "actual_group_rows": group_rows,
                "expected_retrieval_rows": 72,
                "actual_retrieval_rows": retrieval_rows,
                "group_source_table": "task2_group_leaderboard.csv",
                "retrieval_source_table": "task2_retrieval_leaderboard.csv",
                "edist_rows_unranked": bool(edist_unranked),
            },
            "counterexamples": [] if benchmark_pass else benchmark_summary_long.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    uce_rows = benchmark_summary_long.loc[benchmark_summary_long["representation"] == "uce", "caution_codes"].astype(str)
    uce_note_pass = uce_rows.str.contains("UCE_ATTRITION_PRESENT", regex=False).all()
    assertions.append(
        {
            "name": "uce_attrition_note_preserved",
            "pass": bool(uce_note_pass),
            "details": {
                "rules": [
                    "UCE rows must preserve UCE_ATTRITION_PRESENT caution code",
                    "S4 and S5 audited UCE invalid DRUG row count remains 384",
                ],
                "uce_row_count": int((benchmark_summary_long["representation"] == "uce").sum()),
            },
            "counterexamples": [] if uce_note_pass else uce_rows.head(MAX_COUNTEREXAMPLES).tolist(),
        }
    )

    outputs_root_pass = all(path.resolve().is_relative_to(stage_dir.resolve()) for path in [
        stage_dir / "task2_group_concordance_long.csv",
        stage_dir / "task2_group_leaderboard.csv",
        stage_dir / "task2_retrieval_leaderboard.csv",
        stage_dir / "task2_benchmark_summary_long.csv",
    ])
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
    s4_bundle = load_stage_bundle(
        (project_root / args.s4_stage_dir).resolve() if not args.s4_stage_dir.is_absolute() else args.s4_stage_dir.resolve(),
        S4_REQUIRED_FILES,
        ("task2_group_concordance.csv", "task2_group_attrition.csv"),
    )
    s5_bundle = load_stage_bundle(
        (project_root / args.s5_stage_dir).resolve() if not args.s5_stage_dir.is_absolute() else args.s5_stage_dir.resolve(),
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
    validate_cross_stage_contracts(s4_rep_details, s5_rep_details)

    s4_concordance = s4_bundle.tables["task2_group_concordance.csv"]
    s4_attrition = s4_bundle.tables["task2_group_attrition.csv"]
    s5_summary = s5_bundle.tables["task2_retrieval_summary.csv"]

    group_long = build_group_concordance_long(s4_concordance)
    group_leaderboard = build_group_leaderboard(s4_concordance, s4_attrition, s4_rep_details)
    retrieval_leaderboard = build_retrieval_leaderboard(s5_summary)
    benchmark_summary_long = build_benchmark_summary_long(
        group_leaderboard=group_leaderboard,
        retrieval_leaderboard=retrieval_leaderboard,
        s4_run_id=str(s4_bundle.run_manifest["run_id"]),
        s5_run_id=str(s5_bundle.run_manifest["run_id"]),
    )

    if len(group_long) != 459:
        raise ValueError(f"task2_group_concordance_long.csv row count mismatch: {len(group_long)}")
    if len(group_leaderboard) != 9:
        raise ValueError(f"task2_group_leaderboard.csv row count mismatch: {len(group_leaderboard)}")
    if len(retrieval_leaderboard) != 18:
        raise ValueError(f"task2_retrieval_leaderboard.csv row count mismatch: {len(retrieval_leaderboard)}")
    if len(benchmark_summary_long) != EXPECTED_BENCHMARK_ROWS:
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
        non_blocking_notes.append(
            "upstream_git_head_mismatch: S4 and S5 audited runs were generated from different git heads"
        )

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
        raise ValueError(f"Internal Script6 assertion failure: {failed_assertions[:MAX_COUNTEREXAMPLES]}")

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts" / "s6_task2_result_synthesis.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": utc_now_iso(),
        "config": {
            "seed": seed,
            "runs_dir": str(runs_dir),
            "representation_order": list(REPRESENTATION_ORDER),
            "canonical_target_order": list(TARGET_ORDER),
            "target_membership_source": TARGET_MEMBERSHIP_SOURCE,
            "s4_stage_dir": str(s4_bundle.stage_dir),
            "s5_stage_dir": str(s5_bundle.stage_dir),
            "s4_source_run_id": s4_bundle.run_manifest["run_id"],
            "s5_source_run_id": s5_bundle.run_manifest["run_id"],
            "benchmark_summary_expected_rows": EXPECTED_BENCHMARK_ROWS,
            "task2_retrieval_per_query_required": False,
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
            "n_group_concordance_long_rows": int(len(group_long)),
            "n_group_leaderboard_rows": int(len(group_leaderboard)),
            "n_retrieval_leaderboard_rows": int(len(retrieval_leaderboard)),
            "n_benchmark_summary_long_rows": int(len(benchmark_summary_long)),
            "expected_benchmark_summary_long_rows": EXPECTED_BENCHMARK_ROWS,
            "group_rows_in_benchmark_summary": int((benchmark_summary_long["analysis_family"] == "group_concordance").sum()),
            "retrieval_rows_in_benchmark_summary": int((benchmark_summary_long["analysis_family"] == "retrieval").sum()),
            "non_blocking_provenance_notes": non_blocking_notes,
        },
    }
    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})

    manifest_payload = {
        "stage": STAGE,
        "files": [],
    }
    for path in [
        audit_assertions_path,
        benchmark_summary_long_path,
        group_long_path,
        group_leaderboard_path,
        retrieval_leaderboard_path,
        run_manifest_path,
    ]:
        manifest_payload["files"].append(
            {
                "relative_path": path.relative_to(stage_dir).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": compute_sha256(path),
            }
        )
    write_json(manifest_path, manifest_payload)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)
