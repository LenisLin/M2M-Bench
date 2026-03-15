#!/usr/bin/env python3
# SCRIPT_HEADER_CONTRACT
# Script: scripts/s4_task2_group_concordance_multisource.py
# Purpose: Compute corrected Task2 multisource group concordance per eligible
#   dataset/cell_line/target_token/representation from the audited
#   data/task2_snapshot_v2 snapshot only.
# Inputs:
#   - data/task2_snapshot_v2/snapshot_manifest.json
#   - data/task2_snapshot_v2/task2_pairs_coverage.csv
#   - data/task2_snapshot_v2/representation_availability_registry.csv
#   - data/task2_snapshot_v2/lincs/task2_lincs_pairs.csv
#   - data/task2_snapshot_v2/lincs/derived/{delta_meta.csv,task2_row_membership.parquet,gene_delta.npy,pathway_delta.npy}
#   - data/task2_snapshot_v2/scperturb_k562/Common_Targets_K562.csv
#   - data/task2_snapshot_v2/scperturb_k562/derived/{delta_meta.csv,task2_row_membership.parquet,gene_delta.npy,pathway_delta.npy}
#   - data/task2_snapshot_v2/scperturb_k562/fm/<model>/{fm_delta.npy,fm_delta_meta.csv,delta_operator_policy.json}
# Outputs:
#   - runs/<run_id>/s4_task2_group_concordance_multisource/task2_group_concordance.csv
#   - runs/<run_id>/s4_task2_group_concordance_multisource/task2_group_attrition.csv
#   - runs/<run_id>/s4_task2_group_concordance_multisource/run_manifest.json
#   - runs/<run_id>/s4_task2_group_concordance_multisource/audit_assertions.json
#   - runs/<run_id>/s4_task2_group_concordance_multisource/manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s4_task2_group_concordance_multisource/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::compute.workers
# Execution:
#   - python scripts/s4_task2_group_concordance_multisource.py --run-id <run_id> --seed 619 --workers 8
# Failure Modes:
#   - Snapshot root is not data/task2_snapshot_v2 -> exit non-zero
#   - Missing required snapshot files -> exit non-zero
#   - Coverage / registry / schema / key / scope drift -> exit non-zero
#   - Representation validation failure or forbidden 1HAE analysis label -> exit non-zero
#   - Supported analysis grid row count mismatch -> exit non-zero
# Last Updated: 2026-03-10

"""
Inputs:
- corrected Task2 multisource snapshot under data/task2_snapshot_v2/

Outputs:
- task2_group_concordance.csv
- task2_group_attrition.csv
- run_manifest.json
- audit_assertions.json
- manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- EDIST_MAX_N = 256
- Dataset order = LINCS, scPerturb
- Representation order = Gene, Pathway, scgpt, geneformer, scbert, scfoundation, uce, state, tahoe-x1
- Canonical membership source = delta_meta.target_tokens

Attrition rules:
- Representation valid_mask drops are recorded before metric computation.
- If either side has zero valid rows after filtering, all metrics are NA and reason=missing_one_side.
- E-distance is NA unless n_chem_sub >= 2 and n_gen_sub >= 2 after deterministic subsampling.
- not_applicable_scope is a registry-only policy and produces no core rows and no attrition rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml

STAGE = "s4_task2_group_concordance_multisource"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v2")

GLOBAL_SEED = 619
EDIST_MAX_N = 256
MASK_BATCH_ROWS = 256
MAX_COUNTEREXAMPLES = 5

DATASET_ORDER: Tuple[str, str] = ("LINCS", "scPerturb")
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
SIDE_ORDER: Tuple[str, str] = ("Chemical", "Genetic")
METRIC_ORDER: Tuple[str, ...] = ("ALL", "cosine_centroid", "pcc_centroid", "edist_biascorr")

STATUS_AVAILABLE = "available"
STATUS_NOT_APPLICABLE = "not_applicable_scope"
FORBIDDEN_ANALYSIS_CELL_LINE = "1HAE"

FROZEN_AVAILABLE_REPRESENTATIONS: Mapping[str, Tuple[str, ...]] = {
    "LINCS": ("Gene", "Pathway"),
    "scPerturb": REPRESENTATION_ORDER,
}
EXPECTED_AVAILABILITY_REASON: Mapping[str, str] = {
    STATUS_AVAILABLE: "representation_supported_for_dataset_scope",
    STATUS_NOT_APPLICABLE: "representation_not_supported_for_dataset_scope",
}

CONCORDANCE_COLUMNS = [
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

ATTRITION_COLUMNS = [
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

_WORKER_ARRAY: Optional[np.ndarray] = None
_WORKER_REPRESENTATION: Optional[str] = None
_PROCESS_POOL_PERMISSION_DENIED = False


@dataclass(frozen=True)
class DatasetBundle:
    dataset: str
    subtree_relpath: str
    subtree_manifest_relpath: str
    delta_meta_relpath: str
    membership_relpath: str
    mech_order_relpath: str
    common_targets_relpath: Optional[str] = None


@dataclass(frozen=True)
class RepresentationSpec:
    dataset: str
    name: str
    family: str
    array_relpath: str
    meta_relpath: str
    expected_dim: int
    policy_relpath: Optional[str] = None


@dataclass(frozen=True)
class RepresentationPayload:
    spec: RepresentationSpec
    array_path: Path
    meta_path: Path
    valid_mask: np.ndarray
    policy_path: Optional[Path]
    n_valid: int
    n_invalid: int
    invalid_by_side: Mapping[str, int]


@dataclass(frozen=True)
class GroupTask:
    dataset: str
    cell_line: str
    target_token: str
    representation: str
    group_id: str
    n_chem_total: int
    n_gen_total: int
    chem_used_row_ids: np.ndarray
    gen_used_row_ids: np.ndarray
    chem_sub_row_ids: np.ndarray
    gen_sub_row_ids: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S4 corrected Task2 multisource group concordance")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Override config.compute.workers",
    )
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


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=True)
        handle.write("\n")


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def ensure_required_columns(frame: pd.DataFrame, required: Sequence[str], name: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns for {name}: {missing}")


def init_global_rng(seed: int) -> np.random.Generator:
    random.seed(seed)
    np.random.seed(seed)
    return np.random.default_rng(seed)


def parse_bool_series(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(dtype=bool, copy=False)

    normalized = series.astype("string").str.strip().str.lower()
    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
        "y": True,
        "n": False,
        "t": True,
        "f": False,
    }
    unknown = normalized[~normalized.isin(mapping.keys()) & normalized.notna()]
    if not unknown.empty:
        examples = unknown.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"Unsupported boolean values: {examples}")
    return normalized.map(mapping).fillna(False).to_numpy(dtype=bool, copy=False)


def parse_target_tokens_field(value: object) -> Tuple[str, ...]:
    if pd.isna(value):
        raise ValueError("target_tokens contains NA")
    raw = str(value).strip()
    if not raw:
        raise ValueError("target_tokens contains empty string")
    tokens = tuple(token.strip() for token in raw.split(";") if token.strip())
    if not tokens:
        raise ValueError("target_tokens produced zero tokens")
    if len(tokens) != len(set(tokens)):
        raise ValueError(f"target_tokens contains duplicates: {raw!r}")
    normalized = ";".join(tokens)
    if normalized != raw:
        raise ValueError(f"target_tokens is not normalized: raw={raw!r}, normalized={normalized!r}")
    return tokens


def normalize_text_series(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("").str.strip()


def fail_on_forbidden_cell_line(frame: pd.DataFrame, columns: Sequence[str], name: str) -> None:
    for column in columns:
        if column not in frame.columns:
            continue
        values = normalize_text_series(frame[column])
        bad = frame.loc[values.eq(FORBIDDEN_ANALYSIS_CELL_LINE), column]
        if not bad.empty:
            examples = bad.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
            raise ValueError(f"{name} contains forbidden analysis cell_line {FORBIDDEN_ANALYSIS_CELL_LINE}: {examples}")


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> Tuple[float, str]:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0.0:
        return float("nan"), "cosine_zero_norm_centroid"
    return float(np.dot(a, b) / denom), ""


def pearson_corr(a: np.ndarray, b: np.ndarray) -> Tuple[float, str]:
    a_center = a - np.mean(a)
    b_center = b - np.mean(b)
    denom = float(np.linalg.norm(a_center) * np.linalg.norm(b_center))
    if denom == 0.0:
        return float("nan"), "pcc_zero_variance_centroid"
    return float(np.dot(a_center, b_center) / denom), ""


def squared_distance_matrix(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x_norm = np.sum(x * x, axis=1, keepdims=True)
    y_norm = np.sum(y * y, axis=1, keepdims=True).T
    d2 = x_norm + y_norm - 2.0 * (x @ y.T)
    return np.maximum(d2, 0.0)


def energy_distance_biascorr(x: np.ndarray, y: np.ndarray) -> float:
    n = int(x.shape[0])
    m = int(y.shape[0])
    if n < 2 or m < 2:
        return float("nan")

    dxy = squared_distance_matrix(x, y)
    dxx = squared_distance_matrix(x, x)
    dyy = squared_distance_matrix(y, y)

    delta_xy = float(np.sum(dxy) / (n * m))
    sum_x = float(np.sum(dxx) - np.trace(dxx))
    sum_y = float(np.sum(dyy) - np.trace(dyy))
    sigma_x = sum_x / (n * (n - 1))
    sigma_y = sum_y / (m * (m - 1))
    return float(2.0 * delta_xy - sigma_x - sigma_y)


def build_dataset_bundles(snapshot_root: Path) -> Mapping[str, DatasetBundle]:
    return {
        "LINCS": DatasetBundle(
            dataset="LINCS",
            subtree_relpath="lincs",
            subtree_manifest_relpath="lincs/subtree_manifest.json",
            delta_meta_relpath="lincs/derived/delta_meta.csv",
            membership_relpath="lincs/derived/task2_row_membership.parquet",
            mech_order_relpath="lincs/task2_lincs_pairs.csv",
        ),
        "scPerturb": DatasetBundle(
            dataset="scPerturb",
            subtree_relpath="scperturb_k562",
            subtree_manifest_relpath="scperturb_k562/subtree_manifest.json",
            delta_meta_relpath="scperturb_k562/derived/delta_meta.csv",
            membership_relpath="scperturb_k562/derived/task2_row_membership.parquet",
            mech_order_relpath="scperturb_k562/Common_Targets_K562.csv",
            common_targets_relpath="scperturb_k562/Common_Targets_K562.csv",
        ),
    }


def build_physical_representation_specs(snapshot_root: Path) -> Mapping[str, Mapping[str, RepresentationSpec]]:
    return {
        "LINCS": {
            "Gene": RepresentationSpec(
                dataset="LINCS",
                name="Gene",
                family="derived",
                array_relpath=str((snapshot_root / "lincs/derived/gene_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "lincs/derived/delta_meta.csv").relative_to(snapshot_root)),
                expected_dim=2477,
            ),
            "Pathway": RepresentationSpec(
                dataset="LINCS",
                name="Pathway",
                family="derived",
                array_relpath=str((snapshot_root / "lincs/derived/pathway_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "lincs/derived/delta_meta.csv").relative_to(snapshot_root)),
                expected_dim=50,
            ),
        },
        "scPerturb": {
            "Gene": RepresentationSpec(
                dataset="scPerturb",
                name="Gene",
                family="derived",
                array_relpath=str((snapshot_root / "scperturb_k562/derived/gene_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/derived/delta_meta.csv").relative_to(snapshot_root)),
                expected_dim=8363,
            ),
            "Pathway": RepresentationSpec(
                dataset="scPerturb",
                name="Pathway",
                family="derived",
                array_relpath=str((snapshot_root / "scperturb_k562/derived/pathway_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/derived/delta_meta.csv").relative_to(snapshot_root)),
                expected_dim=50,
            ),
            "scgpt": RepresentationSpec(
                dataset="scPerturb",
                name="scgpt",
                family="fm",
                array_relpath=str((snapshot_root / "scperturb_k562/fm/scgpt/fm_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/fm/scgpt/fm_delta_meta.csv").relative_to(snapshot_root)),
                policy_relpath=str((snapshot_root / "scperturb_k562/fm/scgpt/delta_operator_policy.json").relative_to(snapshot_root)),
                expected_dim=512,
            ),
            "geneformer": RepresentationSpec(
                dataset="scPerturb",
                name="geneformer",
                family="fm",
                array_relpath=str((snapshot_root / "scperturb_k562/fm/geneformer/fm_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/fm/geneformer/fm_delta_meta.csv").relative_to(snapshot_root)),
                policy_relpath=str((snapshot_root / "scperturb_k562/fm/geneformer/delta_operator_policy.json").relative_to(snapshot_root)),
                expected_dim=1152,
            ),
            "scbert": RepresentationSpec(
                dataset="scPerturb",
                name="scbert",
                family="fm",
                array_relpath=str((snapshot_root / "scperturb_k562/fm/scbert/fm_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/fm/scbert/fm_delta_meta.csv").relative_to(snapshot_root)),
                policy_relpath=str((snapshot_root / "scperturb_k562/fm/scbert/delta_operator_policy.json").relative_to(snapshot_root)),
                expected_dim=200,
            ),
            "scfoundation": RepresentationSpec(
                dataset="scPerturb",
                name="scfoundation",
                family="fm",
                array_relpath=str((snapshot_root / "scperturb_k562/fm/scfoundation/fm_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/fm/scfoundation/fm_delta_meta.csv").relative_to(snapshot_root)),
                policy_relpath=str((snapshot_root / "scperturb_k562/fm/scfoundation/delta_operator_policy.json").relative_to(snapshot_root)),
                expected_dim=3072,
            ),
            "uce": RepresentationSpec(
                dataset="scPerturb",
                name="uce",
                family="fm",
                array_relpath=str((snapshot_root / "scperturb_k562/fm/uce/fm_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/fm/uce/fm_delta_meta.csv").relative_to(snapshot_root)),
                policy_relpath=str((snapshot_root / "scperturb_k562/fm/uce/delta_operator_policy.json").relative_to(snapshot_root)),
                expected_dim=1280,
            ),
            "state": RepresentationSpec(
                dataset="scPerturb",
                name="state",
                family="fm",
                array_relpath=str((snapshot_root / "scperturb_k562/fm/state/fm_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/fm/state/fm_delta_meta.csv").relative_to(snapshot_root)),
                policy_relpath=str((snapshot_root / "scperturb_k562/fm/state/delta_operator_policy.json").relative_to(snapshot_root)),
                expected_dim=2058,
            ),
            "tahoe-x1": RepresentationSpec(
                dataset="scPerturb",
                name="tahoe-x1",
                family="fm",
                array_relpath=str((snapshot_root / "scperturb_k562/fm/tahoe-x1/fm_delta.npy").relative_to(snapshot_root)),
                meta_relpath=str((snapshot_root / "scperturb_k562/fm/tahoe-x1/fm_delta_meta.csv").relative_to(snapshot_root)),
                policy_relpath=str((snapshot_root / "scperturb_k562/fm/tahoe-x1/delta_operator_policy.json").relative_to(snapshot_root)),
                expected_dim=512,
            ),
        },
    }


def assert_unique_key(frame: pd.DataFrame, columns: Sequence[str], name: str) -> None:
    if frame.empty:
        return
    dup_mask = frame.duplicated(columns, keep=False)
    if bool(dup_mask.any()):
        examples = frame.loc[dup_mask, list(columns)].head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"{name} contains duplicate keys on {list(columns)}: {examples}")


def load_snapshot_manifest(snapshot_root: Path) -> Tuple[Mapping[str, Any], List[Path]]:
    manifest_path = snapshot_root / "snapshot_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing snapshot manifest: {manifest_path}")
    manifest = read_json(manifest_path)
    if manifest.get("analysis_key_fields") != ["dataset", "cell_line", "target_token"]:
        raise ValueError("snapshot_manifest.json analysis_key_fields must equal ['dataset', 'cell_line', 'target_token']")
    datasets = manifest.get("datasets", {})
    if set(datasets.keys()) != set(DATASET_ORDER):
        raise ValueError(f"snapshot_manifest.json datasets must equal {list(DATASET_ORDER)}")
    for dataset in DATASET_ORDER:
        dataset_info = datasets[dataset]
        if dataset_info.get("representations") != list(FROZEN_AVAILABLE_REPRESENTATIONS[dataset]):
            raise ValueError(f"snapshot_manifest.json representations drift for dataset={dataset}")
        cell_lines = dataset_info.get("cell_lines", [])
        if FORBIDDEN_ANALYSIS_CELL_LINE in cell_lines:
            raise ValueError(f"snapshot_manifest.json contains forbidden analysis cell_line {FORBIDDEN_ANALYSIS_CELL_LINE}")
    return manifest, [manifest_path]


def load_shared_contract_tables(snapshot_root: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[Path]]:
    coverage_path = snapshot_root / "task2_pairs_coverage.csv"
    registry_path = snapshot_root / "representation_availability_registry.csv"
    if not coverage_path.is_file():
        raise FileNotFoundError(f"Missing task2_pairs_coverage.csv: {coverage_path}")
    if not registry_path.is_file():
        raise FileNotFoundError(f"Missing representation_availability_registry.csv: {registry_path}")

    coverage = pd.read_csv(coverage_path)
    ensure_required_columns(
        coverage,
        [
            "dataset",
            "cell_line",
            "target_token",
            "n_chem_instances",
            "n_gen_instances",
            "is_eligible_bool",
            "eligibility_reason",
            "source_tag",
        ],
        "task2_pairs_coverage.csv",
    )
    coverage["dataset"] = normalize_text_series(coverage["dataset"])
    coverage["cell_line"] = normalize_text_series(coverage["cell_line"])
    coverage["target_token"] = normalize_text_series(coverage["target_token"])
    coverage["n_chem_instances"] = pd.to_numeric(coverage["n_chem_instances"], errors="raise").astype(np.int64)
    coverage["n_gen_instances"] = pd.to_numeric(coverage["n_gen_instances"], errors="raise").astype(np.int64)
    coverage["is_eligible_bool"] = parse_bool_series(coverage["is_eligible_bool"])
    fail_on_forbidden_cell_line(coverage, ["cell_line"], "task2_pairs_coverage.csv")
    bad_datasets = coverage.loc[~coverage["dataset"].isin(DATASET_ORDER), "dataset"]
    if not bad_datasets.empty:
        examples = bad_datasets.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"task2_pairs_coverage.csv contains unsupported datasets: {examples}")
    eligible = coverage.loc[coverage["is_eligible_bool"]].copy()
    assert_unique_key(eligible, ["dataset", "cell_line", "target_token"], "eligible task2_pairs_coverage.csv")
    bad_eligible = eligible.loc[
        eligible["n_chem_instances"].lt(1) | eligible["n_gen_instances"].lt(1),
        ["dataset", "cell_line", "target_token", "n_chem_instances", "n_gen_instances"],
    ]
    if not bad_eligible.empty:
        examples = bad_eligible.head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"Eligible coverage rows must have both sides present: {examples}")

    registry = pd.read_csv(registry_path)
    ensure_required_columns(
        registry,
        [
            "dataset",
            "cell_line",
            "representation",
            "availability_status",
            "availability_reason",
        ],
        "representation_availability_registry.csv",
    )
    registry["dataset"] = normalize_text_series(registry["dataset"])
    registry["cell_line"] = normalize_text_series(registry["cell_line"])
    registry["representation"] = normalize_text_series(registry["representation"])
    registry["availability_status"] = normalize_text_series(registry["availability_status"])
    registry["availability_reason"] = normalize_text_series(registry["availability_reason"])
    fail_on_forbidden_cell_line(registry, ["cell_line"], "representation_availability_registry.csv")
    assert_unique_key(registry, ["dataset", "cell_line", "representation"], "representation_availability_registry.csv")

    bad_registry_dataset = registry.loc[~registry["dataset"].isin(DATASET_ORDER), "dataset"]
    if not bad_registry_dataset.empty:
        examples = bad_registry_dataset.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"representation_availability_registry.csv contains unsupported datasets: {examples}")
    bad_registry_rep = registry.loc[~registry["representation"].isin(REPRESENTATION_ORDER), "representation"]
    if not bad_registry_rep.empty:
        examples = bad_registry_rep.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"representation_availability_registry.csv contains unsupported representations: {examples}")
    bad_status = registry.loc[
        ~registry["availability_status"].isin([STATUS_AVAILABLE, STATUS_NOT_APPLICABLE]),
        "availability_status",
    ]
    if not bad_status.empty:
        examples = bad_status.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"representation_availability_registry.csv contains unsupported availability_status values: {examples}")

    eligible_cell_lines = eligible[["dataset", "cell_line"]].drop_duplicates().reset_index(drop=True)
    expected_registry_rows = len(eligible_cell_lines) * len(REPRESENTATION_ORDER)
    registry_for_eligible = eligible_cell_lines.merge(
        registry,
        on=["dataset", "cell_line"],
        how="left",
        sort=False,
        indicator=True,
    )
    missing_pairs = registry_for_eligible.loc[registry_for_eligible["_merge"].ne("both"), ["dataset", "cell_line"]]
    if not missing_pairs.empty:
        examples = missing_pairs.drop_duplicates().head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"Missing availability rows for eligible dataset/cell_line pairs: {examples}")
    registry_for_eligible = registry_for_eligible.drop(columns=["_merge"])
    if len(registry_for_eligible) != expected_registry_rows:
        raise ValueError(
            "representation_availability_registry.csv must contain exactly one row per "
            f"eligible dataset/cell_line/representation; expected={expected_registry_rows}, got={len(registry_for_eligible)}"
        )

    for dataset in DATASET_ORDER:
        dataset_registry = registry_for_eligible.loc[registry_for_eligible["dataset"].eq(dataset)].copy()
        allowed_reps = set(FROZEN_AVAILABLE_REPRESENTATIONS[dataset])
        for row in dataset_registry.itertuples(index=False):
            expected_status = STATUS_AVAILABLE if row.representation in allowed_reps else STATUS_NOT_APPLICABLE
            if row.availability_status != expected_status:
                raise ValueError(
                    "representation_availability_registry.csv scope policy drift for "
                    f"dataset={dataset}, cell_line={row.cell_line}, representation={row.representation}, "
                    f"expected_status={expected_status}, got={row.availability_status}"
                )
            expected_reason = EXPECTED_AVAILABILITY_REASON[expected_status]
            if row.availability_reason != expected_reason:
                raise ValueError(
                    "representation_availability_registry.csv availability_reason drift for "
                    f"dataset={dataset}, cell_line={row.cell_line}, representation={row.representation}, "
                    f"expected_reason={expected_reason}, got={row.availability_reason}"
                )

    available_scope = registry.loc[registry["availability_status"].eq(STATUS_AVAILABLE), [
        "dataset",
        "cell_line",
        "representation",
    ]].copy()
    supported_grid = eligible[[
        "dataset",
        "cell_line",
        "target_token",
        "n_chem_instances",
        "n_gen_instances",
    ]].merge(
        available_scope,
        on=["dataset", "cell_line"],
        how="inner",
        sort=False,
    )
    assert_unique_key(
        supported_grid,
        ["dataset", "cell_line", "target_token", "representation"],
        "supported analysis grid",
    )

    return coverage, eligible, supported_grid, [coverage_path, registry_path]


def load_common_targets(common_targets_path: Path) -> Tuple[List[str], List[Path]]:
    if not common_targets_path.is_file():
        raise FileNotFoundError(f"Missing Common_Targets_K562.csv: {common_targets_path}")
    common_targets = pd.read_csv(common_targets_path)
    ensure_required_columns(common_targets, ["Common_Target"], "Common_Targets_K562.csv")
    target_order: List[str] = []
    seen = set()
    for value in common_targets["Common_Target"].tolist():
        if pd.isna(value):
            raise ValueError("Common_Targets_K562.csv contains NA target")
        token = str(value).strip()
        if not token:
            raise ValueError("Common_Targets_K562.csv contains empty target")
        if token in seen:
            raise ValueError(f"Common_Targets_K562.csv contains duplicate target: {token}")
        seen.add(token)
        target_order.append(token)
    if not target_order:
        raise ValueError("Common_Targets_K562.csv produced empty target order")
    return target_order, [common_targets_path]


def load_scperturb_mech_order(
    common_targets_path: Path,
    eligible_coverage: pd.DataFrame,
) -> Tuple[List[Tuple[str, str, str]], List[Path]]:
    target_order, input_paths = load_common_targets(common_targets_path)
    sc_eligible = eligible_coverage.loc[eligible_coverage["dataset"].eq("scPerturb")].copy()
    if sc_eligible.empty:
        raise ValueError("Eligible scPerturb coverage cannot be empty")
    fail_on_forbidden_cell_line(sc_eligible, ["cell_line"], "eligible scPerturb coverage")
    cell_lines = sc_eligible["cell_line"].drop_duplicates().tolist()
    if cell_lines != ["K562"]:
        raise ValueError(f"scPerturb eligible coverage must contain only cell_line=K562, got={cell_lines}")
    eligible_target_set = set(sc_eligible["target_token"].tolist())
    ordered_targets = [token for token in target_order if token in eligible_target_set]
    if ordered_targets != sc_eligible["target_token"].tolist():
        raise ValueError("scPerturb eligible coverage target order must match Common_Targets_K562.csv file order")
    return [("scPerturb", "K562", token) for token in ordered_targets], input_paths


def load_lincs_mech_order(
    lincs_pairs_path: Path,
    eligible_coverage: pd.DataFrame,
) -> Tuple[List[Tuple[str, str, str]], List[Path]]:
    if not lincs_pairs_path.is_file():
        raise FileNotFoundError(f"Missing task2_lincs_pairs.csv: {lincs_pairs_path}")
    lincs_pairs = pd.read_csv(lincs_pairs_path)
    ensure_required_columns(
        lincs_pairs,
        [
            "dataset",
            "cell_line",
            "target_token",
            "n_chem_instances",
            "n_gen_instances",
            "is_eligible_bool",
            "eligibility_reason",
            "source_tag",
        ],
        "task2_lincs_pairs.csv",
    )
    lincs_pairs["dataset"] = normalize_text_series(lincs_pairs["dataset"])
    lincs_pairs["cell_line"] = normalize_text_series(lincs_pairs["cell_line"])
    lincs_pairs["target_token"] = normalize_text_series(lincs_pairs["target_token"])
    lincs_pairs["n_chem_instances"] = pd.to_numeric(lincs_pairs["n_chem_instances"], errors="raise").astype(np.int64)
    lincs_pairs["n_gen_instances"] = pd.to_numeric(lincs_pairs["n_gen_instances"], errors="raise").astype(np.int64)
    lincs_pairs["is_eligible_bool"] = parse_bool_series(lincs_pairs["is_eligible_bool"])
    fail_on_forbidden_cell_line(lincs_pairs, ["cell_line"], "task2_lincs_pairs.csv")
    assert_unique_key(lincs_pairs, ["dataset", "cell_line", "target_token"], "task2_lincs_pairs.csv")
    if set(lincs_pairs["dataset"].drop_duplicates().tolist()) != {"LINCS"}:
        raise ValueError("task2_lincs_pairs.csv must contain only dataset=LINCS")
    if not bool(lincs_pairs["is_eligible_bool"].all()):
        raise ValueError("task2_lincs_pairs.csv must contain only eligible rows")

    expected = eligible_coverage.loc[eligible_coverage["dataset"].eq("LINCS"), lincs_pairs.columns.tolist()].reset_index(drop=True)
    observed = lincs_pairs.reset_index(drop=True)
    if not observed.equals(expected):
        raise ValueError("task2_lincs_pairs.csv must match eligible LINCS rows from task2_pairs_coverage.csv exactly")

    mech_order = [(row.dataset, row.cell_line, row.target_token) for row in lincs_pairs.itertuples(index=False)]
    return mech_order, [lincs_pairs_path]


def load_delta_meta(
    *,
    dataset: str,
    delta_meta_path: Path,
    eligible_cell_lines: Sequence[str],
) -> Tuple[pd.DataFrame, List[Path]]:
    if not delta_meta_path.is_file():
        raise FileNotFoundError(f"Missing delta_meta.csv for dataset={dataset}: {delta_meta_path}")
    delta_meta = pd.read_csv(delta_meta_path)
    required = [
        "row_id",
        "treated_cell_id",
        "perturbation_class",
        "cell_line",
        "target_raw",
        "time",
        "dose_value",
        "specificity_tier",
        "n_controls_used",
        "dataset_side",
        "target_tokens",
        "seed",
    ]
    if dataset == "LINCS":
        required.extend(["source_row_index", "cell_line_raw", "pert_type_raw", "paired_control_id", "sig_id"])
    ensure_required_columns(delta_meta, required, f"{dataset} delta_meta.csv")

    delta_meta["row_id"] = pd.to_numeric(delta_meta["row_id"], errors="raise").astype(np.int64)
    delta_meta["treated_cell_id"] = delta_meta["treated_cell_id"].astype(str)
    delta_meta["perturbation_class"] = normalize_text_series(delta_meta["perturbation_class"])
    delta_meta["cell_line"] = normalize_text_series(delta_meta["cell_line"])
    delta_meta["target_raw"] = delta_meta["target_raw"].astype("string")
    delta_meta["n_controls_used"] = pd.to_numeric(delta_meta["n_controls_used"], errors="raise").astype(np.int64)
    delta_meta["dataset_side"] = normalize_text_series(delta_meta["dataset_side"])
    delta_meta["seed"] = pd.to_numeric(delta_meta["seed"], errors="raise").astype(np.int64)
    delta_meta = delta_meta.sort_values("row_id", kind="mergesort").reset_index(drop=True)

    fail_on_forbidden_cell_line(delta_meta, ["cell_line"], f"{dataset} delta_meta.csv")
    row_ids = delta_meta["row_id"].to_numpy(dtype=np.int64, copy=False)
    if not np.array_equal(row_ids, np.arange(len(delta_meta), dtype=np.int64)):
        raise ValueError(f"{dataset} delta_meta.row_id must be contiguous 0..N-1")
    if delta_meta["treated_cell_id"].duplicated().any():
        dupes = delta_meta.loc[delta_meta["treated_cell_id"].duplicated(), "treated_cell_id"].head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"{dataset} delta_meta treated_cell_id must be unique, examples={dupes}")
    if not delta_meta["seed"].eq(GLOBAL_SEED).all():
        raise ValueError(f"{dataset} delta_meta seed must equal {GLOBAL_SEED}")
    if delta_meta["n_controls_used"].lt(1).any():
        examples = delta_meta.loc[
            delta_meta["n_controls_used"].lt(1),
            ["row_id", "n_controls_used"],
        ].head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"{dataset} delta_meta n_controls_used must be >=1, examples={examples}")
    bad_class = delta_meta.loc[
        ~delta_meta["perturbation_class"].isin(SIDE_ORDER),
        "perturbation_class",
    ]
    if not bad_class.empty:
        examples = bad_class.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"{dataset} delta_meta contains unsupported perturbation_class values: {examples}")

    eligible_cell_line_set = set(eligible_cell_lines)
    delta_cell_line_set = set(delta_meta["cell_line"].drop_duplicates().tolist())
    if delta_cell_line_set != eligible_cell_line_set:
        raise ValueError(
            f"{dataset} delta_meta cell_line set must equal eligible coverage cell_line set: "
            f"expected={sorted(eligible_cell_line_set)}, got={sorted(delta_cell_line_set)}"
        )

    parsed_tokens: List[Tuple[str, ...]] = [parse_target_tokens_field(value) for value in delta_meta["target_tokens"].tolist()]
    delta_meta["parsed_target_tokens"] = parsed_tokens
    bad_genetic = delta_meta.loc[
        delta_meta["perturbation_class"].eq("Genetic")
        & delta_meta["parsed_target_tokens"].map(len).ne(1),
        ["row_id", "target_tokens"],
    ]
    if not bad_genetic.empty:
        examples = bad_genetic.head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"{dataset} Genetic rows must have exactly one target token, examples={examples}")

    expected_dataset_side = {"LINCS"} if dataset == "LINCS" else {"CRISPR", "DRUG"}
    actual_dataset_side = set(delta_meta["dataset_side"].drop_duplicates().tolist())
    if actual_dataset_side != expected_dataset_side:
        raise ValueError(
            f"{dataset} delta_meta dataset_side values drift: expected={sorted(expected_dataset_side)}, got={sorted(actual_dataset_side)}"
        )

    return delta_meta, [delta_meta_path]


def build_membership_from_delta_meta(dataset: str, delta_meta: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for row in delta_meta.itertuples(index=False):
        tokens = row.parsed_target_tokens
        base = {
            "row_id": int(row.row_id),
            "treated_cell_id": row.treated_cell_id,
            "dataset": dataset,
            "cell_line": row.cell_line,
            "perturbation_class": row.perturbation_class,
            "target_raw": row.target_raw,
        }
        if dataset == "LINCS":
            base["source_row_index"] = int(row.source_row_index)
            base["cell_line_raw"] = row.cell_line_raw
        for token in tokens if row.perturbation_class == "Chemical" else (tokens[0],):
            out_row = dict(base)
            out_row["target_token"] = token
            rows.append(out_row)

    membership = pd.DataFrame(rows)
    if dataset == "LINCS":
        membership = membership[
            [
                "row_id",
                "source_row_index",
                "treated_cell_id",
                "cell_line_raw",
                "dataset",
                "cell_line",
                "perturbation_class",
                "target_raw",
                "target_token",
            ]
        ]
    else:
        membership = membership[
            [
                "row_id",
                "treated_cell_id",
                "dataset",
                "cell_line",
                "perturbation_class",
                "target_raw",
                "target_token",
            ]
        ]
    membership = membership.sort_values(["row_id", "target_token"], kind="mergesort").reset_index(drop=True)
    return membership


def load_membership_parity_table(dataset: str, membership_path: Path) -> Tuple[pd.DataFrame, List[Path]]:
    if not membership_path.is_file():
        raise FileNotFoundError(f"Missing task2_row_membership.parquet for dataset={dataset}: {membership_path}")
    membership = pd.read_parquet(membership_path)
    required = ["row_id", "treated_cell_id", "dataset", "cell_line", "perturbation_class", "target_raw", "target_token"]
    if dataset == "LINCS":
        required.extend(["source_row_index", "cell_line_raw"])
    ensure_required_columns(membership, required, f"{dataset} task2_row_membership.parquet")
    membership["row_id"] = pd.to_numeric(membership["row_id"], errors="raise").astype(np.int64)
    membership["treated_cell_id"] = membership["treated_cell_id"].astype(str)
    membership["dataset"] = normalize_text_series(membership["dataset"])
    membership["cell_line"] = normalize_text_series(membership["cell_line"])
    membership["perturbation_class"] = normalize_text_series(membership["perturbation_class"])
    membership["target_token"] = normalize_text_series(membership["target_token"])
    fail_on_forbidden_cell_line(membership, ["cell_line"], f"{dataset} task2_row_membership.parquet")
    if dataset == "LINCS":
        membership["source_row_index"] = pd.to_numeric(membership["source_row_index"], errors="raise").astype(np.int64)
        membership["cell_line_raw"] = membership["cell_line_raw"].astype(str)
    membership = membership.sort_values(["row_id", "target_token"], kind="mergesort").reset_index(drop=True)
    return membership, [membership_path]


def validate_membership_parity(
    *,
    dataset: str,
    exploded_membership: pd.DataFrame,
    parity_membership: pd.DataFrame,
) -> None:
    compare_columns = parity_membership.columns.tolist()
    left = exploded_membership[compare_columns].reset_index(drop=True)
    right = parity_membership[compare_columns].reset_index(drop=True)
    if len(left) != len(right):
        raise ValueError(
            f"{dataset} task2_row_membership.parquet parity row-count mismatch: exploded={len(left)}, parity={len(right)}"
        )

    mismatch_rows: List[Dict[str, object]] = []
    for idx in range(len(left)):
        row_match = True
        for column in compare_columns:
            left_value = left.iloc[idx][column]
            right_value = right.iloc[idx][column]
            if pd.isna(left_value) and pd.isna(right_value):
                continue
            if left_value != right_value:
                row_match = False
                break
        if not row_match:
            mismatch_rows.append(
                {
                    "index": int(idx),
                    "exploded": left.iloc[idx].to_dict(),
                    "parity": right.iloc[idx].to_dict(),
                }
            )
            if len(mismatch_rows) >= MAX_COUNTEREXAMPLES:
                break

    if mismatch_rows:
        raise ValueError(f"{dataset} task2_row_membership.parquet parity check failed: {mismatch_rows}")


def build_base_memberships(
    *,
    dataset: str,
    mech_order: Sequence[Tuple[str, str, str]],
    membership: pd.DataFrame,
    eligible_coverage_subset: pd.DataFrame,
) -> Tuple[Dict[Tuple[str, str], np.ndarray], Dict[Tuple[str, str], np.ndarray]]:
    if dataset == "LINCS":
        eligible_keys = {(cell_line, target_token) for _, cell_line, target_token in mech_order}
    else:
        eligible_keys = {(cell_line, target_token) for _, cell_line, target_token in mech_order}

    chem_rows: MutableMapping[Tuple[str, str], List[int]] = {key: [] for key in eligible_keys}
    gen_rows: MutableMapping[Tuple[str, str], List[int]] = {key: [] for key in eligible_keys}
    for row in membership.itertuples(index=False):
        key = (row.cell_line, row.target_token)
        if key not in eligible_keys:
            continue
        if row.perturbation_class == "Chemical":
            chem_rows[key].append(int(row.row_id))
        else:
            gen_rows[key].append(int(row.row_id))

    chem_arrays = {key: np.asarray(values, dtype=np.int64) for key, values in chem_rows.items()}
    gen_arrays = {key: np.asarray(values, dtype=np.int64) for key, values in gen_rows.items()}

    coverage_keyed = eligible_coverage_subset.set_index(["cell_line", "target_token"])
    for _, cell_line, target_token in mech_order:
        key = (cell_line, target_token)
        expected_chem = int(coverage_keyed.loc[key, "n_chem_instances"])
        expected_gen = int(coverage_keyed.loc[key, "n_gen_instances"])
        actual_chem = int(chem_arrays[key].size)
        actual_gen = int(gen_arrays[key].size)
        if actual_chem != expected_chem or actual_gen != expected_gen:
            raise ValueError(
                f"{dataset} membership count drift for cell_line={cell_line}, target_token={target_token}: "
                f"expected=({expected_chem},{expected_gen}) got=({actual_chem},{actual_gen})"
            )

    return chem_arrays, gen_arrays


def batched_mask_contract_check(arr: np.ndarray, valid_mask: np.ndarray, batch_rows: int) -> List[str]:
    errors: List[str] = []
    n_rows = int(arr.shape[0])
    for start in range(0, n_rows, batch_rows):
        end = min(start + batch_rows, n_rows)
        block = np.asarray(arr[start:end])
        block_valid = valid_mask[start:end]

        all_finite = np.isfinite(block).all(axis=1)
        all_nan = np.isnan(block).all(axis=1)

        if bool(block_valid.any()):
            bad_valid = np.flatnonzero(block_valid & ~all_finite)
            if bad_valid.size > 0:
                errors.append(
                    f"valid rows contain non-finite values; first_row_id={start + int(bad_valid[0])}"
                )
                break

        invalid_mask = ~block_valid
        if bool(invalid_mask.any()):
            bad_invalid = np.flatnonzero(invalid_mask & ~all_nan)
            if bad_invalid.size > 0:
                errors.append(
                    f"invalid rows are not all-NaN; first_row_id={start + int(bad_invalid[0])}"
                )
                break
    return errors


def validate_representation_inputs(
    *,
    snapshot_root: Path,
    spec: RepresentationSpec,
    delta_meta: pd.DataFrame,
) -> Tuple[RepresentationPayload, List[Path], Dict[str, Any]]:
    array_path = snapshot_root / spec.array_relpath
    meta_path = snapshot_root / spec.meta_relpath
    required_inputs = [array_path, meta_path]

    if not array_path.is_file():
        raise FileNotFoundError(f"Missing required array for {spec.dataset}/{spec.name}: {array_path}")
    if not meta_path.is_file():
        raise FileNotFoundError(f"Missing required metadata for {spec.dataset}/{spec.name}: {meta_path}")

    arr = np.load(array_path, mmap_mode="r")
    if arr.ndim != 2:
        raise ValueError(f"{spec.dataset}/{spec.name} array must be 2D, got shape={tuple(arr.shape)}")
    n_rows = int(len(delta_meta))
    if int(arr.shape[0]) != n_rows:
        raise ValueError(
            f"{spec.dataset}/{spec.name} row-count mismatch: array_rows={arr.shape[0]}, delta_meta_rows={n_rows}"
        )
    if int(arr.shape[1]) != int(spec.expected_dim):
        raise ValueError(
            f"{spec.dataset}/{spec.name} dimension mismatch: expected={spec.expected_dim}, got={arr.shape[1]}"
        )

    if spec.family == "derived":
        valid_mask = np.ones(n_rows, dtype=bool)
        invalid_by_side: Dict[str, int] = {}
    else:
        meta = pd.read_csv(meta_path)
        ensure_required_columns(
            meta,
            ["row_id", "treated_cell_id", "valid_mask", "n_controls_used", "dataset_side", "model_name", "seed"],
            f"{spec.dataset}/{spec.name} fm_delta_meta.csv",
        )
        meta = meta.sort_values("row_id", kind="mergesort").reset_index(drop=True)
        meta["row_id"] = pd.to_numeric(meta["row_id"], errors="raise").astype(np.int64)
        meta["treated_cell_id"] = meta["treated_cell_id"].astype(str)
        meta["n_controls_used"] = pd.to_numeric(meta["n_controls_used"], errors="raise").astype(np.int64)
        meta["dataset_side"] = normalize_text_series(meta["dataset_side"])
        meta["model_name"] = normalize_text_series(meta["model_name"])
        meta["seed"] = pd.to_numeric(meta["seed"], errors="raise").astype(np.int64)

        if not np.array_equal(meta["row_id"].to_numpy(dtype=np.int64, copy=False), delta_meta["row_id"].to_numpy(dtype=np.int64, copy=False)):
            raise ValueError(f"{spec.dataset}/{spec.name} fm_delta_meta row_id mismatch against delta_meta.csv")
        if not np.array_equal(meta["treated_cell_id"].to_numpy(dtype=object), delta_meta["treated_cell_id"].to_numpy(dtype=object)):
            raise ValueError(f"{spec.dataset}/{spec.name} treated_cell_id mismatch against delta_meta.csv")
        if not np.array_equal(
            meta["n_controls_used"].to_numpy(dtype=np.int64, copy=False),
            delta_meta["n_controls_used"].to_numpy(dtype=np.int64, copy=False),
        ):
            raise ValueError(f"{spec.dataset}/{spec.name} n_controls_used mismatch against delta_meta.csv")
        if not np.array_equal(meta["dataset_side"].to_numpy(dtype=object), delta_meta["dataset_side"].to_numpy(dtype=object)):
            raise ValueError(f"{spec.dataset}/{spec.name} dataset_side mismatch against delta_meta.csv")
        if not meta["model_name"].eq(spec.name).all():
            bad_models = meta.loc[~meta["model_name"].eq(spec.name), "model_name"].drop_duplicates().tolist()
            raise ValueError(f"{spec.dataset}/{spec.name} model_name mismatch in fm_delta_meta.csv: {bad_models[:MAX_COUNTEREXAMPLES]}")
        if not meta["seed"].eq(GLOBAL_SEED).all():
            raise ValueError(f"{spec.dataset}/{spec.name} fm_delta_meta seed must equal {GLOBAL_SEED}")

        valid_mask = parse_bool_series(meta["valid_mask"])
        invalid_meta = meta.loc[~valid_mask, "dataset_side"]
        invalid_by_side = {str(key): int(value) for key, value in invalid_meta.value_counts().to_dict().items()}

    contract_errors = batched_mask_contract_check(arr, valid_mask, MASK_BATCH_ROWS)
    if contract_errors:
        raise ValueError(f"{spec.dataset}/{spec.name} finite-mask contract violation: {contract_errors[0]}")

    policy_path: Optional[Path] = None
    if spec.policy_relpath is not None:
        candidate = snapshot_root / spec.policy_relpath
        if not candidate.is_file():
            raise FileNotFoundError(f"Missing required policy file for {spec.dataset}/{spec.name}: {candidate}")
        policy_path = candidate
        required_inputs.append(candidate)

    del arr

    details = {
        "dataset": spec.dataset,
        "representation": spec.name,
        "family": spec.family,
        "array_path": str(array_path),
        "meta_path": str(meta_path),
        "policy_path": str(policy_path) if policy_path is not None else "",
        "policy_present": bool(policy_path is not None),
        "n_rows": n_rows,
        "dim": int(spec.expected_dim),
        "n_valid": int(valid_mask.sum()),
        "n_invalid": int((~valid_mask).sum()),
        "invalid_by_side": invalid_by_side,
    }

    return (
        RepresentationPayload(
            spec=spec,
            array_path=array_path,
            meta_path=meta_path,
            valid_mask=valid_mask,
            policy_path=policy_path,
            n_valid=int(valid_mask.sum()),
            n_invalid=int((~valid_mask).sum()),
            invalid_by_side=invalid_by_side,
        ),
        required_inputs,
        details,
    )


def deterministic_subsample_row_ids(
    row_ids: np.ndarray,
    max_n: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, int]:
    if row_ids.size <= max_n:
        return row_ids.copy(), 0
    sorted_row_ids = np.sort(row_ids.astype(np.int64, copy=False), kind="mergesort")
    perm = rng.permutation(sorted_row_ids.size)
    return sorted_row_ids[perm[:max_n]], 1


def make_attrition_row(
    *,
    dataset: str,
    cell_line: str,
    target_token: str,
    representation: str,
    metric_name: str,
    reason: str,
    n_chem_before: int,
    n_chem_after: int,
    n_gen_before: int,
    n_gen_after: int,
    notes: str,
) -> Dict[str, object]:
    return {
        "dataset": dataset,
        "cell_line": cell_line,
        "target_token": target_token,
        "representation": representation,
        "metric_name": metric_name,
        "reason": reason,
        "n_chem_before": int(n_chem_before),
        "n_chem_after": int(n_chem_after),
        "n_gen_before": int(n_gen_before),
        "n_gen_after": int(n_gen_after),
        "n_chem_removed": int(n_chem_before - n_chem_after),
        "n_gen_removed": int(n_gen_before - n_gen_after),
        "notes": notes,
    }


def prepare_group_tasks_for_representation(
    *,
    dataset: str,
    representation: str,
    mech_keys: Sequence[Tuple[str, str]],
    chem_rows_by_key: Mapping[Tuple[str, str], np.ndarray],
    gen_rows_by_key: Mapping[Tuple[str, str], np.ndarray],
    valid_mask: np.ndarray,
    rng: np.random.Generator,
) -> Tuple[List[GroupTask], List[Dict[str, object]], int]:
    tasks: List[GroupTask] = []
    attrition_rows: List[Dict[str, object]] = []
    random_draws = 0

    for cell_line, target_token in mech_keys:
        key = (cell_line, target_token)
        chem_total = chem_rows_by_key[key]
        gen_total = gen_rows_by_key[key]

        chem_used = chem_total[valid_mask[chem_total]] if chem_total.size > 0 else np.empty((0,), dtype=np.int64)
        gen_used = gen_total[valid_mask[gen_total]] if gen_total.size > 0 else np.empty((0,), dtype=np.int64)

        if chem_used.size != chem_total.size or gen_used.size != gen_total.size:
            attrition_rows.append(
                make_attrition_row(
                    dataset=dataset,
                    cell_line=cell_line,
                    target_token=target_token,
                    representation=representation,
                    metric_name="ALL",
                    reason="representation_valid_mask_drop",
                    n_chem_before=int(chem_total.size),
                    n_chem_after=int(chem_used.size),
                    n_gen_before=int(gen_total.size),
                    n_gen_after=int(gen_used.size),
                    notes="Invalid rows excluded via representation valid_mask before metric computation",
                )
            )

        chem_sub, chem_draw = deterministic_subsample_row_ids(chem_used, EDIST_MAX_N, rng)
        gen_sub, gen_draw = deterministic_subsample_row_ids(gen_used, EDIST_MAX_N, rng)
        random_draws += chem_draw + gen_draw

        tasks.append(
            GroupTask(
                dataset=dataset,
                cell_line=cell_line,
                target_token=target_token,
                representation=representation,
                group_id=f"{dataset}||{cell_line}||{target_token}",
                n_chem_total=int(chem_total.size),
                n_gen_total=int(gen_total.size),
                chem_used_row_ids=chem_used.astype(np.int64, copy=False),
                gen_used_row_ids=gen_used.astype(np.int64, copy=False),
                chem_sub_row_ids=chem_sub.astype(np.int64, copy=False),
                gen_sub_row_ids=gen_sub.astype(np.int64, copy=False),
            )
        )

    return tasks, attrition_rows, int(random_draws)


def init_worker(array_path: Path, expected_shape: Tuple[int, int], representation: str) -> None:
    global _WORKER_ARRAY, _WORKER_REPRESENTATION
    _WORKER_ARRAY = np.load(array_path, mmap_mode="r")
    _WORKER_REPRESENTATION = representation
    if _WORKER_ARRAY.shape != expected_shape:
        raise ValueError(
            f"Worker shape mismatch for {representation}: expected={expected_shape}, got={_WORKER_ARRAY.shape}"
        )


def worker_load_rows(row_ids: np.ndarray) -> np.ndarray:
    if _WORKER_ARRAY is None:
        raise RuntimeError("Worker array is not initialized")
    if row_ids.size == 0:
        return np.empty((0, int(_WORKER_ARRAY.shape[1])), dtype=np.float64)
    mat = _WORKER_ARRAY[row_ids.astype(np.int64, copy=False)]
    out = np.asarray(mat, dtype=np.float64)
    if not np.isfinite(out).all():
        raise ValueError(
            f"Non-finite values reached metric computation for representation={_WORKER_REPRESENTATION}"
        )
    return out


def worker_group_metric(task: GroupTask) -> Dict[str, Any]:
    n_chem_used = int(task.chem_used_row_ids.size)
    n_gen_used = int(task.gen_used_row_ids.size)
    n_chem_sub = int(task.chem_sub_row_ids.size)
    n_gen_sub = int(task.gen_sub_row_ids.size)

    row: Dict[str, object] = {
        "dataset": task.dataset,
        "cell_line": task.cell_line,
        "target_token": task.target_token,
        "group_id": task.group_id,
        "representation": task.representation,
        "n_chem_instances_total": int(task.n_chem_total),
        "n_gen_instances_total": int(task.n_gen_total),
        "n_chem_instances_used": n_chem_used,
        "n_gen_instances_used": n_gen_used,
        "n_chem_sub": n_chem_sub,
        "n_gen_sub": n_gen_sub,
        "cosine_centroid": np.nan,
        "pcc_centroid": np.nan,
        "edist_biascorr": np.nan,
        "cosine_valid_bool": False,
        "pcc_valid_bool": False,
        "edist_valid_bool": False,
        "cosine_na_reason": "",
        "pcc_na_reason": "",
        "edist_na_reason": "",
    }

    events: List[Dict[str, object]] = []
    if n_chem_used == 0 or n_gen_used == 0:
        row["cosine_na_reason"] = "missing_one_side"
        row["pcc_na_reason"] = "missing_one_side"
        row["edist_na_reason"] = "missing_one_side"
        events.append(
            make_attrition_row(
                dataset=task.dataset,
                cell_line=task.cell_line,
                target_token=task.target_token,
                representation=task.representation,
                metric_name="ALL",
                reason="missing_one_side",
                n_chem_before=n_chem_used,
                n_chem_after=n_chem_used,
                n_gen_before=n_gen_used,
                n_gen_after=n_gen_used,
                notes="Chemical and Genetic cohorts must both have >=1 valid rows after valid_mask filtering",
            )
        )
        return {"row": row, "events": events}

    chem_used = worker_load_rows(task.chem_used_row_ids)
    gen_used = worker_load_rows(task.gen_used_row_ids)

    chem_centroid = np.mean(chem_used, axis=0, dtype=np.float64)
    gen_centroid = np.mean(gen_used, axis=0, dtype=np.float64)

    cosine, cosine_reason = cosine_similarity(chem_centroid, gen_centroid)
    row["cosine_centroid"] = float(cosine) if np.isfinite(cosine) else np.nan
    row["cosine_valid_bool"] = bool(np.isfinite(cosine))
    row["cosine_na_reason"] = cosine_reason
    if cosine_reason:
        events.append(
            make_attrition_row(
                dataset=task.dataset,
                cell_line=task.cell_line,
                target_token=task.target_token,
                representation=task.representation,
                metric_name="cosine_centroid",
                reason=cosine_reason,
                n_chem_before=n_chem_used,
                n_chem_after=n_chem_used,
                n_gen_before=n_gen_used,
                n_gen_after=n_gen_used,
                notes="Centroid cosine undefined because at least one centroid has zero norm",
            )
        )

    pcc, pcc_reason = pearson_corr(chem_centroid, gen_centroid)
    row["pcc_centroid"] = float(pcc) if np.isfinite(pcc) else np.nan
    row["pcc_valid_bool"] = bool(np.isfinite(pcc))
    row["pcc_na_reason"] = pcc_reason
    if pcc_reason:
        events.append(
            make_attrition_row(
                dataset=task.dataset,
                cell_line=task.cell_line,
                target_token=task.target_token,
                representation=task.representation,
                metric_name="pcc_centroid",
                reason=pcc_reason,
                n_chem_before=n_chem_used,
                n_chem_after=n_chem_used,
                n_gen_before=n_gen_used,
                n_gen_after=n_gen_used,
                notes="Centroid Pearson correlation undefined because at least one centered centroid has zero variance",
            )
        )

    if n_chem_sub < 2 or n_gen_sub < 2:
        row["edist_na_reason"] = "edist_insufficient_cells"
        events.append(
            make_attrition_row(
                dataset=task.dataset,
                cell_line=task.cell_line,
                target_token=task.target_token,
                representation=task.representation,
                metric_name="edist_biascorr",
                reason="edist_insufficient_cells",
                n_chem_before=n_chem_used,
                n_chem_after=n_chem_sub,
                n_gen_before=n_gen_used,
                n_gen_after=n_gen_sub,
                notes="Bias-corrected energy distance requires >=2 rows per side after deterministic subsampling",
            )
        )
        return {"row": row, "events": events}

    chem_sub = worker_load_rows(task.chem_sub_row_ids)
    gen_sub = worker_load_rows(task.gen_sub_row_ids)
    edist = energy_distance_biascorr(chem_sub, gen_sub)
    if np.isfinite(edist):
        row["edist_biascorr"] = float(edist)
        row["edist_valid_bool"] = True
    else:
        row["edist_na_reason"] = "edist_non_finite"
        events.append(
            make_attrition_row(
                dataset=task.dataset,
                cell_line=task.cell_line,
                target_token=task.target_token,
                representation=task.representation,
                metric_name="edist_biascorr",
                reason="edist_non_finite",
                n_chem_before=n_chem_sub,
                n_chem_after=n_chem_sub,
                n_gen_before=n_gen_sub,
                n_gen_after=n_gen_sub,
                notes="Bias-corrected energy distance returned non-finite value despite finite inputs",
            )
        )

    return {"row": row, "events": events}


def run_group_tasks(
    *,
    tasks: Sequence[GroupTask],
    array_path: Path,
    expected_shape: Tuple[int, int],
    worker_name: str,
    workers: int,
) -> Tuple[List[Dict[str, Any]], str]:
    global _PROCESS_POOL_PERMISSION_DENIED

    if not tasks:
        return [], "none"

    max_workers = min(workers, max(1, len(tasks)))
    if max_workers == 1:
        init_worker(array_path, expected_shape, worker_name)
        return [worker_group_metric(task) for task in tasks], "serial"

    if not _PROCESS_POOL_PERMISSION_DENIED:
        try:
            with ProcessPoolExecutor(
                max_workers=max_workers,
                initializer=init_worker,
                initargs=(array_path, expected_shape, worker_name),
            ) as executor:
                return list(executor.map(worker_group_metric, tasks, chunksize=32)), "process"
        except PermissionError:
            _PROCESS_POOL_PERMISSION_DENIED = True
        except OSError as exc:
            if getattr(exc, "errno", None) != 13:
                raise
            _PROCESS_POOL_PERMISSION_DENIED = True

    with ThreadPoolExecutor(
        max_workers=max_workers,
        initializer=init_worker,
        initargs=(array_path, expected_shape, worker_name),
    ) as executor:
        return list(executor.map(worker_group_metric, tasks, chunksize=32)), "thread"


def canonicalize_concordance_frame(
    frame: pd.DataFrame,
    mech_order: Sequence[Tuple[str, str, str]],
    representation_order: Sequence[str],
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=CONCORDANCE_COLUMNS)
    mech_index = {key: idx for idx, key in enumerate(mech_order)}
    rep_index = {name: idx for idx, name in enumerate(representation_order)}
    out = frame.copy()
    out["_mech_idx"] = list(zip(out["dataset"], out["cell_line"], out["target_token"]))
    missing_keys = [key for key in out["_mech_idx"].tolist() if key not in mech_index]
    if missing_keys:
        raise ValueError(f"Concordance frame contains unsupported mech keys: {missing_keys[:MAX_COUNTEREXAMPLES]}")
    out["_mech_idx"] = out["_mech_idx"].map(mech_index)
    out["_rep_idx"] = out["representation"].map(rep_index)
    out = out.sort_values(["_mech_idx", "_rep_idx"], kind="mergesort").drop(columns=["_mech_idx", "_rep_idx"])
    return out.reset_index(drop=True)[CONCORDANCE_COLUMNS]


def canonicalize_attrition_frame(
    frame: pd.DataFrame,
    mech_order: Sequence[Tuple[str, str, str]],
    representation_order: Sequence[str],
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=ATTRITION_COLUMNS)
    mech_index = {key: idx for idx, key in enumerate(mech_order)}
    rep_index = {name: idx for idx, name in enumerate(representation_order)}
    metric_index = {name: idx for idx, name in enumerate(METRIC_ORDER)}
    out = frame.copy()
    out["_mech_idx"] = list(zip(out["dataset"], out["cell_line"], out["target_token"]))
    missing_keys = [key for key in out["_mech_idx"].tolist() if key not in mech_index]
    if missing_keys:
        raise ValueError(f"Attrition frame contains unsupported mech keys: {missing_keys[:MAX_COUNTEREXAMPLES]}")
    out["_mech_idx"] = out["_mech_idx"].map(mech_index)
    out["_rep_idx"] = out["representation"].map(rep_index)
    out["_metric_idx"] = out["metric_name"].map(metric_index).fillna(len(METRIC_ORDER))
    out = out.sort_values(["_mech_idx", "_rep_idx", "_metric_idx", "reason"], kind="mergesort").drop(
        columns=["_mech_idx", "_rep_idx", "_metric_idx"]
    )
    return out.reset_index(drop=True)[ATTRITION_COLUMNS]


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


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()

    config_path = (project_root / CONFIG_PATH).resolve()
    if not config_path.is_file():
        print(f"[ERROR] Missing config file: {config_path}", file=sys.stderr)
        return 2

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    configured_seed = int(config["project"]["seed"])
    seed = configured_seed if args.seed is None else int(args.seed)
    if seed != GLOBAL_SEED:
        print(f"[ERROR] Seed must be GLOBAL_SEED={GLOBAL_SEED}, got {seed}", file=sys.stderr)
        return 3

    workers = int(args.workers) if args.workers is not None else int(config["compute"]["workers"])
    if workers < 1:
        print(f"[ERROR] --workers must be >=1, got {workers}", file=sys.stderr)
        return 4

    runs_dir = resolve_config_path(project_root, str(config["paths"]["runs_dir"]))
    task2_snapshot = (project_root / EXPECTED_TASK2_SNAPSHOT).resolve()
    expected_snapshot = (project_root / EXPECTED_TASK2_SNAPSHOT).resolve()
    if task2_snapshot != expected_snapshot:
        print("[ERROR] Snapshot isolation check failed unexpectedly.", file=sys.stderr)
        return 5
    if not task2_snapshot.is_dir():
        print(f"[ERROR] Missing corrected Task2 snapshot root: {task2_snapshot}", file=sys.stderr)
        return 6

    stage_dir = runs_dir / args.run_id / STAGE
    stage_dir.mkdir(parents=True, exist_ok=True)

    started_at = utc_now_iso()
    rng = init_global_rng(GLOBAL_SEED)

    assertions: List[Dict[str, object]] = []
    input_paths: List[Path] = []

    assertions.append(
        {
            "name": "global_seed_locked",
            "pass": True,
            "details": {
                "seed": GLOBAL_SEED,
                "edist_max_n": EDIST_MAX_N,
                "dataset_order": list(DATASET_ORDER),
                "representation_order": list(REPRESENTATION_ORDER),
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "task2_snapshot_v2_isolation",
            "pass": True,
            "details": {
                "rules": ["S4 reads exclusively from data/task2_snapshot_v2/"],
                "task2_snapshot": str(task2_snapshot),
            },
            "counterexamples": [],
        }
    )

    try:
        snapshot_manifest, shared_manifest_inputs = load_snapshot_manifest(task2_snapshot)
        coverage, eligible_coverage, supported_grid, shared_table_inputs = load_shared_contract_tables(task2_snapshot)
    except Exception as exc:
        assertions.append(
            {
                "name": "shared_snapshot_contracts_ready",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed loading shared snapshot contracts: {exc}", file=sys.stderr)
        return 7

    input_paths.extend(shared_manifest_inputs)
    input_paths.extend(shared_table_inputs)

    assertions.append(
        {
            "name": "dataset_order_frozen",
            "pass": True,
            "details": {"dataset_order": list(DATASET_ORDER)},
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "snapshot_manifest_contract",
            "pass": True,
            "details": {
                "analysis_key_fields": snapshot_manifest["analysis_key_fields"],
                "datasets": {
                    dataset: snapshot_manifest["datasets"][dataset]
                    for dataset in DATASET_ORDER
                },
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "eligible_coverage_contract",
            "pass": True,
            "details": {
                "n_eligible_rows": int(len(eligible_coverage)),
                "eligible_by_dataset": {
                    str(key): int(value)
                    for key, value in eligible_coverage.groupby("dataset", sort=False).size().to_dict().items()
                },
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "availability_registry_scope_policy",
            "pass": True,
            "details": {
                "available_representations": {
                    dataset: list(FROZEN_AVAILABLE_REPRESENTATIONS[dataset])
                    for dataset in DATASET_ORDER
                },
            },
            "counterexamples": [],
        }
    )

    expected_supported_rows = int(len(supported_grid))
    supported_breakdown = (
        supported_grid.groupby(["dataset", "representation"], sort=False).size().reset_index(name="n_rows")
    )
    assertions.append(
        {
            "name": "supported_analysis_grid_contract",
            "pass": True,
            "details": {
                "rules": [
                    "expected_s4_rows := count(inner_join(eligible_mech_keys, available_rep_scope, on=['dataset','cell_line']))",
                ],
                "expected_s4_rows": expected_supported_rows,
                "breakdown": supported_breakdown.to_dict(orient="records"),
            },
            "counterexamples": [],
        }
    )

    dataset_bundles = build_dataset_bundles(task2_snapshot)
    physical_specs = build_physical_representation_specs(task2_snapshot)

    mech_order_all: List[Tuple[str, str, str]] = []
    concordance_rows: List[Dict[str, object]] = []
    attrition_rows: List[Dict[str, object]] = []
    representation_details: List[Dict[str, Any]] = []
    total_random_draws = 0
    parallel_backends_used: List[str] = []

    try:
        for dataset in DATASET_ORDER:
            bundle = dataset_bundles[dataset]
            bundle_inputs = [
                task2_snapshot / bundle.subtree_manifest_relpath,
                task2_snapshot / bundle.delta_meta_relpath,
                task2_snapshot / bundle.membership_relpath,
                task2_snapshot / bundle.mech_order_relpath,
            ]
            for path in bundle_inputs:
                if not path.is_file():
                    raise FileNotFoundError(f"Missing required dataset input for {dataset}: {path}")
            input_paths.extend(bundle_inputs)

            subtree_manifest = read_json(task2_snapshot / bundle.subtree_manifest_relpath)
            if subtree_manifest.get("analysis_key_fields") != ["dataset", "cell_line", "target_token"]:
                raise ValueError(f"{dataset} subtree_manifest analysis_key_fields drift")
            if subtree_manifest.get("representation_scope") != list(FROZEN_AVAILABLE_REPRESENTATIONS[dataset]):
                raise ValueError(f"{dataset} subtree_manifest representation_scope drift")

            eligible_subset = eligible_coverage.loc[eligible_coverage["dataset"].eq(dataset)].copy().reset_index(drop=True)
            eligible_cell_lines = eligible_subset["cell_line"].drop_duplicates().tolist()

            if dataset == "LINCS":
                mech_order_dataset, order_inputs = load_lincs_mech_order(task2_snapshot / bundle.mech_order_relpath, eligible_coverage)
            else:
                mech_order_dataset, order_inputs = load_scperturb_mech_order(task2_snapshot / bundle.common_targets_relpath, eligible_coverage)
            input_paths.extend(order_inputs)
            mech_order_all.extend(mech_order_dataset)

            delta_meta, delta_inputs = load_delta_meta(
                dataset=dataset,
                delta_meta_path=task2_snapshot / bundle.delta_meta_relpath,
                eligible_cell_lines=eligible_cell_lines,
            )
            input_paths.extend(delta_inputs)
            fail_on_forbidden_cell_line(delta_meta, ["cell_line"], f"{dataset} delta_meta.csv")

            exploded_membership = build_membership_from_delta_meta(dataset, delta_meta)
            parity_membership, parity_inputs = load_membership_parity_table(dataset, task2_snapshot / bundle.membership_relpath)
            input_paths.extend(parity_inputs)
            validate_membership_parity(
                dataset=dataset,
                exploded_membership=exploded_membership,
                parity_membership=parity_membership,
            )

            chem_rows_by_key, gen_rows_by_key = build_base_memberships(
                dataset=dataset,
                mech_order=mech_order_dataset,
                membership=exploded_membership,
                eligible_coverage_subset=eligible_subset,
            )

            dataset_supported_grid = supported_grid.loc[supported_grid["dataset"].eq(dataset)].copy()
            dataset_representation_order = [
                representation
                for representation in REPRESENTATION_ORDER
                if representation in dataset_supported_grid["representation"].drop_duplicates().tolist()
            ]

            for representation in dataset_representation_order:
                rep_spec = physical_specs.get(dataset, {}).get(representation)
                if rep_spec is None:
                    raise ValueError(
                        f"Availability registry says available but no physical representation spec exists for dataset={dataset}, representation={representation}"
                    )

                payload, rep_inputs, rep_detail = validate_representation_inputs(
                    snapshot_root=task2_snapshot,
                    spec=rep_spec,
                    delta_meta=delta_meta,
                )
                representation_details.append(rep_detail)
                input_paths.extend(rep_inputs)

                rep_key_set = set(
                    tuple(row)
                    for row in dataset_supported_grid.loc[
                        dataset_supported_grid["representation"].eq(representation),
                        ["cell_line", "target_token"],
                    ].itertuples(index=False, name=None)
                )
                mech_keys_for_rep = [
                    (cell_line, target_token)
                    for dset, cell_line, target_token in mech_order_dataset
                    if (cell_line, target_token) in rep_key_set
                ]

                rep_tasks, rep_attrition, random_draws = prepare_group_tasks_for_representation(
                    dataset=dataset,
                    representation=representation,
                    mech_keys=mech_keys_for_rep,
                    chem_rows_by_key=chem_rows_by_key,
                    gen_rows_by_key=gen_rows_by_key,
                    valid_mask=payload.valid_mask,
                    rng=rng,
                )
                attrition_rows.extend(rep_attrition)
                total_random_draws += random_draws

                rep_results, backend = run_group_tasks(
                    tasks=rep_tasks,
                    array_path=payload.array_path,
                    expected_shape=(len(delta_meta), rep_spec.expected_dim),
                    worker_name=f"{dataset}/{representation}",
                    workers=workers,
                )
                parallel_backends_used.append(backend)
                for result in rep_results:
                    concordance_rows.append(result["row"])
                    attrition_rows.extend(result["events"])
    except Exception as exc:
        assertions.append(
            {
                "name": "representation_processing_fail_fast",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Representation processing failed: {exc}", file=sys.stderr)
        return 8

    assertions.append(
        {
            "name": "deterministic_parallelism_setup",
            "pass": True,
            "details": {
                "rules": [
                    "Process one dataset/representation at a time",
                    "Main process precomputes all EDIST subsamples from seed=619",
                    "Traversal order for deterministic subsampling is dataset order, then representation order, then canonical mech-key order, then side order Chemical -> Genetic",
                    "Worker initializer opens one read-only memmap per process",
                    "Workers never call RNG",
                ],
                "workers": workers,
                "n_rng_draws_main_process": int(total_random_draws),
                "side_order": list(SIDE_ORDER),
                "parallel_backends_used": sorted(set(parallel_backends_used)),
            },
            "counterexamples": [],
        }
    )

    concordance_df = canonicalize_concordance_frame(pd.DataFrame(concordance_rows), mech_order_all, REPRESENTATION_ORDER)
    attrition_df = canonicalize_attrition_frame(pd.DataFrame(attrition_rows), mech_order_all, REPRESENTATION_ORDER)
    assert_unique_key(
        concordance_df,
        ["dataset", "cell_line", "target_token", "representation"],
        "task2_group_concordance.csv",
    )
    assert_unique_key(
        attrition_df,
        ["dataset", "cell_line", "target_token", "representation", "metric_name", "reason"],
        "task2_group_attrition.csv",
    )

    actual_rows = int(len(concordance_df))
    grid_pass = actual_rows == expected_supported_rows
    assertions.append(
        {
            "name": "full_supported_grid_materialized",
            "pass": bool(grid_pass),
            "details": {
                "expected_s4_rows": expected_supported_rows,
                "actual_s4_rows": actual_rows,
            },
            "counterexamples": [] if grid_pass else [{"expected_s4_rows": expected_supported_rows, "actual_s4_rows": actual_rows}],
        }
    )
    if not grid_pass:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Concordance row count does not match supported analysis grid.", file=sys.stderr)
        return 9

    available_scope_keys = set(
        tuple(row)
        for row in supported_grid[["dataset", "cell_line", "target_token", "representation"]].itertuples(index=False, name=None)
    )
    output_scope_keys = set(
        tuple(row)
        for row in concordance_df[["dataset", "cell_line", "target_token", "representation"]].itertuples(index=False, name=None)
    )
    missing_scope = available_scope_keys - output_scope_keys
    extra_scope = output_scope_keys - available_scope_keys
    no_not_applicable_scope_rows = not missing_scope and not extra_scope
    assertions.append(
        {
            "name": "no_not_applicable_scope_rows_materialized",
            "pass": bool(no_not_applicable_scope_rows),
            "details": {
                "supported_grid_rows": int(len(available_scope_keys)),
                "output_rows": int(len(output_scope_keys)),
            },
            "counterexamples": (
                [{"missing_supported_key": key} for key in list(missing_scope)[:MAX_COUNTEREXAMPLES]]
                + [{"unexpected_output_key": key} for key in list(extra_scope)[:MAX_COUNTEREXAMPLES]]
            ),
        }
    )
    if not no_not_applicable_scope_rows:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Output rows drifted from supported analysis grid.", file=sys.stderr)
        return 10

    shape_contract_pass = True
    for detail in representation_details:
        spec = physical_specs[detail["dataset"]][detail["representation"]]
        if int(detail["dim"]) != int(spec.expected_dim):
            shape_contract_pass = False
            break
    assertions.append(
        {
            "name": "physical_representation_contract",
            "pass": bool(shape_contract_pass),
            "details": {"representations": representation_details},
            "counterexamples": [] if shape_contract_pass else representation_details[:MAX_COUNTEREXAMPLES],
        }
    )

    no_rank_columns = all(column not in concordance_df.columns for column in ["rank", "rank_by_mean_edist_biascorr"])
    assertions.append(
        {
            "name": "edist_row_evidence_only",
            "pass": bool(no_rank_columns),
            "details": {
                "rules": [
                    "raw edist_biascorr is emitted as row-level evidence only",
                    "S4 must not add ranking or leaderboard semantics for edist_biascorr",
                ],
            },
            "counterexamples": [] if no_rank_columns else [{"columns": concordance_df.columns.tolist()}],
        }
    )

    concordance_path = stage_dir / "task2_group_concordance.csv"
    attrition_path = stage_dir / "task2_group_attrition.csv"
    write_csv(concordance_df, concordance_path)
    write_csv(attrition_df, attrition_path)

    core_outputs = [concordance_path, attrition_path]
    output_routing_pass = all(path.resolve().is_relative_to(stage_dir.resolve()) for path in core_outputs)
    assertions.append(
        {
            "name": "output_routing_isolated",
            "pass": bool(output_routing_pass),
            "details": {"stage_dir": str(stage_dir)},
            "counterexamples": []
            if output_routing_pass
            else [{"bad_path": str(path)} for path in core_outputs if not path.resolve().is_relative_to(stage_dir.resolve())],
        }
    )

    bad_inputs = [
        str(path.absolute())
        for path in sorted(set(input_paths), key=lambda p: str(p))
        if not path.absolute().is_relative_to(task2_snapshot.absolute())
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "task2_snapshot": str(task2_snapshot.absolute()),
                "n_inputs": int(len(set(input_paths))),
            },
            "counterexamples": [{"bad_input": p} for p in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    completed_at = utc_now_iso()
    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    manifest_path = stage_dir / "manifest.json"

    output_paths = [str(path.resolve()) for path in core_outputs + [run_manifest_path, audit_assertions_path, manifest_path]]
    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s4_task2_group_concordance_multisource.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "edist_max_n": EDIST_MAX_N,
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
            "workers": workers,
            "workers_source": "cli" if args.workers is not None else "config.compute.workers",
            "dataset_order": list(DATASET_ORDER),
            "representation_order": list(REPRESENTATION_ORDER),
            "target_membership_source": "delta_meta.target_tokens",
            "expected_rows_from_supported_analysis_grid": int(expected_supported_rows),
            "rng_draws_main_process": int(total_random_draws),
            "parallel_backends_used": sorted(set(parallel_backends_used)),
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths), key=lambda p: str(p))],
        "outputs": output_paths,
        "summary": {
            "n_concordance_rows": int(len(concordance_df)),
            "n_attrition_rows": int(len(attrition_df)),
            "eligible_mech_keys_by_dataset": {
                str(key): int(value)
                for key, value in eligible_coverage.groupby("dataset", sort=False).size().to_dict().items()
            },
            "supported_grid_breakdown": supported_breakdown.to_dict(orient="records"),
            "representation_details": representation_details,
        },
    }

    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})
    write_json(manifest_path, {"stage": STAGE, "files": build_stage_manifest(stage_dir)})
    return 0


if __name__ == "__main__":
    sys.exit(main())
