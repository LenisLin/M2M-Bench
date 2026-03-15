# SCRIPT_HEADER_CONTRACT
# Script: scripts/s4_task2_group_concordance.py
# Legacy note:
#   - Preserved historical scPerturb-K562 S4 artifact.
#   - Corrected multisource Task2 S4 semantics live in docs/contracts/task2_spec.md and scripts/s4_task2_group_concordance_multisource.py.
# Purpose: Compute legacy/interim Task2 K562 group-level cross-mechanism concordance for Chemical vs Genetic cohorts.
# Inputs:
#   - Task2 Snapshot K562: config/config.yaml::paths.task2_snapshot
#     - k562/Common_Targets_K562.csv
#     - k562/derived/delta_meta.csv
#     - k562/derived/{gene_delta.npy,pathway_delta.npy}
#     - k562/fm/<representation>/{fm_delta.npy,fm_delta_meta.csv}
#     - k562/fm/<representation>/delta_operator_policy.json (optional recorded artifact only)
# Outputs:
#   - task2_group_concordance.csv: runs/<run_id>/s4_task2_group_concordance/
#   - task2_group_attrition.csv: runs/<run_id>/s4_task2_group_concordance/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s4_task2_group_concordance/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::compute.workers
# Execution:
#   - python scripts/s4_task2_group_concordance.py --run-id <run_id> --seed 619 --workers 8
# Failure Modes:
#   - Missing required Task2 inputs -> exit non-zero
#   - delta_meta/fm alignment violation -> exit non-zero
#   - Shape or finite-mask contract violation -> exit non-zero
# Last Updated: 2026-03-13

"""
Inputs:
- data/task2_snapshot_v1/k562/Common_Targets_K562.csv
- data/task2_snapshot_v1/k562/derived/delta_meta.csv
- data/task2_snapshot_v1/k562/derived/gene_delta.npy
- data/task2_snapshot_v1/k562/derived/pathway_delta.npy
- data/task2_snapshot_v1/k562/fm/<representation>/fm_delta.npy
- data/task2_snapshot_v1/k562/fm/<representation>/fm_delta_meta.csv
- data/task2_snapshot_v1/k562/fm/<representation>/delta_operator_policy.json (optional)

Outputs:
- runs/<run_id>/s4_task2_group_concordance/task2_group_concordance.csv
- runs/<run_id>/s4_task2_group_concordance/task2_group_attrition.csv
- runs/<run_id>/s4_task2_group_concordance/run_manifest.json
- runs/<run_id>/s4_task2_group_concordance/audit_assertions.json
- runs/<run_id>/s4_task2_group_concordance/manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- EDIST_MAX_N = 256
- Canonical target order = file order from Common_Targets_K562.csv
- Canonical representation order = Gene, Pathway, scgpt, geneformer, scbert, scfoundation, uce, state, tahoe-x1

Attrition rules:
- Representation valid_mask drops are recorded before metric computation.
- If either side has zero valid rows after filtering, all metrics are NA and reason=missing_one_side.
- E-distance is NA unless n_chem_sub >= 2 and n_gen_sub >= 2 after deterministic subsampling.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml

STAGE = "s4_task2_group_concordance"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
EXPECTED_CELL_LINE = "K562"
EXPECTED_DATASET = "scPerturb"

GLOBAL_SEED = 619
EDIST_MAX_N = 256
MASK_BATCH_ROWS = 256
MAX_COUNTEREXAMPLES = 5

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


@dataclass(frozen=True)
class RepresentationSpec:
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


@dataclass(frozen=True)
class GroupTask:
    representation: str
    target_token: str
    group_id: str
    n_chem_total: int
    n_gen_total: int
    chem_used_row_ids: np.ndarray
    gen_used_row_ids: np.ndarray
    chem_sub_row_ids: np.ndarray
    gen_sub_row_ids: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S4 Task2 K562 group concordance")
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


def build_representation_registry(snapshot_root: Path) -> Tuple[RepresentationSpec, ...]:
    k562_root = snapshot_root / "k562"
    return (
        RepresentationSpec(
            name="Gene",
            family="derived",
            array_relpath=str((k562_root / "derived/gene_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "derived/delta_meta.csv").relative_to(snapshot_root)),
            expected_dim=8363,
        ),
        RepresentationSpec(
            name="Pathway",
            family="derived",
            array_relpath=str((k562_root / "derived/pathway_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "derived/delta_meta.csv").relative_to(snapshot_root)),
            expected_dim=50,
        ),
        RepresentationSpec(
            name="scgpt",
            family="fm",
            array_relpath=str((k562_root / "fm/scgpt/fm_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "fm/scgpt/fm_delta_meta.csv").relative_to(snapshot_root)),
            policy_relpath=str((k562_root / "fm/scgpt/delta_operator_policy.json").relative_to(snapshot_root)),
            expected_dim=512,
        ),
        RepresentationSpec(
            name="geneformer",
            family="fm",
            array_relpath=str((k562_root / "fm/geneformer/fm_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "fm/geneformer/fm_delta_meta.csv").relative_to(snapshot_root)),
            policy_relpath=str((k562_root / "fm/geneformer/delta_operator_policy.json").relative_to(snapshot_root)),
            expected_dim=1152,
        ),
        RepresentationSpec(
            name="scbert",
            family="fm",
            array_relpath=str((k562_root / "fm/scbert/fm_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "fm/scbert/fm_delta_meta.csv").relative_to(snapshot_root)),
            policy_relpath=str((k562_root / "fm/scbert/delta_operator_policy.json").relative_to(snapshot_root)),
            expected_dim=200,
        ),
        RepresentationSpec(
            name="scfoundation",
            family="fm",
            array_relpath=str((k562_root / "fm/scfoundation/fm_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "fm/scfoundation/fm_delta_meta.csv").relative_to(snapshot_root)),
            policy_relpath=str((k562_root / "fm/scfoundation/delta_operator_policy.json").relative_to(snapshot_root)),
            expected_dim=3072,
        ),
        RepresentationSpec(
            name="uce",
            family="fm",
            array_relpath=str((k562_root / "fm/uce/fm_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "fm/uce/fm_delta_meta.csv").relative_to(snapshot_root)),
            policy_relpath=str((k562_root / "fm/uce/delta_operator_policy.json").relative_to(snapshot_root)),
            expected_dim=1280,
        ),
        RepresentationSpec(
            name="state",
            family="fm",
            array_relpath=str((k562_root / "fm/state/fm_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "fm/state/fm_delta_meta.csv").relative_to(snapshot_root)),
            policy_relpath=str((k562_root / "fm/state/delta_operator_policy.json").relative_to(snapshot_root)),
            expected_dim=2058,
        ),
        RepresentationSpec(
            name="tahoe-x1",
            family="fm",
            array_relpath=str((k562_root / "fm/tahoe-x1/fm_delta.npy").relative_to(snapshot_root)),
            meta_relpath=str((k562_root / "fm/tahoe-x1/fm_delta_meta.csv").relative_to(snapshot_root)),
            policy_relpath=str((k562_root / "fm/tahoe-x1/delta_operator_policy.json").relative_to(snapshot_root)),
            expected_dim=512,
        ),
    )


def load_common_target_order(common_targets_path: Path) -> List[str]:
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
    return target_order


def load_delta_meta(delta_meta_path: Path, target_order: Sequence[str]) -> Tuple[pd.DataFrame, List[Tuple[str, int, int]]]:
    if not delta_meta_path.is_file():
        raise FileNotFoundError(f"Missing delta_meta.csv: {delta_meta_path}")
    delta_meta = pd.read_csv(delta_meta_path)
    ensure_required_columns(
        delta_meta,
        [
            "row_id",
            "treated_cell_id",
            "perturbation_class",
            "cell_line",
            "n_controls_used",
            "dataset_side",
            "target_tokens",
        ],
        "delta_meta.csv",
    )

    delta_meta["row_id"] = pd.to_numeric(delta_meta["row_id"], errors="raise").astype(np.int64)
    delta_meta["n_controls_used"] = pd.to_numeric(delta_meta["n_controls_used"], errors="raise").astype(np.int64)
    delta_meta["treated_cell_id"] = delta_meta["treated_cell_id"].astype(str)
    delta_meta["perturbation_class"] = delta_meta["perturbation_class"].astype(str).str.strip()
    delta_meta["dataset_side"] = delta_meta["dataset_side"].astype(str).str.strip().str.upper()
    delta_meta["cell_line"] = delta_meta["cell_line"].astype(str).str.strip()
    delta_meta = delta_meta.sort_values("row_id", kind="mergesort").reset_index(drop=True)

    row_ids = delta_meta["row_id"].to_numpy(dtype=np.int64, copy=False)
    expected_row_ids = np.arange(len(delta_meta), dtype=np.int64)
    if not np.array_equal(row_ids, expected_row_ids):
        raise ValueError("delta_meta.row_id must be contiguous 0..N-1")
    if delta_meta["cell_line"].nunique() != 1 or str(delta_meta["cell_line"].iloc[0]) != EXPECTED_CELL_LINE:
        raise ValueError("delta_meta.csv must contain only cell_line=K562")

    bad_side = delta_meta.loc[~delta_meta["dataset_side"].isin(["CRISPR", "DRUG"]), "dataset_side"]
    if not bad_side.empty:
        examples = bad_side.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"delta_meta.csv contains unsupported dataset_side values: {examples}")

    parsed_tokens: List[Tuple[str, ...]] = [parse_target_tokens_field(value) for value in delta_meta["target_tokens"].tolist()]
    token_universe = {token for tokens in parsed_tokens for token in tokens}
    if token_universe != set(target_order):
        raise ValueError(
            "Atomic token set from delta_meta.target_tokens does not match Common_Targets_K562.csv"
        )

    delta_meta["parsed_target_tokens"] = parsed_tokens
    crispr_bad = delta_meta.loc[
        (delta_meta["dataset_side"] == "CRISPR")
        & (delta_meta["parsed_target_tokens"].map(len) != 1),
        ["row_id", "target_tokens"],
    ]
    if not crispr_bad.empty:
        raise ValueError(
            "CRISPR rows in delta_meta.csv must have exactly one target token, "
            f"examples={crispr_bad.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )

    coverage_before_filter: List[Tuple[str, int, int]] = []
    for target_token in target_order:
        n_chem = 0
        n_gen = 0
        for row in delta_meta.itertuples(index=False):
            tokens = row.parsed_target_tokens
            if row.dataset_side == "DRUG" and target_token in tokens:
                n_chem += 1
            elif row.dataset_side == "CRISPR" and tokens[0] == target_token:
                n_gen += 1
        coverage_before_filter.append((target_token, n_chem, n_gen))

    return delta_meta, coverage_before_filter


def build_base_memberships(
    delta_meta: pd.DataFrame,
    target_order: Sequence[str],
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    chem_rows: MutableMapping[str, List[int]] = {target: [] for target in target_order}
    gen_rows: MutableMapping[str, List[int]] = {target: [] for target in target_order}

    for row in delta_meta.itertuples(index=False):
        row_id = int(row.row_id)
        tokens = row.parsed_target_tokens
        if row.dataset_side == "DRUG":
            for token in tokens:
                chem_rows[token].append(row_id)
        else:
            gen_rows[tokens[0]].append(row_id)

    chem_arrays = {
        target: np.asarray(chem_rows[target], dtype=np.int64)
        for target in target_order
    }
    gen_arrays = {
        target: np.asarray(gen_rows[target], dtype=np.int64)
        for target in target_order
    }
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
        raise FileNotFoundError(f"Missing required array for {spec.name}: {array_path}")
    if not meta_path.is_file():
        raise FileNotFoundError(f"Missing required metadata for {spec.name}: {meta_path}")

    arr = np.load(array_path, mmap_mode="r")
    if arr.ndim != 2:
        raise ValueError(f"{spec.name} array must be 2D, got shape={tuple(arr.shape)}")
    n_rows = int(len(delta_meta))
    if int(arr.shape[0]) != n_rows:
        raise ValueError(
            f"{spec.name} row-count mismatch: array_rows={arr.shape[0]}, delta_meta_rows={n_rows}"
        )
    if int(arr.shape[1]) != int(spec.expected_dim):
        raise ValueError(
            f"{spec.name} dimension mismatch: expected={spec.expected_dim}, got={arr.shape[1]}"
        )

    meta = pd.read_csv(meta_path)
    if spec.family == "derived":
        truth_row_ids = delta_meta["row_id"].to_numpy(dtype=np.int64, copy=False)
        meta_row_ids = pd.to_numeric(meta["row_id"], errors="raise").to_numpy(dtype=np.int64, copy=False)
        if not np.array_equal(meta_row_ids, truth_row_ids):
            raise ValueError(f"{spec.name} metadata row_id mismatch against delta_meta.csv")
        valid_mask = np.ones(n_rows, dtype=bool)
        invalid_by_side: Dict[str, int] = {}
    else:
        ensure_required_columns(
            meta,
            ["row_id", "treated_cell_id", "valid_mask", "n_controls_used", "dataset_side"],
            f"{spec.name} fm_delta_meta.csv",
        )
        meta = meta.sort_values("row_id", kind="mergesort").reset_index(drop=True)
        meta["row_id"] = pd.to_numeric(meta["row_id"], errors="raise").astype(np.int64)
        meta["n_controls_used"] = pd.to_numeric(meta["n_controls_used"], errors="raise").astype(np.int64)
        meta["treated_cell_id"] = meta["treated_cell_id"].astype(str)
        meta["dataset_side"] = meta["dataset_side"].astype(str).str.strip().str.upper()

        truth_row_ids = delta_meta["row_id"].to_numpy(dtype=np.int64, copy=False)
        if not np.array_equal(meta["row_id"].to_numpy(dtype=np.int64, copy=False), truth_row_ids):
            raise ValueError(f"{spec.name} fm_delta_meta row_id mismatch against delta_meta.csv")
        if not np.array_equal(meta["treated_cell_id"].to_numpy(dtype=object), delta_meta["treated_cell_id"].to_numpy(dtype=object)):
            raise ValueError(f"{spec.name} treated_cell_id mismatch against delta_meta.csv")
        if not np.array_equal(meta["dataset_side"].to_numpy(dtype=object), delta_meta["dataset_side"].to_numpy(dtype=object)):
            raise ValueError(f"{spec.name} dataset_side mismatch against delta_meta.csv")
        if not np.array_equal(
            meta["n_controls_used"].to_numpy(dtype=np.int64, copy=False),
            delta_meta["n_controls_used"].to_numpy(dtype=np.int64, copy=False),
        ):
            raise ValueError(f"{spec.name} n_controls_used mismatch against delta_meta.csv")
        if "model_name" in meta.columns:
            bad_models = meta.loc[meta["model_name"].astype(str) != spec.name, "model_name"]
            if not bad_models.empty:
                examples = bad_models.drop_duplicates().head(MAX_COUNTEREXAMPLES).tolist()
                raise ValueError(f"{spec.name} model_name mismatch in fm_delta_meta.csv: {examples}")
        valid_mask = parse_bool_series(meta["valid_mask"])
        invalid_meta = meta.loc[~valid_mask, "dataset_side"]
        invalid_by_side = {str(key): int(value) for key, value in invalid_meta.value_counts().to_dict().items()}

    contract_errors = batched_mask_contract_check(arr, valid_mask, MASK_BATCH_ROWS)
    if contract_errors:
        raise ValueError(f"{spec.name} finite-mask contract violation: {contract_errors[0]}")

    policy_path: Optional[Path] = None
    if spec.policy_relpath is not None:
        candidate = snapshot_root / spec.policy_relpath
        if candidate.is_file():
            policy_path = candidate

    del arr

    rep_inputs = required_inputs.copy()
    if policy_path is not None:
        rep_inputs.append(policy_path)

    details = {
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
        ),
        rep_inputs,
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
        "dataset": EXPECTED_DATASET,
        "cell_line": EXPECTED_CELL_LINE,
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
    representation: str,
    target_order: Sequence[str],
    chem_rows_by_target: Mapping[str, np.ndarray],
    gen_rows_by_target: Mapping[str, np.ndarray],
    valid_mask: np.ndarray,
    rng: np.random.Generator,
) -> Tuple[List[GroupTask], List[Dict[str, object]], int]:
    tasks: List[GroupTask] = []
    attrition_rows: List[Dict[str, object]] = []
    random_draws = 0

    for target_token in target_order:
        chem_total = chem_rows_by_target[target_token]
        gen_total = gen_rows_by_target[target_token]

        chem_used = chem_total[valid_mask[chem_total]] if chem_total.size > 0 else np.empty((0,), dtype=np.int64)
        gen_used = gen_total[valid_mask[gen_total]] if gen_total.size > 0 else np.empty((0,), dtype=np.int64)

        if chem_used.size != chem_total.size or gen_used.size != gen_total.size:
            attrition_rows.append(
                make_attrition_row(
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

        chem_sub, used_chem_draw = deterministic_subsample_row_ids(chem_used, EDIST_MAX_N, rng)
        gen_sub, used_gen_draw = deterministic_subsample_row_ids(gen_used, EDIST_MAX_N, rng)
        random_draws += used_chem_draw + used_gen_draw

        tasks.append(
            GroupTask(
                representation=representation,
                target_token=target_token,
                group_id=f"{EXPECTED_DATASET}||{EXPECTED_CELL_LINE}||{target_token}",
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
        "dataset": EXPECTED_DATASET,
        "cell_line": EXPECTED_CELL_LINE,
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


def canonicalize_concordance_frame(
    frame: pd.DataFrame,
    target_order: Sequence[str],
    representation_order: Sequence[str],
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=CONCORDANCE_COLUMNS)
    rep_index = {name: idx for idx, name in enumerate(representation_order)}
    target_index = {name: idx for idx, name in enumerate(target_order)}
    out = frame.copy()
    out["_rep_idx"] = out["representation"].map(rep_index)
    out["_target_idx"] = out["target_token"].map(target_index)
    out = out.sort_values(["_rep_idx", "_target_idx"], kind="mergesort").drop(columns=["_rep_idx", "_target_idx"])
    return out.reset_index(drop=True)[CONCORDANCE_COLUMNS]


def canonicalize_attrition_frame(
    frame: pd.DataFrame,
    target_order: Sequence[str],
    representation_order: Sequence[str],
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=ATTRITION_COLUMNS)
    rep_index = {name: idx for idx, name in enumerate(representation_order)}
    target_index = {name: idx for idx, name in enumerate(target_order)}
    metric_index = {name: idx for idx, name in enumerate(METRIC_ORDER)}
    out = frame.copy()
    out["_rep_idx"] = out["representation"].map(rep_index)
    out["_target_idx"] = out["target_token"].map(target_index)
    out["_metric_idx"] = out["metric_name"].map(metric_index).fillna(len(METRIC_ORDER))
    out = out.sort_values(["_rep_idx", "_target_idx", "_metric_idx", "reason"], kind="mergesort").drop(
        columns=["_rep_idx", "_target_idx", "_metric_idx"]
    )
    return out.reset_index(drop=True)[ATTRITION_COLUMNS]


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

    task2_snapshot = resolve_config_path(project_root, str(config["paths"]["task2_snapshot"]))
    runs_dir = resolve_config_path(project_root, str(config["paths"]["runs_dir"]))

    expected_snapshot = (project_root / EXPECTED_TASK2_SNAPSHOT).resolve()
    if task2_snapshot != expected_snapshot:
        print(
            "[ERROR] Data isolation violation: config.paths.task2_snapshot must resolve "
            f"to {expected_snapshot}, got {task2_snapshot}",
            file=sys.stderr,
        )
        return 5

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
                "rules": ["GLOBAL_SEED must equal 619", "EDIST_MAX_N must equal 256"],
                "seed": GLOBAL_SEED,
                "edist_max_n": EDIST_MAX_N,
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "task2_snapshot_isolation",
            "pass": True,
            "details": {
                "rules": [
                    "S4 reads exclusively from data/task2_snapshot_v1/",
                    "config.paths.task2_snapshot must resolve to data/task2_snapshot_v1",
                ],
                "task2_snapshot": str(task2_snapshot),
            },
            "counterexamples": [],
        }
    )

    try:
        target_order = load_common_target_order(task2_snapshot / "k562/Common_Targets_K562.csv")
        delta_meta, coverage_before_filter = load_delta_meta(
            task2_snapshot / "k562/derived/delta_meta.csv",
            target_order,
        )
        chem_rows_by_target, gen_rows_by_target = build_base_memberships(delta_meta, target_order)
    except Exception as exc:
        assertions.append(
            {
                "name": "task2_group_inputs_ready",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed loading Task2 inputs: {exc}", file=sys.stderr)
        return 6

    input_paths.extend(
        [
            task2_snapshot / "k562/Common_Targets_K562.csv",
            task2_snapshot / "k562/derived/delta_meta.csv",
        ]
    )

    registry = build_representation_registry(task2_snapshot)
    registry_names = [spec.name for spec in registry]
    registry_pass = tuple(registry_names) == REPRESENTATION_ORDER
    assertions.append(
        {
            "name": "representation_registry_frozen",
            "pass": bool(registry_pass),
            "details": {
                "rules": ["Representation registry must be the reviewed frozen 9-entry list"],
                "representation_order": registry_names,
            },
            "counterexamples": [] if registry_pass else [{"representation_order": registry_names}],
        }
    )
    if not registry_pass:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Representation registry mismatch.", file=sys.stderr)
        return 7

    assertions.append(
        {
            "name": "canonical_target_order_from_common_targets_file",
            "pass": True,
            "details": {
                "rules": [
                    "Canonical target order must be file order from Common_Targets_K562.csv",
                    "No re-sorting by alphabetical order or dict/set order",
                ],
                "target_order": target_order,
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "target_membership_source_target_tokens",
            "pass": True,
            "details": {
                "rules": [
                    "S4 consumes delta_meta.csv.target_tokens as the canonical normalized atomic membership source",
                    "No runtime re-parsing from target_raw is used in this implementation",
                ],
                "column_present": "target_tokens" in delta_meta.columns,
                "n_rows": int(len(delta_meta)),
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "delta_meta_row_universe_contract",
            "pass": True,
            "details": {
                "rules": [
                    "delta_meta.row_id must equal 0..N-1",
                    "delta_meta cell_line must be K562",
                    "Atomic target token set must equal Common_Targets_K562.csv",
                ],
                "n_rows": int(len(delta_meta)),
                "cell_line": EXPECTED_CELL_LINE,
            },
            "counterexamples": [],
        }
    )

    current_snapshot_pass = all(n_chem > 0 and n_gen > 0 for _, n_chem, n_gen in coverage_before_filter)
    assertions.append(
        {
            "name": "current_snapshot_pre_filter_side_coverage_observation",
            "pass": bool(current_snapshot_pass),
            "details": {
                "rules": [
                    "This is a current-snapshot observation only; it is not a universal S4 fail-fast contract",
                ],
                "coverage": [
                    {
                        "target_token": target,
                        "n_chem_instances_total": int(n_chem),
                        "n_gen_instances_total": int(n_gen),
                    }
                    for target, n_chem, n_gen in coverage_before_filter
                ],
            },
            "counterexamples": [
                {
                    "target_token": target,
                    "n_chem_instances_total": int(n_chem),
                    "n_gen_instances_total": int(n_gen),
                }
                for target, n_chem, n_gen in coverage_before_filter
                if n_chem == 0 or n_gen == 0
            ][:MAX_COUNTEREXAMPLES],
        }
    )

    concordance_rows: List[Dict[str, object]] = []
    attrition_rows: List[Dict[str, object]] = []
    representation_details: List[Dict[str, Any]] = []
    total_random_draws = 0

    try:
        for spec in registry:
            payload, rep_inputs, rep_detail = validate_representation_inputs(
                snapshot_root=task2_snapshot,
                spec=spec,
                delta_meta=delta_meta,
            )
            representation_details.append(rep_detail)
            input_paths.extend(rep_inputs)

            rep_tasks, rep_attrition, random_draws = prepare_group_tasks_for_representation(
                representation=spec.name,
                target_order=target_order,
                chem_rows_by_target=chem_rows_by_target,
                gen_rows_by_target=gen_rows_by_target,
                valid_mask=payload.valid_mask,
                rng=rng,
            )
            attrition_rows.extend(rep_attrition)
            total_random_draws += random_draws

            max_workers = min(workers, len(rep_tasks))
            with ProcessPoolExecutor(
                max_workers=max_workers,
                initializer=init_worker,
                initargs=(payload.array_path, (len(delta_meta), spec.expected_dim), spec.name),
            ) as executor:
                futures = [executor.submit(worker_group_metric, task) for task in rep_tasks]
                for future in as_completed(futures):
                    result = future.result()
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
                    "Process one representation at a time",
                "Main process precomputes all EDIST subsamples from seed=619",
                "Traversal order for deterministic subsampling is representation registry order, then canonical target file order, then side order Chemical -> Genetic",
                "Worker initializer opens one read-only memmap per process",
                "Workers never call RNG",
            ],
            "workers": workers,
            "n_representations": int(len(registry)),
            "n_group_tasks": int(len(registry) * len(target_order)),
            "n_rng_draws_main_process": int(total_random_draws),
            "side_order": list(SIDE_ORDER),
        },
        "counterexamples": [],
    }
    )

    concordance_df = canonicalize_concordance_frame(pd.DataFrame(concordance_rows), target_order, REPRESENTATION_ORDER)
    attrition_df = canonicalize_attrition_frame(pd.DataFrame(attrition_rows), target_order, REPRESENTATION_ORDER)

    expected_rows = len(REPRESENTATION_ORDER) * len(target_order)
    grid_pass = int(len(concordance_df)) == int(expected_rows)
    assertions.append(
        {
            "name": "full_representation_target_grid_materialized",
            "pass": bool(grid_pass),
            "details": {
                "rules": [
                    "One task2_group_concordance row must exist for every representation x canonical target pair",
                ],
                "expected_rows": int(expected_rows),
                "actual_rows": int(len(concordance_df)),
            },
            "counterexamples": [] if grid_pass else [{"actual_rows": int(len(concordance_df))}],
        }
    )
    if not grid_pass:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Concordance grid row count mismatch.", file=sys.stderr)
        return 9

    rep_detail_by_name = {detail["representation"]: detail for detail in representation_details}
    shape_contract_pass = all(
        int(detail["n_rows"]) == int(len(delta_meta))
        and int(detail["dim"]) == int(next(spec.expected_dim for spec in registry if spec.name == detail["representation"]))
        for detail in representation_details
    )
    assertions.append(
        {
            "name": "representation_shape_contract",
            "pass": bool(shape_contract_pass),
            "details": {"representations": representation_details},
            "counterexamples": [] if shape_contract_pass else representation_details[:MAX_COUNTEREXAMPLES],
        }
    )

    uce_details = rep_detail_by_name.get("uce", {})
    uce_invalid_ok = (
        int(uce_details.get("n_invalid", -1)) == 384
        and dict(uce_details.get("invalid_by_side", {})) == {"DRUG": 384}
    )
    assertions.append(
        {
            "name": "uce_invalid_rows_current_snapshot_observation",
            "pass": bool(uce_invalid_ok),
            "details": {
                "rules": [
                "Current audited local snapshot is expected to preserve 384 invalid UCE rows",
                    "Current audited local snapshot is expected to place all invalid UCE rows on the DRUG side",
                    "This is a baseline observation, not a universal S4 fail-fast rule",
                ],
                "uce_n_invalid": int(uce_details.get("n_invalid", -1)),
                "uce_invalid_by_side": uce_details.get("invalid_by_side", {}),
            },
            "counterexamples": [] if uce_invalid_ok else [uce_details],
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
            "details": {
                "rules": [f"All outputs must be written under {stage_dir}"],
                "stage_dir": str(stage_dir),
            },
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
                "rules": ["All loaded S4 inputs must be under data/task2_snapshot_v1 lexical namespace"],
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
        "script_path": str((project_root / "scripts/s4_task2_group_concordance.py").resolve()),
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
            "representation_order": list(REPRESENTATION_ORDER),
            "canonical_target_order": target_order,
            "target_membership_source": "delta_meta.target_tokens",
            "fm_delta_operator_policy_json": "optional_recorded_artifact_only",
            "rng_draws_main_process": int(total_random_draws),
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths), key=lambda p: str(p))],
        "outputs": output_paths,
        "summary": {
            "n_row_universe": int(len(delta_meta)),
            "n_targets": int(len(target_order)),
            "n_representations": int(len(REPRESENTATION_ORDER)),
            "n_concordance_rows": int(len(concordance_df)),
            "n_attrition_rows": int(len(attrition_df)),
            "representation_details": representation_details,
            "current_snapshot_pre_filter_side_coverage": [
                {
                    "target_token": target,
                    "n_chem_instances_total": int(n_chem),
                    "n_gen_instances_total": int(n_gen),
                }
                for target, n_chem, n_gen in coverage_before_filter
            ],
        },
    }

    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})

    manifest_entries: List[Dict[str, object]] = []
    for file_path in sorted(stage_dir.iterdir()):
        if file_path.is_file() and file_path.name != "manifest.json":
            manifest_entries.append(
                {
                    "relative_path": file_path.name,
                    "size_bytes": int(file_path.stat().st_size),
                    "sha256": compute_sha256(file_path),
                }
            )
    write_json(manifest_path, {"stage": STAGE, "files": manifest_entries})

    return 0


if __name__ == "__main__":
    sys.exit(main())
