# SCRIPT_HEADER_CONTRACT
# Script: scripts/s2_task1_cross_metrics.py
# Purpose: Compute Task1 cross matched-pair group and retrieval metrics with deterministic parallelism.
# Inputs:
#   - Task1 Snapshot: data/task1_snapshot_v1/
#   - Cross contract: data/task1_snapshot_v1/cross_contract/cross-pairs-genetic-contract
# Outputs:
#   - task1_cross_retrieval_per_query.parquet: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_retrieval_summary.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_chance_identity_check.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_leaderboard_long.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_alignment_proof.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_attrition.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_group_cross.parquet: runs/<run_id>/s2_task1_cross_metrics/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s2_task1_cross_metrics/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task1_snapshot
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::compute.workers
# Execution:
#   - python scripts/s2_task1_cross_metrics.py --run-id <run_id> --seed 619 --workers 20
# Failure Modes:
#   - Cross-contract metadata drift -> exit non-zero
#   - Chance identity check tolerance > 1e-12 -> exit non-zero
# Last Updated: 2026-03-05

"""
Inputs:
- task1 snapshot files under data/task1_snapshot_v1/ only:
  - lincs/lincs-engine1-meta.csv
  - lincs/lincs-engine1-gene-delta.pt
  - pathway/hallmark-w-2477x50.npy
  - pathway/lincs-pathway-policy.json
  - scperturb_delta/scperturb-crispr-delta-meta.csv
  - scperturb_delta/scperturb-crispr-gene-delta.npy
  - scperturb_delta/scperturb-crispr-pathway-delta.npy
  - cross_contract/cross-pairs-genetic-contract

Outputs:
- task1_cross_retrieval_per_query.parquet
- task1_cross_retrieval_summary.csv
- task1_cross_chance_identity_check.csv
- task1_cross_leaderboard_long.csv
- task1_cross_alignment_proof.csv
- task1_cross_attrition.csv
- task1_group_cross.parquet
- run_manifest.json
- audit_assertions.json
- manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- EDIST_MAX_N = 256
- HIT_K_LIST = (1, 5, 10)
- CHANCE_IDENTITY_TOL = 1e-12
- cross_alignment_contract = "global_idx_lincs + sc_delta_row_idx"
- LOO policy (cross) = disjoint_no_leakage

Deterministic parallelism:
- All randomized subsample indices are generated once on main process from GLOBAL_SEED.
- Workers never call RNG.
- Task payloads are lightweight (indices / ids); workers use global read-only subset payload.
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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import torch
import yaml

STAGE = "s2_task1_cross_metrics"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK1_SNAPSHOT = Path("data/task1_snapshot_v1")

GLOBAL_SEED = 619
EDIST_MAX_N = 256
HIT_K_LIST: Tuple[int, ...] = (1, 5, 10)
CHANCE_IDENTITY_TOL = 1e-12
MAX_COUNTEREXAMPLES = 5
QUERY_CHUNK_SIZE = 128
CROSS_ALIGNMENT_CONTRACT = "global_idx_lincs + sc_delta_row_idx"
CROSS_DIRECTIONS: Tuple[str, str] = ("LINCS_to_scPerturb", "scPerturb_to_LINCS")
GROUP_DATASET_OR_DIRECTION = "LINCS_to_scPerturb"

SUMMARY_COLUMNS = [
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

CHANCE_COLUMNS = [
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

ATTRITION_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "reason",
    "n_dropped",
    "n_total_before",
    "notes",
]

LEADERBOARD_COLUMNS = [
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

GROUP_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "group_id",
    "cell_line",
    "target_token",
    "n_L",
    "n_S",
    "n_L_sub",
    "n_S_sub",
    "cosine",
    "pcc",
    "edist",
    "cross_alignment_contract",
]

# Process worker global payload
_WORKER_PAYLOAD: Dict[str, Dict[str, Any]] = {}


class VectorStore:
    def load_rows(self, row_indices: np.ndarray) -> np.ndarray:
        raise NotImplementedError


class TorchMatrixStore(VectorStore):
    def __init__(self, tensor: torch.Tensor) -> None:
        self.tensor = tensor

    def load_rows(self, row_indices: np.ndarray) -> np.ndarray:
        if row_indices.size == 0:
            width = int(self.tensor.shape[1])
            return np.empty((0, width), dtype=np.float64)
        idx_t = torch.as_tensor(row_indices.astype(np.int64, copy=False), dtype=torch.long)
        mat = self.tensor.index_select(0, idx_t)
        return mat.detach().cpu().numpy().astype(np.float64, copy=False)


class NpyMatrixStore(VectorStore):
    def __init__(self, array_path: Path) -> None:
        self.arr = np.load(array_path, mmap_mode="r")

    def load_rows(self, row_indices: np.ndarray) -> np.ndarray:
        if row_indices.size == 0:
            width = int(self.arr.shape[1])
            return np.empty((0, width), dtype=np.float64)
        mat = self.arr[row_indices.astype(np.int64, copy=False)]
        return np.asarray(mat, dtype=np.float64)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S2 Task1 cross metrics")
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


def normalize_text(series: pd.Series) -> pd.Series:
    return series.fillna("NA").astype(str).str.strip()


def normalize_scalar_text(value: object) -> str:
    return str(value).strip()


def normalize_perturbation_type(raw_value: object) -> str:
    value = str(raw_value).strip().lower()
    if value == "drug":
        return "Chemical"
    return "Genetic"


def ensure_required_columns(frame: pd.DataFrame, required: Sequence[str], name: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns for {name}: {missing}")


def init_global_rng(seed: int) -> np.random.Generator:
    random.seed(seed)
    np.random.seed(seed)
    return np.random.default_rng(seed)


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(a, b) / denom)


def pearson_corr(a: np.ndarray, b: np.ndarray) -> float:
    a_center = a - np.mean(a)
    b_center = b - np.mean(b)
    denom = float(np.linalg.norm(a_center) * np.linalg.norm(b_center))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(a_center, b_center) / denom)


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


def expected_mrr_for_single_positive(n_gallery: int) -> float:
    ranks = np.arange(1, n_gallery + 1, dtype=np.float64)
    return float(np.mean(1.0 / ranks))


def expected_hitk_for_single_positive(n_gallery: int, k: int) -> float:
    return float(min(k, n_gallery) / n_gallery)


def hash_gallery_group_ids(group_ids: Sequence[str]) -> str:
    material = "||".join(group_ids)
    return hashlib.sha256(material.encode("utf-8")).hexdigest()


def make_group_id(perturbation_type: str, cell_line: str, target: str) -> str:
    return f"{perturbation_type}||{cell_line}||{target}"


def record_attrition(
    bucket: MutableMapping[Tuple[str, str, str, str, str], Dict[str, object]],
    scope: str,
    dataset_or_direction: str,
    perturbation_type: str,
    representation: str,
    reason: str,
    n_dropped: int,
    n_total_before: int,
    notes: str,
    allow_zero: bool = False,
) -> None:
    if n_dropped < 0:
        return
    if n_dropped == 0 and not allow_zero:
        return

    agg_key = (
        scope,
        dataset_or_direction,
        perturbation_type,
        representation,
        reason,
    )
    if agg_key not in bucket:
        bucket[agg_key] = {
            "scope": scope,
            "dataset_or_direction": dataset_or_direction,
            "perturbation_type": perturbation_type,
            "representation": representation,
            "reason": reason,
            "n_dropped": 0,
            "n_total_before": 0,
            "notes": notes,
        }

    bucket[agg_key]["n_dropped"] = int(bucket[agg_key]["n_dropped"]) + int(n_dropped)
    bucket[agg_key]["n_total_before"] = int(bucket[agg_key]["n_total_before"]) + int(n_total_before)


def per_query_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("scope", pa.string()),
            pa.field("dataset_or_direction", pa.string()),
            pa.field("perturbation_type", pa.string()),
            pa.field("representation", pa.string()),
            pa.field("query_uid", pa.string()),
            pa.field("cell_line", pa.string()),
            pa.field("target_token", pa.string()),
            pa.field("N_gallery", pa.int64()),
            pa.field("N_gallery_max", pa.int64()),
            pa.field("gallery_group_ids_hash", pa.string()),
            pa.field("m_pos", pa.int64()),
            pa.field("rank_true", pa.int64()),
            pa.field("mrr_raw", pa.float64()),
            pa.field("hit1_raw", pa.float64()),
            pa.field("hit5_raw", pa.float64()),
            pa.field("hit10_raw", pa.float64()),
            pa.field("expected_mrr_chance", pa.float64()),
            pa.field("expected_hit1_chance", pa.float64()),
            pa.field("expected_hit5_chance", pa.float64()),
            pa.field("expected_hit10_chance", pa.float64()),
            pa.field("mrr_corrected", pa.float64()),
            pa.field("hit1_corrected", pa.float64()),
            pa.field("hit5_corrected", pa.float64()),
            pa.field("hit10_corrected", pa.float64()),
            pa.field("loo_policy", pa.string()),
            pa.field("cross_alignment_contract", pa.string()),
        ]
    )


def group_cross_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("scope", pa.string()),
            pa.field("dataset_or_direction", pa.string()),
            pa.field("perturbation_type", pa.string()),
            pa.field("representation", pa.string()),
            pa.field("group_id", pa.string()),
            pa.field("cell_line", pa.string()),
            pa.field("target_token", pa.string()),
            pa.field("n_L", pa.int64()),
            pa.field("n_S", pa.int64()),
            pa.field("n_L_sub", pa.int64()),
            pa.field("n_S_sub", pa.int64()),
            pa.field("cosine", pa.float64()),
            pa.field("pcc", pa.float64()),
            pa.field("edist", pa.float64()),
            pa.field("cross_alignment_contract", pa.string()),
        ]
    )


def canonicalize_per_query_frame(frame: pd.DataFrame) -> pd.DataFrame:
    int_cols = ["N_gallery", "N_gallery_max", "m_pos", "rank_true"]
    float_cols = [
        "mrr_raw",
        "hit1_raw",
        "hit5_raw",
        "hit10_raw",
        "expected_mrr_chance",
        "expected_hit1_chance",
        "expected_hit5_chance",
        "expected_hit10_chance",
        "mrr_corrected",
        "hit1_corrected",
        "hit5_corrected",
        "hit10_corrected",
    ]
    str_cols = [
        "scope",
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "query_uid",
        "cell_line",
        "target_token",
        "gallery_group_ids_hash",
        "loo_policy",
        "cross_alignment_contract",
    ]

    out = frame.copy()
    for col in int_cols:
        out[col] = out[col].astype(np.int64)
    for col in float_cols:
        out[col] = out[col].astype(np.float64)
    for col in str_cols:
        out[col] = out[col].fillna("NA").astype(str)
    return out


def canonicalize_group_frame(frame: pd.DataFrame) -> pd.DataFrame:
    int_cols = ["n_L", "n_S", "n_L_sub", "n_S_sub"]
    float_cols = ["cosine", "pcc", "edist"]
    str_cols = [
        "scope",
        "dataset_or_direction",
        "perturbation_type",
        "representation",
        "group_id",
        "cell_line",
        "target_token",
        "cross_alignment_contract",
    ]

    out = frame.copy()
    for col in int_cols:
        out[col] = out[col].astype(np.int64)
    for col in float_cols:
        out[col] = out[col].astype(np.float64)
    for col in str_cols:
        out[col] = out[col].fillna("NA").astype(str)
    return out


def write_parquet_with_schema(
    frame: pd.DataFrame,
    path: Path,
    schema: pa.Schema,
    canonicalizer,
) -> None:
    if frame.empty:
        empty_data = {field.name: pa.array([], type=field.type) for field in schema}
        empty_table = pa.Table.from_pydict(empty_data, schema=schema)
        pq.write_table(empty_table, path)
        return

    out = canonicalizer(frame)
    table = pa.Table.from_pandas(out, schema=schema, preserve_index=False)
    pq.write_table(table, path)


def locate_cross_contract_file(snapshot_root: Path) -> Path:
    cross_root = snapshot_root / "cross_contract"
    candidates = [
        cross_root / "cross-pairs-genetic-contract",
        cross_root / "cross_pairs_genetic_contract",
        cross_root / "cross_pairs_genetic_contract.csv",
        cross_root / "cross-pairs-genetic-contract.csv",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Missing cross_pairs_genetic_contract file under data/task1_snapshot_v1/cross_contract/"
    )


def load_lincs_inputs(
    snapshot_root: Path,
) -> Tuple[pd.DataFrame, VectorStore, np.ndarray, Dict[str, object], List[Path]]:
    lincs_meta_path = snapshot_root / "lincs/lincs-engine1-meta.csv"
    lincs_gene_pt_path = snapshot_root / "lincs/lincs-engine1-gene-delta.pt"
    pathway_w_path = snapshot_root / "pathway/hallmark-w-2477x50.npy"
    pathway_policy_path = snapshot_root / "pathway/lincs-pathway-policy.json"

    lincs_meta = pd.read_csv(lincs_meta_path)
    ensure_required_columns(
        lincs_meta,
        ["pert_type", "target", "cell_line"],
        "lincs-engine1-meta.csv",
    )

    lincs_payload = torch.load(lincs_gene_pt_path, map_location="cpu")
    if not isinstance(lincs_payload, dict) or "y_delta_gene" not in lincs_payload:
        raise ValueError("lincs-engine1-gene-delta.pt must be dict with y_delta_gene")

    y_delta_gene = lincs_payload["y_delta_gene"]
    if not torch.is_tensor(y_delta_gene):
        raise ValueError("y_delta_gene must be torch.Tensor")
    if int(y_delta_gene.shape[0]) != int(len(lincs_meta)):
        raise ValueError(
            f"LINCS row mismatch: tensor={int(y_delta_gene.shape[0])}, meta={len(lincs_meta)}"
        )

    projection_w = np.load(pathway_w_path)
    if projection_w.shape != (2477, 50):
        raise ValueError(f"hallmark-w-2477x50.npy must have shape (2477,50), got {projection_w.shape}")

    with pathway_policy_path.open("r", encoding="utf-8") as handle:
        pathway_policy = json.load(handle)

    mode = str(pathway_policy.get("mode", "")).strip()
    if mode != "project_on_load":
        raise ValueError(f"LINCS pathway policy mode must be project_on_load, got {mode!r}")

    actual_w_sha256 = compute_sha256(pathway_w_path)
    policy_w_sha256 = pathway_policy.get("W_sha256")
    if policy_w_sha256 and str(policy_w_sha256).strip() != actual_w_sha256:
        raise ValueError(
            "hallmark W sha256 mismatch with pathway policy: "
            f"policy={policy_w_sha256}, actual={actual_w_sha256}"
        )

    input_paths = [lincs_meta_path, lincs_gene_pt_path, pathway_w_path, pathway_policy_path]
    return (
        lincs_meta,
        TorchMatrixStore(y_delta_gene),
        np.asarray(projection_w, dtype=np.float64),
        pathway_policy,
        input_paths,
    )


def load_scperturb_crispr_inputs(
    snapshot_root: Path,
) -> Tuple[pd.DataFrame, VectorStore, VectorStore, List[Path]]:
    crispr_meta_path = snapshot_root / "scperturb_delta/scperturb-crispr-delta-meta.csv"
    crispr_gene_path = snapshot_root / "scperturb_delta/scperturb-crispr-gene-delta.npy"
    crispr_pathway_path = snapshot_root / "scperturb_delta/scperturb-crispr-pathway-delta.npy"

    crispr_meta = pd.read_csv(crispr_meta_path)
    ensure_required_columns(
        crispr_meta,
        ["delta_row_idx", "target_std", "cell_std", "pert_type"],
        "scperturb-crispr-delta-meta.csv",
    )
    if crispr_meta["delta_row_idx"].duplicated().any():
        dup = (
            crispr_meta.loc[crispr_meta["delta_row_idx"].duplicated(), "delta_row_idx"]
            .head(MAX_COUNTEREXAMPLES)
            .tolist()
        )
        raise ValueError(f"CRISPR delta_row_idx contains duplicates, examples={dup}")

    crispr_gene = NpyMatrixStore(crispr_gene_path)
    crispr_pathway = NpyMatrixStore(crispr_pathway_path)

    idx = crispr_meta["delta_row_idx"].to_numpy(dtype=np.int64)
    if idx.size > 0:
        if np.min(idx) < 0 or np.max(idx) >= int(crispr_gene.arr.shape[0]):
            raise ValueError("CRISPR delta_row_idx out of gene array bounds")
        if np.min(idx) < 0 or np.max(idx) >= int(crispr_pathway.arr.shape[0]):
            raise ValueError("CRISPR delta_row_idx out of pathway array bounds")

    input_paths = [crispr_meta_path, crispr_gene_path, crispr_pathway_path]
    return crispr_meta, crispr_gene, crispr_pathway, input_paths


def validate_cross_contract(
    contract_df: pd.DataFrame,
    lincs_meta: pd.DataFrame,
    sc_crispr_meta: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    required_columns = [
        "expected_cell_line",
        "expected_target",
        "global_idx_lincs",
        "sc_delta_row_idx",
    ]
    ensure_required_columns(contract_df, required_columns, "cross_pairs_genetic_contract")

    contract = contract_df.copy()
    contract["global_idx_lincs"] = pd.to_numeric(contract["global_idx_lincs"], errors="raise").astype(np.int64)
    contract["sc_delta_row_idx"] = pd.to_numeric(contract["sc_delta_row_idx"], errors="raise").astype(np.int64)

    l_idx = contract["global_idx_lincs"].to_numpy(dtype=np.int64)
    s_idx = contract["sc_delta_row_idx"].to_numpy(dtype=np.int64)

    if l_idx.size > 0:
        if int(np.min(l_idx)) < 0 or int(np.max(l_idx)) >= int(len(lincs_meta)):
            raise ValueError("cross contract global_idx_lincs is out of LINCS metadata bounds")

    sc_lookup = sc_crispr_meta.set_index("delta_row_idx", drop=False)
    missing_sc = [int(v) for v in s_idx.tolist() if int(v) not in sc_lookup.index]
    if missing_sc:
        raise ValueError(
            "cross contract sc_delta_row_idx not found in CRISPR metadata, "
            f"examples={missing_sc[:MAX_COUNTEREXAMPLES]}"
        )

    l_sub = lincs_meta.iloc[l_idx][["cell_line", "target", "pert_type"]].reset_index(drop=True)
    s_sub = sc_lookup.loc[s_idx, ["cell_std", "target_std", "pert_type"]].reset_index(drop=True)

    expected_cell = normalize_text(contract["expected_cell_line"]).to_numpy()
    expected_target = normalize_text(contract["expected_target"]).to_numpy()
    l_cell = normalize_text(l_sub["cell_line"]).to_numpy()
    l_target = normalize_text(l_sub["target"]).to_numpy()
    s_cell = normalize_text(s_sub["cell_std"]).to_numpy()
    s_target = normalize_text(s_sub["target_std"]).to_numpy()

    l_cell_mismatch = l_cell != expected_cell
    l_target_mismatch = l_target != expected_target
    s_cell_mismatch = s_cell != expected_cell
    s_target_mismatch = s_target != expected_target

    if l_cell_mismatch.any() or l_target_mismatch.any() or s_cell_mismatch.any() or s_target_mismatch.any():
        mismatch_idx = np.where(
            l_cell_mismatch | l_target_mismatch | s_cell_mismatch | s_target_mismatch
        )[0][:MAX_COUNTEREXAMPLES]

        examples: List[Dict[str, object]] = []
        for i in mismatch_idx.tolist():
            examples.append(
                {
                    "row": int(i),
                    "expected_cell_line": expected_cell[i],
                    "expected_target": expected_target[i],
                    "lincs_cell_line": l_cell[i],
                    "lincs_target": l_target[i],
                    "sc_cell_line": s_cell[i],
                    "sc_target": s_target[i],
                }
            )
        raise ValueError(
            "Cross-contract validation failed (Section 1.2 metadata drift). "
            f"counterexamples={examples}"
        )

    l_pert = np.asarray([normalize_perturbation_type(v) for v in l_sub["pert_type"]], dtype=object)
    s_pert = np.asarray([normalize_perturbation_type(v) for v in s_sub["pert_type"]], dtype=object)
    pert_mismatch = l_pert != s_pert
    if pert_mismatch.any():
        idx_bad = np.where(pert_mismatch)[0][:MAX_COUNTEREXAMPLES]
        examples = [
            {"row": int(i), "lincs_perturbation_type": l_pert[i], "sc_perturbation_type": s_pert[i]}
            for i in idx_bad.tolist()
        ]
        raise ValueError(
            "Cross-contract validation failed (perturbation_type mismatch between sides). "
            f"counterexamples={examples}"
        )

    validated_pairs = pd.DataFrame(
        {
            "perturbation_type": l_pert,
            "cell_line": expected_cell,
            "target": expected_target,
            "global_idx_lincs": l_idx.astype(np.int64),
            "sc_delta_row_idx": s_idx.astype(np.int64),
        }
    )

    details = {
        "n_rows": int(len(validated_pairs)),
        "n_unique_lincs_indices": int(validated_pairs["global_idx_lincs"].nunique()),
        "n_unique_sc_indices": int(validated_pairs["sc_delta_row_idx"].nunique()),
        "n_unique_perturbation_types": int(validated_pairs["perturbation_type"].nunique()),
    }
    return validated_pairs, details


def build_worker_payload(
    validated_pairs: pd.DataFrame,
    lincs_store: VectorStore,
    sc_gene_store: VectorStore,
    sc_pathway_store: VectorStore,
    pathway_w: np.ndarray,
    eligible_perturbation_types: Sequence[str],
) -> Dict[str, Dict[str, Any]]:
    payload: Dict[str, Dict[str, Any]] = {}
    for perturbation_type in sorted(eligible_perturbation_types):
        pair_sub = (
            validated_pairs.loc[validated_pairs["perturbation_type"] == perturbation_type]
            .copy()
            .sort_values(
                ["cell_line", "target", "global_idx_lincs", "sc_delta_row_idx"],
                kind="mergesort",
            )
            .reset_index(drop=True)
        )
        if pair_sub.empty:
            continue

        l_rows = pair_sub["global_idx_lincs"].to_numpy(dtype=np.int64)
        s_rows = pair_sub["sc_delta_row_idx"].to_numpy(dtype=np.int64)

        l_gene = lincs_store.load_rows(l_rows)
        l_pathway = np.asarray(l_gene @ pathway_w, dtype=np.float64)
        s_gene = sc_gene_store.load_rows(s_rows)
        s_pathway = sc_pathway_store.load_rows(s_rows)

        row_cell_line = normalize_text(pair_sub["cell_line"]).to_numpy()
        row_target = normalize_text(pair_sub["target"]).to_numpy()
        row_group_ids = np.asarray(
            [
                make_group_id(perturbation_type, str(c), str(t))
                for c, t in zip(row_cell_line.tolist(), row_target.tolist())
            ],
            dtype=object,
        )

        l_uids = l_rows.astype(str)
        s_uids = s_rows.astype(str)

        group_to_all_positions: Dict[str, List[int]] = {}
        for i, group_id in enumerate(row_group_ids.tolist()):
            group_to_all_positions.setdefault(str(group_id), []).append(i)

        group_ids = sorted(group_to_all_positions.keys())

        reps = {
            "Gene": (l_gene, s_gene),
            "Pathway": (l_pathway, s_pathway),
        }

        for representation, (l_vectors_raw, s_vectors_raw) in reps.items():
            l_vectors = np.asarray(l_vectors_raw, dtype=np.float64)
            s_vectors = np.asarray(s_vectors_raw, dtype=np.float64)

            l_finite = np.isfinite(l_vectors).all(axis=1)
            s_finite = np.isfinite(s_vectors).all(axis=1)

            group_to_l_pos: Dict[str, np.ndarray] = {}
            group_to_s_pos: Dict[str, np.ndarray] = {}
            group_to_cell_target: Dict[str, Tuple[str, str]] = {}

            for group_id in group_ids:
                pos_all = np.asarray(group_to_all_positions[group_id], dtype=np.int64)
                group_to_l_pos[group_id] = pos_all[l_finite[pos_all]]
                group_to_s_pos[group_id] = pos_all[s_finite[pos_all]]
                first = int(pos_all[0])
                group_to_cell_target[group_id] = (
                    str(row_cell_line[first]),
                    str(row_target[first]),
                )

            l_gallery_group_ids = sorted([g for g in group_ids if group_to_l_pos[g].size > 0])
            s_gallery_group_ids = sorted([g for g in group_ids if group_to_s_pos[g].size > 0])

            dim = int(l_vectors.shape[1])
            if l_gallery_group_ids:
                l_centroids = np.stack(
                    [np.mean(l_vectors[group_to_l_pos[g]], axis=0) for g in l_gallery_group_ids],
                    axis=0,
                ).astype(np.float64, copy=False)
            else:
                l_centroids = np.empty((0, dim), dtype=np.float64)

            if s_gallery_group_ids:
                s_centroids = np.stack(
                    [np.mean(s_vectors[group_to_s_pos[g]], axis=0) for g in s_gallery_group_ids],
                    axis=0,
                ).astype(np.float64, copy=False)
            else:
                s_centroids = np.empty((0, dim), dtype=np.float64)

            dataset_id = f"{perturbation_type}|{representation}"
            payload[dataset_id] = {
                "scope": "cross",
                "dataset_id": dataset_id,
                "perturbation_type": perturbation_type,
                "representation": representation,
                "lincs_vectors": l_vectors,
                "sc_vectors": s_vectors,
                "lincs_uids": l_uids,
                "sc_uids": s_uids,
                "row_cell_line": row_cell_line.astype(object),
                "row_target": row_target.astype(object),
                "row_group_ids": row_group_ids,
                "group_ids": group_ids,
                "group_to_lincs_pos": group_to_l_pos,
                "group_to_sc_pos": group_to_s_pos,
                "group_to_cell_target": group_to_cell_target,
                "lincs_gallery_group_ids": l_gallery_group_ids,
                "sc_gallery_group_ids": s_gallery_group_ids,
                "lincs_centroids": l_centroids,
                "sc_centroids": s_centroids,
                "lincs_centroid_norms": np.sum(l_centroids * l_centroids, axis=1, dtype=np.float64),
                "sc_centroid_norms": np.sum(s_centroids * s_centroids, axis=1, dtype=np.float64),
                "group_to_lincs_gallery_idx": {
                    g: i for i, g in enumerate(l_gallery_group_ids)
                },
                "group_to_sc_gallery_idx": {
                    g: i for i, g in enumerate(s_gallery_group_ids)
                },
                "gallery_hash_l2s": hash_gallery_group_ids(s_gallery_group_ids),
                "gallery_hash_s2l": hash_gallery_group_ids(l_gallery_group_ids),
                "n_pairs": int(len(pair_sub)),
            }

    return payload


def deterministic_subsample_positions(
    positions: np.ndarray,
    uid_array: np.ndarray,
    max_n: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, int]:
    if positions.size <= max_n:
        return positions.copy(), 0

    ordered = np.argsort(uid_array[positions].astype(str), kind="mergesort")
    sorted_pos = positions[ordered]
    perm = rng.permutation(sorted_pos.size)
    return sorted_pos[perm[:max_n]], 1


def precompute_group_tasks(
    payload: Mapping[str, Mapping[str, Any]],
    rng: np.random.Generator,
) -> Tuple[List[Tuple[str, str, np.ndarray, np.ndarray]], int]:
    tasks: List[Tuple[str, str, np.ndarray, np.ndarray]] = []
    random_draws = 0

    for dataset_id in sorted(payload.keys()):
        rep = payload[dataset_id]
        group_ids = sorted(rep["group_ids"])
        l_uids = rep["lincs_uids"]
        s_uids = rep["sc_uids"]

        for group_id in group_ids:
            l_pos = rep["group_to_lincs_pos"][group_id]
            s_pos = rep["group_to_sc_pos"][group_id]

            l_sub, used_l = deterministic_subsample_positions(l_pos, l_uids, EDIST_MAX_N, rng)
            s_sub, used_s = deterministic_subsample_positions(s_pos, s_uids, EDIST_MAX_N, rng)
            random_draws += used_l + used_s

            tasks.append(
                (
                    str(dataset_id),
                    str(group_id),
                    l_sub.astype(np.int64, copy=False),
                    s_sub.astype(np.int64, copy=False),
                )
            )
    return tasks, int(random_draws)


def precompute_retrieval_tasks(
    payload: Mapping[str, Mapping[str, Any]],
    chunk_size: int,
) -> List[Tuple[str, str, np.ndarray]]:
    tasks: List[Tuple[str, str, np.ndarray]] = []
    chunk_size = max(1, int(chunk_size))

    for dataset_id in sorted(payload.keys()):
        rep = payload[dataset_id]

        for direction in CROSS_DIRECTIONS:
            query_uids = rep["lincs_uids"] if direction == "LINCS_to_scPerturb" else rep["sc_uids"]
            order = np.argsort(query_uids.astype(str), kind="mergesort")
            for start in range(0, order.size, chunk_size):
                chunk = order[start : start + chunk_size].astype(np.int64, copy=False)
                tasks.append((str(dataset_id), direction, chunk))
    return tasks


def init_worker(payload: Dict[str, Dict[str, Any]]) -> None:
    global _WORKER_PAYLOAD
    _WORKER_PAYLOAD = payload


def worker_group_metric(task: Tuple[str, str, np.ndarray, np.ndarray]) -> Dict[str, Any]:
    dataset_id, group_id, l_sub_pos, s_sub_pos = task
    rep = _WORKER_PAYLOAD[dataset_id]

    l_pos = rep["group_to_lincs_pos"][group_id]
    s_pos = rep["group_to_sc_pos"][group_id]
    l_vectors = rep["lincs_vectors"]
    s_vectors = rep["sc_vectors"]
    cell_line, target_token = rep["group_to_cell_target"][group_id]

    n_l = int(l_pos.size)
    n_s = int(s_pos.size)
    n_l_sub = int(l_sub_pos.size)
    n_s_sub = int(s_sub_pos.size)

    row = {
        "scope": "cross",
        "dataset_or_direction": GROUP_DATASET_OR_DIRECTION,
        "perturbation_type": rep["perturbation_type"],
        "representation": rep["representation"],
        "group_id": group_id,
        "cell_line": cell_line,
        "target_token": target_token,
        "n_L": n_l,
        "n_S": n_s,
        "n_L_sub": n_l_sub,
        "n_S_sub": n_s_sub,
        "cosine": np.nan,
        "pcc": np.nan,
        "edist": np.nan,
        "cross_alignment_contract": CROSS_ALIGNMENT_CONTRACT,
    }

    events: List[Dict[str, Any]] = []
    if n_l == 0 or n_s == 0:
        events.append(
            {
                "reason": "missing_one_side",
                "n_dropped": 1,
                "n_total_before": 1,
                "notes": "Cross group requires finite vectors on both sides",
            }
        )
        return {"row": row, "events": events}

    centroid_l = np.mean(l_vectors[l_pos], axis=0)
    centroid_s = np.mean(s_vectors[s_pos], axis=0)
    row["cosine"] = float(cosine_similarity(centroid_l, centroid_s))
    row["pcc"] = float(pearson_corr(centroid_l, centroid_s))

    if n_l_sub < 2 or n_s_sub < 2:
        events.append(
            {
                "reason": "edist_insufficient_cells",
                "n_dropped": 1,
                "n_total_before": 1,
                "notes": "Energy distance requires >=2 per side after deterministic subsample",
            }
        )
        return {"row": row, "events": events}

    edist = energy_distance_biascorr(l_vectors[l_sub_pos], s_vectors[s_sub_pos])
    row["edist"] = float(edist) if np.isfinite(edist) else np.nan
    return {"row": row, "events": events}


def worker_retrieval_chunk(task: Tuple[str, str, np.ndarray]) -> Dict[str, Any]:
    dataset_id, direction, query_positions = task
    rep = _WORKER_PAYLOAD[dataset_id]

    if direction == "LINCS_to_scPerturb":
        query_vectors = rep["lincs_vectors"]
        query_uids = rep["lincs_uids"]
        gallery_centroids = rep["sc_centroids"]
        gallery_norms = rep["sc_centroid_norms"]
        group_to_true = rep["group_to_sc_gallery_idx"]
        gallery_hash = rep["gallery_hash_l2s"]
    else:
        query_vectors = rep["sc_vectors"]
        query_uids = rep["sc_uids"]
        gallery_centroids = rep["lincs_centroids"]
        gallery_norms = rep["lincs_centroid_norms"]
        group_to_true = rep["group_to_lincs_gallery_idx"]
        gallery_hash = rep["gallery_hash_s2l"]

    row_cells = rep["row_cell_line"]
    row_targets = rep["row_target"]
    row_groups = rep["row_group_ids"]

    n_gallery = int(gallery_centroids.shape[0])
    expected_mrr = expected_mrr_for_single_positive(n_gallery) if n_gallery > 0 else np.nan
    expected_hit = {
        k: expected_hitk_for_single_positive(n_gallery, k) if n_gallery > 0 else np.nan
        for k in HIT_K_LIST
    }

    n_total = int(query_positions.size)
    n_valid = 0
    n_missing_metric = 0
    n_mpos0 = 0

    sum_mrr_raw = 0.0
    sum_exp_mrr = 0.0
    sum_mrr_corr = 0.0

    sum_hit_raw: Dict[int, float] = {k: 0.0 for k in HIT_K_LIST}
    sum_hit_exp: Dict[int, float] = {k: 0.0 for k in HIT_K_LIST}
    sum_hit_corr: Dict[int, float] = {k: 0.0 for k in HIT_K_LIST}

    rows: List[Dict[str, Any]] = []

    for pos in query_positions.tolist():
        q = query_vectors[int(pos)]
        if not np.isfinite(q).all():
            n_missing_metric += 1
            continue

        group_id = str(row_groups[int(pos)])
        true_idx = group_to_true.get(group_id, None)
        if true_idx is None or n_gallery <= 0:
            n_mpos0 += 1
            continue

        q_norm = float(np.dot(q, q))
        scores = -(q_norm + gallery_norms - 2.0 * (gallery_centroids @ q))
        score_true = float(scores[int(true_idx)])
        rank_true = int(1 + np.sum(scores > score_true))

        mrr_raw = 1.0 / float(rank_true)
        hit_raw = {k: 1.0 if rank_true <= k else 0.0 for k in HIT_K_LIST}
        mrr_corr = float(mrr_raw - expected_mrr)
        hit_corr = {k: float(hit_raw[k] - expected_hit[k]) for k in HIT_K_LIST}

        n_valid += 1
        sum_mrr_raw += mrr_raw
        sum_exp_mrr += float(expected_mrr)
        sum_mrr_corr += mrr_corr
        for k in HIT_K_LIST:
            sum_hit_raw[k] += hit_raw[k]
            sum_hit_exp[k] += float(expected_hit[k])
            sum_hit_corr[k] += hit_corr[k]

        rows.append(
            {
                "scope": "cross",
                "dataset_or_direction": direction,
                "perturbation_type": rep["perturbation_type"],
                "representation": rep["representation"],
                "query_uid": str(query_uids[int(pos)]),
                "cell_line": str(row_cells[int(pos)]),
                "target_token": str(row_targets[int(pos)]),
                "N_gallery": int(n_gallery),
                "N_gallery_max": int(n_gallery),
                "gallery_group_ids_hash": gallery_hash,
                "m_pos": 1,
                "rank_true": int(rank_true),
                "mrr_raw": float(mrr_raw),
                "hit1_raw": float(hit_raw[1]),
                "hit5_raw": float(hit_raw[5]),
                "hit10_raw": float(hit_raw[10]),
                "expected_mrr_chance": float(expected_mrr),
                "expected_hit1_chance": float(expected_hit[1]),
                "expected_hit5_chance": float(expected_hit[5]),
                "expected_hit10_chance": float(expected_hit[10]),
                "mrr_corrected": float(mrr_corr),
                "hit1_corrected": float(hit_corr[1]),
                "hit5_corrected": float(hit_corr[5]),
                "hit10_corrected": float(hit_corr[10]),
                "loo_policy": "disjoint_no_leakage",
                "cross_alignment_contract": CROSS_ALIGNMENT_CONTRACT,
            }
        )

    return {
        "scope": "cross",
        "dataset_or_direction": direction,
        "perturbation_type": rep["perturbation_type"],
        "representation": rep["representation"],
        "n_total": int(n_total),
        "n_valid": int(n_valid),
        "n_missing_metric": int(n_missing_metric),
        "n_mpos0": int(n_mpos0),
        "n_gallery": int(n_gallery),
        "sum_mrr_raw": float(sum_mrr_raw),
        "sum_exp_mrr": float(sum_exp_mrr),
        "sum_mrr_corr": float(sum_mrr_corr),
        "sum_hit1_raw": float(sum_hit_raw[1]),
        "sum_hit5_raw": float(sum_hit_raw[5]),
        "sum_hit10_raw": float(sum_hit_raw[10]),
        "sum_hit1_exp": float(sum_hit_exp[1]),
        "sum_hit5_exp": float(sum_hit_exp[5]),
        "sum_hit10_exp": float(sum_hit_exp[10]),
        "sum_hit1_corr": float(sum_hit_corr[1]),
        "sum_hit5_corr": float(sum_hit_corr[5]),
        "sum_hit10_corr": float(sum_hit_corr[10]),
        "rows": rows,
    }


def build_leaderboard_rows(
    summary_rows: Sequence[Mapping[str, object]],
    group_rows: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    leaderboard: List[Dict[str, object]] = []

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        retrieval_metric_map = {
            "mean_mrr_corrected": "mean_mrr_corrected",
            "mean_hit1_corrected": "mean_hit1_corrected",
            "mean_hit5_corrected": "mean_hit5_corrected",
            "mean_hit10_corrected": "mean_hit10_corrected",
        }
        for _, row in summary_df.iterrows():
            for metric_name, src in retrieval_metric_map.items():
                leaderboard.append(
                    {
                        "scope": row["scope"],
                        "dataset_or_direction": row["dataset_or_direction"],
                        "perturbation_type": row["perturbation_type"],
                        "representation": row["representation"],
                        "metric_name": metric_name,
                        "metric_value": row[src],
                        "n_total": int(row["n_total"]),
                        "n_valid": int(row["n_valid"]),
                        "n_excluded": int(row["n_excluded_missing_metric_or_mpos0"]),
                        "N_gallery_max": int(row["N_gallery_max"]),
                        "cross_alignment_contract": CROSS_ALIGNMENT_CONTRACT,
                    }
                )

    group_df = pd.DataFrame(group_rows)
    if not group_df.empty:
        metric_map = {
            "mean_cosine_centroid": "cosine",
            "mean_pcc_centroid": "pcc",
            "mean_edist_biascorr": "edist",
        }
        key_cols = ["scope", "dataset_or_direction", "perturbation_type", "representation"]
        for key_values, sub in group_df.groupby(key_cols, sort=False):
            key_map = dict(zip(key_cols, key_values))
            n_total = int(len(sub))
            for metric_name, metric_col in metric_map.items():
                finite = sub[metric_col].dropna()
                n_valid = int(len(finite))
                n_excluded = int(n_total - n_valid)
                metric_value = float(finite.mean()) if n_valid > 0 else np.nan
                leaderboard.append(
                    {
                        "scope": key_map["scope"],
                        "dataset_or_direction": key_map["dataset_or_direction"],
                        "perturbation_type": key_map["perturbation_type"],
                        "representation": key_map["representation"],
                        "metric_name": metric_name,
                        "metric_value": metric_value,
                        "n_total": n_total,
                        "n_valid": n_valid,
                        "n_excluded": n_excluded,
                        "N_gallery_max": 0,
                        "cross_alignment_contract": CROSS_ALIGNMENT_CONTRACT,
                    }
                )

    return leaderboard


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
        return 12

    task1_snapshot = (project_root / config["paths"]["task1_snapshot"]).resolve()
    runs_dir = (project_root / config["paths"]["runs_dir"]).resolve()

    expected_snapshot = (project_root / EXPECTED_TASK1_SNAPSHOT).resolve()
    if task1_snapshot != expected_snapshot:
        print(
            "[ERROR] Data isolation violation: config.paths.task1_snapshot must resolve "
            f"to {expected_snapshot}, got {task1_snapshot}",
            file=sys.stderr,
        )
        return 4

    if not task1_snapshot.is_dir():
        print(f"[ERROR] Missing task1 snapshot directory: {task1_snapshot}", file=sys.stderr)
        return 5

    stage_dir = runs_dir / args.run_id / STAGE
    stage_dir.mkdir(parents=True, exist_ok=True)

    started_at = utc_now_iso()
    rng = init_global_rng(GLOBAL_SEED)

    assertions: List[Dict[str, object]] = []
    input_paths: List[Path] = []
    attrition_bucket: Dict[Tuple[str, str, str, str, str], Dict[str, object]] = {}

    assertions.append(
        {
            "name": "task1_snapshot_isolation",
            "pass": True,
            "details": {
                "rules": [
                    "S2 reads exclusively from data/task1_snapshot_v1/",
                    "config.paths.task1_snapshot must equal data/task1_snapshot_v1",
                ],
                "task1_snapshot": str(task1_snapshot),
            },
            "counterexamples": [],
        }
    )

    try:
        lincs_meta, lincs_store, pathway_w, pathway_policy, lincs_paths = load_lincs_inputs(task1_snapshot)
        crispr_meta, crispr_gene_store, crispr_pathway_store, crispr_paths = load_scperturb_crispr_inputs(
            task1_snapshot
        )
        cross_contract_path = locate_cross_contract_file(task1_snapshot)
        cross_contract_df = pd.read_csv(cross_contract_path)

        input_paths.extend(lincs_paths)
        input_paths.extend(crispr_paths)
        input_paths.append(cross_contract_path)

        summary_path = task1_snapshot / "cross_contract/cross-pairs-alignment-summary.json"
        if summary_path.is_file():
            input_paths.append(summary_path)

        validated_pairs, validation_details = validate_cross_contract(
            cross_contract_df, lincs_meta, crispr_meta
        )
    except Exception as exc:
        assertions.append(
            {
                "name": "cross_contract_validation_fail_fast",
                "pass": False,
                "details": {
                    "rules": [
                        "Section 1.2: expected_cell_line/expected_target must match LINCS and scPerturb metadata exactly",
                        "Cross perturbation_type must match across sides",
                    ],
                    "error": str(exc),
                },
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Cross-contract fail-fast: {exc}", file=sys.stderr)
        return 6

    assertions.append(
        {
            "name": "cross_contract_validation_fail_fast",
            "pass": True,
            "details": {
                "rows_validated": int(validation_details["n_rows"]),
                "n_unique_lincs_indices": int(validation_details["n_unique_lincs_indices"]),
                "n_unique_sc_indices": int(validation_details["n_unique_sc_indices"]),
                "cross_alignment_contract": CROSS_ALIGNMENT_CONTRACT,
            },
            "counterexamples": [],
        }
    )

    mode_ok = str(pathway_policy.get("mode", "")).strip() == "project_on_load"
    assertions.append(
        {
            "name": "lincs_pathway_project_on_load_policy",
            "pass": mode_ok,
            "details": {
                "rules": [
                    "LINCS pathway must use project_on_load policy",
                ],
                "mode": pathway_policy.get("mode", "NA"),
                "W_path": pathway_policy.get("W_path", "NA"),
            },
            "counterexamples": []
            if mode_ok
            else [{"mode": pathway_policy.get("mode", "NA")}],
        }
    )
    if not mode_ok:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] LINCS pathway policy mode must be project_on_load", file=sys.stderr)
        return 6

    proof_rows: List[Dict[str, object]] = []
    eligible_perturbation_types: List[str] = []
    for perturbation_type in ["Chemical", "Genetic"]:
        subset = validated_pairs.loc[validated_pairs["perturbation_type"] == perturbation_type]
        n_matched_keys = int(
            subset[["cell_line", "target", "perturbation_type"]].drop_duplicates().shape[0]
        )
        eligible = n_matched_keys >= 5
        excluded_reason = "" if eligible else "matched_keys_lt5"

        if eligible:
            eligible_perturbation_types.append(perturbation_type)
        else:
            record_attrition(
                attrition_bucket,
                scope="cross",
                dataset_or_direction="cross",
                perturbation_type=perturbation_type,
                representation="ALL",
                reason="cross_matched_keys_lt_5",
                n_dropped=n_matched_keys,
                n_total_before=n_matched_keys,
                notes="Cross eligibility gate: n_matched_keys < 5",
                allow_zero=True,
            )

        proof_rows.append(
            {
                "cross_alignment_contract": CROSS_ALIGNMENT_CONTRACT,
                "perturbation_type": perturbation_type,
                "n_matched_keys": n_matched_keys,
                "eligible_bool": bool(eligible),
                "excluded_reason": excluded_reason,
            }
        )

    alignment_proof_df = pd.DataFrame(proof_rows).sort_values(
        ["perturbation_type"], kind="mergesort"
    )
    alignment_proof_path = stage_dir / "task1_cross_alignment_proof.csv"
    write_csv(alignment_proof_df, alignment_proof_path)

    assertions.append(
        {
            "name": "cross_eligibility_gate",
            "pass": True,
            "details": {
                "rules": [
                    "Per perturbation_type: eligible iff n_matched_keys >= 5",
                    "Ineligible perturbation_type must be excluded from metrics",
                ],
                "proof_rows": alignment_proof_df.to_dict(orient="records"),
            },
            "counterexamples": [],
        }
    )

    worker_payload: Dict[str, Dict[str, Any]] = {}
    try:
        if eligible_perturbation_types:
            worker_payload = build_worker_payload(
                validated_pairs=validated_pairs,
                lincs_store=lincs_store,
                sc_gene_store=crispr_gene_store,
                sc_pathway_store=crispr_pathway_store,
                pathway_w=pathway_w,
                eligible_perturbation_types=eligible_perturbation_types,
            )
    except Exception as exc:
        assertions.append(
            {
                "name": "cross_representation_payload_build",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed building cross representation payload: {exc}", file=sys.stderr)
        return 6

    group_tasks, random_draws = precompute_group_tasks(worker_payload, rng)
    retrieval_tasks = precompute_retrieval_tasks(worker_payload, QUERY_CHUNK_SIZE)

    assertions.append(
        {
            "name": "deterministic_parallelism_setup",
            "pass": True,
            "details": {
                "rules": [
                    "Only main process RNG (seed=619) precomputes all randomized subsample indices",
                    "Workers do not call RNG",
                    "Lightweight task payloads (indices/ids) only",
                ],
                "workers": workers,
                "n_cohorts": int(len(worker_payload)),
                "n_group_tasks": int(len(group_tasks)),
                "n_retrieval_tasks": int(len(retrieval_tasks)),
                "n_rng_draws_main_process": int(random_draws),
            },
            "counterexamples": [],
        }
    )

    group_rows: List[Dict[str, object]] = []
    per_query_rows: List[Dict[str, object]] = []

    summary_acc: Dict[Tuple[str, str, str, str], Dict[str, object]] = {}

    try:
        if worker_payload and (group_tasks or retrieval_tasks):
            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=init_worker,
                initargs=(worker_payload,),
            ) as executor:
                group_futures = [executor.submit(worker_group_metric, task) for task in group_tasks]
                for future in as_completed(group_futures):
                    result = future.result()
                    row = result["row"]
                    events = result["events"]
                    group_rows.append(row)

                    for event in events:
                        record_attrition(
                            attrition_bucket,
                            scope="cross",
                            dataset_or_direction=GROUP_DATASET_OR_DIRECTION,
                            perturbation_type=str(row["perturbation_type"]),
                            representation=str(row["representation"]),
                            reason=str(event["reason"]),
                            n_dropped=int(event["n_dropped"]),
                            n_total_before=int(event["n_total_before"]),
                            notes=str(event["notes"]),
                        )

                retrieval_futures = [
                    executor.submit(worker_retrieval_chunk, task) for task in retrieval_tasks
                ]
                for future in as_completed(retrieval_futures):
                    result = future.result()
                    key = (
                        str(result["scope"]),
                        str(result["dataset_or_direction"]),
                        str(result["perturbation_type"]),
                        str(result["representation"]),
                    )
                    acc = summary_acc.get(key)
                    if acc is None:
                        acc = {
                            "scope": key[0],
                            "dataset_or_direction": key[1],
                            "perturbation_type": key[2],
                            "representation": key[3],
                            "n_total": 0,
                            "n_valid": 0,
                            "N_gallery_max": 0,
                            "sum_mrr_raw": 0.0,
                            "sum_exp_mrr": 0.0,
                            "sum_mrr_corr": 0.0,
                            "sum_hit1_raw": 0.0,
                            "sum_hit5_raw": 0.0,
                            "sum_hit10_raw": 0.0,
                            "sum_hit1_exp": 0.0,
                            "sum_hit5_exp": 0.0,
                            "sum_hit10_exp": 0.0,
                            "sum_hit1_corr": 0.0,
                            "sum_hit5_corr": 0.0,
                            "sum_hit10_corr": 0.0,
                        }
                        summary_acc[key] = acc

                    acc["n_total"] = int(acc["n_total"]) + int(result["n_total"])
                    acc["n_valid"] = int(acc["n_valid"]) + int(result["n_valid"])
                    acc["N_gallery_max"] = max(int(acc["N_gallery_max"]), int(result["n_gallery"]))

                    acc["sum_mrr_raw"] = float(acc["sum_mrr_raw"]) + float(result["sum_mrr_raw"])
                    acc["sum_exp_mrr"] = float(acc["sum_exp_mrr"]) + float(result["sum_exp_mrr"])
                    acc["sum_mrr_corr"] = float(acc["sum_mrr_corr"]) + float(result["sum_mrr_corr"])

                    acc["sum_hit1_raw"] = float(acc["sum_hit1_raw"]) + float(result["sum_hit1_raw"])
                    acc["sum_hit5_raw"] = float(acc["sum_hit5_raw"]) + float(result["sum_hit5_raw"])
                    acc["sum_hit10_raw"] = float(acc["sum_hit10_raw"]) + float(result["sum_hit10_raw"])
                    acc["sum_hit1_exp"] = float(acc["sum_hit1_exp"]) + float(result["sum_hit1_exp"])
                    acc["sum_hit5_exp"] = float(acc["sum_hit5_exp"]) + float(result["sum_hit5_exp"])
                    acc["sum_hit10_exp"] = float(acc["sum_hit10_exp"]) + float(result["sum_hit10_exp"])
                    acc["sum_hit1_corr"] = float(acc["sum_hit1_corr"]) + float(result["sum_hit1_corr"])
                    acc["sum_hit5_corr"] = float(acc["sum_hit5_corr"]) + float(result["sum_hit5_corr"])
                    acc["sum_hit10_corr"] = float(acc["sum_hit10_corr"]) + float(result["sum_hit10_corr"])

                    missing_metric = int(result["n_missing_metric"])
                    m_pos0 = int(result["n_mpos0"])
                    if missing_metric > 0:
                        record_attrition(
                            attrition_bucket,
                            scope="cross",
                            dataset_or_direction=str(result["dataset_or_direction"]),
                            perturbation_type=str(result["perturbation_type"]),
                            representation=str(result["representation"]),
                            reason="missing_metric",
                            n_dropped=missing_metric,
                            n_total_before=int(result["n_total"]),
                            notes="Query excluded due to non-finite vector",
                        )
                    if m_pos0 > 0:
                        record_attrition(
                            attrition_bucket,
                            scope="cross",
                            dataset_or_direction=str(result["dataset_or_direction"]),
                            perturbation_type=str(result["perturbation_type"]),
                            representation=str(result["representation"]),
                            reason="m_pos0_after_intersection",
                            n_dropped=m_pos0,
                            n_total_before=int(result["n_total"]),
                            notes="No positive group in disjoint cross gallery",
                        )

                    per_query_rows.extend(result["rows"])
    except Exception as exc:
        assertions.append(
            {
                "name": "parallel_compute",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Parallel compute failed: {exc}", file=sys.stderr)
        return 13

    summary_rows: List[Dict[str, object]] = []
    chance_rows: List[Dict[str, object]] = []

    for key in sorted(summary_acc.keys()):
        acc = summary_acc[key]
        n_total = int(acc["n_total"])
        n_valid = int(acc["n_valid"])
        n_excluded = int(n_total - n_valid)
        n_gallery_max = int(acc["N_gallery_max"])

        if n_valid > 0:
            mean_mrr_raw = float(acc["sum_mrr_raw"]) / n_valid
            mean_exp_mrr = float(acc["sum_exp_mrr"]) / n_valid
            mean_mrr_corr = float(acc["sum_mrr_corr"]) / n_valid

            mean_hit1_raw = float(acc["sum_hit1_raw"]) / n_valid
            mean_hit5_raw = float(acc["sum_hit5_raw"]) / n_valid
            mean_hit10_raw = float(acc["sum_hit10_raw"]) / n_valid

            mean_hit1_exp = float(acc["sum_hit1_exp"]) / n_valid
            mean_hit5_exp = float(acc["sum_hit5_exp"]) / n_valid
            mean_hit10_exp = float(acc["sum_hit10_exp"]) / n_valid

            mean_hit1_corr = float(acc["sum_hit1_corr"]) / n_valid
            mean_hit5_corr = float(acc["sum_hit5_corr"]) / n_valid
            mean_hit10_corr = float(acc["sum_hit10_corr"]) / n_valid

            delta_mrr = float(mean_mrr_corr - (mean_mrr_raw - mean_exp_mrr))
            delta_hit1 = float(mean_hit1_corr - (mean_hit1_raw - mean_hit1_exp))
            delta_hit5 = float(mean_hit5_corr - (mean_hit5_raw - mean_hit5_exp))
            delta_hit10 = float(mean_hit10_corr - (mean_hit10_raw - mean_hit10_exp))
        else:
            mean_mrr_raw = np.nan
            mean_exp_mrr = np.nan
            mean_mrr_corr = np.nan
            mean_hit1_raw = np.nan
            mean_hit5_raw = np.nan
            mean_hit10_raw = np.nan
            mean_hit1_exp = np.nan
            mean_hit5_exp = np.nan
            mean_hit10_exp = np.nan
            mean_hit1_corr = np.nan
            mean_hit5_corr = np.nan
            mean_hit10_corr = np.nan
            delta_mrr = np.nan
            delta_hit1 = np.nan
            delta_hit5 = np.nan
            delta_hit10 = np.nan

        summary_rows.append(
            {
                "scope": acc["scope"],
                "dataset_or_direction": acc["dataset_or_direction"],
                "perturbation_type": acc["perturbation_type"],
                "representation": acc["representation"],
                "n_total": n_total,
                "n_valid": n_valid,
                "n_excluded_missing_metric_or_mpos0": n_excluded,
                "N_gallery_max": n_gallery_max,
                "mean_mrr_raw": mean_mrr_raw,
                "mean_expected_mrr_chance": mean_exp_mrr,
                "mean_mrr_corrected": mean_mrr_corr,
                "mean_hit1_raw": mean_hit1_raw,
                "mean_expected_hit1_chance": mean_hit1_exp,
                "mean_hit1_corrected": mean_hit1_corr,
                "mean_hit5_raw": mean_hit5_raw,
                "mean_expected_hit5_chance": mean_hit5_exp,
                "mean_hit5_corrected": mean_hit5_corr,
                "mean_hit10_raw": mean_hit10_raw,
                "mean_expected_hit10_chance": mean_hit10_exp,
                "mean_hit10_corrected": mean_hit10_corr,
            }
        )

        chance_rows.append(
            {
                "scope": acc["scope"],
                "dataset_or_direction": acc["dataset_or_direction"],
                "perturbation_type": acc["perturbation_type"],
                "representation": acc["representation"],
                "delta_mrr": delta_mrr,
                "abs_delta_mrr": abs(delta_mrr) if np.isfinite(delta_mrr) else np.nan,
                "delta_hit1": delta_hit1,
                "abs_delta_hit1": abs(delta_hit1) if np.isfinite(delta_hit1) else np.nan,
                "delta_hit5": delta_hit5,
                "abs_delta_hit5": abs(delta_hit5) if np.isfinite(delta_hit5) else np.nan,
                "delta_hit10": delta_hit10,
                "abs_delta_hit10": abs(delta_hit10) if np.isfinite(delta_hit10) else np.nan,
            }
        )

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows).sort_values(
            ["scope", "dataset_or_direction", "perturbation_type", "representation"],
            kind="mergesort",
        )
    else:
        summary_df = pd.DataFrame(columns=SUMMARY_COLUMNS)

    if chance_rows:
        chance_df = pd.DataFrame(chance_rows).sort_values(
            ["scope", "dataset_or_direction", "perturbation_type", "representation"],
            kind="mergesort",
        )
    else:
        chance_df = pd.DataFrame(columns=CHANCE_COLUMNS)

    if group_rows:
        group_df = pd.DataFrame(group_rows).sort_values(
            ["scope", "dataset_or_direction", "perturbation_type", "representation", "group_id"],
            kind="mergesort",
        )
    else:
        group_df = pd.DataFrame(columns=GROUP_COLUMNS)

    if per_query_rows:
        per_query_df = pd.DataFrame(per_query_rows).sort_values(
            ["scope", "dataset_or_direction", "perturbation_type", "representation", "query_uid"],
            kind="mergesort",
        )
    else:
        per_query_df = pd.DataFrame(columns=[f.name for f in per_query_schema()])

    attrition_df = pd.DataFrame(list(attrition_bucket.values()))
    if attrition_df.empty:
        attrition_df = pd.DataFrame(columns=ATTRITION_COLUMNS)
    else:
        attrition_df = attrition_df.sort_values(
            ["scope", "dataset_or_direction", "perturbation_type", "representation", "reason"],
            kind="mergesort",
        )

    leaderboard_rows = build_leaderboard_rows(summary_rows=summary_rows, group_rows=group_rows)
    if leaderboard_rows:
        leaderboard_df = pd.DataFrame(leaderboard_rows).sort_values(
            ["scope", "dataset_or_direction", "perturbation_type", "representation", "metric_name"],
            kind="mergesort",
        )
    else:
        leaderboard_df = pd.DataFrame(columns=LEADERBOARD_COLUMNS)

    per_query_path = stage_dir / "task1_cross_retrieval_per_query.parquet"
    summary_path = stage_dir / "task1_cross_retrieval_summary.csv"
    chance_path = stage_dir / "task1_cross_chance_identity_check.csv"
    leaderboard_path = stage_dir / "task1_cross_leaderboard_long.csv"
    attrition_path = stage_dir / "task1_cross_attrition.csv"
    group_cross_path = stage_dir / "task1_group_cross.parquet"

    write_parquet_with_schema(
        per_query_df,
        per_query_path,
        per_query_schema(),
        canonicalize_per_query_frame,
    )
    write_csv(summary_df, summary_path)
    write_csv(chance_df, chance_path)
    write_csv(leaderboard_df, leaderboard_path)
    write_csv(attrition_df, attrition_path)
    write_parquet_with_schema(
        group_df,
        group_cross_path,
        group_cross_schema(),
        canonicalize_group_frame,
    )

    denominator_failures = summary_df.loc[
        summary_df["n_total"]
        != (summary_df["n_valid"] + summary_df["n_excluded_missing_metric_or_mpos0"])
    ]
    denominator_pass = denominator_failures.empty
    assertions.append(
        {
            "name": "denominator_conservation",
            "pass": bool(denominator_pass),
            "details": {
                "rules": [
                    "For each cohort: n_total == n_valid + n_excluded_missing_metric_or_mpos0",
                ],
                "n_checked": int(len(summary_df)),
            },
            "counterexamples": denominator_failures.head(MAX_COUNTEREXAMPLES).to_dict(
                orient="records"
            ),
        }
    )

    tol_cols = ["abs_delta_mrr", "abs_delta_hit1", "abs_delta_hit5", "abs_delta_hit10"]
    chance_finite = chance_df.dropna(subset=tol_cols, how="any")
    if chance_finite.empty:
        chance_violations = pd.DataFrame(columns=chance_df.columns)
    else:
        chance_violations = chance_finite.loc[
            (chance_finite["abs_delta_mrr"] > CHANCE_IDENTITY_TOL)
            | (chance_finite["abs_delta_hit1"] > CHANCE_IDENTITY_TOL)
            | (chance_finite["abs_delta_hit5"] > CHANCE_IDENTITY_TOL)
            | (chance_finite["abs_delta_hit10"] > CHANCE_IDENTITY_TOL)
        ]

    chance_pass = chance_violations.empty
    assertions.append(
        {
            "name": "chance_identity_tolerance",
            "pass": bool(chance_pass),
            "details": {
                "rules": [
                    "For each finite cohort: abs(delta_mrr) <= 1e-12",
                    "For each finite cohort: abs(delta_hit1/hit5/hit10) <= 1e-12",
                ],
                "tolerance": CHANCE_IDENTITY_TOL,
                "n_checked": int(len(chance_finite)),
            },
            "counterexamples": chance_violations.head(MAX_COUNTEREXAMPLES).to_dict(
                orient="records"
            ),
        }
    )

    core_outputs = [
        per_query_path,
        summary_path,
        chance_path,
        leaderboard_path,
        alignment_proof_path,
        attrition_path,
        group_cross_path,
    ]
    output_routing_pass = all(
        path.resolve().is_relative_to(stage_dir.resolve()) for path in core_outputs
    )
    assertions.append(
        {
            "name": "output_routing_isolated",
            "pass": bool(output_routing_pass),
            "details": {
                "rules": [
                    f"All outputs must be written under {stage_dir}",
                ],
                "stage_dir": str(stage_dir),
            },
            "counterexamples": []
            if output_routing_pass
            else [
                {"bad_path": str(path)}
                for path in core_outputs
                if not path.resolve().is_relative_to(stage_dir.resolve())
            ],
        }
    )

    bad_inputs = [
        str(path)
        for path in input_paths
        if not path.absolute().is_relative_to(task1_snapshot.absolute())
    ]
    assertions.append(
        {
            "name": "input_path_data_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "rules": [
                    "All loaded inputs must be under data/task1_snapshot_v1/ lexical namespace",
                ],
                "task1_snapshot": str(task1_snapshot),
                "n_inputs": int(len(input_paths)),
            },
            "counterexamples": [{"bad_input": p} for p in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    completed_at = utc_now_iso()

    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    manifest_path = stage_dir / "manifest.json"

    output_paths = [
        str(path.resolve())
        for path in core_outputs + [run_manifest_path, audit_assertions_path, manifest_path]
    ]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s2_task1_cross_metrics.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "edist_max_n": EDIST_MAX_N,
            "hit_k_list": list(HIT_K_LIST),
            "chance_identity_tol": CHANCE_IDENTITY_TOL,
            "cross_alignment_contract": CROSS_ALIGNMENT_CONTRACT,
            "task1_snapshot": str(task1_snapshot),
            "runs_dir": str(runs_dir),
            "workers": workers,
            "workers_source": "cli" if args.workers is not None else "config.compute.workers",
            "query_chunk_size": QUERY_CHUNK_SIZE,
            "rng_draws_main_process": int(random_draws),
        },
        "inputs": [str(path) for path in sorted(input_paths)],
        "outputs": output_paths,
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

    if not denominator_pass:
        print("[ERROR] denominator conservation assertion failed", file=sys.stderr)
        return 7

    if not output_routing_pass:
        print("[ERROR] output routing assertion failed", file=sys.stderr)
        return 9

    if not chance_pass:
        print(
            f"[ERROR] chance identity check failed: tolerance={CHANCE_IDENTITY_TOL}",
            file=sys.stderr,
        )
        return 10

    if bad_inputs:
        print("[ERROR] input path data isolation assertion failed", file=sys.stderr)
        return 11

    return 0


if __name__ == "__main__":
    sys.exit(main())