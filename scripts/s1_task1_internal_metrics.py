# SCRIPT_HEADER_CONTRACT
# Script: scripts/s1_task1_internal_metrics.py
# Purpose: 计算 Task1 internal 组级一致性与检索指标 (Modality Concordance)
# Inputs:
#   - Task1 Snapshot: data/task1_snapshot_v1/
# Outputs:
#   - task1_retrieval_per_query.parquet: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_retrieval_summary.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_chance_identity_check.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_leaderboard_long.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_attrition.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s1_task1_internal_metrics/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task1_snapshot
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s1_task1_internal_metrics.py --run-id <run_id> --seed 619
# Failure Modes:
#   - LOO strict_recompute violated -> exit non-zero
#   - Chance identity check tolerance > 1e-12 -> exit non-zero
# Last Updated: 2026-03-03

"""
Inputs:
- task1 snapshot files under data/task1_snapshot_v1/ only:
  - lincs/lincs-engine1-meta.csv
  - lincs/lincs-engine1-gene-delta.pt
  - scperturb_delta/*.csv, *.npy
  - fm_delta/<model>/*.csv, *.npy
  - pathway/hallmark-w-2477x50.npy
  - pathway/lincs-pathway-policy.json

Outputs:
- task1_retrieval_per_query.parquet
- task1_retrieval_summary.csv
- task1_chance_identity_check.csv
- task1_leaderboard_long.csv
- task1_attrition.csv
- run_manifest.json
- audit_assertions.json
- manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- EDIST_MAX_N = 256
- HIT_K_LIST = (1, 5, 10)
- CHANCE_IDENTITY_TOL = 1e-12
- LOO policy = strict_recompute
- LINCS pathway projection policy = project_on_load using hallmark-w-2477x50.npy

Attrition rules:
- invalid_delta
- split_half_requires_n>=4
- edist_requires_n>=2_per_split
- loo_requires_count>1
- missing_metric
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import torch
import yaml

STAGE = "s1_task1_internal_metrics"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK1_SNAPSHOT = Path("data/task1_snapshot_v1")

GLOBAL_SEED = 619
EDIST_MAX_N = 256
HIT_K_LIST: Tuple[int, ...] = (1, 5, 10)
CHANCE_IDENTITY_TOL = 1e-12
MAX_COUNTEREXAMPLES = 5


@dataclass(frozen=True)
class CohortKey:
    scope: str
    dataset_or_direction: str
    perturbation_type: str
    representation: str


@dataclass
class CohortData:
    key: CohortKey
    meta: pd.DataFrame
    vector_store: "VectorStore"


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


class ProjectedPathwayStore(VectorStore):
    """
    Project pathway on load from gene delta using frozen W.
    This enforces: pathway_delta = gene_delta @ hallmark_W_2477x50.
    """

    def __init__(self, gene_store: VectorStore, projection_w: np.ndarray) -> None:
        self.gene_store = gene_store
        self.w = np.asarray(projection_w, dtype=np.float64)

    def load_rows(self, row_indices: np.ndarray) -> np.ndarray:
        gene = self.gene_store.load_rows(row_indices)
        if gene.size == 0:
            return np.empty((0, self.w.shape[1]), dtype=np.float64)
        return gene @ self.w


class PerQueryParquetSink:
    def __init__(self, output_path: Path, schema: pa.Schema) -> None:
        self.output_path = output_path
        self.schema = schema
        self.writer: Optional[pq.ParquetWriter] = None

    def write(self, frame: pd.DataFrame) -> None:
        if frame.empty:
            return
        table = pa.Table.from_pandas(frame, schema=self.schema, preserve_index=False)
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)
        self.writer.write_table(table)

    def close(self) -> None:
        if self.writer is not None:
            self.writer.close()
            return

        empty_data = {
            field.name: pa.array([], type=field.type) for field in self.schema
        }
        empty_table = pa.Table.from_pydict(empty_data, schema=self.schema)
        pq.write_table(empty_table, self.output_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S1 Task1 internal metrics")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--n-gallery-max",
        type=int,
        default=None,
        help="Optional deterministic cap for candidate centroids in retrieval",
    )
    parser.add_argument(
        "--query-batch-size",
        type=int,
        default=256,
        help="Batch size for query scoring",
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


def normalize_perturbation_type(raw_value: str) -> str:
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


def parse_bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    mapped = series.fillna(False).astype(str).str.strip().str.lower().map(
        {
            "true": True,
            "1": True,
            "yes": True,
            "y": True,
            "false": False,
            "0": False,
            "no": False,
            "n": False,
        }
    )
    return mapped.fillna(False).astype(bool)


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


def deterministic_subsample_by_uid(
    vectors: np.ndarray,
    uids: np.ndarray,
    max_n: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, int]:
    if vectors.shape[0] <= max_n:
        return vectors, int(vectors.shape[0])

    uid_order = np.argsort(uids.astype(str), kind="mergesort")
    vec_sorted = vectors[uid_order]
    perm = rng.permutation(vec_sorted.shape[0])
    picked = perm[:max_n]
    return vec_sorted[picked], int(max_n)


def expected_mrr_for_single_positive(n_gallery: int) -> float:
    ranks = np.arange(1, n_gallery + 1, dtype=np.float64)
    return float(np.mean(1.0 / ranks))


def expected_hitk_for_single_positive(n_gallery: int, k: int) -> float:
    return float(min(k, n_gallery) / n_gallery)


def hash_gallery_group_ids(group_ids: Sequence[str]) -> str:
    material = "||".join(group_ids)
    return hashlib.sha256(material.encode("utf-8")).hexdigest()


def make_group_id(cell_line: str, target: str) -> str:
    return f"{cell_line}||{target}"


def record_attrition(
    bucket: MutableMapping[Tuple[str, str, str, str, str], Dict[str, object]],
    key: CohortKey,
    reason: str,
    n_dropped: int,
    n_total_before: int,
    notes: str,
) -> None:
    if n_dropped <= 0:
        return

    agg_key = (
        key.scope,
        key.dataset_or_direction,
        key.perturbation_type,
        key.representation,
        reason,
    )
    if agg_key not in bucket:
        bucket[agg_key] = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
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


def load_lincs_inputs(snapshot_root: Path) -> Tuple[pd.DataFrame, VectorStore, np.ndarray, Dict[str, object], List[Path]]:
    lincs_meta_path = snapshot_root / "lincs/lincs-engine1-meta.csv"
    lincs_gene_pt_path = snapshot_root / "lincs/lincs-engine1-gene-delta.pt"
    pathway_w_path = snapshot_root / "pathway/hallmark-w-2477x50.npy"
    pathway_policy_path = snapshot_root / "pathway/lincs-pathway-policy.json"

    lincs_meta = pd.read_csv(lincs_meta_path)
    ensure_required_columns(
        lincs_meta,
        ["pert_type", "target", "cell_line", "dose_val", "time_val"],
        "lincs-engine1-meta.csv",
    )

    lincs_payload = torch.load(lincs_gene_pt_path, map_location="cpu")
    if not isinstance(lincs_payload, dict) or "y_delta_gene" not in lincs_payload:
        raise ValueError("lincs-engine1-gene-delta.pt must be a dict containing y_delta_gene")

    y_delta_gene = lincs_payload["y_delta_gene"]
    if not torch.is_tensor(y_delta_gene):
        raise ValueError("lincs y_delta_gene must be a torch.Tensor")

    if int(y_delta_gene.shape[0]) != int(len(lincs_meta)):
        raise ValueError(
            f"LINCS row mismatch: tensor rows={int(y_delta_gene.shape[0])}, meta rows={len(lincs_meta)}"
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

    lincs_store = TorchMatrixStore(y_delta_gene)
    input_paths = [lincs_meta_path, lincs_gene_pt_path, pathway_w_path, pathway_policy_path]
    return lincs_meta, lincs_store, np.asarray(projection_w, dtype=np.float64), pathway_policy, input_paths


def load_scperturb_inputs(snapshot_root: Path) -> Tuple[Dict[str, pd.DataFrame], Dict[str, VectorStore], Dict[str, VectorStore], List[Path]]:
    drug_meta_path = snapshot_root / "scperturb_delta/scperturb-drug-delta-meta.csv"
    crispr_meta_path = snapshot_root / "scperturb_delta/scperturb-crispr-delta-meta.csv"

    drug_gene_path = snapshot_root / "scperturb_delta/scperturb-drug-gene-delta.npy"
    drug_pathway_path = snapshot_root / "scperturb_delta/scperturb-drug-pathway-delta.npy"
    crispr_gene_path = snapshot_root / "scperturb_delta/scperturb-crispr-gene-delta.npy"
    crispr_pathway_path = snapshot_root / "scperturb_delta/scperturb-crispr-pathway-delta.npy"

    required_meta_cols = ["delta_row_idx", "target_std", "cell_std", "dose_val", "time_val", "pert_type"]
    drug_meta = pd.read_csv(drug_meta_path)
    crispr_meta = pd.read_csv(crispr_meta_path)
    ensure_required_columns(drug_meta, required_meta_cols, "scperturb-drug-delta-meta.csv")
    ensure_required_columns(crispr_meta, required_meta_cols, "scperturb-crispr-delta-meta.csv")

    for name, frame in [("drug", drug_meta), ("crispr", crispr_meta)]:
        if frame["delta_row_idx"].duplicated().any():
            dup = frame.loc[frame["delta_row_idx"].duplicated(), "delta_row_idx"].head(MAX_COUNTEREXAMPLES)
            raise ValueError(f"{name} delta_row_idx contains duplicates, examples={dup.tolist()}")

    drug_gene = NpyMatrixStore(drug_gene_path)
    drug_pathway = NpyMatrixStore(drug_pathway_path)
    crispr_gene = NpyMatrixStore(crispr_gene_path)
    crispr_pathway = NpyMatrixStore(crispr_pathway_path)

    # Row-range checks
    drug_idx = drug_meta["delta_row_idx"].to_numpy(dtype=np.int64)
    crispr_idx = crispr_meta["delta_row_idx"].to_numpy(dtype=np.int64)

    if np.min(drug_idx) < 0 or np.max(drug_idx) >= int(drug_gene.arr.shape[0]):
        raise ValueError("drug delta_row_idx out of gene array bounds")
    if np.min(drug_idx) < 0 or np.max(drug_idx) >= int(drug_pathway.arr.shape[0]):
        raise ValueError("drug delta_row_idx out of pathway array bounds")
    if np.min(crispr_idx) < 0 or np.max(crispr_idx) >= int(crispr_gene.arr.shape[0]):
        raise ValueError("crispr delta_row_idx out of gene array bounds")
    if np.min(crispr_idx) < 0 or np.max(crispr_idx) >= int(crispr_pathway.arr.shape[0]):
        raise ValueError("crispr delta_row_idx out of pathway array bounds")

    meta_map = {
        "Chemical": drug_meta,
        "Genetic": crispr_meta,
    }
    gene_store_map = {
        "Chemical": drug_gene,
        "Genetic": crispr_gene,
    }
    pathway_store_map = {
        "Chemical": drug_pathway,
        "Genetic": crispr_pathway,
    }
    input_paths = [
        drug_meta_path,
        crispr_meta_path,
        drug_gene_path,
        drug_pathway_path,
        crispr_gene_path,
        crispr_pathway_path,
    ]
    return meta_map, gene_store_map, pathway_store_map, input_paths


def load_fm_inputs(
    snapshot_root: Path,
    sc_meta_map: Mapping[str, pd.DataFrame],
) -> Tuple[List[CohortData], List[Path], List[Dict[str, object]]]:
    fm_root = snapshot_root / "fm_delta"
    cohort_list: List[CohortData] = []
    input_paths: List[Path] = []
    policy_rows: List[Dict[str, object]] = []

    if not fm_root.exists():
        return cohort_list, input_paths, policy_rows

    sc_row_sets: Dict[str, set] = {}
    for perturbation_type, frame in sc_meta_map.items():
        sc_row_sets[perturbation_type] = set(frame["delta_row_idx"].astype(int).tolist())

    for model_dir in sorted([p for p in fm_root.iterdir() if p.is_dir()]):
        model = model_dir.name
        policy_path = model_dir / "delta_operator_policy.json"
        policy_info: Dict[str, object] = {
            "model": model,
            "policy_found": policy_path.exists(),
            "policy_path": str(policy_path),
            "operator_type": "NA",
            "eps": "NA",
        }

        if policy_path.exists():
            with policy_path.open("r", encoding="utf-8") as handle:
                policy = json.load(handle)
            policy_info["operator_type"] = str(policy.get("operator_type", "NA"))
            policy_info["eps"] = policy.get("eps", "NA")
            input_paths.append(policy_path)

        policy_rows.append(policy_info)

        for perturbation_type in ["Chemical", "Genetic"]:
            meta_path = model_dir / f"{perturbation_type}_delta_meta.csv"
            delta_path = model_dir / f"{perturbation_type}_delta_cell.npy"
            if not meta_path.is_file() or not delta_path.is_file():
                raise ValueError(f"Missing FM files for model={model}, perturbation_type={perturbation_type}")

            meta = pd.read_csv(meta_path)
            required_cols = ["row_in_modality_matrix", "target_std", "cell_std", "dose_val", "time_val", "pert_type"]
            ensure_required_columns(meta, required_cols, f"{model}/{perturbation_type}_delta_meta.csv")

            if meta["row_in_modality_matrix"].duplicated().any():
                dup = (
                    meta.loc[
                        meta["row_in_modality_matrix"].duplicated(),
                        "row_in_modality_matrix",
                    ]
                    .head(MAX_COUNTEREXAMPLES)
                    .tolist()
                )
                raise ValueError(
                    f"FM row_in_modality_matrix has duplicates for model={model}, perturbation_type={perturbation_type}, examples={dup}"
                )

            rows = meta["row_in_modality_matrix"].astype(np.int64)
            if not set(rows.tolist()).issubset(sc_row_sets[perturbation_type]):
                raise ValueError(
                    f"FM rows must be subset of sc delta_row_idx for {model}/{perturbation_type}"
                )

            store = NpyMatrixStore(delta_path)
            if int(np.min(rows)) < 0 or int(np.max(rows)) >= int(store.arr.shape[0]):
                raise ValueError(f"FM row_in_modality_matrix out of bounds for {model}/{perturbation_type}")

            normalized = pd.DataFrame(
                {
                    "cell_line": normalize_text(meta["cell_std"]),
                    "target": normalize_text(meta["target_std"]),
                    "dose_val": meta["dose_val"],
                    "time_val": meta["time_val"],
                    "canonical_query_uid": meta["row_in_modality_matrix"].astype(np.int64).astype(str),
                    "delta_row_index": meta["row_in_modality_matrix"].astype(np.int64),
                    "delta_valid_bool": np.ones(len(meta), dtype=bool),
                }
            )

            valid_col = None
            for candidate in ["valid_mask", "delta_valid_bool"]:
                if candidate in meta.columns:
                    valid_col = candidate
                    break
            if valid_col is not None:
                normalized["delta_valid_bool"] = parse_bool_series(meta[valid_col])

            cohort_list.append(
                CohortData(
                    key=CohortKey(
                        scope="internal",
                        dataset_or_direction="scPerturb",
                        perturbation_type=perturbation_type,
                        representation=f"FM:{model}",
                    ),
                    meta=normalized,
                    vector_store=store,
                )
            )

            input_paths.extend([meta_path, delta_path])

    return cohort_list, input_paths, policy_rows


def build_cohorts(snapshot_root: Path) -> Tuple[List[CohortData], List[Path], Dict[str, object]]:
    lincs_meta, lincs_gene_store, pathway_w, pathway_policy, lincs_paths = load_lincs_inputs(snapshot_root)
    sc_meta_map, sc_gene_map, sc_pathway_map, sc_paths = load_scperturb_inputs(snapshot_root)

    cohorts: List[CohortData] = []
    input_paths: List[Path] = []

    # LINCS: gene + pathway(project_on_load), no FM.
    lincs_pert = lincs_meta["pert_type"].map(normalize_perturbation_type)
    lincs_base = pd.DataFrame(
        {
            "cell_line": normalize_text(lincs_meta["cell_line"]),
            "target": normalize_text(lincs_meta["target"]),
            "dose_val": lincs_meta["dose_val"],
            "time_val": lincs_meta["time_val"],
            "perturbation_type": lincs_pert,
            "delta_row_index": np.arange(len(lincs_meta), dtype=np.int64),
        }
    )
    lincs_base["canonical_query_uid"] = lincs_base["delta_row_index"].astype(str)
    lincs_base["delta_valid_bool"] = True

    lincs_pathway_store = ProjectedPathwayStore(lincs_gene_store, pathway_w)

    for perturbation_type in ["Chemical", "Genetic"]:
        subset = lincs_base.loc[lincs_base["perturbation_type"] == perturbation_type].copy()
        meta = subset[
            [
                "cell_line",
                "target",
                "dose_val",
                "time_val",
                "canonical_query_uid",
                "delta_row_index",
                "delta_valid_bool",
            ]
        ].reset_index(drop=True)

        cohorts.append(
            CohortData(
                key=CohortKey(
                    scope="internal",
                    dataset_or_direction="LINCS",
                    perturbation_type=perturbation_type,
                    representation="Gene",
                ),
                meta=meta.copy(),
                vector_store=lincs_gene_store,
            )
        )
        cohorts.append(
            CohortData(
                key=CohortKey(
                    scope="internal",
                    dataset_or_direction="LINCS",
                    perturbation_type=perturbation_type,
                    representation="Pathway",
                ),
                meta=meta.copy(),
                vector_store=lincs_pathway_store,
            )
        )

    # scPerturb: gene + pathway.
    for perturbation_type in ["Chemical", "Genetic"]:
        src = sc_meta_map[perturbation_type]
        normalized = pd.DataFrame(
            {
                "cell_line": normalize_text(src["cell_std"]),
                "target": normalize_text(src["target_std"]),
                "dose_val": src["dose_val"],
                "time_val": src["time_val"],
                "canonical_query_uid": src["delta_row_idx"].astype(np.int64).astype(str),
                "delta_row_index": src["delta_row_idx"].astype(np.int64),
                "delta_valid_bool": np.ones(len(src), dtype=bool),
            }
        )

        cohorts.append(
            CohortData(
                key=CohortKey(
                    scope="internal",
                    dataset_or_direction="scPerturb",
                    perturbation_type=perturbation_type,
                    representation="Gene",
                ),
                meta=normalized.copy(),
                vector_store=sc_gene_map[perturbation_type],
            )
        )
        cohorts.append(
            CohortData(
                key=CohortKey(
                    scope="internal",
                    dataset_or_direction="scPerturb",
                    perturbation_type=perturbation_type,
                    representation="Pathway",
                ),
                meta=normalized.copy(),
                vector_store=sc_pathway_map[perturbation_type],
            )
        )

    # scPerturb FM cohorts.
    fm_cohorts, fm_paths, fm_policy_rows = load_fm_inputs(snapshot_root, sc_meta_map)
    cohorts.extend(fm_cohorts)

    input_paths.extend(lincs_paths)
    input_paths.extend(sc_paths)
    input_paths.extend(fm_paths)

    return (
        sorted(
            cohorts,
            key=lambda c: (
                c.key.scope,
                c.key.dataset_or_direction,
                c.key.perturbation_type,
                c.key.representation,
            ),
        ),
        sorted(set(input_paths)),
        {
            "pathway_policy": pathway_policy,
            "fm_policy_rows": fm_policy_rows,
        },
    )


def process_cohort(
    cohort: CohortData,
    rng: np.random.Generator,
    n_gallery_max_cap: Optional[int],
    query_batch_size: int,
    per_query_sink: PerQueryParquetSink,
    attrition_bucket: MutableMapping[Tuple[str, str, str, str, str], Dict[str, object]],
) -> Tuple[Dict[str, object], Dict[str, object], List[Dict[str, object]], Dict[str, int]]:
    key = cohort.key
    meta = cohort.meta.copy().reset_index(drop=True)

    n_total = int(len(meta))
    if n_total == 0:
        empty_summary = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
            "n_total": 0,
            "n_valid": 0,
            "n_excluded_missing_metric_or_mpos0": 0,
            "N_gallery_max": 0,
            "mean_mrr_raw": np.nan,
            "mean_expected_mrr_chance": np.nan,
            "mean_mrr_corrected": np.nan,
            "mean_hit1_raw": np.nan,
            "mean_expected_hit1_chance": np.nan,
            "mean_hit1_corrected": np.nan,
            "mean_hit5_raw": np.nan,
            "mean_expected_hit5_chance": np.nan,
            "mean_hit5_corrected": np.nan,
            "mean_hit10_raw": np.nan,
            "mean_expected_hit10_chance": np.nan,
            "mean_hit10_corrected": np.nan,
        }
        chance_row = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
            "delta_mrr": np.nan,
            "abs_delta_mrr": np.nan,
            "delta_hit1": np.nan,
            "abs_delta_hit1": np.nan,
            "delta_hit5": np.nan,
            "abs_delta_hit5": np.nan,
            "delta_hit10": np.nan,
            "abs_delta_hit10": np.nan,
        }
        return empty_summary, chance_row, [], {"loo_violation": 0}

    meta["delta_valid_bool"] = meta["delta_valid_bool"].fillna(False).astype(bool)
    mask_by_policy = meta["delta_valid_bool"].to_numpy()
    dropped_by_policy = int((~mask_by_policy).sum())
    record_attrition(
        attrition_bucket,
        key,
        "invalid_delta",
        n_dropped=dropped_by_policy,
        n_total_before=n_total,
        notes="Excluded by delta_valid_bool policy",
    )

    meta_work = meta.loc[mask_by_policy].copy().reset_index(drop=True)

    if meta_work.empty:
        summary = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
            "n_total": n_total,
            "n_valid": 0,
            "n_excluded_missing_metric_or_mpos0": n_total,
            "N_gallery_max": 0,
            "mean_mrr_raw": np.nan,
            "mean_expected_mrr_chance": np.nan,
            "mean_mrr_corrected": np.nan,
            "mean_hit1_raw": np.nan,
            "mean_expected_hit1_chance": np.nan,
            "mean_hit1_corrected": np.nan,
            "mean_hit5_raw": np.nan,
            "mean_expected_hit5_chance": np.nan,
            "mean_hit5_corrected": np.nan,
            "mean_hit10_raw": np.nan,
            "mean_expected_hit10_chance": np.nan,
            "mean_hit10_corrected": np.nan,
        }
        chance = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
            "delta_mrr": np.nan,
            "abs_delta_mrr": np.nan,
            "delta_hit1": np.nan,
            "abs_delta_hit1": np.nan,
            "delta_hit5": np.nan,
            "abs_delta_hit5": np.nan,
            "delta_hit10": np.nan,
            "abs_delta_hit10": np.nan,
        }
        return summary, chance, [], {"loo_violation": 0}

    meta_work["group_id"] = [
        make_group_id(cell_line=row.cell_line, target=row.target)
        for row in meta_work.itertuples(index=False)
    ]

    group_to_positions: Dict[str, np.ndarray] = {}
    for group_id, idx in meta_work.groupby("group_id", sort=False).groups.items():
        group_to_positions[group_id] = np.asarray(sorted(idx), dtype=np.int64)

    sorted_group_ids = sorted(group_to_positions.keys())

    # First pass: enforce finite vectors, collect group stats, compute split-half group metrics.
    group_positions_valid: Dict[str, np.ndarray] = {}
    group_sum: Dict[str, np.ndarray] = {}
    group_count: Dict[str, int] = {}
    group_metrics: List[Dict[str, object]] = []

    for group_id in sorted_group_ids:
        pos_all = group_to_positions[group_id]
        row_idx = meta_work.loc[pos_all, "delta_row_index"].to_numpy(dtype=np.int64)
        vectors = cohort.vector_store.load_rows(row_idx)
        finite_mask = np.isfinite(vectors).all(axis=1)

        dropped_nonfinite = int((~finite_mask).sum())
        record_attrition(
            attrition_bucket,
            key,
            "invalid_delta",
            n_dropped=dropped_nonfinite,
            n_total_before=int(len(pos_all)),
            notes="Excluded due to non-finite values",
        )

        if not finite_mask.any():
            continue

        pos_valid = pos_all[finite_mask]
        vec_valid = vectors[finite_mask]
        uid_valid = meta_work.loc[pos_valid, "canonical_query_uid"].astype(str).to_numpy()

        group_positions_valid[group_id] = pos_valid
        group_sum[group_id] = np.sum(vec_valid, axis=0, dtype=np.float64)
        group_count[group_id] = int(vec_valid.shape[0])

        n_group = int(vec_valid.shape[0])
        metric_row = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
            "group_id": group_id,
            "n_total": n_group,
            "n_A": np.nan,
            "n_B": np.nan,
            "n_A_sub": np.nan,
            "n_B_sub": np.nan,
            "cosine": np.nan,
            "pcc": np.nan,
            "edist": np.nan,
        }

        if n_group < 4:
            record_attrition(
                attrition_bucket,
                key,
                "split_half_requires_n>=4",
                n_dropped=1,
                n_total_before=n_group,
                notes="Group-level split-half requires n>=4",
            )
            group_metrics.append(metric_row)
            continue

        sort_idx = np.argsort(uid_valid, kind="mergesort")
        vec_sorted = vec_valid[sort_idx]
        uid_sorted = uid_valid[sort_idx]

        perm = rng.permutation(n_group)
        half = n_group // 2
        a_idx = perm[:half]
        b_idx = perm[half:]

        a_vec = vec_sorted[a_idx]
        b_vec = vec_sorted[b_idx]
        a_uid = uid_sorted[a_idx]
        b_uid = uid_sorted[b_idx]

        centroid_a = np.mean(a_vec, axis=0)
        centroid_b = np.mean(b_vec, axis=0)

        cosine = cosine_similarity(centroid_a, centroid_b)
        pcc = pearson_corr(centroid_a, centroid_b)

        if len(a_vec) < 2 or len(b_vec) < 2:
            record_attrition(
                attrition_bucket,
                key,
                "edist_requires_n>=2_per_split",
                n_dropped=1,
                n_total_before=n_group,
                notes="E-distance split requires at least 2 samples per side",
            )
            edist = np.nan
            n_a_sub = int(len(a_vec))
            n_b_sub = int(len(b_vec))
        else:
            a_sub, n_a_sub = deterministic_subsample_by_uid(a_vec, a_uid, EDIST_MAX_N, rng)
            b_sub, n_b_sub = deterministic_subsample_by_uid(b_vec, b_uid, EDIST_MAX_N, rng)
            if len(a_sub) < 2 or len(b_sub) < 2:
                record_attrition(
                    attrition_bucket,
                    key,
                    "edist_requires_n>=2_per_split",
                    n_dropped=1,
                    n_total_before=n_group,
                    notes="E-distance after subsample requires at least 2 samples per side",
                )
                edist = np.nan
            else:
                edist = energy_distance_biascorr(a_sub, b_sub)

        metric_row.update(
            {
                "n_A": int(len(a_vec)),
                "n_B": int(len(b_vec)),
                "n_A_sub": int(n_a_sub),
                "n_B_sub": int(n_b_sub),
                "cosine": float(cosine) if np.isfinite(cosine) else np.nan,
                "pcc": float(pcc) if np.isfinite(pcc) else np.nan,
                "edist": float(edist) if np.isfinite(edist) else np.nan,
            }
        )
        group_metrics.append(metric_row)

    valid_group_ids = sorted(group_sum.keys())

    if not valid_group_ids:
        summary = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
            "n_total": n_total,
            "n_valid": 0,
            "n_excluded_missing_metric_or_mpos0": n_total,
            "N_gallery_max": 0,
            "mean_mrr_raw": np.nan,
            "mean_expected_mrr_chance": np.nan,
            "mean_mrr_corrected": np.nan,
            "mean_hit1_raw": np.nan,
            "mean_expected_hit1_chance": np.nan,
            "mean_hit1_corrected": np.nan,
            "mean_hit5_raw": np.nan,
            "mean_expected_hit5_chance": np.nan,
            "mean_hit5_corrected": np.nan,
            "mean_hit10_raw": np.nan,
            "mean_expected_hit10_chance": np.nan,
            "mean_hit10_corrected": np.nan,
        }
        chance = {
            "scope": key.scope,
            "dataset_or_direction": key.dataset_or_direction,
            "perturbation_type": key.perturbation_type,
            "representation": key.representation,
            "delta_mrr": np.nan,
            "abs_delta_mrr": np.nan,
            "delta_hit1": np.nan,
            "abs_delta_hit1": np.nan,
            "delta_hit5": np.nan,
            "abs_delta_hit5": np.nan,
            "delta_hit10": np.nan,
            "abs_delta_hit10": np.nan,
        }
        return summary, chance, group_metrics, {"loo_violation": 0}

    centroid_matrix = np.stack(
        [group_sum[group_id] / float(group_count[group_id]) for group_id in valid_group_ids],
        axis=0,
    )
    centroid_norms = np.sum(centroid_matrix * centroid_matrix, axis=1)
    group_to_centroid_idx = {group_id: i for i, group_id in enumerate(valid_group_ids)}

    cap = n_gallery_max_cap
    if cap is None or cap <= 0 or cap >= len(valid_group_ids):
        base_gallery_ids = list(valid_group_ids)
    else:
        base_gallery_ids = list(valid_group_ids[:cap])

    base_gallery_idx = np.asarray([group_to_centroid_idx[g] for g in base_gallery_ids], dtype=np.int64)
    base_centroids = centroid_matrix[base_gallery_idx]
    base_centroid_norms = centroid_norms[base_gallery_idx]
    base_gallery_hash = hash_gallery_group_ids(base_gallery_ids)

    # Retrieval accumulators.
    n_valid_queries = 0
    max_n_gallery = 0
    sum_mrr_raw = 0.0
    sum_exp_mrr = 0.0
    sum_mrr_corr = 0.0
    sum_hit_raw: Dict[int, float] = {k: 0.0 for k in HIT_K_LIST}
    sum_hit_exp: Dict[int, float] = {k: 0.0 for k in HIT_K_LIST}
    sum_hit_corr: Dict[int, float] = {k: 0.0 for k in HIT_K_LIST}

    loo_violations = 0

    for group_id in valid_group_ids:
        count_g = int(group_count[group_id])
        pos_g = group_positions_valid[group_id]

        if count_g <= 1:
            record_attrition(
                attrition_bucket,
                key,
                "loo_requires_count>1",
                n_dropped=int(count_g),
                n_total_before=int(count_g),
                notes="Strict LOO requires group count > 1",
            )
            continue

        if group_id in base_gallery_ids:
            gallery_ids = base_gallery_ids
            gallery_centroids = base_centroids
            gallery_norms = base_centroid_norms
            gallery_hash = base_gallery_hash
            true_gallery_pos = gallery_ids.index(group_id)
        else:
            # Deterministic inclusion policy when capped.
            gallery_ids = list(base_gallery_ids[:-1]) + [group_id]
            gallery_idx = np.asarray([group_to_centroid_idx[g] for g in gallery_ids], dtype=np.int64)
            gallery_centroids = centroid_matrix[gallery_idx]
            gallery_norms = centroid_norms[gallery_idx]
            gallery_hash = hash_gallery_group_ids(gallery_ids)
            true_gallery_pos = len(gallery_ids) - 1

        n_gallery = int(len(gallery_ids))
        max_n_gallery = max(max_n_gallery, n_gallery)
        expected_mrr = expected_mrr_for_single_positive(n_gallery)
        expected_hit = {k: expected_hitk_for_single_positive(n_gallery, k) for k in HIT_K_LIST}

        # Query vectors for this group.
        row_idx = meta_work.loc[pos_g, "delta_row_index"].to_numpy(dtype=np.int64)
        q_matrix_all = cohort.vector_store.load_rows(row_idx)
        finite_q = np.isfinite(q_matrix_all).all(axis=1)

        dropped_missing = int((~finite_q).sum())
        record_attrition(
            attrition_bucket,
            key,
            "missing_metric",
            n_dropped=dropped_missing,
            n_total_before=int(len(pos_g)),
            notes="Query excluded due to non-finite vector during retrieval",
        )

        if not finite_q.any():
            continue

        q_matrix = q_matrix_all[finite_q]
        pos_valid = pos_g[finite_q]

        if count_g <= 1:
            loo_violations += int(len(pos_valid))
            continue

        group_sum_vec = group_sum[group_id]

        for start in range(0, q_matrix.shape[0], query_batch_size):
            end = min(start + query_batch_size, q_matrix.shape[0])
            q_batch = q_matrix[start:end]
            pos_batch = pos_valid[start:end]

            q_norm = np.sum(q_batch * q_batch, axis=1)
            scores = -(
                q_norm[:, None] + gallery_norms[None, :] - 2.0 * (q_batch @ gallery_centroids.T)
            )

            true_centroid = (group_sum_vec[None, :] - q_batch) / float(count_g - 1)
            true_norm = np.sum(true_centroid * true_centroid, axis=1)
            score_true = -(q_norm + true_norm - 2.0 * np.sum(q_batch * true_centroid, axis=1))

            scores[:, true_gallery_pos] = score_true
            rank_true = 1 + np.sum(scores > score_true[:, None], axis=1)

            mrr_raw = 1.0 / rank_true.astype(np.float64)
            hit_raw = {
                k: (rank_true <= k).astype(np.float64)
                for k in HIT_K_LIST
            }
            mrr_corr = mrr_raw - expected_mrr
            hit_corr = {k: hit_raw[k] - expected_hit[k] for k in HIT_K_LIST}

            n_batch = int(len(rank_true))
            n_valid_queries += n_batch
            sum_mrr_raw += float(np.sum(mrr_raw))
            sum_exp_mrr += float(expected_mrr * n_batch)
            sum_mrr_corr += float(np.sum(mrr_corr))
            for k in HIT_K_LIST:
                sum_hit_raw[k] += float(np.sum(hit_raw[k]))
                sum_hit_exp[k] += float(expected_hit[k] * n_batch)
                sum_hit_corr[k] += float(np.sum(hit_corr[k]))

            meta_batch = meta_work.loc[pos_batch, ["canonical_query_uid", "cell_line", "target"]].copy()
            frame = pd.DataFrame(
                {
                    "scope": key.scope,
                    "dataset_or_direction": key.dataset_or_direction,
                    "perturbation_type": key.perturbation_type,
                    "representation": key.representation,
                    "query_uid": meta_batch["canonical_query_uid"].astype(str).to_numpy(),
                    "cell_line": meta_batch["cell_line"].astype(str).to_numpy(),
                    "target_token": meta_batch["target"].astype(str).to_numpy(),
                    "N_gallery": np.full(n_batch, n_gallery, dtype=np.int64),
                    "N_gallery_max": np.full(n_batch, n_gallery, dtype=np.int64),
                    "gallery_group_ids_hash": np.full(n_batch, gallery_hash, dtype=object),
                    "m_pos": np.ones(n_batch, dtype=np.int64),
                    "rank_true": rank_true.astype(np.int64),
                    "mrr_raw": mrr_raw.astype(np.float64),
                    "hit1_raw": hit_raw[1].astype(np.float64),
                    "hit5_raw": hit_raw[5].astype(np.float64),
                    "hit10_raw": hit_raw[10].astype(np.float64),
                    "expected_mrr_chance": np.full(n_batch, expected_mrr, dtype=np.float64),
                    "expected_hit1_chance": np.full(n_batch, expected_hit[1], dtype=np.float64),
                    "expected_hit5_chance": np.full(n_batch, expected_hit[5], dtype=np.float64),
                    "expected_hit10_chance": np.full(n_batch, expected_hit[10], dtype=np.float64),
                    "mrr_corrected": mrr_corr.astype(np.float64),
                    "hit1_corrected": hit_corr[1].astype(np.float64),
                    "hit5_corrected": hit_corr[5].astype(np.float64),
                    "hit10_corrected": hit_corr[10].astype(np.float64),
                    "loo_policy": np.full(n_batch, "strict_recompute", dtype=object),
                    "cross_alignment_contract": np.full(n_batch, "NA", dtype=object),
                }
            )
            per_query_sink.write(canonicalize_per_query_frame(frame))

    n_excluded = int(n_total - n_valid_queries)

    if n_valid_queries > 0:
        mean_mrr_raw = sum_mrr_raw / n_valid_queries
        mean_expected_mrr = sum_exp_mrr / n_valid_queries
        mean_mrr_corr = sum_mrr_corr / n_valid_queries

        mean_hit_raw = {k: sum_hit_raw[k] / n_valid_queries for k in HIT_K_LIST}
        mean_hit_exp = {k: sum_hit_exp[k] / n_valid_queries for k in HIT_K_LIST}
        mean_hit_corr = {k: sum_hit_corr[k] / n_valid_queries for k in HIT_K_LIST}

        delta_mrr = mean_mrr_corr - (mean_mrr_raw - mean_expected_mrr)
        delta_hit = {
            k: mean_hit_corr[k] - (mean_hit_raw[k] - mean_hit_exp[k])
            for k in HIT_K_LIST
        }
    else:
        mean_mrr_raw = np.nan
        mean_expected_mrr = np.nan
        mean_mrr_corr = np.nan
        mean_hit_raw = {k: np.nan for k in HIT_K_LIST}
        mean_hit_exp = {k: np.nan for k in HIT_K_LIST}
        mean_hit_corr = {k: np.nan for k in HIT_K_LIST}
        delta_mrr = np.nan
        delta_hit = {k: np.nan for k in HIT_K_LIST}

    summary = {
        "scope": key.scope,
        "dataset_or_direction": key.dataset_or_direction,
        "perturbation_type": key.perturbation_type,
        "representation": key.representation,
        "n_total": int(n_total),
        "n_valid": int(n_valid_queries),
        "n_excluded_missing_metric_or_mpos0": int(n_excluded),
        "N_gallery_max": int(max_n_gallery),
        "mean_mrr_raw": float(mean_mrr_raw) if np.isfinite(mean_mrr_raw) else np.nan,
        "mean_expected_mrr_chance": float(mean_expected_mrr) if np.isfinite(mean_expected_mrr) else np.nan,
        "mean_mrr_corrected": float(mean_mrr_corr) if np.isfinite(mean_mrr_corr) else np.nan,
        "mean_hit1_raw": float(mean_hit_raw[1]) if np.isfinite(mean_hit_raw[1]) else np.nan,
        "mean_expected_hit1_chance": float(mean_hit_exp[1]) if np.isfinite(mean_hit_exp[1]) else np.nan,
        "mean_hit1_corrected": float(mean_hit_corr[1]) if np.isfinite(mean_hit_corr[1]) else np.nan,
        "mean_hit5_raw": float(mean_hit_raw[5]) if np.isfinite(mean_hit_raw[5]) else np.nan,
        "mean_expected_hit5_chance": float(mean_hit_exp[5]) if np.isfinite(mean_hit_exp[5]) else np.nan,
        "mean_hit5_corrected": float(mean_hit_corr[5]) if np.isfinite(mean_hit_corr[5]) else np.nan,
        "mean_hit10_raw": float(mean_hit_raw[10]) if np.isfinite(mean_hit_raw[10]) else np.nan,
        "mean_expected_hit10_chance": float(mean_hit_exp[10]) if np.isfinite(mean_hit_exp[10]) else np.nan,
        "mean_hit10_corrected": float(mean_hit_corr[10]) if np.isfinite(mean_hit_corr[10]) else np.nan,
    }

    chance = {
        "scope": key.scope,
        "dataset_or_direction": key.dataset_or_direction,
        "perturbation_type": key.perturbation_type,
        "representation": key.representation,
        "delta_mrr": float(delta_mrr) if np.isfinite(delta_mrr) else np.nan,
        "abs_delta_mrr": float(abs(delta_mrr)) if np.isfinite(delta_mrr) else np.nan,
        "delta_hit1": float(delta_hit[1]) if np.isfinite(delta_hit[1]) else np.nan,
        "abs_delta_hit1": float(abs(delta_hit[1])) if np.isfinite(delta_hit[1]) else np.nan,
        "delta_hit5": float(delta_hit[5]) if np.isfinite(delta_hit[5]) else np.nan,
        "abs_delta_hit5": float(abs(delta_hit[5])) if np.isfinite(delta_hit[5]) else np.nan,
        "delta_hit10": float(delta_hit[10]) if np.isfinite(delta_hit[10]) else np.nan,
        "abs_delta_hit10": float(abs(delta_hit[10])) if np.isfinite(delta_hit[10]) else np.nan,
    }

    counters = {"loo_violation": loo_violations}
    return summary, chance, group_metrics, counters


def build_leaderboard_rows(
    summary_rows: Sequence[Mapping[str, object]],
    group_metric_rows: Sequence[Mapping[str, object]],
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
                        "cross_alignment_contract": "NA",
                    }
                )

    group_df = pd.DataFrame(group_metric_rows)
    if group_df.empty:
        return leaderboard

    metric_map = {
        "mean_cosine_centroid": "cosine",
        "mean_pcc_centroid": "pcc",
        "mean_edist_biascorr": "edist",
    }

    group_keys = ["scope", "dataset_or_direction", "perturbation_type", "representation"]
    for key_values, sub in group_df.groupby(group_keys, sort=False):
        n_total = int(len(sub))
        key_dict = dict(zip(group_keys, key_values))

        for metric_name, metric_col in metric_map.items():
            finite = sub[metric_col].dropna()
            n_valid = int(len(finite))
            n_excluded = int(n_total - n_valid)
            value = float(finite.mean()) if n_valid > 0 else np.nan

            leaderboard.append(
                {
                    "scope": key_dict["scope"],
                    "dataset_or_direction": key_dict["dataset_or_direction"],
                    "perturbation_type": key_dict["perturbation_type"],
                    "representation": key_dict["representation"],
                    "metric_name": metric_name,
                    "metric_value": value,
                    "n_total": n_total,
                    "n_valid": n_valid,
                    "n_excluded": n_excluded,
                    "N_gallery_max": 0,
                    "cross_alignment_contract": "NA",
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
        print(
            f"[ERROR] Seed must be GLOBAL_SEED={GLOBAL_SEED}, got {seed}",
            file=sys.stderr,
        )
        return 3

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
    assertions.append(
        {
            "name": "task1_snapshot_isolation",
            "pass": True,
            "details": {
                "rules": [
                    "S1 reads exclusively from data/task1_snapshot_v1/",
                    "config.paths.task1_snapshot must equal data/task1_snapshot_v1",
                ],
                "task1_snapshot": str(task1_snapshot),
            },
            "counterexamples": [],
        }
    )

    try:
        cohorts, input_paths, policy_info = build_cohorts(task1_snapshot)
    except Exception as exc:
        assertions.append(
            {
                "name": "input_loading_and_contract_validation",
                "pass": False,
                "details": {
                    "rules": [
                        "All required snapshot files/columns/shape contracts must pass before compute",
                    ],
                    "error": str(exc),
                },
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        return 6

    pathway_policy = policy_info["pathway_policy"]
    fm_policy_rows = policy_info["fm_policy_rows"]

    assertions.append(
        {
            "name": "lincs_pathway_project_on_load_policy",
            "pass": str(pathway_policy.get("mode", "")).strip() == "project_on_load",
            "details": {
                "rules": [
                    "LINCS pathway must be derived dynamically: pathway = gene @ W",
                    "Pathway policy mode must equal project_on_load",
                ],
                "mode": pathway_policy.get("mode", "NA"),
                "W_path": pathway_policy.get("W_path", "NA"),
            },
            "counterexamples": []
            if str(pathway_policy.get("mode", "")).strip() == "project_on_load"
            else [{"mode": pathway_policy.get("mode", "NA")}],
        }
    )

    # Deterministic compute over all cohorts.
    per_query_path = stage_dir / "task1_retrieval_per_query.parquet"
    sink = PerQueryParquetSink(per_query_path, per_query_schema())

    attrition_bucket: Dict[Tuple[str, str, str, str, str], Dict[str, object]] = {}
    summary_rows: List[Dict[str, object]] = []
    chance_rows: List[Dict[str, object]] = []
    group_metric_rows: List[Dict[str, object]] = []

    loo_violation_total = 0

    for cohort in cohorts:
        summary, chance, group_metrics, counters = process_cohort(
            cohort=cohort,
            rng=rng,
            n_gallery_max_cap=args.n_gallery_max,
            query_batch_size=max(1, int(args.query_batch_size)),
            per_query_sink=sink,
            attrition_bucket=attrition_bucket,
        )
        summary_rows.append(summary)
        chance_rows.append(chance)
        group_metric_rows.extend(group_metrics)
        loo_violation_total += int(counters.get("loo_violation", 0))

    sink.close()

    # Build final tabular outputs.
    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["scope", "dataset_or_direction", "perturbation_type", "representation"],
        kind="mergesort",
    )

    chance_df = pd.DataFrame(chance_rows).sort_values(
        ["scope", "dataset_or_direction", "perturbation_type", "representation"],
        kind="mergesort",
    )

    attrition_df = pd.DataFrame(list(attrition_bucket.values()))
    if attrition_df.empty:
        attrition_df = pd.DataFrame(
            columns=[
                "scope",
                "dataset_or_direction",
                "perturbation_type",
                "representation",
                "reason",
                "n_dropped",
                "n_total_before",
                "notes",
            ]
        )
    else:
        attrition_df = attrition_df.sort_values(
            ["scope", "dataset_or_direction", "perturbation_type", "representation", "reason"],
            kind="mergesort",
        )

    leaderboard_df = pd.DataFrame(
        build_leaderboard_rows(summary_rows=summary_rows, group_metric_rows=group_metric_rows)
    )
    if leaderboard_df.empty:
        leaderboard_df = pd.DataFrame(
            columns=[
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
        )
    else:
        leaderboard_df = leaderboard_df.sort_values(
            [
                "scope",
                "dataset_or_direction",
                "perturbation_type",
                "representation",
                "metric_name",
            ],
            kind="mergesort",
        )

    # Persist canonical outputs.
    retrieval_summary_path = stage_dir / "task1_retrieval_summary.csv"
    chance_check_path = stage_dir / "task1_chance_identity_check.csv"
    leaderboard_path = stage_dir / "task1_leaderboard_long.csv"
    attrition_path = stage_dir / "task1_attrition.csv"

    write_csv(summary_df, retrieval_summary_path)
    write_csv(chance_df, chance_check_path)
    write_csv(leaderboard_df, leaderboard_path)
    write_csv(attrition_df, attrition_path)

    # Assertions
    denominator_failures = summary_df.loc[
        summary_df["n_total"] != (
            summary_df["n_valid"] + summary_df["n_excluded_missing_metric_or_mpos0"]
        )
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
            "counterexamples": denominator_failures.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    assertions.append(
        {
            "name": "loo_strict_recompute_policy",
            "pass": loo_violation_total == 0,
            "details": {
                "rules": [
                    "Internal retrieval must use strict LOO true centroid recomputation",
                    "Queries from groups with count <= 1 must be excluded",
                ],
                "loo_violations": int(loo_violation_total),
            },
            "counterexamples": [] if loo_violation_total == 0 else [{"loo_violations": loo_violation_total}],
        }
    )

    # Chance identity fail-fast.
    tol_cols = [
        "abs_delta_mrr",
        "abs_delta_hit1",
        "abs_delta_hit5",
        "abs_delta_hit10",
    ]
    chance_finite = chance_df.dropna(subset=tol_cols, how="any")
    if chance_finite.empty:
        chance_violation = pd.DataFrame(columns=chance_df.columns)
    else:
        chance_violation = chance_finite.loc[
            (chance_finite["abs_delta_mrr"] > CHANCE_IDENTITY_TOL)
            | (chance_finite["abs_delta_hit1"] > CHANCE_IDENTITY_TOL)
            | (chance_finite["abs_delta_hit5"] > CHANCE_IDENTITY_TOL)
            | (chance_finite["abs_delta_hit10"] > CHANCE_IDENTITY_TOL)
        ]

    chance_pass = chance_violation.empty
    assertions.append(
        {
            "name": "chance_identity_tolerance",
            "pass": bool(chance_pass),
            "details": {
                "rules": [
                    "For each cohort: |delta_mrr| <= 1e-12",
                    "For each cohort: |delta_hit1|, |delta_hit5|, |delta_hit10| <= 1e-12",
                ],
                "tolerance": CHANCE_IDENTITY_TOL,
                "n_checked": int(len(chance_finite)),
            },
            "counterexamples": chance_violation.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    # Output routing check.
    core_outputs = [
        per_query_path,
        retrieval_summary_path,
        chance_check_path,
        leaderboard_path,
        attrition_path,
    ]
    output_routing_pass = all(path.resolve().is_relative_to(stage_dir.resolve()) for path in core_outputs)
    assertions.append(
        {
            "name": "output_routing_isolated",
            "pass": bool(output_routing_pass),
            "details": {
                "rules": [
                    f"All task1_*.csv and *.parquet outputs must be inside {stage_dir}",
                ],
                "stage_dir": str(stage_dir),
            },
            "counterexamples": []
            if output_routing_pass
            else [{"bad_path": str(path)} for path in core_outputs if not path.resolve().is_relative_to(stage_dir.resolve())],
        }
    )

    # Input path routing (lexical path rooted at snapshot path).
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
                    "All loaded inputs must be under data/task1_snapshot_v1/ path namespace",
                ],
                "task1_snapshot": str(task1_snapshot),
                "n_inputs": len(input_paths),
            },
            "counterexamples": [{"bad_input": p} for p in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    # Optional FM policy evidence (soft assertion).
    fm_policy_missing = [r for r in fm_policy_rows if not r["policy_found"]]
    assertions.append(
        {
            "name": "fm_policy_presence_soft",
            "pass": len(fm_policy_missing) == 0,
            "details": {
                "rules": [
                    "FM model directories should provide delta_operator_policy.json when available",
                    "Missing policy does not block S1; finite fallback is enforced",
                ],
                "n_models": len(fm_policy_rows),
            },
            "counterexamples": fm_policy_missing[:MAX_COUNTEREXAMPLES],
        }
    )

    completed_at = utc_now_iso()

    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    manifest_path = stage_dir / "manifest.json"

    output_paths = [
        str(path.resolve())
        for path in core_outputs
        + [run_manifest_path, audit_assertions_path, manifest_path]
    ]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s1_task1_internal_metrics.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "edist_max_n": EDIST_MAX_N,
            "hit_k_list": list(HIT_K_LIST),
            "chance_identity_tol": CHANCE_IDENTITY_TOL,
            "task1_snapshot": str(task1_snapshot),
            "runs_dir": str(runs_dir),
            "n_gallery_max_cap": args.n_gallery_max,
            "query_batch_size": int(max(1, args.query_batch_size)),
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

    if loo_violation_total != 0:
        print("[ERROR] LOO strict_recompute assertion failed", file=sys.stderr)
        return 8

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
