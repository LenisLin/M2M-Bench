# SCRIPT_HEADER_CONTRACT
# Script: scripts/s3_build_task2_snapshot.py
# Purpose: Build deterministic Task2 K562 core snapshot artifacts (pairing + Gene/Pathway deltas) without FM extraction.
# Inputs:
#   - Raw K562 Data: config/config.yaml::paths.raw_k562_dir
#   - Task1 Pathway Assets: config/config.yaml::paths.task1_snapshot/pathway + lincs-gene-alignment.csv
# Outputs:
#   - Snapshot K562 base files:
#     <paths.task2_snapshot>/k562/{CRISPR_counts.pt,CRISPR_meta.csv,Drug_counts.pt,Drug_meta.csv,Common_Targets_K562.csv,shared_var_names.csv}
#   - pair_list.parquet: <paths.task2_snapshot>/k562/derived/pair_list.parquet
#   - delta_meta.csv: <paths.task2_snapshot>/k562/derived/delta_meta.csv
#   - gene_delta.npy: <paths.task2_snapshot>/k562/derived/gene_delta.npy
#   - pathway_delta.npy: <paths.task2_snapshot>/k562/derived/pathway_delta.npy
#   - task2_post_build_inventory.csv: runs/<run_id>/s3_build_task2_snapshot/
#   - task2_pairs_coverage.csv: runs/<run_id>/s3_build_task2_snapshot/
#   - k562_attrition.csv: runs/<run_id>/s3_build_task2_snapshot/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates/updates snapshot directory under config.paths.task2_snapshot
#   - Creates isolated run directory: runs/<run_id>/s3_build_task2_snapshot/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.raw_k562_dir
#   - config/config.yaml::paths.task1_snapshot
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s3_build_task2_snapshot.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Seed != 619 -> exit non-zero
#   - Missing required raw inputs/pathway assets -> exit non-zero
#   - No valid treated rows after pairing -> exit non-zero
#   - Pathway projection contract invalid -> exit non-zero
# Last Updated: 2026-03-05

"""
Inputs:
- Raw K562 files under config.paths.raw_k562_dir:
  - CRISPR_counts.pt
  - CRISPR_meta.csv
  - Drug_counts.pt
  - Drug_meta.csv
  - Common_Targets_K562.csv
  - shared_var_names.csv
- Task1 frozen pathway projection assets:
  - pathway/hallmark-w-2477x50.npy
  - pathway/lincs-pathway-policy.json
  - lincs/lincs-gene-alignment.csv

Outputs:
- data/task2_snapshot_v1/k562/{raw files...}
- data/task2_snapshot_v1/k562/derived/pair_list.parquet
- data/task2_snapshot_v1/k562/derived/delta_meta.csv
- data/task2_snapshot_v1/k562/derived/gene_delta.npy
- data/task2_snapshot_v1/k562/derived/pathway_delta.npy
- runs/<run_id>/s3_build_task2_snapshot/task2_post_build_inventory.csv
- runs/<run_id>/s3_build_task2_snapshot/task2_pairs_coverage.csv
- runs/<run_id>/s3_build_task2_snapshot/k562_attrition.csv
- runs/<run_id>/s3_build_task2_snapshot/run_manifest.json
- runs/<run_id>/s3_build_task2_snapshot/audit_assertions.json
- runs/<run_id>/s3_build_task2_snapshot/manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- K_CONTROLS_MAX = 50
- MIN_CONTROLS_REQUIRED = 1
- Cell-line is fixed to K562

Attrition rules:
- pool_size == 0 => soft-fail, drop treated row, record reason=no_controls_available
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import torch
import yaml

STAGE = "s3_build_task2_snapshot"
CONFIG_PATH = Path("config/config.yaml")

GLOBAL_SEED = 619
K_CONTROLS_MAX = 50
MIN_CONTROLS_REQUIRED = 1
MAX_COUNTEREXAMPLES = 5
PAIR_PARQUET_BATCH_ROWS = 200_000
PATHWAY_BATCH_ROWS = 512
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
EXPECTED_CELL_LINE = "K562"


@dataclass
class TreatedRecord:
    row_id: int
    treated_cell_id: str
    dataset_side: str
    perturbation_class: str
    treated_index: int
    control_indices: np.ndarray
    n_controls_used: int
    target_raw: str
    target_tokens: Tuple[str, ...]
    time_value: float
    dose_value: float
    specificity_tier: str


class PairListWriter:
    def __init__(self, output_path: Path, schema: pa.Schema, batch_rows: int) -> None:
        self.output_path = output_path
        self.schema = schema
        self.batch_rows = max(1, int(batch_rows))
        self.writer: Optional[pq.ParquetWriter] = None
        self.rows_written = 0
        self.buffer: Dict[str, List[Any]] = {field.name: [] for field in schema}

    def append_row(
        self,
        *,
        treated_cell_id: str,
        control_cell_id: str,
        control_rank: int,
        n_controls_used: int,
        dataset_side: str,
        perturbation_class: str,
        cell_line: str,
        target_raw: str,
        target_tokens: str,
        time_value: float,
        dose_value: float,
        specificity_tier: str,
        seed: int,
    ) -> None:
        self.buffer["treated_cell_id"].append(treated_cell_id)
        self.buffer["control_cell_id"].append(control_cell_id)
        self.buffer["control_rank"].append(int(control_rank))
        self.buffer["n_controls_used"].append(int(n_controls_used))
        self.buffer["dataset_side"].append(dataset_side)
        self.buffer["perturbation_class"].append(perturbation_class)
        self.buffer["cell_line"].append(cell_line)
        self.buffer["target_raw"].append(target_raw)
        self.buffer["target_tokens"].append(target_tokens)
        self.buffer["time"].append(float(time_value))
        self.buffer["dose_value"].append(float(dose_value))
        self.buffer["specificity_tier"].append(specificity_tier)
        self.buffer["seed"].append(int(seed))

        if len(self.buffer["treated_cell_id"]) >= self.batch_rows:
            self.flush()

    def flush(self) -> None:
        n = len(self.buffer["treated_cell_id"])
        if n == 0:
            return

        table = pa.Table.from_pydict(self.buffer, schema=self.schema)
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)
        self.writer.write_table(table)
        self.rows_written += n
        self.buffer = {field.name: [] for field in self.schema}

    def close(self) -> None:
        self.flush()
        if self.writer is None:
            empty = pa.Table.from_pydict(
                {field.name: pa.array([], type=field.type) for field in self.schema},
                schema=self.schema,
            )
            pq.write_table(empty, self.output_path)
        else:
            self.writer.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S3 Task2 K562 core snapshot builder")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
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


def normalize_scalar(value: object, default: str = "NA") -> str:
    if pd.isna(value):
        return default
    text = str(value).strip()
    return text if text else default


def to_float_or_nan(value: object) -> float:
    if pd.isna(value):
        return float("nan")
    return float(value)


def tokenize_target(raw_target: object) -> Tuple[str, ...]:
    raw = normalize_scalar(raw_target, default="NA")
    tokens: List[str] = []
    for chunk in raw.split(";"):
        token = chunk.strip()
        if token:
            tokens.append(token)
    if not tokens:
        return ("NA",)
    return tuple(sorted(set(tokens)))


def join_tokens(tokens: Sequence[str]) -> str:
    return ";".join(tokens)


def init_global_rng(seed: int) -> np.random.Generator:
    random.seed(seed)
    np.random.seed(seed)
    return np.random.default_rng(seed)


def pair_list_schema() -> pa.Schema:
    return pa.schema(
        [
            ("treated_cell_id", pa.string()),
            ("control_cell_id", pa.string()),
            ("control_rank", pa.int64()),
            ("n_controls_used", pa.int64()),
            ("dataset_side", pa.string()),
            ("perturbation_class", pa.string()),
            ("cell_line", pa.string()),
            ("target_raw", pa.string()),
            ("target_tokens", pa.string()),
            ("time", pa.float64()),
            ("dose_value", pa.float64()),
            ("specificity_tier", pa.string()),
            ("seed", pa.int64()),
        ]
    )


def link_or_copy(src: Path, dst: Path) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.is_file() and dst.stat().st_size == src.stat().st_size:
            return "reuse"
        dst.unlink()

    try:
        os.link(src, dst)
        return "hardlink"
    except OSError:
        shutil.copy2(src, dst)
        return "copy"


def select_controls(
    pool_indices: np.ndarray,
    control_ids: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    pool_size = int(pool_indices.size)
    if pool_size == 0:
        return np.empty((0,), dtype=np.int64)

    if pool_size <= K_CONTROLS_MAX:
        return pool_indices.copy()

    sampled_positions = rng.choice(pool_size, size=K_CONTROLS_MAX, replace=False)
    sampled = pool_indices[np.asarray(sampled_positions, dtype=np.int64)]
    order = np.argsort(control_ids[sampled], kind="mergesort")
    return sampled[order]


def append_attrition(
    rows: List[Dict[str, object]],
    *,
    dataset_side: str,
    perturbation_class: str,
    treated_cell_id: str,
    target_raw: str,
    time_value: float,
    dose_value: float,
    pool_size: int,
    reason: str,
) -> None:
    rows.append(
        {
            "dataset_side": dataset_side,
            "perturbation_class": perturbation_class,
            "treated_cell_id": treated_cell_id,
            "target_raw": target_raw,
            "time": time_value,
            "dose_value": dose_value,
            "pool_size": int(pool_size),
            "reason": reason,
        }
    )


def build_crispr_pairings(
    crispr_meta: pd.DataFrame,
    rng: np.random.Generator,
    pair_writer: PairListWriter,
    attrition_rows: List[Dict[str, object]],
    start_row_id: int,
) -> Tuple[List[TreatedRecord], int, Dict[str, int]]:
    id_col = "cell_barcode"
    cell_ids = crispr_meta[id_col].map(normalize_scalar).to_numpy(dtype=object)
    benchmark = crispr_meta["benchmark_group"].map(normalize_scalar)
    treated_indices = np.flatnonzero(benchmark.eq("Test").to_numpy())
    control_indices = np.flatnonzero(benchmark.eq("Control").to_numpy())

    if control_indices.size > 0:
        control_order = np.argsort(cell_ids[control_indices], kind="mergesort")
        control_pool = control_indices[control_order].astype(np.int64, copy=False)
    else:
        control_pool = np.empty((0,), dtype=np.int64)

    treated_order = np.argsort(cell_ids[treated_indices], kind="mergesort")
    treated_sorted = treated_indices[treated_order].astype(np.int64, copy=False)

    records: List[TreatedRecord] = []
    row_id = int(start_row_id)
    n_pairs_written = 0

    for treated_idx in treated_sorted.tolist():
        row = crispr_meta.iloc[int(treated_idx)]
        treated_cell_id = normalize_scalar(row[id_col])
        target_raw = normalize_scalar(row["clean_target_mapped"])
        target_tokens = tokenize_target(target_raw)
        specificity_tier = normalize_scalar(row["specificity_tier"])
        time_value = float("nan")
        dose_value = float("nan")

        pool_size = int(control_pool.size)
        if pool_size < MIN_CONTROLS_REQUIRED:
            append_attrition(
                attrition_rows,
                dataset_side="CRISPR",
                perturbation_class="Genetic",
                treated_cell_id=treated_cell_id,
                target_raw=target_raw,
                time_value=time_value,
                dose_value=dose_value,
                pool_size=pool_size,
                reason="no_controls_available",
            )
            continue

        selected = select_controls(control_pool, cell_ids, rng)
        n_controls_used = int(selected.size)

        records.append(
            TreatedRecord(
                row_id=row_id,
                treated_cell_id=treated_cell_id,
                dataset_side="CRISPR",
                perturbation_class="Genetic",
                treated_index=int(treated_idx),
                control_indices=selected.astype(np.int64, copy=False),
                n_controls_used=n_controls_used,
                target_raw=target_raw,
                target_tokens=target_tokens,
                time_value=time_value,
                dose_value=dose_value,
                specificity_tier=specificity_tier,
            )
        )

        token_text = join_tokens(target_tokens)
        for rank, control_idx in enumerate(selected.tolist(), start=1):
            pair_writer.append_row(
                treated_cell_id=treated_cell_id,
                control_cell_id=str(cell_ids[int(control_idx)]),
                control_rank=rank,
                n_controls_used=n_controls_used,
                dataset_side="CRISPR",
                perturbation_class="Genetic",
                cell_line=EXPECTED_CELL_LINE,
                target_raw=target_raw,
                target_tokens=token_text,
                time_value=time_value,
                dose_value=dose_value,
                specificity_tier=specificity_tier,
                seed=GLOBAL_SEED,
            )
            n_pairs_written += 1

        row_id += 1

    summary = {
        "n_treated_total": int(treated_sorted.size),
        "n_treated_kept": int(len(records)),
        "n_treated_dropped": int(treated_sorted.size - len(records)),
        "n_pairs_written": int(n_pairs_written),
    }
    return records, row_id, summary


def build_drug_pairings(
    drug_meta: pd.DataFrame,
    rng: np.random.Generator,
    pair_writer: PairListWriter,
    attrition_rows: List[Dict[str, object]],
    start_row_id: int,
) -> Tuple[List[TreatedRecord], int, Dict[str, int]]:
    id_col = "Unnamed: 0" if "Unnamed: 0" in drug_meta.columns else "cell_id"
    if id_col not in drug_meta.columns:
        raise ValueError("Drug_meta.csv must include either 'Unnamed: 0' or 'cell_id'")

    cell_ids = drug_meta[id_col].map(normalize_scalar).to_numpy(dtype=object)
    benchmark = drug_meta["benchmark_group"].map(normalize_scalar)
    times = pd.to_numeric(drug_meta["time"], errors="coerce")

    treated_indices = np.flatnonzero(benchmark.eq("Test").to_numpy())
    control_indices = np.flatnonzero(benchmark.eq("Control").to_numpy())

    controls_by_time: Dict[float, np.ndarray] = {}
    if control_indices.size > 0:
        for time_value in sorted(times.iloc[control_indices].dropna().astype(float).unique().tolist()):
            subset = control_indices[np.where(times.iloc[control_indices].to_numpy() == time_value)[0]]
            order = np.argsort(cell_ids[subset], kind="mergesort")
            controls_by_time[float(time_value)] = subset[order].astype(np.int64, copy=False)

    treated_order = np.argsort(cell_ids[treated_indices], kind="mergesort")
    treated_sorted = treated_indices[treated_order].astype(np.int64, copy=False)

    records: List[TreatedRecord] = []
    row_id = int(start_row_id)
    n_pairs_written = 0

    for treated_idx in treated_sorted.tolist():
        row = drug_meta.iloc[int(treated_idx)]
        treated_cell_id = normalize_scalar(row[id_col])
        target_raw = normalize_scalar(row["clean_target_mapped"])
        target_tokens = tokenize_target(target_raw)
        specificity_tier = normalize_scalar(row["specificity_tier"])
        time_value = to_float_or_nan(row["time"])
        dose_value = to_float_or_nan(row["dose_value"])

        pool = controls_by_time.get(float(time_value), np.empty((0,), dtype=np.int64))
        pool_size = int(pool.size)

        if pool_size < MIN_CONTROLS_REQUIRED:
            append_attrition(
                attrition_rows,
                dataset_side="DRUG",
                perturbation_class="Chemical",
                treated_cell_id=treated_cell_id,
                target_raw=target_raw,
                time_value=time_value,
                dose_value=dose_value,
                pool_size=pool_size,
                reason="no_controls_available",
            )
            continue

        selected = select_controls(pool, cell_ids, rng)
        n_controls_used = int(selected.size)

        records.append(
            TreatedRecord(
                row_id=row_id,
                treated_cell_id=treated_cell_id,
                dataset_side="DRUG",
                perturbation_class="Chemical",
                treated_index=int(treated_idx),
                control_indices=selected.astype(np.int64, copy=False),
                n_controls_used=n_controls_used,
                target_raw=target_raw,
                target_tokens=target_tokens,
                time_value=time_value,
                dose_value=dose_value,
                specificity_tier=specificity_tier,
            )
        )

        token_text = join_tokens(target_tokens)
        for rank, control_idx in enumerate(selected.tolist(), start=1):
            pair_writer.append_row(
                treated_cell_id=treated_cell_id,
                control_cell_id=str(cell_ids[int(control_idx)]),
                control_rank=rank,
                n_controls_used=n_controls_used,
                dataset_side="DRUG",
                perturbation_class="Chemical",
                cell_line=EXPECTED_CELL_LINE,
                target_raw=target_raw,
                target_tokens=token_text,
                time_value=time_value,
                dose_value=dose_value,
                specificity_tier=specificity_tier,
                seed=GLOBAL_SEED,
            )
            n_pairs_written += 1

        row_id += 1

    summary = {
        "n_treated_total": int(treated_sorted.size),
        "n_treated_kept": int(len(records)),
        "n_treated_dropped": int(treated_sorted.size - len(records)),
        "n_pairs_written": int(n_pairs_written),
    }
    return records, row_id, summary


def check_pairing_contract(records: Sequence[TreatedRecord]) -> List[Dict[str, object]]:
    violations: List[Dict[str, object]] = []
    for record in records:
        control_indices = record.control_indices
        unique_controls = np.unique(control_indices)
        if record.n_controls_used != int(control_indices.size):
            violations.append(
                {
                    "row_id": record.row_id,
                    "treated_cell_id": record.treated_cell_id,
                    "violation": "n_controls_used_mismatch",
                    "n_controls_used": record.n_controls_used,
                    "actual": int(control_indices.size),
                }
            )
        if record.n_controls_used > K_CONTROLS_MAX:
            violations.append(
                {
                    "row_id": record.row_id,
                    "treated_cell_id": record.treated_cell_id,
                    "violation": "n_controls_used_gt_k",
                    "n_controls_used": record.n_controls_used,
                    "k_controls_max": K_CONTROLS_MAX,
                }
            )
        if unique_controls.size != control_indices.size:
            violations.append(
                {
                    "row_id": record.row_id,
                    "treated_cell_id": record.treated_cell_id,
                    "violation": "controls_not_unique",
                    "n_controls_used": int(control_indices.size),
                    "n_unique_controls": int(unique_controls.size),
                }
            )
        if len(violations) >= MAX_COUNTEREXAMPLES:
            break
    return violations


def compute_gene_delta_and_meta(
    records: Sequence[TreatedRecord],
    crispr_counts: np.ndarray,
    drug_counts: np.ndarray,
    gene_delta_path: Path,
) -> pd.DataFrame:
    if len(records) == 0:
        raise ValueError("No treated records available after pairing.")

    if crispr_counts.shape[1] != drug_counts.shape[1]:
        raise ValueError("CRISPR and Drug counts must have same number of genes.")

    n_rows = len(records)
    n_genes = int(crispr_counts.shape[1])
    gene_out = np.lib.format.open_memmap(
        gene_delta_path,
        mode="w+",
        dtype=np.float32,
        shape=(n_rows, n_genes),
    )

    delta_meta_rows: List[Dict[str, object]] = []
    records_sorted = sorted(records, key=lambda r: r.row_id)

    for record in records_sorted:
        source = crispr_counts if record.dataset_side == "CRISPR" else drug_counts
        ctrl_matrix = source[record.control_indices]
        if ctrl_matrix.ndim != 2 or ctrl_matrix.shape[0] < MIN_CONTROLS_REQUIRED:
            raise ValueError(
                f"Invalid control matrix for row_id={record.row_id}, "
                f"treated_cell_id={record.treated_cell_id}"
            )
        ctrl_mean = np.mean(ctrl_matrix, axis=0, dtype=np.float32)
        delta = source[record.treated_index] - ctrl_mean
        gene_out[record.row_id] = delta.astype(np.float32, copy=False)

        delta_meta_rows.append(
            {
                "row_id": int(record.row_id),
                "treated_cell_id": record.treated_cell_id,
                "perturbation_class": record.perturbation_class,
                "cell_line": EXPECTED_CELL_LINE,
                "target_raw": record.target_raw,
                "time": float(record.time_value),
                "dose_value": float(record.dose_value),
                "specificity_tier": record.specificity_tier,
                "n_controls_used": int(record.n_controls_used),
                "dataset_side": record.dataset_side,
                "target_tokens": join_tokens(record.target_tokens),
                "seed": GLOBAL_SEED,
            }
        )

    gene_out.flush()
    del gene_out

    delta_meta = pd.DataFrame(delta_meta_rows).sort_values("row_id", kind="mergesort")
    return delta_meta.reset_index(drop=True)


def load_pathway_contract(
    task1_snapshot: Path,
    shared_var_names: pd.Series,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object], List[Path]]:
    w_path = task1_snapshot / "pathway/hallmark-w-2477x50.npy"
    policy_path = task1_snapshot / "pathway/lincs-pathway-policy.json"
    alignment_path = task1_snapshot / "lincs/lincs-gene-alignment.csv"

    for path in [w_path, policy_path, alignment_path]:
        if not path.is_file():
            raise FileNotFoundError(f"Missing required pathway asset: {path}")

    with policy_path.open("r", encoding="utf-8") as handle:
        policy = json.load(handle)

    mode = str(policy.get("mode", "")).strip()
    if mode != "project_on_load":
        raise ValueError(f"Pathway policy mode must be project_on_load, got {mode!r}")

    w = np.load(w_path).astype(np.float32, copy=False)
    if w.shape != (2477, 50):
        raise ValueError(f"hallmark-w-2477x50.npy must have shape (2477,50), got {w.shape}")

    policy_sha = str(policy.get("W_sha256", "")).strip()
    actual_sha = compute_sha256(w_path)
    if policy_sha and policy_sha != actual_sha:
        raise ValueError(f"W sha mismatch: policy={policy_sha}, actual={actual_sha}")

    alignment = pd.read_csv(alignment_path)
    ensure_required_columns(alignment, ["gene_symbol", "local_idx"], "lincs-gene-alignment.csv")
    alignment["local_idx"] = pd.to_numeric(alignment["local_idx"], errors="raise").astype(np.int64)

    if alignment["local_idx"].duplicated().any():
        dup = alignment.loc[alignment["local_idx"].duplicated(), "local_idx"].head(MAX_COUNTEREXAMPLES)
        raise ValueError(f"Alignment local_idx has duplicates, examples={dup.tolist()}")

    if alignment["local_idx"].min() < 0 or alignment["local_idx"].max() >= 2477:
        raise ValueError("Alignment local_idx must be in [0, 2476]")

    shared = shared_var_names.fillna("").astype(str).str.strip().str.upper()
    shared_lookup = {gene: idx for idx, gene in enumerate(shared.tolist())}

    index_map = np.full((2477,), -1, dtype=np.int64)
    alignment_sorted = alignment.sort_values("local_idx", kind="mergesort")
    for _, row in alignment_sorted.iterrows():
        local_idx = int(row["local_idx"])
        symbol = str(row["gene_symbol"]).strip().upper()
        if symbol in shared_lookup:
            index_map[local_idx] = int(shared_lookup[symbol])

    n_mapped = int(np.sum(index_map >= 0))
    if n_mapped < 1:
        raise ValueError("No pathway-alignment genes mapped to K562 shared_var_names.")

    details = {
        "mode": mode,
        "w_path": str(w_path),
        "alignment_path": str(alignment_path),
        "w_sha256": actual_sha,
        "n_w_rows": int(w.shape[0]),
        "n_w_cols": int(w.shape[1]),
        "n_mapped_genes": int(n_mapped),
        "n_unmapped_genes": int(w.shape[0] - n_mapped),
    }
    return w, index_map, details, [w_path, policy_path, alignment_path]


def project_pathway_delta(
    gene_delta_path: Path,
    pathway_delta_path: Path,
    w: np.ndarray,
    index_map: np.ndarray,
) -> Dict[str, int]:
    gene_mmap = np.load(gene_delta_path, mmap_mode="r")
    if gene_mmap.ndim != 2:
        raise ValueError(f"gene_delta.npy must be 2D, got shape={gene_mmap.shape}")

    n_rows = int(gene_mmap.shape[0])
    pathway_dim = int(w.shape[1])

    out = np.lib.format.open_memmap(
        pathway_delta_path,
        mode="w+",
        dtype=np.float32,
        shape=(n_rows, pathway_dim),
    )

    valid_local = np.where(index_map >= 0)[0].astype(np.int64, copy=False)
    valid_source = index_map[valid_local].astype(np.int64, copy=False)

    for start in range(0, n_rows, PATHWAY_BATCH_ROWS):
        end = min(start + PATHWAY_BATCH_ROWS, n_rows)
        block = np.asarray(gene_mmap[start:end], dtype=np.float32)
        projected_input = np.zeros((block.shape[0], w.shape[0]), dtype=np.float32)
        if valid_local.size > 0:
            projected_input[:, valid_local] = block[:, valid_source]
        out[start:end] = projected_input @ w

    out.flush()
    del out

    return {
        "n_rows": n_rows,
        "pathway_dim": pathway_dim,
        "n_mapped_genes": int(valid_local.size),
    }


def build_task2_pairs_coverage(
    records: Sequence[TreatedRecord],
    common_targets: Sequence[str],
) -> pd.DataFrame:
    chem_counts: MutableMapping[str, int] = {}
    gen_counts: MutableMapping[str, int] = {}

    for record in records:
        tokens = [t for t in record.target_tokens if t != "NA"]
        if record.perturbation_class == "Chemical":
            for token in tokens:
                chem_counts[token] = chem_counts.get(token, 0) + 1
        else:
            if tokens:
                token = tokens[0]
                gen_counts[token] = gen_counts.get(token, 0) + 1

    all_tokens = sorted(set(common_targets) | set(chem_counts.keys()) | set(gen_counts.keys()))
    rows: List[Dict[str, object]] = []
    for token in all_tokens:
        n_chem = int(chem_counts.get(token, 0))
        n_gen = int(gen_counts.get(token, 0))
        rows.append(
            {
                "dataset": "scPerturb",
                "cell_line": EXPECTED_CELL_LINE,
                "target_token": token,
                "n_chem_instances": n_chem,
                "n_gen_instances": n_gen,
                "is_eligible_bool": bool(n_chem > 0 and n_gen > 0),
                "source_tag": "K562 common targets",
            }
        )

    return pd.DataFrame(rows).sort_values("target_token", kind="mergesort").reset_index(drop=True)


def build_post_build_inventory(
    crispr_summary: Mapping[str, int],
    drug_summary: Mapping[str, int],
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for representation in ["Gene", "Pathway"]:
        for perturbation_type, summary in [
            ("Genetic", crispr_summary),
            ("Chemical", drug_summary),
        ]:
            drop = int(summary["n_treated_dropped"])
            rows.append(
                {
                    "dataset": "scPerturb",
                    "perturbation_type": perturbation_type,
                    "cell_line": EXPECTED_CELL_LINE,
                    "representation": representation,
                    "stage": "post_build",
                    "n_total_instances": int(summary["n_treated_total"]),
                    "n_valid_instances": int(summary["n_treated_kept"]),
                    "n_dropped": drop,
                    "drop_reason_breakdown": json.dumps(
                        {"no_controls_available": drop},
                        ensure_ascii=True,
                        sort_keys=True,
                    ),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["perturbation_type", "representation"], kind="mergesort"
    ).reset_index(drop=True)


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

    raw_k562_dir = resolve_config_path(project_root, str(config["paths"]["raw_k562_dir"]))
    task1_snapshot = resolve_config_path(project_root, str(config["paths"]["task1_snapshot"]))
    task2_snapshot = resolve_config_path(project_root, str(config["paths"]["task2_snapshot"]))
    runs_dir = resolve_config_path(project_root, str(config["paths"]["runs_dir"]))

    expected_task2_snapshot = (project_root / EXPECTED_TASK2_SNAPSHOT).resolve()
    if task2_snapshot != expected_task2_snapshot:
        print(
            "[ERROR] Data isolation violation: config.paths.task2_snapshot must resolve "
            f"to {expected_task2_snapshot}, got {task2_snapshot}",
            file=sys.stderr,
        )
        return 4

    stage_dir = runs_dir / args.run_id / STAGE
    stage_dir.mkdir(parents=True, exist_ok=True)

    started_at = utc_now_iso()
    rng = init_global_rng(GLOBAL_SEED)

    assertions: List[Dict[str, object]] = []
    input_paths: List[Path] = []
    snapshot_outputs: List[Path] = []
    stage_outputs: List[Path] = []

    assertions.append(
        {
            "name": "seed_locked_global_619",
            "pass": True,
            "details": {
                "rules": ["GLOBAL_SEED must be fixed to 619 for deterministic control sampling."],
                "seed": GLOBAL_SEED,
            },
            "counterexamples": [],
        }
    )

    if not raw_k562_dir.is_dir():
        assertions.append(
            {
                "name": "raw_k562_directory_exists",
                "pass": False,
                "details": {
                    "rules": ["config.paths.raw_k562_dir must exist and be a directory."],
                    "raw_k562_dir": str(raw_k562_dir),
                },
                "counterexamples": [{"missing_path": str(raw_k562_dir)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Missing raw_k562_dir: {raw_k562_dir}", file=sys.stderr)
        return 5

    required_raw_files = [
        "CRISPR_counts.pt",
        "CRISPR_meta.csv",
        "Drug_counts.pt",
        "Drug_meta.csv",
        "Common_Targets_K562.csv",
        "shared_var_names.csv",
    ]
    missing_raw = [name for name in required_raw_files if not (raw_k562_dir / name).is_file()]
    if missing_raw:
        assertions.append(
            {
                "name": "raw_k562_required_files_present",
                "pass": False,
                "details": {
                    "rules": ["All required raw K562 files must exist."],
                    "required_files": required_raw_files,
                },
                "counterexamples": [{"missing_file": name} for name in missing_raw[:MAX_COUNTEREXAMPLES]],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Missing required raw inputs: {missing_raw}", file=sys.stderr)
        return 6

    assertions.append(
        {
            "name": "raw_k562_required_files_present",
            "pass": True,
            "details": {
                "rules": ["All required raw K562 files must exist."],
                "required_files": required_raw_files,
                "raw_k562_dir": str(raw_k562_dir),
            },
            "counterexamples": [],
        }
    )

    raw_paths = {name: (raw_k562_dir / name).resolve() for name in required_raw_files}
    input_paths.extend(raw_paths.values())

    try:
        crispr_meta = pd.read_csv(raw_paths["CRISPR_meta.csv"])
        drug_meta = pd.read_csv(raw_paths["Drug_meta.csv"])
        common_targets_df = pd.read_csv(raw_paths["Common_Targets_K562.csv"])
        shared_var_names = pd.read_csv(raw_paths["shared_var_names.csv"])

        crispr_counts_t = torch.load(raw_paths["CRISPR_counts.pt"], map_location="cpu")
        drug_counts_t = torch.load(raw_paths["Drug_counts.pt"], map_location="cpu")
    except Exception as exc:
        assertions.append(
            {
                "name": "raw_input_loading",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed loading raw inputs: {exc}", file=sys.stderr)
        return 7

    required_crispr_cols = [
        "cell_barcode",
        "benchmark_group",
        "clean_target_mapped",
        "specificity_tier",
    ]
    required_drug_cols = [
        "benchmark_group",
        "clean_target_mapped",
        "specificity_tier",
        "time",
        "dose_value",
    ]
    ensure_required_columns(crispr_meta, required_crispr_cols, "CRISPR_meta.csv")
    ensure_required_columns(drug_meta, required_drug_cols, "Drug_meta.csv")
    if "Unnamed: 0" not in drug_meta.columns and "cell_id" not in drug_meta.columns:
        raise ValueError("Drug_meta.csv missing cell-id column: need 'Unnamed: 0' or 'cell_id'")

    if "Common_Target" not in common_targets_df.columns:
        raise ValueError("Common_Targets_K562.csv must include column 'Common_Target'")
    if "gene_symbol" not in shared_var_names.columns:
        raise ValueError("shared_var_names.csv must include column 'gene_symbol'")

    if not torch.is_tensor(crispr_counts_t) or not torch.is_tensor(drug_counts_t):
        raise ValueError("CRISPR_counts.pt and Drug_counts.pt must be torch.Tensor")

    crispr_counts = crispr_counts_t.detach().cpu().numpy().astype(np.float32, copy=False)
    drug_counts = drug_counts_t.detach().cpu().numpy().astype(np.float32, copy=False)

    if crispr_counts.shape[0] != len(crispr_meta):
        raise ValueError(
            f"CRISPR row mismatch: counts={crispr_counts.shape[0]}, meta={len(crispr_meta)}"
        )
    if drug_counts.shape[0] != len(drug_meta):
        raise ValueError(f"Drug row mismatch: counts={drug_counts.shape[0]}, meta={len(drug_meta)}")
    if crispr_counts.shape[1] != drug_counts.shape[1]:
        raise ValueError("CRISPR and Drug counts must share gene dimension.")
    if int(shared_var_names.shape[0]) != int(crispr_counts.shape[1]):
        raise ValueError(
            "shared_var_names length must equal count gene dimension, got "
            f"{shared_var_names.shape[0]} vs {crispr_counts.shape[1]}"
        )

    assertions.append(
        {
            "name": "k562_raw_alignment",
            "pass": True,
            "details": {
                "rules": [
                    "CRISPR_counts rows == CRISPR_meta rows",
                    "Drug_counts rows == Drug_meta rows",
                    "shared_var_names rows == gene dimension",
                ],
                "crispr_rows": int(crispr_counts.shape[0]),
                "drug_rows": int(drug_counts.shape[0]),
                "gene_dim": int(crispr_counts.shape[1]),
            },
            "counterexamples": [],
        }
    )

    # Materialize raw K562 files into task2 snapshot.
    snapshot_k562_dir = task2_snapshot / "k562"
    snapshot_derived_dir = snapshot_k562_dir / "derived"
    snapshot_derived_dir.mkdir(parents=True, exist_ok=True)

    materialization_rows: List[Dict[str, object]] = []
    for name in required_raw_files:
        src = raw_paths[name]
        dst = snapshot_k562_dir / name
        mode = link_or_copy(src, dst)
        snapshot_outputs.append(dst.resolve())
        materialization_rows.append({"file": name, "mode": mode, "dst": str(dst.resolve())})

    pair_list_path = (snapshot_derived_dir / "pair_list.parquet").resolve()
    delta_meta_path = (snapshot_derived_dir / "delta_meta.csv").resolve()
    gene_delta_path = (snapshot_derived_dir / "gene_delta.npy").resolve()
    pathway_delta_path = (snapshot_derived_dir / "pathway_delta.npy").resolve()

    pair_writer = PairListWriter(pair_list_path, pair_list_schema(), PAIR_PARQUET_BATCH_ROWS)
    attrition_rows: List[Dict[str, object]] = []
    records: List[TreatedRecord] = []

    crispr_records, next_row_id, crispr_summary = build_crispr_pairings(
        crispr_meta=crispr_meta,
        rng=rng,
        pair_writer=pair_writer,
        attrition_rows=attrition_rows,
        start_row_id=0,
    )
    drug_records, next_row_id, drug_summary = build_drug_pairings(
        drug_meta=drug_meta,
        rng=rng,
        pair_writer=pair_writer,
        attrition_rows=attrition_rows,
        start_row_id=next_row_id,
    )
    records.extend(crispr_records)
    records.extend(drug_records)

    pair_writer.close()
    snapshot_outputs.append(pair_list_path)

    if len(records) == 0:
        assertions.append(
            {
                "name": "non_empty_pairing",
                "pass": False,
                "details": {
                    "rules": ["At least one treated row must survive pairing."],
                    "crispr_summary": crispr_summary,
                    "drug_summary": drug_summary,
                },
                "counterexamples": [],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] No treated rows survived pairing.", file=sys.stderr)
        return 8

    pairing_violations = check_pairing_contract(records)
    pairing_contract_pass = len(pairing_violations) == 0
    assertions.append(
        {
            "name": "pairing_contract_n_controls_used_min_pool_50_no_replacement",
            "pass": pairing_contract_pass,
            "details": {
                "rules": [
                    "n_controls_used == min(pool_size, 50)",
                    "controls per treated are unique (without replacement)",
                    "n_controls_used <= 50",
                ],
                "k_controls_max": K_CONTROLS_MAX,
                "n_records_checked": len(records),
            },
            "counterexamples": pairing_violations[:MAX_COUNTEREXAMPLES],
        }
    )

    try:
        delta_meta = compute_gene_delta_and_meta(
            records=records,
            crispr_counts=crispr_counts,
            drug_counts=drug_counts,
            gene_delta_path=gene_delta_path,
        )
    except Exception as exc:
        assertions.append(
            {
                "name": "gene_delta_build",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed building gene deltas: {exc}", file=sys.stderr)
        return 9

    write_csv(delta_meta, delta_meta_path)
    snapshot_outputs.extend([gene_delta_path, delta_meta_path])

    try:
        w, index_map, pathway_details, pathway_input_paths = load_pathway_contract(
            task1_snapshot=task1_snapshot,
            shared_var_names=shared_var_names["gene_symbol"],
        )
        input_paths.extend(pathway_input_paths)

        pathway_proj_stats = project_pathway_delta(
            gene_delta_path=gene_delta_path,
            pathway_delta_path=pathway_delta_path,
            w=w,
            index_map=index_map,
        )
    except Exception as exc:
        assertions.append(
            {
                "name": "pathway_projection_contract",
                "pass": False,
                "details": {
                    "rules": [
                        "Use frozen project_on_load policy",
                        "W shape must be (2477, 50)",
                        "At least one pathway-alignment gene must map to K562 shared_var_names",
                    ],
                    "error": str(exc),
                },
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed building pathway deltas: {exc}", file=sys.stderr)
        return 10

    assertions.append(
        {
            "name": "pathway_projection_contract",
            "pass": True,
            "details": {
                "rules": [
                    "Use frozen project_on_load policy",
                    "W shape must be (2477, 50)",
                    "At least one pathway-alignment gene must map to K562 shared_var_names",
                ],
                "pathway_details": pathway_details,
                "projection_stats": pathway_proj_stats,
            },
            "counterexamples": [],
        }
    )
    snapshot_outputs.append(pathway_delta_path)

    common_targets = (
        common_targets_df["Common_Target"]
        .fillna("")
        .astype(str)
        .str.strip()
        .replace("", np.nan)
        .dropna()
        .tolist()
    )
    common_targets = sorted(set(common_targets))

    coverage_df = build_task2_pairs_coverage(records=records, common_targets=common_targets)
    inventory_df = build_post_build_inventory(crispr_summary=crispr_summary, drug_summary=drug_summary)

    attrition_df = pd.DataFrame(attrition_rows)
    if attrition_df.empty:
        attrition_df = pd.DataFrame(
            columns=[
                "dataset_side",
                "perturbation_class",
                "treated_cell_id",
                "target_raw",
                "time",
                "dose_value",
                "pool_size",
                "reason",
            ]
        )
    else:
        attrition_df = attrition_df.sort_values(
            ["dataset_side", "treated_cell_id"], kind="mergesort"
        ).reset_index(drop=True)

    coverage_path = (stage_dir / "task2_pairs_coverage.csv").resolve()
    inventory_path = (stage_dir / "task2_post_build_inventory.csv").resolve()
    attrition_path = (stage_dir / "k562_attrition.csv").resolve()

    write_csv(coverage_df, coverage_path)
    write_csv(inventory_df, inventory_path)
    write_csv(attrition_df, attrition_path)

    stage_outputs.extend([coverage_path, inventory_path, attrition_path])

    output_routing_snapshot_pass = all(
        path.resolve().is_relative_to(task2_snapshot.resolve()) for path in snapshot_outputs
    )
    output_routing_stage_pass = all(
        path.resolve().is_relative_to(stage_dir.resolve()) for path in stage_outputs
    )
    assertions.append(
        {
            "name": "output_routing_isolated",
            "pass": bool(output_routing_snapshot_pass and output_routing_stage_pass),
            "details": {
                "rules": [
                    "Snapshot artifacts must be under config.paths.task2_snapshot",
                    "Run artifacts must be under runs/<run_id>/s3_build_task2_snapshot",
                ],
                "task2_snapshot": str(task2_snapshot),
                "stage_dir": str(stage_dir),
            },
            "counterexamples": []
            if output_routing_snapshot_pass and output_routing_stage_pass
            else [
                {"bad_path": str(path)}
                for path in snapshot_outputs + stage_outputs
                if not (
                    path.resolve().is_relative_to(task2_snapshot.resolve())
                    or path.resolve().is_relative_to(stage_dir.resolve())
                )
            ][:MAX_COUNTEREXAMPLES],
        }
    )

    allowed_input_roots = [raw_k562_dir.resolve(), task1_snapshot.resolve()]
    bad_inputs = [
        str(path)
        for path in sorted(set(input_paths))
        if not any(path.absolute().is_relative_to(root) for root in allowed_input_roots)
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "rules": [
                    "K562 biological inputs are read from config.paths.raw_k562_dir",
                    "Pathway projection assets are read from config.paths.task1_snapshot",
                ],
                "raw_k562_dir": str(raw_k562_dir),
                "task1_snapshot": str(task1_snapshot),
                "n_inputs": len(sorted(set(input_paths))),
            },
            "counterexamples": [{"bad_input": p} for p in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    completed_at = utc_now_iso()

    run_manifest_path = (stage_dir / "run_manifest.json").resolve()
    audit_assertions_path = (stage_dir / "audit_assertions.json").resolve()
    manifest_path = (stage_dir / "manifest.json").resolve()

    output_paths = [
        str(path.resolve())
        for path in (
            snapshot_outputs
            + stage_outputs
            + [run_manifest_path, audit_assertions_path, manifest_path]
        )
    ]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s3_build_task2_snapshot.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "k_controls_max": K_CONTROLS_MAX,
            "min_controls_required": MIN_CONTROLS_REQUIRED,
            "raw_k562_dir": str(raw_k562_dir),
            "task1_snapshot": str(task1_snapshot),
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths))],
        "outputs": output_paths,
        "pairing_summary": {
            "crispr": crispr_summary,
            "drug": drug_summary,
            "n_treated_total": int(crispr_summary["n_treated_total"] + drug_summary["n_treated_total"]),
            "n_treated_kept": int(crispr_summary["n_treated_kept"] + drug_summary["n_treated_kept"]),
            "n_pairs_written": int(crispr_summary["n_pairs_written"] + drug_summary["n_pairs_written"]),
        },
        "pathway_projection": pathway_details,
        "raw_materialization": materialization_rows,
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

    if not pairing_contract_pass:
        print("[ERROR] Pairing contract assertion failed.", file=sys.stderr)
        return 11

    if not (output_routing_snapshot_pass and output_routing_stage_pass):
        print("[ERROR] Output routing assertion failed.", file=sys.stderr)
        return 12

    if bad_inputs:
        print("[ERROR] Input path isolation assertion failed.", file=sys.stderr)
        return 13

    return 0


if __name__ == "__main__":
    sys.exit(main())