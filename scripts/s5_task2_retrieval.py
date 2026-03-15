# SCRIPT_HEADER_CONTRACT
# Script: scripts/s5_task2_retrieval.py
# Legacy note:
#   - Preserved historical scPerturb-K562 S5 artifact with context-based G2C semantics.
#   - Corrected multisource Task2 S5 target-level semantics live in docs/contracts/task2_spec.md, docs/contracts/output-schemas.md, and scripts/s5_task2_retrieval_multisource.py.
# Purpose: 计算 legacy/interim Task2 实例级双向检索指标 (C2G, G2C) 及其 multi-positive chance-correction
# Inputs:
#   - Task2 Snapshot K562: config/config.yaml::paths.task2_snapshot
#     - k562/Common_Targets_K562.csv
#     - k562/derived/delta_meta.csv
#     - k562/Drug_meta.csv
#     - k562/derived/{gene_delta.npy,pathway_delta.npy}
#     - k562/fm/<representation>/{fm_delta.npy,fm_delta_meta.csv}
#     - k562/fm/<representation>/delta_operator_policy.json (optional recorded artifact only)
# Outputs:
#   - task2_retrieval_per_query.parquet: runs/<run_id>/s5_task2_retrieval/
#   - task2_retrieval_summary.csv: runs/<run_id>/s5_task2_retrieval/
#   - task2_retrieval_attrition.csv: runs/<run_id>/s5_task2_retrieval/
#   - task2_chance_identity_check.csv: runs/<run_id>/s5_task2_retrieval/
#   - task2_retrieval_summary_long.csv: runs/<run_id>/s5_task2_retrieval/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s5_task2_retrieval/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s5_task2_retrieval.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Missing required Task2 inputs -> exit non-zero
#   - Row-alignment / join / shape / finite-mask contract violation -> exit non-zero
#   - Chance identity check tolerance > 1e-12 -> exit non-zero
# Last Updated: 2026-03-13

"""
Inputs:
- data/task2_snapshot_v1/k562/Common_Targets_K562.csv
- data/task2_snapshot_v1/k562/derived/delta_meta.csv
- data/task2_snapshot_v1/k562/Drug_meta.csv
- data/task2_snapshot_v1/k562/derived/gene_delta.npy
- data/task2_snapshot_v1/k562/derived/pathway_delta.npy
- data/task2_snapshot_v1/k562/fm/<representation>/fm_delta.npy
- data/task2_snapshot_v1/k562/fm/<representation>/fm_delta_meta.csv
- data/task2_snapshot_v1/k562/fm/<representation>/delta_operator_policy.json (optional)

Outputs:
- task2_retrieval_per_query.parquet
- task2_retrieval_summary.csv
- task2_retrieval_attrition.csv
- task2_chance_identity_check.csv
- task2_retrieval_summary_long.csv
- run_manifest.json
- audit_assertions.json
- manifest.json

Frozen constants:
- GLOBAL_SEED = 619
- HIT_K_LIST = (1, 5, 10)
- CHANCE_IDENTITY_TOL = 1e-12
- QUERY_BATCH_ROWS = 1024
- score(q, c) = -||q - c||^2
- rank_true = 1 + count(score_all > best_positive_score)
- Canonical target order = file order from Common_Targets_K562.csv
- Canonical representation order = Gene, Pathway, scgpt, geneformer, scbert, scfoundation, uce, state, tahoe-x1

Attrition rules:
- Representation valid_mask is applied before query selection and gallery construction.
- m_pos == 0 after valid intersection => soft-fail exclusion with reason=m_pos0_after_intersection.
- Empty gallery after valid filtering => soft-fail summary row with NA metrics and attrition, not global termination.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

STAGE = "s5_task2_retrieval"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
EXPECTED_CELL_LINE = "K562"
EXPECTED_DATASET = "scPerturb"

GLOBAL_SEED = 619
HIT_K_LIST: Tuple[int, ...] = (1, 5, 10)
CHANCE_IDENTITY_TOL = 1e-12
QUERY_BATCH_ROWS = 1024
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

DIRECTION_ORDER: Tuple[str, str] = ("C2G", "G2C")
SUMMARY_LONG_METRICS: Tuple[str, ...] = (
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
)

C2G_GALLERY_DEFINITION = "C2G_genetic_target_centroid_gallery_v1"
G2C_GALLERY_DEFINITION = "G2C_chemical_context_centroid_gallery_v1"
C2G_POS_DEFINITION = "C2G_tokens_intersect_gallery_targets_v1"
G2C_POS_DEFINITION = "G2C_target_in_context_tokens_v1"

SUMMARY_COLUMNS = [
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

CHANCE_COLUMNS = [
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

ATTRITION_COLUMNS = [
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

SUMMARY_LONG_COLUMNS = [
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
    invalid_by_side: Dict[str, int]


@dataclass(frozen=True)
class GalleryPayload:
    gallery_ids: List[str]
    member_row_ids: List[np.ndarray]
    centroids: np.ndarray
    norms: np.ndarray
    gallery_hash: str


@dataclass
class DirectionAccumulator:
    dataset: str
    cell_line: str
    direction: str
    representation: str
    gallery_definition_id: str
    pos_definition_id: str
    n_total: int
    n_gallery: int
    n_valid: int = 0
    n_mpos0: int = 0
    sum_mrr_raw: float = 0.0
    sum_expected_mrr: float = 0.0
    sum_mrr_corrected: float = 0.0
    sum_hit_raw: Dict[int, float] = field(default_factory=lambda: {k: 0.0 for k in HIT_K_LIST})
    sum_hit_expected: Dict[int, float] = field(default_factory=lambda: {k: 0.0 for k in HIT_K_LIST})
    sum_hit_corrected: Dict[int, float] = field(default_factory=lambda: {k: 0.0 for k in HIT_K_LIST})
    m_pos_values: List[int] = field(default_factory=list)

    def add_valid(
        self,
        *,
        m_pos: int,
        mrr_raw: float,
        expected_mrr: float,
        mrr_corrected: float,
        hit_raw: Mapping[int, float],
        hit_expected: Mapping[int, float],
        hit_corrected: Mapping[int, float],
    ) -> None:
        self.n_valid += 1
        self.sum_mrr_raw += float(mrr_raw)
        self.sum_expected_mrr += float(expected_mrr)
        self.sum_mrr_corrected += float(mrr_corrected)
        self.m_pos_values.append(int(m_pos))
        for k in HIT_K_LIST:
            self.sum_hit_raw[k] += float(hit_raw[k])
            self.sum_hit_expected[k] += float(hit_expected[k])
            self.sum_hit_corrected[k] += float(hit_corrected[k])


class PerQueryParquetSink:
    def __init__(self, output_path: Path, schema: pa.Schema) -> None:
        self.output_path = output_path
        self.schema = schema
        self.writer: Optional[pq.ParquetWriter] = None

    def write(self, frame: pd.DataFrame) -> None:
        if frame.empty:
            return
        table = pa.Table.from_pandas(
            canonicalize_per_query_frame(frame),
            schema=self.schema,
            preserve_index=False,
        )
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)
        self.writer.write_table(table)

    def close(self) -> None:
        if self.writer is not None:
            self.writer.close()
            return
        empty_data = {field.name: pa.array([], type=field.type) for field in self.schema}
        empty_table = pa.Table.from_pydict(empty_data, schema=self.schema)
        pq.write_table(empty_table, self.output_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S5 Task2 opposite-mechanism retrieval")
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


def init_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def normalize_scalar_text(value: object, default: str = "NA") -> str:
    if pd.isna(value):
        return default
    text = str(value).strip()
    return text if text else default


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


def hash_ordered_ids(ids: Sequence[str]) -> str:
    material = "||".join(ids)
    return hashlib.sha256(material.encode("utf-8")).hexdigest()


def make_context_id(cell_line: str, perturbation_raw: str, time_value: float, dose_value: float) -> str:
    payload = [str(cell_line), str(perturbation_raw), float(time_value), float(dose_value)]
    return json.dumps(payload, ensure_ascii=True, separators=(",", ":"))


def per_query_schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("dataset", pa.string()),
            pa.field("cell_line", pa.string()),
            pa.field("direction", pa.string()),
            pa.field("representation", pa.string()),
            pa.field("query_row_id", pa.int64()),
            pa.field("query_uid", pa.string()),
            pa.field("treated_cell_id", pa.string()),
            pa.field("query_target_raw", pa.string()),
            pa.field("query_target_tokens", pa.string()),
            pa.field("query_target_token", pa.string()),
            pa.field("query_n_targets", pa.int64()),
            pa.field("query_time", pa.float64()),
            pa.field("query_dose_value", pa.float64()),
            pa.field("query_specificity_tier", pa.string()),
            pa.field("gallery_definition_id", pa.string()),
            pa.field("pos_definition_id", pa.string()),
            pa.field("gallery_ids_hash", pa.string()),
            pa.field("N_gallery", pa.int64()),
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
        ]
    )


def canonicalize_per_query_frame(frame: pd.DataFrame) -> pd.DataFrame:
    int_cols = ["query_row_id", "query_n_targets", "N_gallery", "m_pos", "rank_true"]
    float_cols = [
        "query_time",
        "query_dose_value",
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
        "dataset",
        "cell_line",
        "direction",
        "representation",
        "query_uid",
        "treated_cell_id",
        "query_target_raw",
        "query_target_tokens",
        "query_target_token",
        "query_specificity_tier",
        "gallery_definition_id",
        "pos_definition_id",
        "gallery_ids_hash",
    ]

    out = frame.copy()
    for col in int_cols:
        out[col] = out[col].astype(np.int64)
    for col in float_cols:
        out[col] = out[col].astype(np.float64)
    for col in str_cols:
        out[col] = out[col].fillna("NA").astype(str)
    return out


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
        token = normalize_scalar_text(value, default="")
        if not token:
            raise ValueError("Common_Targets_K562.csv contains empty target")
        if token in seen:
            raise ValueError(f"Common_Targets_K562.csv contains duplicate target: {token}")
        seen.add(token)
        target_order.append(token)

    if not target_order:
        raise ValueError("Common_Targets_K562.csv produced empty target order")
    return target_order


def load_delta_meta(delta_meta_path: Path, target_order: Sequence[str]) -> pd.DataFrame:
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
            "target_raw",
            "time",
            "dose_value",
            "specificity_tier",
            "n_controls_used",
            "dataset_side",
            "target_tokens",
        ],
        "delta_meta.csv",
    )

    delta_meta["row_id"] = pd.to_numeric(delta_meta["row_id"], errors="raise").astype(np.int64)
    delta_meta["n_controls_used"] = pd.to_numeric(delta_meta["n_controls_used"], errors="raise").astype(np.int64)
    delta_meta["time"] = pd.to_numeric(delta_meta["time"], errors="coerce").astype(np.float64)
    delta_meta["dose_value"] = pd.to_numeric(delta_meta["dose_value"], errors="coerce").astype(np.float64)
    delta_meta["treated_cell_id"] = delta_meta["treated_cell_id"].astype(str)
    delta_meta["perturbation_class"] = delta_meta["perturbation_class"].astype(str).str.strip()
    delta_meta["cell_line"] = delta_meta["cell_line"].astype(str).str.strip()
    delta_meta["target_raw"] = delta_meta["target_raw"].astype(str).str.strip()
    delta_meta["dataset_side"] = delta_meta["dataset_side"].astype(str).str.strip().str.upper()
    delta_meta["specificity_tier"] = delta_meta["specificity_tier"].map(normalize_scalar_text)
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

    parsed_tokens = [parse_target_tokens_field(value) for value in delta_meta["target_tokens"].tolist()]
    token_universe = {token for tokens in parsed_tokens for token in tokens}
    if token_universe != set(target_order):
        raise ValueError(
            "Atomic token set from delta_meta.target_tokens does not match Common_Targets_K562.csv"
        )

    delta_meta["parsed_target_tokens"] = parsed_tokens
    delta_meta["n_targets"] = delta_meta["parsed_target_tokens"].map(len).astype(np.int64)

    crispr_bad = delta_meta.loc[
        (delta_meta["dataset_side"] == "CRISPR")
        & (delta_meta["n_targets"] != 1),
        ["row_id", "target_tokens"],
    ]
    if not crispr_bad.empty:
        raise ValueError(
            "CRISPR rows in delta_meta.csv must have exactly one target token, "
            f"examples={crispr_bad.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )

    return delta_meta


def load_drug_context_frame(drug_meta_path: Path, delta_meta: pd.DataFrame) -> pd.DataFrame:
    if not drug_meta_path.is_file():
        raise FileNotFoundError(f"Missing Drug_meta.csv: {drug_meta_path}")
    drug_meta = pd.read_csv(drug_meta_path)
    if "treated_cell_id" not in drug_meta.columns and "Unnamed: 0" in drug_meta.columns:
        drug_meta = drug_meta.rename(columns={"Unnamed: 0": "treated_cell_id"})

    ensure_required_columns(
        drug_meta,
        ["treated_cell_id", "perturbation_raw", "time", "dose_value", "clean_target_mapped"],
        "Drug_meta.csv",
    )

    drug_meta["treated_cell_id"] = drug_meta["treated_cell_id"].astype(str)
    drug_meta["perturbation_raw"] = drug_meta["perturbation_raw"].map(lambda value: normalize_scalar_text(value, default=""))
    drug_meta["clean_target_mapped"] = drug_meta["clean_target_mapped"].map(normalize_scalar_text)
    drug_meta["time"] = pd.to_numeric(drug_meta["time"], errors="coerce").astype(np.float64)
    drug_meta["dose_value"] = pd.to_numeric(drug_meta["dose_value"], errors="coerce").astype(np.float64)

    if drug_meta["treated_cell_id"].duplicated().any():
        dup = drug_meta.loc[drug_meta["treated_cell_id"].duplicated(), "treated_cell_id"].head(MAX_COUNTEREXAMPLES)
        raise ValueError(f"Drug_meta.csv treated_cell_id has duplicates, examples={dup.tolist()}")
    if (drug_meta["perturbation_raw"] == "").any():
        examples = drug_meta.loc[drug_meta["perturbation_raw"] == "", "treated_cell_id"].head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"Drug_meta.csv contains empty perturbation_raw, examples={examples}")

    drug_rows = delta_meta.loc[delta_meta["dataset_side"] == "DRUG"].copy()
    merged = drug_rows.merge(
        drug_meta[["treated_cell_id", "perturbation_raw", "time", "dose_value", "clean_target_mapped"]],
        on="treated_cell_id",
        how="left",
        validate="one_to_one",
        sort=False,
        suffixes=("_delta", "_raw"),
    )

    if merged["perturbation_raw"].isna().any():
        examples = merged.loc[merged["perturbation_raw"].isna(), "treated_cell_id"].head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"Drug_meta.csv join missing perturbation_raw for treated_cell_id examples={examples}")

    target_match = merged["clean_target_mapped"].astype(str) == merged["target_raw"].astype(str)
    if not bool(target_match.all()):
        examples = merged.loc[~target_match, ["treated_cell_id", "clean_target_mapped", "target_raw"]]
        raise ValueError(
            "Drug_meta.csv clean_target_mapped mismatch against delta_meta target_raw, "
            f"examples={examples.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )

    time_match = np.isclose(
        merged["time_delta"].to_numpy(dtype=np.float64),
        merged["time_raw"].to_numpy(dtype=np.float64),
        equal_nan=True,
    )
    dose_match = np.isclose(
        merged["dose_value_delta"].to_numpy(dtype=np.float64),
        merged["dose_value_raw"].to_numpy(dtype=np.float64),
        equal_nan=True,
    )

    merged = merged.rename(
        columns={
            "time_delta": "delta_time",
            "dose_value_delta": "delta_dose_value",
        }
    )

    if not bool(time_match.all()):
        examples = merged.loc[~time_match, ["treated_cell_id", "delta_time"]]
        raise ValueError(
            "Drug_meta.csv time mismatch against delta_meta, "
            f"examples={examples.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )
    if not bool(dose_match.all()):
        examples = merged.loc[~dose_match, ["treated_cell_id", "delta_dose_value"]]
        raise ValueError(
            "Drug_meta.csv dose_value mismatch against delta_meta, "
            f"examples={examples.head(MAX_COUNTEREXAMPLES).to_dict(orient='records')}"
        )

    context_ids = [
        make_context_id(
            EXPECTED_CELL_LINE,
            str(perturbation_raw),
            float(time_value),
            float(dose_value),
        )
        for perturbation_raw, time_value, dose_value in zip(
            merged["perturbation_raw"].tolist(),
            merged["delta_time"].tolist(),
            merged["delta_dose_value"].tolist(),
        )
    ]

    out = pd.DataFrame(
        {
            "row_id": merged["row_id"].to_numpy(dtype=np.int64, copy=False),
            "treated_cell_id": merged["treated_cell_id"].astype(str).to_numpy(),
            "perturbation_raw": merged["perturbation_raw"].astype(str).to_numpy(),
            "time": merged["delta_time"].to_numpy(dtype=np.float64, copy=False),
            "dose_value": merged["delta_dose_value"].to_numpy(dtype=np.float64, copy=False),
            "context_id": np.asarray(context_ids, dtype=object),
        }
    )
    return out.sort_values(
        ["perturbation_raw", "time", "dose_value", "row_id"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_genetic_rows_by_target(
    delta_meta: pd.DataFrame,
    target_order: Sequence[str],
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    chem_rows: List[int] = []
    gen_rows: List[int] = []
    gen_by_target: MutableMapping[str, List[int]] = {target: [] for target in target_order}

    for row in delta_meta.itertuples(index=False):
        row_id = int(row.row_id)
        if row.dataset_side == "DRUG":
            chem_rows.append(row_id)
        else:
            gen_rows.append(row_id)
            gen_by_target[row.parsed_target_tokens[0]].append(row_id)

    return (
        np.asarray(chem_rows, dtype=np.int64),
        np.asarray(gen_rows, dtype=np.int64),
        {target: np.asarray(rows, dtype=np.int64) for target, rows in gen_by_target.items()},
    )


def build_all_context_items(drug_context_frame: pd.DataFrame) -> List[Tuple[str, np.ndarray]]:
    context_items: List[Tuple[str, np.ndarray]] = []
    for context_id, sub in drug_context_frame.groupby("context_id", sort=False):
        row_ids = sub["row_id"].to_numpy(dtype=np.int64, copy=False)
        context_items.append((str(context_id), row_ids))
    return context_items


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
    truth_row_ids = delta_meta["row_id"].to_numpy(dtype=np.int64, copy=False)

    if spec.family == "derived":
        ensure_required_columns(meta, ["row_id", "treated_cell_id"], f"{spec.name} metadata")
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

        if not np.array_equal(meta["row_id"].to_numpy(dtype=np.int64, copy=False), truth_row_ids):
            raise ValueError(f"{spec.name} fm_delta_meta row_id mismatch against delta_meta.csv")
        if not np.array_equal(
            meta["treated_cell_id"].to_numpy(dtype=object),
            delta_meta["treated_cell_id"].to_numpy(dtype=object),
        ):
            raise ValueError(f"{spec.name} treated_cell_id mismatch against delta_meta.csv")
        if not np.array_equal(
            meta["dataset_side"].to_numpy(dtype=object),
            delta_meta["dataset_side"].to_numpy(dtype=object),
        ):
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
            invalid_by_side=invalid_by_side,
        ),
        rep_inputs,
        details,
    )


def build_centroids(
    arr: np.ndarray,
    ordered_groups: Sequence[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray]:
    if not ordered_groups:
        return np.empty((0, int(arr.shape[1])), dtype=np.float64), np.empty((0,), dtype=np.float64)

    centroids: List[np.ndarray] = []
    for row_ids in ordered_groups:
        block = np.asarray(arr[row_ids.astype(np.int64, copy=False)], dtype=np.float64)
        if block.ndim != 2 or block.shape[0] == 0:
            raise ValueError("Gallery centroid block must contain at least one row")
        centroid = np.mean(block, axis=0, dtype=np.float64)
        if not bool(np.isfinite(centroid).all()):
            raise ValueError("Gallery centroid contains non-finite values after valid_mask filtering")
        centroids.append(centroid)

    centroid_matrix = np.stack(centroids, axis=0).astype(np.float64, copy=False)
    norms = np.sum(centroid_matrix * centroid_matrix, axis=1, dtype=np.float64)
    return centroid_matrix, norms


def build_exact_chance_lookup(n_gallery: int) -> Dict[int, Dict[str, float]]:
    lookup: Dict[int, Dict[str, float]] = {}
    if n_gallery <= 0:
        return lookup

    for m_pos in range(1, n_gallery + 1):
        total = math.comb(n_gallery, m_pos)
        max_rank = n_gallery - m_pos + 1
        probs: List[float] = []
        for rank in range(1, max_rank + 1):
            probs.append(math.comb(n_gallery - rank, m_pos - 1) / total)

        expected_mrr = float(
            sum(prob / float(rank) for rank, prob in enumerate(probs, start=1))
        )
        expected_hit = {
            k: float(sum(probs[: min(int(k), len(probs))]))
            for k in HIT_K_LIST
        }
        lookup[m_pos] = {
            "expected_mrr": expected_mrr,
            "expected_hit1": expected_hit[1],
            "expected_hit5": expected_hit[5],
            "expected_hit10": expected_hit[10],
        }
    return lookup


def make_attrition_row(
    *,
    direction: str,
    representation: str,
    reason: str,
    n_query_rows_before: int,
    n_query_rows_after: int,
    n_gallery_member_rows_before: int,
    n_gallery_member_rows_after: int,
    n_gallery_items_before: int,
    n_gallery_items_after: int,
    notes: str,
) -> Dict[str, object]:
    return {
        "dataset": EXPECTED_DATASET,
        "cell_line": EXPECTED_CELL_LINE,
        "direction": direction,
        "representation": representation,
        "reason": reason,
        "n_query_rows_before": int(n_query_rows_before),
        "n_query_rows_after": int(n_query_rows_after),
        "n_query_rows_removed": int(n_query_rows_before - n_query_rows_after),
        "n_gallery_member_rows_before": int(n_gallery_member_rows_before),
        "n_gallery_member_rows_after": int(n_gallery_member_rows_after),
        "n_gallery_member_rows_removed": int(n_gallery_member_rows_before - n_gallery_member_rows_after),
        "n_gallery_items_before": int(n_gallery_items_before),
        "n_gallery_items_after": int(n_gallery_items_after),
        "n_gallery_items_removed": int(n_gallery_items_before - n_gallery_items_after),
        "notes": notes,
    }


def maybe_append_attrition(
    rows: List[Dict[str, object]],
    payload: Mapping[str, object],
) -> None:
    if int(payload["n_query_rows_removed"]) <= 0 and int(payload["n_gallery_member_rows_removed"]) <= 0 and int(payload["n_gallery_items_removed"]) <= 0:
        return
    rows.append(dict(payload))


def make_per_query_row(
    *,
    direction: str,
    representation: str,
    row_id: int,
    gallery_definition_id: str,
    pos_definition_id: str,
    gallery_ids_hash: str,
    n_gallery: int,
    m_pos: int,
    rank_true: int,
    expected: Mapping[str, float],
    target_raw: Sequence[str],
    target_tokens_text: Sequence[str],
    query_n_targets: Sequence[int],
    query_time: np.ndarray,
    query_dose_value: np.ndarray,
    specificity_tier: Sequence[str],
    treated_cell_id: Sequence[str],
    parsed_target_tokens: Sequence[Tuple[str, ...]],
) -> Dict[str, object]:
    mrr_raw = 1.0 / float(rank_true)
    hit_raw = {k: float(rank_true <= k) for k in HIT_K_LIST}
    hit_expected = {
        1: float(expected["expected_hit1"]),
        5: float(expected["expected_hit5"]),
        10: float(expected["expected_hit10"]),
    }
    hit_corrected = {k: float(hit_raw[k] - hit_expected[k]) for k in HIT_K_LIST}
    mrr_corrected = float(mrr_raw - float(expected["expected_mrr"]))
    query_target_token = parsed_target_tokens[row_id][0] if direction == "G2C" else "NA"

    return {
        "dataset": EXPECTED_DATASET,
        "cell_line": EXPECTED_CELL_LINE,
        "direction": direction,
        "representation": representation,
        "query_row_id": int(row_id),
        "query_uid": str(row_id),
        "treated_cell_id": str(treated_cell_id[row_id]),
        "query_target_raw": str(target_raw[row_id]),
        "query_target_tokens": str(target_tokens_text[row_id]),
        "query_target_token": query_target_token,
        "query_n_targets": int(query_n_targets[row_id]),
        "query_time": float(query_time[row_id]),
        "query_dose_value": float(query_dose_value[row_id]),
        "query_specificity_tier": str(specificity_tier[row_id]),
        "gallery_definition_id": gallery_definition_id,
        "pos_definition_id": pos_definition_id,
        "gallery_ids_hash": gallery_ids_hash,
        "N_gallery": int(n_gallery),
        "m_pos": int(m_pos),
        "rank_true": int(rank_true),
        "mrr_raw": float(mrr_raw),
        "hit1_raw": float(hit_raw[1]),
        "hit5_raw": float(hit_raw[5]),
        "hit10_raw": float(hit_raw[10]),
        "expected_mrr_chance": float(expected["expected_mrr"]),
        "expected_hit1_chance": float(expected["expected_hit1"]),
        "expected_hit5_chance": float(expected["expected_hit5"]),
        "expected_hit10_chance": float(expected["expected_hit10"]),
        "mrr_corrected": float(mrr_corrected),
        "hit1_corrected": float(hit_corrected[1]),
        "hit5_corrected": float(hit_corrected[5]),
        "hit10_corrected": float(hit_corrected[10]),
    }


def build_summary_row(acc: DirectionAccumulator) -> Dict[str, object]:
    n_excluded = int(acc.n_total - acc.n_valid)
    if acc.n_valid > 0:
        m_pos_arr = np.asarray(acc.m_pos_values, dtype=np.float64)
        n_gallery_mean = float(acc.n_gallery)
        mean_mrr_raw = acc.sum_mrr_raw / acc.n_valid
        mean_expected_mrr = acc.sum_expected_mrr / acc.n_valid
        mean_mrr_corrected = acc.sum_mrr_corrected / acc.n_valid
        mean_hit_raw = {k: acc.sum_hit_raw[k] / acc.n_valid for k in HIT_K_LIST}
        mean_hit_expected = {k: acc.sum_hit_expected[k] / acc.n_valid for k in HIT_K_LIST}
        mean_hit_corrected = {k: acc.sum_hit_corrected[k] / acc.n_valid for k in HIT_K_LIST}
        m_pos_mean = float(np.mean(m_pos_arr))
        m_pos_p50 = float(np.percentile(m_pos_arr, 50))
        m_pos_p90 = float(np.percentile(m_pos_arr, 90))
    else:
        n_gallery_mean = float(acc.n_gallery) if acc.n_gallery > 0 else 0.0
        mean_mrr_raw = np.nan
        mean_expected_mrr = np.nan
        mean_mrr_corrected = np.nan
        mean_hit_raw = {k: np.nan for k in HIT_K_LIST}
        mean_hit_expected = {k: np.nan for k in HIT_K_LIST}
        mean_hit_corrected = {k: np.nan for k in HIT_K_LIST}
        m_pos_mean = np.nan
        m_pos_p50 = np.nan
        m_pos_p90 = np.nan

    return {
        "dataset": acc.dataset,
        "cell_line": acc.cell_line,
        "direction": acc.direction,
        "representation": acc.representation,
        "gallery_definition_id": acc.gallery_definition_id,
        "pos_definition_id": acc.pos_definition_id,
        "n_total": int(acc.n_total),
        "n_valid": int(acc.n_valid),
        "n_excluded_missing_metric_or_mpos0": int(n_excluded),
        "N_gallery_mean": float(n_gallery_mean),
        "N_gallery_max": int(acc.n_gallery),
        "m_pos_mean": float(m_pos_mean) if np.isfinite(m_pos_mean) else np.nan,
        "m_pos_p50": float(m_pos_p50) if np.isfinite(m_pos_p50) else np.nan,
        "m_pos_p90": float(m_pos_p90) if np.isfinite(m_pos_p90) else np.nan,
        "mean_mrr_raw": float(mean_mrr_raw) if np.isfinite(mean_mrr_raw) else np.nan,
        "mean_expected_mrr_chance": float(mean_expected_mrr) if np.isfinite(mean_expected_mrr) else np.nan,
        "mean_mrr_corrected": float(mean_mrr_corrected) if np.isfinite(mean_mrr_corrected) else np.nan,
        "mean_hit1_raw": float(mean_hit_raw[1]) if np.isfinite(mean_hit_raw[1]) else np.nan,
        "mean_expected_hit1_chance": float(mean_hit_expected[1]) if np.isfinite(mean_hit_expected[1]) else np.nan,
        "mean_hit1_corrected": float(mean_hit_corrected[1]) if np.isfinite(mean_hit_corrected[1]) else np.nan,
        "mean_hit5_raw": float(mean_hit_raw[5]) if np.isfinite(mean_hit_raw[5]) else np.nan,
        "mean_expected_hit5_chance": float(mean_hit_expected[5]) if np.isfinite(mean_hit_expected[5]) else np.nan,
        "mean_hit5_corrected": float(mean_hit_corrected[5]) if np.isfinite(mean_hit_corrected[5]) else np.nan,
        "mean_hit10_raw": float(mean_hit_raw[10]) if np.isfinite(mean_hit_raw[10]) else np.nan,
        "mean_expected_hit10_chance": float(mean_hit_expected[10]) if np.isfinite(mean_hit_expected[10]) else np.nan,
        "mean_hit10_corrected": float(mean_hit_corrected[10]) if np.isfinite(mean_hit_corrected[10]) else np.nan,
    }


def build_chance_row(summary_row: Mapping[str, object]) -> Dict[str, object]:
    if int(summary_row["n_valid"]) > 0:
        delta_mrr = float(summary_row["mean_mrr_corrected"]) - (
            float(summary_row["mean_mrr_raw"]) - float(summary_row["mean_expected_mrr_chance"])
        )
        delta_hit1 = float(summary_row["mean_hit1_corrected"]) - (
            float(summary_row["mean_hit1_raw"]) - float(summary_row["mean_expected_hit1_chance"])
        )
        delta_hit5 = float(summary_row["mean_hit5_corrected"]) - (
            float(summary_row["mean_hit5_raw"]) - float(summary_row["mean_expected_hit5_chance"])
        )
        delta_hit10 = float(summary_row["mean_hit10_corrected"]) - (
            float(summary_row["mean_hit10_raw"]) - float(summary_row["mean_expected_hit10_chance"])
        )
    else:
        delta_mrr = np.nan
        delta_hit1 = np.nan
        delta_hit5 = np.nan
        delta_hit10 = np.nan

    return {
        "dataset": summary_row["dataset"],
        "cell_line": summary_row["cell_line"],
        "direction": summary_row["direction"],
        "representation": summary_row["representation"],
        "delta_mrr": float(delta_mrr) if np.isfinite(delta_mrr) else np.nan,
        "abs_delta_mrr": float(abs(delta_mrr)) if np.isfinite(delta_mrr) else np.nan,
        "delta_hit1": float(delta_hit1) if np.isfinite(delta_hit1) else np.nan,
        "abs_delta_hit1": float(abs(delta_hit1)) if np.isfinite(delta_hit1) else np.nan,
        "delta_hit5": float(delta_hit5) if np.isfinite(delta_hit5) else np.nan,
        "abs_delta_hit5": float(abs(delta_hit5)) if np.isfinite(delta_hit5) else np.nan,
        "delta_hit10": float(delta_hit10) if np.isfinite(delta_hit10) else np.nan,
        "abs_delta_hit10": float(abs(delta_hit10)) if np.isfinite(delta_hit10) else np.nan,
    }


def sort_summary_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)
    order_map_dir = {name: i for i, name in enumerate(DIRECTION_ORDER)}
    order_map_rep = {name: i for i, name in enumerate(REPRESENTATION_ORDER)}
    out = frame.copy()
    out["_direction_order"] = out["direction"].map(order_map_dir)
    out["_representation_order"] = out["representation"].map(order_map_rep)
    out = out.sort_values(["_direction_order", "_representation_order"], kind="mergesort")
    out = out.drop(columns=["_direction_order", "_representation_order"])
    return out[SUMMARY_COLUMNS].reset_index(drop=True)


def sort_chance_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=CHANCE_COLUMNS)
    order_map_dir = {name: i for i, name in enumerate(DIRECTION_ORDER)}
    order_map_rep = {name: i for i, name in enumerate(REPRESENTATION_ORDER)}
    out = frame.copy()
    out["_direction_order"] = out["direction"].map(order_map_dir)
    out["_representation_order"] = out["representation"].map(order_map_rep)
    out = out.sort_values(["_direction_order", "_representation_order"], kind="mergesort")
    out = out.drop(columns=["_direction_order", "_representation_order"])
    return out[CHANCE_COLUMNS].reset_index(drop=True)


def sort_attrition_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=ATTRITION_COLUMNS)
    order_map_dir = {name: i for i, name in enumerate(DIRECTION_ORDER)}
    order_map_rep = {name: i for i, name in enumerate(REPRESENTATION_ORDER)}
    out = frame.copy()
    out["_direction_order"] = out["direction"].map(order_map_dir)
    out["_representation_order"] = out["representation"].map(order_map_rep)
    out = out.sort_values(
        ["_direction_order", "_representation_order", "reason"],
        kind="mergesort",
    )
    out = out.drop(columns=["_direction_order", "_representation_order"])
    return out[ATTRITION_COLUMNS].reset_index(drop=True)


def build_summary_long_rows(summary_rows: Sequence[Mapping[str, object]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for row in summary_rows:
        for metric_name in SUMMARY_LONG_METRICS:
            rows.append(
                {
                    "direction": row["direction"],
                    "dataset": row["dataset"],
                    "representation": row["representation"],
                    "metric_name": metric_name,
                    "metric_value": row[metric_name],
                    "n_valid": row["n_valid"],
                    "N_gallery_mean": row["N_gallery_mean"],
                    "N_gallery_max": row["N_gallery_max"],
                    "m_pos_mean": row["m_pos_mean"],
                    "m_pos_p50": row["m_pos_p50"],
                    "m_pos_p90": row["m_pos_p90"],
                }
            )
    return rows


def sort_summary_long_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=SUMMARY_LONG_COLUMNS)
    order_map_dir = {name: i for i, name in enumerate(DIRECTION_ORDER)}
    order_map_rep = {name: i for i, name in enumerate(REPRESENTATION_ORDER)}
    order_map_metric = {name: i for i, name in enumerate(SUMMARY_LONG_METRICS)}
    out = frame.copy()
    out["_direction_order"] = out["direction"].map(order_map_dir)
    out["_representation_order"] = out["representation"].map(order_map_rep)
    out["_metric_order"] = out["metric_name"].map(order_map_metric)
    out = out.sort_values(
        ["_direction_order", "_representation_order", "_metric_order"],
        kind="mergesort",
    )
    out = out.drop(columns=["_direction_order", "_representation_order", "_metric_order"])
    return out[SUMMARY_LONG_COLUMNS].reset_index(drop=True)


def run_c2g_direction(
    *,
    arr: np.ndarray,
    representation: str,
    valid_mask: np.ndarray,
    chem_all_row_ids: np.ndarray,
    gen_all_row_ids: np.ndarray,
    gen_rows_by_target: Mapping[str, np.ndarray],
    target_order: Sequence[str],
    target_raw: Sequence[str],
    target_tokens_text: Sequence[str],
    query_n_targets: Sequence[int],
    query_time: np.ndarray,
    query_dose_value: np.ndarray,
    specificity_tier: Sequence[str],
    treated_cell_id: Sequence[str],
    parsed_target_tokens: Sequence[Tuple[str, ...]],
    sink: PerQueryParquetSink,
) -> Tuple[Dict[str, object], Dict[str, object], List[Dict[str, object]], Dict[str, object]]:
    attrition_rows: List[Dict[str, object]] = []

    query_before = chem_all_row_ids
    query_after = query_before[valid_mask[query_before]] if query_before.size > 0 else np.empty((0,), dtype=np.int64)

    gallery_member_before = gen_all_row_ids
    gallery_member_after = gallery_member_before[valid_mask[gallery_member_before]] if gallery_member_before.size > 0 else np.empty((0,), dtype=np.int64)

    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="C2G",
            representation=representation,
            reason="representation_valid_mask_drop_query",
            n_query_rows_before=int(query_before.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_member_before.size),
            n_gallery_member_rows_after=int(gallery_member_before.size),
            n_gallery_items_before=int(len(target_order)),
            n_gallery_items_after=int(len(target_order)),
            notes="Query rows filtered by representation valid_mask before retrieval",
        ),
    )
    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="C2G",
            representation=representation,
            reason="representation_valid_mask_drop_gallery_member",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_member_before.size),
            n_gallery_member_rows_after=int(gallery_member_after.size),
            n_gallery_items_before=int(len(target_order)),
            n_gallery_items_after=int(len(target_order)),
            notes="Genetic gallery member rows filtered by representation valid_mask before centroid construction",
        ),
    )

    ordered_gallery_ids: List[str] = []
    ordered_member_rows: List[np.ndarray] = []
    for target in target_order:
        rows = gen_rows_by_target[target]
        valid_rows = rows[valid_mask[rows]] if rows.size > 0 else np.empty((0,), dtype=np.int64)
        if valid_rows.size > 0:
            ordered_gallery_ids.append(target)
            ordered_member_rows.append(valid_rows)

    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="C2G",
            representation=representation,
            reason="gallery_item_empty_after_filter",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_member_after.size),
            n_gallery_member_rows_after=int(gallery_member_after.size),
            n_gallery_items_before=int(len(target_order)),
            n_gallery_items_after=int(len(ordered_gallery_ids)),
            notes="Target centroids with zero valid genetic members are removed from the gallery",
        ),
    )

    gallery_payload = GalleryPayload(
        gallery_ids=ordered_gallery_ids,
        member_row_ids=ordered_member_rows,
        centroids=np.empty((0, int(arr.shape[1])), dtype=np.float64),
        norms=np.empty((0,), dtype=np.float64),
        gallery_hash=hash_ordered_ids(ordered_gallery_ids),
    )
    if ordered_member_rows:
        centroids, norms = build_centroids(arr, ordered_member_rows)
        gallery_payload = GalleryPayload(
            gallery_ids=ordered_gallery_ids,
            member_row_ids=ordered_member_rows,
            centroids=centroids,
            norms=norms,
            gallery_hash=hash_ordered_ids(ordered_gallery_ids),
        )

    acc = DirectionAccumulator(
        dataset=EXPECTED_DATASET,
        cell_line=EXPECTED_CELL_LINE,
        direction="C2G",
        representation=representation,
        gallery_definition_id=C2G_GALLERY_DEFINITION,
        pos_definition_id=C2G_POS_DEFINITION,
        n_total=int(query_after.size),
        n_gallery=int(len(gallery_payload.gallery_ids)),
    )

    if gallery_payload.gallery_ids:
        target_to_gallery_idx = {target: idx for idx, target in enumerate(gallery_payload.gallery_ids)}
        chance_lookup = build_exact_chance_lookup(len(gallery_payload.gallery_ids))
        for start in range(0, query_after.size, QUERY_BATCH_ROWS):
            end = min(start + QUERY_BATCH_ROWS, query_after.size)
            row_batch = query_after[start:end]
            q_batch = np.asarray(arr[row_batch.astype(np.int64, copy=False)], dtype=np.float64)
            if not bool(np.isfinite(q_batch).all()):
                raise ValueError(f"{representation}/C2G query batch contains non-finite values after valid_mask filtering")
            q_norm = np.sum(q_batch * q_batch, axis=1, dtype=np.float64)
            scores = -(q_norm[:, None] + gallery_payload.norms[None, :] - 2.0 * (q_batch @ gallery_payload.centroids.T))

            per_query_rows: List[Dict[str, object]] = []
            for batch_pos, row_id in enumerate(row_batch.tolist()):
                pos_idx = [target_to_gallery_idx[token] for token in parsed_target_tokens[row_id] if token in target_to_gallery_idx]
                if not pos_idx:
                    acc.n_mpos0 += 1
                    continue
                best_positive_score = float(np.max(scores[batch_pos, pos_idx]))
                rank_true = int(1 + np.sum(scores[batch_pos] > best_positive_score))
                expected = chance_lookup[len(pos_idx)]
                row = make_per_query_row(
                    direction="C2G",
                    representation=representation,
                    row_id=int(row_id),
                    gallery_definition_id=C2G_GALLERY_DEFINITION,
                    pos_definition_id=C2G_POS_DEFINITION,
                    gallery_ids_hash=gallery_payload.gallery_hash,
                    n_gallery=len(gallery_payload.gallery_ids),
                    m_pos=len(pos_idx),
                    rank_true=rank_true,
                    expected=expected,
                    target_raw=target_raw,
                    target_tokens_text=target_tokens_text,
                    query_n_targets=query_n_targets,
                    query_time=query_time,
                    query_dose_value=query_dose_value,
                    specificity_tier=specificity_tier,
                    treated_cell_id=treated_cell_id,
                    parsed_target_tokens=parsed_target_tokens,
                )
                acc.add_valid(
                    m_pos=len(pos_idx),
                    mrr_raw=row["mrr_raw"],
                    expected_mrr=row["expected_mrr_chance"],
                    mrr_corrected=row["mrr_corrected"],
                    hit_raw={1: row["hit1_raw"], 5: row["hit5_raw"], 10: row["hit10_raw"]},
                    hit_expected={
                        1: row["expected_hit1_chance"],
                        5: row["expected_hit5_chance"],
                        10: row["expected_hit10_chance"],
                    },
                    hit_corrected={1: row["hit1_corrected"], 5: row["hit5_corrected"], 10: row["hit10_corrected"]},
                )
                per_query_rows.append(row)
            if per_query_rows:
                sink.write(pd.DataFrame(per_query_rows))

    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="C2G",
            representation=representation,
            reason="m_pos0_after_intersection",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size - acc.n_mpos0),
            n_gallery_member_rows_before=int(gallery_member_after.size),
            n_gallery_member_rows_after=int(gallery_member_after.size),
            n_gallery_items_before=int(len(gallery_payload.gallery_ids)),
            n_gallery_items_after=int(len(gallery_payload.gallery_ids)),
            notes="Chemical query lost all positive genetic target centroids after valid intersection",
        ),
    )

    summary_row = build_summary_row(acc)
    chance_row = build_chance_row(summary_row)
    readiness = {
        "direction": "C2G",
        "representation": representation,
        "n_query_rows_before": int(query_before.size),
        "n_query_rows_after_valid_mask": int(query_after.size),
        "n_gallery_member_rows_before": int(gallery_member_before.size),
        "n_gallery_member_rows_after_valid_mask": int(gallery_member_after.size),
        "n_gallery_items_after_filter": int(len(gallery_payload.gallery_ids)),
        "n_mpos0": int(acc.n_mpos0),
        "n_valid": int(acc.n_valid),
        "gallery_ids": list(gallery_payload.gallery_ids),
    }
    return summary_row, chance_row, attrition_rows, readiness


def run_g2c_direction(
    *,
    arr: np.ndarray,
    representation: str,
    valid_mask: np.ndarray,
    gen_all_row_ids: np.ndarray,
    chem_all_row_ids: np.ndarray,
    all_context_items: Sequence[Tuple[str, np.ndarray]],
    target_order: Sequence[str],
    target_raw: Sequence[str],
    target_tokens_text: Sequence[str],
    query_n_targets: Sequence[int],
    query_time: np.ndarray,
    query_dose_value: np.ndarray,
    specificity_tier: Sequence[str],
    treated_cell_id: Sequence[str],
    parsed_target_tokens: Sequence[Tuple[str, ...]],
    sink: PerQueryParquetSink,
) -> Tuple[Dict[str, object], Dict[str, object], List[Dict[str, object]], Dict[str, object]]:
    attrition_rows: List[Dict[str, object]] = []

    query_before = gen_all_row_ids
    query_after = query_before[valid_mask[query_before]] if query_before.size > 0 else np.empty((0,), dtype=np.int64)

    gallery_member_before = chem_all_row_ids
    gallery_member_after = gallery_member_before[valid_mask[gallery_member_before]] if gallery_member_before.size > 0 else np.empty((0,), dtype=np.int64)

    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="G2C",
            representation=representation,
            reason="representation_valid_mask_drop_query",
            n_query_rows_before=int(query_before.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_member_before.size),
            n_gallery_member_rows_after=int(gallery_member_before.size),
            n_gallery_items_before=int(len(all_context_items)),
            n_gallery_items_after=int(len(all_context_items)),
            notes="Genetic query rows filtered by representation valid_mask before retrieval",
        ),
    )
    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="G2C",
            representation=representation,
            reason="representation_valid_mask_drop_gallery_member",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_member_before.size),
            n_gallery_member_rows_after=int(gallery_member_after.size),
            n_gallery_items_before=int(len(all_context_items)),
            n_gallery_items_after=int(len(all_context_items)),
            notes="Chemical context member rows filtered by representation valid_mask before centroid construction",
        ),
    )

    ordered_gallery_ids: List[str] = []
    ordered_member_rows: List[np.ndarray] = []
    token_to_context_indices: MutableMapping[str, List[int]] = {target: [] for target in target_order}

    for context_id, row_ids in all_context_items:
        valid_rows = row_ids[valid_mask[row_ids]] if row_ids.size > 0 else np.empty((0,), dtype=np.int64)
        if valid_rows.size == 0:
            continue
        ordered_gallery_ids.append(str(context_id))
        ordered_member_rows.append(valid_rows)
        token_union = {
            token
            for row_id in valid_rows.tolist()
            for token in parsed_target_tokens[row_id]
        }
        ordered_tokens = [target for target in target_order if target in token_union]
        gallery_idx = len(ordered_gallery_ids) - 1
        for token in ordered_tokens:
            token_to_context_indices[token].append(gallery_idx)

    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="G2C",
            representation=representation,
            reason="gallery_item_empty_after_filter",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_member_after.size),
            n_gallery_member_rows_after=int(gallery_member_after.size),
            n_gallery_items_before=int(len(all_context_items)),
            n_gallery_items_after=int(len(ordered_gallery_ids)),
            notes="Chemical contexts with zero valid members are removed from the gallery",
        ),
    )

    gallery_payload = GalleryPayload(
        gallery_ids=ordered_gallery_ids,
        member_row_ids=ordered_member_rows,
        centroids=np.empty((0, int(arr.shape[1])), dtype=np.float64),
        norms=np.empty((0,), dtype=np.float64),
        gallery_hash=hash_ordered_ids(ordered_gallery_ids),
    )
    if ordered_member_rows:
        centroids, norms = build_centroids(arr, ordered_member_rows)
        gallery_payload = GalleryPayload(
            gallery_ids=ordered_gallery_ids,
            member_row_ids=ordered_member_rows,
            centroids=centroids,
            norms=norms,
            gallery_hash=hash_ordered_ids(ordered_gallery_ids),
        )

    acc = DirectionAccumulator(
        dataset=EXPECTED_DATASET,
        cell_line=EXPECTED_CELL_LINE,
        direction="G2C",
        representation=representation,
        gallery_definition_id=G2C_GALLERY_DEFINITION,
        pos_definition_id=G2C_POS_DEFINITION,
        n_total=int(query_after.size),
        n_gallery=int(len(gallery_payload.gallery_ids)),
    )

    if gallery_payload.gallery_ids:
        chance_lookup = build_exact_chance_lookup(len(gallery_payload.gallery_ids))
        for start in range(0, query_after.size, QUERY_BATCH_ROWS):
            end = min(start + QUERY_BATCH_ROWS, query_after.size)
            row_batch = query_after[start:end]
            q_batch = np.asarray(arr[row_batch.astype(np.int64, copy=False)], dtype=np.float64)
            if not bool(np.isfinite(q_batch).all()):
                raise ValueError(f"{representation}/G2C query batch contains non-finite values after valid_mask filtering")
            q_norm = np.sum(q_batch * q_batch, axis=1, dtype=np.float64)
            scores = -(q_norm[:, None] + gallery_payload.norms[None, :] - 2.0 * (q_batch @ gallery_payload.centroids.T))

            per_query_rows: List[Dict[str, object]] = []
            for batch_pos, row_id in enumerate(row_batch.tolist()):
                query_token = parsed_target_tokens[row_id][0]
                pos_idx = token_to_context_indices.get(query_token, [])
                if not pos_idx:
                    acc.n_mpos0 += 1
                    continue
                best_positive_score = float(np.max(scores[batch_pos, pos_idx]))
                rank_true = int(1 + np.sum(scores[batch_pos] > best_positive_score))
                expected = chance_lookup[len(pos_idx)]
                row = make_per_query_row(
                    direction="G2C",
                    representation=representation,
                    row_id=int(row_id),
                    gallery_definition_id=G2C_GALLERY_DEFINITION,
                    pos_definition_id=G2C_POS_DEFINITION,
                    gallery_ids_hash=gallery_payload.gallery_hash,
                    n_gallery=len(gallery_payload.gallery_ids),
                    m_pos=len(pos_idx),
                    rank_true=rank_true,
                    expected=expected,
                    target_raw=target_raw,
                    target_tokens_text=target_tokens_text,
                    query_n_targets=query_n_targets,
                    query_time=query_time,
                    query_dose_value=query_dose_value,
                    specificity_tier=specificity_tier,
                    treated_cell_id=treated_cell_id,
                    parsed_target_tokens=parsed_target_tokens,
                )
                acc.add_valid(
                    m_pos=len(pos_idx),
                    mrr_raw=row["mrr_raw"],
                    expected_mrr=row["expected_mrr_chance"],
                    mrr_corrected=row["mrr_corrected"],
                    hit_raw={1: row["hit1_raw"], 5: row["hit5_raw"], 10: row["hit10_raw"]},
                    hit_expected={
                        1: row["expected_hit1_chance"],
                        5: row["expected_hit5_chance"],
                        10: row["expected_hit10_chance"],
                    },
                    hit_corrected={1: row["hit1_corrected"], 5: row["hit5_corrected"], 10: row["hit10_corrected"]},
                )
                per_query_rows.append(row)
            if per_query_rows:
                sink.write(pd.DataFrame(per_query_rows))

    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            direction="G2C",
            representation=representation,
            reason="m_pos0_after_intersection",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size - acc.n_mpos0),
            n_gallery_member_rows_before=int(gallery_member_after.size),
            n_gallery_member_rows_after=int(gallery_member_after.size),
            n_gallery_items_before=int(len(gallery_payload.gallery_ids)),
            n_gallery_items_after=int(len(gallery_payload.gallery_ids)),
            notes="Genetic query lost all positive chemical contexts after valid intersection",
        ),
    )

    summary_row = build_summary_row(acc)
    chance_row = build_chance_row(summary_row)
    readiness = {
        "direction": "G2C",
        "representation": representation,
        "n_query_rows_before": int(query_before.size),
        "n_query_rows_after_valid_mask": int(query_after.size),
        "n_gallery_member_rows_before": int(gallery_member_before.size),
        "n_gallery_member_rows_after_valid_mask": int(gallery_member_after.size),
        "n_gallery_items_after_filter": int(len(gallery_payload.gallery_ids)),
        "n_mpos0": int(acc.n_mpos0),
        "n_valid": int(acc.n_valid),
        "gallery_ids_count": int(len(gallery_payload.gallery_ids)),
    }
    return summary_row, chance_row, attrition_rows, readiness


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

    init_global_seed(seed)

    task2_snapshot = resolve_config_path(project_root, str(config["paths"]["task2_snapshot"]))
    runs_dir = resolve_config_path(project_root, str(config["paths"]["runs_dir"]))

    expected_snapshot = (project_root / EXPECTED_TASK2_SNAPSHOT).resolve()
    if task2_snapshot != expected_snapshot:
        print(
            "[ERROR] Data isolation violation: config.paths.task2_snapshot must resolve "
            f"to {expected_snapshot}, got {task2_snapshot}",
            file=sys.stderr,
        )
        return 4

    stage_dir = runs_dir / args.run_id / STAGE
    stage_dir.mkdir(parents=True, exist_ok=True)

    started_at = utc_now_iso()
    assertions: List[Dict[str, object]] = []
    input_paths: List[Path] = []

    assertions.append(
        {
            "name": "global_seed_locked",
            "pass": True,
            "details": {
                "rules": [
                    "GLOBAL_SEED must equal 619",
                    "HIT_K_LIST must equal [1,5,10]",
                    "QUERY_BATCH_ROWS is deterministic and fixed",
                ],
                "seed": GLOBAL_SEED,
                "hit_k_list": list(HIT_K_LIST),
                "query_batch_rows": QUERY_BATCH_ROWS,
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
                    "S5 reads exclusively from data/task2_snapshot_v1/",
                    "config.paths.task2_snapshot must resolve to data/task2_snapshot_v1",
                ],
                "task2_snapshot": str(task2_snapshot),
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "scoring_kernel_locked",
            "pass": True,
            "details": {
                "rules": [
                    "score(q,c) = -||q-c||^2",
                    "rank_true = 1 + count(score_all > best_positive_score)",
                ],
                "kernel": "negative_squared_euclidean",
                "tie_rule": "strict_greater_than",
            },
            "counterexamples": [],
        }
    )

    try:
        target_order = load_common_target_order(task2_snapshot / "k562/Common_Targets_K562.csv")
        delta_meta = load_delta_meta(task2_snapshot / "k562/derived/delta_meta.csv", target_order)
        drug_context_frame = load_drug_context_frame(task2_snapshot / "k562/Drug_meta.csv", delta_meta)
        chem_all_row_ids, gen_all_row_ids, gen_rows_by_target = build_genetic_rows_by_target(delta_meta, target_order)
        all_context_items = build_all_context_items(drug_context_frame)
    except Exception as exc:
        assertions.append(
            {
                "name": "task2_retrieval_inputs_ready",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed loading Task2 retrieval inputs: {exc}", file=sys.stderr)
        return 5

    input_paths.extend(
        [
            task2_snapshot / "k562/Common_Targets_K562.csv",
            task2_snapshot / "k562/derived/delta_meta.csv",
            task2_snapshot / "k562/Drug_meta.csv",
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
        return 6

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
                    "S5 consumes delta_meta.csv.target_tokens as the canonical normalized membership source",
                    "No runtime re-parsing from target_raw is used for positives",
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
    assertions.append(
        {
            "name": "drug_perturbation_raw_join_complete",
            "pass": True,
            "details": {
                "rules": [
                    "Every DRUG delta row must map to one Drug_meta perturbation_raw",
                    "G2C context key is (cell_line, perturbation_raw, time, dose_value)",
                ],
                "n_drug_rows": int(len(drug_context_frame)),
                "n_contexts": int(len(all_context_items)),
            },
            "counterexamples": [],
        }
    )

    target_raw = delta_meta["target_raw"].astype(str).tolist()
    target_tokens_text = delta_meta["target_tokens"].astype(str).tolist()
    query_n_targets = delta_meta["n_targets"].to_numpy(dtype=np.int64, copy=False)
    query_time = delta_meta["time"].to_numpy(dtype=np.float64, copy=False)
    query_dose_value = delta_meta["dose_value"].to_numpy(dtype=np.float64, copy=False)
    specificity_tier = delta_meta["specificity_tier"].astype(str).tolist()
    treated_cell_id = delta_meta["treated_cell_id"].astype(str).tolist()
    parsed_target_tokens = delta_meta["parsed_target_tokens"].tolist()

    per_query_path = stage_dir / "task2_retrieval_per_query.parquet"
    sink = PerQueryParquetSink(per_query_path, per_query_schema())

    summary_rows: List[Dict[str, object]] = []
    chance_rows: List[Dict[str, object]] = []
    attrition_rows: List[Dict[str, object]] = []
    readiness_rows: List[Dict[str, object]] = []
    representation_details: List[Dict[str, Any]] = []

    try:
        for spec in registry:
            payload, rep_inputs, rep_detail = validate_representation_inputs(
                snapshot_root=task2_snapshot,
                spec=spec,
                delta_meta=delta_meta,
            )
            representation_details.append(rep_detail)
            input_paths.extend(rep_inputs)

            arr = np.load(payload.array_path, mmap_mode="r")

            c2g_summary, c2g_chance, c2g_attrition, c2g_readiness = run_c2g_direction(
                arr=arr,
                representation=spec.name,
                valid_mask=payload.valid_mask,
                chem_all_row_ids=chem_all_row_ids,
                gen_all_row_ids=gen_all_row_ids,
                gen_rows_by_target=gen_rows_by_target,
                target_order=target_order,
                target_raw=target_raw,
                target_tokens_text=target_tokens_text,
                query_n_targets=query_n_targets,
                query_time=query_time,
                query_dose_value=query_dose_value,
                specificity_tier=specificity_tier,
                treated_cell_id=treated_cell_id,
                parsed_target_tokens=parsed_target_tokens,
                sink=sink,
            )
            summary_rows.append(c2g_summary)
            chance_rows.append(c2g_chance)
            attrition_rows.extend(c2g_attrition)
            readiness_rows.append(c2g_readiness)

            g2c_summary, g2c_chance, g2c_attrition, g2c_readiness = run_g2c_direction(
                arr=arr,
                representation=spec.name,
                valid_mask=payload.valid_mask,
                gen_all_row_ids=gen_all_row_ids,
                chem_all_row_ids=chem_all_row_ids,
                all_context_items=all_context_items,
                target_order=target_order,
                target_raw=target_raw,
                target_tokens_text=target_tokens_text,
                query_n_targets=query_n_targets,
                query_time=query_time,
                query_dose_value=query_dose_value,
                specificity_tier=specificity_tier,
                treated_cell_id=treated_cell_id,
                parsed_target_tokens=parsed_target_tokens,
                sink=sink,
            )
            summary_rows.append(g2c_summary)
            chance_rows.append(g2c_chance)
            attrition_rows.extend(g2c_attrition)
            readiness_rows.append(g2c_readiness)

            del arr
    except Exception as exc:
        sink.close()
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
        return 7

    sink.close()

    summary_df = sort_summary_frame(pd.DataFrame(summary_rows))
    chance_df = sort_chance_frame(pd.DataFrame(chance_rows))
    attrition_df = sort_attrition_frame(pd.DataFrame(attrition_rows))
    summary_long_df = sort_summary_long_frame(pd.DataFrame(build_summary_long_rows(summary_rows)))

    summary_path = stage_dir / "task2_retrieval_summary.csv"
    chance_path = stage_dir / "task2_chance_identity_check.csv"
    attrition_path = stage_dir / "task2_retrieval_attrition.csv"
    summary_long_path = stage_dir / "task2_retrieval_summary_long.csv"

    write_csv(summary_df, summary_path)
    write_csv(chance_df, chance_path)
    write_csv(attrition_df, attrition_path)
    write_csv(summary_long_df, summary_long_path)

    direction_nonempty_pass = all(
        int(row["n_query_rows_after_valid_mask"]) > 0 and int(row["n_gallery_items_after_filter"]) > 0
        for row in readiness_rows
    )
    assertions.append(
        {
            "name": "current_snapshot_direction_nonempty_observation",
            "pass": bool(direction_nonempty_pass),
            "details": {
                "rules": [
                    "Current local snapshot is expected to leave both directions non-empty for every representation",
                    "This is a current-snapshot observation, not a universal fail-fast rule",
                ],
                "rows": readiness_rows,
            },
            "counterexamples": [
                row
                for row in readiness_rows
                if int(row["n_query_rows_after_valid_mask"]) == 0 or int(row["n_gallery_items_after_filter"]) == 0
            ][:MAX_COUNTEREXAMPLES],
        }
    )

    mpos0_pass = all(int(row["n_mpos0"]) == 0 for row in readiness_rows)
    assertions.append(
        {
            "name": "current_snapshot_no_mpos0_observation",
            "pass": bool(mpos0_pass),
            "details": {
                "rules": [
                    "Current local readiness audit expects every valid query to have at least one positive opposite-side gallery item",
                    "This is an observation and the implementation still soft-fails m_pos0 if drift occurs",
                ],
                "rows": readiness_rows,
            },
            "counterexamples": [row for row in readiness_rows if int(row["n_mpos0"]) > 0][:MAX_COUNTEREXAMPLES],
        }
    )

    shape_contract_pass = all(
        int(detail["n_rows"]) == int(len(delta_meta))
        and int(detail["dim"]) == int(next(spec.expected_dim for spec in registry if spec.name == detail["representation"]))
        for detail in representation_details
    )
    assertions.append(
        {
            "name": "representation_shape_and_finite_mask_contract",
            "pass": bool(shape_contract_pass),
            "details": {"representations": representation_details},
            "counterexamples": [] if shape_contract_pass else representation_details[:MAX_COUNTEREXAMPLES],
        }
    )

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
                    "For each direction/representation: n_total == n_valid + n_excluded_missing_metric_or_mpos0",
                ],
                "n_checked": int(len(summary_df)),
            },
            "counterexamples": denominator_failures.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
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
                    "For each finite summary row: |delta_mrr| <= 1e-12",
                    "For each finite summary row: |delta_hit1|, |delta_hit5|, |delta_hit10| <= 1e-12",
                ],
                "tolerance": CHANCE_IDENTITY_TOL,
                "n_checked": int(len(chance_finite)),
            },
            "counterexamples": chance_violations.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    core_outputs = [
        per_query_path,
        summary_path,
        attrition_path,
        chance_path,
        summary_long_path,
    ]
    output_routing_pass = all(path.resolve().is_relative_to(stage_dir.resolve()) for path in core_outputs)
    assertions.append(
        {
            "name": "output_routing_isolated",
            "pass": bool(output_routing_pass),
            "details": {
                "rules": [f"All S5 outputs must be written under {stage_dir}"],
                "stage_dir": str(stage_dir),
            },
            "counterexamples": []
            if output_routing_pass
            else [{"bad_path": str(path)} for path in core_outputs if not path.resolve().is_relative_to(stage_dir.resolve())],
        }
    )

    bad_inputs = [
        str(path.resolve())
        for path in sorted(set(input_paths), key=lambda p: str(p))
        if not path.resolve().is_relative_to(task2_snapshot.resolve())
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "rules": ["All loaded S5 inputs must be under data/task2_snapshot_v1 lexical namespace"],
                "task2_snapshot": str(task2_snapshot.resolve()),
                "n_inputs": int(len(set(input_paths))),
            },
            "counterexamples": [{"bad_input": path} for path in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    completed_at = utc_now_iso()
    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    manifest_path = stage_dir / "manifest.json"

    output_paths = [
        str(path.resolve())
        for path in [
            per_query_path,
            summary_path,
            attrition_path,
            chance_path,
            summary_long_path,
            run_manifest_path,
            audit_assertions_path,
            manifest_path,
        ]
    ]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s5_task2_retrieval.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "hit_k_list": list(HIT_K_LIST),
            "chance_identity_tol": CHANCE_IDENTITY_TOL,
            "query_batch_rows": QUERY_BATCH_ROWS,
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
            "representation_order": list(REPRESENTATION_ORDER),
            "direction_order": list(DIRECTION_ORDER),
            "canonical_target_order": target_order,
            "target_membership_source": "delta_meta.target_tokens",
            "scoring_kernel": "negative_squared_euclidean",
            "tie_rule": "strict_greater_than",
            "gallery_definition_ids": {
                "C2G": C2G_GALLERY_DEFINITION,
                "G2C": G2C_GALLERY_DEFINITION,
            },
            "pos_definition_ids": {
                "C2G": C2G_POS_DEFINITION,
                "G2C": G2C_POS_DEFINITION,
            },
            "rng_draws_main_process": 0,
            "s3_provenance_note": "Known .absolute() vs .resolve() inconsistency remains non-blocking for S5",
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths), key=lambda p: str(p))],
        "outputs": output_paths,
        "summary": {
            "n_row_universe": int(len(delta_meta)),
            "n_targets": int(len(target_order)),
            "n_representations": int(len(REPRESENTATION_ORDER)),
            "n_summary_rows": int(len(summary_df)),
            "n_chance_rows": int(len(chance_df)),
            "n_attrition_rows": int(len(attrition_df)),
            "n_summary_long_rows": int(len(summary_long_df)),
            "n_per_query_rows": int(summary_df["n_valid"].sum()) if not summary_df.empty else 0,
            "representation_details": representation_details,
            "readiness_rows": readiness_rows,
        },
    }

    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})

    manifest_entries: List[Dict[str, object]] = []
    for file_path in [
        per_query_path,
        summary_path,
        attrition_path,
        chance_path,
        summary_long_path,
        run_manifest_path,
        audit_assertions_path,
    ]:
        manifest_entries.append(
            {
                "relative_path": file_path.name,
                "size_bytes": int(file_path.stat().st_size),
                "sha256": compute_sha256(file_path),
            }
        )
    write_json(manifest_path, {"stage": STAGE, "files": manifest_entries})

    if not chance_pass:
        print(
            f"[ERROR] chance identity check failed: tolerance={CHANCE_IDENTITY_TOL}",
            file=sys.stderr,
        )
        return 8

    return 0


if __name__ == "__main__":
    sys.exit(main())
