#!/usr/bin/env python3
# SCRIPT_HEADER_CONTRACT
# Script: scripts/s5_task2_retrieval_multisource.py
# Purpose: Compute corrected Task2 multisource target-level retrieval metrics
#   under the approved S5 contract from data/task2_snapshot_v2 only.
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
#   - runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet
#   - runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_summary.csv
#   - runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_summary_long.csv
#   - runs/<run_id>/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv
#   - runs/<run_id>/s5_task2_retrieval_multisource/task2_chance_identity_check.csv
#   - runs/<run_id>/s5_task2_retrieval_multisource/run_manifest.json
#   - runs/<run_id>/s5_task2_retrieval_multisource/audit_assertions.json
#   - runs/<run_id>/s5_task2_retrieval_multisource/manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s5_task2_retrieval_multisource/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s5_task2_retrieval_multisource.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Snapshot root is not data/task2_snapshot_v2 -> exit non-zero
#   - Missing required snapshot files -> exit non-zero
#   - Coverage / registry / schema / key / scope drift -> exit non-zero
#   - Representation validation failure or forbidden 1HAE analysis label -> exit non-zero
#   - Summary-grid / denominator / chance-identity contract violation -> exit non-zero
# Last Updated: 2026-03-11

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from s4_task2_group_concordance_multisource import (
    DATASET_ORDER,
    MAX_COUNTEREXAMPLES,
    REPRESENTATION_ORDER,
    STATUS_AVAILABLE,
    RepresentationPayload,
    RepresentationSpec,
    assert_unique_key,
    build_dataset_bundles,
    build_physical_representation_specs,
    build_stage_manifest,
    ensure_required_columns,
    fail_on_forbidden_cell_line,
    load_common_targets,
    load_delta_meta,
    load_membership_parity_table,
    load_shared_contract_tables,
    load_snapshot_manifest,
    normalize_text_series,
    resolve_config_path,
    safe_git_head,
    validate_membership_parity,
    validate_representation_inputs,
    write_csv,
    write_json,
    build_membership_from_delta_meta,
)

STAGE = "s5_task2_retrieval_multisource"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v2")

GLOBAL_SEED = 619
HIT_K_LIST: Tuple[int, ...] = (1, 5, 10)
CHANCE_IDENTITY_TOL = 1e-12
QUERY_BATCH_ROWS = 1024

DIRECTION_ORDER: Tuple[str, str] = ("C2G", "G2C")
SUMMARY_LONG_METRICS: Tuple[str, ...] = (
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
)

C2G_GALLERY_DEFINITION = "C2G_genetic_target_centroid_gallery_v1"
G2C_GALLERY_DEFINITION = "G2C_chemical_target_centroid_gallery_v1"
C2G_POS_DEFINITION = "C2G_tokens_intersect_genetic_target_gallery_v1"
G2C_POS_DEFINITION = "G2C_atomic_token_matches_chemical_target_gallery_v1"

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


@dataclass(frozen=True)
class GalleryPayload:
    gallery_ids: List[str]
    member_row_ids: List[np.ndarray]
    centroids: np.ndarray
    norms: np.ndarray
    gallery_hash: str
    source_target_count: int
    gallery_member_rows_before: np.ndarray
    gallery_member_rows_after: np.ndarray


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
    parser = argparse.ArgumentParser(description="S5 corrected Task2 multisource target-level retrieval")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    return parser.parse_args()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def init_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def hash_ordered_ids(ids: Sequence[str]) -> str:
    material = "||".join(ids)
    return hashlib.sha256(material.encode("utf-8")).hexdigest()


def unique_union(arrays: Sequence[np.ndarray]) -> np.ndarray:
    non_empty = [array.astype(np.int64, copy=False) for array in arrays if array.size > 0]
    if not non_empty:
        return np.empty((0,), dtype=np.int64)
    return np.unique(np.concatenate(non_empty))


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


@lru_cache(maxsize=None)
def exact_chance_metrics(n_gallery: int, m_pos: int) -> Tuple[float, float, float, float]:
    if n_gallery <= 0:
        raise ValueError("n_gallery must be >=1 for chance correction")
    if m_pos < 1 or m_pos > n_gallery:
        raise ValueError(f"m_pos must satisfy 1 <= m_pos <= n_gallery, got m_pos={m_pos}, n_gallery={n_gallery}")

    total = math.comb(n_gallery, m_pos)
    max_rank = n_gallery - m_pos + 1
    expected_mrr = 0.0
    expected_hit = {k: 0.0 for k in HIT_K_LIST}
    for rank in range(1, max_rank + 1):
        prob = math.comb(n_gallery - rank, m_pos - 1) / total
        expected_mrr += prob / float(rank)
        for k in HIT_K_LIST:
            if rank <= k:
                expected_hit[k] += prob
    return (
        float(expected_mrr),
        float(expected_hit[1]),
        float(expected_hit[5]),
        float(expected_hit[10]),
    )


def build_centroids(arr: np.ndarray, ordered_groups: Sequence[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
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


def build_cell_line_order_map(snapshot_manifest: Mapping[str, Any]) -> Mapping[str, Mapping[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    for dataset in DATASET_ORDER:
        cell_lines = snapshot_manifest["datasets"][dataset]["cell_lines"]
        out[dataset] = {cell_line: idx for idx, cell_line in enumerate(cell_lines)}
    return out


def load_available_rep_scope(
    snapshot_root: Path,
    eligible_coverage: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[Path]]:
    registry_path = snapshot_root / "representation_availability_registry.csv"
    if not registry_path.is_file():
        raise FileNotFoundError(f"Missing representation_availability_registry.csv: {registry_path}")

    registry = pd.read_csv(registry_path)
    ensure_required_columns(
        registry,
        ["dataset", "cell_line", "representation", "availability_status", "availability_reason"],
        "representation_availability_registry.csv",
    )
    registry["dataset"] = normalize_text_series(registry["dataset"])
    registry["cell_line"] = normalize_text_series(registry["cell_line"])
    registry["representation"] = normalize_text_series(registry["representation"])
    registry["availability_status"] = normalize_text_series(registry["availability_status"])
    registry["availability_reason"] = normalize_text_series(registry["availability_reason"])
    fail_on_forbidden_cell_line(registry, ["cell_line"], "representation_availability_registry.csv")
    assert_unique_key(registry, ["dataset", "cell_line", "representation"], "representation_availability_registry.csv")

    eligible_pairs = eligible_coverage[["dataset", "cell_line"]].drop_duplicates().reset_index(drop=True)
    scoped_registry = eligible_pairs.merge(
        registry,
        on=["dataset", "cell_line"],
        how="left",
        sort=False,
        indicator=True,
    )
    missing_pairs = scoped_registry.loc[scoped_registry["_merge"].ne("both"), ["dataset", "cell_line"]]
    if not missing_pairs.empty:
        examples = missing_pairs.drop_duplicates().head(MAX_COUNTEREXAMPLES).to_dict(orient="records")
        raise ValueError(f"Missing availability rows for eligible dataset/cell_line pairs: {examples}")
    scoped_registry = scoped_registry.drop(columns=["_merge"])

    available = scoped_registry.loc[
        scoped_registry["availability_status"].eq(STATUS_AVAILABLE),
        ["dataset", "cell_line", "representation"],
    ].copy()
    assert_unique_key(available, ["dataset", "cell_line", "representation"], "available representation scope")
    return registry, available.reset_index(drop=True), [registry_path]


def build_source_target_order_by_slice(
    snapshot_root: Path,
    eligible_coverage: pd.DataFrame,
) -> Tuple[Mapping[Tuple[str, str], List[str]], List[Path]]:
    source_order: Dict[Tuple[str, str], List[str]] = {}
    input_paths: List[Path] = []

    lincs_eligible = eligible_coverage.loc[eligible_coverage["dataset"].eq("LINCS")].copy().reset_index(drop=True)
    assert_unique_key(lincs_eligible, ["dataset", "cell_line", "target_token"], "eligible LINCS coverage")
    for cell_line in lincs_eligible["cell_line"].drop_duplicates().tolist():
        ordered_targets = lincs_eligible.loc[lincs_eligible["cell_line"].eq(cell_line), "target_token"].tolist()
        if not ordered_targets:
            raise ValueError(f"LINCS eligible target order is empty for cell_line={cell_line}")
        if len(ordered_targets) != len(set(ordered_targets)):
            raise ValueError(f"LINCS eligible target order contains duplicates for cell_line={cell_line}")
        source_order[("LINCS", cell_line)] = ordered_targets

    common_targets_path = snapshot_root / "scperturb_k562/Common_Targets_K562.csv"
    common_targets, common_target_inputs = load_common_targets(common_targets_path)
    input_paths.extend(common_target_inputs)
    sc_eligible = eligible_coverage.loc[eligible_coverage["dataset"].eq("scPerturb")].copy().reset_index(drop=True)
    if sc_eligible.empty:
        raise ValueError("Eligible scPerturb coverage cannot be empty")
    sc_cell_lines = sc_eligible["cell_line"].drop_duplicates().tolist()
    if sc_cell_lines != ["K562"]:
        raise ValueError(f"scPerturb eligible coverage must contain only K562, got={sc_cell_lines}")
    eligible_target_set = set(sc_eligible["target_token"].tolist())
    ordered_targets = [target for target in common_targets if target in eligible_target_set]
    if set(ordered_targets) != eligible_target_set:
        raise ValueError("scPerturb Common_Targets_K562.csv order does not cover eligible targets exactly")
    if len(ordered_targets) != len(set(ordered_targets)):
        raise ValueError("scPerturb Common_Targets_K562.csv filtered order contains duplicates")
    source_order[("scPerturb", "K562")] = ordered_targets
    return source_order, input_paths


def build_rows_by_slice_and_target(
    *,
    dataset: str,
    membership: pd.DataFrame,
    eligible_coverage_subset: pd.DataFrame,
    source_target_order_by_slice: Mapping[Tuple[str, str], List[str]],
) -> Tuple[
    Mapping[str, np.ndarray],
    Mapping[str, np.ndarray],
    Mapping[Tuple[str, str], np.ndarray],
    Mapping[Tuple[str, str], np.ndarray],
]:
    chem_query_rows: Dict[str, np.ndarray] = {}
    gen_query_rows: Dict[str, np.ndarray] = {}
    chem_target_rows: Dict[Tuple[str, str], np.ndarray] = {}
    gen_target_rows: Dict[Tuple[str, str], np.ndarray] = {}

    rows_by_side_target: MutableMapping[Tuple[str, str, str], List[int]] = {}
    for (dataset_name, cell_line), target_order in source_target_order_by_slice.items():
        if dataset_name != dataset:
            continue
        target_set = set(target_order)
        for target in target_order:
            rows_by_side_target[(cell_line, target, "Chemical")] = []
            rows_by_side_target[(cell_line, target, "Genetic")] = []
        subset = membership.loc[membership["cell_line"].eq(cell_line) & membership["target_token"].isin(target_set)].copy()
        for row in subset.itertuples(index=False):
            rows_by_side_target[(row.cell_line, row.target_token, row.perturbation_class)].append(int(row.row_id))

        chem_arrays: List[np.ndarray] = []
        gen_arrays: List[np.ndarray] = []
        for target in target_order:
            chem_key = (cell_line, target, "Chemical")
            gen_key = (cell_line, target, "Genetic")
            chem_arr = np.unique(np.asarray(rows_by_side_target[chem_key], dtype=np.int64))
            gen_arr = np.unique(np.asarray(rows_by_side_target[gen_key], dtype=np.int64))
            chem_target_rows[(cell_line, target)] = chem_arr
            gen_target_rows[(cell_line, target)] = gen_arr
            if chem_arr.size > 0:
                chem_arrays.append(chem_arr)
            if gen_arr.size > 0:
                gen_arrays.append(gen_arr)

        chem_query_rows[cell_line] = unique_union(chem_arrays)
        gen_query_rows[cell_line] = unique_union(gen_arrays)

    coverage_keyed = eligible_coverage_subset.set_index(["cell_line", "target_token"])
    for (dataset_name, cell_line), target_order in source_target_order_by_slice.items():
        if dataset_name != dataset:
            continue
        for target in target_order:
            key = (cell_line, target)
            expected_chem = int(coverage_keyed.loc[key, "n_chem_instances"])
            expected_gen = int(coverage_keyed.loc[key, "n_gen_instances"])
            actual_chem = int(chem_target_rows[key].size)
            actual_gen = int(gen_target_rows[key].size)
            if actual_chem != expected_chem or actual_gen != expected_gen:
                raise ValueError(
                    f"{dataset} membership count drift for cell_line={cell_line}, target_token={target}: "
                    f"expected=({expected_chem},{expected_gen}) got=({actual_chem},{actual_gen})"
                )

    return chem_query_rows, gen_query_rows, chem_target_rows, gen_target_rows


def build_target_gallery(
    *,
    arr: np.ndarray,
    target_order: Sequence[str],
    target_rows_by_key: Mapping[Tuple[str, str], np.ndarray],
    cell_line: str,
    valid_mask: np.ndarray,
) -> GalleryPayload:
    member_rows_before_list: List[np.ndarray] = []
    member_rows_after_list: List[np.ndarray] = []
    ordered_gallery_ids: List[str] = []
    ordered_member_rows: List[np.ndarray] = []

    for target in target_order:
        rows_before = target_rows_by_key[(cell_line, target)]
        member_rows_before_list.append(rows_before)
        valid_rows = rows_before[valid_mask[rows_before]] if rows_before.size > 0 else np.empty((0,), dtype=np.int64)
        if valid_rows.size > 0:
            ordered_gallery_ids.append(target)
            ordered_member_rows.append(valid_rows)
            member_rows_after_list.append(valid_rows)

    gallery_member_rows_before = unique_union(member_rows_before_list)
    gallery_member_rows_after = unique_union(member_rows_after_list)
    centroids, norms = build_centroids(arr, ordered_member_rows)
    return GalleryPayload(
        gallery_ids=ordered_gallery_ids,
        member_row_ids=ordered_member_rows,
        centroids=centroids,
        norms=norms,
        gallery_hash=hash_ordered_ids(ordered_gallery_ids),
        source_target_count=int(len(target_order)),
        gallery_member_rows_before=gallery_member_rows_before,
        gallery_member_rows_after=gallery_member_rows_after,
    )


def make_attrition_row(
    *,
    dataset: str,
    cell_line: str,
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
        "dataset": dataset,
        "cell_line": cell_line,
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


def maybe_append_attrition(rows: List[Dict[str, object]], payload: Mapping[str, object]) -> None:
    if (
        int(payload["n_query_rows_removed"]) <= 0
        and int(payload["n_gallery_member_rows_removed"]) <= 0
        and int(payload["n_gallery_items_removed"]) <= 0
    ):
        return
    rows.append(dict(payload))


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


def make_per_query_row(
    *,
    dataset: str,
    cell_line: str,
    direction: str,
    representation: str,
    row_id: int,
    gallery_definition_id: str,
    pos_definition_id: str,
    gallery_ids_hash: str,
    n_gallery: int,
    m_pos: int,
    rank_true: int,
    expected: Tuple[float, float, float, float],
    target_raw: Sequence[str],
    target_tokens_text: Sequence[str],
    query_n_targets: np.ndarray,
    query_time: np.ndarray,
    query_dose_value: np.ndarray,
    treated_cell_id: Sequence[str],
    query_target_token: str,
) -> Dict[str, object]:
    expected_mrr, expected_hit1, expected_hit5, expected_hit10 = expected
    mrr_raw = 1.0 / float(rank_true)
    hit_raw = {k: float(rank_true <= k) for k in HIT_K_LIST}
    hit_expected = {1: expected_hit1, 5: expected_hit5, 10: expected_hit10}
    hit_corrected = {k: float(hit_raw[k] - hit_expected[k]) for k in HIT_K_LIST}
    mrr_corrected = float(mrr_raw - expected_mrr)

    return {
        "dataset": dataset,
        "cell_line": cell_line,
        "direction": direction,
        "representation": representation,
        "query_row_id": int(row_id),
        "query_uid": f"{dataset}||{cell_line}||{row_id}",
        "treated_cell_id": str(treated_cell_id[row_id]),
        "query_target_raw": str(target_raw[row_id]),
        "query_target_tokens": str(target_tokens_text[row_id]),
        "query_target_token": str(query_target_token),
        "query_n_targets": int(query_n_targets[row_id]),
        "query_time": float(query_time[row_id]),
        "query_dose_value": float(query_dose_value[row_id]),
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
        "expected_mrr_chance": float(expected_mrr),
        "expected_hit1_chance": float(expected_hit1),
        "expected_hit5_chance": float(expected_hit5),
        "expected_hit10_chance": float(expected_hit10),
        "mrr_corrected": float(mrr_corrected),
        "hit1_corrected": float(hit_corrected[1]),
        "hit5_corrected": float(hit_corrected[5]),
        "hit10_corrected": float(hit_corrected[10]),
    }


def build_order_maps(
    cell_line_order_map: Mapping[str, Mapping[str, int]],
) -> Tuple[Mapping[str, int], Mapping[Tuple[str, str], int], Mapping[str, int], Mapping[str, int], Mapping[str, int]]:
    dataset_order = {name: idx for idx, name in enumerate(DATASET_ORDER)}
    cell_line_order = {
        (dataset, cell_line): idx
        for dataset, mapping in cell_line_order_map.items()
        for cell_line, idx in mapping.items()
    }
    direction_order = {name: idx for idx, name in enumerate(DIRECTION_ORDER)}
    representation_order = {name: idx for idx, name in enumerate(REPRESENTATION_ORDER)}
    metric_order = {name: idx for idx, name in enumerate(SUMMARY_LONG_METRICS)}
    return dataset_order, cell_line_order, direction_order, representation_order, metric_order


def sort_summary_frame(frame: pd.DataFrame, cell_line_order_map: Mapping[str, Mapping[str, int]]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)
    dataset_order, cell_line_order, direction_order, representation_order, _ = build_order_maps(cell_line_order_map)
    out = frame.copy()
    out["_dataset_order"] = out["dataset"].map(dataset_order)
    out["_cell_line_order"] = [
        cell_line_order[(dataset, cell_line)]
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_direction_order"] = out["direction"].map(direction_order)
    out["_representation_order"] = out["representation"].map(representation_order)
    out = out.sort_values(
        ["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order"],
        kind="mergesort",
    )
    out = out.drop(columns=["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order"])
    return out[SUMMARY_COLUMNS].reset_index(drop=True)


def sort_chance_frame(frame: pd.DataFrame, cell_line_order_map: Mapping[str, Mapping[str, int]]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=CHANCE_COLUMNS)
    dataset_order, cell_line_order, direction_order, representation_order, _ = build_order_maps(cell_line_order_map)
    out = frame.copy()
    out["_dataset_order"] = out["dataset"].map(dataset_order)
    out["_cell_line_order"] = [
        cell_line_order[(dataset, cell_line)]
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_direction_order"] = out["direction"].map(direction_order)
    out["_representation_order"] = out["representation"].map(representation_order)
    out = out.sort_values(
        ["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order"],
        kind="mergesort",
    )
    out = out.drop(columns=["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order"])
    return out[CHANCE_COLUMNS].reset_index(drop=True)


def sort_attrition_frame(frame: pd.DataFrame, cell_line_order_map: Mapping[str, Mapping[str, int]]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=ATTRITION_COLUMNS)
    dataset_order, cell_line_order, direction_order, representation_order, _ = build_order_maps(cell_line_order_map)
    out = frame.copy()
    out["_dataset_order"] = out["dataset"].map(dataset_order)
    out["_cell_line_order"] = [
        cell_line_order[(dataset, cell_line)]
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_direction_order"] = out["direction"].map(direction_order)
    out["_representation_order"] = out["representation"].map(representation_order)
    out = out.sort_values(
        ["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order", "reason"],
        kind="mergesort",
    )
    out = out.drop(columns=["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order"])
    return out[ATTRITION_COLUMNS].reset_index(drop=True)


def build_summary_long_rows(summary_rows: Sequence[Mapping[str, object]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for row in summary_rows:
        for metric_name in SUMMARY_LONG_METRICS:
            rows.append(
                {
                    "dataset": row["dataset"],
                    "cell_line": row["cell_line"],
                    "direction": row["direction"],
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


def sort_summary_long_frame(frame: pd.DataFrame, cell_line_order_map: Mapping[str, Mapping[str, int]]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=SUMMARY_LONG_COLUMNS)
    dataset_order, cell_line_order, direction_order, representation_order, metric_order = build_order_maps(cell_line_order_map)
    out = frame.copy()
    out["_dataset_order"] = out["dataset"].map(dataset_order)
    out["_cell_line_order"] = [
        cell_line_order[(dataset, cell_line)]
        for dataset, cell_line in zip(out["dataset"].tolist(), out["cell_line"].tolist())
    ]
    out["_direction_order"] = out["direction"].map(direction_order)
    out["_representation_order"] = out["representation"].map(representation_order)
    out["_metric_order"] = out["metric_name"].map(metric_order)
    out = out.sort_values(
        ["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order", "_metric_order"],
        kind="mergesort",
    )
    out = out.drop(
        columns=["_dataset_order", "_cell_line_order", "_direction_order", "_representation_order", "_metric_order"]
    )
    return out[SUMMARY_LONG_COLUMNS].reset_index(drop=True)


def build_expected_summary_grid(available_rep_scope: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for row in available_rep_scope.itertuples(index=False):
        for direction in DIRECTION_ORDER:
            rows.append(
                {
                    "dataset": row.dataset,
                    "cell_line": row.cell_line,
                    "direction": direction,
                    "representation": row.representation,
                }
            )
    expected = pd.DataFrame(rows)
    assert_unique_key(expected, ["dataset", "cell_line", "direction", "representation"], "expected summary grid")
    return expected


def run_direction_for_slice(
    *,
    arr: np.ndarray,
    dataset: str,
    cell_line: str,
    direction: str,
    representation: str,
    valid_mask: np.ndarray,
    target_order: Sequence[str],
    query_row_ids: np.ndarray,
    gallery_target_rows_by_key: Mapping[Tuple[str, str], np.ndarray],
    target_raw: Sequence[str],
    target_tokens_text: Sequence[str],
    query_n_targets: np.ndarray,
    query_time: np.ndarray,
    query_dose_value: np.ndarray,
    treated_cell_id: Sequence[str],
    parsed_target_tokens: Sequence[Tuple[str, ...]],
    sink: PerQueryParquetSink,
) -> Tuple[Dict[str, object], Dict[str, object], List[Dict[str, object]], Dict[str, object]]:
    if direction not in DIRECTION_ORDER:
        raise ValueError(f"Unsupported direction: {direction}")

    gallery_definition_id = C2G_GALLERY_DEFINITION if direction == "C2G" else G2C_GALLERY_DEFINITION
    pos_definition_id = C2G_POS_DEFINITION if direction == "C2G" else G2C_POS_DEFINITION

    query_before = query_row_ids.astype(np.int64, copy=False)
    query_after = query_before[valid_mask[query_before]] if query_before.size > 0 else np.empty((0,), dtype=np.int64)
    gallery_payload = build_target_gallery(
        arr=arr,
        target_order=target_order,
        target_rows_by_key=gallery_target_rows_by_key,
        cell_line=cell_line,
        valid_mask=valid_mask,
    )

    attrition_rows: List[Dict[str, object]] = []
    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            dataset=dataset,
            cell_line=cell_line,
            direction=direction,
            representation=representation,
            reason="representation_valid_mask_drop_query",
            n_query_rows_before=int(query_before.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_payload.gallery_member_rows_before.size),
            n_gallery_member_rows_after=int(gallery_payload.gallery_member_rows_before.size),
            n_gallery_items_before=int(gallery_payload.source_target_count),
            n_gallery_items_after=int(gallery_payload.source_target_count),
            notes="Query rows filtered by representation valid_mask before retrieval",
        ),
    )
    maybe_append_attrition(
        attrition_rows,
        make_attrition_row(
            dataset=dataset,
            cell_line=cell_line,
            direction=direction,
            representation=representation,
            reason="representation_valid_mask_drop_gallery_member",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size),
            n_gallery_member_rows_before=int(gallery_payload.gallery_member_rows_before.size),
            n_gallery_member_rows_after=int(gallery_payload.gallery_member_rows_after.size),
            n_gallery_items_before=int(gallery_payload.source_target_count),
            n_gallery_items_after=int(gallery_payload.source_target_count),
            notes="Gallery member rows filtered by representation valid_mask before centroid construction",
        ),
    )

    acc = DirectionAccumulator(
        dataset=dataset,
        cell_line=cell_line,
        direction=direction,
        representation=representation,
        gallery_definition_id=gallery_definition_id,
        pos_definition_id=pos_definition_id,
        n_total=int(query_after.size),
        n_gallery=int(len(gallery_payload.gallery_ids)),
    )

    if not gallery_payload.gallery_ids:
        maybe_append_attrition(
            attrition_rows,
            make_attrition_row(
                dataset=dataset,
                cell_line=cell_line,
                direction=direction,
                representation=representation,
                reason="gallery_target_empty_after_filter",
                n_query_rows_before=int(query_after.size),
                n_query_rows_after=int(query_after.size),
                n_gallery_member_rows_before=int(gallery_payload.gallery_member_rows_before.size),
                n_gallery_member_rows_after=int(gallery_payload.gallery_member_rows_after.size),
                n_gallery_items_before=int(gallery_payload.source_target_count),
                n_gallery_items_after=0,
                notes="No valid gallery target centroid remains after validity filtering",
            ),
        )
        summary_row = build_summary_row(acc)
        chance_row = build_chance_row(summary_row)
        readiness = {
            "dataset": dataset,
            "cell_line": cell_line,
            "direction": direction,
            "representation": representation,
            "n_query_rows_before": int(query_before.size),
            "n_query_rows_after_valid_mask": int(query_after.size),
            "n_gallery_member_rows_before": int(gallery_payload.gallery_member_rows_before.size),
            "n_gallery_member_rows_after_valid_mask": int(gallery_payload.gallery_member_rows_after.size),
            "n_gallery_items_after_filter": 0,
            "n_mpos0": 0,
            "n_valid": 0,
        }
        return summary_row, chance_row, attrition_rows, readiness

    target_to_gallery_idx = {target: idx for idx, target in enumerate(gallery_payload.gallery_ids)}
    for start in range(0, query_after.size, QUERY_BATCH_ROWS):
        end = min(start + QUERY_BATCH_ROWS, query_after.size)
        row_batch = query_after[start:end]
        q_batch = np.asarray(arr[row_batch.astype(np.int64, copy=False)], dtype=np.float64)
        if not bool(np.isfinite(q_batch).all()):
            raise ValueError(
                f"{dataset}/{cell_line}/{representation}/{direction} query batch contains non-finite values after valid_mask"
            )
        q_norm = np.sum(q_batch * q_batch, axis=1, dtype=np.float64)
        scores = -(q_norm[:, None] + gallery_payload.norms[None, :] - 2.0 * (q_batch @ gallery_payload.centroids.T))

        per_query_rows: List[Dict[str, object]] = []
        for batch_pos, row_id in enumerate(row_batch.tolist()):
            if direction == "C2G":
                pos_idx = [target_to_gallery_idx[token] for token in parsed_target_tokens[row_id] if token in target_to_gallery_idx]
                query_target_token = "NA"
            else:
                query_target_token = parsed_target_tokens[row_id][0]
                pos_idx = [target_to_gallery_idx[query_target_token]] if query_target_token in target_to_gallery_idx else []
            if not pos_idx:
                acc.n_mpos0 += 1
                continue
            best_positive_score = float(np.max(scores[batch_pos, pos_idx]))
            rank_true = int(1 + np.sum(scores[batch_pos] > best_positive_score))
            expected = exact_chance_metrics(len(gallery_payload.gallery_ids), len(pos_idx))
            row = make_per_query_row(
                dataset=dataset,
                cell_line=cell_line,
                direction=direction,
                representation=representation,
                row_id=int(row_id),
                gallery_definition_id=gallery_definition_id,
                pos_definition_id=pos_definition_id,
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
                treated_cell_id=treated_cell_id,
                query_target_token=query_target_token,
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
            dataset=dataset,
            cell_line=cell_line,
            direction=direction,
            representation=representation,
            reason="m_pos0_after_intersection",
            n_query_rows_before=int(query_after.size),
            n_query_rows_after=int(query_after.size - acc.n_mpos0),
            n_gallery_member_rows_before=int(gallery_payload.gallery_member_rows_after.size),
            n_gallery_member_rows_after=int(gallery_payload.gallery_member_rows_after.size),
            n_gallery_items_before=int(len(gallery_payload.gallery_ids)),
            n_gallery_items_after=int(len(gallery_payload.gallery_ids)),
            notes="Non-empty gallery remained, but query positives were absent after valid intersection",
        ),
    )

    summary_row = build_summary_row(acc)
    chance_row = build_chance_row(summary_row)
    readiness = {
        "dataset": dataset,
        "cell_line": cell_line,
        "direction": direction,
        "representation": representation,
        "n_query_rows_before": int(query_before.size),
        "n_query_rows_after_valid_mask": int(query_after.size),
        "n_gallery_member_rows_before": int(gallery_payload.gallery_member_rows_before.size),
        "n_gallery_member_rows_after_valid_mask": int(gallery_payload.gallery_member_rows_after.size),
        "n_gallery_items_after_filter": int(len(gallery_payload.gallery_ids)),
        "n_mpos0": int(acc.n_mpos0),
        "n_valid": int(acc.n_valid),
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

    runs_dir = resolve_config_path(project_root, str(config["paths"]["runs_dir"]))
    task2_snapshot = (project_root / EXPECTED_TASK2_SNAPSHOT).resolve()
    expected_snapshot = (project_root / EXPECTED_TASK2_SNAPSHOT).resolve()
    if task2_snapshot != expected_snapshot:
        print("[ERROR] Snapshot isolation check failed unexpectedly.", file=sys.stderr)
        return 4
    if not task2_snapshot.is_dir():
        print(f"[ERROR] Missing corrected Task2 snapshot root: {task2_snapshot}", file=sys.stderr)
        return 5

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
                "seed": GLOBAL_SEED,
                "hit_k_list": list(HIT_K_LIST),
                "chance_identity_tol": CHANCE_IDENTITY_TOL,
                "query_batch_rows": QUERY_BATCH_ROWS,
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "task2_snapshot_v2_isolation",
            "pass": True,
            "details": {
                "rules": ["S5 reads exclusively from data/task2_snapshot_v2/"],
                "task2_snapshot": str(task2_snapshot),
            },
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "dataset_direction_representation_order_frozen",
            "pass": True,
            "details": {
                "dataset_order": list(DATASET_ORDER),
                "direction_order": list(DIRECTION_ORDER),
                "representation_order": list(REPRESENTATION_ORDER),
            },
            "counterexamples": [],
        }
    )

    try:
        snapshot_manifest, shared_manifest_inputs = load_snapshot_manifest(task2_snapshot)
        coverage, eligible_coverage, _supported_grid, shared_table_inputs = load_shared_contract_tables(task2_snapshot)
        registry, available_rep_scope, registry_inputs = load_available_rep_scope(task2_snapshot, eligible_coverage)
        source_target_order_by_slice, source_order_inputs = build_source_target_order_by_slice(task2_snapshot, eligible_coverage)
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
        return 6

    input_paths.extend(shared_manifest_inputs)
    input_paths.extend(shared_table_inputs)
    input_paths.extend(registry_inputs)
    input_paths.extend(source_order_inputs)

    cell_line_order_map = build_cell_line_order_map(snapshot_manifest)
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
            "name": "availability_registry_scope_policy",
            "pass": True,
            "details": {
                "available_rep_scope_rows": int(len(available_rep_scope)),
                "available_by_dataset": {
                    str(dataset): int(count)
                    for dataset, count in available_rep_scope.groupby("dataset", sort=False).size().to_dict().items()
                },
            },
            "counterexamples": [],
        }
    )

    dataset_bundles = build_dataset_bundles(task2_snapshot)
    physical_specs = build_physical_representation_specs(task2_snapshot)

    contexts: Dict[str, Dict[str, Any]] = {}
    representation_details: List[Dict[str, Any]] = []
    membership_details: List[Dict[str, object]] = []

    try:
        for dataset in DATASET_ORDER:
            bundle = dataset_bundles[dataset]
            subtree_manifest_path = task2_snapshot / bundle.subtree_manifest_relpath
            if not subtree_manifest_path.is_file():
                raise FileNotFoundError(f"Missing required dataset input for {dataset}: {subtree_manifest_path}")
            input_paths.append(subtree_manifest_path)
            subtree_manifest = json.loads(subtree_manifest_path.read_text(encoding="utf-8"))
            if subtree_manifest.get("analysis_key_fields") != ["dataset", "cell_line", "target_token"]:
                raise ValueError(f"{dataset} subtree_manifest analysis_key_fields drift")
            if subtree_manifest.get("representation_scope") != list(snapshot_manifest["datasets"][dataset]["representations"]):
                raise ValueError(f"{dataset} subtree_manifest representation_scope drift")

            eligible_subset = eligible_coverage.loc[eligible_coverage["dataset"].eq(dataset)].copy().reset_index(drop=True)
            eligible_cell_lines = eligible_subset["cell_line"].drop_duplicates().tolist()

            delta_meta, delta_inputs = load_delta_meta(
                dataset=dataset,
                delta_meta_path=task2_snapshot / bundle.delta_meta_relpath,
                eligible_cell_lines=eligible_cell_lines,
            )
            input_paths.extend(delta_inputs)
            fail_on_forbidden_cell_line(delta_meta, ["cell_line"], f"{dataset} delta_meta.csv")

            membership = build_membership_from_delta_meta(dataset, delta_meta)
            parity_membership, membership_inputs = load_membership_parity_table(dataset, task2_snapshot / bundle.membership_relpath)
            input_paths.extend(membership_inputs)
            validate_membership_parity(
                dataset=dataset,
                exploded_membership=membership,
                parity_membership=parity_membership,
            )

            (
                chem_query_rows_by_cell_line,
                gen_query_rows_by_cell_line,
                chem_target_rows_by_key,
                gen_target_rows_by_key,
            ) = build_rows_by_slice_and_target(
                dataset=dataset,
                membership=membership,
                eligible_coverage_subset=eligible_subset,
                source_target_order_by_slice=source_target_order_by_slice,
            )

            contexts[dataset] = {
                "delta_meta": delta_meta,
                "membership": membership,
                "chem_query_rows_by_cell_line": chem_query_rows_by_cell_line,
                "gen_query_rows_by_cell_line": gen_query_rows_by_cell_line,
                "chem_target_rows_by_key": chem_target_rows_by_key,
                "gen_target_rows_by_key": gen_target_rows_by_key,
            }
            membership_details.append(
                {
                    "dataset": dataset,
                    "delta_meta_rows": int(len(delta_meta)),
                    "exploded_membership_rows": int(len(membership)),
                    "cell_lines": eligible_cell_lines,
                }
            )
    except Exception as exc:
        assertions.append(
            {
                "name": "dataset_inputs_ready",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed loading dataset inputs: {exc}", file=sys.stderr)
        return 7

    no_forbidden_1hae = (
        not coverage["cell_line"].eq("1HAE").any()
        and not registry["cell_line"].eq("1HAE").any()
        and all(not contexts[dataset]["delta_meta"]["cell_line"].eq("1HAE").any() for dataset in DATASET_ORDER)
    )
    assertions.append(
        {
            "name": "no_forbidden_1HAE_in_analysis_inputs",
            "pass": bool(no_forbidden_1hae),
            "details": {
                "rules": ["Normalized Task2 analysis-facing inputs must not contain cell_line=1HAE"],
            },
            "counterexamples": [] if no_forbidden_1hae else [{"cell_line": "1HAE"}],
        }
    )
    assertions.append(
        {
            "name": "membership_parity_from_delta_meta",
            "pass": True,
            "details": {"datasets": membership_details},
            "counterexamples": [],
        }
    )

    per_query_path = stage_dir / "task2_retrieval_per_query.parquet"
    sink = PerQueryParquetSink(per_query_path, per_query_schema())

    summary_rows: List[Dict[str, object]] = []
    chance_rows: List[Dict[str, object]] = []
    attrition_rows: List[Dict[str, object]] = []
    readiness_rows: List[Dict[str, object]] = []

    try:
        for dataset in DATASET_ORDER:
            delta_meta = contexts[dataset]["delta_meta"]
            target_raw = delta_meta["target_raw"].astype(str).tolist()
            target_tokens_text = delta_meta["target_tokens"].astype(str).tolist()
            query_n_targets = delta_meta["parsed_target_tokens"].map(len).to_numpy(dtype=np.int64, copy=False)
            query_time = pd.to_numeric(delta_meta["time"], errors="coerce").to_numpy(dtype=np.float64, copy=False)
            query_dose_value = pd.to_numeric(delta_meta["dose_value"], errors="coerce").to_numpy(dtype=np.float64, copy=False)
            treated_cell_id = delta_meta["treated_cell_id"].astype(str).tolist()
            parsed_target_tokens = delta_meta["parsed_target_tokens"].tolist()

            for representation in REPRESENTATION_ORDER:
                available_subset = available_rep_scope.loc[
                    available_rep_scope["dataset"].eq(dataset) & available_rep_scope["representation"].eq(representation)
                ].copy()
                if available_subset.empty:
                    continue

                spec = physical_specs[dataset][representation]
                payload, rep_inputs, rep_detail = validate_representation_inputs(
                    snapshot_root=task2_snapshot,
                    spec=spec,
                    delta_meta=delta_meta,
                )
                representation_details.append(rep_detail)
                input_paths.extend(rep_inputs)

                arr = np.load(payload.array_path, mmap_mode="r")
                cell_lines = available_subset["cell_line"].tolist()
                cell_lines = sorted(cell_lines, key=lambda cell_line: cell_line_order_map[dataset][cell_line])
                for cell_line in cell_lines:
                    target_order = source_target_order_by_slice[(dataset, cell_line)]
                    c2g_summary, c2g_chance, c2g_attrition, c2g_readiness = run_direction_for_slice(
                        arr=arr,
                        dataset=dataset,
                        cell_line=cell_line,
                        direction="C2G",
                        representation=representation,
                        valid_mask=payload.valid_mask,
                        target_order=target_order,
                        query_row_ids=contexts[dataset]["chem_query_rows_by_cell_line"][cell_line],
                        gallery_target_rows_by_key=contexts[dataset]["gen_target_rows_by_key"],
                        target_raw=target_raw,
                        target_tokens_text=target_tokens_text,
                        query_n_targets=query_n_targets,
                        query_time=query_time,
                        query_dose_value=query_dose_value,
                        treated_cell_id=treated_cell_id,
                        parsed_target_tokens=parsed_target_tokens,
                        sink=sink,
                    )
                    summary_rows.append(c2g_summary)
                    chance_rows.append(c2g_chance)
                    attrition_rows.extend(c2g_attrition)
                    readiness_rows.append(c2g_readiness)

                    g2c_summary, g2c_chance, g2c_attrition, g2c_readiness = run_direction_for_slice(
                        arr=arr,
                        dataset=dataset,
                        cell_line=cell_line,
                        direction="G2C",
                        representation=representation,
                        valid_mask=payload.valid_mask,
                        target_order=target_order,
                        query_row_ids=contexts[dataset]["gen_query_rows_by_cell_line"][cell_line],
                        gallery_target_rows_by_key=contexts[dataset]["chem_target_rows_by_key"],
                        target_raw=target_raw,
                        target_tokens_text=target_tokens_text,
                        query_n_targets=query_n_targets,
                        query_time=query_time,
                        query_dose_value=query_dose_value,
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
        return 8

    sink.close()

    summary_df = sort_summary_frame(pd.DataFrame(summary_rows), cell_line_order_map)
    chance_df = sort_chance_frame(pd.DataFrame(chance_rows), cell_line_order_map)
    attrition_df = sort_attrition_frame(pd.DataFrame(attrition_rows), cell_line_order_map)
    summary_long_df = sort_summary_long_frame(pd.DataFrame(build_summary_long_rows(summary_rows)), cell_line_order_map)

    summary_path = stage_dir / "task2_retrieval_summary.csv"
    chance_path = stage_dir / "task2_chance_identity_check.csv"
    attrition_path = stage_dir / "task2_retrieval_attrition.csv"
    summary_long_path = stage_dir / "task2_retrieval_summary_long.csv"

    write_csv(summary_df, summary_path)
    write_csv(chance_df, chance_path)
    write_csv(attrition_df, attrition_path)
    write_csv(summary_long_df, summary_long_path)

    expected_summary_grid = build_expected_summary_grid(available_rep_scope)
    actual_summary_keys = summary_df[["dataset", "cell_line", "direction", "representation"]].copy()
    assert_unique_key(actual_summary_keys, ["dataset", "cell_line", "direction", "representation"], "actual summary grid")
    summary_grid_pass = len(actual_summary_keys) == len(expected_summary_grid) and actual_summary_keys.equals(
        sort_summary_frame(
            expected_summary_grid.assign(
                gallery_definition_id="",
                pos_definition_id="",
                n_total=0,
                n_valid=0,
                n_excluded_missing_metric_or_mpos0=0,
                N_gallery_mean=0.0,
                N_gallery_max=0,
                m_pos_mean=np.nan,
                m_pos_p50=np.nan,
                m_pos_p90=np.nan,
                mean_mrr_raw=np.nan,
                mean_expected_mrr_chance=np.nan,
                mean_mrr_corrected=np.nan,
                mean_hit1_raw=np.nan,
                mean_expected_hit1_chance=np.nan,
                mean_hit1_corrected=np.nan,
                mean_hit5_raw=np.nan,
                mean_expected_hit5_chance=np.nan,
                mean_hit5_corrected=np.nan,
                mean_hit10_raw=np.nan,
                mean_expected_hit10_chance=np.nan,
                mean_hit10_corrected=np.nan,
            )[["dataset", "cell_line", "direction", "representation", "gallery_definition_id", "pos_definition_id", "n_total", "n_valid", "n_excluded_missing_metric_or_mpos0", "N_gallery_mean", "N_gallery_max", "m_pos_mean", "m_pos_p50", "m_pos_p90", "mean_mrr_raw", "mean_expected_mrr_chance", "mean_mrr_corrected", "mean_hit1_raw", "mean_expected_hit1_chance", "mean_hit1_corrected", "mean_hit5_raw", "mean_expected_hit5_chance", "mean_hit5_corrected", "mean_hit10_raw", "mean_expected_hit10_chance", "mean_hit10_corrected"]],
            cell_line_order_map,
        )[["dataset", "cell_line", "direction", "representation"]]
    )
    expected_by_dataset = {
        str(dataset): int(count * len(DIRECTION_ORDER))
        for dataset, count in available_rep_scope.groupby("dataset", sort=False).size().to_dict().items()
    }
    actual_by_dataset = {
        str(dataset): int(count)
        for dataset, count in summary_df.groupby("dataset", sort=False).size().to_dict().items()
    }
    assertions.append(
        {
            "name": "summary_grid_expected_vs_actual",
            "pass": bool(summary_grid_pass),
            "details": {
                "expected_summary_rows": int(len(expected_summary_grid)),
                "actual_summary_rows": int(len(summary_df)),
                "expected_by_dataset": expected_by_dataset,
                "actual_by_dataset": actual_by_dataset,
            },
            "counterexamples": []
            if summary_grid_pass
            else actual_summary_keys.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    expected_summary_long_rows = int(len(summary_df) * len(SUMMARY_LONG_METRICS))
    summary_long_pass = int(len(summary_long_df)) == expected_summary_long_rows
    assert_unique_key(
        summary_long_df,
        ["dataset", "cell_line", "direction", "representation", "metric_name"],
        "task2_retrieval_summary_long.csv",
    )
    assertions.append(
        {
            "name": "summary_long_expected_vs_actual",
            "pass": bool(summary_long_pass),
            "details": {
                "expected_summary_long_rows": expected_summary_long_rows,
                "actual_summary_long_rows": int(len(summary_long_df)),
            },
            "counterexamples": []
            if summary_long_pass
            else summary_long_df.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    not_applicable = registry.loc[
        registry["availability_status"].ne(STATUS_AVAILABLE),
        ["dataset", "cell_line", "representation"],
    ].drop_duplicates()
    if summary_df.empty:
        bad_scope_rows = pd.DataFrame(columns=["dataset", "cell_line", "representation"])
    else:
        bad_scope_rows = summary_df[["dataset", "cell_line", "representation"]].merge(
            not_applicable,
            on=["dataset", "cell_line", "representation"],
            how="inner",
            sort=False,
        )
    assertions.append(
        {
            "name": "no_not_applicable_scope_rows_materialized",
            "pass": bad_scope_rows.empty,
            "details": {
                "summary_rows": int(len(summary_df)),
                "attrition_rows": int(len(attrition_df)),
            },
            "counterexamples": bad_scope_rows.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    empty_gallery_pairs = set(
        tuple(row)
        for row in attrition_df.loc[
            attrition_df["reason"].eq("gallery_target_empty_after_filter"),
            ["dataset", "cell_line", "direction", "representation"],
        ].itertuples(index=False, name=None)
    )
    mpos0_pairs = set(
        tuple(row)
        for row in attrition_df.loc[
            attrition_df["reason"].eq("m_pos0_after_intersection"),
            ["dataset", "cell_line", "direction", "representation"],
        ].itertuples(index=False, name=None)
    )
    empty_gallery_distinct_pass = len(empty_gallery_pairs & mpos0_pairs) == 0
    assertions.append(
        {
            "name": "empty_gallery_and_mpos0_semantics_distinct",
            "pass": bool(empty_gallery_distinct_pass),
            "details": {
                "empty_gallery_rows": int(len(empty_gallery_pairs)),
                "mpos0_rows": int(len(mpos0_pairs)),
            },
            "counterexamples": [
                {
                    "dataset": row[0],
                    "cell_line": row[1],
                    "direction": row[2],
                    "representation": row[3],
                }
                for row in list(empty_gallery_pairs & mpos0_pairs)[:MAX_COUNTEREXAMPLES]
            ],
        }
    )

    denominator_failures = summary_df.loc[
        summary_df["n_total"] != (
            summary_df["n_valid"] + summary_df["n_excluded_missing_metric_or_mpos0"]
        )
    ]
    assertions.append(
        {
            "name": "denominator_conservation",
            "pass": denominator_failures.empty,
            "details": {"n_checked": int(len(summary_df))},
            "counterexamples": denominator_failures.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    tol_cols = ["abs_delta_mrr", "abs_delta_hit1", "abs_delta_hit5", "abs_delta_hit10"]
    chance_finite = chance_df.dropna(subset=tol_cols, how="any")
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
                "tolerance": CHANCE_IDENTITY_TOL,
                "n_checked": int(len(chance_finite)),
            },
            "counterexamples": chance_violations.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )

    core_outputs = [per_query_path, summary_path, summary_long_path, attrition_path, chance_path]
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
        str(path.resolve())
        for path in sorted(set(input_paths), key=lambda p: str(p))
        if not path.resolve().is_relative_to(task2_snapshot.resolve())
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "task2_snapshot": str(task2_snapshot),
                "n_inputs": int(len(set(input_paths))),
            },
            "counterexamples": [{"bad_input": path} for path in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    per_query_rows = int(summary_df["n_valid"].sum()) if not summary_df.empty else 0
    completed_at = utc_now_iso()
    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    manifest_path = stage_dir / "manifest.json"

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s5_task2_retrieval_multisource.py").resolve()),
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
            "dataset_order": list(DATASET_ORDER),
            "direction_order": list(DIRECTION_ORDER),
            "representation_order": list(REPRESENTATION_ORDER),
            "target_membership_source": "delta_meta.target_tokens",
            "gallery_definition_ids": {
                "C2G": C2G_GALLERY_DEFINITION,
                "G2C": G2C_GALLERY_DEFINITION,
            },
            "pos_definition_ids": {
                "C2G": C2G_POS_DEFINITION,
                "G2C": G2C_POS_DEFINITION,
            },
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths), key=lambda p: str(p))],
        "outputs": [
            str(path.resolve())
            for path in [
                per_query_path,
                summary_path,
                summary_long_path,
                attrition_path,
                chance_path,
                run_manifest_path,
                audit_assertions_path,
                manifest_path,
            ]
        ],
        "summary": {
            "n_summary_rows": int(len(summary_df)),
            "n_summary_long_rows": int(len(summary_long_df)),
            "n_per_query_rows": per_query_rows,
            "n_attrition_rows": int(len(attrition_df)),
            "n_chance_rows": int(len(chance_df)),
            "expected_summary_rows": int(len(expected_summary_grid)),
            "expected_summary_long_rows": expected_summary_long_rows,
            "expected_by_dataset": expected_by_dataset,
            "actual_by_dataset": actual_by_dataset,
            "representation_details": representation_details,
            "readiness_rows": readiness_rows,
        },
    }

    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})
    write_json(manifest_path, {"stage": STAGE, "files": build_stage_manifest(stage_dir)})

    if not summary_grid_pass:
        print("[ERROR] Summary grid mismatch against input-derived expected grid.", file=sys.stderr)
        return 9
    if not summary_long_pass:
        print("[ERROR] Summary-long row count mismatch.", file=sys.stderr)
        return 10
    if not empty_gallery_distinct_pass:
        print("[ERROR] Empty-gallery and m_pos0 attrition semantics were mixed.", file=sys.stderr)
        return 11
    if not denominator_failures.empty:
        print("[ERROR] Denominator conservation failed.", file=sys.stderr)
        return 12
    if not chance_pass:
        print(f"[ERROR] chance identity check failed: tolerance={CHANCE_IDENTITY_TOL}", file=sys.stderr)
        return 13

    return 0


if __name__ == "__main__":
    sys.exit(main())
