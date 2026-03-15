# SCRIPT_HEADER_CONTRACT
# Script: scripts/s3_build_task2_multisource_snapshot.py
# Purpose: Build corrected Task2 multisource snapshot artifacts under data/task2_snapshot_v2
#   by reusing the audited legacy scPerturb-K562 Task2 snapshot and materializing
#   downstream-friendly LINCS Gene/Pathway derived assets from Task1 LINCS sources.
# Inputs:
#   - Legacy scPerturb Task2 snapshot: config/config.yaml::paths.task2_snapshot (must resolve to data/task2_snapshot_v1)
#   - Task1 LINCS assets:
#       config/config.yaml::paths.task1_snapshot/lincs/lincs-engine1-meta.csv
#       config/config.yaml::paths.task1_snapshot/lincs/lincs-engine1-gene-delta.pt
#       config/config.yaml::paths.task1_snapshot/lincs/lincs-gene-alignment.csv
#       config/config.yaml::paths.task1_snapshot/pathway/hallmark-w-2477x50.npy
#       config/config.yaml::paths.task1_snapshot/pathway/lincs-pathway-policy.json
# Outputs:
#   - Shared snapshot artifacts under --output-root (default data/task2_snapshot_v2):
#       snapshot_manifest.json
#       source_lineage_manifest.json
#       task2_pairs_coverage.csv
#       representation_availability_registry.csv
#       pathway/{hallmark-w-2477x50.npy,lincs-pathway-policy.json}
#   - scPerturb subtree:
#       scperturb_k562/{raw files, derived files, fm files, subtree_manifest.json}
#       scperturb_k562/derived/task2_row_membership.parquet
#   - LINCS subtree:
#       lincs/{lincs-engine1-meta.csv,lincs-gene-alignment.csv,task2_lincs_pairs.csv,subtree_manifest.json}
#       lincs/derived/{gene_delta.npy,delta_meta.csv,pathway_delta.npy,task2_row_membership.parquet}
#   - Optional/helper run outputs:
#       runs/<run_id>/s3_build_task2_multisource_snapshot/task2_post_build_inventory.csv
#   - AVCP Artifacts:
#       runs/<run_id>/s3_build_task2_multisource_snapshot/{run_manifest.json,audit_assertions.json,manifest.json}
# Side Effects:
#   - Creates/updates snapshot directory under --output-root
#   - Creates isolated run directory: runs/<run_id>/s3_build_task2_multisource_snapshot/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task1_snapshot
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s3_build_task2_multisource_snapshot.py --run-id <run_id> --output-root data/task2_snapshot_v2 --seed 619
# Failure Modes:
#   - Seed != 619 -> exit non-zero
#   - config.paths.task2_snapshot does not resolve to legacy data/task2_snapshot_v1 -> exit non-zero
#   - Missing required legacy scPerturb or Task1 LINCS assets -> exit non-zero
#   - LINCS source tensor/meta row mismatch or missing required source keys -> exit non-zero
#   - No eligible LINCS Task2 cohorts found -> exit non-zero
#   - Output routing outside --output-root or runs/<run_id>/... -> exit non-zero
# Last Updated: 2026-03-10

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
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import torch
import yaml

STAGE = "s3_build_task2_multisource_snapshot"
CONFIG_PATH = Path("config/config.yaml")

GLOBAL_SEED = 619
MAX_COUNTEREXAMPLES = 5
LINCS_BATCH_ROWS = 4096
EXPECTED_LEGACY_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
DEFAULT_OUTPUT_ROOT = Path("data/task2_snapshot_v2")
EXPECTED_TASK1_SNAPSHOT = Path("data/task1_snapshot_v1")

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
LINCS_REPRESENTATIONS = {"Gene", "Pathway"}
SCPERTURB_K562_REPRESENTATIONS = set(REPRESENTATION_ORDER)
LINCS_ALLOWED_PERT_TYPES = {"drug", "sh", "crispr", "oe"}
LINCS_CELL_LINE_ALIAS_MAP = {
    "1HAE": "HA1E",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S3 corrected Task2 multisource snapshot builder")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=True)
        handle.write("\n")


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(frame, preserve_index=False)
    pq.write_table(table, path)


def ensure_required_columns(frame: pd.DataFrame, required: Sequence[str], name: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns for {name}: {missing}")


def normalize_scalar(value: object, default: str = "NA") -> str:
    if pd.isna(value):
        return default
    text = str(value).strip()
    return text if text else default


def normalize_required_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_lincs_cell_line(value: object) -> str:
    raw = normalize_required_text(value)
    return LINCS_CELL_LINE_ALIAS_MAP.get(raw, raw)


def to_float_or_nan(value: object) -> float:
    if pd.isna(value):
        return float("nan")
    return float(value)


def init_global_rng(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def join_tokens(tokens: Sequence[str]) -> str:
    return ";".join(tokens)


def tokenize_scperturb_target(raw_target: object) -> Tuple[str, ...]:
    raw = normalize_required_text(raw_target)
    tokens: List[str] = []
    for chunk in raw.split(";"):
        token = chunk.strip()
        if token:
            tokens.append(token)
    return tuple(sorted(set(tokens)))


def tokenize_lincs_target(raw_target: object) -> Tuple[str, ...]:
    raw = normalize_required_text(raw_target)
    if not raw:
        return tuple()
    normalized = raw.replace(";", "|")
    tokens = [chunk.strip() for chunk in normalized.split("|") if chunk.strip()]
    return tuple(sorted(set(tokens)))


def link_or_copy(src: Path, dst: Path) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        if dst.is_file() and src.is_file() and dst.stat().st_size == src.stat().st_size:
            return "reuse"
        dst.unlink()
    try:
        os.link(src, dst)
        return "hardlink"
    except OSError:
        shutil.copy2(src, dst)
        return "copy"


def symlink_or_reuse(src: Path, dst: Path) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    src_resolved = src.resolve()
    if dst.exists() or dst.is_symlink():
        if dst.is_symlink() and dst.resolve() == src_resolved:
            return "reuse"
        dst.unlink()
    os.symlink(src_resolved, dst)
    return "symlink"


def local_output_path(path: Path) -> str:
    return str(path.absolute())


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


def legacy_scperturb_required_relpaths() -> List[Path]:
    relpaths = [
        Path("k562/CRISPR_counts.pt"),
        Path("k562/CRISPR_meta.csv"),
        Path("k562/Drug_counts.pt"),
        Path("k562/Drug_meta.csv"),
        Path("k562/Common_Targets_K562.csv"),
        Path("k562/shared_var_names.csv"),
        Path("k562/derived/pair_list.parquet"),
        Path("k562/derived/delta_meta.csv"),
        Path("k562/derived/gene_delta.npy"),
        Path("k562/derived/pathway_delta.npy"),
    ]
    for model in REPRESENTATION_ORDER[2:]:
        relpaths.extend(
            [
                Path(f"k562/fm/{model}/fm_delta.npy"),
                Path(f"k562/fm/{model}/fm_delta_meta.csv"),
                Path(f"k562/fm/{model}/delta_operator_policy.json"),
            ]
        )
    return relpaths


def validate_legacy_scperturb_snapshot(legacy_root: Path) -> Dict[str, object]:
    missing = [str(rel) for rel in legacy_scperturb_required_relpaths() if not (legacy_root / rel).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing legacy scPerturb Task2 files: {missing[:MAX_COUNTEREXAMPLES]}")

    delta_meta_path = legacy_root / "k562/derived/delta_meta.csv"
    gene_delta_path = legacy_root / "k562/derived/gene_delta.npy"
    pathway_delta_path = legacy_root / "k562/derived/pathway_delta.npy"

    delta_meta = pd.read_csv(delta_meta_path)
    ensure_required_columns(
        delta_meta,
        [
            "row_id",
            "treated_cell_id",
            "perturbation_class",
            "cell_line",
            "target_raw",
            "target_tokens",
        ],
        "legacy scPerturb delta_meta.csv",
    )
    if delta_meta["row_id"].duplicated().any():
        raise ValueError("Legacy scPerturb delta_meta row_id must be unique.")

    gene_shape = np.load(gene_delta_path, mmap_mode="r").shape
    pathway_shape = np.load(pathway_delta_path, mmap_mode="r").shape
    if int(gene_shape[0]) != len(delta_meta):
        raise ValueError(
            f"Legacy scPerturb gene delta row mismatch: delta_meta={len(delta_meta)} gene_delta={gene_shape[0]}"
        )
    if int(pathway_shape[0]) != len(delta_meta):
        raise ValueError(
            f"Legacy scPerturb pathway delta row mismatch: delta_meta={len(delta_meta)} pathway_delta={pathway_shape[0]}"
        )

    fm_summaries: Dict[str, Dict[str, int]] = {}
    for model in REPRESENTATION_ORDER[2:]:
        fm_meta_path = legacy_root / f"k562/fm/{model}/fm_delta_meta.csv"
        fm_delta_path = legacy_root / f"k562/fm/{model}/fm_delta.npy"
        fm_meta = pd.read_csv(fm_meta_path)
        ensure_required_columns(
            fm_meta,
            ["row_id", "treated_cell_id", "valid_mask", "invalid_reason", "model_name"],
            f"legacy scPerturb fm_delta_meta.csv::{model}",
        )
        if len(fm_meta) != len(delta_meta):
            raise ValueError(
                f"Legacy scPerturb FM meta row mismatch for {model}: fm_meta={len(fm_meta)} delta_meta={len(delta_meta)}"
            )
        fm_shape = np.load(fm_delta_path, mmap_mode="r").shape
        if int(fm_shape[0]) != len(delta_meta):
            raise ValueError(
                f"Legacy scPerturb FM delta row mismatch for {model}: fm_delta={fm_shape[0]} delta_meta={len(delta_meta)}"
            )
        fm_summaries[model] = {
            "n_total_rows": int(len(fm_meta)),
            "n_valid_rows": int(pd.Series(fm_meta["valid_mask"]).astype(bool).sum()),
            "n_invalid_rows": int(len(fm_meta) - pd.Series(fm_meta["valid_mask"]).astype(bool).sum()),
        }

    common_targets_path = legacy_root / "k562/Common_Targets_K562.csv"
    common_targets_df = pd.read_csv(common_targets_path)
    if "Common_Target" not in common_targets_df.columns:
        raise ValueError("Legacy Common_Targets_K562.csv must contain Common_Target")
    common_targets = (
        common_targets_df["Common_Target"]
        .fillna("")
        .astype(str)
        .str.strip()
        .replace("", np.nan)
        .dropna()
        .tolist()
    )
    if not common_targets:
        raise ValueError("Legacy Common_Targets_K562.csv cannot be empty.")

    return {
        "delta_meta_rows": int(len(delta_meta)),
        "gene_dim": int(gene_shape[1]),
        "pathway_dim": int(pathway_shape[1]),
        "common_targets": common_targets,
        "fm_summaries": fm_summaries,
    }


def materialize_legacy_scperturb_subtree(
    legacy_root: Path,
    output_root: Path,
) -> Tuple[List[Path], List[Dict[str, object]]]:
    outputs: List[Path] = []
    rows: List[Dict[str, object]] = []
    for rel in legacy_scperturb_required_relpaths():
        src = legacy_root / rel
        mapped_rel = Path("scperturb_k562") / rel.relative_to("k562")
        dst = output_root / mapped_rel
        mode = link_or_copy(src, dst)
        outputs.append(dst.resolve())
        rows.append(
            {
                "dataset": "scPerturb",
                "source": str(src.resolve()),
                "destination": str(dst.resolve()),
                "mode": mode,
            }
        )
    return outputs, rows


def build_scperturb_membership(
    legacy_delta_meta_path: Path,
    common_targets: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    delta_meta = pd.read_csv(legacy_delta_meta_path)
    ensure_required_columns(
        delta_meta,
        [
            "row_id",
            "treated_cell_id",
            "perturbation_class",
            "cell_line",
            "target_raw",
            "target_tokens",
        ],
        "legacy scPerturb delta_meta.csv",
    )
    membership = (
        delta_meta.assign(
            target_token=delta_meta["target_tokens"].map(tokenize_scperturb_target),
            dataset="scPerturb",
        )
        .explode("target_token")
        .rename(columns={"target_token": "target_token"})
    )
    membership["target_token"] = membership["target_token"].fillna("").astype(str)
    membership = membership.loc[membership["target_token"].ne("")].copy()
    membership = membership[
        ["row_id", "treated_cell_id", "dataset", "cell_line", "perturbation_class", "target_raw", "target_token"]
    ].sort_values(["row_id", "target_token"], kind="mergesort")
    membership.reset_index(drop=True, inplace=True)

    chem_counts = (
        membership.loc[membership["perturbation_class"].eq("Chemical")]
        .groupby("target_token", sort=False)
        .size()
        .to_dict()
    )
    gen_counts = (
        membership.loc[membership["perturbation_class"].eq("Genetic")]
        .groupby("target_token", sort=False)
        .size()
        .to_dict()
    )
    rows: List[Dict[str, object]] = []
    for token in common_targets:
        n_chem = int(chem_counts.get(token, 0))
        n_gen = int(gen_counts.get(token, 0))
        rows.append(
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target_token": token,
                "n_chem_instances": n_chem,
                "n_gen_instances": n_gen,
                "is_eligible_bool": bool(n_chem > 0 and n_gen > 0),
                "eligibility_reason": "chemical_and_genetic_present"
                if (n_chem > 0 and n_gen > 0)
                else "missing_chemical_or_genetic",
                "source_tag": "legacy_task2_v1_k562",
            }
        )
    coverage = pd.DataFrame(rows)
    return membership, coverage


def load_lincs_tensor_bundle(bundle_path: Path, expected_rows: int) -> Tuple[torch.Tensor, torch.Tensor]:
    bundle = torch.load(bundle_path, map_location="cpu", mmap=True)
    if not isinstance(bundle, dict):
        raise ValueError("LINCS tensor bundle must be a dict.")
    for key in ["y_delta_gene", "y_delta_pathway"]:
        if key not in bundle:
            raise ValueError(f"LINCS tensor bundle missing key {key!r}.")
        if not torch.is_tensor(bundle[key]):
            raise ValueError(f"LINCS tensor bundle key {key!r} must be a tensor.")
    gene = bundle["y_delta_gene"]
    pathway = bundle["y_delta_pathway"]
    if int(gene.shape[0]) != expected_rows:
        raise ValueError(
            f"LINCS gene tensor row mismatch: meta={expected_rows}, gene_tensor={int(gene.shape[0])}"
        )
    if int(pathway.shape[0]) != expected_rows:
        raise ValueError(
            f"LINCS pathway tensor row mismatch: meta={expected_rows}, pathway_tensor={int(pathway.shape[0])}"
        )
    if int(gene.shape[1]) != 2477:
        raise ValueError(f"LINCS gene tensor width must be 2477, got {int(gene.shape[1])}")
    if int(pathway.shape[1]) != 50:
        raise ValueError(f"LINCS pathway tensor width must be 50, got {int(pathway.shape[1])}")
    return gene, pathway


def normalize_lincs_metadata(lincs_meta: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    required_cols = [
        "unique_id",
        "cell_line",
        "pert_type",
        "target",
        "dose_val",
        "time_val",
        "paired_control_id",
        "sig_id",
    ]
    ensure_required_columns(lincs_meta, required_cols, "lincs-engine1-meta.csv")

    meta = lincs_meta.copy()
    meta["source_row_index"] = np.arange(len(meta), dtype=np.int64)
    meta["treated_cell_id"] = meta["unique_id"].map(normalize_required_text)
    fallback_mask = meta["treated_cell_id"].eq("")
    if fallback_mask.any():
        meta.loc[fallback_mask, "treated_cell_id"] = meta.loc[fallback_mask, "source_row_index"].map(
            lambda idx: f"LINCS_ROW_{int(idx)}"
        )

    meta["cell_line_raw"] = meta["cell_line"].map(normalize_required_text)
    meta["cell_line_norm"] = meta["cell_line_raw"].map(normalize_lincs_cell_line)
    meta["target_raw"] = meta["target"].map(normalize_required_text)
    meta["pert_type_raw"] = meta["pert_type"].map(normalize_required_text)
    meta["pert_type_norm"] = meta["pert_type_raw"].str.lower()
    meta["time"] = pd.to_numeric(meta["time_val"], errors="coerce")
    meta["dose_value"] = pd.to_numeric(meta["dose_val"], errors="coerce")
    meta["paired_control_id_norm"] = meta["paired_control_id"].map(normalize_required_text)
    meta["sig_id_norm"] = meta["sig_id"].map(normalize_required_text)

    exclusion_reason = pd.Series(pd.NA, index=meta.index, dtype="object")
    exclusion_reason = exclusion_reason.mask(meta["pert_type_raw"].eq(""), "missing_pert_type")
    unsupported_mask = exclusion_reason.isna() & (~meta["pert_type_norm"].isin(LINCS_ALLOWED_PERT_TYPES))
    exclusion_reason = exclusion_reason.mask(unsupported_mask, "unsupported_pert_type")
    exclusion_reason = exclusion_reason.mask(
        exclusion_reason.isna() & meta["cell_line_norm"].eq(""),
        "missing_cell_line",
    )
    token_lists = meta["target_raw"].map(tokenize_lincs_target)
    empty_target_mask = token_lists.map(len).eq(0)
    exclusion_reason = exclusion_reason.mask(
        exclusion_reason.isna() & empty_target_mask,
        "missing_target",
    )

    meta["target_tokens_all"] = token_lists
    meta["perturbation_class"] = meta["pert_type_norm"].map(
        {"drug": "Chemical", "sh": "Genetic", "crispr": "Genetic", "oe": "Genetic"}
    )
    meta["exclusion_reason"] = exclusion_reason

    excluded = meta.loc[meta["exclusion_reason"].notna(), [
        "source_row_index",
        "treated_cell_id",
        "cell_line_raw",
        "cell_line_norm",
        "pert_type_raw",
        "target_raw",
        "exclusion_reason",
    ]].rename(columns={"cell_line_norm": "cell_line"})

    supported = meta.loc[meta["exclusion_reason"].isna(), [
        "source_row_index",
        "treated_cell_id",
        "cell_line_raw",
        "cell_line_norm",
        "perturbation_class",
        "target_raw",
        "target_tokens_all",
        "time",
        "dose_value",
        "pert_type_raw",
        "paired_control_id_norm",
        "sig_id_norm",
    ]].rename(
        columns={
            "cell_line_norm": "cell_line",
            "paired_control_id_norm": "paired_control_id",
            "sig_id_norm": "sig_id",
        }
    )

    if supported["treated_cell_id"].duplicated().any():
        dupes = supported.loc[supported["treated_cell_id"].duplicated(), "treated_cell_id"].head(
            MAX_COUNTEREXAMPLES
        )
        raise ValueError(f"LINCS treated_cell_id must be unique, examples={dupes.tolist()}")

    return supported.reset_index(drop=True), excluded.reset_index(drop=True)


def build_lincs_task2_tables(
    supported_rows: pd.DataFrame,
    excluded_rows: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    full_membership = (
        supported_rows[["source_row_index", "treated_cell_id", "cell_line", "perturbation_class", "target_raw", "target_tokens_all"]]
        .explode("target_tokens_all")
        .rename(columns={"target_tokens_all": "target_token"})
    )
    full_membership["target_token"] = full_membership["target_token"].fillna("").astype(str)
    full_membership = full_membership.loc[full_membership["target_token"].ne("")].copy()

    counts = (
        full_membership.groupby(["cell_line", "target_token", "perturbation_class"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    if "Chemical" not in counts.columns:
        counts["Chemical"] = 0
    if "Genetic" not in counts.columns:
        counts["Genetic"] = 0

    counts["dataset"] = "LINCS"
    counts["n_chem_instances"] = counts["Chemical"].astype(np.int64)
    counts["n_gen_instances"] = counts["Genetic"].astype(np.int64)
    counts["is_eligible_bool"] = counts["n_chem_instances"].gt(0) & counts["n_gen_instances"].gt(0)
    counts["eligibility_reason"] = np.where(
        counts["is_eligible_bool"], "chemical_and_genetic_present", "missing_chemical_or_genetic"
    )
    counts["source_tag"] = "task1_lincs_supported_rows"
    coverage = counts[
        [
            "dataset",
            "cell_line",
            "target_token",
            "n_chem_instances",
            "n_gen_instances",
            "is_eligible_bool",
            "eligibility_reason",
            "source_tag",
        ]
    ].sort_values(["cell_line", "target_token"], kind="mergesort")
    coverage.reset_index(drop=True, inplace=True)

    eligible_pairs = coverage.loc[coverage["is_eligible_bool"]].copy()
    if eligible_pairs.empty:
        raise ValueError("No eligible LINCS Task2 cohorts found after supported-row normalization.")
    eligible_lookup: Dict[str, set[str]] = {
        cell_line: set(group["target_token"].tolist())
        for cell_line, group in eligible_pairs.groupby("cell_line", sort=False)
    }

    eligible_token_lists = [
        tuple(token for token in tokens if token in eligible_lookup.get(cell_line, set()))
        for cell_line, tokens in zip(supported_rows["cell_line"], supported_rows["target_tokens_all"])
    ]
    included_mask = pd.Series([len(tokens) > 0 for tokens in eligible_token_lists], index=supported_rows.index)

    cohort_excluded = supported_rows.loc[~included_mask, [
        "source_row_index",
        "treated_cell_id",
        "cell_line",
        "perturbation_class",
        "target_raw",
    ]].copy()
    cohort_excluded["exclusion_reason"] = "cohort_not_task2_eligible"

    included = supported_rows.loc[included_mask].copy()
    included["eligible_tokens"] = [eligible_token_lists[idx] for idx in included.index]
    included = included.sort_values("source_row_index", kind="mergesort").reset_index(drop=True)
    included["row_id"] = np.arange(len(included), dtype=np.int64)
    included["target_tokens"] = included["eligible_tokens"].map(join_tokens)
    included["specificity_tier"] = "NA"
    included["n_controls_used"] = 1
    included["dataset_side"] = "LINCS"
    included["seed"] = GLOBAL_SEED

    delta_meta = included[
        [
            "row_id",
            "source_row_index",
            "treated_cell_id",
            "cell_line_raw",
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
            "pert_type_raw",
            "paired_control_id",
            "sig_id",
        ]
    ].copy()

    membership = (
        included[
            [
                "row_id",
                "source_row_index",
                "treated_cell_id",
                "cell_line_raw",
                "cell_line",
                "perturbation_class",
                "target_raw",
                "eligible_tokens",
            ]
        ]
        .explode("eligible_tokens")
        .rename(columns={"eligible_tokens": "target_token"})
    )
    membership["dataset"] = "LINCS"
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
    ].sort_values(["row_id", "target_token"], kind="mergesort")
    membership.reset_index(drop=True, inplace=True)

    exclusion_summary = pd.concat([excluded_rows, cohort_excluded], ignore_index=True, sort=False)
    reason_counts = exclusion_summary["exclusion_reason"].fillna("NA").value_counts().sort_index().to_dict()
    included_by_class = included["perturbation_class"].value_counts().sort_index().to_dict()
    supported_by_class = supported_rows["perturbation_class"].value_counts().sort_index().to_dict()
    summary = {
        "n_source_rows_total": int(len(supported_rows) + len(excluded_rows)),
        "n_source_rows_supported": int(len(supported_rows)),
        "n_source_rows_excluded_pre_support": int(len(excluded_rows)),
        "n_rows_included": int(len(included)),
        "n_rows_excluded_post_support": int(len(cohort_excluded)),
        "excluded_by_reason": {str(k): int(v) for k, v in reason_counts.items()},
        "supported_rows_by_perturbation_class": {str(k): int(v) for k, v in supported_by_class.items()},
        "included_rows_by_perturbation_class": {str(k): int(v) for k, v in included_by_class.items()},
        "n_eligible_cohorts": int(len(eligible_pairs)),
        "eligible_cell_lines": sorted(eligible_lookup.keys()),
    }

    return delta_meta, membership, coverage, eligible_pairs.reset_index(drop=True), exclusion_summary, summary


def materialize_lincs_subset_tensor(
    source_tensor: torch.Tensor,
    source_row_indices: np.ndarray,
    output_path: Path,
) -> Tuple[int, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = int(source_row_indices.size)
    width = int(source_tensor.shape[1])
    out = np.lib.format.open_memmap(output_path, mode="w+", dtype=np.float32, shape=(n_rows, width))
    for start in range(0, n_rows, LINCS_BATCH_ROWS):
        end = min(start + LINCS_BATCH_ROWS, n_rows)
        idx_t = torch.as_tensor(source_row_indices[start:end], dtype=torch.long)
        block = source_tensor.index_select(0, idx_t).detach().cpu().numpy().astype(np.float32, copy=False)
        out[start:end] = block
    out.flush()
    del out
    return n_rows, width


def build_representation_availability_registry(coverage_df: pd.DataFrame) -> pd.DataFrame:
    eligible_pairs = coverage_df.loc[coverage_df["is_eligible_bool"]].copy()
    cohort_rows: List[Dict[str, object]] = []
    for dataset, dataset_frame in eligible_pairs.groupby("dataset", sort=False):
        cell_lines = sorted(dataset_frame["cell_line"].dropna().astype(str).unique().tolist())
        allowed = LINCS_REPRESENTATIONS if dataset == "LINCS" else SCPERTURB_K562_REPRESENTATIONS
        for cell_line in cell_lines:
            for representation in REPRESENTATION_ORDER:
                if representation in allowed:
                    status = "available"
                    reason = "representation_supported_for_dataset_scope"
                else:
                    status = "not_applicable_scope"
                    reason = "representation_not_supported_for_dataset_scope"
                cohort_rows.append(
                    {
                        "dataset": dataset,
                        "cell_line": cell_line,
                        "representation": representation,
                        "availability_status": status,
                        "availability_reason": reason,
                    }
                )
    registry = pd.DataFrame(cohort_rows)
    dataset_order = {"LINCS": 0, "scPerturb": 1}
    rep_order = {name: idx for idx, name in enumerate(REPRESENTATION_ORDER)}
    registry["_dataset_order"] = registry["dataset"].map(dataset_order)
    registry["_rep_order"] = registry["representation"].map(rep_order)
    registry = registry.sort_values(
        ["_dataset_order", "cell_line", "_rep_order"],
        kind="mergesort",
    ).drop(columns=["_dataset_order", "_rep_order"])
    registry.reset_index(drop=True, inplace=True)
    return registry


def build_post_build_inventory(
    scperturb_delta_meta_path: Path,
    scperturb_fm_root: Path,
    lincs_supported_rows: pd.DataFrame,
    lincs_included_delta_meta: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []

    sc_meta = pd.read_csv(scperturb_delta_meta_path)
    sc_counts_total = (
        sc_meta.groupby("perturbation_class", sort=False)
        .size()
        .rename(index={"Chemical": "Chemical", "Genetic": "Genetic"})
        .to_dict()
    )
    for representation in ["Gene", "Pathway"]:
        for perturbation_type in ["Chemical", "Genetic"]:
            total = int(sc_counts_total.get(perturbation_type, 0))
            rows.append(
                {
                    "dataset": "scPerturb",
                    "perturbation_type": perturbation_type,
                    "cell_line": "K562",
                    "representation": representation,
                    "stage": "post_build",
                    "n_total_instances": total,
                    "n_valid_instances": total,
                    "n_dropped": 0,
                    "drop_reason_breakdown": json.dumps({}, sort_keys=True, ensure_ascii=True),
                }
            )

    class_row_ids: Dict[str, pd.Series] = {
        label: sc_meta.loc[sc_meta["perturbation_class"].eq(label), "row_id"].astype(np.int64)
        for label in ["Chemical", "Genetic"]
    }
    for model in REPRESENTATION_ORDER[2:]:
        fm_meta = pd.read_csv(scperturb_fm_root / model / "fm_delta_meta.csv")
        valid_mask = pd.Series(fm_meta["valid_mask"]).astype(bool)
        for perturbation_type in ["Chemical", "Genetic"]:
            row_ids = class_row_ids.get(perturbation_type, pd.Series(dtype=np.int64))
            subset = fm_meta.loc[fm_meta["row_id"].isin(row_ids.tolist())].copy()
            total = int(len(subset))
            valid = int(pd.Series(subset["valid_mask"]).astype(bool).sum())
            invalid_counts = (
                subset.loc[~pd.Series(subset["valid_mask"]).astype(bool), "invalid_reason"]
                .fillna("missing_invalid_reason")
                .astype(str)
                .value_counts()
                .sort_index()
                .to_dict()
            )
            rows.append(
                {
                    "dataset": "scPerturb",
                    "perturbation_type": perturbation_type,
                    "cell_line": "K562",
                    "representation": model,
                    "stage": "post_build",
                    "n_total_instances": total,
                    "n_valid_instances": valid,
                    "n_dropped": int(total - valid),
                    "drop_reason_breakdown": json.dumps(
                        {str(k): int(v) for k, v in invalid_counts.items()},
                        sort_keys=True,
                        ensure_ascii=True,
                    ),
                }
            )

    lincs_supported = lincs_supported_rows.groupby(["cell_line", "perturbation_class"], sort=False).size()
    lincs_valid = lincs_included_delta_meta.groupby(["cell_line", "perturbation_class"], sort=False).size()
    for (cell_line, perturbation_type), total in lincs_supported.items():
        valid = int(lincs_valid.get((cell_line, perturbation_type), 0))
        drop_counts = {}
        if int(total) - valid > 0:
            drop_counts["cohort_not_task2_eligible"] = int(total) - valid
        for representation in ["Gene", "Pathway"]:
            rows.append(
                {
                    "dataset": "LINCS",
                    "perturbation_type": perturbation_type,
                    "cell_line": cell_line,
                    "representation": representation,
                    "stage": "post_build",
                    "n_total_instances": int(total),
                    "n_valid_instances": valid,
                    "n_dropped": int(total) - valid,
                    "drop_reason_breakdown": json.dumps(drop_counts, sort_keys=True, ensure_ascii=True),
                }
            )

    inventory = pd.DataFrame(rows)
    dataset_order = {"LINCS": 0, "scPerturb": 1}
    rep_order = {name: idx for idx, name in enumerate(REPRESENTATION_ORDER)}
    inventory["_dataset_order"] = inventory["dataset"].map(dataset_order)
    inventory["_rep_order"] = inventory["representation"].map(rep_order)
    inventory = inventory.sort_values(
        ["_dataset_order", "cell_line", "perturbation_type", "_rep_order"],
        kind="mergesort",
    ).drop(columns=["_dataset_order", "_rep_order"])
    inventory.reset_index(drop=True, inplace=True)
    return inventory


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

    init_global_rng(seed)

    task1_snapshot = resolve_config_path(project_root, str(config["paths"]["task1_snapshot"]))
    legacy_task2_root = resolve_config_path(project_root, str(config["paths"]["task2_snapshot"]))
    runs_dir = resolve_config_path(project_root, str(config["paths"]["runs_dir"]))
    output_root = resolve_config_path(project_root, str(args.output_root))

    expected_legacy_root = (project_root / EXPECTED_LEGACY_TASK2_SNAPSHOT).resolve()
    if legacy_task2_root != expected_legacy_root:
        print(
            "[ERROR] Legacy Task2 input mismatch: config.paths.task2_snapshot must resolve "
            f"to {expected_legacy_root}, got {legacy_task2_root}",
            file=sys.stderr,
        )
        return 4
    if output_root == legacy_task2_root:
        print(
            f"[ERROR] --output-root must not equal legacy Task2 root {legacy_task2_root}",
            file=sys.stderr,
        )
        return 5

    stage_dir = runs_dir / args.run_id / STAGE
    stage_dir.mkdir(parents=True, exist_ok=True)
    started_at = utc_now_iso()

    assertions: List[Dict[str, object]] = []
    input_paths: List[Path] = []
    output_paths: List[Path] = []
    lineage_rows: List[Dict[str, object]] = []

    assertions.append(
        {
            "name": "seed_locked_global_619",
            "pass": True,
            "details": {"seed": GLOBAL_SEED},
            "counterexamples": [],
        }
    )
    assertions.append(
        {
            "name": "output_root_cli_override_policy",
            "pass": True,
            "details": {
                "rules": [
                    "Corrected Task2 v2 writes to explicit --output-root",
                    "config.paths.task2_snapshot remains the legacy v1 input root during migration",
                ],
                "legacy_task2_root": str(legacy_task2_root),
                "output_root": str(output_root),
            },
            "counterexamples": [],
        }
    )

    try:
        legacy_summary = validate_legacy_scperturb_snapshot(legacy_task2_root)
    except Exception as exc:
        assertions.append(
            {
                "name": "legacy_scperturb_snapshot_ready",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Legacy scPerturb snapshot validation failed: {exc}", file=sys.stderr)
        return 6

    assertions.append(
        {
            "name": "legacy_scperturb_snapshot_ready",
            "pass": True,
            "details": {
                "delta_meta_rows": legacy_summary["delta_meta_rows"],
                "gene_dim": legacy_summary["gene_dim"],
                "pathway_dim": legacy_summary["pathway_dim"],
                "n_common_targets": len(legacy_summary["common_targets"]),
            },
            "counterexamples": [],
        }
    )

    required_lincs_paths = {
        "meta_csv": task1_snapshot / "lincs/lincs-engine1-meta.csv",
        "tensor_pt": task1_snapshot / "lincs/lincs-engine1-gene-delta.pt",
        "alignment_csv": task1_snapshot / "lincs/lincs-gene-alignment.csv",
        "pathway_w": task1_snapshot / "pathway/hallmark-w-2477x50.npy",
        "pathway_policy": task1_snapshot / "pathway/lincs-pathway-policy.json",
    }
    missing_lincs = [name for name, path in required_lincs_paths.items() if not path.is_file()]
    if missing_lincs:
        assertions.append(
            {
                "name": "task1_lincs_assets_present",
                "pass": False,
                "details": {"missing_assets": missing_lincs},
                "counterexamples": [{"missing_asset": name} for name in missing_lincs[:MAX_COUNTEREXAMPLES]],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Missing required Task1 LINCS assets: {missing_lincs}", file=sys.stderr)
        return 7

    input_paths.extend(required_lincs_paths.values())
    input_paths.extend((legacy_task2_root / rel).resolve() for rel in legacy_scperturb_required_relpaths())

    assertions.append(
        {
            "name": "task1_lincs_assets_present",
            "pass": True,
            "details": {name: str(path.resolve()) for name, path in required_lincs_paths.items()},
            "counterexamples": [],
        }
    )

    try:
        lincs_meta = pd.read_csv(required_lincs_paths["meta_csv"])
        gene_tensor, pathway_tensor = load_lincs_tensor_bundle(
            required_lincs_paths["tensor_pt"], expected_rows=len(lincs_meta)
        )
    except Exception as exc:
        assertions.append(
            {
                "name": "lincs_source_alignment",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed loading LINCS sources: {exc}", file=sys.stderr)
        return 8

    assertions.append(
        {
            "name": "lincs_source_alignment",
            "pass": True,
            "details": {
                "meta_rows": int(len(lincs_meta)),
                "gene_tensor_shape": [int(gene_tensor.shape[0]), int(gene_tensor.shape[1])],
                "pathway_tensor_shape": [int(pathway_tensor.shape[0]), int(pathway_tensor.shape[1])],
            },
            "counterexamples": [],
        }
    )

    try:
        supported_lincs_rows, excluded_lincs_rows = normalize_lincs_metadata(lincs_meta)
        (
            lincs_delta_meta,
            lincs_membership,
            lincs_coverage,
            lincs_pairs,
            lincs_exclusion_summary,
            lincs_summary,
        ) = build_lincs_task2_tables(supported_lincs_rows, excluded_lincs_rows)
    except Exception as exc:
        assertions.append(
            {
                "name": "lincs_task2_normalization_contract",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed normalizing LINCS Task2 rows: {exc}", file=sys.stderr)
        return 9

    lincs_cell_line_trace = pd.concat(
        [
            supported_lincs_rows[["cell_line_raw", "cell_line"]],
            excluded_lincs_rows[["cell_line_raw", "cell_line"]],
        ],
        ignore_index=True,
    )
    alias_rows = lincs_cell_line_trace.loc[
        lincs_cell_line_trace["cell_line_raw"] != lincs_cell_line_trace["cell_line"]
    ].copy()
    alias_counts = (
        alias_rows.groupby(["cell_line_raw", "cell_line"], sort=False)
        .size()
        .reset_index(name="n_rows")
    )
    assertions.append(
        {
            "name": "lincs_cell_line_normalization_rule",
            "pass": True,
            "details": {
                "rule_name": "explicit_task2_v2_lincs_cell_line_alias_map",
                "alias_map": LINCS_CELL_LINE_ALIAS_MAP,
                "trace_fields": ["cell_line_raw", "cell_line"],
                "n_rows_normalized": int(len(alias_rows)),
                "normalized_pairs": alias_counts.to_dict("records"),
            },
            "counterexamples": [],
        }
    )

    unsupported_count = int(
        excluded_lincs_rows["exclusion_reason"].eq("unsupported_pert_type").sum()
        if not excluded_lincs_rows.empty
        else 0
    )
    assertions.append(
        {
            "name": "lincs_perturbation_mapping_contract",
            "pass": True,
            "details": {
                "mapping": {
                    "drug": "Chemical",
                    "sh": "Genetic",
                    "crispr": "Genetic",
                    "oe": "Genetic",
                },
                "unsupported_pert_type_rows": unsupported_count,
                "pre_support_exclusions": lincs_summary["n_source_rows_excluded_pre_support"],
            },
            "counterexamples": []
            if excluded_lincs_rows.empty
            else excluded_lincs_rows.head(MAX_COUNTEREXAMPLES).to_dict("records"),
        }
    )

    if lincs_delta_meta["row_id"].duplicated().any():
        assertions.append(
            {
                "name": "lincs_row_identity_unique",
                "pass": False,
                "details": {"rules": ["LINCS delta_meta row_id must be unique."]},
                "counterexamples": lincs_delta_meta.loc[
                    lincs_delta_meta["row_id"].duplicated(), ["row_id", "source_row_index", "treated_cell_id"]
                ].head(MAX_COUNTEREXAMPLES).to_dict("records"),
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] LINCS row identity uniqueness failed.", file=sys.stderr)
        return 10

    assertions.append(
        {
            "name": "lincs_row_identity_unique",
            "pass": True,
            "details": {
                "row_identity_fields": ["row_id", "source_row_index", "treated_cell_id"],
                "n_rows": int(len(lincs_delta_meta)),
            },
            "counterexamples": [],
        }
    )

    snapshot_root = output_root
    snapshot_root.mkdir(parents=True, exist_ok=True)
    shared_pathway_dir = snapshot_root / "pathway"
    lincs_root = snapshot_root / "lincs"
    lincs_derived_dir = lincs_root / "derived"
    scperturb_root = snapshot_root / "scperturb_k562"
    scperturb_derived_dir = scperturb_root / "derived"

    shared_outputs: List[Path] = []
    lincs_outputs: List[Path] = []
    scperturb_outputs: List[Path] = []

    for name in ["pathway_w", "pathway_policy"]:
        src = required_lincs_paths[name]
        dst = shared_pathway_dir / src.name
        mode = link_or_copy(src.resolve(), dst)
        shared_outputs.append(dst.resolve())
        lineage_rows.append(
            {
                "dataset": "shared",
                "source": str(src.resolve()),
                "destination": str(dst.resolve()),
                "mode": mode,
            }
        )

    for name in ["meta_csv", "alignment_csv"]:
        src = required_lincs_paths[name]
        dst = lincs_root / src.name
        mode = link_or_copy(src.resolve(), dst)
        lincs_outputs.append(dst.resolve())
        lineage_rows.append(
            {
                "dataset": "LINCS",
                "source": str(src.resolve()),
                "destination": str(dst.resolve()),
                "mode": mode,
            }
        )

    tensor_pt_dst = lincs_root / required_lincs_paths["tensor_pt"].name
    tensor_mode = symlink_or_reuse(required_lincs_paths["tensor_pt"], tensor_pt_dst)
    lincs_outputs.append(tensor_pt_dst.absolute())
    lineage_rows.append(
        {
            "dataset": "LINCS",
            "source": str(required_lincs_paths["tensor_pt"].resolve()),
            "destination": local_output_path(tensor_pt_dst),
            "mode": tensor_mode,
        }
    )

    try:
        sc_outputs, sc_lineage_rows = materialize_legacy_scperturb_subtree(
            legacy_root=legacy_task2_root,
            output_root=snapshot_root,
        )
    except Exception as exc:
        assertions.append(
            {
                "name": "scperturb_subtree_materialization",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed materializing legacy scPerturb subtree: {exc}", file=sys.stderr)
        return 11

    scperturb_outputs.extend(sc_outputs)
    lineage_rows.extend(sc_lineage_rows)
    assertions.append(
        {
            "name": "scperturb_subtree_materialization",
            "pass": True,
            "details": {
                "n_files": len(sc_lineage_rows),
                "destination_root": str(scperturb_root.resolve()),
            },
            "counterexamples": [],
        }
    )

    sc_membership, sc_coverage = build_scperturb_membership(
        legacy_delta_meta_path=scperturb_derived_dir / "delta_meta.csv",
        common_targets=legacy_summary["common_targets"],
    )
    sc_membership_path = (scperturb_derived_dir / "task2_row_membership.parquet").resolve()
    write_parquet(sc_membership, sc_membership_path)
    scperturb_outputs.append(sc_membership_path)

    lincs_source_indices = lincs_delta_meta["source_row_index"].to_numpy(dtype=np.int64, copy=False)
    lincs_gene_delta_path = (lincs_derived_dir / "gene_delta.npy").resolve()
    lincs_pathway_delta_path = (lincs_derived_dir / "pathway_delta.npy").resolve()
    lincs_delta_meta_path = (lincs_derived_dir / "delta_meta.csv").resolve()
    lincs_membership_path = (lincs_derived_dir / "task2_row_membership.parquet").resolve()

    try:
        gene_stats = materialize_lincs_subset_tensor(gene_tensor, lincs_source_indices, lincs_gene_delta_path)
        pathway_stats = materialize_lincs_subset_tensor(
            pathway_tensor, lincs_source_indices, lincs_pathway_delta_path
        )
    except Exception as exc:
        assertions.append(
            {
                "name": "lincs_derived_materialization",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed materializing LINCS derived arrays: {exc}", file=sys.stderr)
        return 12

    write_csv(lincs_delta_meta, lincs_delta_meta_path)
    write_parquet(lincs_membership, lincs_membership_path)
    lincs_outputs.extend(
        [
            lincs_gene_delta_path,
            lincs_pathway_delta_path,
            lincs_delta_meta_path,
            lincs_membership_path,
        ]
    )
    assertions.append(
        {
            "name": "lincs_derived_materialization",
            "pass": True,
            "details": {
                "gene_delta_shape": [int(gene_stats[0]), int(gene_stats[1])],
                "pathway_delta_shape": [int(pathway_stats[0]), int(pathway_stats[1])],
                "delta_meta_rows": int(len(lincs_delta_meta)),
            },
            "counterexamples": [],
        }
    )

    lincs_pairs_path = (lincs_root / "task2_lincs_pairs.csv").resolve()
    write_csv(lincs_pairs, lincs_pairs_path)
    lincs_outputs.append(lincs_pairs_path)

    lincs_subtree_manifest_path = (lincs_root / "subtree_manifest.json").resolve()
    write_json(
        lincs_subtree_manifest_path,
        {
            "dataset": "LINCS",
            "analysis_key_fields": ["dataset", "cell_line", "target_token"],
            "row_identity_fields": ["row_id", "source_row_index", "treated_cell_id"],
            "cell_line_trace_fields": ["cell_line_raw", "cell_line"],
            "representation_scope": ["Gene", "Pathway"],
            "perturbation_mapping": {
                "drug": "Chemical",
                "sh": "Genetic",
                "crispr": "Genetic",
                "oe": "Genetic",
            },
            "cell_line_normalization": {
                "rule_name": "explicit_task2_v2_lincs_cell_line_alias_map",
                "alias_map": LINCS_CELL_LINE_ALIAS_MAP,
            },
            "source_assets": {
                "meta_csv": str((lincs_root / "lincs-engine1-meta.csv").resolve()),
                "alignment_csv": str((lincs_root / "lincs-gene-alignment.csv").resolve()),
                "tensor_pt": local_output_path(tensor_pt_dst),
            },
            "summary": lincs_summary,
            "derived_outputs": {
                "gene_delta": str(lincs_gene_delta_path),
                "pathway_delta": str(lincs_pathway_delta_path),
                "delta_meta": str(lincs_delta_meta_path),
                "membership": str(lincs_membership_path),
                "task2_lincs_pairs": str(lincs_pairs_path),
            },
        },
    )
    lincs_outputs.append(lincs_subtree_manifest_path)

    scperturb_subtree_manifest_path = (scperturb_root / "subtree_manifest.json").resolve()
    write_json(
        scperturb_subtree_manifest_path,
        {
            "dataset": "scPerturb",
            "cell_line": "K562",
            "analysis_key_fields": ["dataset", "cell_line", "target_token"],
            "row_identity_fields": ["row_id", "treated_cell_id"],
            "representation_scope": list(REPRESENTATION_ORDER),
            "source_root": str((legacy_task2_root / "k562").resolve()),
            "common_targets": legacy_summary["common_targets"],
            "summary": {
                "delta_meta_rows": legacy_summary["delta_meta_rows"],
                "gene_dim": legacy_summary["gene_dim"],
                "pathway_dim": legacy_summary["pathway_dim"],
                "fm": legacy_summary["fm_summaries"],
            },
            "membership_output": str(sc_membership_path),
        },
    )
    scperturb_outputs.append(scperturb_subtree_manifest_path)

    combined_coverage = pd.concat([lincs_coverage, sc_coverage], ignore_index=True, sort=False)
    dataset_order = {"LINCS": 0, "scPerturb": 1}
    sc_target_order = {token: idx for idx, token in enumerate(legacy_summary["common_targets"])}
    combined_coverage["_dataset_order"] = combined_coverage["dataset"].map(dataset_order)
    combined_coverage["_target_order"] = combined_coverage.apply(
        lambda row: sc_target_order.get(str(row["target_token"]), 10_000)
        if row["dataset"] == "scPerturb"
        else -1,
        axis=1,
    )
    combined_coverage = combined_coverage.sort_values(
        ["_dataset_order", "cell_line", "_target_order", "target_token"],
        kind="mergesort",
    ).drop(columns=["_dataset_order", "_target_order"])
    combined_coverage.reset_index(drop=True, inplace=True)

    coverage_path = (snapshot_root / "task2_pairs_coverage.csv").resolve()
    registry = build_representation_availability_registry(combined_coverage)
    registry_path = (snapshot_root / "representation_availability_registry.csv").resolve()
    write_csv(combined_coverage, coverage_path)
    write_csv(registry, registry_path)
    shared_outputs.extend([coverage_path, registry_path])

    inventory = build_post_build_inventory(
        scperturb_delta_meta_path=scperturb_derived_dir / "delta_meta.csv",
        scperturb_fm_root=scperturb_root / "fm",
        lincs_supported_rows=supported_lincs_rows,
        lincs_included_delta_meta=lincs_delta_meta,
    )
    inventory_path = (stage_dir / "task2_post_build_inventory.csv").resolve()
    write_csv(inventory, inventory_path)

    lingering_alias_rows: List[Dict[str, object]] = []
    for output_name, frame in [
        ("task2_pairs_coverage.csv", combined_coverage.loc[combined_coverage["dataset"].eq("LINCS")]),
        ("lincs/derived/delta_meta.csv", lincs_delta_meta),
        ("representation_availability_registry.csv", registry.loc[registry["dataset"].eq("LINCS")]),
        ("task2_post_build_inventory.csv", inventory.loc[inventory["dataset"].eq("LINCS")]),
    ]:
        for raw_alias in LINCS_CELL_LINE_ALIAS_MAP:
            alias_frame = frame.loc[frame["cell_line"].astype(str).eq(raw_alias)]
            if alias_frame.empty:
                continue
            lingering_alias_rows.append(
                {
                    "output": output_name,
                    "cell_line": raw_alias,
                    "n_rows": int(len(alias_frame)),
                }
            )
    assertions.append(
        {
            "name": "lincs_cell_line_normalization_applied",
            "pass": len(lingering_alias_rows) == 0,
            "details": {
                "checked_outputs": [
                    "task2_pairs_coverage.csv",
                    "lincs/derived/delta_meta.csv",
                    "representation_availability_registry.csv",
                    "task2_post_build_inventory.csv",
                ],
                "forbidden_normalized_labels": sorted(LINCS_CELL_LINE_ALIAS_MAP.keys()),
            },
            "counterexamples": lingering_alias_rows[:MAX_COUNTEREXAMPLES],
        }
    )
    if lingering_alias_rows:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] LINCS cell-line normalization was not fully applied.", file=sys.stderr)
        return 13

    snapshot_manifest_path = (snapshot_root / "snapshot_manifest.json").resolve()
    source_lineage_manifest_path = (snapshot_root / "source_lineage_manifest.json").resolve()

    write_json(
        snapshot_manifest_path,
        {
            "snapshot_root": str(snapshot_root),
            "stage": STAGE,
            "analysis_key_fields": ["dataset", "cell_line", "target_token"],
            "row_identity_contract": {
                "scPerturb": ["row_id", "treated_cell_id"],
                "LINCS": ["row_id", "source_row_index", "treated_cell_id"],
            },
            "cell_line_normalization": {
                "LINCS": {
                    "rule_name": "explicit_task2_v2_lincs_cell_line_alias_map",
                    "alias_map": LINCS_CELL_LINE_ALIAS_MAP,
                    "trace_fields": ["cell_line_raw", "cell_line"],
                }
            },
            "target_membership_source": "derived/delta_meta.csv::target_tokens",
            "datasets": {
                "scPerturb": {
                    "cell_lines": ["K562"],
                    "subtree": str(scperturb_root.resolve()),
                    "representations": list(REPRESENTATION_ORDER),
                },
                "LINCS": {
                    "cell_lines": sorted(lincs_pairs["cell_line"].drop_duplicates().astype(str).tolist()),
                    "subtree": str(lincs_root.resolve()),
                    "representations": ["Gene", "Pathway"],
                },
            },
            "counts": {
                "scperturb_k562_rows": int(legacy_summary["delta_meta_rows"]),
                "lincs_rows": int(len(lincs_delta_meta)),
                "total_coverage_rows": int(len(combined_coverage)),
                "eligible_lincs_cohorts": int(len(lincs_pairs)),
            },
            "shared_artifacts": {
                "coverage": str(coverage_path),
                "representation_registry": str(registry_path),
                "pathway_dir": str(shared_pathway_dir.resolve()),
            },
        },
    )
    write_json(
        source_lineage_manifest_path,
        {
            "legacy_task2_snapshot": str(legacy_task2_root),
            "task1_snapshot": str(task1_snapshot),
            "cell_line_normalization": {
                "LINCS": {
                    "rule_name": "explicit_task2_v2_lincs_cell_line_alias_map",
                    "alias_map": LINCS_CELL_LINE_ALIAS_MAP,
                }
            },
            "materialization": lineage_rows,
        },
    )
    shared_outputs.extend([snapshot_manifest_path, source_lineage_manifest_path])

    output_paths.extend(shared_outputs)
    output_paths.extend(lincs_outputs)
    output_paths.extend(scperturb_outputs)
    output_paths.append(inventory_path)

    output_routing_snapshot_pass = all(
        path.absolute().is_relative_to(snapshot_root) or path.resolve().is_relative_to(snapshot_root)
        for path in output_paths
        if path != inventory_path
    )
    output_routing_stage_pass = inventory_path.absolute().is_relative_to(stage_dir.resolve())
    assertions.append(
        {
            "name": "output_routing_isolated",
            "pass": bool(output_routing_snapshot_pass and output_routing_stage_pass),
            "details": {
                "snapshot_root": str(snapshot_root),
                "stage_dir": str(stage_dir),
            },
            "counterexamples": []
            if output_routing_snapshot_pass and output_routing_stage_pass
            else [{"bad_output": str(path)} for path in output_paths[:MAX_COUNTEREXAMPLES]],
        }
    )

    logical_allowed_input_roots = [task1_snapshot.resolve(), legacy_task2_root.resolve()]
    resolved_allowed_input_roots = sorted(
        {
            path.resolve().parent
            for path in input_paths
            if not any(path.resolve().is_relative_to(root) for root in logical_allowed_input_roots)
        },
        key=lambda path: str(path),
    )
    allowed_input_roots = logical_allowed_input_roots + resolved_allowed_input_roots
    bad_inputs = [
        str(path)
        for path in sorted(set(input_paths))
        if not any(
            path.absolute().is_relative_to(root) or path.resolve().is_relative_to(root)
            for root in allowed_input_roots
        )
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "rules": [
                    "Legacy scPerturb reuse inputs are read from the legacy Task2 snapshot root.",
                    "Task1-derived LINCS inputs may be consumed via logical task1_snapshot paths or their resolved external LINCS_Processed source roots.",
                ],
                "allowed_input_logical_roots": [str(path) for path in logical_allowed_input_roots],
                "allowed_input_resolved_roots": [str(path) for path in resolved_allowed_input_roots],
                "n_inputs": len(sorted(set(input_paths))),
            },
            "counterexamples": [{"bad_input": path} for path in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    eligible_lines = sorted(lincs_pairs["cell_line"].drop_duplicates().astype(str).tolist())
    if registry.empty or set(eligible_lines) - set(
        registry.loc[registry["dataset"].eq("LINCS"), "cell_line"].drop_duplicates().astype(str).tolist()
    ):
        assertions.append(
            {
                "name": "representation_availability_policy",
                "pass": False,
                "details": {"error": "Representation registry missing eligible LINCS cell lines."},
                "counterexamples": [],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Representation availability registry is incomplete.", file=sys.stderr)
        return 14

    assertions.append(
        {
            "name": "representation_availability_policy",
            "pass": True,
            "details": {
                "lincs_allowed": sorted(LINCS_REPRESENTATIONS),
                "scperturb_k562_allowed": list(REPRESENTATION_ORDER),
                "eligible_lincs_cell_lines": eligible_lines,
            },
            "counterexamples": [],
        }
    )

    completed_at = utc_now_iso()
    run_manifest_path = (stage_dir / "run_manifest.json").resolve()
    audit_assertions_path = (stage_dir / "audit_assertions.json").resolve()
    manifest_path = (stage_dir / "manifest.json").resolve()

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s3_build_task2_multisource_snapshot.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "task1_snapshot": str(task1_snapshot),
            "legacy_task2_snapshot": str(legacy_task2_root),
            "output_root": str(snapshot_root),
            "runs_dir": str(runs_dir),
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths))],
        "outputs": [
            local_output_path(path)
            for path in output_paths + [run_manifest_path, audit_assertions_path, manifest_path]
        ],
        "scperturb_reuse_summary": {
            "source_root": str((legacy_task2_root / "k562").resolve()),
            "destination_root": str(scperturb_root.resolve()),
            "delta_meta_rows": legacy_summary["delta_meta_rows"],
            "common_targets": legacy_summary["common_targets"],
            "fm": legacy_summary["fm_summaries"],
        },
        "lincs_build_summary": {
            **lincs_summary,
            "cell_line_normalization": {
                "rule_name": "explicit_task2_v2_lincs_cell_line_alias_map",
                "alias_map": LINCS_CELL_LINE_ALIAS_MAP,
                "trace_fields": ["cell_line_raw", "cell_line"],
            },
            "gene_delta_shape": [int(gene_stats[0]), int(gene_stats[1])],
            "pathway_delta_shape": [int(pathway_stats[0]), int(pathway_stats[1])],
        },
    }

    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})
    write_json(manifest_path, {"stage": STAGE, "files": build_stage_manifest(stage_dir)})

    if not (output_routing_snapshot_pass and output_routing_stage_pass):
        print("[ERROR] Output routing assertion failed.", file=sys.stderr)
        return 15
    if bad_inputs:
        print("[ERROR] Input path isolation assertion failed.", file=sys.stderr)
        return 16

    return 0


if __name__ == "__main__":
    sys.exit(main())
