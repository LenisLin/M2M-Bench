# SCRIPT_HEADER_CONTRACT
# Script: scripts/fm_extractors/extract_tahoex1.py
# Purpose: Build Task2 K562 Tahoe-x1 FM deltas using pair_list with strict row alignment to delta_meta.
# Inputs:
#   - Task2 Snapshot K562: config/config.yaml::paths.task2_snapshot
#     - k562/derived/pair_list.parquet
#     - k562/derived/delta_meta.csv
#     - k562/{CRISPR_counts.pt,CRISPR_meta.csv,Drug_counts.pt,Drug_meta.csv,shared_var_names.csv}
#   - Tahoe-x1 model/code dir:
#     - config/config.yaml::fm_extractors.{tahoex1|tahoe_x1|tahoe-x1}.model_dir (optional)
#     - CLI override: --model-dir
# Outputs:
#   - fm_delta.npy: <paths.task2_snapshot>/k562/fm/tahoe-x1/
#   - fm_delta_meta.csv: <paths.task2_snapshot>/k562/fm/tahoe-x1/
#   - delta_operator_policy.json: <paths.task2_snapshot>/k562/fm/tahoe-x1/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
#     under runs/<run_id>/extract_tahoex1/
# Side Effects:
#   - Creates FM output directory under task2 snapshot
#   - Creates isolated run directory: runs/<run_id>/extract_tahoex1/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::fm_extractors.{tahoex1|tahoe_x1|tahoe-x1} (optional)
# Execution:
#   - python scripts/fm_extractors/extract_tahoex1.py --run-id <run_id> [--model-dir <dir>]
# Failure Modes:
#   - Missing required snapshot files -> exit non-zero
#   - Missing/invalid Tahoe-x1 assets -> exit non-zero
#   - Complete embedding failure for a side -> exit non-zero
# Last Updated: 2026-03-05

"""
Task2 Tahoe-x1 FM extractor with strict row-preservation contract:

1) Output fm_delta.npy must have exactly N rows where N = len(delta_meta).
2) Row i always corresponds to delta_meta row_id i (contiguous 0..N-1 required).
3) Any per-row embedding failure does NOT drop rows:
   - fm_delta[row_id, :] = NaN
   - fm_delta_meta.valid_mask[row_id] = False
4) Output routing is fixed:
   - Snapshot outputs -> data/task2_snapshot_v1/k562/fm/tahoe-x1/
   - Run artifacts -> runs/<run_id>/extract_tahoex1/
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import pickle
import random
import shlex
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import yaml

STAGE = "extract_tahoex1"
MODEL_NAME = "tahoe-x1"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
GLOBAL_SEED = 619
MAX_COUNTEREXAMPLES = 5
MIN_MAPPED_GENES = 500
DEFAULT_DOCKER_IMAGE = "ghcr.io/tahoebio/tahoe-x1:latest"
DEFAULT_DOCKER_SHM_SIZE = "64g"
CONTAINER_CODE_ROOT = PurePosixPath("/workspace/tahoe-x1")
CONTAINER_MODEL_DIR = PurePosixPath("/workspace/model-package")
CONTAINER_JOB_DIR = PurePosixPath("/workspace/job")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Task2 K562 Tahoe-x1 FM delta extractor")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="Override Tahoe-x1 model/code directory configured in config/config.yaml",
    )
    parser.add_argument("--python-bin", type=str, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--cell-chunk-size", type=int, default=None)
    parser.add_argument("--seq-len-dataset", type=int, default=None)
    parser.add_argument("--num-workers", type=int, default=None)
    parser.add_argument("--prefetch-factor", type=int, default=None)
    parser.add_argument("--model-size", type=str, default=None, choices=["70m", "1b", "3b"])
    parser.add_argument("--tx-model-name", type=str, default=None)
    parser.add_argument("--gene-map-pkl", type=Path, default=None)
    parser.add_argument("--cell-type-key", type=str, default=None)
    parser.add_argument("--gene-id-key", type=str, default=None)
    parser.add_argument("--default-cell-type", type=str, default=None)
    parser.add_argument("--docker-image", type=str, default=None)
    parser.add_argument("--docker-shm-size", type=str, default=None)
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


def init_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def choose_side_id_column(meta: pd.DataFrame, side: str) -> str:
    if side == "CRISPR":
        candidates = ["cell_barcode", "cell_id", "Unnamed: 0"]
    else:
        candidates = ["Unnamed: 0", "cell_id", "cell_barcode"]
    for column in candidates:
        if column in meta.columns:
            return column
    raise ValueError(f"Unable to locate cell-id column for side={side}")


def is_oom_text(text: str) -> bool:
    t = text.lower()
    return (
        "out of memory" in t
        or "cuda out of memory" in t
        or ("cuda" in t and "memory" in t)
        or "cublas" in t
    )


def pick_tahoex1_cfg(config: Mapping[str, object]) -> Tuple[Mapping[str, object], str]:
    fm_cfg_root = config.get("fm_extractors", {})
    if not isinstance(fm_cfg_root, Mapping):
        return {}, "missing"

    key_candidates = [MODEL_NAME, "tahoex1", "tahoe_x1", "tahoe-x1"]
    for key in key_candidates:
        node = fm_cfg_root.get(key)
        if isinstance(node, Mapping):
            return node, f"fm_extractors.{key}"
    return {}, "missing"


def has_tahoex1_package_files(path: Path) -> bool:
    required = ["model_config.yml", "collator_config.yml", "vocab.json", "best-model.pt"]
    return path.is_dir() and all((path / name).is_file() for name in required)


def resolve_tahoex1_assets(
    *,
    project_root: Path,
    model_dir: Path,
    cfg: Mapping[str, object],
    args: argparse.Namespace,
    model_size: str,
) -> Tuple[Path, Path, Path, Path, Dict[str, str]]:
    sources: Dict[str, str] = {}
    root = model_dir.resolve()

    predict_candidates: List[Tuple[Path, str]] = []
    if "predict_script" in cfg and cfg["predict_script"] is not None:
        predict_candidates.append(
            (resolve_config_path(project_root, str(cfg["predict_script"])), "config.predict_script")
        )
    predict_candidates.extend(
        [
            (root / "scripts/inference/predict_embeddings.py", "model_dir/scripts/inference/predict_embeddings.py"),
            (root / "predict_embeddings.py", "model_dir/predict_embeddings.py"),
            (root.parent / "tahoe-x1/scripts/inference/predict_embeddings.py", "model_dir_parent/tahoe-x1/scripts/inference/predict_embeddings.py"),
            (root.parent / "Tahoe-x1/scripts/inference/predict_embeddings.py", "model_dir_parent/Tahoe-x1/scripts/inference/predict_embeddings.py"),
        ]
    )
    for anc in [root] + list(root.parents)[:4]:
        predict_candidates.extend(
            [
                (anc / "benchmark/tahoe-x1/scripts/inference/predict_embeddings.py", f"{anc}/benchmark/tahoe-x1/scripts/inference/predict_embeddings.py"),
                (anc / "M2MBench/benchmark/tahoe-x1/scripts/inference/predict_embeddings.py", f"{anc}/M2MBench/benchmark/tahoe-x1/scripts/inference/predict_embeddings.py"),
            ]
        )

    predict_script: Optional[Path] = None
    for p, source in predict_candidates:
        if p.is_file():
            predict_script = p.resolve()
            sources["predict_script"] = source
            break
    if predict_script is None:
        raise RuntimeError("Cannot resolve Tahoe-x1 predict_embeddings.py from --model-dir")

    code_root = predict_script.parent.parent.parent.resolve()

    package_candidates: List[Tuple[Path, str]] = []
    if "model_package_dir" in cfg and cfg["model_package_dir"] is not None:
        package_candidates.append(
            (resolve_config_path(project_root, str(cfg["model_package_dir"])), "config.model_package_dir")
        )
    package_candidates.extend(
        [
            (root, "model_dir"),
            (root / f"{model_size}-model", f"model_dir/{model_size}-model"),
            (root.parent / f"{model_size}-model", f"model_dir_parent/{model_size}-model"),
            (root.parent / "Tahoe-x1" / f"{model_size}-model", f"model_dir_parent/Tahoe-x1/{model_size}-model"),
            (root.parent / "tahoe-x1" / f"{model_size}-model", f"model_dir_parent/tahoe-x1/{model_size}-model"),
            (code_root / f"{model_size}-model", f"code_root/{model_size}-model"),
            (code_root.parent / "Tahoe-x1" / f"{model_size}-model", f"code_root_parent/Tahoe-x1/{model_size}-model"),
            (code_root.parent / "tahoe-x1" / f"{model_size}-model", f"code_root_parent/tahoe-x1/{model_size}-model"),
        ]
    )

    model_package_dir: Optional[Path] = None
    for p, source in package_candidates:
        if has_tahoex1_package_files(p):
            model_package_dir = p.resolve()
            sources["model_package_dir"] = source
            break
    if model_package_dir is None:
        raise RuntimeError(
            "Cannot resolve Tahoe-x1 package dir with model_config.yml/collator_config.yml/vocab.json/best-model.pt"
        )

    gene_map_candidates: List[Tuple[Path, str]] = []
    if args.gene_map_pkl is not None:
        gene_map_candidates.append((args.gene_map_pkl.resolve(), "cli(--gene-map-pkl)"))
    if "gene_map_pkl" in cfg and cfg["gene_map_pkl"] is not None:
        gene_map_candidates.append(
            (resolve_config_path(project_root, str(cfg["gene_map_pkl"])), "config.gene_map_pkl")
        )
    gene_map_candidates.extend(
        [
            (
                root.parent / "Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
                "model_dir_parent/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
            ),
            (
                code_root.parent / "Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
                "code_root_parent/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
            ),
            (
                project_root.parent / "M2MBench/benchmark/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
                "project_root_parent/M2MBench/benchmark/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
            ),
        ]
    )
    for anc in [root] + list(root.parents)[:4]:
        gene_map_candidates.extend(
            [
                (
                    anc / "benchmark/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
                    f"{anc}/benchmark/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
                ),
                (
                    anc / "Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
                    f"{anc}/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl",
                ),
            ]
        )

    gene_map_pkl: Optional[Path] = None
    for p, source in gene_map_candidates:
        if p.is_file():
            gene_map_pkl = p.resolve()
            sources["gene_map_pkl"] = source
            break
    if gene_map_pkl is None:
        raise RuntimeError(
            "Cannot resolve gene_name_id_dict_gc104M.pkl. Provide --gene-map-pkl or config gene_map_pkl."
        )

    return code_root, predict_script, model_package_dir, gene_map_pkl, sources


def load_symbol_to_ensembl(pkl_path: Path) -> Dict[str, str]:
    with pkl_path.open("rb") as handle:
        obj = pickle.load(handle)
    if not isinstance(obj, dict):
        raise TypeError(f"Expected dict in {pkl_path}, got {type(obj)}")
    output: Dict[str, str] = {}
    for k, v in obj.items():
        if k is None or v is None:
            continue
        key = str(k).strip()
        val = str(v).strip()
        if key and val:
            output[key] = val
    return output


def build_tahoex1_gene_mapper(
    shared_genes: Sequence[str],
    symbol2ens: Mapping[str, str],
) -> Tuple[np.ndarray, List[str], List[str]]:
    src_cols: List[int] = []
    mapped_symbols: List[str] = []
    mapped_ensembl: List[str] = []

    for i, symbol in enumerate(shared_genes):
        sym = str(symbol)
        ens = symbol2ens.get(sym)
        if ens is None:
            continue
        src_cols.append(int(i))
        mapped_symbols.append(sym)
        mapped_ensembl.append(str(ens))

    src_arr = np.asarray(src_cols, dtype=np.int64)
    if src_arr.size < MIN_MAPPED_GENES:
        raise RuntimeError(
            f"Too few genes mapped to Ensembl IDs: {src_arr.size} < {MIN_MAPPED_GENES}."
        )

    return src_arr, mapped_symbols, mapped_ensembl


def write_tahoex1_input_h5ad(
    *,
    counts_dense: np.ndarray,
    mapped_gene_symbols: Sequence[str],
    mapped_ensembl: Sequence[str],
    gene_id_key: str,
    cell_ids: Sequence[str],
    cell_type_key: str,
    default_cell_type: str,
    out_h5ad: Path,
) -> None:
    try:
        import anndata as ad  # type: ignore
        from scipy import sparse  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("anndata/scipy is required for Tahoe-x1 extractor") from exc

    x = counts_dense.astype(np.float32, copy=False)
    x_csr = sparse.csr_matrix(x)

    obs = pd.DataFrame(index=pd.Index([str(cid) for cid in cell_ids], name="cell_id"))
    obs[cell_type_key] = [str(default_cell_type)] * len(cell_ids)

    var = pd.DataFrame(index=pd.Index([str(g) for g in mapped_gene_symbols], name="gene_symbol"))
    var[gene_id_key] = [str(e) for e in mapped_ensembl]

    adata = ad.AnnData(X=x_csr, obs=obs, var=var)
    adata.var_names_make_unique()
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)


def write_tahoex1_predict_yaml(
    *,
    cfg_yaml: Path,
    model_package_dir: Path,
    adata_input: Path,
    adata_output: Path,
    tx_model_name: str,
    cell_type_key: str,
    gene_id_key: str,
    batch_size: int,
    seq_len_dataset: int,
    num_workers: int,
    prefetch_factor: int,
) -> None:
    payload = {
        "model_name": tx_model_name,
        "paths": {
            "model_dir": str(model_package_dir),
            "adata_input": str(adata_input),
            "adata_output": str(adata_output),
        },
        "data": {
            "cell_type_key": str(cell_type_key),
            "gene_id_key": str(gene_id_key),
        },
        "predict": {
            "batch_size": int(batch_size),
            "seq_len_dataset": int(seq_len_dataset),
            "num_workers": int(num_workers),
            "prefetch_factor": int(prefetch_factor),
            "return_gene_embeddings": False,
            "use_chem_inf": False,
        },
    }

    cfg_yaml.parent.mkdir(parents=True, exist_ok=True)
    with cfg_yaml.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def build_container_job_path(host_path: Path) -> PurePosixPath:
    return CONTAINER_JOB_DIR / host_path.name


def run_tahoex1_once(
    *,
    python_bin: str,
    predict_script: Path,
    model_root: Path,
    model_package_dir: Path,
    docker_image: str,
    docker_shm_size: str,
    job_dir: Path,
    cfg_yaml: Path,
    batch_size: int,
) -> Tuple[bool, str]:
    container_cfg_yaml = build_container_job_path(cfg_yaml)
    shell_cmd = (
        "python -m pip install -e . --no-deps >/tmp/tahoex1_pip_install.log 2>&1"
        f" && {shlex.quote(python_bin)} scripts/inference/predict_embeddings.py "
        f"{shlex.quote(str(container_cfg_yaml))} --predict.batch_size={int(batch_size)}"
    )
    cmd = [
        "docker",
        "run",
        "--rm",
        "--gpus",
        "all",
        "--shm-size",
        str(docker_shm_size),
        "-v",
        f"{model_root}:{CONTAINER_CODE_ROOT}",
        "-v",
        f"{model_package_dir}:{CONTAINER_MODEL_DIR}",
        "-v",
        f"{job_dir}:{CONTAINER_JOB_DIR}",
        "-w",
        str(CONTAINER_CODE_ROOT),
        str(docker_image),
        "/bin/bash",
        "-lc",
        shell_cmd,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)
    if result.returncode == 0:
        return True, ""

    msg = (result.stderr or "") + "\n" + (result.stdout or "")
    return False, msg[-3000:]


def read_h5ad(path: Path) -> Any:
    try:
        import scanpy as sc  # type: ignore

        return sc.read_h5ad(str(path))
    except Exception:
        try:
            import anndata as ad  # type: ignore

            return ad.read_h5ad(str(path))
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError("Unable to read h5ad (scanpy/anndata missing).") from exc


def extract_tahoex1_aligned(
    *,
    output_h5ad: Path,
    expected_ids: Sequence[str],
    tx_model_name: str,
) -> Tuple[np.ndarray, List[str], str]:
    adata = read_h5ad(output_h5ad)
    keys = list(adata.obsm_keys())

    chosen = None
    if tx_model_name in keys:
        chosen = tx_model_name
    else:
        candidates = [k for k in keys if "tx1" in str(k).lower() or "tahoe" in str(k).lower()]
        if candidates:
            chosen = candidates[0]

    if chosen is None:
        raise KeyError(f"Cannot find Tahoe-x1 embedding key. expected={tx_model_name}, available={keys}")

    emb = np.asarray(adata.obsm[chosen], dtype=np.float32)
    if emb.ndim != 2:
        raise RuntimeError(f"Tahoe-x1 embedding matrix must be 2D, got shape={emb.shape}")

    obs_ids = adata.obs_names.astype(str).tolist()
    emb_df = pd.DataFrame(emb, index=pd.Index(obs_ids, name="cell_id"))
    aligned = emb_df.reindex([str(x) for x in expected_ids])

    missing_mask = aligned.isna().any(axis=1)
    missing_ids = aligned.index[missing_mask].astype(str).tolist()
    matrix = aligned.to_numpy(dtype=np.float32, copy=False)
    return matrix, missing_ids, chosen


def run_tahoex1_segment_with_retry(
    *,
    python_bin: str,
    predict_script: Path,
    model_root: Path,
    model_package_dir: Path,
    docker_image: str,
    docker_shm_size: str,
    job_dir: Path,
    cfg_yaml: Path,
    output_h5ad: Path,
    initial_batch_size: int,
    tx_model_name: str,
    expected_ids: Sequence[str],
) -> Tuple[np.ndarray, int, Optional[str], Optional[str]]:
    batch_size = max(1, int(initial_batch_size))
    last_err: Optional[str] = None

    while True:
        ok, msg = run_tahoex1_once(
            python_bin=python_bin,
            predict_script=predict_script,
            model_root=model_root,
            model_package_dir=model_package_dir,
            docker_image=docker_image,
            docker_shm_size=docker_shm_size,
            job_dir=job_dir,
            cfg_yaml=cfg_yaml,
            batch_size=batch_size,
        )
        if ok:
            matrix, _missing, used_key = extract_tahoex1_aligned(
                output_h5ad=output_h5ad,
                expected_ids=expected_ids,
                tx_model_name=tx_model_name,
            )
            return matrix, batch_size, used_key, None

        last_err = msg
        if batch_size <= 1 or not is_oom_text(msg):
            return np.empty((0, 0), dtype=np.float32), batch_size, None, last_err

        batch_size = max(1, batch_size // 2)
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


def extract_side_embeddings_tahoex1(
    *,
    side: str,
    meta_path: Path,
    counts_path: Path,
    required_cell_ids: Sequence[str],
    src_cols: np.ndarray,
    mapped_gene_symbols: Sequence[str],
    mapped_ensembl: Sequence[str],
    python_bin: str,
    predict_script: Path,
    model_root: Path,
    model_package_dir: Path,
    docker_image: str,
    docker_shm_size: str,
    tx_model_name: str,
    batch_size: int,
    seq_len_dataset: int,
    num_workers: int,
    prefetch_factor: int,
    cell_chunk_size: int,
    cell_type_key: str,
    gene_id_key: str,
    default_cell_type: str,
    stage_dir: Path,
) -> Tuple[Dict[str, np.ndarray], Sequence[str], int, Dict[str, int], List[str]]:
    meta = pd.read_csv(meta_path)
    id_col = choose_side_id_column(meta, side)
    meta_ids = meta[id_col].map(normalize_scalar).astype(str)

    if meta_ids.duplicated().any():
        dup = meta_ids.loc[meta_ids.duplicated()].head(MAX_COUNTEREXAMPLES).tolist()
        raise ValueError(f"{side} metadata id column has duplicates, examples={dup}")

    id_to_row: Dict[str, int] = {cell_id: int(i) for i, cell_id in enumerate(meta_ids.tolist())}

    counts_t = torch.load(counts_path, map_location="cpu")
    if not torch.is_tensor(counts_t):
        raise ValueError(f"{side} counts file must be torch.Tensor, got {type(counts_t)}")
    if int(counts_t.shape[0]) != len(meta):
        raise ValueError(
            f"{side} counts/meta row mismatch: counts={int(counts_t.shape[0])}, meta={len(meta)}"
        )

    required_unique = sorted(set(str(x) for x in required_cell_ids))
    missing_ids = [cell_id for cell_id in required_unique if cell_id not in id_to_row]
    present_ids = [cell_id for cell_id in required_unique if cell_id in id_to_row]

    present_rows = np.asarray([id_to_row[cell_id] for cell_id in present_ids], dtype=np.int64)
    order = np.argsort(present_rows, kind="mergesort")
    ordered_rows = present_rows[order]
    ordered_ids = [present_ids[int(i)] for i in order.tolist()]

    vectors_by_cell: Dict[str, np.ndarray] = {}
    invalid_cells: set[str] = set(missing_ids)
    emb_dim: Optional[int] = None
    chunk_errors: List[str] = []

    stats = {
        "n_required_cells": int(len(required_unique)),
        "n_missing_in_meta": int(len(missing_ids)),
        "n_present_in_meta": int(len(ordered_ids)),
        "n_chunks_total": 0,
        "n_chunks_failed": 0,
        "n_chunks_succeeded": 0,
        "n_rows_filtered_zero_after_mapping": 0,
    }

    chunk_n = max(1, int(cell_chunk_size))
    for start in range(0, len(ordered_rows), chunk_n):
        end = min(start + chunk_n, len(ordered_rows))
        chunk_rows = ordered_rows[start:end]
        chunk_ids = ordered_ids[start:end]
        stats["n_chunks_total"] += 1

        queue: List[np.ndarray] = [np.arange(len(chunk_rows), dtype=np.int64)]
        segment_i = 0
        chunk_success = True

        while queue:
            seg_local = queue.pop(0)
            seg_ids = [chunk_ids[int(i)] for i in seg_local.tolist()]
            seg_rows = chunk_rows[seg_local]

            seg_tensor = torch.as_tensor(seg_rows, dtype=torch.long)
            seg_counts = counts_t.index_select(0, seg_tensor).detach().cpu().numpy().astype(np.float32, copy=False)

            try:
                seg_mapped = seg_counts[:, src_cols]
            except Exception as exc:  # noqa: BLE001
                seg_mapped = None
                map_exc = exc
            else:
                map_exc = None

            if seg_mapped is None:
                if len(seg_local) > 1:
                    mid = len(seg_local) // 2
                    queue.insert(0, seg_local[mid:])
                    queue.insert(0, seg_local[:mid])
                    chunk_errors.append(
                        f"{side} chunk[{start}:{end}] segment split len={len(seg_local)} due to map error: {map_exc}"
                    )
                else:
                    invalid_cells.add(seg_ids[0])
                    chunk_errors.append(
                        f"{side} chunk[{start}:{end}] segment single-row map failure cell={seg_ids[0]} error: {map_exc}"
                    )
                    chunk_success = False
                del seg_counts
                gc.collect()
                continue

            row_nonzero = (seg_mapped > 0).sum(axis=1) > 0
            keep_idx = np.where(row_nonzero)[0].astype(np.int64)
            drop_idx = np.where(~row_nonzero)[0].astype(np.int64)

            if len(drop_idx) > 0:
                stats["n_rows_filtered_zero_after_mapping"] += int(len(drop_idx))
                for i in drop_idx.tolist():
                    invalid_cells.add(seg_ids[int(i)])

            if len(keep_idx) == 0:
                chunk_errors.append(
                    f"{side} chunk[{start}:{end}] segment all rows zero after Tahoe gene mapping"
                )
                chunk_success = False
                del seg_counts
                del seg_mapped
                gc.collect()
                continue

            kept_ids = [seg_ids[int(i)] for i in keep_idx.tolist()]
            kept_counts = seg_mapped[keep_idx]

            with tempfile.TemporaryDirectory(prefix=f"tahoex1_{side}_{start}_{end}_seg{segment_i:04d}_", dir=str(stage_dir)) as tmp:
                tmp_dir = Path(tmp)
                input_h5ad = tmp_dir / f"chunk_{side}_{start}_{end}_seg{segment_i:04d}.h5ad"
                output_h5ad = tmp_dir / f"chunk_{side}_{start}_{end}_seg{segment_i:04d}_tx1.h5ad"
                cfg_yaml = tmp_dir / f"chunk_{side}_{start}_{end}_seg{segment_i:04d}.yaml"
                segment_i += 1
                container_input_h5ad = build_container_job_path(input_h5ad)
                container_output_h5ad = build_container_job_path(output_h5ad)
                container_cfg_yaml = build_container_job_path(cfg_yaml)

                try:
                    write_tahoex1_input_h5ad(
                        counts_dense=kept_counts,
                        mapped_gene_symbols=mapped_gene_symbols,
                        mapped_ensembl=mapped_ensembl,
                        gene_id_key=gene_id_key,
                        cell_ids=kept_ids,
                        cell_type_key=cell_type_key,
                        default_cell_type=default_cell_type,
                        out_h5ad=input_h5ad,
                    )

                    write_tahoex1_predict_yaml(
                        cfg_yaml=cfg_yaml,
                        model_package_dir=Path(str(CONTAINER_MODEL_DIR)),
                        adata_input=Path(str(container_input_h5ad)),
                        adata_output=Path(str(container_output_h5ad)),
                        tx_model_name=tx_model_name,
                        cell_type_key=cell_type_key,
                        gene_id_key=gene_id_key,
                        batch_size=batch_size,
                        seq_len_dataset=seq_len_dataset,
                        num_workers=num_workers,
                        prefetch_factor=prefetch_factor,
                    )

                    emb_matrix, used_bs, _used_key, err = run_tahoex1_segment_with_retry(
                        python_bin=python_bin,
                        predict_script=predict_script,
                        model_root=model_root,
                        model_package_dir=model_package_dir,
                        docker_image=docker_image,
                        docker_shm_size=docker_shm_size,
                        job_dir=tmp_dir,
                        cfg_yaml=cfg_yaml,
                        output_h5ad=output_h5ad,
                        initial_batch_size=batch_size,
                        tx_model_name=tx_model_name,
                        expected_ids=kept_ids,
                    )

                    if emb_matrix.size == 0:
                        raise RuntimeError(f"Tahoe-x1 segment failed after retry (batch_size_final={used_bs}): {err}")

                    if emb_dim is None:
                        emb_dim = int(emb_matrix.shape[1])
                    elif int(emb_matrix.shape[1]) != emb_dim:
                        raise RuntimeError(
                            f"Embedding dim mismatch: got={emb_matrix.shape[1]} expected={emb_dim}"
                        )

                    for cid, vec in zip(kept_ids, emb_matrix):
                        if not np.isfinite(vec).all():
                            invalid_cells.add(cid)
                            continue
                        vectors_by_cell[cid] = np.asarray(vec, dtype=np.float32)

                except Exception as exc:  # noqa: BLE001
                    if len(seg_local) > 1:
                        mid = len(seg_local) // 2
                        queue.insert(0, seg_local[mid:])
                        queue.insert(0, seg_local[:mid])
                        chunk_errors.append(
                            f"{side} chunk[{start}:{end}] segment split len={len(seg_local)} due to error: {exc}"
                        )
                    else:
                        invalid_cells.add(seg_ids[0])
                        chunk_errors.append(
                            f"{side} chunk[{start}:{end}] segment single-row failure cell={seg_ids[0]} error: {exc}"
                        )
                        chunk_success = False

            del seg_counts
            del seg_mapped
            del kept_counts
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

        if chunk_success:
            stats["n_chunks_succeeded"] += 1
        else:
            stats["n_chunks_failed"] += 1

    del counts_t

    if emb_dim is None:
        raise RuntimeError(f"{side} embedding failed for all chunks; cannot determine embedding dimension.")

    return vectors_by_cell, sorted(invalid_cells), emb_dim, stats, chunk_errors


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
    snapshot_outputs: List[Path] = []
    stage_outputs: List[Path] = []

    init_global_seed(GLOBAL_SEED)
    assertions.append(
        {
            "name": "global_seed_locked",
            "pass": True,
            "details": {"rules": ["GLOBAL_SEED must equal 619"], "seed": GLOBAL_SEED},
            "counterexamples": [],
        }
    )

    fm_cfg, fm_cfg_source = pick_tahoex1_cfg(config)

    model_dir: Optional[Path]
    if args.model_dir is not None:
        model_dir = args.model_dir.resolve()
        model_dir_source = "cli(--model-dir)"
    elif "model_dir" in fm_cfg and fm_cfg["model_dir"] is not None:
        model_dir = resolve_config_path(project_root, str(fm_cfg["model_dir"]))
        model_dir_source = f"{fm_cfg_source}.model_dir"
    else:
        model_dir = None
        model_dir_source = "missing"

    if model_dir is None or not model_dir.is_dir():
        assertions.append(
            {
                "name": "tahoex1_model_dir_resolved",
                "pass": False,
                "details": {
                    "rules": [
                        "Tahoe-x1 model dir must be provided by --model-dir or config.fm_extractors.{tahoex1|tahoe_x1|tahoe-x1}.model_dir"
                    ],
                    "model_dir_source": model_dir_source,
                    "model_dir": str(model_dir) if model_dir else "NA",
                },
                "counterexamples": [{"model_dir": str(model_dir) if model_dir else "NA"}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Missing valid Tahoe-x1 model directory.", file=sys.stderr)
        return 5

    model_size = str(args.model_size) if args.model_size is not None else str(fm_cfg.get("model_size", "70m"))
    if model_size not in {"70m", "1b", "3b"}:
        print(f"[ERROR] Unsupported --model-size={model_size}", file=sys.stderr)
        return 6

    python_bin = str(args.python_bin) if args.python_bin is not None else str(fm_cfg.get("python_bin", "python"))
    batch_size = int(args.batch_size) if args.batch_size is not None else int(fm_cfg.get("batch_size", 8))
    cell_chunk_size = (
        int(args.cell_chunk_size) if args.cell_chunk_size is not None else int(fm_cfg.get("cell_chunk_size", 2048))
    )
    seq_len_dataset = (
        int(args.seq_len_dataset) if args.seq_len_dataset is not None else int(fm_cfg.get("seq_len_dataset", 2048))
    )
    num_workers = int(args.num_workers) if args.num_workers is not None else int(fm_cfg.get("num_workers", 8))
    prefetch_factor = (
        int(args.prefetch_factor) if args.prefetch_factor is not None else int(fm_cfg.get("prefetch_factor", 48))
    )
    tx_model_name = (
        str(args.tx_model_name)
        if args.tx_model_name is not None
        else str(fm_cfg.get("tx_model_name", f"Tx1-{model_size}"))
    )
    cell_type_key = str(args.cell_type_key) if args.cell_type_key is not None else str(fm_cfg.get("cell_type_key", "cell_type"))
    gene_id_key = str(args.gene_id_key) if args.gene_id_key is not None else str(fm_cfg.get("gene_id_key", "ensembl_id"))
    default_cell_type = (
        str(args.default_cell_type)
        if args.default_cell_type is not None
        else str(fm_cfg.get("default_cell_type", "unknown"))
    )
    docker_image = str(args.docker_image) if args.docker_image is not None else str(fm_cfg.get("docker_image", DEFAULT_DOCKER_IMAGE))
    docker_shm_size = (
        str(args.docker_shm_size) if args.docker_shm_size is not None else str(fm_cfg.get("docker_shm_size", DEFAULT_DOCKER_SHM_SIZE))
    )

    try:
        model_root, predict_script, model_package_dir, gene_map_pkl, asset_sources = resolve_tahoex1_assets(
            project_root=project_root,
            model_dir=model_dir,
            cfg=fm_cfg,
            args=args,
            model_size=model_size,
        )
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "tahoex1_assets_resolved",
                "pass": False,
                "details": {"error": str(exc), "model_dir": str(model_dir), "model_size": model_size},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 7

    assertions.append(
        {
            "name": "tahoex1_assets_resolved",
            "pass": True,
            "details": {
                "model_dir": str(model_dir),
                "model_dir_source": model_dir_source,
                "fm_cfg_source": fm_cfg_source,
                "model_root": str(model_root),
                "predict_script": str(predict_script),
                "model_package_dir": str(model_package_dir),
                "gene_map_pkl": str(gene_map_pkl),
                "asset_sources": asset_sources,
                "python_bin": python_bin,
                "docker_image": docker_image,
                "docker_shm_size": docker_shm_size,
                "model_size": model_size,
                "tx_model_name": tx_model_name,
                "batch_size": int(batch_size),
                "cell_chunk_size": int(cell_chunk_size),
                "seq_len_dataset": int(seq_len_dataset),
                "num_workers": int(num_workers),
                "prefetch_factor": int(prefetch_factor),
                "cell_type_key": cell_type_key,
                "gene_id_key": gene_id_key,
                "default_cell_type": default_cell_type,
            },
            "counterexamples": [],
        }
    )

    k562_dir = task2_snapshot / "k562"
    derived_dir = k562_dir / "derived"
    pair_path = (derived_dir / "pair_list.parquet").resolve()
    delta_meta_path = (derived_dir / "delta_meta.csv").resolve()
    shared_var_path = (k562_dir / "shared_var_names.csv").resolve()
    crispr_meta_path = (k562_dir / "CRISPR_meta.csv").resolve()
    drug_meta_path = (k562_dir / "Drug_meta.csv").resolve()
    crispr_counts_path = (k562_dir / "CRISPR_counts.pt").resolve()
    drug_counts_path = (k562_dir / "Drug_counts.pt").resolve()

    required_inputs = [
        pair_path,
        delta_meta_path,
        shared_var_path,
        crispr_meta_path,
        drug_meta_path,
        crispr_counts_path,
        drug_counts_path,
        predict_script,
        model_package_dir,
        gene_map_pkl,
    ]
    missing_inputs = [path for path in required_inputs if not path.exists()]
    if missing_inputs:
        assertions.append(
            {
                "name": "task2_snapshot_inputs_present",
                "pass": False,
                "details": {"rules": ["All required inputs must exist."]},
                "counterexamples": [{"missing_input": str(p)} for p in missing_inputs[:MAX_COUNTEREXAMPLES]],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Missing required inputs: {missing_inputs}", file=sys.stderr)
        return 8

    input_paths.extend(
        [
            pair_path,
            delta_meta_path,
            shared_var_path,
            crispr_meta_path,
            drug_meta_path,
            crispr_counts_path,
            drug_counts_path,
            model_dir,
            model_root,
            predict_script,
            model_package_dir,
            gene_map_pkl,
        ]
    )

    assertions.append(
        {
            "name": "task2_snapshot_inputs_present",
            "pass": True,
            "details": {
                "rules": ["All required task2 snapshot inputs must exist."],
                "k562_dir": str(k562_dir.resolve()),
            },
            "counterexamples": [],
        }
    )

    try:
        pair_df = pd.read_parquet(
            pair_path,
            columns=[
                "treated_cell_id",
                "control_cell_id",
                "control_rank",
                "n_controls_used",
                "dataset_side",
            ],
        )
        delta_meta = pd.read_csv(delta_meta_path)
        shared_var = pd.read_csv(shared_var_path)
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "snapshot_loading",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Failed loading snapshot inputs: {exc}", file=sys.stderr)
        return 9

    ensure_required_columns(
        pair_df,
        ["treated_cell_id", "control_cell_id", "control_rank", "n_controls_used", "dataset_side"],
        "pair_list.parquet",
    )
    ensure_required_columns(delta_meta, ["row_id", "treated_cell_id", "n_controls_used"], "delta_meta.csv")

    if "dataset_side" not in delta_meta.columns:
        if "perturbation_class" not in delta_meta.columns:
            raise ValueError("delta_meta.csv requires dataset_side or perturbation_class")
        delta_meta["dataset_side"] = delta_meta["perturbation_class"].map(
            lambda x: "CRISPR" if normalize_scalar(x).lower() == "genetic" else "DRUG"
        )

    if "gene_symbol" not in shared_var.columns:
        raise ValueError("shared_var_names.csv must include gene_symbol column")
    shared_genes = shared_var["gene_symbol"].fillna("").astype(str).tolist()

    pair_df["dataset_side"] = pair_df["dataset_side"].astype(str).str.strip().str.upper()
    pair_df["treated_cell_id"] = pair_df["treated_cell_id"].map(normalize_scalar)
    pair_df["control_cell_id"] = pair_df["control_cell_id"].map(normalize_scalar)
    pair_df["control_rank"] = pd.to_numeric(pair_df["control_rank"], errors="raise").astype(np.int64)
    pair_df["n_controls_used"] = pd.to_numeric(pair_df["n_controls_used"], errors="raise").astype(np.int64)

    delta_meta["row_id"] = pd.to_numeric(delta_meta["row_id"], errors="raise").astype(np.int64)
    delta_meta["treated_cell_id"] = delta_meta["treated_cell_id"].map(normalize_scalar)
    delta_meta["n_controls_used"] = pd.to_numeric(delta_meta["n_controls_used"], errors="raise").astype(np.int64)
    delta_meta["dataset_side"] = delta_meta["dataset_side"].astype(str).str.strip().str.upper()
    delta_meta = delta_meta.sort_values("row_id", kind="mergesort").reset_index(drop=True)

    n_rows = int(len(delta_meta))
    expected_row_ids = np.arange(n_rows, dtype=np.int64)
    row_ids = delta_meta["row_id"].to_numpy(dtype=np.int64)
    row_alignment_input_ok = np.array_equal(row_ids, expected_row_ids)
    assertions.append(
        {
            "name": "delta_meta_row_id_contiguous",
            "pass": bool(row_alignment_input_ok),
            "details": {
                "rules": ["delta_meta.row_id must be exactly 0..N-1 for strict row alignment contract"],
                "n_rows": n_rows,
            },
            "counterexamples": []
            if row_alignment_input_ok
            else delta_meta.head(MAX_COUNTEREXAMPLES).to_dict(orient="records"),
        }
    )
    if not row_alignment_input_ok:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] delta_meta.row_id is not contiguous 0..N-1.", file=sys.stderr)
        return 10

    pair_sorted = pair_df.sort_values(["dataset_side", "treated_cell_id", "control_rank"], kind="mergesort")
    controls_by_key: Dict[Tuple[str, str], List[str]] = {}
    pair_contract_violations: List[Dict[str, object]] = []

    for (side, treated_cell_id), grp in pair_sorted.groupby(["dataset_side", "treated_cell_id"], sort=False):
        controls = grp["control_cell_id"].astype(str).tolist()
        n_used_unique = grp["n_controls_used"].unique()
        if len(n_used_unique) != 1:
            pair_contract_violations.append(
                {
                    "dataset_side": side,
                    "treated_cell_id": treated_cell_id,
                    "violation": "n_controls_used_inconsistent_within_group",
                }
            )
            if len(pair_contract_violations) >= MAX_COUNTEREXAMPLES:
                break
            continue
        n_used = int(n_used_unique[0])
        if len(controls) != n_used:
            pair_contract_violations.append(
                {
                    "dataset_side": side,
                    "treated_cell_id": treated_cell_id,
                    "violation": "group_size_ne_n_controls_used",
                    "n_controls_used": n_used,
                    "group_size": len(controls),
                }
            )
            if len(pair_contract_violations) >= MAX_COUNTEREXAMPLES:
                break
            continue
        controls_by_key[(str(side), str(treated_cell_id))] = controls

    pair_contract_pass = len(pair_contract_violations) == 0
    assertions.append(
        {
            "name": "pair_list_group_contract",
            "pass": pair_contract_pass,
            "details": {
                "rules": [
                    "Each (dataset_side, treated_cell_id) group must have consistent n_controls_used",
                    "Group row count must equal n_controls_used",
                ],
                "n_groups": int(len(controls_by_key)),
            },
            "counterexamples": pair_contract_violations[:MAX_COUNTEREXAMPLES],
        }
    )
    if not pair_contract_pass:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] pair_list grouping contract failed.", file=sys.stderr)
        return 11

    required_by_side: Dict[str, set[str]] = {"CRISPR": set(), "DRUG": set()}
    missing_pair_keys = 0
    for row in delta_meta.itertuples(index=False):
        side = str(row.dataset_side)
        treated_id = str(row.treated_cell_id)
        key = (side, treated_id)
        required_by_side.setdefault(side, set()).add(treated_id)
        controls = controls_by_key.get(key)
        if controls is None:
            missing_pair_keys += 1
            continue
        required_by_side[side].update(controls)

    assertions.append(
        {
            "name": "delta_rows_have_pair_keys",
            "pass": missing_pair_keys == 0,
            "details": {
                "rules": ["Every delta_meta row must resolve to a pair_list group key"],
                "missing_pair_keys": int(missing_pair_keys),
            },
            "counterexamples": []
            if missing_pair_keys == 0
            else [{"missing_pair_keys": int(missing_pair_keys)}],
        }
    )

    try:
        symbol2ens = load_symbol_to_ensembl(gene_map_pkl)
        src_cols, mapped_gene_symbols, mapped_ensembl = build_tahoex1_gene_mapper(shared_genes, symbol2ens)
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "tahoex1_gene_mapping_ready",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 12

    assertions.append(
        {
            "name": "tahoex1_gene_mapping_ready",
            "pass": True,
            "details": {
                "n_input_genes": int(len(shared_genes)),
                "n_mapped_genes": int(len(src_cols)),
                "min_mapped_genes_required": int(MIN_MAPPED_GENES),
                "gene_map_pkl": str(gene_map_pkl),
            },
            "counterexamples": [],
        }
    )

    side_vectors: Dict[str, Dict[str, np.ndarray]] = {}
    side_stats: Dict[str, Dict[str, int]] = {}
    side_chunk_errors: Dict[str, List[str]] = {}
    emb_dim: Optional[int] = None

    side_specs = [
        ("CRISPR", crispr_meta_path, crispr_counts_path),
        ("DRUG", drug_meta_path, drug_counts_path),
    ]

    for side, meta_path, counts_path in side_specs:
        try:
            vectors, _invalid_cells, side_dim, stats, chunk_errors = extract_side_embeddings_tahoex1(
                side=side,
                meta_path=meta_path,
                counts_path=counts_path,
                required_cell_ids=sorted(required_by_side.get(side, set())),
                src_cols=src_cols,
                mapped_gene_symbols=mapped_gene_symbols,
                mapped_ensembl=mapped_ensembl,
                python_bin=python_bin,
                predict_script=predict_script,
                model_root=model_root,
                model_package_dir=model_package_dir,
                docker_image=docker_image,
                docker_shm_size=docker_shm_size,
                tx_model_name=tx_model_name,
                batch_size=batch_size,
                seq_len_dataset=seq_len_dataset,
                num_workers=num_workers,
                prefetch_factor=prefetch_factor,
                cell_chunk_size=cell_chunk_size,
                cell_type_key=cell_type_key,
                gene_id_key=gene_id_key,
                default_cell_type=default_cell_type,
                stage_dir=stage_dir,
            )
        except Exception as exc:  # noqa: BLE001
            assertions.append(
                {
                    "name": f"tahoex1_embedding_side_{side.lower()}",
                    "pass": False,
                    "details": {"error": str(exc)},
                    "counterexamples": [{"error": str(exc)}],
                }
            )
            write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
            print(f"[ERROR] Failed embedding side={side}: {exc}", file=sys.stderr)
            return 13

        side_vectors[side] = vectors
        side_stats[side] = stats
        side_chunk_errors[side] = chunk_errors[:MAX_COUNTEREXAMPLES]

        if emb_dim is None:
            emb_dim = side_dim
        elif side_dim != emb_dim:
            assertions.append(
                {
                    "name": "embedding_dimension_consistent_across_sides",
                    "pass": False,
                    "details": {"expected_dim": emb_dim, "side": side, "side_dim": side_dim},
                    "counterexamples": [{"side": side, "side_dim": side_dim}],
                }
            )
            write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
            print("[ERROR] CRISPR/DRUG embedding dimensions do not match.", file=sys.stderr)
            return 14

    assert emb_dim is not None
    assertions.append(
        {
            "name": "embedding_dimension_consistent_across_sides",
            "pass": True,
            "details": {
                "rules": ["CRISPR and DRUG embeddings must share one latent dimension"],
                "embedding_dim": int(emb_dim),
                "side_stats": side_stats,
                "side_chunk_errors_topk": side_chunk_errors,
            },
            "counterexamples": [],
        }
    )

    fm_delta = np.full((n_rows, emb_dim), np.nan, dtype=np.float32)
    valid_mask = np.zeros(n_rows, dtype=bool)
    invalid_reason = np.full(n_rows, "", dtype=object)

    for row in delta_meta.itertuples(index=False):
        row_id = int(row.row_id)
        side = str(row.dataset_side)
        treated_id = str(row.treated_cell_id)
        key = (side, treated_id)

        controls = controls_by_key.get(key)
        if controls is None:
            invalid_reason[row_id] = "missing_pair_key"
            continue

        expected_n_controls = int(row.n_controls_used)
        if len(controls) != expected_n_controls:
            invalid_reason[row_id] = "n_controls_mismatch"
            continue
        if expected_n_controls < 1:
            invalid_reason[row_id] = "no_controls"
            continue

        treated_vec = side_vectors.get(side, {}).get(treated_id)
        if treated_vec is None or not np.isfinite(treated_vec).all():
            invalid_reason[row_id] = "treated_embedding_missing_or_invalid"
            continue

        ctrl_vecs: List[np.ndarray] = []
        missing_ctrl = False
        for control_id in controls:
            cvec = side_vectors.get(side, {}).get(control_id)
            if cvec is None or not np.isfinite(cvec).all():
                missing_ctrl = True
                break
            ctrl_vecs.append(cvec)
        if missing_ctrl:
            invalid_reason[row_id] = "control_embedding_missing_or_invalid"
            continue

        ctrl_mean = np.mean(np.stack(ctrl_vecs, axis=0), axis=0, dtype=np.float32)
        delta_vec = treated_vec - ctrl_mean
        if not np.isfinite(delta_vec).all():
            invalid_reason[row_id] = "delta_non_finite"
            continue

        fm_delta[row_id] = delta_vec.astype(np.float32, copy=False)
        valid_mask[row_id] = True

    valid_finite_ok = bool(np.isfinite(fm_delta[valid_mask]).all()) if bool(valid_mask.any()) else True
    invalid_nan_ok = bool(np.isnan(fm_delta[~valid_mask]).all()) if bool((~valid_mask).any()) else True

    assertions.append(
        {
            "name": "row_alignment_preserved_never_drop_rows",
            "pass": bool(
                fm_delta.shape[0] == n_rows
                and int(valid_mask.shape[0]) == n_rows
                and valid_finite_ok
                and invalid_nan_ok
            ),
            "details": {
                "rules": [
                    "fm_delta rows must equal len(delta_meta)",
                    "valid rows must be finite",
                    "invalid rows must remain all-NaN",
                ],
                "n_rows_delta_meta": n_rows,
                "n_rows_fm_delta": int(fm_delta.shape[0]),
                "n_valid": int(valid_mask.sum()),
                "n_invalid": int((~valid_mask).sum()),
            },
            "counterexamples": []
            if valid_finite_ok and invalid_nan_ok
            else [{"valid_finite_ok": valid_finite_ok, "invalid_nan_ok": invalid_nan_ok}],
        }
    )

    fm_root = (k562_dir / f"fm/{MODEL_NAME}").resolve()
    fm_root.mkdir(parents=True, exist_ok=True)
    fm_delta_path = (fm_root / "fm_delta.npy").resolve()
    fm_meta_path = (fm_root / "fm_delta_meta.csv").resolve()
    policy_path = (fm_root / "delta_operator_policy.json").resolve()

    np.save(fm_delta_path, fm_delta)

    fm_meta = pd.DataFrame(
        {
            "row_id": delta_meta["row_id"].to_numpy(dtype=np.int64),
            "treated_cell_id": delta_meta["treated_cell_id"].astype(str).to_numpy(),
            "valid_mask": valid_mask.astype(bool),
            "n_controls_used": delta_meta["n_controls_used"].to_numpy(dtype=np.int64),
            "dataset_side": delta_meta["dataset_side"].astype(str).to_numpy(),
            "invalid_reason": invalid_reason.astype(str),
            "model_name": np.full(n_rows, MODEL_NAME, dtype=object),
            "seed": np.full(n_rows, GLOBAL_SEED, dtype=np.int64),
        }
    )
    write_csv(fm_meta, fm_meta_path)

    policy_payload = {
        "model_name": MODEL_NAME,
        "operator_type": "euclidean",
        "eps": 1e-12,
        "delta_definition": "treated_minus_mean_controls",
        "row_alignment_contract": {
            "reference_table": "k562/derived/delta_meta.csv",
            "rule": "fm_delta[row_id] aligned to delta_meta.row_id with no row drops",
            "invalid_row_fill": "all_nan_and_valid_mask_false",
        },
        "seed": GLOBAL_SEED,
        "embedding_config": {
            "model_dir": str(model_dir),
            "model_root": str(model_root),
            "predict_script": str(predict_script),
            "model_package_dir": str(model_package_dir),
            "gene_map_pkl": str(gene_map_pkl),
            "python_bin": python_bin,
            "docker_image": docker_image,
            "docker_shm_size": docker_shm_size,
            "model_size": model_size,
            "tx_model_name": tx_model_name,
            "batch_size": int(batch_size),
            "cell_chunk_size": int(cell_chunk_size),
            "seq_len_dataset": int(seq_len_dataset),
            "num_workers": int(num_workers),
            "prefetch_factor": int(prefetch_factor),
            "cell_type_key": cell_type_key,
            "gene_id_key": gene_id_key,
            "default_cell_type": default_cell_type,
            "n_mapped_genes": int(len(src_cols)),
        },
    }
    write_json(policy_path, policy_payload)

    snapshot_outputs.extend([fm_delta_path, fm_meta_path, policy_path])

    output_routing_snapshot_pass = all(
        path.resolve().is_relative_to(task2_snapshot.resolve()) for path in snapshot_outputs
    )
    assertions.append(
        {
            "name": "snapshot_output_routing",
            "pass": bool(output_routing_snapshot_pass),
            "details": {
                "rules": ["FM outputs must be written under config.paths.task2_snapshot"],
                "task2_snapshot": str(task2_snapshot),
                "fm_root": str(fm_root),
            },
            "counterexamples": []
            if output_routing_snapshot_pass
            else [
                {"bad_output": str(path)}
                for path in snapshot_outputs
                if not path.resolve().is_relative_to(task2_snapshot.resolve())
            ][:MAX_COUNTEREXAMPLES],
        }
    )

    # IMPORTANT: use .absolute() (not .resolve()) so symlinked model assets keep logical path isolation.
    allowed_roots = [
        task2_snapshot.absolute(),
        model_dir.absolute(),
        model_root.absolute(),
        predict_script.parent.absolute(),
        model_package_dir.absolute(),
        gene_map_pkl.parent.absolute(),
    ]
    bad_inputs = [
        str(path.absolute())
        for path in sorted(set(input_paths), key=lambda p: str(p))
        if not any(path.absolute().is_relative_to(root) for root in allowed_roots)
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "rules": [
                    "Extractor reads only from task2 snapshot inputs and resolved model assets",
                    "Path isolation uses .absolute() to preserve symlinked model-root semantics",
                ],
                "task2_snapshot": str(task2_snapshot.absolute()),
                "model_dir": str(model_dir.absolute()),
                "n_inputs": int(len(set(input_paths))),
            },
            "counterexamples": [{"bad_input": p} for p in bad_inputs[:MAX_COUNTEREXAMPLES]],
        }
    )

    completed_at = utc_now_iso()
    run_manifest_path = (stage_dir / "run_manifest.json").resolve()
    audit_assertions_path = (stage_dir / "audit_assertions.json").resolve()
    manifest_path = (stage_dir / "manifest.json").resolve()
    stage_outputs.extend([run_manifest_path, audit_assertions_path, manifest_path])

    output_paths = [str(path.resolve()) for path in (snapshot_outputs + stage_outputs)]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/fm_extractors/extract_tahoex1.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
            "model_dir": str(model_dir),
            "model_dir_source": model_dir_source,
            "fm_cfg_source": fm_cfg_source,
            "model_root": str(model_root),
            "predict_script": str(predict_script),
            "model_package_dir": str(model_package_dir),
            "gene_map_pkl": str(gene_map_pkl),
            "python_bin": python_bin,
            "docker_image": docker_image,
            "docker_shm_size": docker_shm_size,
            "model_size": model_size,
            "tx_model_name": tx_model_name,
            "batch_size": int(batch_size),
            "cell_chunk_size": int(cell_chunk_size),
            "seq_len_dataset": int(seq_len_dataset),
            "num_workers": int(num_workers),
            "prefetch_factor": int(prefetch_factor),
            "cell_type_key": cell_type_key,
            "gene_id_key": gene_id_key,
            "default_cell_type": default_cell_type,
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths), key=lambda p: str(p))],
        "outputs": output_paths,
        "summary": {
            "n_rows": n_rows,
            "embedding_dim": int(emb_dim),
            "n_valid": int(valid_mask.sum()),
            "n_invalid": int((~valid_mask).sum()),
            "invalid_reason_top": pd.Series(invalid_reason[invalid_reason != ""]).value_counts().head(10).to_dict(),
            "side_stats": side_stats,
            "n_mapped_genes": int(len(src_cols)),
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

    output_routing_stage_pass = all(path.resolve().is_relative_to(stage_dir.resolve()) for path in stage_outputs)
    if not output_routing_stage_pass:
        print("[ERROR] Stage output routing assertion failed.", file=sys.stderr)
        return 15

    if not output_routing_snapshot_pass:
        print("[ERROR] Snapshot output routing assertion failed.", file=sys.stderr)
        return 16

    if bad_inputs:
        print("[ERROR] Input path isolation assertion failed.", file=sys.stderr)
        return 17

    if not (valid_finite_ok and invalid_nan_ok):
        print("[ERROR] Row alignment finite/NaN contract assertion failed.", file=sys.stderr)
        return 18

    return 0


if __name__ == "__main__":
    sys.exit(main())
