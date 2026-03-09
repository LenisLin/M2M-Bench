# SCRIPT_HEADER_CONTRACT
# Script: scripts/fm_extractors/extract_scfoundation.py
# Purpose: Build Task2 K562 scFoundation FM deltas using pair_list with strict row alignment to delta_meta.
# Inputs:
#   - Task2 Snapshot K562: config/config.yaml::paths.task2_snapshot
#     - k562/derived/pair_list.parquet
#     - k562/derived/delta_meta.csv
#     - k562/{CRISPR_counts.pt,CRISPR_meta.csv,Drug_counts.pt,Drug_meta.csv,shared_var_names.csv}
#   - scFoundation repo/model dir:
#     - config/config.yaml::fm_extractors.scfoundation.model_dir (optional)
#     - CLI override: --model-dir
# Outputs:
#   - fm_delta.npy: <paths.task2_snapshot>/k562/fm/scfoundation/
#   - fm_delta_meta.csv: <paths.task2_snapshot>/k562/fm/scfoundation/
#   - delta_operator_policy.json: <paths.task2_snapshot>/k562/fm/scfoundation/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
#     under runs/<run_id>/extract_scfoundation/
# Side Effects:
#   - Creates FM output directory under task2 snapshot
#   - Creates isolated run directory: runs/<run_id>/extract_scfoundation/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::fm_extractors.scfoundation (optional)
# Execution:
#   - python scripts/fm_extractors/extract_scfoundation.py --run-id <run_id> [--model-dir <repo_dir>]
# Failure Modes:
#   - Missing required snapshot files -> exit non-zero
#   - Missing/invalid scFoundation assets -> exit non-zero
#   - Complete embedding failure for a side -> exit non-zero
# Last Updated: 2026-03-05

"""
Task2 scFoundation FM extractor with strict row-preservation contract:

1) Output fm_delta.npy must have exactly N rows where N = len(delta_meta).
2) Row i always corresponds to delta_meta row_id i (contiguous 0..N-1 required).
3) Any per-row embedding failure does NOT drop rows:
   - fm_delta[row_id, :] = NaN
   - fm_delta_meta.valid_mask[row_id] = False
4) Output routing is fixed:
   - Snapshot outputs -> data/task2_snapshot_v1/k562/fm/scfoundation/
   - Run artifacts -> runs/<run_id>/extract_scfoundation/
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import random
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import yaml

STAGE = "extract_scfoundation"
MODEL_NAME = "scfoundation"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
GLOBAL_SEED = 619
MAX_COUNTEREXAMPLES = 5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Task2 K562 scFoundation FM delta extractor")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="Override scFoundation root configured in config/config.yaml",
    )
    parser.add_argument("--python-bin", type=str, default=None, help="Python executable for get_embedding.py")
    parser.add_argument("--version", type=str, default=None)
    parser.add_argument("--tgthighres", type=str, default=None)
    parser.add_argument("--pool-type", type=str, default=None)
    parser.add_argument("--pre-normalized", type=str, default=None)
    parser.add_argument("--cell-chunk-size", type=int, default=None)
    parser.add_argument("--min-positive-tokens", type=int, default=None)
    parser.add_argument("--to-upper", action="store_true")
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


def load_scfoundation_gene_list(path: Path, to_upper: bool) -> List[str]:
    try:
        df = pd.read_csv(path, sep="\t")
        if "gene_name" in df.columns:
            genes = df["gene_name"].astype(str).tolist()
        else:
            genes = df.iloc[:, 0].astype(str).tolist()
    except Exception:
        df = pd.read_csv(path, sep="\t", header=None)
        genes = df.iloc[:, 0].astype(str).tolist()
    if to_upper:
        genes = [g.upper() for g in genes]
    return genes


def build_gene_mapper(shared_genes: Sequence[str], sf_genes: Sequence[str], to_upper: bool) -> Tuple[np.ndarray, np.ndarray]:
    if to_upper:
        src_names = [str(g).upper() for g in shared_genes]
        tgt_names = [str(g).upper() for g in sf_genes]
    else:
        src_names = [str(g) for g in shared_genes]
        tgt_names = [str(g) for g in sf_genes]

    sf_index: Dict[str, int] = {g: i for i, g in enumerate(tgt_names)}
    src_cols: List[int] = []
    tgt_cols: List[int] = []
    for j, gene in enumerate(src_names):
        idx = sf_index.get(gene)
        if idx is None:
            continue
        src_cols.append(int(j))
        tgt_cols.append(int(idx))

    src_arr = np.asarray(src_cols, dtype=np.int64)
    tgt_arr = np.asarray(tgt_cols, dtype=np.int64)
    if src_arr.size == 0:
        raise RuntimeError("No overlap between shared_var_names genes and scFoundation gene list.")
    return src_arr, tgt_arr


def align_counts_to_scfoundation(
    counts_dense: np.ndarray,
    src_cols: np.ndarray,
    tgt_cols: np.ndarray,
    n_sf_genes: int,
) -> Any:
    try:
        from scipy import sparse  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("scipy is required for scFoundation alignment.") from exc

    x = counts_dense.astype(np.float32, copy=False)
    x_csr = sparse.csr_matrix(x)
    x_sub = x_csr[:, src_cols].tocoo()
    mapped_cols = tgt_cols[x_sub.col]
    x_full = sparse.coo_matrix(
        (x_sub.data, (x_sub.row, mapped_cols)),
        shape=(x.shape[0], int(n_sf_genes)),
        dtype=np.float32,
    ).tocsr()
    return x_full


def _safe_float(text: str, default: float = 0.0) -> float:
    try:
        return float(text)
    except Exception:
        return float(default)


def _row_nonzero_feature_count(x_csr: Any) -> np.ndarray:
    try:
        from scipy import sparse  # type: ignore

        if sparse.isspmatrix_csr(x_csr):
            return np.diff(x_csr.indptr).astype(np.int32, copy=False)
        if sparse.issparse(x_csr):
            x_csr = x_csr.tocsr()
            return np.diff(x_csr.indptr).astype(np.int32, copy=False)
    except Exception:
        pass
    return np.asarray((x_csr > 0).sum(axis=1)).reshape(-1).astype(np.int32, copy=False)


def estimate_positive_token_count(
    *,
    row_sum: np.ndarray,
    row_nnz: np.ndarray,
    tgthighres: str,
) -> np.ndarray:
    extra2 = (row_sum > 1.0).astype(np.int32, copy=False)
    extra1 = np.zeros_like(extra2, dtype=np.int32)

    mode = str(tgthighres).strip()[:1].lower()
    value = _safe_float(str(tgthighres).strip()[1:] if len(str(tgthighres).strip()) > 1 else "0", 0.0)
    if mode == "f":
        if value > 0:
            threshold = 1.0 / float(value)
            extra1 = (row_sum > threshold).astype(np.int32, copy=False)
    elif mode == "a":
        threshold = float(np.power(10.0, -float(value)))
        extra1 = (row_sum > threshold).astype(np.int32, copy=False)
    elif mode == "t":
        extra1 = np.full(row_sum.shape, 1 if value > 0 else 0, dtype=np.int32)

    return row_nnz.astype(np.int32, copy=False) + extra1 + extra2


def select_inference_rows(
    *,
    x_csr: Any,
    tgthighres: str,
    min_positive_tokens: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    row_sum = np.asarray(x_csr.sum(axis=1)).reshape(-1).astype(np.float32, copy=False)
    row_nnz = _row_nonzero_feature_count(x_csr)
    pos_tokens = estimate_positive_token_count(row_sum=row_sum, row_nnz=row_nnz, tgthighres=tgthighres)

    nonzero_mask = row_sum > 0
    token_safe_mask = pos_tokens >= int(min_positive_tokens)
    run_mask = nonzero_mask & token_safe_mask
    run_idx = np.where(run_mask)[0].astype(np.int64, copy=False)
    return run_idx, row_sum, row_nnz, pos_tokens


def write_h5ad_for_scfoundation(
    *,
    x_csr: Any,
    gene_names: Sequence[str],
    cell_ids: Sequence[str],
    out_h5ad: Path,
) -> None:
    try:
        import anndata as ad  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("anndata is required for scFoundation extractor.") from exc

    obs = pd.DataFrame(index=pd.Index([str(c) for c in cell_ids], name="cell_id"))
    obs["total_count"] = np.asarray(x_csr.sum(axis=1)).reshape(-1)
    var = pd.DataFrame(index=pd.Index([str(g) for g in gene_names], name="gene_symbol"))

    adata = ad.AnnData(X=x_csr, obs=obs, var=var)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)


def run_scfoundation_once(
    *,
    python_bin: str,
    get_embedding_py: Path,
    task_name: str,
    data_path: Path,
    save_path: Path,
    version: str,
    tgthighres: str,
    pool_type: str,
    pre_normalized: str,
) -> Tuple[bool, str]:
    save_path.mkdir(parents=True, exist_ok=True)

    cmd = [
        python_bin,
        str(get_embedding_py),
        "--task_name",
        task_name,
        "--input_type",
        "singlecell",
        "--output_type",
        "cell",
        "--pool_type",
        pool_type,
        "--tgthighres",
        tgthighres,
        "--data_path",
        str(data_path),
        "--save_path",
        str(save_path),
        "--pre_normalized",
        pre_normalized,
        "--version",
        version,
    ]

    result = subprocess.run(
        cmd,
        cwd=str(get_embedding_py.parent),
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        return True, ""

    msg = (result.stderr or "") + "\n" + (result.stdout or "")
    return False, msg[-2000:]


def pick_latest_embedding_file(save_path: Path, task_name: str) -> Path:
    patterns = [
        f"{task_name}*cell_embedding*.npy",
        f"*{task_name}*cell_embedding*.npy",
        f"{task_name}*.npy",
    ]
    candidates: List[Path] = []
    for p in patterns:
        candidates.extend(save_path.glob(p))
    unique = sorted(set(candidates), key=lambda x: x.stat().st_mtime, reverse=True)
    if not unique:
        raise FileNotFoundError(f"No scFoundation embedding npy found under {save_path} for task={task_name}")
    return unique[0]


def extract_segment_embedding(
    *,
    x_csr_segment: Any,
    segment_ids: Sequence[str],
    sf_genes: Sequence[str],
    side: str,
    chunk_tag: str,
    seg_tag: str,
    stage_dir: Path,
    python_bin: str,
    get_embedding_py: Path,
    version: str,
    tgthighres: str,
    pool_type: str,
    pre_normalized: str,
) -> Tuple[np.ndarray, Optional[str]]:
    with tempfile.TemporaryDirectory(prefix=f"scfoundation_{side}_{chunk_tag}_{seg_tag}_", dir=str(stage_dir)) as tmp:
        tmp_dir = Path(tmp)
        h5ad_path = tmp_dir / f"input_{side}_{chunk_tag}_{seg_tag}.h5ad"
        save_path = tmp_dir / "out"
        task_name = f"m2m_{MODEL_NAME}_{side}_{chunk_tag}_{seg_tag}"

        write_h5ad_for_scfoundation(
            x_csr=x_csr_segment,
            gene_names=sf_genes,
            cell_ids=segment_ids,
            out_h5ad=h5ad_path,
        )

        ok, err = run_scfoundation_once(
            python_bin=python_bin,
            get_embedding_py=get_embedding_py,
            task_name=task_name,
            data_path=h5ad_path,
            save_path=save_path,
            version=version,
            tgthighres=tgthighres,
            pool_type=pool_type,
            pre_normalized=pre_normalized,
        )
        if not ok:
            raise RuntimeError(f"get_embedding failed: {err}")

        emb_file = pick_latest_embedding_file(save_path, task_name)
        emb = np.load(emb_file)
        emb = np.asarray(emb, dtype=np.float32)
        if emb.ndim != 2:
            raise RuntimeError(f"scFoundation embedding must be 2D, got shape={emb.shape}")
        return emb, None


def extract_side_embeddings_scfoundation(
    *,
    side: str,
    meta_path: Path,
    counts_path: Path,
    required_cell_ids: Sequence[str],
    src_cols: np.ndarray,
    tgt_cols: np.ndarray,
    sf_genes: Sequence[str],
    python_bin: str,
    get_embedding_py: Path,
    version: str,
    tgthighres: str,
    pool_type: str,
    pre_normalized: str,
    cell_chunk_size: int,
    min_positive_tokens: int,
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
        "n_rows_nonzero": 0,
        "n_rows_run": 0,
        "n_rows_token_filtered": 0,
    }

    chunk_n = max(1, int(cell_chunk_size))
    for start in range(0, len(ordered_rows), chunk_n):
        end = min(start + chunk_n, len(ordered_rows))
        chunk_rows = ordered_rows[start:end]
        chunk_ids = ordered_ids[start:end]
        chunk_tag = f"{start}_{end}"
        stats["n_chunks_total"] += 1

        row_tensor = torch.as_tensor(chunk_rows, dtype=torch.long)
        chunk_dense = counts_t.index_select(0, row_tensor).detach().cpu().numpy().astype(np.float32, copy=False)

        try:
            x_aligned = align_counts_to_scfoundation(
                counts_dense=chunk_dense,
                src_cols=src_cols,
                tgt_cols=tgt_cols,
                n_sf_genes=len(sf_genes),
            )
        except Exception as exc:  # noqa: BLE001
            stats["n_chunks_failed"] += 1
            invalid_cells.update(chunk_ids)
            chunk_errors.append(f"{side} chunk[{start}:{end}] align_error: {exc}")
            continue

        run_idx, row_sum, _row_nnz, _pos_tokens = select_inference_rows(
            x_csr=x_aligned,
            tgthighres=tgthighres,
            min_positive_tokens=int(min_positive_tokens),
        )

        n_nonzero = int((row_sum > 0).sum())
        n_run = int(len(run_idx))
        stats["n_rows_nonzero"] += n_nonzero
        stats["n_rows_run"] += n_run
        stats["n_rows_token_filtered"] += int(max(0, n_nonzero - n_run))

        if n_run == 0:
            stats["n_chunks_failed"] += 1
            invalid_cells.update(chunk_ids)
            chunk_errors.append(f"{side} chunk[{start}:{end}] no rows passed token filter")
            continue

        run_idx_arr = np.asarray(run_idx, dtype=np.int64)
        queue: List[np.ndarray] = [run_idx_arr]
        segment_id = 0
        chunk_success = True

        while queue:
            seg = queue.pop(0)
            seg_ids = [chunk_ids[int(i)] for i in seg.tolist()]
            seg_x = x_aligned[seg]
            seg_tag = f"seg{segment_id:04d}"
            segment_id += 1

            try:
                emb_arr, _ = extract_segment_embedding(
                    x_csr_segment=seg_x,
                    segment_ids=seg_ids,
                    sf_genes=sf_genes,
                    side=side,
                    chunk_tag=chunk_tag,
                    seg_tag=seg_tag,
                    stage_dir=stage_dir,
                    python_bin=python_bin,
                    get_embedding_py=get_embedding_py,
                    version=version,
                    tgthighres=tgthighres,
                    pool_type=pool_type,
                    pre_normalized=pre_normalized,
                )
            except Exception as exc:  # noqa: BLE001
                if len(seg) > 1:
                    mid = len(seg) // 2
                    queue.insert(0, seg[mid:])
                    queue.insert(0, seg[:mid])
                    chunk_errors.append(
                        f"{side} chunk[{start}:{end}] segment split len={len(seg)} due to error: {exc}"
                    )
                    continue

                invalid_cells.add(seg_ids[0])
                chunk_errors.append(
                    f"{side} chunk[{start}:{end}] segment single-row failure cell={seg_ids[0]} error: {exc}"
                )
                chunk_success = False
                continue

            if emb_arr.shape[0] < 1:
                for cid in seg_ids:
                    invalid_cells.add(cid)
                chunk_errors.append(f"{side} chunk[{start}:{end}] segment empty embedding output")
                chunk_success = False
                continue

            if emb_dim is None:
                emb_dim = int(emb_arr.shape[1])
            elif int(emb_arr.shape[1]) != emb_dim:
                for cid in seg_ids:
                    invalid_cells.add(cid)
                chunk_errors.append(
                    f"{side} chunk[{start}:{end}] dim mismatch: got={emb_arr.shape[1]} expected={emb_dim}"
                )
                chunk_success = False
                continue

            n_copy = min(len(seg_ids), int(emb_arr.shape[0]))
            for i in range(n_copy):
                vec = emb_arr[i]
                cid = seg_ids[i]
                if not np.isfinite(vec).all():
                    invalid_cells.add(cid)
                    continue
                vectors_by_cell[cid] = np.asarray(vec, dtype=np.float32)

            if n_copy < len(seg_ids):
                for i in range(n_copy, len(seg_ids)):
                    invalid_cells.add(seg_ids[i])
                chunk_errors.append(
                    f"{side} chunk[{start}:{end}] segment row mismatch emb_rows={emb_arr.shape[0]} expected={len(seg_ids)}"
                )
                chunk_success = False

        if chunk_success:
            stats["n_chunks_succeeded"] += 1
        else:
            stats["n_chunks_failed"] += 1

        del chunk_dense
        del x_aligned
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    del counts_t

    if emb_dim is None:
        raise RuntimeError(f"{side} embedding failed for all chunks; cannot determine embedding dimension.")

    return vectors_by_cell, sorted(invalid_cells), emb_dim, stats, chunk_errors


def resolve_scfoundation_assets(
    *,
    project_root: Path,
    model_dir: Path,
) -> Tuple[Path, Path, Path, Dict[str, str]]:
    sources: Dict[str, str] = {}

    model_dir_abs = model_dir.resolve()
    candidates: List[Tuple[Path, str]] = [
        (model_dir_abs / "model/get_embedding.py", "model_dir/model/get_embedding.py"),
        (model_dir_abs / "get_embedding.py", "model_dir/get_embedding.py"),
    ]

    get_embedding_py: Optional[Path] = None
    for p, source in candidates:
        if p.is_file():
            get_embedding_py = p.resolve()
            sources["get_embedding_py"] = source
            break

    if get_embedding_py is None:
        raise RuntimeError("Cannot resolve scFoundation get_embedding.py from --model-dir")

    scfoundation_root = get_embedding_py.parent.parent if get_embedding_py.parent.name == "model" else get_embedding_py.parent

    gene_list_candidates: List[Tuple[Path, str]] = [
        (scfoundation_root / "OS_scRNA_gene_index.19264.tsv", "root/OS_scRNA_gene_index.19264.tsv"),
        (scfoundation_root / "model/OS_scRNA_gene_index.19264.tsv", "root/model/OS_scRNA_gene_index.19264.tsv"),
        (get_embedding_py.parent / "OS_scRNA_gene_index.19264.tsv", "get_embedding_dir/OS_scRNA_gene_index.19264.tsv"),
    ]

    gene_list_path: Optional[Path] = None
    for p, source in gene_list_candidates:
        if p.is_file():
            gene_list_path = p.resolve()
            sources["gene_list_path"] = source
            break

    if gene_list_path is None:
        raise RuntimeError("Cannot resolve OS_scRNA_gene_index.19264.tsv from --model-dir")

    return scfoundation_root.resolve(), get_embedding_py.resolve(), gene_list_path.resolve(), sources


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

    fm_cfg_root = config.get("fm_extractors", {})
    fm_cfg: Mapping[str, object]
    if isinstance(fm_cfg_root, Mapping):
        node = fm_cfg_root.get(MODEL_NAME, {})
        fm_cfg = node if isinstance(node, Mapping) else {}
    else:
        fm_cfg = {}

    model_dir: Optional[Path]
    if args.model_dir is not None:
        model_dir = args.model_dir.resolve()
        model_dir_source = "cli(--model-dir)"
    elif "model_dir" in fm_cfg and fm_cfg["model_dir"] is not None:
        model_dir = resolve_config_path(project_root, str(fm_cfg["model_dir"]))
        model_dir_source = "config.fm_extractors.scfoundation.model_dir"
    else:
        model_dir = None
        model_dir_source = "missing"

    if model_dir is None or not model_dir.is_dir():
        assertions.append(
            {
                "name": "scfoundation_model_dir_resolved",
                "pass": False,
                "details": {
                    "rules": [
                        "scFoundation model dir must be provided by --model-dir or config.fm_extractors.scfoundation.model_dir"
                    ],
                    "model_dir_source": model_dir_source,
                    "model_dir": str(model_dir) if model_dir else "NA",
                },
                "counterexamples": [{"model_dir": str(model_dir) if model_dir else "NA"}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Missing valid scFoundation model directory.", file=sys.stderr)
        return 5

    try:
        scfoundation_root, get_embedding_py, gene_list_path, asset_sources = resolve_scfoundation_assets(
            project_root=project_root,
            model_dir=model_dir,
        )
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "scfoundation_assets_resolved",
                "pass": False,
                "details": {"error": str(exc), "model_dir": str(model_dir)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 6

    python_bin = args.python_bin if args.python_bin else str(fm_cfg.get("python_bin", "python"))
    version = str(args.version) if args.version is not None else str(fm_cfg.get("version", "rde"))
    tgthighres = str(args.tgthighres) if args.tgthighres is not None else str(fm_cfg.get("tgthighres", "f1"))
    pool_type = str(args.pool_type) if args.pool_type is not None else str(fm_cfg.get("pool_type", "all"))
    pre_normalized = (
        str(args.pre_normalized) if args.pre_normalized is not None else str(fm_cfg.get("pre_normalized", "F"))
    )
    cell_chunk_size = (
        int(args.cell_chunk_size) if args.cell_chunk_size is not None else int(fm_cfg.get("cell_chunk_size", 8000))
    )
    min_positive_tokens = (
        int(args.min_positive_tokens)
        if args.min_positive_tokens is not None
        else int(fm_cfg.get("min_positive_tokens", 2))
    )
    to_upper = bool(args.to_upper or bool(fm_cfg.get("to_upper", False)))

    assertions.append(
        {
            "name": "scfoundation_assets_resolved",
            "pass": True,
            "details": {
                "model_dir": str(model_dir),
                "model_dir_source": model_dir_source,
                "scfoundation_root": str(scfoundation_root),
                "get_embedding_py": str(get_embedding_py),
                "gene_list_path": str(gene_list_path),
                "asset_sources": asset_sources,
                "python_bin": python_bin,
                "version": version,
                "tgthighres": tgthighres,
                "pool_type": pool_type,
                "pre_normalized": pre_normalized,
                "cell_chunk_size": int(cell_chunk_size),
                "min_positive_tokens": int(min_positive_tokens),
                "to_upper": bool(to_upper),
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
        get_embedding_py,
        gene_list_path,
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
        return 7

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
            scfoundation_root,
            get_embedding_py,
            gene_list_path,
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
        return 8

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
        return 9

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
        return 10

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
        sf_genes = load_scfoundation_gene_list(gene_list_path, to_upper=to_upper)
        src_cols, tgt_cols = build_gene_mapper(shared_genes=shared_genes, sf_genes=sf_genes, to_upper=to_upper)
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "scfoundation_gene_mapping_ready",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 11

    assertions.append(
        {
            "name": "scfoundation_gene_mapping_ready",
            "pass": True,
            "details": {
                "n_input_genes": int(len(shared_genes)),
                "n_scfoundation_genes": int(len(sf_genes)),
                "n_overlap_genes": int(len(src_cols)),
                "to_upper": bool(to_upper),
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
            vectors, _invalid_cells, side_dim, stats, chunk_errors = extract_side_embeddings_scfoundation(
                side=side,
                meta_path=meta_path,
                counts_path=counts_path,
                required_cell_ids=sorted(required_by_side.get(side, set())),
                src_cols=src_cols,
                tgt_cols=tgt_cols,
                sf_genes=sf_genes,
                python_bin=python_bin,
                get_embedding_py=get_embedding_py,
                version=version,
                tgthighres=tgthighres,
                pool_type=pool_type,
                pre_normalized=pre_normalized,
                cell_chunk_size=cell_chunk_size,
                min_positive_tokens=min_positive_tokens,
                stage_dir=stage_dir,
            )
        except Exception as exc:  # noqa: BLE001
            assertions.append(
                {
                    "name": f"scfoundation_embedding_side_{side.lower()}",
                    "pass": False,
                    "details": {"error": str(exc)},
                    "counterexamples": [{"error": str(exc)}],
                }
            )
            write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
            print(f"[ERROR] Failed embedding side={side}: {exc}", file=sys.stderr)
            return 12

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
            return 13

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
            "scfoundation_root": str(scfoundation_root),
            "get_embedding_py": str(get_embedding_py),
            "gene_list_path": str(gene_list_path),
            "python_bin": python_bin,
            "version": version,
            "tgthighres": tgthighres,
            "pool_type": pool_type,
            "pre_normalized": pre_normalized,
            "cell_chunk_size": int(cell_chunk_size),
            "min_positive_tokens": int(min_positive_tokens),
            "to_upper": bool(to_upper),
            "n_overlap_genes": int(len(src_cols)),
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
        scfoundation_root.absolute(),
        get_embedding_py.parent.absolute(),
        gene_list_path.parent.absolute(),
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
        "script_path": str((project_root / "scripts/fm_extractors/extract_scfoundation.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
            "model_dir": str(model_dir),
            "model_dir_source": model_dir_source,
            "scfoundation_root": str(scfoundation_root),
            "get_embedding_py": str(get_embedding_py),
            "gene_list_path": str(gene_list_path),
            "python_bin": python_bin,
            "version": version,
            "tgthighres": tgthighres,
            "pool_type": pool_type,
            "pre_normalized": pre_normalized,
            "cell_chunk_size": int(cell_chunk_size),
            "min_positive_tokens": int(min_positive_tokens),
            "to_upper": bool(to_upper),
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
            "n_overlap_genes": int(len(src_cols)),
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
        return 14

    if not output_routing_snapshot_pass:
        print("[ERROR] Snapshot output routing assertion failed.", file=sys.stderr)
        return 15

    if bad_inputs:
        print("[ERROR] Input path isolation assertion failed.", file=sys.stderr)
        return 16

    if not (valid_finite_ok and invalid_nan_ok):
        print("[ERROR] Row alignment finite/NaN contract assertion failed.", file=sys.stderr)
        return 17

    return 0


if __name__ == "__main__":
    sys.exit(main())
