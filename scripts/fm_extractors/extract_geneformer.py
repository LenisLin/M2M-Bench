# SCRIPT_HEADER_CONTRACT
# Script: scripts/fm_extractors/extract_geneformer.py
# Purpose: Build Task2 K562 Geneformer FM deltas using pair_list with strict row alignment to delta_meta.
# Inputs:
#   - Task2 Snapshot K562: config/config.yaml::paths.task2_snapshot
#     - k562/derived/pair_list.parquet
#     - k562/derived/delta_meta.csv
#     - k562/{CRISPR_counts.pt,CRISPR_meta.csv,Drug_counts.pt,Drug_meta.csv,shared_var_names.csv}
#   - Geneformer model/code dir:
#     - config/config.yaml::fm_extractors.geneformer.model_dir (optional)
#     - CLI override: --model-dir
# Outputs:
#   - fm_delta.npy: <paths.task2_snapshot>/k562/fm/geneformer/
#   - fm_delta_meta.csv: <paths.task2_snapshot>/k562/fm/geneformer/
#   - delta_operator_policy.json: <paths.task2_snapshot>/k562/fm/geneformer/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
#     under runs/<run_id>/extract_geneformer/
# Side Effects:
#   - Creates FM output directory under task2 snapshot
#   - Creates isolated run directory: runs/<run_id>/extract_geneformer/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::fm_extractors.geneformer (optional)
# Execution:
#   - python scripts/fm_extractors/extract_geneformer.py --run-id <run_id> [--model-dir <dir>]
# Failure Modes:
#   - Missing required snapshot files -> exit non-zero
#   - Missing/invalid Geneformer model assets -> exit non-zero
#   - Complete embedding failure for a side -> exit non-zero
# Last Updated: 2026-03-05

"""
Task2 Geneformer FM extractor with strict row-preservation contract:

1) Output fm_delta.npy must have exactly N rows where N = len(delta_meta).
2) Row i always corresponds to delta_meta row_id i (contiguous 0..N-1 required).
3) Any per-row embedding failure does NOT drop rows:
   - fm_delta[row_id, :] = NaN
   - fm_delta_meta.valid_mask[row_id] = False
4) Output routing is fixed:
   - Snapshot outputs -> data/task2_snapshot_v1/k562/fm/geneformer/
   - Run artifacts -> runs/<run_id>/extract_geneformer/
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pickle
import random
import shutil
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

STAGE = "extract_geneformer"
MODEL_NAME = "geneformer"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
GLOBAL_SEED = 619
MAX_COUNTEREXAMPLES = 5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Task2 K562 Geneformer FM delta extractor")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="Override Geneformer code/model root configured in config/config.yaml",
    )
    parser.add_argument(
        "--gene-dict-pkl",
        type=Path,
        default=None,
        help="Optional override for Geneformer gene_name_id_dict_gc104M.pkl",
    )
    parser.add_argument("--nproc", type=int, default=None)
    parser.add_argument("--cell-chunk-size", type=int, default=None)
    parser.add_argument("--tokenizer-chunk-size", type=int, default=None)
    parser.add_argument("--forward-batch-size", type=int, default=None)
    parser.add_argument("--emb-mode", type=str, default=None, choices=["cls", "cell"])
    parser.add_argument("--emb-layer", type=int, default=None)
    parser.add_argument(
        "--use-h5ad-index",
        action="store_true",
        help="Use AnnData var index as input key for Geneformer mapping (legacy-compatible).",
    )
    parser.add_argument(
        "--force-retokenize",
        action="store_true",
        help="If set, regenerate tokenized .dataset files for each chunk.",
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


def is_oom_error(exc: BaseException) -> bool:
    text = str(exc).lower()
    return "out of memory" in text or ("cuda" in text and "memory" in text)


def choose_side_id_column(meta: pd.DataFrame, side: str) -> str:
    if side == "CRISPR":
        candidates = ["cell_barcode", "cell_id", "Unnamed: 0"]
    else:
        candidates = ["Unnamed: 0", "cell_id", "cell_barcode"]
    for column in candidates:
        if column in meta.columns:
            return column
    raise ValueError(f"Unable to locate cell-id column for side={side}")


def looks_like_hf_model_dir(path: Path) -> bool:
    return path.is_dir() and (path / "config.json").is_file() and (
        (path / "model.safetensors").is_file() or (path / "pytorch_model.bin").is_file()
    )


def first_existing_path(paths: Sequence[Path], expect: str) -> Optional[Path]:
    for path in paths:
        if expect == "file" and path.is_file():
            return path.resolve()
        if expect == "dir" and path.is_dir():
            return path.resolve()
    return None


def choose_geneformer_classes(code_root: Optional[Path]) -> Tuple[Any, Any, str]:
    try:
        from geneformer import EmbExtractor, TranscriptomeTokenizer  # type: ignore

        return TranscriptomeTokenizer, EmbExtractor, "python_env"
    except Exception as first_exc:  # noqa: BLE001
        if code_root is not None and (code_root / "geneformer").is_dir():
            sys.path.insert(0, str(code_root))
            try:
                from geneformer import EmbExtractor, TranscriptomeTokenizer  # type: ignore

                return TranscriptomeTokenizer, EmbExtractor, f"sys.path:{code_root}"
            except Exception as second_exc:  # noqa: BLE001
                raise RuntimeError(
                    "Unable to import geneformer. Activate the Geneformer env or provide --model-dir with usable source code."
                ) from second_exc
        raise RuntimeError(
            "Unable to import geneformer. Activate the Geneformer env or provide --model-dir with usable source code."
        ) from first_exc


def resolve_geneformer_assets(
    *,
    project_root: Path,
    model_dir: Path,
    cfg: Mapping[str, object],
    gene_dict_override: Optional[Path],
) -> Tuple[Optional[Path], Path, Path, Dict[str, str]]:
    sources: Dict[str, str] = {}

    code_root: Optional[Path] = None
    if (model_dir / "geneformer").is_dir():
        code_root = model_dir.resolve()
        sources["code_root"] = "model_dir"
    elif (model_dir.parent / "geneformer").is_dir():
        code_root = model_dir.parent.resolve()
        sources["code_root"] = "model_dir_parent"

    pretrained_candidates: List[Tuple[Path, str]] = []
    if "pretrained_model_dir" in cfg and cfg["pretrained_model_dir"] is not None:
        p = resolve_config_path(project_root, str(cfg["pretrained_model_dir"]))
        pretrained_candidates.append((p, "config.fm_extractors.geneformer.pretrained_model_dir"))
    pretrained_candidates.extend(
        [
            (model_dir.resolve(), "model_dir"),
            ((model_dir / "Geneformer-V2-104M").resolve(), "model_dir/Geneformer-V2-104M"),
            ((model_dir / "Geneformer-V2-316M").resolve(), "model_dir/Geneformer-V2-316M"),
            ((model_dir / "Geneformer-V1-10M").resolve(), "model_dir/Geneformer-V1-10M"),
            ((model_dir.parent / "Geneformer-V2-104M").resolve(), "model_dir_parent/Geneformer-V2-104M"),
        ]
    )

    pretrained_model_dir: Optional[Path] = None
    for candidate, source in pretrained_candidates:
        if looks_like_hf_model_dir(candidate):
            pretrained_model_dir = candidate.resolve()
            sources["pretrained_model_dir"] = source
            break
    if pretrained_model_dir is None:
        raise RuntimeError(
            "Cannot resolve Geneformer pretrained model directory. Expected config.json + model weights."
        )

    gene_dict_candidates: List[Tuple[Path, str]] = []
    if gene_dict_override is not None:
        gene_dict_candidates.append((gene_dict_override.resolve(), "cli(--gene-dict-pkl)"))
    if "gene_mapping_file" in cfg and cfg["gene_mapping_file"] is not None:
        p = resolve_config_path(project_root, str(cfg["gene_mapping_file"]))
        gene_dict_candidates.append((p, "config.fm_extractors.geneformer.gene_mapping_file"))

    roots_for_gene_dict = [
        model_dir.resolve(),
        pretrained_model_dir.resolve(),
        pretrained_model_dir.parent.resolve(),
    ]
    if code_root is not None:
        roots_for_gene_dict.append(code_root.resolve())

    for root in roots_for_gene_dict:
        gene_dict_candidates.extend(
            [
                (root / "geneformer/gene_name_id_dict_gc104M.pkl", f"{root}/geneformer"),
                (root / "build/lib/geneformer/gene_name_id_dict_gc104M.pkl", f"{root}/build/lib/geneformer"),
                (root / "gene_name_id_dict_gc104M.pkl", f"{root}"),
            ]
        )

    gene_mapping_file: Optional[Path] = None
    for candidate, source in gene_dict_candidates:
        if candidate.is_file():
            gene_mapping_file = candidate.resolve()
            sources["gene_mapping_file"] = source
            break

    if gene_mapping_file is None:
        raise RuntimeError(
            "Cannot resolve gene_name_id_dict_gc104M.pkl. Provide --gene-dict-pkl or config.fm_extractors.geneformer.gene_mapping_file."
        )

    return code_root, pretrained_model_dir, gene_mapping_file, sources


def load_gene_symbol_to_ensembl(path: Path) -> Dict[str, str]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict in {path}, got {type(payload)}")
    result: Dict[str, str] = {}
    for key, value in payload.items():
        if key is None or value is None:
            continue
        k = str(key).strip()
        v = str(value).strip()
        if not k or not v:
            continue
        result[k] = v
    return result


def build_gene_mapping(shared_gene_symbols: Sequence[str], symbol2ens: Mapping[str, str]) -> Tuple[np.ndarray, List[str], List[str]]:
    keep_cols: List[int] = []
    kept_symbols: List[str] = []
    kept_ens: List[str] = []

    for j, sym in enumerate(shared_gene_symbols):
        sym_clean = str(sym).strip()
        if not sym_clean:
            continue
        ens = symbol2ens.get(sym_clean)
        if ens is None:
            ens = symbol2ens.get(sym_clean.upper())
        if ens is None:
            continue
        ens_clean = str(ens).strip()
        if not ens_clean.startswith("ENSG"):
            continue
        keep_cols.append(int(j))
        kept_symbols.append(sym_clean)
        kept_ens.append(ens_clean)

    keep_cols_arr = np.asarray(keep_cols, dtype=np.int64)
    if keep_cols_arr.size == 0:
        raise RuntimeError("No genes can be mapped to Ensembl IDs using the provided dictionary.")
    return keep_cols_arr, kept_symbols, kept_ens


def make_h5ad_from_chunk(
    *,
    counts_t: torch.Tensor,
    row_index: np.ndarray,
    chunk_ids: Sequence[str],
    meta_chunk: pd.DataFrame,
    keep_cols: np.ndarray,
    kept_symbols: Sequence[str],
    kept_ens: Sequence[str],
    out_h5ad: Path,
) -> None:
    try:
        import anndata as ad  # type: ignore
        from scipy import sparse  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "Unable to import anndata/scipy. Run this script in the Geneformer environment."
        ) from exc

    row_tensor = torch.as_tensor(row_index, dtype=torch.long)
    col_tensor = torch.as_tensor(keep_cols, dtype=torch.long)

    chunk_counts = counts_t.index_select(0, row_tensor).index_select(1, col_tensor)
    x = chunk_counts.detach().cpu().numpy().astype(np.float32, copy=False)
    x_sp = sparse.csr_matrix(x, dtype=np.float32)

    obs = meta_chunk.copy().reset_index(drop=True)
    obs["cell_id"] = [str(xid) for xid in chunk_ids]
    if "target" not in obs.columns:
        obs["target"] = "NA"
    if "benchmark_group" not in obs.columns:
        obs["benchmark_group"] = "NA"
    if "specificity_tier" not in obs.columns:
        obs["specificity_tier"] = "NA"
    if "clean_target_mapped" not in obs.columns:
        obs["clean_target_mapped"] = "NA"

    obs["n_counts"] = np.asarray(x_sp.sum(axis=1)).reshape(-1).astype(np.float32)
    obs["filter_pass"] = np.ones((x_sp.shape[0],), dtype=np.int64)

    var = pd.DataFrame(index=pd.Index([str(g) for g in kept_symbols], name="gene_symbol"))
    var["ensembl_id"] = [str(e) for e in kept_ens]

    adata = ad.AnnData(X=x_sp, obs=obs, var=var)
    adata.obs_names = pd.Index([str(xid) for xid in chunk_ids], name="cell_id")
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)


def tokenize_h5ad_to_dataset_singlefile(
    *,
    transcriptome_tokenizer_cls: Any,
    h5ad_path: Path,
    token_out_dir: Path,
    output_prefix: str,
    nproc: int,
    tokenizer_chunk_size: int,
    use_h5ad_index: bool,
    force_retokenize: bool,
) -> Path:
    token_out_dir.mkdir(parents=True, exist_ok=True)
    dataset_path = token_out_dir / f"{output_prefix}.dataset"
    if dataset_path.exists() and force_retokenize:
        shutil.rmtree(dataset_path, ignore_errors=True)

    if dataset_path.exists() and not force_retokenize:
        return dataset_path.resolve()

    tokenizer = transcriptome_tokenizer_cls(
        model_version="V2",
        nproc=int(nproc),
        chunk_size=int(tokenizer_chunk_size),
        use_h5ad_index=bool(use_h5ad_index),
        custom_attr_name_dict={
            "cell_id": "cell_id",
            "target": "target",
            "benchmark_group": "benchmark_group",
            "specificity_tier": "specificity_tier",
            "clean_target_mapped": "clean_target_mapped",
        },
    )

    tmp_dir = Path(tempfile.mkdtemp(prefix=f"geneformer_tokenize_{output_prefix}_"))
    try:
        tmp_h5ad = tmp_dir / h5ad_path.name
        shutil.copy2(h5ad_path, tmp_h5ad)
        tokenizer.tokenize_data(
            data_directory=str(tmp_dir),
            output_directory=str(token_out_dir),
            output_prefix=output_prefix,
            use_generator=True,
            file_format="h5ad",
        )
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    if dataset_path.exists():
        return dataset_path.resolve()

    candidates = sorted(
        [
            p
            for p in token_out_dir.iterdir()
            if p.name.startswith(output_prefix) and p.name.endswith(".dataset")
        ]
    )
    if not candidates:
        raise FileNotFoundError(f"Geneformer tokenization finished but no .dataset found in {token_out_dir}")
    return candidates[0].resolve()


def extract_embeddings_from_dataset(
    *,
    emb_extractor_cls: Any,
    pretrained_model_dir: Path,
    dataset_path: Path,
    output_dir: Path,
    output_prefix: str,
    forward_batch_size: int,
    nproc: int,
    emb_mode: str,
    emb_layer: int,
) -> Tuple[pd.DataFrame, np.ndarray]:
    output_dir.mkdir(parents=True, exist_ok=True)

    embex = emb_extractor_cls(
        model_type="Pretrained",
        num_classes=0,
        emb_mode=emb_mode,
        max_ncells=None,
        emb_layer=int(emb_layer),
        emb_label=["cell_id"],
        forward_batch_size=int(forward_batch_size),
        nproc=int(nproc),
        model_version="V2",
    )

    embs_df, embs_tensor = embex.extract_embs(
        model_directory=str(pretrained_model_dir),
        input_data_file=str(dataset_path),
        output_directory=str(output_dir),
        output_prefix=output_prefix,
        output_torch_embs=True,
    )

    arr = embs_tensor.detach().cpu().numpy().astype(np.float32, copy=False)
    if arr.ndim != 2:
        raise RuntimeError(f"Geneformer embedding tensor must be 2D, got shape={arr.shape}")
    return embs_df, arr


def align_embeddings_to_chunk_ids(
    *,
    embs_df: pd.DataFrame,
    emb_arr: np.ndarray,
    chunk_ids: Sequence[str],
) -> Tuple[Dict[str, np.ndarray], List[str], Optional[str]]:
    if emb_arr.shape[0] == 0:
        return {}, list(chunk_ids), "empty_embedding_rows"

    cell_id_col: Optional[str] = None
    if "cell_id" in embs_df.columns:
        cell_id_col = "cell_id"
    else:
        for candidate in ["cell", "barcode", "obs_name", "obs_names", "Unnamed: 0"]:
            if candidate in embs_df.columns:
                cell_id_col = candidate
                break

    vectors_by_cell: Dict[str, np.ndarray] = {}
    missing_cells: List[str] = []

    if cell_id_col is not None:
        emb_ids = embs_df[cell_id_col].astype(str).tolist()
        if len(emb_ids) != emb_arr.shape[0]:
            return {}, list(chunk_ids), f"id_length_mismatch:{len(emb_ids)}!={emb_arr.shape[0]}"

        first_pos: Dict[str, int] = {}
        for idx, cell_id in enumerate(emb_ids):
            if cell_id not in first_pos:
                first_pos[cell_id] = idx

        for cell_id in chunk_ids:
            pos = first_pos.get(str(cell_id))
            if pos is None:
                missing_cells.append(str(cell_id))
                continue
            vec = emb_arr[pos]
            if not np.isfinite(vec).all():
                missing_cells.append(str(cell_id))
                continue
            vectors_by_cell[str(cell_id)] = np.asarray(vec, dtype=np.float32)

        return vectors_by_cell, missing_cells, None

    # last-resort positional fallback (legacy behavior)
    n_copy = min(len(chunk_ids), emb_arr.shape[0])
    for i in range(n_copy):
        vec = emb_arr[i]
        if np.isfinite(vec).all():
            vectors_by_cell[str(chunk_ids[i])] = np.asarray(vec, dtype=np.float32)
        else:
            missing_cells.append(str(chunk_ids[i]))
    for i in range(n_copy, len(chunk_ids)):
        missing_cells.append(str(chunk_ids[i]))
    return vectors_by_cell, missing_cells, "position_fallback_without_cell_id"


def embed_chunk_with_retry(
    *,
    emb_extractor_cls: Any,
    pretrained_model_dir: Path,
    dataset_path: Path,
    output_dir: Path,
    output_prefix: str,
    initial_forward_batch_size: int,
    nproc: int,
    emb_mode: str,
    emb_layer: int,
) -> Tuple[Optional[pd.DataFrame], Optional[np.ndarray], int, Optional[str]]:
    batch_size = max(1, int(initial_forward_batch_size))
    last_error: Optional[str] = None

    while True:
        try:
            embs_df, arr = extract_embeddings_from_dataset(
                emb_extractor_cls=emb_extractor_cls,
                pretrained_model_dir=pretrained_model_dir,
                dataset_path=dataset_path,
                output_dir=output_dir,
                output_prefix=output_prefix,
                forward_batch_size=batch_size,
                nproc=nproc,
                emb_mode=emb_mode,
                emb_layer=emb_layer,
            )
            return embs_df, arr, batch_size, None
        except Exception as exc:  # noqa: BLE001
            last_error = str(exc)
            if batch_size <= 1 or not is_oom_error(exc):
                return None, None, batch_size, last_error
            batch_size = max(1, batch_size // 2)
            if torch.cuda.is_available():
                torch.cuda.empty_cache()


def extract_side_embeddings_geneformer(
    *,
    side: str,
    meta_path: Path,
    counts_path: Path,
    required_cell_ids: Sequence[str],
    keep_cols: np.ndarray,
    kept_symbols: Sequence[str],
    kept_ens: Sequence[str],
    transcriptome_tokenizer_cls: Any,
    emb_extractor_cls: Any,
    pretrained_model_dir: Path,
    nproc: int,
    cell_chunk_size: int,
    tokenizer_chunk_size: int,
    forward_batch_size: int,
    emb_mode: str,
    emb_layer: int,
    use_h5ad_index: bool,
    force_retokenize: bool,
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
    }

    chunk_n = max(1, int(cell_chunk_size))
    for start in range(0, len(ordered_rows), chunk_n):
        end = min(start + chunk_n, len(ordered_rows))
        chunk_rows = ordered_rows[start:end]
        chunk_ids = ordered_ids[start:end]
        stats["n_chunks_total"] += 1

        chunk_meta = meta.iloc[chunk_rows].copy().reset_index(drop=True)

        chunk_prefix = f"{side.lower()}_{start}_{end}"
        with tempfile.TemporaryDirectory(prefix=f"geneformer_{chunk_prefix}_", dir=str(stage_dir)) as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            h5ad_path = tmp_dir / f"{chunk_prefix}.h5ad"
            token_dir = tmp_dir / "tokenized"
            emb_out_dir = tmp_dir / "embeddings"

            try:
                make_h5ad_from_chunk(
                    counts_t=counts_t,
                    row_index=chunk_rows,
                    chunk_ids=chunk_ids,
                    meta_chunk=chunk_meta,
                    keep_cols=keep_cols,
                    kept_symbols=kept_symbols,
                    kept_ens=kept_ens,
                    out_h5ad=h5ad_path,
                )

                dataset_path = tokenize_h5ad_to_dataset_singlefile(
                    transcriptome_tokenizer_cls=transcriptome_tokenizer_cls,
                    h5ad_path=h5ad_path,
                    token_out_dir=token_dir,
                    output_prefix=chunk_prefix,
                    nproc=nproc,
                    tokenizer_chunk_size=tokenizer_chunk_size,
                    use_h5ad_index=use_h5ad_index,
                    force_retokenize=force_retokenize,
                )

                embs_df, emb_arr, used_bs, error_msg = embed_chunk_with_retry(
                    emb_extractor_cls=emb_extractor_cls,
                    pretrained_model_dir=pretrained_model_dir,
                    dataset_path=dataset_path,
                    output_dir=emb_out_dir,
                    output_prefix=chunk_prefix,
                    initial_forward_batch_size=forward_batch_size,
                    nproc=nproc,
                    emb_mode=emb_mode,
                    emb_layer=emb_layer,
                )

                if emb_arr is None or embs_df is None:
                    stats["n_chunks_failed"] += 1
                    invalid_cells.update(chunk_ids)
                    chunk_errors.append(
                        f"{side} chunk[{start}:{end}] failed after retry "
                        f"(forward_batch_size_final={used_bs}): {error_msg}"
                    )
                    continue

                if emb_dim is None:
                    emb_dim = int(emb_arr.shape[1])
                elif int(emb_arr.shape[1]) != emb_dim:
                    stats["n_chunks_failed"] += 1
                    invalid_cells.update(chunk_ids)
                    chunk_errors.append(
                        f"{side} chunk[{start}:{end}] dim mismatch: got={emb_arr.shape[1]}, expected={emb_dim}"
                    )
                    continue

                aligned, missing_chunk_cells, align_note = align_embeddings_to_chunk_ids(
                    embs_df=embs_df,
                    emb_arr=emb_arr,
                    chunk_ids=chunk_ids,
                )

                if align_note is not None:
                    chunk_errors.append(f"{side} chunk[{start}:{end}] alignment_note={align_note}")

                for cid in missing_chunk_cells:
                    invalid_cells.add(str(cid))

                for cid, vec in aligned.items():
                    if not np.isfinite(vec).all():
                        invalid_cells.add(str(cid))
                        continue
                    vectors_by_cell[str(cid)] = np.asarray(vec, dtype=np.float32)

                stats["n_chunks_succeeded"] += 1

            except Exception as exc:  # noqa: BLE001
                stats["n_chunks_failed"] += 1
                invalid_cells.update(chunk_ids)
                chunk_errors.append(f"{side} chunk[{start}:{end}] exception: {exc}")

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    del counts_t

    if emb_dim is None:
        raise RuntimeError(
            f"{side} embedding failed for all chunks; cannot determine embedding dimension."
        )

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
        model_dir_source = "config.fm_extractors.geneformer.model_dir"
    else:
        model_dir = None
        model_dir_source = "missing"

    if model_dir is None or not model_dir.is_dir():
        assertions.append(
            {
                "name": "geneformer_model_dir_resolved",
                "pass": False,
                "details": {
                    "rules": [
                        "Geneformer model/code dir must be provided by --model-dir or config.fm_extractors.geneformer.model_dir"
                    ],
                    "model_dir_source": model_dir_source,
                    "model_dir": str(model_dir) if model_dir else "NA",
                },
                "counterexamples": [{"model_dir": str(model_dir) if model_dir else "NA"}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Missing valid Geneformer model directory.", file=sys.stderr)
        return 5

    try:
        code_root, pretrained_model_dir, gene_mapping_file, asset_sources = resolve_geneformer_assets(
            project_root=project_root,
            model_dir=model_dir,
            cfg=fm_cfg,
            gene_dict_override=args.gene_dict_pkl,
        )
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "geneformer_model_assets_resolved",
                "pass": False,
                "details": {"error": str(exc), "model_dir": str(model_dir)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 6

    nproc_cfg = int(fm_cfg.get("nproc", 8))
    cell_chunk_cfg = int(fm_cfg.get("cell_chunk_size", 2048))
    tokenizer_chunk_cfg = int(fm_cfg.get("tokenizer_chunk_size", 512))
    fwd_bs_cfg = int(fm_cfg.get("forward_batch_size", 32))
    emb_mode_cfg = str(fm_cfg.get("emb_mode", "cls"))
    emb_layer_cfg = int(fm_cfg.get("emb_layer", -1))
    use_h5ad_index_cfg = bool(fm_cfg.get("use_h5ad_index", True))

    nproc = int(args.nproc) if args.nproc is not None else nproc_cfg
    cell_chunk_size = int(args.cell_chunk_size) if args.cell_chunk_size is not None else cell_chunk_cfg
    tokenizer_chunk_size = (
        int(args.tokenizer_chunk_size) if args.tokenizer_chunk_size is not None else tokenizer_chunk_cfg
    )
    forward_batch_size = (
        int(args.forward_batch_size) if args.forward_batch_size is not None else fwd_bs_cfg
    )
    emb_mode = str(args.emb_mode) if args.emb_mode is not None else emb_mode_cfg
    emb_layer = int(args.emb_layer) if args.emb_layer is not None else emb_layer_cfg
    use_h5ad_index = bool(args.use_h5ad_index or use_h5ad_index_cfg)

    assertions.append(
        {
            "name": "geneformer_model_assets_resolved",
            "pass": True,
            "details": {
                "rules": [
                    "Geneformer assets must be resolvable from --model-dir/config without hardcoded absolute paths"
                ],
                "model_dir": str(model_dir),
                "model_dir_source": model_dir_source,
                "code_root": str(code_root) if code_root else "python_env",
                "pretrained_model_dir": str(pretrained_model_dir),
                "gene_mapping_file": str(gene_mapping_file),
                "asset_sources": asset_sources,
                "nproc": int(nproc),
                "cell_chunk_size": int(cell_chunk_size),
                "tokenizer_chunk_size": int(tokenizer_chunk_size),
                "forward_batch_size": int(forward_batch_size),
                "emb_mode": emb_mode,
                "emb_layer": int(emb_layer),
                "use_h5ad_index": bool(use_h5ad_index),
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
        pretrained_model_dir,
        gene_mapping_file,
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
            pretrained_model_dir,
            gene_mapping_file,
        ]
    )
    if code_root is not None:
        input_paths.append(code_root)

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

    shared_gene_symbols = shared_var["gene_symbol"].fillna("").astype(str).tolist()

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
        transcriptome_tokenizer_cls, emb_extractor_cls, import_source = choose_geneformer_classes(code_root)
        symbol2ens = load_gene_symbol_to_ensembl(gene_mapping_file)
        keep_cols, kept_symbols, kept_ens = build_gene_mapping(shared_gene_symbols, symbol2ens)
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "geneformer_import_and_mapping",
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
            "name": "geneformer_import_and_mapping",
            "pass": True,
            "details": {
                "import_source": import_source,
                "n_input_genes": int(len(shared_gene_symbols)),
                "n_mappable_genes": int(len(keep_cols)),
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
            vectors, _invalid_cells, side_dim, stats, chunk_errors = extract_side_embeddings_geneformer(
                side=side,
                meta_path=meta_path,
                counts_path=counts_path,
                required_cell_ids=sorted(required_by_side.get(side, set())),
                keep_cols=keep_cols,
                kept_symbols=kept_symbols,
                kept_ens=kept_ens,
                transcriptome_tokenizer_cls=transcriptome_tokenizer_cls,
                emb_extractor_cls=emb_extractor_cls,
                pretrained_model_dir=pretrained_model_dir,
                nproc=nproc,
                cell_chunk_size=cell_chunk_size,
                tokenizer_chunk_size=tokenizer_chunk_size,
                forward_batch_size=forward_batch_size,
                emb_mode=emb_mode,
                emb_layer=emb_layer,
                use_h5ad_index=use_h5ad_index,
                force_retokenize=args.force_retokenize,
                stage_dir=stage_dir,
            )
        except Exception as exc:  # noqa: BLE001
            assertions.append(
                {
                    "name": f"geneformer_embedding_side_{side.lower()}",
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
            "code_root": str(code_root) if code_root else "python_env",
            "pretrained_model_dir": str(pretrained_model_dir),
            "gene_mapping_file": str(gene_mapping_file),
            "nproc": int(nproc),
            "cell_chunk_size": int(cell_chunk_size),
            "tokenizer_chunk_size": int(tokenizer_chunk_size),
            "forward_batch_size": int(forward_batch_size),
            "emb_mode": emb_mode,
            "emb_layer": int(emb_layer),
            "use_h5ad_index": bool(use_h5ad_index),
            "force_retokenize": bool(args.force_retokenize),
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

    allowed_roots = [task2_snapshot.resolve()]
    for p in [model_dir, pretrained_model_dir, gene_mapping_file.parent]:
        if p is not None:
            allowed_roots.append(p.resolve())
    if code_root is not None:
        allowed_roots.append(code_root.resolve())

    bad_inputs = [
        str(path.resolve())
        for path in sorted(set(input_paths))
        if not any(path.resolve().is_relative_to(root) for root in allowed_roots)
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "rules": [
                    "Extractor reads only from task2 snapshot inputs and resolved model assets"
                ],
                "task2_snapshot": str(task2_snapshot),
                "model_dir": str(model_dir),
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
        "script_path": str((project_root / "scripts/fm_extractors/extract_geneformer.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
            "model_dir": str(model_dir),
            "model_dir_source": model_dir_source,
            "code_root": str(code_root) if code_root else "python_env",
            "pretrained_model_dir": str(pretrained_model_dir),
            "gene_mapping_file": str(gene_mapping_file),
            "nproc": int(nproc),
            "cell_chunk_size": int(cell_chunk_size),
            "tokenizer_chunk_size": int(tokenizer_chunk_size),
            "forward_batch_size": int(forward_batch_size),
            "emb_mode": emb_mode,
            "emb_layer": int(emb_layer),
            "use_h5ad_index": bool(use_h5ad_index),
            "force_retokenize": bool(args.force_retokenize),
        },
        "inputs": [str(path.resolve()) for path in sorted(set(input_paths))],
        "outputs": output_paths,
        "summary": {
            "n_rows": n_rows,
            "embedding_dim": int(emb_dim),
            "n_valid": int(valid_mask.sum()),
            "n_invalid": int((~valid_mask).sum()),
            "invalid_reason_top": pd.Series(invalid_reason[invalid_reason != ""]).value_counts().head(10).to_dict(),
            "side_stats": side_stats,
            "mappable_gene_count": int(len(keep_cols)),
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
