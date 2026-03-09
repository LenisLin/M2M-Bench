# SCRIPT_HEADER_CONTRACT
# Script: scripts/fm_extractors/extract_uce.py
# Purpose: Build Task2 K562 UCE FM deltas using pair_list with strict row alignment to delta_meta.
# Inputs:
#   - Task2 Snapshot K562: config/config.yaml::paths.task2_snapshot
#     - k562/derived/pair_list.parquet
#     - k562/derived/delta_meta.csv
#     - k562/{CRISPR_counts.pt,CRISPR_meta.csv,Drug_counts.pt,Drug_meta.csv,shared_var_names.csv}
#   - UCE repo/model dir:
#     - config/config.yaml::fm_extractors.uce.model_dir (optional)
#     - CLI override: --model-dir
# Outputs:
#   - fm_delta.npy: <paths.task2_snapshot>/k562/fm/uce/
#   - fm_delta_meta.csv: <paths.task2_snapshot>/k562/fm/uce/
#   - delta_operator_policy.json: <paths.task2_snapshot>/k562/fm/uce/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
#     under runs/<run_id>/extract_uce/
# Side Effects:
#   - Creates FM output directory under task2 snapshot
#   - Creates isolated run directory: runs/<run_id>/extract_uce/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::fm_extractors.uce (optional)
# Execution:
#   - python scripts/fm_extractors/extract_uce.py --run-id <run_id> [--model-dir <repo_dir>]
# Failure Modes:
#   - Missing required snapshot files -> exit non-zero
#   - Missing/invalid UCE assets -> exit non-zero
#   - Complete embedding failure for a side -> exit non-zero
# Last Updated: 2026-03-05

"""
Task2 UCE FM extractor with strict row-preservation contract:

1) Output fm_delta.npy must have exactly N rows where N = len(delta_meta).
2) Row i always corresponds to delta_meta row_id i (contiguous 0..N-1 required).
3) Any per-row embedding failure does NOT drop rows:
   - fm_delta[row_id, :] = NaN
   - fm_delta_meta.valid_mask[row_id] = False
4) Output routing is fixed:
   - Snapshot outputs -> data/task2_snapshot_v1/k562/fm/uce/
   - Run artifacts -> runs/<run_id>/extract_uce/
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

STAGE = "extract_uce"
MODEL_NAME = "uce"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
GLOBAL_SEED = 619
MAX_COUNTEREXAMPLES = 5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Task2 K562 UCE FM delta extractor")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="Override UCE root configured in config/config.yaml",
    )
    parser.add_argument("--python-bin", type=str, default=None, help="Python executable for eval_single_anndata.py")
    parser.add_argument("--species", type=str, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--cell-chunk-size", type=int, default=None)
    parser.add_argument("--nlayers", type=int, default=None)
    parser.add_argument("--model-loc", type=Path, default=None)
    parser.add_argument("--token-file", type=Path, default=None)
    parser.add_argument("--protein-embeddings-dir", type=Path, default=None)
    parser.add_argument("--spec-chrom-csv", type=Path, default=None)
    parser.add_argument("--offset-pkl", type=Path, default=None)
    parser.add_argument("--sample-size", type=int, default=None)
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


def resolve_uce_assets(
    *,
    project_root: Path,
    model_dir: Path,
    cfg: Mapping[str, object],
    args: argparse.Namespace,
) -> Tuple[Path, Path, Path, Path, Path, Path, Dict[str, str]]:
    sources: Dict[str, str] = {}
    root = model_dir.resolve()

    eval_script_candidates = [
        (root / "eval_single_anndata.py", "model_dir/eval_single_anndata.py"),
    ]
    eval_script: Optional[Path] = None
    for p, source in eval_script_candidates:
        if p.is_file():
            eval_script = p.resolve()
            sources["eval_script"] = source
            break
    if eval_script is None:
        raise RuntimeError("Cannot resolve eval_single_anndata.py from --model-dir")

    if args.model_loc is not None:
        model_loc = args.model_loc.resolve()
        sources["model_loc"] = "cli(--model-loc)"
    elif "model_loc" in cfg and cfg["model_loc"] is not None:
        model_loc = resolve_config_path(project_root, str(cfg["model_loc"]))
        sources["model_loc"] = "config.fm_extractors.uce.model_loc"
    else:
        model_loc = (root / "model_files/4layer_model.torch").resolve()
        sources["model_loc"] = "model_dir/model_files/4layer_model.torch"

    if args.token_file is not None:
        token_file = args.token_file.resolve()
        sources["token_file"] = "cli(--token-file)"
    elif "token_file" in cfg and cfg["token_file"] is not None:
        token_file = resolve_config_path(project_root, str(cfg["token_file"]))
        sources["token_file"] = "config.fm_extractors.uce.token_file"
    else:
        token_file = (root / "model_files/all_tokens.torch").resolve()
        sources["token_file"] = "model_dir/model_files/all_tokens.torch"

    if args.protein_embeddings_dir is not None:
        protein_embeddings_dir = args.protein_embeddings_dir.resolve()
        sources["protein_embeddings_dir"] = "cli(--protein-embeddings-dir)"
    elif "protein_embeddings_dir" in cfg and cfg["protein_embeddings_dir"] is not None:
        protein_embeddings_dir = resolve_config_path(project_root, str(cfg["protein_embeddings_dir"]))
        sources["protein_embeddings_dir"] = "config.fm_extractors.uce.protein_embeddings_dir"
    else:
        protein_embeddings_dir = (root / "model_files/protein_embeddings").resolve()
        sources["protein_embeddings_dir"] = "model_dir/model_files/protein_embeddings"

    if args.spec_chrom_csv is not None:
        spec_chrom_csv = args.spec_chrom_csv.resolve()
        sources["spec_chrom_csv"] = "cli(--spec-chrom-csv)"
    elif "spec_chrom_csv" in cfg and cfg["spec_chrom_csv"] is not None:
        spec_chrom_csv = resolve_config_path(project_root, str(cfg["spec_chrom_csv"]))
        sources["spec_chrom_csv"] = "config.fm_extractors.uce.spec_chrom_csv"
    else:
        spec_chrom_csv = (root / "model_files/species_chrom.csv").resolve()
        sources["spec_chrom_csv"] = "model_dir/model_files/species_chrom.csv"

    if args.offset_pkl is not None:
        offset_pkl = args.offset_pkl.resolve()
        sources["offset_pkl"] = "cli(--offset-pkl)"
    elif "offset_pkl" in cfg and cfg["offset_pkl"] is not None:
        offset_pkl = resolve_config_path(project_root, str(cfg["offset_pkl"]))
        sources["offset_pkl"] = "config.fm_extractors.uce.offset_pkl"
    else:
        offset_pkl = (root / "model_files/species_offsets.pkl").resolve()
        sources["offset_pkl"] = "model_dir/model_files/species_offsets.pkl"

    return (
        root,
        eval_script,
        model_loc,
        token_file,
        protein_embeddings_dir,
        spec_chrom_csv,
        offset_pkl,
        sources,
    )


def write_uce_input_h5ad(
    *,
    counts_dense: np.ndarray,
    gene_symbols: Sequence[str],
    cell_ids: Sequence[str],
    out_h5ad: Path,
) -> None:
    try:
        import anndata as ad  # type: ignore
        from scipy import sparse  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("anndata/scipy is required for UCE extractor") from exc

    x = counts_dense.astype(np.float32, copy=False)
    x_csr = sparse.csr_matrix(x)

    obs = pd.DataFrame(index=pd.Index([str(cid) for cid in cell_ids], name="cell_id"))
    var = pd.DataFrame(index=pd.Index([str(g) for g in gene_symbols], name="gene_symbol"))

    adata = ad.AnnData(X=x_csr, obs=obs, var=var)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)


def run_uce_once(
    *,
    python_bin: str,
    eval_script: Path,
    model_root: Path,
    input_h5ad: Path,
    work_dir: Path,
    species: str,
    model_loc: Path,
    batch_size: int,
    nlayers: int,
    sample_size: int,
    token_file: Path,
    protein_embeddings_dir: Path,
    spec_chrom_csv: Path,
    offset_pkl: Path,
) -> Tuple[bool, str]:
    work_dir.mkdir(parents=True, exist_ok=True)
    work_dir_arg = str(work_dir) + ("" if str(work_dir).endswith("/") else "/")

    cmd = [
        python_bin,
        str(eval_script),
        "--adata_path",
        str(input_h5ad),
        "--dir",
        work_dir_arg,
        "--species",
        str(species),
        "--model_loc",
        str(model_loc),
        "--batch_size",
        str(int(batch_size)),
        "--nlayers",
        str(int(nlayers)),
        "--sample_size",
        str(int(sample_size)),
        "--token_file",
        str(token_file),
        "--protein_embeddings_dir",
        str(protein_embeddings_dir),
        "--spec_chrom_csv_path",
        str(spec_chrom_csv),
        "--offset_pkl_path",
        str(offset_pkl),
    ]

    result = subprocess.run(
        cmd,
        cwd=str(model_root),
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        return True, ""

    msg = (result.stderr or "") + "\n" + (result.stdout or "")
    return False, msg[-3000:]


def find_uce_output_h5ad(work_dir: Path, input_h5ad: Path) -> Path:
    base = input_h5ad.stem
    expected = work_dir / f"{base}_uce_adata.h5ad"
    if expected.is_file():
        return expected.resolve()

    candidates = sorted(work_dir.glob("*_uce_adata.h5ad"), key=lambda p: p.stat().st_mtime, reverse=True)
    if candidates:
        return candidates[0].resolve()

    raise FileNotFoundError(f"No UCE output h5ad found under {work_dir}")


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


def extract_uce_aligned(
    *,
    output_h5ad: Path,
    expected_ids: Sequence[str],
    obsm_key: str = "X_uce",
) -> Tuple[np.ndarray, List[str]]:
    adata = read_h5ad(output_h5ad)

    if obsm_key not in adata.obsm_keys():
        raise KeyError(f"Expected embeddings in adata.obsm['{obsm_key}'], keys={list(adata.obsm_keys())}")

    emb = np.asarray(adata.obsm[obsm_key], dtype=np.float32)
    if emb.ndim != 2:
        raise RuntimeError(f"UCE embedding matrix must be 2D, got shape={emb.shape}")

    obs_ids = adata.obs_names.astype(str).tolist()
    emb_df = pd.DataFrame(emb, index=pd.Index(obs_ids, name="cell_id"))
    aligned = emb_df.reindex([str(x) for x in expected_ids])

    missing_mask = aligned.isna().any(axis=1)
    missing_ids = aligned.index[missing_mask].astype(str).tolist()
    matrix = aligned.to_numpy(dtype=np.float32, copy=False)
    return matrix, missing_ids


def run_uce_segment_with_retry(
    *,
    python_bin: str,
    eval_script: Path,
    model_root: Path,
    input_h5ad: Path,
    work_dir: Path,
    species: str,
    model_loc: Path,
    initial_batch_size: int,
    nlayers: int,
    sample_size: int,
    token_file: Path,
    protein_embeddings_dir: Path,
    spec_chrom_csv: Path,
    offset_pkl: Path,
    expected_ids: Sequence[str],
) -> Tuple[np.ndarray, int, Optional[str]]:
    batch_size = max(1, int(initial_batch_size))
    last_err: Optional[str] = None

    while True:
        ok, msg = run_uce_once(
            python_bin=python_bin,
            eval_script=eval_script,
            model_root=model_root,
            input_h5ad=input_h5ad,
            work_dir=work_dir,
            species=species,
            model_loc=model_loc,
            batch_size=batch_size,
            nlayers=nlayers,
            sample_size=sample_size,
            token_file=token_file,
            protein_embeddings_dir=protein_embeddings_dir,
            spec_chrom_csv=spec_chrom_csv,
            offset_pkl=offset_pkl,
        )
        if ok:
            output_h5ad = find_uce_output_h5ad(work_dir, input_h5ad)
            matrix, _ = extract_uce_aligned(output_h5ad=output_h5ad, expected_ids=expected_ids)
            return matrix, batch_size, None

        last_err = msg
        if batch_size <= 1 or not is_oom_text(msg):
            return np.empty((0, 0), dtype=np.float32), batch_size, last_err

        batch_size = max(1, batch_size // 2)
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


def extract_side_embeddings_uce(
    *,
    side: str,
    meta_path: Path,
    counts_path: Path,
    required_cell_ids: Sequence[str],
    gene_symbols: Sequence[str],
    python_bin: str,
    eval_script: Path,
    model_root: Path,
    species: str,
    model_loc: Path,
    batch_size: int,
    nlayers: int,
    sample_size: int,
    token_file: Path,
    protein_embeddings_dir: Path,
    spec_chrom_csv: Path,
    offset_pkl: Path,
    cell_chunk_size: int,
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

        queue: List[np.ndarray] = [np.arange(len(chunk_rows), dtype=np.int64)]
        segment_i = 0
        chunk_success = True

        while queue:
            seg_local = queue.pop(0)
            seg_ids = [chunk_ids[int(i)] for i in seg_local.tolist()]
            seg_rows = chunk_rows[seg_local]

            seg_tensor = torch.as_tensor(seg_rows, dtype=torch.long)
            seg_counts = counts_t.index_select(0, seg_tensor).detach().cpu().numpy().astype(np.float32, copy=False)

            with tempfile.TemporaryDirectory(prefix=f"uce_{side}_{start}_{end}_seg{segment_i:04d}_", dir=str(stage_dir)) as tmp:
                tmp_dir = Path(tmp)
                input_h5ad = tmp_dir / f"chunk_{side}_{start}_{end}_seg{segment_i:04d}.h5ad"
                work_dir = tmp_dir / "uce_work"
                segment_i += 1

                try:
                    write_uce_input_h5ad(
                        counts_dense=seg_counts,
                        gene_symbols=gene_symbols,
                        cell_ids=seg_ids,
                        out_h5ad=input_h5ad,
                    )

                    emb_matrix, used_bs, err = run_uce_segment_with_retry(
                        python_bin=python_bin,
                        eval_script=eval_script,
                        model_root=model_root,
                        input_h5ad=input_h5ad,
                        work_dir=work_dir,
                        species=species,
                        model_loc=model_loc,
                        initial_batch_size=batch_size,
                        nlayers=nlayers,
                        sample_size=sample_size,
                        token_file=token_file,
                        protein_embeddings_dir=protein_embeddings_dir,
                        spec_chrom_csv=spec_chrom_csv,
                        offset_pkl=offset_pkl,
                        expected_ids=seg_ids,
                    )

                    if emb_matrix.size == 0:
                        raise RuntimeError(f"UCE segment failed after retry (batch_size_final={used_bs}): {err}")

                    if emb_dim is None:
                        emb_dim = int(emb_matrix.shape[1])
                    elif int(emb_matrix.shape[1]) != emb_dim:
                        raise RuntimeError(
                            f"Embedding dim mismatch: got={emb_matrix.shape[1]} expected={emb_dim}"
                        )

                    for cid, vec in zip(seg_ids, emb_matrix):
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
        model_dir_source = "config.fm_extractors.uce.model_dir"
    else:
        model_dir = None
        model_dir_source = "missing"

    if model_dir is None or not model_dir.is_dir():
        assertions.append(
            {
                "name": "uce_model_dir_resolved",
                "pass": False,
                "details": {
                    "rules": ["UCE model dir must be provided by --model-dir or config.fm_extractors.uce.model_dir"],
                    "model_dir_source": model_dir_source,
                    "model_dir": str(model_dir) if model_dir else "NA",
                },
                "counterexamples": [{"model_dir": str(model_dir) if model_dir else "NA"}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Missing valid UCE model directory.", file=sys.stderr)
        return 5

    try:
        (
            model_root,
            eval_script,
            model_loc,
            token_file,
            protein_embeddings_dir,
            spec_chrom_csv,
            offset_pkl,
            asset_sources,
        ) = resolve_uce_assets(project_root=project_root, model_dir=model_dir, cfg=fm_cfg, args=args)
    except Exception as exc:  # noqa: BLE001
        assertions.append(
            {
                "name": "uce_assets_resolved",
                "pass": False,
                "details": {"error": str(exc), "model_dir": str(model_dir)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 6

    python_bin = args.python_bin if args.python_bin else str(fm_cfg.get("python_bin", "python"))
    species = str(args.species) if args.species is not None else str(fm_cfg.get("species", "human"))
    batch_size = int(args.batch_size) if args.batch_size is not None else int(fm_cfg.get("batch_size", 24))
    cell_chunk_size = (
        int(args.cell_chunk_size) if args.cell_chunk_size is not None else int(fm_cfg.get("cell_chunk_size", 2048))
    )
    nlayers = int(args.nlayers) if args.nlayers is not None else int(fm_cfg.get("nlayers", 4))
    sample_size = int(args.sample_size) if args.sample_size is not None else int(fm_cfg.get("sample_size", 1024))

    assertions.append(
        {
            "name": "uce_assets_resolved",
            "pass": True,
            "details": {
                "model_dir": str(model_dir),
                "model_dir_source": model_dir_source,
                "model_root": str(model_root),
                "eval_script": str(eval_script),
                "model_loc": str(model_loc),
                "token_file": str(token_file),
                "protein_embeddings_dir": str(protein_embeddings_dir),
                "spec_chrom_csv": str(spec_chrom_csv),
                "offset_pkl": str(offset_pkl),
                "asset_sources": asset_sources,
                "python_bin": python_bin,
                "species": species,
                "batch_size": int(batch_size),
                "cell_chunk_size": int(cell_chunk_size),
                "nlayers": int(nlayers),
                "sample_size": int(sample_size),
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
        eval_script,
        model_loc,
        token_file,
        protein_embeddings_dir,
        spec_chrom_csv,
        offset_pkl,
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
            model_root,
            eval_script,
            model_loc,
            token_file,
            protein_embeddings_dir,
            spec_chrom_csv,
            offset_pkl,
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
    gene_symbols = shared_var["gene_symbol"].fillna("").astype(str).tolist()

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
            vectors, _invalid_cells, side_dim, stats, chunk_errors = extract_side_embeddings_uce(
                side=side,
                meta_path=meta_path,
                counts_path=counts_path,
                required_cell_ids=sorted(required_by_side.get(side, set())),
                gene_symbols=gene_symbols,
                python_bin=python_bin,
                eval_script=eval_script,
                model_root=model_root,
                species=species,
                model_loc=model_loc,
                batch_size=batch_size,
                nlayers=nlayers,
                sample_size=sample_size,
                token_file=token_file,
                protein_embeddings_dir=protein_embeddings_dir,
                spec_chrom_csv=spec_chrom_csv,
                offset_pkl=offset_pkl,
                cell_chunk_size=cell_chunk_size,
                stage_dir=stage_dir,
            )
        except Exception as exc:  # noqa: BLE001
            assertions.append(
                {
                    "name": f"uce_embedding_side_{side.lower()}",
                    "pass": False,
                    "details": {"error": str(exc)},
                    "counterexamples": [{"error": str(exc)}],
                }
            )
            write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
            print(f"[ERROR] Failed embedding side={side}: {exc}", file=sys.stderr)
            return 11

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
            return 12

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
            "eval_script": str(eval_script),
            "model_loc": str(model_loc),
            "token_file": str(token_file),
            "protein_embeddings_dir": str(protein_embeddings_dir),
            "spec_chrom_csv": str(spec_chrom_csv),
            "offset_pkl": str(offset_pkl),
            "python_bin": python_bin,
            "species": species,
            "batch_size": int(batch_size),
            "cell_chunk_size": int(cell_chunk_size),
            "nlayers": int(nlayers),
            "sample_size": int(sample_size),
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
        eval_script.parent.absolute(),
        model_loc.parent.absolute(),
        token_file.parent.absolute(),
        protein_embeddings_dir.absolute(),
        spec_chrom_csv.parent.absolute(),
        offset_pkl.parent.absolute(),
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
        "script_path": str((project_root / "scripts/fm_extractors/extract_uce.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
            "model_dir": str(model_dir),
            "model_dir_source": model_dir_source,
            "model_root": str(model_root),
            "eval_script": str(eval_script),
            "model_loc": str(model_loc),
            "token_file": str(token_file),
            "protein_embeddings_dir": str(protein_embeddings_dir),
            "spec_chrom_csv": str(spec_chrom_csv),
            "offset_pkl": str(offset_pkl),
            "python_bin": python_bin,
            "species": species,
            "batch_size": int(batch_size),
            "cell_chunk_size": int(cell_chunk_size),
            "nlayers": int(nlayers),
            "sample_size": int(sample_size),
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
        return 13

    if not output_routing_snapshot_pass:
        print("[ERROR] Snapshot output routing assertion failed.", file=sys.stderr)
        return 14

    if bad_inputs:
        print("[ERROR] Input path isolation assertion failed.", file=sys.stderr)
        return 15

    if not (valid_finite_ok and invalid_nan_ok):
        print("[ERROR] Row alignment finite/NaN contract assertion failed.", file=sys.stderr)
        return 16

    return 0


if __name__ == "__main__":
    sys.exit(main())
