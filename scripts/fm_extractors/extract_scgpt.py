# SCRIPT_HEADER_CONTRACT
# Script: scripts/fm_extractors/extract_scgpt.py
# Purpose: Build Task2 K562 scGPT FM deltas using pair_list with strict row alignment to delta_meta.
# Inputs:
#   - Task2 Snapshot K562: config/config.yaml::paths.task2_snapshot
#     - k562/derived/pair_list.parquet
#     - k562/derived/delta_meta.csv
#     - k562/{CRISPR_counts.pt,CRISPR_meta.csv,Drug_counts.pt,Drug_meta.csv,shared_var_names.csv}
#   - scGPT model checkpoint dir:
#     - config/config.yaml::fm_extractors.scgpt.model_dir (optional)
#     - CLI override: --model-dir
# Outputs:
#   - fm_delta.npy: <paths.task2_snapshot>/k562/fm/scgpt/
#   - fm_delta_meta.csv: <paths.task2_snapshot>/k562/fm/scgpt/
#   - delta_operator_policy.json: <paths.task2_snapshot>/k562/fm/scgpt/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
#     under runs/<run_id>/extract_scgpt/
# Side Effects:
#   - Creates FM output directory under task2 snapshot
#   - Creates isolated run directory: runs/<run_id>/extract_scgpt/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task2_snapshot
#   - config/config.yaml::paths.runs_dir
#   - config/config.yaml::fm_extractors.scgpt (optional)
# Execution:
#   - python scripts/fm_extractors/extract_scgpt.py --run-id <run_id> [--model-dir <checkpoint_dir>]
# Failure Modes:
#   - Missing required snapshot files -> exit non-zero
#   - Missing/invalid scGPT model directory -> exit non-zero
#   - Complete embedding failure for a side -> exit non-zero
# Last Updated: 2026-03-05

"""
Task2 scGPT FM extractor with strict row-preservation contract:

1) Output fm_delta.npy must have exactly N rows where N = len(delta_meta).
2) Row i always corresponds to delta_meta row_id i (contiguous 0..N-1 required).
3) Any per-row embedding failure does NOT drop rows:
   - fm_delta[row_id, :] = NaN
   - fm_delta_meta.valid_mask[row_id] = False
4) Output routing is fixed:
   - Snapshot outputs -> data/task2_snapshot_v1/k562/fm/scgpt/
   - Run artifacts -> runs/<run_id>/extract_scgpt/
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import os
import random
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import yaml

STAGE = "extract_scgpt"
CONFIG_PATH = Path("config/config.yaml")
EXPECTED_TASK2_SNAPSHOT = Path("data/task2_snapshot_v1")
GLOBAL_SEED = 619
MAX_COUNTEREXAMPLES = 5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Task2 K562 scGPT FM delta extractor")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="Override scGPT checkpoint directory configured in config/config.yaml",
    )
    parser.add_argument("--device", type=str, default=None, help="Override device, e.g. cuda or cpu")
    parser.add_argument("--batch-size", type=int, default=None, help="Override scGPT batch size")
    parser.add_argument("--chunk-size", type=int, default=None, help="Override cell chunk size")
    parser.add_argument("--max-length", type=int, default=None, help="Override max sequence length")
    parser.add_argument(
        "--use-fast-transformer",
        action="store_true",
        help="Enable scGPT fast transformer path if supported by installed version",
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


def choose_scgpt_embed_fn() -> Any:
    try:
        import scgpt as scg  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "Unable to import scgpt. Run this script in the dedicated scGPT conda environment."
        ) from exc

    if hasattr(scg, "tasks") and hasattr(scg.tasks, "embed_data"):
        return scg.tasks.embed_data

    try:
        from scgpt.tasks.cell_emb import embed_data  # type: ignore

        return embed_data
    except Exception as exc:
        raise RuntimeError("Unable to locate scGPT embed_data function.") from exc


def call_scgpt_embed_data(
    *,
    embed_fn: Any,
    adata: Any,
    model_dir: Path,
    max_length: int,
    batch_size: int,
    device: str,
    use_fast_transformer: bool,
) -> Any:
    sig = inspect.signature(embed_fn)
    params = sig.parameters

    kwargs: Dict[str, Any] = {}
    pos_args: List[Any] = []

    if "adata_or_file" in params:
        kwargs["adata_or_file"] = adata
    elif "adata" in params:
        kwargs["adata"] = adata
    else:
        pos_args.append(adata)

    if "model_dir" in params:
        kwargs["model_dir"] = str(model_dir)
    elif "model_path" in params:
        kwargs["model_path"] = str(model_dir)
    else:
        pos_args.append(str(model_dir))

    if "gene_col" in params:
        kwargs["gene_col"] = "gene_symbol"
    if "cell_type_key" in params:
        kwargs["cell_type_key"] = "cell_type"
    if "max_length" in params:
        kwargs["max_length"] = int(max_length)
    if "batch_size" in params:
        kwargs["batch_size"] = int(batch_size)
    if "device" in params:
        kwargs["device"] = device
    if "return_new_adata" in params:
        kwargs["return_new_adata"] = True
    if "use_fast_transformer" in params:
        kwargs["use_fast_transformer"] = bool(use_fast_transformer)

    return embed_fn(*pos_args, **kwargs)


def to_dense_2d_matrix(output_obj: Any) -> np.ndarray:
    x = output_obj.X if hasattr(output_obj, "X") else output_obj
    try:
        import scipy.sparse as sp  # type: ignore

        if sp.issparse(x):
            x = x.toarray()
    except Exception:
        pass

    arr = np.asarray(x, dtype=np.float32)
    if arr.ndim != 2:
        raise RuntimeError(f"scGPT embedding output must be 2D, got shape={arr.shape}")
    return arr


def is_oom_error(exc: BaseException) -> bool:
    text = str(exc).lower()
    return "out of memory" in text or ("cuda" in text and "memory" in text)


def embed_chunk_with_retry(
    *,
    embed_fn: Any,
    adata_chunk: Any,
    model_dir: Path,
    max_length: int,
    initial_batch_size: int,
    device: str,
    use_fast_transformer: bool,
) -> Tuple[Optional[np.ndarray], int, Optional[str]]:
    batch_size = max(1, int(initial_batch_size))
    last_error: Optional[str] = None

    while True:
        try:
            output_obj = call_scgpt_embed_data(
                embed_fn=embed_fn,
                adata=adata_chunk,
                model_dir=model_dir,
                max_length=max_length,
                batch_size=batch_size,
                device=device,
                use_fast_transformer=use_fast_transformer,
            )
            matrix = to_dense_2d_matrix(output_obj)
            return matrix, batch_size, None
        except Exception as exc:  # noqa: BLE001
            last_error = str(exc)
            if batch_size <= 1 or not is_oom_error(exc):
                return None, batch_size, last_error
            batch_size = max(1, batch_size // 2)
            if torch.cuda.is_available():
                torch.cuda.empty_cache()


def choose_side_id_column(meta: pd.DataFrame, side: str) -> str:
    if side == "CRISPR":
        candidates = ["cell_barcode", "cell_id", "Unnamed: 0"]
    else:
        candidates = ["Unnamed: 0", "cell_id", "cell_barcode"]
    for column in candidates:
        if column in meta.columns:
            return column
    raise ValueError(f"Unable to locate cell-id column for side={side}")


def build_adata_chunk(
    *,
    chunk_counts: np.ndarray,
    chunk_meta: pd.DataFrame,
    gene_symbols: Sequence[str],
) -> Any:
    try:
        import anndata as ad  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "Unable to import anndata. Run this script in the dedicated scGPT conda environment."
        ) from exc

    obs = chunk_meta.copy()
    if "cell_type" not in obs.columns:
        obs["cell_type"] = "NA"

    adata = ad.AnnData(X=chunk_counts.astype(np.float32, copy=False), obs=obs)
    adata.var["gene_symbol"] = [str(g) for g in gene_symbols]
    return adata


def extract_side_embeddings(
    *,
    side: str,
    meta_path: Path,
    counts_path: Path,
    required_cell_ids: Sequence[str],
    gene_symbols: Sequence[str],
    embed_fn: Any,
    model_dir: Path,
    max_length: int,
    batch_size: int,
    chunk_size: int,
    device: str,
    use_fast_transformer: bool,
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

    for start in range(0, len(ordered_rows), max(1, int(chunk_size))):
        end = min(start + max(1, int(chunk_size)), len(ordered_rows))
        chunk_rows = ordered_rows[start:end]
        chunk_ids = ordered_ids[start:end]
        stats["n_chunks_total"] += 1

        row_tensor = torch.as_tensor(chunk_rows, dtype=torch.long)
        chunk_counts = (
            counts_t.index_select(0, row_tensor)
            .detach()
            .cpu()
            .numpy()
            .astype(np.float32, copy=False)
        )
        chunk_meta = meta.iloc[chunk_rows].copy().reset_index(drop=True)
        chunk_meta["cell_id"] = chunk_ids

        adata_chunk = build_adata_chunk(
            chunk_counts=chunk_counts,
            chunk_meta=chunk_meta,
            gene_symbols=gene_symbols,
        )
        matrix, used_bs, error_msg = embed_chunk_with_retry(
            embed_fn=embed_fn,
            adata_chunk=adata_chunk,
            model_dir=model_dir,
            max_length=max_length,
            initial_batch_size=batch_size,
            device=device,
            use_fast_transformer=use_fast_transformer,
        )

        if matrix is None:
            stats["n_chunks_failed"] += 1
            invalid_cells.update(chunk_ids)
            msg = (
                f"{side} chunk[{start}:{end}] failed after retry "
                f"(batch_size_final={used_bs}): {error_msg}"
            )
            chunk_errors.append(msg)
            continue

        if matrix.shape[0] != len(chunk_ids):
            stats["n_chunks_failed"] += 1
            invalid_cells.update(chunk_ids)
            chunk_errors.append(
                f"{side} chunk[{start}:{end}] row mismatch: emb_rows={matrix.shape[0]}, "
                f"chunk_rows={len(chunk_ids)}"
            )
            continue

        if emb_dim is None:
            emb_dim = int(matrix.shape[1])
        elif int(matrix.shape[1]) != emb_dim:
            stats["n_chunks_failed"] += 1
            invalid_cells.update(chunk_ids)
            chunk_errors.append(
                f"{side} chunk[{start}:{end}] dim mismatch: got={matrix.shape[1]}, expected={emb_dim}"
            )
            continue

        stats["n_chunks_succeeded"] += 1
        for cell_id, vec in zip(chunk_ids, matrix):
            if not np.isfinite(vec).all():
                invalid_cells.add(cell_id)
                continue
            vectors_by_cell[cell_id] = np.asarray(vec, dtype=np.float32)

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    if emb_dim is None:
        raise RuntimeError(
            f"{side} embedding failed for all chunks; cannot determine embedding dimension."
        )

    return vectors_by_cell, sorted(invalid_cells), emb_dim, stats, chunk_errors


def load_max_length(model_dir: Path, cfg: Mapping[str, object]) -> Tuple[int, Optional[Path], str]:
    if "max_length" in cfg and cfg["max_length"] is not None:
        return int(cfg["max_length"]), None, "config.fm_extractors.scgpt.max_length"

    args_json = model_dir / "args.json"
    if args_json.is_file():
        with args_json.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        value = payload.get("max_seq_len", payload.get("max_length", 1200))
        return int(value), args_json.resolve(), "model_dir/args.json"

    return 1200, None, "default_1200"


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
    try:
        from scgpt.utils import set_seed as scgpt_set_seed  # type: ignore

        scgpt_set_seed(GLOBAL_SEED)
    except Exception:
        pass

    assertions.append(
        {
            "name": "global_seed_locked",
            "pass": True,
            "details": {"rules": ["GLOBAL_SEED must equal 619"], "seed": GLOBAL_SEED},
            "counterexamples": [],
        }
    )

    scgpt_cfg_raw = config.get("fm_extractors", {})
    scgpt_cfg: Mapping[str, object]
    if isinstance(scgpt_cfg_raw, Mapping):
        node = scgpt_cfg_raw.get("scgpt", {})
        scgpt_cfg = node if isinstance(node, Mapping) else {}
    else:
        scgpt_cfg = {}

    model_dir: Optional[Path]
    if args.model_dir is not None:
        model_dir = args.model_dir.resolve()
        model_dir_source = "cli(--model-dir)"
    elif "model_dir" in scgpt_cfg and scgpt_cfg["model_dir"] is not None:
        model_dir = resolve_config_path(project_root, str(scgpt_cfg["model_dir"]))
        model_dir_source = "config.fm_extractors.scgpt.model_dir"
    else:
        model_dir = None
        model_dir_source = "missing"

    if model_dir is None or not model_dir.is_dir():
        assertions.append(
            {
                "name": "scgpt_model_dir_resolved",
                "pass": False,
                "details": {
                    "rules": [
                        "scGPT checkpoint dir must be provided by --model-dir or config.fm_extractors.scgpt.model_dir"
                    ],
                    "model_dir_source": model_dir_source,
                    "model_dir": str(model_dir) if model_dir else "NA",
                },
                "counterexamples": [{"model_dir": str(model_dir) if model_dir else "NA"}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print("[ERROR] Missing valid scGPT model directory.", file=sys.stderr)
        return 5

    max_length_cfg, args_json_path, max_length_source = load_max_length(model_dir, scgpt_cfg)
    max_length = int(args.max_length) if args.max_length is not None else int(max_length_cfg)
    device_cfg = str(scgpt_cfg.get("device", "auto"))
    batch_size_cfg = int(scgpt_cfg.get("batch_size", 96))
    chunk_size_cfg = int(scgpt_cfg.get("chunk_size", 2048))
    use_fast_cfg = bool(scgpt_cfg.get("use_fast_transformer", False))

    device = args.device if args.device else device_cfg
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"

    batch_size = int(args.batch_size) if args.batch_size is not None else batch_size_cfg
    chunk_size = int(args.chunk_size) if args.chunk_size is not None else chunk_size_cfg
    use_fast_transformer = bool(args.use_fast_transformer or use_fast_cfg)

    assertions.append(
        {
            "name": "scgpt_model_dir_resolved",
            "pass": True,
            "details": {
                "rules": [
                    "scGPT checkpoint dir must be provided by --model-dir or config.fm_extractors.scgpt.model_dir"
                ],
                "model_dir": str(model_dir),
                "model_dir_source": model_dir_source,
                "max_length": int(max_length),
                "max_length_source": max_length_source,
                "batch_size": int(batch_size),
                "chunk_size": int(chunk_size),
                "device": device,
                "use_fast_transformer": bool(use_fast_transformer),
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
    ]
    missing_inputs = [path for path in required_inputs if not path.is_file()]
    if missing_inputs:
        assertions.append(
            {
                "name": "task2_snapshot_inputs_present",
                "pass": False,
                "details": {"rules": ["All required task2 snapshot inputs must exist."]},
                "counterexamples": [{"missing_input": str(p)} for p in missing_inputs[:MAX_COUNTEREXAMPLES]],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] Missing required inputs: {missing_inputs}", file=sys.stderr)
        return 6

    input_paths.extend(required_inputs)
    input_paths.append(model_dir)
    if args_json_path is not None and args_json_path.is_file():
        input_paths.append(args_json_path)

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
    except Exception as exc:
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
        return 7

    ensure_required_columns(pair_df, ["treated_cell_id", "control_cell_id", "control_rank", "n_controls_used", "dataset_side"], "pair_list.parquet")
    ensure_required_columns(delta_meta, ["row_id", "treated_cell_id", "n_controls_used"], "delta_meta.csv")
    if "dataset_side" not in delta_meta.columns:
        if "perturbation_class" not in delta_meta.columns:
            raise ValueError("delta_meta.csv requires dataset_side or perturbation_class")
        mapped = delta_meta["perturbation_class"].map(
            lambda x: "CRISPR" if normalize_scalar(x).lower() == "genetic" else "DRUG"
        )
        delta_meta["dataset_side"] = mapped

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
        return 8

    pair_sorted = pair_df.sort_values(
        ["dataset_side", "treated_cell_id", "control_rank"], kind="mergesort"
    )
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
        return 9

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
        embed_fn = choose_scgpt_embed_fn()
    except Exception as exc:
        assertions.append(
            {
                "name": "scgpt_import_and_embed_function",
                "pass": False,
                "details": {"error": str(exc)},
                "counterexamples": [{"error": str(exc)}],
            }
        )
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 10

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
            vectors, _invalid_cells, side_dim, stats, chunk_errors = extract_side_embeddings(
                side=side,
                meta_path=meta_path,
                counts_path=counts_path,
                required_cell_ids=sorted(required_by_side.get(side, set())),
                gene_symbols=gene_symbols,
                embed_fn=embed_fn,
                model_dir=model_dir,
                max_length=max_length,
                batch_size=batch_size,
                chunk_size=chunk_size,
                device=device,
                use_fast_transformer=use_fast_transformer,
            )
        except Exception as exc:
            assertions.append(
                {
                    "name": f"scgpt_embedding_side_{side.lower()}",
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
            else [
                {
                    "valid_finite_ok": valid_finite_ok,
                    "invalid_nan_ok": invalid_nan_ok,
                }
            ],
        }
    )

    fm_root = (k562_dir / "fm/scgpt").resolve()
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
            "model_name": np.full(n_rows, "scgpt", dtype=object),
            "seed": np.full(n_rows, GLOBAL_SEED, dtype=np.int64),
        }
    )
    write_csv(fm_meta, fm_meta_path)

    policy_payload = {
        "model_name": "scgpt",
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
            "max_length": int(max_length),
            "batch_size": int(batch_size),
            "chunk_size": int(chunk_size),
            "device": device,
            "use_fast_transformer": bool(use_fast_transformer),
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

    allowed_roots = [task2_snapshot.resolve(), model_dir.resolve()]
    bad_inputs = [
        str(path.resolve())
        for path in sorted(set(input_paths))
        if not any(path.absolute().is_relative_to(root) for root in allowed_roots)
    ]
    assertions.append(
        {
            "name": "input_path_isolation",
            "pass": len(bad_inputs) == 0,
            "details": {
                "rules": [
                    "Extractor reads only from task2 snapshot inputs and specified model_dir"
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

    output_paths = [
        str(path.resolve())
        for path in (snapshot_outputs + stage_outputs)
    ]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/fm_extractors/extract_scgpt.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": GLOBAL_SEED,
            "task2_snapshot": str(task2_snapshot),
            "runs_dir": str(runs_dir),
            "model_dir": str(model_dir),
            "model_dir_source": model_dir_source,
            "max_length": int(max_length),
            "max_length_source": max_length_source,
            "batch_size": int(batch_size),
            "chunk_size": int(chunk_size),
            "device": device,
            "use_fast_transformer": bool(use_fast_transformer),
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

    output_routing_stage_pass = all(
        path.resolve().is_relative_to(stage_dir.resolve()) for path in stage_outputs
    )
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
