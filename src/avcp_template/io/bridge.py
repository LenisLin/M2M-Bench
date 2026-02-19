from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from ..utils.logging import get_logger
from ..validation.contracts import DataContractError, validate_for_r_export


def _load_config(config_path: str | Path) -> dict[str, Any]:
    resolved_path = Path(config_path).resolve()
    if not resolved_path.exists():
        raise FileNotFoundError(f"Config file not found: {resolved_path}")
    with resolved_path.open("r", encoding="utf-8") as file_obj:
        cfg = yaml.safe_load(file_obj) or {}
    return cfg


def _resolve_output_dir(cfg: dict[str, Any], config_path: Path) -> Path:
    configured_path = Path(cfg.get("paths", {}).get("interim_viz_dir", "data/interim_viz"))
    if configured_path.is_absolute():
        return configured_path

    if config_path.as_posix().endswith("config/config.yaml"):
        repo_root = config_path.resolve().parents[1]
    else:
        repo_root = Path.cwd()

    return repo_root / configured_path


def _ensure_primary_key(df: pd.DataFrame, primary_key: str) -> pd.DataFrame:
    logger = get_logger("avcp_template.io.bridge")
    if primary_key in df.columns:
        return df

    logger.info(
        "Primary key '%s' missing; creating from index via reset_index().",
        primary_key,
    )

    try:
        df_with_key = df.reset_index(drop=False, names=[primary_key])
    except TypeError:
        df_with_key = df.reset_index(drop=False)
        if "index" in df_with_key.columns and primary_key not in df_with_key.columns:
            df_with_key = df_with_key.rename(columns={"index": primary_key})

    if primary_key not in df_with_key.columns:
        df_with_key.insert(0, primary_key, range(df_with_key.shape[0]))

    return df_with_key


def save_for_r(
    df: pd.DataFrame,
    filename: str,
    config_path: str | Path = "config/config.yaml",
    primary_key: str = "row_id",
    provenance_script: str = "unknown",
    git_commit: str = "unknown",
) -> tuple[Path, Path]:
    """
    Save a DataFrame for R visualization with an enforced meta sidecar.

    - Output directory: config.paths.interim_viz_dir (default: data/interim_viz)
    - Formats: .parquet or .csv only
    - Enforces explicit primary key column; creates from index if missing
    - Writes <stem>_meta.json with schema + provenance

    Returns:
        (data_path, meta_path)
    """
    logger = get_logger("avcp_template.io.bridge")
    resolved_config_path = Path(config_path).resolve()

    cfg = _load_config(resolved_config_path)
    out_dir = _resolve_output_dir(cfg=cfg, config_path=resolved_config_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / filename
    suffix = out_path.suffix.lower()
    if suffix not in {".parquet", ".csv"}:
        raise DataContractError("filename must end with .parquet or .csv")

    df = _ensure_primary_key(df=df, primary_key=primary_key)
    validate_for_r_export(df, primary_key=primary_key)

    logger.info("Saving for R: rows=%s cols=%s -> %s", df.shape[0], df.shape[1], out_path)
    if suffix == ".parquet":
        df.to_parquet(out_path, index=False)
    else:
        df.to_csv(out_path, index=False)

    meta = {
        "file": out_path.name,
        "primary_key": primary_key,
        "columns": {col: str(dtype) for col, dtype in df.dtypes.items()},
        "provenance": {
            "script": provenance_script,
            "git_commit": git_commit,
            "config": str(resolved_config_path),
        },
    }
    meta_path = out_path.with_name(f"{out_path.stem}_meta.json")
    with meta_path.open("w", encoding="utf-8") as file_obj:
        json.dump(meta, file_obj, indent=2, ensure_ascii=False)

    logger.info("Wrote meta sidecar: %s", meta_path)
    return out_path, meta_path
