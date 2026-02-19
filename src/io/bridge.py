from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import pandas as pd
import yaml

from src.utils.logging import get_logger
from src.validation.contracts import DataContractError, validate_for_r_export


def _load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    return cfg


def save_for_r(
    df: pd.DataFrame,
    filename: str,
    config_path: Union[str, Path] = "config/config.yaml",
    primary_key: str = "row_id",
    provenance_script: str = "unknown",
    git_commit: str = "unknown",
) -> Tuple[Path, Path]:
    """
    Save a DataFrame for R visualization with an enforced meta sidecar.

    - Output directory: config.paths.interim_viz_dir (default: data/interim_viz)
    - Formats: .parquet or .csv only
    - Enforces explicit primary key column; creates from index if missing
    - Writes <stem>_meta.json with schema + provenance

    Returns:
        (data_path, meta_path)
    """
    logger = get_logger("src.io.bridge")

    cfg = _load_config(config_path)
    out_dir = Path(cfg.get("paths", {}).get("interim_viz_dir", "data/interim_viz"))
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / filename
    suffix = out_path.suffix.lower()
    if suffix not in {".parquet", ".csv"}:
        raise DataContractError("filename must end with .parquet or .csv")

    # Ensure explicit primary key
    if primary_key not in df.columns:
        logger.info(f"Primary key '{primary_key}' missing; creating from index via reset_index().")
        df = df.reset_index(drop=False)
        # If reset_index created a column named "index", rename it to primary_key
        if "index" in df.columns and primary_key not in df.columns:
            df = df.rename(columns={"index": primary_key})
        # If still missing, create sequential id
        if primary_key not in df.columns:
            df.insert(0, primary_key, range(df.shape[0]))

    validate_for_r_export(df, primary_key=primary_key)

    # Persist
    logger.info(f"Saving for R: rows={df.shape[0]} cols={df.shape[1]} -> {out_path}")
    if suffix == ".parquet":
        df.to_parquet(out_path, index=False)
    else:
        df.to_csv(out_path, index=False)

    # Sidecar meta
    meta = {
        "file": out_path.name,
        "primary_key": primary_key,
        "columns": {col: str(dtype) for col, dtype in df.dtypes.items()},
        "provenance": {
            "script": provenance_script,
            "git_commit": git_commit,
            "config": str(Path(config_path)),
        },
    }
    meta_path = out_path.with_name(f"{out_path.stem}_meta.json")
    with meta_path.open("w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    logger.info(f"Wrote meta sidecar: {meta_path}")
    return out_path, meta_path
