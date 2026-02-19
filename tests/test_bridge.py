from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from avcp_template.io.bridge import save_for_r
from avcp_template.validation.contracts import DataContractError


def test_save_for_r_writes_data_and_meta(tmp_path: Path) -> None:
    cfg_dir = tmp_path / "config"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = cfg_dir / "config.yaml"
    cfg_path.write_text(
        "paths:\n  interim_viz_dir: " + str((tmp_path / "viz").as_posix()) + "\n",
        encoding="utf-8",
    )

    df = pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]})

    data_path, meta_path = save_for_r(
        df=df,
        filename="example.parquet",
        config_path=cfg_path,
        primary_key="row_id",
        provenance_script="tests/test_bridge.py",
        git_commit="unknown",
    )

    assert data_path.exists()
    assert meta_path.exists()

    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    assert meta["file"] == data_path.name
    assert meta["primary_key"] == "row_id"
    assert "columns" in meta and isinstance(meta["columns"], dict)
    assert "provenance" in meta and isinstance(meta["provenance"], dict)
    assert meta["provenance"]["config"].endswith("config.yaml")


def test_save_for_r_resolves_relative_output_path_from_repo_root(tmp_path: Path) -> None:
    cfg_path = tmp_path / "config" / "config.yaml"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text("paths:\n  interim_viz_dir: data/interim_viz\n", encoding="utf-8")

    df = pd.DataFrame({"row_id": [1, 2], "value": [10, 20]})
    data_path, meta_path = save_for_r(
        df=df,
        filename="example.csv",
        config_path=cfg_path,
        primary_key="row_id",
    )

    expected = tmp_path / "data" / "interim_viz" / "example.csv"
    assert data_path == expected
    assert data_path.exists()
    assert meta_path.exists()


def test_save_for_r_rejects_invalid_extension(tmp_path: Path) -> None:
    cfg_path = tmp_path / "config" / "config.yaml"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text("paths:\n  interim_viz_dir: data/interim_viz\n", encoding="utf-8")

    df = pd.DataFrame({"row_id": [1], "value": [10]})
    with pytest.raises(DataContractError, match="filename must end with .parquet or .csv"):
        save_for_r(df=df, filename="bad.json", config_path=cfg_path, primary_key="row_id")
