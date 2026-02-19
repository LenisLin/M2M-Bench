from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.io.bridge import save_for_r


def test_save_for_r_writes_data_and_meta(tmp_path: Path) -> None:
    # Create an isolated config pointing interim_viz_dir to tmp_path
    cfg_dir = tmp_path / "config"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = cfg_dir / "config.yaml"
    cfg_path.write_text("paths:\n  interim_viz_dir: " + str((tmp_path / "viz").as_posix()) + "\n", encoding="utf-8")

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
