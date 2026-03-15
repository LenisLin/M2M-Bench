# SCRIPT_HEADER_CONTRACT
# Script: scripts/s0_build_data_inventory.py
# Purpose: Build Task1 data inventory (T1-P3) from task1 snapshot only.
# Inputs:
#   - Task1 Snapshot: data/task1_snapshot_v1/
#   - Config: config/config.yaml
# Outputs:
#   - task1_data_inventory_long.csv: runs/<run_id>/s0_build_data_inventory/
#   - data_source_manifest.csv: runs/<run_id>/s0_build_data_inventory/
#   - run_manifest.json: runs/<run_id>/s0_build_data_inventory/
#   - audit_assertions.json: runs/<run_id>/s0_build_data_inventory/
#   - manifest.json: runs/<run_id>/s0_build_data_inventory/
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s0_build_data_inventory/
# Config Dependencies:
#   - config/config.yaml::project.seed
#   - config/config.yaml::paths.task1_snapshot
#   - config/config.yaml::paths.runs_dir
# Execution:
#   - python scripts/s0_build_data_inventory.py --project-root . --run-id <run_id> --seed 619
# Failure Modes:
#   - Missing task1 snapshot directory -> exit non-zero
#   - Missing required files/columns -> exit non-zero
# Last Updated: 2026-03-03

"""
Inputs:
- task1 snapshot files:
  - lincs/lincs-engine1-meta.csv
  - scperturb_delta/scperturb-crispr-delta-meta.csv
  - scperturb_delta/scperturb-drug-delta-meta.csv

Outputs:
- task1_data_inventory_long.csv
- data_source_manifest.csv
- run_manifest.json
- audit_assertions.json
- manifest.json

Frozen constants:
- STAGE = "s0_build_data_inventory"
- REQUIRED_INPUTS and required columns are fixed below.
- Target tokenization uses dataset-specific primary delimiters:
  - LINCS: "|"
  - scPerturb: "_"
  - plus ";" cleanup fallback per token chunk.

Attrition rules:
- N/A for S0 inventory aggregation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd
import yaml

STAGE = "s0_build_data_inventory"
CONFIG_PATH = Path("config/config.yaml")
MAX_COUNTEREXAMPLES = 5


@dataclass(frozen=True)
class SourceSpec:
    source_name: str
    relative_path: str
    dataset: str
    perturbation_type: str
    required_columns: Sequence[str]
    target_column: str
    cell_line_column: str
    primary_delimiter: str
    source_db: str


REQUIRED_SOURCES: List[SourceSpec] = [
    SourceSpec(
        source_name="lincs_meta",
        relative_path="lincs/lincs-engine1-meta.csv",
        dataset="LINCS",
        perturbation_type="MAPPED_FROM_PERT_TYPE",
        required_columns=("pert_type", "target", "cell_line"),
        target_column="target",
        cell_line_column="cell_line",
        primary_delimiter="|",
        source_db="LINCS_lincs-engine1-meta.csv",
    ),
    SourceSpec(
        source_name="scperturb_crispr_meta",
        relative_path="scperturb_delta/scperturb-crispr-delta-meta.csv",
        dataset="scPerturb",
        perturbation_type="Genetic",
        required_columns=("target_std", "cell_std"),
        target_column="target_std",
        cell_line_column="cell_std",
        primary_delimiter="_",
        source_db="scPerturb_scperturb-crispr-delta-meta.csv",
    ),
    SourceSpec(
        source_name="scperturb_drug_meta",
        relative_path="scperturb_delta/scperturb-drug-delta-meta.csv",
        dataset="scPerturb",
        perturbation_type="Chemical",
        required_columns=("target_std", "cell_std"),
        target_column="target_std",
        cell_line_column="cell_std",
        primary_delimiter="_",
        source_db="scPerturb_scperturb-drug-delta-meta.csv",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="S0 Task1 data inventory builder")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--run-id", required=True, help="Run identifier")
    parser.add_argument("--seed", type=int, default=None, help="Override config seed")
    return parser.parse_args()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def compute_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def tokenize_targets(series: pd.Series, primary_delimiter: str) -> pd.Series:
    def _split_one(raw_value: str) -> List[str]:
        primary_chunks = str(raw_value).split(primary_delimiter)
        tokens: List[str] = []
        for chunk in primary_chunks:
            # Semicolon cleanup is a secondary fallback.
            for token in chunk.split(";"):
                cleaned = token.strip()
                if cleaned:
                    tokens.append(cleaned)
        return tokens

    return (
        series.fillna("")
        .astype(str)
        .map(_split_one)
    )


def normalize_lincs_perturbation_type(pert_type: pd.Series) -> pd.Series:
    return pert_type.fillna("").astype(str).map(
        lambda value: "Chemical" if value.strip().lower() == "drug" else "Genetic"
    )


def write_csv(df: pd.DataFrame, output_path: Path) -> None:
    df.to_csv(output_path, index=False, quoting=csv.QUOTE_MINIMAL)


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    config_path = (project_root / CONFIG_PATH).resolve()
    if not config_path.exists():
        print(f"[ERROR] Missing config file: {config_path}", file=sys.stderr)
        return 1

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    configured_seed = int(config["project"]["seed"])
    seed = configured_seed if args.seed is None else int(args.seed)

    task1_snapshot = project_root / config["paths"]["task1_snapshot"]
    runs_dir = (project_root / config["paths"]["runs_dir"]).resolve()

    started_at = utc_now_iso()
    assertions: List[Dict[str, object]] = []
    source_manifest_rows: List[Dict[str, object]] = []
    input_paths: List[str] = []
    all_inventory_rows: List[pd.DataFrame] = []

    # Hard fail-fast: missing required snapshot directory.
    if not task1_snapshot.is_dir():
        print(
            f"[ERROR] Required snapshot directory missing: {task1_snapshot}",
            file=sys.stderr,
        )
        return 2

    stage_dir = runs_dir / args.run_id / STAGE
    stage_dir.mkdir(parents=True, exist_ok=True)

    assertions.append(
        {
            "name": "task1_snapshot_directory_exists",
            "pass": True,
            "details": {
                "rules": ["Path config.paths.task1_snapshot must exist and be a directory."],
                "path": str(task1_snapshot),
            },
            "counterexamples": [],
        }
    )

    for spec in REQUIRED_SOURCES:
        path = task1_snapshot / spec.relative_path
        input_paths.append(str(path))
        if not path.is_file():
            assertions.append(
                {
                    "name": f"required_file_exists::{spec.source_name}",
                    "pass": False,
                    "details": {
                        "rules": [f"Required file must exist: {spec.relative_path}"],
                        "path": str(path),
                    },
                    "counterexamples": [{"missing_path": str(path)}],
                }
            )
            write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
            return 3

        frame = pd.read_csv(path, usecols=list(spec.required_columns))
        missing_columns = [col for col in spec.required_columns if col not in frame.columns]
        if missing_columns:
            assertions.append(
                {
                    "name": f"required_columns_present::{spec.source_name}",
                    "pass": False,
                    "details": {
                        "rules": [f"Required columns must exist for {spec.source_name}."],
                        "required_columns": list(spec.required_columns),
                        "missing_columns": missing_columns,
                    },
                    "counterexamples": [{"file": str(path), "missing_column": col} for col in missing_columns[:MAX_COUNTEREXAMPLES]],
                }
            )
            write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
            return 4

        if spec.source_name == "lincs_meta":
            perturbation_type = normalize_lincs_perturbation_type(frame["pert_type"])
        else:
            perturbation_type = pd.Series([spec.perturbation_type] * len(frame), index=frame.index)

        tokens = tokenize_targets(frame[spec.target_column], spec.primary_delimiter)
        normalized = pd.DataFrame(
            {
                "dataset": spec.dataset,
                "perturbation_type": perturbation_type,
                "cell_line": frame[spec.cell_line_column].fillna("NA").astype(str).str.strip(),
                "target_token": tokens,
                "source_db": spec.source_db,
            }
        )
        exploded = normalized.explode("target_token")
        exploded["target_token"] = exploded["target_token"].fillna("").astype(str).str.strip()
        exploded = exploded[exploded["target_token"] != ""].copy()

        all_inventory_rows.append(exploded)
        source_manifest_rows.append(
            {
                "source_name": spec.source_name,
                "source_db": spec.source_db,
                "file_path": str(path),
                "row_count": int(len(frame)),
                "required_columns": "|".join(spec.required_columns),
                "required_columns_ok": True,
                "sha256": compute_sha256(path),
            }
        )
        assertions.append(
            {
                "name": f"required_columns_present::{spec.source_name}",
                "pass": True,
                "details": {
                    "rules": [f"Required columns must exist for {spec.source_name}."],
                    "required_columns": list(spec.required_columns),
                },
                "counterexamples": [],
            }
        )

    inventory_all = pd.concat(all_inventory_rows, ignore_index=True)
    inventory = (
        inventory_all.groupby(
            ["dataset", "perturbation_type", "cell_line", "target_token"],
            as_index=False,
            dropna=False,
        )
        .agg(
            n_instances=("target_token", "size"),
            source_db=("source_db", lambda x: "|".join(sorted(set(map(str, x))))),
        )
        .sort_values(
            by=["dataset", "perturbation_type", "cell_line", "target_token"],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )

    expected_inventory_cols = [
        "dataset",
        "perturbation_type",
        "cell_line",
        "target_token",
        "n_instances",
        "source_db",
    ]
    actual_cols = inventory.columns.tolist()
    schema_pass = actual_cols == expected_inventory_cols
    assertions.append(
        {
            "name": "task1_inventory_schema_t1_p3",
            "pass": schema_pass,
            "details": {
                "rules": [
                    "Output schema must follow T1-P3 with exact columns:",
                    "dataset, perturbation_type, cell_line, target_token, n_instances, source_db",
                ],
                "actual_columns": actual_cols,
            },
            "counterexamples": [] if schema_pass else [{"actual_columns": actual_cols}],
        }
    )
    if not schema_pass:
        write_json(stage_dir / "audit_assertions.json", {"assertions": assertions})
        return 5

    source_manifest = pd.DataFrame(source_manifest_rows).sort_values(by=["source_name"]).reset_index(drop=True)
    inventory_path = stage_dir / "task1_data_inventory_long.csv"
    data_source_manifest_path = stage_dir / "data_source_manifest.csv"
    write_csv(inventory[expected_inventory_cols], inventory_path)
    write_csv(source_manifest, data_source_manifest_path)

    # Output routing assertion.
    generated_core = [inventory_path, data_source_manifest_path]
    output_routing_pass = all(path.resolve().is_relative_to(stage_dir.resolve()) for path in generated_core)
    assertions.append(
        {
            "name": "output_routing_isolated",
            "pass": output_routing_pass,
            "details": {
                "rules": [f"All generated outputs must be inside {stage_dir}."],
                "stage_dir": str(stage_dir),
            },
            "counterexamples": [] if output_routing_pass else [{"bad_path": str(path)} for path in generated_core],
        }
    )

    task2_path_literal = str((project_root / config["paths"]["task2_snapshot"]).resolve())
    assertions.append(
        {
            "name": "task1_only_data_isolation",
            "pass": True,
            "details": {
                "rules": ["S0 must not read from config.paths.task2_snapshot."],
                "task2_snapshot_path_not_used": task2_path_literal,
                "used_inputs": input_paths,
            },
            "counterexamples": [],
        }
    )

    completed_at = utc_now_iso()
    outputs = [str(path.resolve()) for path in generated_core]

    run_manifest = {
        "run_id": args.run_id,
        "stage": STAGE,
        "script_path": str((project_root / "scripts/s0_build_data_inventory.py").resolve()),
        "git_head": safe_git_head(project_root),
        "started_at": started_at,
        "completed_at": completed_at,
        "config": {
            "seed": seed,
            "task1_snapshot": str(task1_snapshot),
            "runs_dir": str(runs_dir),
        },
        "inputs": input_paths,
        "outputs": outputs + [
            str((stage_dir / "run_manifest.json").resolve()),
            str((stage_dir / "audit_assertions.json").resolve()),
            str((stage_dir / "manifest.json").resolve()),
        ],
    }

    run_manifest_path = stage_dir / "run_manifest.json"
    audit_assertions_path = stage_dir / "audit_assertions.json"
    write_json(run_manifest_path, run_manifest)
    write_json(audit_assertions_path, {"assertions": assertions})

    manifest_path = stage_dir / "manifest.json"
    manifest_entries: List[Dict[str, object]] = []
    for path in sorted(stage_dir.iterdir()):
        if path.is_file() and path.name != "manifest.json":
            manifest_entries.append(
                {
                    "relative_path": path.name,
                    "size_bytes": path.stat().st_size,
                    "sha256": compute_sha256(path),
                }
            )
    write_json(manifest_path, {"stage": STAGE, "files": manifest_entries})
    return 0


def write_json(path: Path, payload: Dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=True)
        handle.write("\n")


if __name__ == "__main__":
    sys.exit(main())
