#!/usr/bin/env python3
"""
Task0 analysis-style pipeline runner.

Why this script exists
----------------------
Compared with thin wrapper scripts, this file provides one readable place for:
  1) Task0 curation execution
  2) optional Task0->Task1 attrition audit
  3) run manifest output

This keeps the workflow easy to follow and easy to modify.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from _script_bootstrap import setup_project_imports

setup_project_imports()

from m2m_bench.task0.curate import Task0Config, run_task0

# Reuse audit function from the dedicated script module.
from task01_attrition_audit import run_audit as run_task01_attrition_audit


@dataclass
class Task0PipelineConfig:
    # Core IO
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    output_dir: str = "./outputs/task0_curated"
    lincs_meta_name: str = "LINCS_Engine1_Meta.csv"
    lincs_tensor_name: str = "LINCS_Engine1_TrainData.pt"
    sc_meta_glob: str = "*Meta*.csv"

    # Execution switches
    run_curation: bool = True
    run_attrition_audit: bool = False
    dry_run: bool = False

    # Optional audit inputs
    task1_matched_pairs_path: str = "./outputs/task1/data/m1_matched_pairs.csv"
    audit_output_dir: str = "./outputs/task1_audit"

    # Selected curation knobs (frequently used)
    save_bundle: bool = True
    save_parquet: bool = True
    save_tensors_separately: bool = True
    strict_load: bool = True
    enable_qc_gates: bool = True
    max_pairs_per_group_level1: int | None = None
    qc_seed: int = 42

    @property
    def pipeline_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task0_pipeline.json")


def _resolve_table_path(path_like: str) -> Path:
    """
    Resolve a csv/parquet input path with fallback to sibling extension.
    """
    path = Path(path_like)
    if path.exists():
        return path

    if path.suffix.lower() == ".parquet":
        fallback = path.with_suffix(".csv")
        if fallback.exists():
            return fallback
    if path.suffix.lower() == ".csv":
        fallback = path.with_suffix(".parquet")
        if fallback.exists():
            return fallback

    raise FileNotFoundError(
        f"Cannot find matched-pairs table. Checked: {path} and sibling csv/parquet fallback."
    )


def _to_task0_config(cfg: Task0PipelineConfig) -> Task0Config:
    """Convert pipeline config to module-native Task0Config."""
    return Task0Config(
        processed_dir=cfg.processed_dir,
        output_dir=cfg.output_dir,
        lincs_meta_name=cfg.lincs_meta_name,
        lincs_tensor_name=cfg.lincs_tensor_name,
        sc_meta_glob=cfg.sc_meta_glob,
        save_bundle=cfg.save_bundle,
        save_parquet=cfg.save_parquet,
        save_tensors_separately=cfg.save_tensors_separately,
        strict_load=cfg.strict_load,
        enable_qc_gates=cfg.enable_qc_gates,
        max_pairs_per_group_level1=cfg.max_pairs_per_group_level1,
        qc_seed=cfg.qc_seed,
    )


def run_task0_pipeline(cfg: Task0PipelineConfig) -> None:
    t0 = time.time()
    summary = {
        "ran_curation": False,
        "ran_attrition_audit": False,
    }

    if cfg.dry_run:
        print("[DRY RUN] Task0 pipeline configuration:")
        print(json.dumps(asdict(cfg), indent=2))
        return

    # ------------------------------------------------------------
    # Step A. Task0 curation
    # ------------------------------------------------------------
    if cfg.run_curation:
        print(">>> [Task0 Pipeline] Step A: run Task0 curation")
        run_task0(_to_task0_config(cfg))
        summary["ran_curation"] = True

    # ------------------------------------------------------------
    # Step B. Optional Task0->Task1 attrition audit
    # ------------------------------------------------------------
    if cfg.run_attrition_audit:
        print(">>> [Task0 Pipeline] Step B: run Task0->Task1 attrition audit")
        bundle_path = Path(cfg.output_dir) / "bundle" / "m2m_task0_bundle.pt"
        unified_meta_path = Path(cfg.output_dir) / "metadata" / "unified_meta.parquet"
        matched_pairs_path = _resolve_table_path(cfg.task1_matched_pairs_path)
        run_task01_attrition_audit(
            bundle_path=bundle_path,
            unified_meta_path=unified_meta_path,
            matched_pairs_path=matched_pairs_path,
            output_dir=Path(cfg.audit_output_dir),
        )
        summary["ran_attrition_audit"] = True

    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": summary,
    }
    Path(cfg.pipeline_manifest_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f">>> [Task0 Pipeline] done. Manifest: {cfg.pipeline_manifest_path}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Task0 readable pipeline runner")

    parser.add_argument("--processed-dir", type=str, default=Task0PipelineConfig.processed_dir)
    parser.add_argument("--output-dir", type=str, default=Task0PipelineConfig.output_dir)
    parser.add_argument("--lincs-meta-name", type=str, default=Task0PipelineConfig.lincs_meta_name)
    parser.add_argument("--lincs-tensor-name", type=str, default=Task0PipelineConfig.lincs_tensor_name)
    parser.add_argument("--sc-meta-glob", type=str, default=Task0PipelineConfig.sc_meta_glob)

    parser.add_argument("--run-curation", action=argparse.BooleanOptionalAction, default=Task0PipelineConfig.run_curation)
    parser.add_argument(
        "--run-attrition-audit",
        action=argparse.BooleanOptionalAction,
        default=Task0PipelineConfig.run_attrition_audit,
    )
    parser.add_argument("--dry-run", action=argparse.BooleanOptionalAction, default=Task0PipelineConfig.dry_run)

    parser.add_argument("--task1-matched-pairs-path", type=str, default=Task0PipelineConfig.task1_matched_pairs_path)
    parser.add_argument("--audit-output-dir", type=str, default=Task0PipelineConfig.audit_output_dir)

    parser.add_argument("--save-bundle", action=argparse.BooleanOptionalAction, default=Task0PipelineConfig.save_bundle)
    parser.add_argument("--save-parquet", action=argparse.BooleanOptionalAction, default=Task0PipelineConfig.save_parquet)
    parser.add_argument(
        "--save-tensors-separately",
        action=argparse.BooleanOptionalAction,
        default=Task0PipelineConfig.save_tensors_separately,
    )
    parser.add_argument("--strict-load", action=argparse.BooleanOptionalAction, default=Task0PipelineConfig.strict_load)
    parser.add_argument(
        "--enable-qc-gates",
        action=argparse.BooleanOptionalAction,
        default=Task0PipelineConfig.enable_qc_gates,
    )
    parser.add_argument("--max-pairs-per-group-level1", type=int, default=Task0PipelineConfig.max_pairs_per_group_level1)
    parser.add_argument("--qc-seed", type=int, default=Task0PipelineConfig.qc_seed)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task0PipelineConfig(
        processed_dir=args.processed_dir,
        output_dir=args.output_dir,
        lincs_meta_name=args.lincs_meta_name,
        lincs_tensor_name=args.lincs_tensor_name,
        sc_meta_glob=args.sc_meta_glob,
        run_curation=args.run_curation,
        run_attrition_audit=args.run_attrition_audit,
        dry_run=args.dry_run,
        task1_matched_pairs_path=args.task1_matched_pairs_path,
        audit_output_dir=args.audit_output_dir,
        save_bundle=args.save_bundle,
        save_parquet=args.save_parquet,
        save_tensors_separately=args.save_tensors_separately,
        strict_load=args.strict_load,
        enable_qc_gates=args.enable_qc_gates,
        max_pairs_per_group_level1=args.max_pairs_per_group_level1,
        qc_seed=args.qc_seed,
    )
    run_task0_pipeline(cfg)


if __name__ == "__main__":
    main()
