#!/usr/bin/env python3
"""
Task1 analysis-style pipeline runner.

This script provides a readable and configurable Task1 execution flow:
  1) Group-wise modality gap
  2) Instance-level retrieval
  3) Confounder analysis
  4) Set-level centroid analysis

Each step can run independently via flags.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from _script_bootstrap import setup_project_imports

setup_project_imports()

from m2m_bench.task1.groupwise_gap import Task1Config, run_task1
from m2m_bench.task1.retrieval_instance import Task1RetrievalConfig, run_task1_retrieval
from m2m_bench.task1.confounder_analysis import Task1ConfounderConfig, run_task1_confounder
from m2m_bench.task1.set_level_centroid import Task1SetLevelConfig, run_task1_set_level


@dataclass
class Task1PipelineConfig:
    # shared paths
    task0_unified_meta_path: str = "./outputs/task0_curated/metadata/unified_meta.parquet"
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    lincs_tensor_name: str = "LINCS_Engine1_TrainData.pt"
    task1_output_dir: str = "./outputs/task1"

    # step switches
    run_groupwise: bool = True
    run_retrieval: bool = True
    run_confounder: bool = True
    run_set_level: bool = True
    dry_run: bool = False

    # groupwise params
    max_pairs_total: int | None = 200000
    n_perm_groupwise: int = 200
    n_bootstrap_groupwise: int = 200
    max_groups_per_type: int = 300
    min_group_n_for_stats: int = 20
    strict_enable: bool = True
    strict_chemical_max_dose_logdiff: float | None = 0.5
    strict_chemical_max_time_absdiff: float | None = 24.0

    # retrieval params
    retrieval_output_dir: str = "./outputs/task1/retrieval"
    retrieval_tracks: str = "gene,path"
    retrieval_directions: str = "LINCS->scPerturb,scPerturb->LINCS"
    retrieval_topk: str = "1,5,10"
    retrieval_batch_size: int = 256
    retrieval_n_perm: int = 200
    retrieval_min_queries_per_group: int = 10
    retrieval_run_balanced_eval: bool = True
    retrieval_balanced_gallery_size: int = 256
    retrieval_balanced_true_per_query: int = 1
    retrieval_balanced_n_repeats: int = 50

    # confounder params
    confounder_output_dir: str = "./outputs/task1/confounder"
    confounder_n_perm: int = 200
    confounder_min_group_n: int = 10
    confounder_max_groups_for_report: int = 100

    # set-level params
    set_level_output_dir: str = "./outputs/task1/retrieval_set_level"
    set_level_min_n_per_source: int = 1
    set_level_n_perm: int = 100
    set_level_n_bootstrap: int = 100

    random_seed: int = 42

    @property
    def pipeline_manifest_path(self) -> str:
        return str(Path(self.task1_output_dir) / "run_manifest_task1_pipeline.json")


def _resolve_with_sibling_fallback(path_like: str) -> Path:
    """
    Resolve a table path with csv/parquet sibling fallback.
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

    raise FileNotFoundError(f"Missing required input: {path} (and sibling csv/parquet fallback not found).")


def _resolve_candidates_path(base_task1_dir: str) -> Path:
    """Prefer parquet; fallback to csv when parquet engine is unavailable elsewhere."""
    parquet_path = Path(base_task1_dir) / "data" / "m1_candidates.parquet"
    return _resolve_with_sibling_fallback(str(parquet_path))


def run_task1_pipeline(cfg: Task1PipelineConfig) -> None:
    t0 = time.time()
    summary = {
        "ran_groupwise": False,
        "ran_retrieval": False,
        "ran_confounder": False,
        "ran_set_level": False,
    }

    if cfg.dry_run:
        print("[DRY RUN] Task1 pipeline configuration:")
        print(json.dumps(asdict(cfg), indent=2))
        return

    # ------------------------------------------------------------
    # Step A. Group-wise modality gap
    # ------------------------------------------------------------
    if cfg.run_groupwise:
        print(">>> [Task1 Pipeline] Step A: group-wise modality gap")
        group_cfg = Task1Config(
            task0_unified_meta_path=cfg.task0_unified_meta_path,
            processed_dir=cfg.processed_dir,
            lincs_tensor_name=cfg.lincs_tensor_name,
            output_dir=cfg.task1_output_dir,
            max_pairs_total=cfg.max_pairs_total,
            min_group_n_for_stats=cfg.min_group_n_for_stats,
            max_groups_per_type=cfg.max_groups_per_type,
            n_perm=cfg.n_perm_groupwise,
            n_bootstrap=cfg.n_bootstrap_groupwise,
            random_seed=cfg.random_seed,
            enable_strict_protocol_subset=cfg.strict_enable,
            strict_chemical_max_dose_logdiff=cfg.strict_chemical_max_dose_logdiff,
            strict_chemical_max_time_absdiff=cfg.strict_chemical_max_time_absdiff,
        )
        run_task1(group_cfg)
        summary["ran_groupwise"] = True

    # ------------------------------------------------------------
    # Step B. Retrieval
    # ------------------------------------------------------------
    if cfg.run_retrieval:
        print(">>> [Task1 Pipeline] Step B: retrieval")
        candidates_path = _resolve_candidates_path(cfg.task1_output_dir)
        retrieval_cfg = Task1RetrievalConfig(
            m1_candidates_path=str(candidates_path),
            processed_dir=cfg.processed_dir,
            lincs_tensor_name=cfg.lincs_tensor_name,
            output_dir=cfg.retrieval_output_dir,
            tracks=cfg.retrieval_tracks,
            directions=cfg.retrieval_directions,
            topk=cfg.retrieval_topk,
            batch_size=cfg.retrieval_batch_size,
            n_perm=cfg.retrieval_n_perm,
            min_queries_per_group=cfg.retrieval_min_queries_per_group,
            random_seed=cfg.random_seed,
            run_balanced_eval=cfg.retrieval_run_balanced_eval,
            balanced_gallery_size=cfg.retrieval_balanced_gallery_size,
            balanced_true_per_query=cfg.retrieval_balanced_true_per_query,
            balanced_n_repeats=cfg.retrieval_balanced_n_repeats,
        )
        run_task1_retrieval(retrieval_cfg)
        summary["ran_retrieval"] = True

    # ------------------------------------------------------------
    # Step C. Confounder analysis
    # ------------------------------------------------------------
    if cfg.run_confounder:
        print(">>> [Task1 Pipeline] Step C: confounder analysis")
        matched_pairs_path = _resolve_with_sibling_fallback(
            str(Path(cfg.task1_output_dir) / "data" / "m1_matched_pairs.parquet")
        )
        per_pair_path = _resolve_with_sibling_fallback(
            str(Path(cfg.task1_output_dir) / "analysis" / "modality_gap_per_pair.csv")
        )
        conf_cfg = Task1ConfounderConfig(
            per_pair_path=str(per_pair_path),
            matched_pairs_path=str(matched_pairs_path),
            output_dir=cfg.confounder_output_dir,
            n_perm=cfg.confounder_n_perm,
            min_group_n=cfg.confounder_min_group_n,
            max_groups_for_report=cfg.confounder_max_groups_for_report,
            random_seed=cfg.random_seed,
        )
        run_task1_confounder(conf_cfg)
        summary["ran_confounder"] = True

    # ------------------------------------------------------------
    # Step D. Set-level centroid analysis
    # ------------------------------------------------------------
    if cfg.run_set_level:
        print(">>> [Task1 Pipeline] Step D: set-level centroid analysis")
        candidates_path = _resolve_candidates_path(cfg.task1_output_dir)
        per_pair_path = _resolve_with_sibling_fallback(
            str(Path(cfg.task1_output_dir) / "analysis" / "modality_gap_per_pair.csv")
        )
        set_cfg = Task1SetLevelConfig(
            m1_candidates_path=str(candidates_path),
            processed_dir=cfg.processed_dir,
            lincs_tensor_name=cfg.lincs_tensor_name,
            output_dir=cfg.set_level_output_dir,
            per_pair_path=str(per_pair_path),
            min_n_per_source=cfg.set_level_min_n_per_source,
            n_perm=cfg.set_level_n_perm,
            n_bootstrap=cfg.set_level_n_bootstrap,
            random_seed=cfg.random_seed,
        )
        run_task1_set_level(set_cfg)
        summary["ran_set_level"] = True

    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": summary,
    }
    Path(cfg.pipeline_manifest_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f">>> [Task1 Pipeline] done. Manifest: {cfg.pipeline_manifest_path}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Task1 readable pipeline runner")

    # shared
    parser.add_argument("--task0-unified-meta-path", type=str, default=Task1PipelineConfig.task0_unified_meta_path)
    parser.add_argument("--processed-dir", type=str, default=Task1PipelineConfig.processed_dir)
    parser.add_argument("--lincs-tensor-name", type=str, default=Task1PipelineConfig.lincs_tensor_name)
    parser.add_argument("--task1-output-dir", type=str, default=Task1PipelineConfig.task1_output_dir)

    # steps
    parser.add_argument("--run-groupwise", action=argparse.BooleanOptionalAction, default=Task1PipelineConfig.run_groupwise)
    parser.add_argument("--run-retrieval", action=argparse.BooleanOptionalAction, default=Task1PipelineConfig.run_retrieval)
    parser.add_argument("--run-confounder", action=argparse.BooleanOptionalAction, default=Task1PipelineConfig.run_confounder)
    parser.add_argument("--run-set-level", action=argparse.BooleanOptionalAction, default=Task1PipelineConfig.run_set_level)
    parser.add_argument("--dry-run", action=argparse.BooleanOptionalAction, default=Task1PipelineConfig.dry_run)

    # groupwise
    parser.add_argument("--max-pairs-total", type=int, default=Task1PipelineConfig.max_pairs_total)
    parser.add_argument("--n-perm-groupwise", type=int, default=Task1PipelineConfig.n_perm_groupwise)
    parser.add_argument("--n-bootstrap-groupwise", type=int, default=Task1PipelineConfig.n_bootstrap_groupwise)
    parser.add_argument("--max-groups-per-type", type=int, default=Task1PipelineConfig.max_groups_per_type)
    parser.add_argument("--min-group-n-for-stats", type=int, default=Task1PipelineConfig.min_group_n_for_stats)
    parser.add_argument(
        "--strict-enable",
        action=argparse.BooleanOptionalAction,
        default=Task1PipelineConfig.strict_enable,
        help="Enable strict protocol subset analysis in Task1 groupwise step.",
    )
    parser.add_argument(
        "--strict-chemical-max-dose-logdiff",
        type=float,
        default=Task1PipelineConfig.strict_chemical_max_dose_logdiff,
        help="Chemical strict threshold for abs(log1p dose diff).",
    )
    parser.add_argument(
        "--strict-chemical-max-time-absdiff",
        type=float,
        default=Task1PipelineConfig.strict_chemical_max_time_absdiff,
        help="Chemical strict threshold for absolute time diff (hours).",
    )

    # retrieval
    parser.add_argument("--retrieval-output-dir", type=str, default=Task1PipelineConfig.retrieval_output_dir)
    parser.add_argument("--retrieval-tracks", type=str, default=Task1PipelineConfig.retrieval_tracks)
    parser.add_argument("--retrieval-directions", type=str, default=Task1PipelineConfig.retrieval_directions)
    parser.add_argument("--retrieval-topk", type=str, default=Task1PipelineConfig.retrieval_topk)
    parser.add_argument("--retrieval-batch-size", type=int, default=Task1PipelineConfig.retrieval_batch_size)
    parser.add_argument("--retrieval-n-perm", type=int, default=Task1PipelineConfig.retrieval_n_perm)
    parser.add_argument(
        "--retrieval-min-queries-per-group",
        type=int,
        default=Task1PipelineConfig.retrieval_min_queries_per_group,
    )
    parser.add_argument(
        "--retrieval-run-balanced-eval",
        action=argparse.BooleanOptionalAction,
        default=Task1PipelineConfig.retrieval_run_balanced_eval,
        help="Compute balanced candidate-space retrieval metrics.",
    )
    parser.add_argument(
        "--retrieval-balanced-gallery-size",
        type=int,
        default=Task1PipelineConfig.retrieval_balanced_gallery_size,
    )
    parser.add_argument(
        "--retrieval-balanced-true-per-query",
        type=int,
        default=Task1PipelineConfig.retrieval_balanced_true_per_query,
    )
    parser.add_argument(
        "--retrieval-balanced-n-repeats",
        type=int,
        default=Task1PipelineConfig.retrieval_balanced_n_repeats,
    )

    # confounder
    parser.add_argument("--confounder-output-dir", type=str, default=Task1PipelineConfig.confounder_output_dir)
    parser.add_argument("--confounder-n-perm", type=int, default=Task1PipelineConfig.confounder_n_perm)
    parser.add_argument("--confounder-min-group-n", type=int, default=Task1PipelineConfig.confounder_min_group_n)
    parser.add_argument(
        "--confounder-max-groups-for-report",
        type=int,
        default=Task1PipelineConfig.confounder_max_groups_for_report,
    )

    # set-level
    parser.add_argument("--set-level-output-dir", type=str, default=Task1PipelineConfig.set_level_output_dir)
    parser.add_argument("--set-level-min-n-per-source", type=int, default=Task1PipelineConfig.set_level_min_n_per_source)
    parser.add_argument("--set-level-n-perm", type=int, default=Task1PipelineConfig.set_level_n_perm)
    parser.add_argument("--set-level-n-bootstrap", type=int, default=Task1PipelineConfig.set_level_n_bootstrap)

    parser.add_argument("--random-seed", type=int, default=Task1PipelineConfig.random_seed)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task1PipelineConfig(
        task0_unified_meta_path=args.task0_unified_meta_path,
        processed_dir=args.processed_dir,
        lincs_tensor_name=args.lincs_tensor_name,
        task1_output_dir=args.task1_output_dir,
        run_groupwise=args.run_groupwise,
        run_retrieval=args.run_retrieval,
        run_confounder=args.run_confounder,
        run_set_level=args.run_set_level,
        dry_run=args.dry_run,
        max_pairs_total=args.max_pairs_total,
        n_perm_groupwise=args.n_perm_groupwise,
        n_bootstrap_groupwise=args.n_bootstrap_groupwise,
        max_groups_per_type=args.max_groups_per_type,
        min_group_n_for_stats=args.min_group_n_for_stats,
        strict_enable=args.strict_enable,
        strict_chemical_max_dose_logdiff=args.strict_chemical_max_dose_logdiff,
        strict_chemical_max_time_absdiff=args.strict_chemical_max_time_absdiff,
        retrieval_output_dir=args.retrieval_output_dir,
        retrieval_tracks=args.retrieval_tracks,
        retrieval_directions=args.retrieval_directions,
        retrieval_topk=args.retrieval_topk,
        retrieval_batch_size=args.retrieval_batch_size,
        retrieval_n_perm=args.retrieval_n_perm,
        retrieval_min_queries_per_group=args.retrieval_min_queries_per_group,
        retrieval_run_balanced_eval=args.retrieval_run_balanced_eval,
        retrieval_balanced_gallery_size=args.retrieval_balanced_gallery_size,
        retrieval_balanced_true_per_query=args.retrieval_balanced_true_per_query,
        retrieval_balanced_n_repeats=args.retrieval_balanced_n_repeats,
        confounder_output_dir=args.confounder_output_dir,
        confounder_n_perm=args.confounder_n_perm,
        confounder_min_group_n=args.confounder_min_group_n,
        confounder_max_groups_for_report=args.confounder_max_groups_for_report,
        set_level_output_dir=args.set_level_output_dir,
        set_level_min_n_per_source=args.set_level_min_n_per_source,
        set_level_n_perm=args.set_level_n_perm,
        set_level_n_bootstrap=args.set_level_n_bootstrap,
        random_seed=args.random_seed,
    )
    run_task1_pipeline(cfg)


if __name__ == "__main__":
    main()
