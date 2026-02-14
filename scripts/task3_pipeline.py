#!/usr/bin/env python3
"""
Task3 readable pipeline runner.

This script wraps legacy Task3 scripts (9/10/11) into one explicit workflow:
  A) prepare deltas
  B) run analysis
  C) prepare R-ready tables
  D) run result audit

Each step is independently switchable for iterative reruns.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from task3_results_audit import Task3AuditConfig, run_task3_audit


@dataclass
class Task3PipelineConfig:
    eval_dir: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562"
    model_base_path: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results"
    deltas_out_dir: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Deltas_K562_CellLevel"
    analysis_out_dir: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Analysis"
    unified_out_dir: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified"
    local_audit_output_dir: str = "./outputs/task3_audit"

    run_prepare_deltas: bool = False
    run_analysis: bool = False
    run_collect_unified: bool = False
    run_audit: bool = True
    dry_run: bool = False

    # Step-9 options
    seed: int = 42
    k_controls: int = 1
    batch_size: int = 1024
    pb_scale: float = 10000.0
    use_zscore: bool = False
    z_eps: float = 1e-6
    z_clip: float = 10.0
    pathway_csv: str = ""
    pathway_name: str = "HALLMARK"

    # Step-10 options
    views: str = "Standard,Systema"
    tracks: str = "ALL"
    do_retrieval: bool = True
    do_pairwise: bool = True
    nan_policy: str = "drop"
    edist_boot: int = 200
    edist_max_n: int = 300

    # Step-11 options
    panels: str = "All,Cleanest,Family,Promiscuous"
    baseline_track: str = "Gene"
    fig2_direction: str = "CRISPR->Drug"
    fig2_eps_ed: float = 1e-8
    fig2_boot_median: int = 500
    fig2_boot_max_n: int = 2000
    fig2_seed: int = 42

    @property
    def pipeline_manifest_path(self) -> str:
        return str(Path(self.local_audit_output_dir).parent / "task3_pipeline_manifest.json")


def _run_cmd(cmd: Sequence[str], cwd: str | None = None, dry_run: bool = False) -> None:
    printable = " ".join(shlex.quote(str(x)) for x in cmd)
    if cwd:
        printable = f"(cd {shlex.quote(cwd)} && {printable})"
    print(f">>> {printable}")
    if dry_run:
        return
    subprocess.run([str(x) for x in cmd], cwd=cwd, check=True)


def run_task3_pipeline(cfg: Task3PipelineConfig) -> None:
    t0 = time.time()
    summary = {
        "ran_prepare_deltas": False,
        "ran_analysis": False,
        "ran_collect_unified": False,
        "ran_audit": False,
    }

    if cfg.dry_run:
        print("[DRY RUN] Task3 pipeline configuration:")
        print(json.dumps(asdict(cfg), indent=2))

    repo_root = Path(__file__).resolve().parents[1]
    legacy_dir = str(repo_root / "Chem2Gene")

    # ------------------------------------------------------------
    # Step A. Prepare deltas (Script 9)
    # ------------------------------------------------------------
    if cfg.run_prepare_deltas:
        cmd = [
            "python",
            "9_Task3_PrepareAnalysis.py",
            "--eval_dir",
            cfg.eval_dir,
            "--base_path",
            cfg.model_base_path,
            "--out_dir",
            cfg.deltas_out_dir,
            "--seed",
            str(cfg.seed),
            "--group_level",
            "cell",
            "--k_controls",
            str(cfg.k_controls),
            "--batch_size",
            str(cfg.batch_size),
            "--pb_scale",
            str(cfg.pb_scale),
            "--z_eps",
            str(cfg.z_eps),
            "--z_clip",
            str(cfg.z_clip),
            "--pathway_name",
            str(cfg.pathway_name),
        ]
        if cfg.use_zscore:
            cmd.append("--use_zscore")
        if str(cfg.pathway_csv).strip():
            cmd.extend(["--pathway_csv", str(cfg.pathway_csv).strip()])
        _run_cmd(cmd, cwd=legacy_dir, dry_run=cfg.dry_run)
        summary["ran_prepare_deltas"] = True

    # ------------------------------------------------------------
    # Step B. Analyze deltas (Script 10)
    # ------------------------------------------------------------
    if cfg.run_analysis:
        cmd = [
            "python",
            "10_Task3_AnalyzeOutputs.py",
            "--task3_out_dir",
            cfg.deltas_out_dir,
            "--out_dir",
            cfg.analysis_out_dir,
            "--views",
            cfg.views,
            "--tracks",
            cfg.tracks,
            "--nan_policy",
            cfg.nan_policy,
            "--edist_boot",
            str(cfg.edist_boot),
            "--edist_max_n",
            str(cfg.edist_max_n),
            "--seed",
            str(cfg.seed),
        ]
        if cfg.do_retrieval:
            cmd.append("--do_retrieval")
        if cfg.do_pairwise:
            cmd.append("--do_pairwise")
        _run_cmd(cmd, cwd=legacy_dir, dry_run=cfg.dry_run)
        summary["ran_analysis"] = True

    # ------------------------------------------------------------
    # Step C. Build unified visualization tables (Script 11)
    # ------------------------------------------------------------
    if cfg.run_collect_unified:
        cmd = [
            "python",
            "11_Task3_OutputCollect.py",
            "--eval_dir",
            cfg.eval_dir,
            "--analysis_dir",
            cfg.analysis_out_dir,
            "--out_dir",
            cfg.unified_out_dir,
            "--panels",
            cfg.panels,
            "--baseline_track",
            cfg.baseline_track,
            "--fig2_direction",
            cfg.fig2_direction,
            "--fig2_eps_ed",
            str(cfg.fig2_eps_ed),
            "--fig2_boot_median",
            str(cfg.fig2_boot_median),
            "--fig2_boot_max_n",
            str(cfg.fig2_boot_max_n),
            "--fig2_seed",
            str(cfg.fig2_seed),
        ]
        _run_cmd(cmd, cwd=legacy_dir, dry_run=cfg.dry_run)
        summary["ran_collect_unified"] = True

    # ------------------------------------------------------------
    # Step D. Audit finished outputs
    # ------------------------------------------------------------
    if cfg.run_audit:
        audit_cfg = Task3AuditConfig(
            deltas_dir=cfg.deltas_out_dir,
            analysis_dir=cfg.analysis_out_dir,
            unified_dir=cfg.unified_out_dir,
            output_dir=cfg.local_audit_output_dir,
        )
        if cfg.dry_run:
            print(">>> [DRY RUN] skip run_task3_audit()")
        else:
            run_task3_audit(audit_cfg)
        summary["ran_audit"] = True

    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": summary,
    }
    Path(cfg.pipeline_manifest_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f">>> [Task3 Pipeline] done. Manifest: {cfg.pipeline_manifest_path}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task3 readable pipeline runner")

    parser.add_argument("--eval-dir", type=str, default=Task3PipelineConfig.eval_dir)
    parser.add_argument("--model-base-path", type=str, default=Task3PipelineConfig.model_base_path)
    parser.add_argument("--deltas-out-dir", type=str, default=Task3PipelineConfig.deltas_out_dir)
    parser.add_argument("--analysis-out-dir", type=str, default=Task3PipelineConfig.analysis_out_dir)
    parser.add_argument("--unified-out-dir", type=str, default=Task3PipelineConfig.unified_out_dir)
    parser.add_argument("--local-audit-output-dir", type=str, default=Task3PipelineConfig.local_audit_output_dir)

    parser.add_argument(
        "--run-prepare-deltas",
        action=argparse.BooleanOptionalAction,
        default=Task3PipelineConfig.run_prepare_deltas,
    )
    parser.add_argument(
        "--run-analysis",
        action=argparse.BooleanOptionalAction,
        default=Task3PipelineConfig.run_analysis,
    )
    parser.add_argument(
        "--run-collect-unified",
        action=argparse.BooleanOptionalAction,
        default=Task3PipelineConfig.run_collect_unified,
    )
    parser.add_argument("--run-audit", action=argparse.BooleanOptionalAction, default=Task3PipelineConfig.run_audit)
    parser.add_argument("--dry-run", action=argparse.BooleanOptionalAction, default=Task3PipelineConfig.dry_run)

    parser.add_argument("--seed", type=int, default=Task3PipelineConfig.seed)
    parser.add_argument("--k-controls", type=int, default=Task3PipelineConfig.k_controls)
    parser.add_argument("--batch-size", type=int, default=Task3PipelineConfig.batch_size)
    parser.add_argument("--pb-scale", type=float, default=Task3PipelineConfig.pb_scale)
    parser.add_argument("--use-zscore", action=argparse.BooleanOptionalAction, default=Task3PipelineConfig.use_zscore)
    parser.add_argument("--z-eps", type=float, default=Task3PipelineConfig.z_eps)
    parser.add_argument("--z-clip", type=float, default=Task3PipelineConfig.z_clip)
    parser.add_argument("--pathway-csv", type=str, default=Task3PipelineConfig.pathway_csv)
    parser.add_argument("--pathway-name", type=str, default=Task3PipelineConfig.pathway_name)

    parser.add_argument("--views", type=str, default=Task3PipelineConfig.views)
    parser.add_argument("--tracks", type=str, default=Task3PipelineConfig.tracks)
    parser.add_argument("--do-retrieval", action=argparse.BooleanOptionalAction, default=Task3PipelineConfig.do_retrieval)
    parser.add_argument("--do-pairwise", action=argparse.BooleanOptionalAction, default=Task3PipelineConfig.do_pairwise)
    parser.add_argument("--nan-policy", type=str, choices=["drop", "zero"], default=Task3PipelineConfig.nan_policy)
    parser.add_argument("--edist-boot", type=int, default=Task3PipelineConfig.edist_boot)
    parser.add_argument("--edist-max-n", type=int, default=Task3PipelineConfig.edist_max_n)

    parser.add_argument("--panels", type=str, default=Task3PipelineConfig.panels)
    parser.add_argument("--baseline-track", type=str, default=Task3PipelineConfig.baseline_track)
    parser.add_argument("--fig2-direction", type=str, default=Task3PipelineConfig.fig2_direction)
    parser.add_argument("--fig2-eps-ed", type=float, default=Task3PipelineConfig.fig2_eps_ed)
    parser.add_argument("--fig2-boot-median", type=int, default=Task3PipelineConfig.fig2_boot_median)
    parser.add_argument("--fig2-boot-max-n", type=int, default=Task3PipelineConfig.fig2_boot_max_n)
    parser.add_argument("--fig2-seed", type=int, default=Task3PipelineConfig.fig2_seed)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task3PipelineConfig(
        eval_dir=args.eval_dir,
        model_base_path=args.model_base_path,
        deltas_out_dir=args.deltas_out_dir,
        analysis_out_dir=args.analysis_out_dir,
        unified_out_dir=args.unified_out_dir,
        local_audit_output_dir=args.local_audit_output_dir,
        run_prepare_deltas=args.run_prepare_deltas,
        run_analysis=args.run_analysis,
        run_collect_unified=args.run_collect_unified,
        run_audit=args.run_audit,
        dry_run=args.dry_run,
        seed=args.seed,
        k_controls=args.k_controls,
        batch_size=args.batch_size,
        pb_scale=args.pb_scale,
        use_zscore=args.use_zscore,
        z_eps=args.z_eps,
        z_clip=args.z_clip,
        pathway_csv=args.pathway_csv,
        pathway_name=args.pathway_name,
        views=args.views,
        tracks=args.tracks,
        do_retrieval=args.do_retrieval,
        do_pairwise=args.do_pairwise,
        nan_policy=args.nan_policy,
        edist_boot=args.edist_boot,
        edist_max_n=args.edist_max_n,
        panels=args.panels,
        baseline_track=args.baseline_track,
        fig2_direction=args.fig2_direction,
        fig2_eps_ed=args.fig2_eps_ed,
        fig2_boot_median=args.fig2_boot_median,
        fig2_boot_max_n=args.fig2_boot_max_n,
        fig2_seed=args.fig2_seed,
    )
    run_task3_pipeline(cfg)


if __name__ == "__main__":
    main()
