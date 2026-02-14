#!/usr/bin/env python3
"""
Task2 result audit (non-destructive, reproducibility-oriented).

This script audits already-generated Task2 CSV outputs and exports:
  - file manifest
  - consistency checks
  - key numerical summaries
  - markdown report

Default input path points to existing NAS Task2 outputs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


REQUIRED_FILES = [
    "Step1_L1_Instance_Tidy.csv",
    "Step1_L1_Context_Aggregated.csv",
    "Step1_Tests.csv",
    "Step2_Context_Labels.csv",
    "Step2_Diagnostics_SourcePairLevel.csv",
    "Step4_CaseStudy_Tracer.csv",
    "Step5_Enrichment_Cells_RobustHigh.csv",
    "Step5_Enrichment_Targets_ProtocolSensitive.csv",
    "Step5_Enrichment_Targets_RobustHigh.csv",
    "Step5_Enrichment_Targets_RobustLow.csv",
    "Step5_Protocol_Correlations_ProtocolSensitive.csv",
]


def _write_md(path: Path, lines: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_audit(task2_dir: Path, output_dir: Path) -> None:
    analysis_dir = output_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------
    # 1) File manifest
    # -----------------------------------------------------------------
    manifest_rows = []
    missing = []
    for fname in REQUIRED_FILES:
        fpath = task2_dir / fname
        exists = fpath.exists()
        size_mb = (fpath.stat().st_size / (1024 * 1024)) if exists else np.nan
        manifest_rows.append(
            {
                "file_name": fname,
                "exists": bool(exists),
                "size_mb": float(size_mb) if np.isfinite(size_mb) else np.nan,
                "path": str(fpath),
            }
        )
        if not exists:
            missing.append(fname)
    file_manifest = pd.DataFrame(manifest_rows)
    file_manifest.to_csv(analysis_dir / "task2_file_manifest.csv", index=False)
    if missing:
        raise FileNotFoundError(f"Missing required Task2 files: {missing}")

    # -----------------------------------------------------------------
    # 2) Load core tables
    # -----------------------------------------------------------------
    step1_inst = pd.read_csv(task2_dir / "Step1_L1_Instance_Tidy.csv")
    step1_ctx = pd.read_csv(task2_dir / "Step1_L1_Context_Aggregated.csv")
    step1_test = pd.read_csv(task2_dir / "Step1_Tests.csv")
    step2_label = pd.read_csv(task2_dir / "Step2_Context_Labels.csv")
    step2_diag = pd.read_csv(task2_dir / "Step2_Diagnostics_SourcePairLevel.csv")
    step4_case = pd.read_csv(task2_dir / "Step4_CaseStudy_Tracer.csv")
    enr_cells = pd.read_csv(task2_dir / "Step5_Enrichment_Cells_RobustHigh.csv")
    enr_ps = pd.read_csv(task2_dir / "Step5_Enrichment_Targets_ProtocolSensitive.csv")
    enr_rh = pd.read_csv(task2_dir / "Step5_Enrichment_Targets_RobustHigh.csv")
    enr_rl = pd.read_csv(task2_dir / "Step5_Enrichment_Targets_RobustLow.csv")
    prot_corr = pd.read_csv(task2_dir / "Step5_Protocol_Correlations_ProtocolSensitive.csv")

    # -----------------------------------------------------------------
    # 3) Consistency checks
    # -----------------------------------------------------------------
    checks = []

    checks.append(
        {
            "check_name": "step1_instance_nonempty",
            "passed": bool(len(step1_inst) > 0),
            "detail": f"n_rows={len(step1_inst)}",
        }
    )
    checks.append(
        {
            "check_name": "step1_context_nonempty",
            "passed": bool(len(step1_ctx) > 0),
            "detail": f"n_rows={len(step1_ctx)}",
        }
    )

    # Re-aggregate from step1 instance and compare with step1 context table size.
    # Step1 context table applies MIN_QUERY_SAMPLES_PER_CONTEXT=3 in legacy pipeline.
    need_ctx_cols = ["Track", "View", "DomainType", "Cell", "Target", "Source_Chem", "Source_Gene"]
    has_cols = all(c in step1_inst.columns for c in need_ctx_cols)
    if has_cols:
        rebuilt = step1_inst.groupby(need_ctx_cols, dropna=False).size().reset_index(name="n_instances")
        rebuild_ctx_n_raw = int(len(rebuilt))
        rebuild_ctx_n_min3 = int((rebuilt["n_instances"] >= 3).sum())
        checks.append(
            {
                "check_name": "step1_context_key_count_match_min3",
                "passed": bool(rebuild_ctx_n_min3 == len(step1_ctx)),
                "detail": (
                    f"from_instance_raw={rebuild_ctx_n_raw}, "
                    f"from_instance_min3={rebuild_ctx_n_min3}, "
                    f"from_context_csv={len(step1_ctx)}"
                ),
            }
        )
    else:
        checks.append(
            {
                "check_name": "step1_context_key_count_match_min3",
                "passed": False,
                "detail": "missing required key columns in Step1_L1_Instance_Tidy.csv",
            }
        )

    # Value-range checks
    for col in ["MRR", "Success_Score"]:
        if col in step1_inst.columns:
            ok = bool(step1_inst[col].dropna().between(0, 1).all())
            checks.append(
                {
                    "check_name": f"step1_{col}_in_0_1",
                    "passed": ok,
                    "detail": f"min={step1_inst[col].min()}, max={step1_inst[col].max()}",
                }
            )

    if "P_Value" in step1_test.columns:
        ok = bool(step1_test["P_Value"].dropna().between(0, 1).all())
        checks.append(
            {
                "check_name": "step1_tests_pvalue_in_0_1",
                "passed": ok,
                "detail": f"min={step1_test['P_Value'].min()}, max={step1_test['P_Value'].max()}",
            }
        )
    if "FDR_BH" in step1_test.columns:
        ok = bool(step1_test["FDR_BH"].dropna().between(0, 1).all())
        checks.append(
            {
                "check_name": "step1_tests_fdr_in_0_1",
                "passed": ok,
                "detail": f"min={step1_test['FDR_BH'].min()}, max={step1_test['FDR_BH'].max()}",
            }
        )

    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(analysis_dir / "task2_consistency_checks.csv", index=False)

    # -----------------------------------------------------------------
    # 4) Key summaries
    # -----------------------------------------------------------------
    summary_rows = []
    summary_rows.append({"metric": "step1_instance_rows", "value": float(len(step1_inst))})
    summary_rows.append({"metric": "step1_context_rows", "value": float(len(step1_ctx))})
    summary_rows.append({"metric": "step2_context_labels_rows", "value": float(len(step2_label))})
    summary_rows.append({"metric": "step2_sourcepair_diag_rows", "value": float(len(step2_diag))})
    summary_rows.append({"metric": "step4_case_rows", "value": float(len(step4_case))})
    summary_rows.append({"metric": "step5_protocol_corr_rows", "value": float(len(prot_corr))})

    if "DomainType" in step1_inst.columns:
        dom = step1_inst["DomainType"].value_counts()
        for k, v in dom.items():
            summary_rows.append({"metric": f"step1_domain_rows::{k}", "value": float(v)})

    if "Scenario" in step1_inst.columns:
        scn = step1_inst["Scenario"].value_counts()
        for k, v in scn.items():
            summary_rows.append({"metric": f"step1_scenario_rows::{k}", "value": float(v)})

    key_summary = pd.DataFrame(summary_rows)
    key_summary.to_csv(analysis_dir / "task2_key_summary.csv", index=False)

    # Primary slice summary (A_L1Base, Chem->Gene, Gene track)
    primary = step1_inst.copy()
    for col, val in [("Scenario", "A_L1Base"), ("Direction", "Chem->Gene"), ("Track", "Gene")]:
        if col in primary.columns:
            primary = primary[primary[col] == val]

    prim_grp_cols = [c for c in ["View", "DomainType"] if c in primary.columns]
    if prim_grp_cols:
        primary_summary = (
            primary.groupby(prim_grp_cols, dropna=False)
            .agg(
                n_instances=("MRR", "count"),
                mean_mrr=("MRR", "mean"),
                mean_success=("Success_Score", "mean"),
                median_rank=("True_Rank", "median"),
            )
            .reset_index()
            .sort_values(["DomainType", "View"])
        )
    else:
        primary_summary = pd.DataFrame()
    primary_summary.to_csv(analysis_dir / "task2_primary_slice_summary.csv", index=False)

    # Step2 labels summary
    label_col = "Performance_Class" if "Performance_Class" in step2_label.columns else None
    if label_col is not None:
        label_counts = (
            step2_label[label_col].value_counts(dropna=False)
            .rename_axis("Performance_Class")
            .reset_index(name="n_contexts")
        )
    else:
        label_counts = pd.DataFrame()
    label_counts.to_csv(analysis_dir / "task2_label_counts.csv", index=False)

    # Enrichment top hits
    enrich_frames = []
    for name, df in [
        ("Cells_RobustHigh", enr_cells),
        ("Targets_ProtocolSensitive", enr_ps),
        ("Targets_RobustHigh", enr_rh),
        ("Targets_RobustLow", enr_rl),
    ]:
        tmp = df.copy()
        tmp["table_name"] = name
        if "FDR_BH" in tmp.columns:
            tmp = tmp.sort_values(["FDR_BH", "P_Value"], ascending=[True, True])
        enrich_frames.append(tmp.head(20))
    enrich_top = pd.concat(enrich_frames, ignore_index=True)
    enrich_top.to_csv(analysis_dir / "task2_enrichment_top_hits.csv", index=False)

    # -----------------------------------------------------------------
    # 5) Markdown report
    # -----------------------------------------------------------------
    lines: List[str] = []
    lines.append("# Task2 Result Audit")
    lines.append("")
    lines.append("## File Completeness")
    lines.append(f"- Required files: {len(REQUIRED_FILES)}")
    lines.append(f"- All files present: {str((~file_manifest['exists'].eq(False)).all())}")
    lines.append("")
    lines.append("## Core Volumes")
    lines.append(f"- Step1 instance rows: {len(step1_inst)}")
    lines.append(f"- Step1 context rows: {len(step1_ctx)}")
    lines.append(f"- Step2 labeled contexts: {len(step2_label)}")
    lines.append("")
    lines.append("## Consistency Checks")
    for _, r in checks_df.iterrows():
        lines.append(f"- {r['check_name']}: {bool(r['passed'])} ({r['detail']})")
    lines.append("")

    if not primary_summary.empty:
        lines.append("## Primary Slice (A_L1Base, Chem->Gene, Gene)")
        for _, r in primary_summary.iterrows():
            lines.append(
                "- "
                + ", ".join([f"{c}={r[c]}" for c in prim_grp_cols])
                + f", n={int(r['n_instances'])}, mean_success={r['mean_success']:.4f}, mean_mrr={r['mean_mrr']:.4f}"
            )

    _write_md(analysis_dir / "task2_audit_report.md", lines)

    run_manifest = {
        "task2_dir": str(task2_dir),
        "summary": {
            "n_required_files": len(REQUIRED_FILES),
            "n_instance_rows": int(len(step1_inst)),
            "n_context_rows": int(len(step1_ctx)),
            "n_labels": int(len(step2_label)),
            "n_checks": int(len(checks_df)),
            "n_checks_passed": int(checks_df["passed"].sum()),
        },
        "outputs": {"analysis_dir": str(analysis_dir)},
    }
    (output_dir / "run_manifest_task2_audit.json").write_text(
        json.dumps(run_manifest, indent=2), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser("Task2 result audit")
    parser.add_argument(
        "--task2-dir",
        type=str,
        default="/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task2_Unified",
        help="Directory that already contains Step1/2/4/5 Task2 CSV outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./outputs/task2_audit",
    )
    args = parser.parse_args()

    run_audit(task2_dir=Path(args.task2_dir), output_dir=Path(args.output_dir))


if __name__ == "__main__":
    main()
