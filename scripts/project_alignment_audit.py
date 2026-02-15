#!/usr/bin/env python3
"""
Project-level alignment audit for M2M-Bench.

Checks whether core manuscript-facing claims remain aligned with generated outputs.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class ProjectAlignmentAuditConfig:
    output_dir: str = "./outputs/project_alignment_audit"

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_project_alignment_audit.json")


def _pass_fail(cond: bool) -> str:
    return "PASS" if bool(cond) else "FAIL"


def run_project_alignment_audit(cfg: ProjectAlignmentAuditConfig) -> None:
    t0 = time.time()
    out = Path(cfg.analysis_dir)
    out.mkdir(parents=True, exist_ok=True)

    checks: list[dict] = []

    # Task1: pair counts
    per_pair = pd.read_csv("./outputs/task1/analysis/modality_gap_per_pair.csv")
    checks.append(
        {
            "check_id": "task1_pair_count",
            "status": _pass_fail(len(per_pair) == 670),
            "value": int(len(per_pair)),
            "expected": "670",
            "detail": "Task1 matched pair count",
        }
    )

    # Task1: strict subset composition
    strict = pd.read_csv("./outputs/task1_reviewer_fixes/analysis/strict_subset_composition.csv")
    strict_true = strict[strict["is_strict_protocol_pair"] == True]  # noqa: E712
    strict_true_mods = sorted(strict_true["modality"].astype(str).tolist())
    checks.append(
        {
            "check_id": "task1_strict_composition",
            "status": _pass_fail(strict_true_mods == ["Genetic"]),
            "value": ",".join(strict_true_mods),
            "expected": "Genetic-only",
            "detail": "Strict protocol subset modality composition",
        }
    )

    # Task1: protocol partial trend sign
    partial = pd.read_csv("./outputs/task1_reviewer_fixes/analysis/protocol_continuous_partial_spearman.csv")
    trow = partial[
        (partial["x_var"] == "time_absdiff")
        & (partial["y_var"] == "cosine_gene")
        & (partial["control_spec"] == "cell_target_plus_other_mismatch_plus_match_score")
    ]
    if trow.empty:
        trow = partial[
            (partial["x_var"] == "time_absdiff")
            & (partial["y_var"] == "cosine_gene")
            & (partial["control_spec"] == "cell_target_plus_other_mismatch")
        ]
    rho = float(trow["partial_spearman_rho"].iloc[0]) if not trow.empty else np.nan
    checks.append(
        {
            "check_id": "task1_time_partial_negative",
            "status": _pass_fail(np.isfinite(rho) and rho < 0),
            "value": rho,
            "expected": "<0",
            "detail": "Deconfounded time mismatch effect should stay negative",
        }
    )

    # Task1: retrieval sensitivity presence
    sens_path = Path("./outputs/task1_reviewer_fixes/retrieval_sensitivity/analysis/retrieval_sensitivity_summary.csv")
    checks.append(
        {
            "check_id": "task1_retrieval_sensitivity_exists",
            "status": _pass_fail(sens_path.exists()),
            "value": int(sens_path.exists()),
            "expected": "1",
            "detail": "Balanced retrieval sensitivity output should exist",
        }
    )

    # Task2: labeled context count and class sum
    step2 = pd.read_csv("./outputs/task2_nodomain/analysis/Step2_Context_Labels_NoDomain.csv")
    checks.append(
        {
            "check_id": "task2_labeled_context_count",
            "status": _pass_fail(len(step2) == 3061),
            "value": int(len(step2)),
            "expected": "3061",
            "detail": "Task2 labeled context row count",
        }
    )

    # Task3: mean_scaled_best definition exists
    score_def = Path("./outputs/task3_meta/analysis/task3_mean_scaled_best_definition.csv")
    checks.append(
        {
            "check_id": "task3_score_definition_exists",
            "status": _pass_fail(score_def.exists()),
            "value": int(score_def.exists()),
            "expected": "1",
            "detail": "Task3 mean_scaled_best definition export",
        }
    )

    # Repro pack coverage
    inv = pd.read_csv("./outputs/reproducibility/analysis/repro_manifest_inventory.csv")
    n_found = int(inv["exists"].sum())
    n_total = int(len(inv))
    checks.append(
        {
            "check_id": "repro_manifest_coverage",
            "status": _pass_fail(n_found == n_total),
            "value": f"{n_found}/{n_total}",
            "expected": "all found",
            "detail": "Repro manifest coverage",
        }
    )

    df = pd.DataFrame(checks)
    df.to_csv(out / "project_alignment_checks.csv", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "n_checks", "value": int(len(df))},
            {"metric": "n_pass", "value": int((df["status"] == "PASS").sum())},
            {"metric": "n_fail", "value": int((df["status"] == "FAIL").sum())},
        ]
    )
    summary.to_csv(out / "project_alignment_summary.csv", index=False)

    report_lines = [
        "# Project Alignment Audit",
        "",
        f"- Generated: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}",
        f"- Checks: {int(len(df))}, PASS: {int((df['status'] == 'PASS').sum())}, FAIL: {int((df['status'] == 'FAIL').sum())}",
        "",
        "## Failed checks",
    ]
    fail_df = df[df["status"] == "FAIL"]
    if fail_df.empty:
        report_lines.append("- None")
    else:
        for _, r in fail_df.iterrows():
            report_lines.append(
                f"- `{r['check_id']}`: value={r['value']} expected={r['expected']} ({r['detail']})"
            )
    (out / "project_alignment_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_checks": int(len(df)),
            "n_pass": int((df["status"] == "PASS").sum()),
            "n_fail": int((df["status"] == "FAIL").sum()),
        },
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[Project Alignment Audit] done. Outputs: {cfg.analysis_dir}")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser("Project alignment audit")
    p.add_argument("--output-dir", type=str, default=ProjectAlignmentAuditConfig.output_dir)
    return p


def main() -> None:
    args = _build_parser().parse_args()
    cfg = ProjectAlignmentAuditConfig(output_dir=args.output_dir)
    run_project_alignment_audit(cfg)


if __name__ == "__main__":
    main()
