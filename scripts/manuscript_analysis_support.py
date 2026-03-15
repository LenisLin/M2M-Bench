#!/usr/bin/env python3
# SCRIPT_HEADER_CONTRACT
# Script: scripts/manuscript_analysis_support.py
# Purpose: Provide a registry and local-input checks for planned manuscript-
#   support analyses without recomputing benchmark-stage metrics.
# Inputs:
#   - existing audited Task1/Task2 outputs and manuscript planning files
# Outputs:
#   - stdout only
# Side Effects:
#   - none
# Failure Modes:
#   - missing required local inputs when --check-inputs is used
# Last Updated: 2026-03-13

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
NAS_RUNS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs")
MANUSCRIPT_PHASE1_ROOT = NAS_RUNS_ROOT / "manuscript_phase1"


@dataclass(frozen=True)
class AnalysisSpec:
    analysis_id: str
    phase: str
    analysis_name: str
    analysis_class: str
    reviewer_attack_vector: str
    local_inputs: tuple[str, ...]
    prerequisite_outputs: tuple[str, ...]
    planned_outputs: tuple[str, ...]
    local_only: bool = False
    conditional: bool = False
    last_resort_upstream: tuple[str, ...] = field(default_factory=tuple)


ANALYSIS_SPECS: tuple[AnalysisSpec, ...] = (
    AnalysisSpec(
        analysis_id="A3",
        phase="Phase 1",
        analysis_name="C2G vs G2C direction-specific robustness audit",
        analysis_class="reviewer_defense",
        reviewer_attack_vector="asymmetry_artifact_challenge",
        local_inputs=(
            "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_chance_identity_check.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv",
        ),
        prerequisite_outputs=(),
        planned_outputs=(str(MANUSCRIPT_PHASE1_ROOT / "a3_direction_robustness/task2_direction_robustness_audit.csv"),),
    ),
    AnalysisSpec(
        analysis_id="B6",
        phase="Phase 1",
        analysis_name="Dose/time covariate-overlap audit",
        analysis_class="high_risk_gate",
        reviewer_attack_vector="dose_time_confounding_challenge",
        local_inputs=(
            "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
            "data/task2_snapshot_v2/lincs/derived/delta_meta.csv",
            "data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv",
        ),
        prerequisite_outputs=(),
        planned_outputs=(str(MANUSCRIPT_PHASE1_ROOT / "b6_dose_time_overlap/task2_dose_time_overlap_audit.csv"),),
    ),
    AnalysisSpec(
        analysis_id="A1",
        phase="Phase 2",
        analysis_name="Task1 internal vs cross common-scope comparison",
        analysis_class="main_bridge",
        reviewer_attack_vector="virtual_cell_scope_challenge",
        local_inputs=(
            "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_per_query.parquet",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_leaderboard_long.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_group_cross.parquet",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_per_query.parquet",
        ),
        prerequisite_outputs=(),
        planned_outputs=(str(MANUSCRIPT_PHASE1_ROOT / "a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv"),),
        last_resort_upstream=("/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_group_internal.parquet",),
    ),
    AnalysisSpec(
        analysis_id="A2",
        phase="Phase 2",
        analysis_name="Task1 vs Task2 shared-slice bridge comparison",
        analysis_class="main_bridge",
        reviewer_attack_vector="virtual_cell_scope_challenge",
        local_inputs=(
            "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv",
        ),
        prerequisite_outputs=(),
        planned_outputs=(str(MANUSCRIPT_PHASE1_ROOT / "a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv"),),
        last_resort_upstream=("/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_group_internal.parquet",),
    ),
    AnalysisSpec(
        analysis_id="A4",
        phase="Phase 3",
        analysis_name="Paired/statistical framework for core claims",
        analysis_class="statistical_packaging",
        reviewer_attack_vector="denominator_fragility_challenge",
        local_inputs=(
            "/mnt/NAS_21T/ProjectData/M2M/runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
        ),
        prerequisite_outputs=(
            str(MANUSCRIPT_PHASE1_ROOT / "a3_direction_robustness/task2_direction_robustness_audit.csv"),
            str(MANUSCRIPT_PHASE1_ROOT / "a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv"),
            str(MANUSCRIPT_PHASE1_ROOT / "a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv"),
        ),
        planned_outputs=("docs/manuscript_core_paired_statistics.csv",),
        last_resort_upstream=("/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_group_internal.parquet",),
    ),
    AnalysisSpec(
        analysis_id="B5",
        phase="Phase 4",
        analysis_name="Target-complexity analysis",
        analysis_class="explanation_layer",
        reviewer_attack_vector="denominator_fragility_challenge",
        local_inputs=(
            "docs/task2_c2g_query_n_targets_sensitivity.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
            "runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv",
        ),
        prerequisite_outputs=(),
        planned_outputs=("docs/task2_n_targets_pairwise_stats.csv",),
    ),
    AnalysisSpec(
        analysis_id="B10",
        phase="Phase 4",
        analysis_name="Dataset / cell-line composition analysis",
        analysis_class="explanation_layer",
        reviewer_attack_vector="survivor_bias_enrichment_challenge",
        local_inputs=(
            "data/task2_snapshot_v2/task2_pairs_coverage.csv",
            "data/task2_snapshot_v2/representation_availability_registry.csv",
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_attrition.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv",
            "runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv",
        ),
        prerequisite_outputs=(),
        planned_outputs=("docs/task2_dataset_cellline_composition_analysis.csv",),
    ),
    AnalysisSpec(
        analysis_id="B9",
        phase="Phase 4",
        analysis_name="Enrichment for high-/low-concordance target sets",
        analysis_class="explanation_layer",
        reviewer_attack_vector="survivor_bias_enrichment_challenge",
        local_inputs=(
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv",
        ),
        prerequisite_outputs=(
            "docs/task2_dataset_cellline_composition_analysis.csv",
            "docs/task2_target_pool_background_registry.csv",
        ),
        planned_outputs=(
            "docs/task2_target_concordance_rankings.csv",
            "docs/task2_high_low_concordance_enrichment.csv",
        ),
    ),
    AnalysisSpec(
        analysis_id="B7",
        phase="Phase 4",
        analysis_name="Matched-subset / controlled dose-time analysis",
        analysis_class="explanation_layer",
        reviewer_attack_vector="dose_time_confounding_challenge",
        local_inputs=(
            "/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
            "/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv",
        ),
        prerequisite_outputs=(str(MANUSCRIPT_PHASE1_ROOT / "b6_dose_time_overlap/task2_dose_time_overlap_audit.csv"),),
        planned_outputs=("docs/task2_dose_time_matched_analysis.csv",),
        conditional=True,
    ),
    AnalysisSpec(
        analysis_id="B8",
        phase="Phase 4",
        analysis_name="Specificity local sensitivity",
        analysis_class="reviewer_defense",
        reviewer_attack_vector="denominator_fragility_challenge",
        local_inputs=(
            "docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
            "data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv",
        ),
        prerequisite_outputs=(),
        planned_outputs=("docs/task2_specificity_local_stats.csv",),
        local_only=True,
    ),
    AnalysisSpec(
        analysis_id="C11",
        phase="Phase 5",
        analysis_name="Local FM comparison within supported scPerturb K562 scope",
        analysis_class="representation_layer",
        reviewer_attack_vector="local_fm_outlier_challenge",
        local_inputs=(
            "runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv",
            "runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv",
            "docs/manuscript_panel_support_ledger.csv",
            "data/task2_snapshot_v2/representation_availability_registry.csv",
        ),
        prerequisite_outputs=(),
        planned_outputs=("docs/task2_fm_k562_local_comparison_summary.csv",),
        local_only=True,
    ),
    AnalysisSpec(
        analysis_id="C12",
        phase="Phase 5",
        analysis_name="Per-target FM robustness / paired tests",
        analysis_class="representation_layer",
        reviewer_attack_vector="local_fm_outlier_challenge",
        local_inputs=(
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
        ),
        prerequisite_outputs=("docs/task2_fm_k562_local_comparison_summary.csv",),
        planned_outputs=(
            "docs/task2_fm_k562_per_target_robustness.csv",
            "docs/task2_fm_k562_paired_tests.csv",
        ),
        local_only=True,
    ),
    AnalysisSpec(
        analysis_id="D15",
        phase="Phase 6",
        analysis_name="Modality-good but mechanism-bad cases",
        analysis_class="case_layer",
        reviewer_attack_vector="case_selection_bias_challenge",
        local_inputs=(
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
        ),
        prerequisite_outputs=(
            str(MANUSCRIPT_PHASE1_ROOT / "a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv"),
            str(MANUSCRIPT_PHASE1_ROOT / "a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv"),
            "docs/manuscript_core_paired_statistics.csv",
        ),
        planned_outputs=("docs/manuscript_case_shortlist.csv",),
    ),
    AnalysisSpec(
        analysis_id="D13",
        phase="Phase 6",
        analysis_name="Representative aligned cases",
        analysis_class="case_layer",
        reviewer_attack_vector="case_selection_bias_challenge",
        local_inputs=(
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
        ),
        prerequisite_outputs=(
            str(MANUSCRIPT_PHASE1_ROOT / "a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv"),
            str(MANUSCRIPT_PHASE1_ROOT / "a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv"),
            "docs/manuscript_core_paired_statistics.csv",
        ),
        planned_outputs=("docs/manuscript_case_shortlist.csv",),
    ),
    AnalysisSpec(
        analysis_id="D14",
        phase="Phase 6",
        analysis_name="Representative mechanism-mismatch cases",
        analysis_class="case_layer",
        reviewer_attack_vector="case_selection_bias_challenge",
        local_inputs=(
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
        ),
        prerequisite_outputs=(
            str(MANUSCRIPT_PHASE1_ROOT / "a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv"),
            str(MANUSCRIPT_PHASE1_ROOT / "a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv"),
            "docs/manuscript_core_paired_statistics.csv",
        ),
        planned_outputs=("docs/manuscript_case_shortlist.csv",),
    ),
    AnalysisSpec(
        analysis_id="D16",
        phase="Phase 6",
        analysis_name="Representation-switch cases if present",
        analysis_class="case_layer",
        reviewer_attack_vector="case_selection_bias_challenge",
        local_inputs=(
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet",
        ),
        prerequisite_outputs=(
            "docs/task2_fm_k562_local_comparison_summary.csv",
            "docs/task2_fm_k562_paired_tests.csv",
            "docs/manuscript_core_paired_statistics.csv",
        ),
        planned_outputs=("docs/manuscript_case_shortlist.csv",),
        conditional=True,
        local_only=True,
    ),
)


def rel_exists(rel_path: str) -> bool:
    return (ROOT / rel_path).is_file()


def list_specs() -> list[dict[str, object]]:
    return [asdict(spec) for spec in ANALYSIS_SPECS]


def iter_check_rows() -> Iterable[dict[str, object]]:
    for spec in ANALYSIS_SPECS:
        local_missing = [path for path in spec.local_inputs if not rel_exists(path)]
        prereq_present = [path for path in spec.prerequisite_outputs if rel_exists(path)]
        prereq_missing = [path for path in spec.prerequisite_outputs if not rel_exists(path)]
        yield {
            "analysis_id": spec.analysis_id,
            "phase": spec.phase,
            "analysis_class": spec.analysis_class,
            "local_inputs_missing": len(local_missing),
            "local_missing_paths": local_missing,
            "prerequisite_outputs_present": len(prereq_present),
            "prerequisite_outputs_missing": len(prereq_missing),
            "prerequisite_missing_paths": prereq_missing,
            "local_only": spec.local_only,
            "conditional": spec.conditional,
            "last_resort_upstream": list(spec.last_resort_upstream),
        }


def print_table(rows: Sequence[dict[str, object]]) -> None:
    headers = [
        "analysis_id",
        "phase",
        "analysis_class",
        "local_inputs_missing",
        "prerequisite_outputs_missing",
        "local_only",
        "conditional",
    ]
    widths = {header: len(header) for header in headers}
    for row in rows:
        for header in headers:
            widths[header] = max(widths[header], len(str(row[header])))
    print(" ".join(header.ljust(widths[header]) for header in headers))
    print(" ".join("-" * widths[header] for header in headers))
    for row in rows:
        print(" ".join(str(row[header]).ljust(widths[header]) for header in headers))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Registry and local-input checks for manuscript-support analyses")
    parser.add_argument("--json", action="store_true", help="emit JSON instead of the default table output")
    parser.add_argument(
        "--check-inputs",
        action="store_true",
        help="check required local upstream inputs; prerequisite outputs are reported but do not count as missing inputs",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.check_inputs:
        rows = list(iter_check_rows())
        if args.json:
            print(json.dumps(rows, indent=2, sort_keys=True))
        else:
            print_table(rows)
        has_missing = any(int(row["local_inputs_missing"]) > 0 for row in rows)
        return 1 if has_missing else 0

    payload = list_specs()
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print_table(
            [
                {
                    "analysis_id": spec.analysis_id,
                    "phase": spec.phase,
                    "analysis_class": spec.analysis_class,
                    "local_inputs_missing": 0,
                    "prerequisite_outputs_missing": len(spec.prerequisite_outputs),
                    "local_only": spec.local_only,
                    "conditional": spec.conditional,
                }
                for spec in ANALYSIS_SPECS
            ]
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
