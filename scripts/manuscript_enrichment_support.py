#!/usr/bin/env python3
# SCRIPT_HEADER_CONTRACT
# Script: scripts/manuscript_enrichment_support.py
# Purpose: Provide a registry and local-input checks for the manuscript-only
#   enrichment and composition layer without recomputing benchmark metrics.
# Inputs:
#   - existing audited Task2 outputs and manuscript planning files
# Outputs:
#   - stdout only
# Side Effects:
#   - none
# Failure Modes:
#   - missing required local benchmark inputs when --check-inputs is used
# Last Updated: 2026-03-13

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class EnrichmentSpec:
    analysis_id: str
    analysis_name: str
    reviewer_attack_vector: str
    local_inputs: tuple[str, ...]
    prerequisite_outputs: tuple[str, ...]
    planned_outputs: tuple[str, ...]
    unresolved_external_assets: tuple[str, ...]
    notes: str


ENRICHMENT_SPECS: tuple[EnrichmentSpec, ...] = (
    EnrichmentSpec(
        analysis_id="B10",
        analysis_name="Dataset / cell-line composition analysis",
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
        unresolved_external_assets=(),
        notes="Interpret enrichment only after or alongside this composition context.",
    ),
    EnrichmentSpec(
        analysis_id="B9",
        analysis_name="Enrichment for high-/low-concordance target sets",
        reviewer_attack_vector="survivor_bias_enrichment_challenge",
        local_inputs=(
            "runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv",
            "runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv",
            "data/task2_snapshot_v2/task2_pairs_coverage.csv",
        ),
        prerequisite_outputs=(
            "docs/task2_dataset_cellline_composition_analysis.csv",
            "docs/task2_target_pool_background_registry.csv",
        ),
        planned_outputs=(
            "docs/task2_target_concordance_rankings.csv",
            "docs/task2_high_low_concordance_enrichment.csv",
        ),
        unresolved_external_assets=(
            "annotation asset for target-family or pathway membership",
            "annotation asset for degree-aware or network context, if used",
        ),
        notes="Matched and evaluable target pools are required backgrounds; never use the whole genome.",
    ),
)


def rel_exists(rel_path: str) -> bool:
    return (ROOT / rel_path).is_file()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Registry and local-input checks for enrichment-layer manuscript support")
    parser.add_argument("--json", action="store_true", help="emit JSON instead of human-readable text")
    parser.add_argument("--check-inputs", action="store_true", help="check required local benchmark inputs")
    return parser.parse_args()


def check_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for spec in ENRICHMENT_SPECS:
        local_missing = [path for path in spec.local_inputs if not rel_exists(path)]
        prereq_missing = [path for path in spec.prerequisite_outputs if not rel_exists(path)]
        rows.append(
            {
                "analysis_id": spec.analysis_id,
                "analysis_name": spec.analysis_name,
                "reviewer_attack_vector": spec.reviewer_attack_vector,
                "local_missing_paths": local_missing,
                "local_inputs_missing": len(local_missing),
                "prerequisite_missing_paths": prereq_missing,
                "prerequisite_outputs_missing": len(prereq_missing),
                "unresolved_external_assets": list(spec.unresolved_external_assets),
                "notes": spec.notes,
            }
        )
    return rows


def print_text(rows: Sequence[dict[str, object]]) -> None:
    for row in rows:
        print(f"{row['analysis_id']}: {row['analysis_name']}")
        print(f"  reviewer_attack_vector: {row['reviewer_attack_vector']}")
        print(f"  local_inputs_missing: {row['local_inputs_missing']}")
        if row["local_missing_paths"]:
            print(f"  local_missing_paths: {', '.join(row['local_missing_paths'])}")
        print(f"  prerequisite_outputs_missing: {row['prerequisite_outputs_missing']}")
        if row["prerequisite_missing_paths"]:
            print(f"  prerequisite_missing_paths: {', '.join(row['prerequisite_missing_paths'])}")
        if row["unresolved_external_assets"]:
            print(f"  unresolved_external_assets: {', '.join(row['unresolved_external_assets'])}")
        print(f"  notes: {row['notes']}")


def main() -> int:
    args = parse_args()
    if args.check_inputs:
        rows = check_rows()
        if args.json:
            print(json.dumps(rows, indent=2, sort_keys=True))
        else:
            print_text(rows)
        return 1 if any(int(row["local_inputs_missing"]) > 0 for row in rows) else 0

    payload = [asdict(spec) for spec in ENRICHMENT_SPECS]
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print_text(check_rows())
    return 0


if __name__ == "__main__":
    sys.exit(main())
