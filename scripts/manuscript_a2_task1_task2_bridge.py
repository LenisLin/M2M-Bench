#!/usr/bin/env python3
"""
Phase 1 manuscript-support skeleton for A2: Task1 to Task2 bridge scaffold.

Analysis purpose:
- Prepare a lawful shared-slice bridge comparison between Task1 and Task2 on
  exact shared keys only.

Why this script exists in the manuscript workflow:
- The manuscript needs a downstream bridge scaffold that makes exact shared-key
  coverage explicit while leaving unresolved bridge semantics under human review.

Lawful input scope:
- reviewed local downstream Task1 outputs
- reviewed local S4 `task2_group_concordance.csv`
- reviewed local S6 `task2_group_leaderboard.csv`
- reviewed local S6 `task2_retrieval_leaderboard.csv`
- optional Task1-side keyed inputs may be passed if reviewed and lawful;
  `task1_group_internal.parquet` is only one possible candidate path

Prohibited behavior:
- do not invent bridge rows
- do not infer metric-family mappings silently
- do not convert unresolved bridge semantics into hidden defaults
- do not mutate benchmark-stage outputs

Intended output semantics:
- one row per lawful bridge attempt after exact shared-key review
- successful rows carry Task1-side and Task2-side values only when both sides
  exist on exact reviewed keys
- unavailable or blocked rows must remain explicit rather than inferred

What remains TODO in later implementation:
- confirm whether any reviewed Task1-side input from the allowed input set can
  expose exact shared keys centered on `dataset, cell_line, target_token,
  representation`
- if not, review whether an additional Task1-side keyed export is lawful
- freeze metric-family bridge semantics between Task1 and Task2
- build the final bridge table only after those review points are settled
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
NAS_RUNS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs")
MANUSCRIPT_PHASE1_ROOT = NAS_RUNS_ROOT / "manuscript_phase1" / "a2_task1_task2_bridge"

DEFAULT_TASK1_LEADERBOARD_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv"
)
DEFAULT_TASK1_GROUP_INTERNAL_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_group_internal.parquet"
)
DEFAULT_TASK2_GROUP_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv"
)
DEFAULT_TASK2_GROUP_LEADERBOARD_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv"
)
DEFAULT_TASK2_RETRIEVAL_LEADERBOARD_PATH = Path(
    "/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv"
)
DEFAULT_OUTPUT_PATH = MANUSCRIPT_PHASE1_ROOT / "task1_task2_shared_scope_bridge_table.csv"

TASK1_LEADERBOARD_REQUIRED_COLUMNS = [
    "scope",
    "dataset_or_direction",
    "perturbation_type",
    "representation",
    "metric_name",
    "metric_value",
    "n_total",
    "n_valid",
    "n_excluded",
    "N_gallery_max",
]
TASK1_GROUP_INTERNAL_MINIMAL_COLUMNS = [
    "dataset_or_direction",
    "cell_line",
    "target_token",
    "representation",
]
TASK2_GROUP_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "target_token",
    "representation",
    "n_chem_instances_used",
    "n_gen_instances_used",
    "cosine_centroid",
    "pcc_centroid",
    "edist_biascorr",
]
TASK2_GROUP_LEADERBOARD_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "representation",
    "mean_cosine_centroid",
    "mean_pcc_centroid",
    "mean_edist_biascorr",
]
TASK2_RETRIEVAL_LEADERBOARD_REQUIRED_COLUMNS = [
    "dataset",
    "cell_line",
    "direction",
    "representation",
    "n_total",
    "n_valid",
    "mean_mrr_corrected",
    "mean_hit1_corrected",
    "mean_hit5_corrected",
    "mean_hit10_corrected",
]
BRIDGE_OUTPUT_COLUMNS = [
    "dataset",
    "cell_line",
    "target_token",
    "representation",
    "task1_metric_family",
    "task1_metric_name",
    "task1_value",
    "task2_analysis_family",
    "task2_metric_name",
    "task2_value",
    "bridge_status",
    "blocked_reason",
    "notes",
]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the A2 scaffold."""
    parser = argparse.ArgumentParser(description="Phase 1 skeleton for A2 Task1 to Task2 bridge.")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--task1-leaderboard-path", type=Path, default=DEFAULT_TASK1_LEADERBOARD_PATH)
    parser.add_argument("--task1-group-internal-path", type=Path, default=DEFAULT_TASK1_GROUP_INTERNAL_PATH)
    parser.add_argument("--task2-group-path", type=Path, default=DEFAULT_TASK2_GROUP_PATH)
    parser.add_argument("--task2-group-leaderboard-path", type=Path, default=DEFAULT_TASK2_GROUP_LEADERBOARD_PATH)
    parser.add_argument("--task2-retrieval-leaderboard-path", type=Path, default=DEFAULT_TASK2_RETRIEVAL_LEADERBOARD_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    return parser.parse_args()


def resolve_path(project_root: Path, raw_path: Path) -> Path:
    """Resolve a CLI path against the requested project root."""
    return raw_path if raw_path.is_absolute() else (project_root / raw_path)


def load_table(path: Path, label: str) -> pd.DataFrame:
    """Load a CSV or parquet table without hidden coercion."""
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if path.suffix == ".csv":
        return pd.read_csv(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"{label} must be CSV or parquet, got: {path}")


def ensure_required_columns(frame: pd.DataFrame, required: Iterable[str], label: str) -> None:
    """Fail fast when a required frozen column is absent."""
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def ensure_unique_keys(frame: pd.DataFrame, key_columns: list[str], label: str) -> None:
    """Reject duplicate rows on exact shared-key inputs."""
    duplicated = frame.duplicated(subset=key_columns, keep=False)
    if duplicated.any():
        raise ValueError(f"{label} has duplicate key rows on {key_columns}")


def load_reviewed_inputs(args: argparse.Namespace, project_root: Path) -> dict[str, pd.DataFrame | None]:
    """Load the reviewed local A2 inputs, keeping the last-resort Task1 group path optional."""
    task1_group_internal_path = resolve_path(project_root, args.task1_group_internal_path)
    task1_group_internal = None
    if task1_group_internal_path.exists():
        task1_group_internal = load_table(task1_group_internal_path, "optional Task1 group internal table")

    return {
        "task1_leaderboard": load_table(
            resolve_path(project_root, args.task1_leaderboard_path),
            "Task1 leaderboard",
        ),
        "task1_group_internal": task1_group_internal,
        "task2_group": load_table(resolve_path(project_root, args.task2_group_path), "S4 Task2 group concordance"),
        "task2_group_leaderboard": load_table(
            resolve_path(project_root, args.task2_group_leaderboard_path),
            "S6 Task2 group leaderboard",
        ),
        "task2_retrieval_leaderboard": load_table(
            resolve_path(project_root, args.task2_retrieval_leaderboard_path),
            "S6 Task2 retrieval leaderboard",
        ),
    }


def validate_inputs(inputs: dict[str, pd.DataFrame | None]) -> None:
    """Run the reviewed schema checks that are safe to freeze now."""
    ensure_required_columns(
        inputs["task1_leaderboard"],
        TASK1_LEADERBOARD_REQUIRED_COLUMNS,
        "Task1 leaderboard",
    )
    ensure_required_columns(inputs["task2_group"], TASK2_GROUP_REQUIRED_COLUMNS, "S4 Task2 group concordance")
    ensure_required_columns(
        inputs["task2_group_leaderboard"],
        TASK2_GROUP_LEADERBOARD_REQUIRED_COLUMNS,
        "S6 Task2 group leaderboard",
    )
    ensure_required_columns(
        inputs["task2_retrieval_leaderboard"],
        TASK2_RETRIEVAL_LEADERBOARD_REQUIRED_COLUMNS,
        "S6 Task2 retrieval leaderboard",
    )

    ensure_unique_keys(
        inputs["task2_group"],
        ["dataset", "cell_line", "target_token", "representation"],
        "S4 Task2 group concordance",
    )

    if inputs["task1_group_internal"] is not None:
        ensure_required_columns(
            inputs["task1_group_internal"],
            TASK1_GROUP_INTERNAL_MINIMAL_COLUMNS,
            "optional Task1 group internal table",
        )


def check_exact_shared_key_coverage(inputs: dict[str, pd.DataFrame | None]) -> pd.DataFrame:
    """
    Check whether exact shared-key coverage can be reviewed from current files.

    This function exists now because the planning docs explicitly require exact
    shared-key checking. The actual coverage logic remains blocked only when no
    reviewed Task1-side keyed input is currently available from the allowed
    input set.
    """
    task1_group_internal = inputs["task1_group_internal"]
    task2_group = inputs["task2_group"]

    if task1_group_internal is None:
        return pd.DataFrame(
            [
                {
                    "coverage_status": "blocked",
                    "blocked_reason": "no_reviewed_task1_side_exact_key_input_available",
                    "notes": (
                        "Current reviewed Task1 inputs loaded by this skeleton "
                        "do not include a keyed table that exposes exact shared "
                        "keys on `dataset, cell_line, target_token, "
                        "representation`. This is a general gating condition, "
                        "not a claim that one specific candidate path must "
                        "exist."
                    ),
                }
            ]
        )

    # TODO: Generalize this branch to whichever reviewed Task1-side keyed input
    # is eventually approved for the bridge. `task1_group_internal.parquet` is
    # only one optional candidate and must not become a hidden default.
    #
    # TODO: Normalize the Task1-side dataset field once the reviewed keyed input
    # schema is frozen.
    #
    # TODO: Intersect exact shared keys on
    # `dataset, cell_line, target_token, representation`.
    #
    # TODO: Report how much Task1-side and Task2-side support remains on the
    # retained exact-key slice before any bridge row is built.
    _ = task2_group
    return pd.DataFrame(columns=["coverage_status", "blocked_reason", "notes"])


def build_bridge_table(
    inputs: dict[str, pd.DataFrame | None],
    coverage: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build the future A2 bridge table.

    The metric-family bridge semantics remain unresolved by the planning docs,
    so this scaffold returns only the output schema.
    """
    _ = inputs
    _ = coverage

    # TODO: Freeze the reviewed Task1-to-Task2 metric-family bridge semantics.
    # The scaffold must not silently map Task1 retrieval metrics onto Task2
    # group or retrieval metrics.
    #
    # TODO: Build bridge rows only on exact shared keys that passed the coverage
    # review. No bridge row may be invented from summary-level coincidences.
    #
    # TODO: Emit `successful`, `unavailable`, or `blocked` row statuses exactly
    # as defined in the Phase 1 script spec.
    return pd.DataFrame(columns=BRIDGE_OUTPUT_COLUMNS)


def write_output(frame: pd.DataFrame, output_path: Path) -> None:
    """Write the manuscript-support bridge table."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)


def main() -> int:
    """Run the A2 skeleton through loading, validation, and placeholder output."""
    args = parse_args()
    project_root = resolve_path(Path.cwd(), args.project_root)
    inputs = load_reviewed_inputs(args, project_root)
    validate_inputs(inputs)
    coverage = check_exact_shared_key_coverage(inputs)
    output = build_bridge_table(inputs, coverage)
    write_output(output, resolve_path(project_root, args.output_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
