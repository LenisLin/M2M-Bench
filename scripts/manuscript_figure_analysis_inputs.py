#!/usr/bin/env python3
"""
Build compact figure-facing analysis inputs from the active manuscript-support outputs.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


DEFAULT_ACTIVE_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active")

# Use a simple conservative threshold so supported slices are not driven by tiny strata.
B6_SUPPORTED_N_QUERIES_THRESHOLD = 5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create figure-facing analysis inputs from active manuscript outputs.")
    parser.add_argument("--active-root", type=Path, default=DEFAULT_ACTIVE_ROOT)
    return parser.parse_args()


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def build_a1_outputs(group_bridge_dir: Path, analysis_dir: Path) -> list[str]:
    bridge_path = group_bridge_dir / "task1_internal_vs_cross_group_bridge.csv"
    require_path(bridge_path, "A1 bridge table")

    bridge = pd.read_csv(bridge_path)
    require_columns(
        bridge,
        {
            "dataset",
            "cell_line",
            "target",
            "metric_name",
            "internal_value",
            "cross_value",
            "cross_minus_internal",
        },
        "A1 bridge table",
    )
    bridge_complete = bridge.loc[
        bridge[["internal_value", "cross_value", "cross_minus_internal"]].notna().all(axis=1)
    ].copy()

    dataset_cellline = (
        bridge_complete.groupby(["metric_name", "dataset", "cell_line"], dropna=False, observed=True, sort=True)
        .agg(
            n_groups=("target", "size"),
            internal_mean=("internal_value", "mean"),
            cross_mean=("cross_value", "mean"),
            delta_mean=("cross_minus_internal", "mean"),
            delta_median=("cross_minus_internal", "median"),
        )
        .reset_index()
    )
    dataset_cellline = dataset_cellline.sort_values(
        ["metric_name", "dataset", "cell_line"],
        kind="mergesort",
    ).reset_index(drop=True)

    metric_overview = (
        bridge_complete.groupby(["metric_name"], dropna=False, observed=True, sort=True)
        .agg(
            n_groups=("target", "size"),
            internal_mean=("internal_value", "mean"),
            cross_mean=("cross_value", "mean"),
            delta_mean=("cross_minus_internal", "mean"),
            delta_median=("cross_minus_internal", "median"),
        )
        .reset_index()
        .rename(
        columns={
            "n_groups": "n_groups_total",
            "internal_mean": "internal_mean_overall",
            "cross_mean": "cross_mean_overall",
            "delta_mean": "delta_mean_overall",
        }
    )
    )
    metric_overview = metric_overview.rename(columns={"delta_median": "delta_median_overall"})
    metric_overview = metric_overview.sort_values(["metric_name"], kind="mergesort").reset_index(drop=True)

    dataset_cellline_path = analysis_dir / "figure2_a1_dataset_cellline_overview.csv"
    metric_overview_path = analysis_dir / "figure2_a1_metric_overview.csv"
    dataset_cellline.to_csv(dataset_cellline_path, index=False)
    metric_overview.to_csv(metric_overview_path, index=False)
    return [metric_overview_path.name, dataset_cellline_path.name]


def build_a2_outputs(group_bridge_dir: Path, analysis_dir: Path) -> list[str]:
    bridge_path = group_bridge_dir / "task1_task2_group_bridge.csv"
    require_path(bridge_path, "A2 bridge table")

    bridge = pd.read_csv(bridge_path)
    require_columns(
        bridge,
        {
            "dataset",
            "cell_line",
            "target",
            "task1_metric_name",
            "task2_metric_name",
            "task1_value",
            "task2_value",
            "task2_minus_task1",
        },
        "A2 bridge table",
    )
    bridge_complete = bridge.loc[
        bridge[["task1_value", "task2_value", "task2_minus_task1"]].notna().all(axis=1)
    ].copy()

    dataset_cellline = (
        bridge_complete.groupby(
            ["task1_metric_name", "task2_metric_name", "dataset", "cell_line"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_groups=("target", "size"),
            task1_mean=("task1_value", "mean"),
            task2_mean=("task2_value", "mean"),
            delta_mean=("task2_minus_task1", "mean"),
            delta_median=("task2_minus_task1", "median"),
        )
        .reset_index()
    )
    dataset_cellline = dataset_cellline.sort_values(
        ["task1_metric_name", "task2_metric_name", "dataset", "cell_line"],
        kind="mergesort",
    ).reset_index(drop=True)

    metricpair_overview = (
        bridge_complete.groupby(
            ["task1_metric_name", "task2_metric_name"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_groups=("target", "size"),
            task1_mean=("task1_value", "mean"),
            task2_mean=("task2_value", "mean"),
            delta_mean=("task2_minus_task1", "mean"),
            delta_median=("task2_minus_task1", "median"),
        )
        .reset_index()
        .rename(
        columns={
            "n_groups": "n_groups_total",
            "task1_mean": "task1_mean_overall",
            "task2_mean": "task2_mean_overall",
            "delta_mean": "delta_mean_overall",
        }
    )
    )
    metricpair_overview = metricpair_overview.rename(columns={"delta_median": "delta_median_overall"})
    metricpair_overview = metricpair_overview.sort_values(
        ["task1_metric_name", "task2_metric_name"],
        kind="mergesort",
    ).reset_index(drop=True)

    dataset_cellline_path = analysis_dir / "figure3_a2_dataset_cellline_overview.csv"
    metricpair_overview_path = analysis_dir / "figure3_a2_metricpair_overview.csv"
    dataset_cellline.to_csv(dataset_cellline_path, index=False)
    metricpair_overview.to_csv(metricpair_overview_path, index=False)
    return [metricpair_overview_path.name, dataset_cellline_path.name]


def build_b6_outputs(sensitivity_dir: Path, analysis_dir: Path) -> list[str]:
    base_path = sensitivity_dir / "task2_c2g_dose_time_sensitivity.csv"
    require_path(base_path, "B6 sensitivity table")

    base = pd.read_csv(
        base_path,
        usecols=[
            "dataset",
            "cell_line",
            "target",
            "retrieval_metric_name",
            "dose_value",
            "time",
            "retrieval_value",
        ],
    )

    dataset_cellline = (
        base.groupby(
            ["dataset", "cell_line", "retrieval_metric_name", "dose_value", "time"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_queries=("target", "size"),
            mean_value=("retrieval_value", "mean"),
            median_value=("retrieval_value", "median"),
        )
        .reset_index()
    )
    dataset_cellline = dataset_cellline.sort_values(
        ["dataset", "cell_line", "retrieval_metric_name", "dose_value", "time"],
        kind="mergesort",
    ).reset_index(drop=True)

    metric_dose_time = (
        base.groupby(
            ["retrieval_metric_name", "dose_value", "time"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .agg(
            n_queries=("target", "size"),
            mean_value=("retrieval_value", "mean"),
            median_value=("retrieval_value", "median"),
        )
        .reset_index()
        .rename(
        columns={
            "n_queries": "n_queries_total",
            "mean_value": "mean_value_overall",
            "median_value": "median_value_overall",
        }
    )
    )
    metric_dose_time = metric_dose_time.sort_values(
        ["retrieval_metric_name", "dose_value", "time"],
        kind="mergesort",
    ).reset_index(drop=True)

    supported_slices = (
        base.groupby(
            ["dataset", "cell_line", "target", "retrieval_metric_name", "dose_value", "time"],
            dropna=False,
            observed=True,
            sort=True,
        )
        .size()
        .rename("n_queries")
        .reset_index()
    )
    supported_slices["support_flag"] = supported_slices["n_queries"].map(
        lambda value: "supported" if int(value) >= B6_SUPPORTED_N_QUERIES_THRESHOLD else "below_threshold"
    )
    supported_slices = supported_slices.sort_values(
        ["dataset", "cell_line", "target", "retrieval_metric_name", "dose_value", "time"],
        kind="mergesort",
    ).reset_index(drop=True)

    metric_dose_time_path = analysis_dir / "figure3_b6_metric_dose_time_overview.csv"
    dataset_cellline_path = analysis_dir / "figure3_b6_dataset_cellline_overview.csv"
    supported_slices_path = analysis_dir / "figure3_b6_supported_slices.csv"
    metric_dose_time.to_csv(metric_dose_time_path, index=False)
    dataset_cellline.to_csv(dataset_cellline_path, index=False)
    supported_slices.to_csv(supported_slices_path, index=False)
    return [metric_dose_time_path.name, dataset_cellline_path.name, supported_slices_path.name]


def write_analysis_note(active_root: Path, analysis_dir: Path, created_files: list[str]) -> str:
    note_path = analysis_dir / "figure_analysis_update.md"
    lines = [
        "# Figure Analysis Update",
        "",
        f"- Relocated active outputs: `{active_root / 'group_bridge'}`, `{active_root / 'sensitivity'}`, `{active_root / 'notes'}`.",
        f"- Figure-facing tables created: {', '.join(f'`{name}`' for name in created_files)}.",
        f"- B6 supported-slice threshold: `n_queries >= {B6_SUPPORTED_N_QUERIES_THRESHOLD}`.",
        "- Blocker for figure drafting: none.",
    ]
    note_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return note_path.name


def write_manifest(active_root: Path, created_files: list[str]) -> str:
    manifest_path = active_root / "manifest.json"
    payload = {
        "batch_name": "manuscript_active",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "files": [
            "group_bridge/task1_internal_vs_cross_group_bridge.csv",
            "group_bridge/task1_internal_vs_cross_group_bridge_summary.csv",
            "group_bridge/task1_task2_group_bridge.csv",
            "group_bridge/task1_task2_group_bridge_summary.csv",
            "sensitivity/task2_c2g_dose_time_sensitivity.csv",
            "sensitivity/task2_c2g_dose_time_sensitivity_summary.csv",
            "notes/manuscript_active_experiment_update.md",
            *[f"analysis/{name}" for name in created_files],
        ],
        "source_scripts": [
            "scripts/manuscript_a1_task1_internal_cross.py",
            "scripts/manuscript_a2_task1_task2_bridge.py",
            "scripts/manuscript_b6_c2g_dose_time_sensitivity.py",
            "scripts/manuscript_figure_analysis_inputs.py",
        ],
        "note": "downstream manuscript-support outputs",
    }
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return manifest_path.name


def main() -> int:
    args = parse_args()
    active_root = args.active_root.resolve()
    group_bridge_dir = active_root / "group_bridge"
    sensitivity_dir = active_root / "sensitivity"
    analysis_dir = active_root / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    created_files: list[str] = []
    created_files.extend(build_a1_outputs(group_bridge_dir, analysis_dir))
    created_files.extend(build_a2_outputs(group_bridge_dir, analysis_dir))
    created_files.extend(build_b6_outputs(sensitivity_dir, analysis_dir))
    created_files.append(write_analysis_note(active_root, analysis_dir, created_files.copy()))
    created_files.append(write_manifest(active_root, created_files.copy()))

    print(f"Wrote {len(created_files)} analysis artifacts under {analysis_dir}")
    for name in created_files:
        print(name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
