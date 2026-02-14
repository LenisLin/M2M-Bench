#!/usr/bin/env python3
"""
Legacy-compatible Task1 group-wise entrypoint.

Use this script for only the group-wise modality-gap step.
For the full readable Task1 workflow (group-wise + retrieval + confounder + set-level),
prefer `scripts/task1_pipeline.py`.
"""

from __future__ import annotations

from _script_bootstrap import setup_project_imports

setup_project_imports()

from m2m_bench.task1.groupwise_gap import main as run_task1_groupwise_cli


if __name__ == "__main__":
    run_task1_groupwise_cli()
