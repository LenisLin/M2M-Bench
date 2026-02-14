#!/usr/bin/env python3
"""
Legacy-compatible Task1 retrieval entrypoint.

Use this script for only the instance-level retrieval step.
For the full readable Task1 workflow, prefer `scripts/task1_pipeline.py`.
"""

from __future__ import annotations

from _script_bootstrap import setup_project_imports

setup_project_imports()

from m2m_bench.task1.retrieval_instance import main as run_task1_retrieval_cli


if __name__ == "__main__":
    run_task1_retrieval_cli()
