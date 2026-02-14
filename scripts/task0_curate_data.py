#!/usr/bin/env python3
"""
Legacy-compatible Task0 entrypoint.

Use this script when you only want Task0 data curation.
For a fuller readable workflow (Task0 + optional Task0->Task1 attrition audit),
prefer `scripts/task0_pipeline.py`.
"""

from __future__ import annotations

from _script_bootstrap import setup_project_imports

setup_project_imports()

from m2m_bench.task0.curate import main as run_task0_cli


if __name__ == "__main__":
    run_task0_cli()
