#!/usr/bin/env python3
"""
Shared bootstrap helper for local script execution.

All analysis scripts in `scripts/` can call `setup_project_imports()` to ensure:
1) `src/` is importable (for `m2m_bench.*` modules),
2) `scripts/` is importable (for script-to-script reuse).
"""

from __future__ import annotations

from pathlib import Path
import sys


def setup_project_imports() -> Path:
    """
    Add project import paths and return repository root path.

    Returns
    -------
    Path
        Repository root (parent of `scripts/`).
    """
    scripts_dir = Path(__file__).resolve().parent
    root_dir = scripts_dir.parent
    src_dir = root_dir / "src"

    for path in (src_dir, scripts_dir):
        path_str = str(path)
        if path_str not in sys.path:
            sys.path.insert(0, path_str)

    return root_dir
