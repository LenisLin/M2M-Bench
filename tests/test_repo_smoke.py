from __future__ import annotations

import py_compile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FM_EXTRACTOR_DIR = REPO_ROOT / "scripts" / "fm_extractors"


def test_fm_extractors_are_present() -> None:
    scripts = sorted(FM_EXTRACTOR_DIR.glob("*.py"))
    assert scripts, "Expected FM extractor scripts to remain in the repo surface."


def test_fm_extractors_compile() -> None:
    scripts = sorted(FM_EXTRACTOR_DIR.glob("*.py"))
    for script in scripts:
        py_compile.compile(str(script), doraise=True)
