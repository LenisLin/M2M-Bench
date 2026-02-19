from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "dev" / "update_changelog.py"


def _run_updater(changelog_path: Path, entry: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--file",
            str(changelog_path),
            "--entry",
            entry,
        ],
        check=False,
        capture_output=True,
        text=True,
    )


def test_insert_unreleased_entry_avoids_duplicate(tmp_path: Path) -> None:
    path = tmp_path / "dev_log.md"
    initial = "# Changelog\n\n## Unreleased\n- feat: thing\n"
    path.write_text(initial, encoding="utf-8")

    result = _run_updater(path, "- feat: thing")
    assert result.returncode == 0, result.stderr
    assert path.read_text(encoding="utf-8") == initial


def test_insert_unreleased_entry_appends_at_end_of_unreleased(tmp_path: Path) -> None:
    path = tmp_path / "dev_log.md"
    initial = "# Changelog\n\n## Unreleased\n- first\n\n## 0.1.0\n- old\n"
    path.write_text(initial, encoding="utf-8")

    result = _run_updater(path, "- second")
    assert result.returncode == 0, result.stderr

    expected = "# Changelog\n\n## Unreleased\n- first\n- second\n\n## 0.1.0\n- old\n"
    assert path.read_text(encoding="utf-8") == expected


def test_insert_unreleased_entry_creates_headers_when_missing(tmp_path: Path) -> None:
    path = tmp_path / "dev_log.md"
    path.write_text("Some intro\n", encoding="utf-8")

    result = _run_updater(path, "feat: init")
    assert result.returncode == 0, result.stderr

    updated = path.read_text(encoding="utf-8")
    assert updated.startswith("# Changelog\n")
    assert "## Unreleased\n" in updated
    assert "- feat: init\n" in updated
    assert updated.endswith("\n")
