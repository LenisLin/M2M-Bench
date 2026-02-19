from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "dev" / "generate_readme.py"


def _write_project_yaml(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "name: demo",
                "title: Demo Title",
                "one_liner: Demo one-liner",
                "domain: Demo domain",
                "stage: Demo stage",
                "owner: Demo owner",
                "license: Demo license",
                "entrypoints: []",
                "datasets: []",
                "outputs: []",
                "components: []",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _run_generator(
    tmp_path: Path,
    project_path: Path,
    template_path: Path,
    readme_path: Path,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(SCRIPT_PATH),
        "--project-file",
        str(project_path),
        "--template-file",
        str(template_path),
        "--readme-file",
        str(readme_path),
    ]
    if check:
        command.append("--check")

    return subprocess.run(command, cwd=tmp_path, check=False, capture_output=True, text=True)


def test_generate_readme_write_mode_updates_marker_block(tmp_path: Path) -> None:
    project_path = tmp_path / "project.yaml"
    _write_project_yaml(project_path)

    docs_dir = tmp_path / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    template_path = docs_dir / "readme.template.md"
    template_path.write_text("## {{ title }}\n\n{{ one_liner }}\n", encoding="utf-8")

    readme_path = tmp_path / "README.md"
    readme_path.write_text(
        "# Demo\n\n<!-- AVCP:README:START -->\nOLD\n<!-- AVCP:README:END -->\n",
        encoding="utf-8",
    )

    result = _run_generator(tmp_path, project_path, template_path, readme_path)
    assert result.returncode == 0, result.stderr

    updated = readme_path.read_text(encoding="utf-8")
    assert "## Demo Title" in updated
    assert "Demo one-liner" in updated
    assert "<!-- AVCP:README:START -->" in updated
    assert "<!-- AVCP:README:END -->" in updated


def test_generate_readme_check_mode_fails_when_readme_is_stale(tmp_path: Path) -> None:
    project_path = tmp_path / "project.yaml"
    _write_project_yaml(project_path)

    docs_dir = tmp_path / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)
    template_path = docs_dir / "readme.template.md"
    template_path.write_text("## {{ title }}\n\n{{ one_liner }}\n", encoding="utf-8")

    readme_path = tmp_path / "README.md"
    readme_path.write_text(
        "# Demo\n\n<!-- AVCP:README:START -->\nSTALE\n<!-- AVCP:README:END -->\n",
        encoding="utf-8",
    )

    result = _run_generator(tmp_path, project_path, template_path, readme_path, check=True)
    assert result.returncode != 0
