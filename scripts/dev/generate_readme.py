from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Any

import yaml
from jinja2 import Environment, StrictUndefined

START_MARKER = "<!-- AVCP:README:START -->"
END_MARKER = "<!-- AVCP:README:END -->"


def _single_trailing_newline(text: str) -> str:
    return text.rstrip("\n") + "\n"


def _load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing project metadata file: {path}")

    parsed = yaml.safe_load(path.read_text(encoding="utf-8"))
    if parsed is None:
        return {}
    if not isinstance(parsed, dict):
        raise ValueError(f"Project metadata must be a mapping: {path}")
    return parsed


def _render_template(template_path: Path, context: dict[str, Any]) -> str:
    if not template_path.exists():
        raise FileNotFoundError(f"Missing README template file: {template_path}")

    environment = Environment(undefined=StrictUndefined, trim_blocks=True, lstrip_blocks=True)
    template = environment.from_string(template_path.read_text(encoding="utf-8"))
    return template.render(**context).strip("\n")


def _build_controlled_block(rendered: str) -> str:
    return f"{START_MARKER}\n{rendered}\n{END_MARKER}"


def _insert_block_near_top(readme_text: str, controlled_block: str) -> str:
    lines = readme_text.splitlines()
    block_lines = controlled_block.splitlines()

    if not lines:
        return _single_trailing_newline(controlled_block)

    h1_index = next((i for i, line in enumerate(lines) if line.startswith("# ")), None)
    if h1_index is None:
        merged = [*block_lines, "", *lines]
        return _single_trailing_newline("\n".join(merged))

    prefix = lines[: h1_index + 1]
    suffix = lines[h1_index + 1 :]

    merged = [*prefix, "", *block_lines]
    if suffix:
        merged.append("")
        merged.extend(suffix)

    return _single_trailing_newline("\n".join(merged))


def update_readme_marked_block(readme_text: str, rendered_block: str) -> str:
    controlled_block = _build_controlled_block(rendered_block)

    has_start = START_MARKER in readme_text
    has_end = END_MARKER in readme_text
    if has_start and has_end:
        pattern = re.compile(
            rf"{re.escape(START_MARKER)}.*?{re.escape(END_MARKER)}",
            re.DOTALL,
        )
        updated = pattern.sub(controlled_block, readme_text, count=1)
        return _single_trailing_newline(updated)

    return _insert_block_near_top(readme_text, controlled_block)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Exit non-zero if README is stale")
    parser.add_argument("--project-file", default="project.yaml")
    parser.add_argument("--template-file", default="docs/readme.template.md")
    parser.add_argument("--readme-file", default="README.md")
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    logger = logging.getLogger("scripts.dev.generate_readme")

    args = _parse_args()
    project_path = Path(args.project_file)
    template_path = Path(args.template_file)
    readme_path = Path(args.readme_file)

    context = _load_yaml(project_path)
    rendered = _render_template(template_path, context)
    current_readme = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""
    expected_readme = update_readme_marked_block(current_readme, rendered)

    normalized_current = _single_trailing_newline(current_readme) if current_readme else ""
    if args.check:
        if normalized_current != expected_readme:
            logger.error("README is out of date. Run: python scripts/dev/generate_readme.py")
            raise SystemExit(1)
        logger.info("README is up to date.")
        return

    readme_path.write_text(expected_readme, encoding="utf-8")
    logger.info("Updated README markers in %s", readme_path)


if __name__ == "__main__":
    main()
