from __future__ import annotations

import argparse
from pathlib import Path

CHANGELOG_HEADER = "# Changelog"
UNRELEASED_HEADER = "## Unreleased"


def _normalize_entry(entry: str) -> str:
    stripped = entry.strip()
    if not stripped.startswith("- "):
        return "- " + stripped
    return stripped


def _render_lines(lines: list[str]) -> str:
    return "\n".join(lines).rstrip("\n") + "\n"


def _find_line(lines: list[str], expected: str) -> int | None:
    for index, line in enumerate(lines):
        if line.strip() == expected:
            return index
    return None


def _ensure_changelog_header(lines: list[str]) -> list[str]:
    if _find_line(lines, CHANGELOG_HEADER) is not None:
        return lines

    if lines:
        return [CHANGELOG_HEADER, "", *lines]
    return [CHANGELOG_HEADER, ""]


def insert_unreleased_entry(md_text: str, entry: str) -> str:
    entry_line = _normalize_entry(entry)
    lines = md_text.splitlines()

    if any(line.strip() == entry_line for line in lines):
        return _render_lines(lines)

    lines = _ensure_changelog_header(lines)

    unreleased_index = _find_line(lines, UNRELEASED_HEADER)
    if unreleased_index is None:
        changelog_index = _find_line(lines, CHANGELOG_HEADER)
        if changelog_index is None:
            raise ValueError("Missing changelog header after normalization.")
        insert_at = changelog_index + 1
        lines[insert_at:insert_at] = ["", UNRELEASED_HEADER]
        unreleased_index = _find_line(lines, UNRELEASED_HEADER)
        if unreleased_index is None:
            raise ValueError("Failed to create Unreleased section.")

    end_index = len(lines)
    for index in range(unreleased_index + 1, len(lines)):
        if lines[index].startswith("## "):
            end_index = index
            break

    insert_at = end_index
    while insert_at > unreleased_index + 1 and lines[insert_at - 1].strip() == "":
        insert_at -= 1

    lines.insert(insert_at, entry_line)
    return _render_lines(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--entry", required=True, help="Changelog bullet, e.g. '- feat(x): ...'")
    parser.add_argument("--file", default="docs/dev_log.md", help="Markdown changelog file path")
    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{CHANGELOG_HEADER}\n", encoding="utf-8")

    existing = path.read_text(encoding="utf-8")
    updated = insert_unreleased_entry(existing, args.entry)
    path.write_text(updated, encoding="utf-8")


if __name__ == "__main__":
    main()
