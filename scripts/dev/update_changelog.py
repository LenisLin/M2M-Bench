from __future__ import annotations

import argparse
from pathlib import Path


UNRELEASED_HEADER = "## Unreleased"


def insert_unreleased_entry(md_text: str, entry: str) -> str:
    entry = entry.strip()
    if not entry.startswith("- "):
        entry = "- " + entry

    lines = md_text.splitlines()

    # Ensure Unreleased exists
    if UNRELEASED_HEADER not in md_text:
        # Insert after first H1 if present, else at top
        out = []
        inserted = False
        for i, line in enumerate(lines):
            out.append(line)
            if (not inserted) and line.startswith("# "):
                out.append("")
                out.append(UNRELEASED_HEADER)
                out.append(entry)
                inserted = True
        if not inserted:
            out = [UNRELEASED_HEADER, entry, ""] + lines
        return "\n".join(out).rstrip() + "\n"

    # Insert directly after "## Unreleased" line, before next header of same or higher level
    out = []
    inserted = False
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        if (line.strip() == UNRELEASED_HEADER) and (not inserted):
            # Insert on next line
            out.append(entry)
            inserted = True
        i += 1

    # If entry already present, do not duplicate (simple check)
    # (Optional: can be enhanced)
    if md_text.find(entry) != -1:
        return md_text if md_text.endswith("\n") else (md_text + "\n")

    return "\n".join(out).rstrip() + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--entry", required=True, help="Changelog bullet, e.g. '- feat(x): ...'")
    parser.add_argument("--file", default="docs/dev_log.md", help="Markdown changelog file path")
    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# Changelog\n\n## Unreleased\n", encoding="utf-8")

    md = path.read_text(encoding="utf-8")
    updated = insert_unreleased_entry(md, args.entry)
    path.write_text(updated, encoding="utf-8")


if __name__ == "__main__":
    main()
