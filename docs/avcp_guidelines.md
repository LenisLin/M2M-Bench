# AVCP Guidelines

## Bridge
**Rule:** Any Python→R handover must be produced via `src/io/bridge.py::save_for_r()` unless explicitly waived in `docs/constraints.md`.

Output rules:
- Location: `paths.interim_viz_dir` from `config/config.yaml` (default: `data/interim_viz`)
- Format: `.parquet` for large, `.csv` for small
- Sidecar: always generate `<stem>_meta.json` with `file`, `primary_key`, `columns`, `provenance`
- No implicit index: always have an explicit primary key column

## Git and SemVer
SemVer applies to release artifacts (repo-level or component-level), not individual files.
Use Conventional Commits: `<type>(<scope>): <message>`.
Bump rules:
- MAJOR: breaking change in api/data contracts
- MINOR: backward-compatible feature
- PATCH: bugfix/perf/docs/tests not changing external contracts

## Changelog
Do NOT append blindly using shell echo.
Preferred:
- Provide unified diff patches, OR
- Use `scripts/dev/update_changelog.py` to insert under `## Unreleased`.
