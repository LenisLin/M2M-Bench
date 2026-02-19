# AVCP Guidelines

## Bridge
**Rule:** Any Python→R handover must be produced via `src/avcp_template/io/bridge.py::save_for_r()` unless explicitly waived in `docs/constraints.md`.

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
- `update_changelog.py` is the canonical way to add entries without breaking markdown structure or duplicating bullets.

## 4.1 Script Header Contract
All new or modified executable scripts under `scripts/` must include this header block at the top of file-level comments.

```text
# SCRIPT_HEADER_CONTRACT
# Script: <repo-relative-path>
# Purpose: <one-line objective>
# Inputs:
#   - <name>: <source/path/type>
# Outputs:
#   - <artifact>: <path/format>
# Side Effects:
#   - <created/modified paths>
# Config Dependencies:
#   - config/config.yaml::<key.path>
# Execution:
#   - python <script> [args]
# Failure Modes:
#   - <condition> -> <behavior/exit code>
# Last Updated: <YYYY-MM-DD>
```

Rules:
- Keep it synchronized with actual script behavior in the same patch.
- Do not hardcode absolute paths; reference `config/config.yaml` keys.
- If a field is not applicable, write `N/A` explicitly.

## 4.2 AI Role Positioning: Objectivity and Evidence
The AI must operate as an objective engineering collaborator.

Behavioral requirements:
- No flattery, appeasement, or persuasive language that is not technically relevant.
- No fabricated outcomes, metrics, code behavior, or experiment conclusions.
- No certainty claims without verifiable evidence.

Conclusion requirements:
- Every non-trivial conclusion must include a numbered evidence list.
- Each evidence item should reference a concrete source:
  - file path + line/section,
  - command output,
  - table/metric artifact,
  - or external citation (if used).
- If evidence is incomplete, explicitly say `Insufficient evidence` and provide the next verification action.

Recommended conclusion format:
1. Conclusion
2. Evidence
3. Confidence / uncertainty
4. Next verification step (if needed)

## README (Derived Artifact)
- `README.md` is a derived artifact generated from `project.yaml` + `docs/readme.template.md`.
- Use `scripts/dev/generate_readme.py` in write mode locally; CI enforces `--check`.
