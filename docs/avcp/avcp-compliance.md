# AVCP Compliance Mapping

## AVCP invariants we enforce
1) Root sentinel: project root must contain .avcp_root.
2) Dry-run first: every CODEX execution prints planned writes before writing.
3) Run isolation: outputs must go to runs/<run_id>/... (no silent overwrites).
4) Evidence-first: every claim must be supported by an artifact (report + manifest + samples).
5) Fail-fast: contract violations cause non-zero exit after writing minimal evidence artifacts.
6) No scope creep: tasks are bound by docs/contracts/task-definitions.md and policy files.
7) E-min pack 与 question pack index 属于 Tier-1 evidence outputs；禁止在没有 snapshot validation 通过的情况下生成。

## Required artifacts per run
- run_manifest.json (inputs, outputs, versions, hashes)
- manifest.json (sha256+size of generated files)
- report.md (fact-only)
- sample evidence files (small CSV/JSONL)
- audit_assertions.json (rules + pass/fail) when metrics are computed

## Human governance
- Any change to task definition / policy requires Human Lead sign-off and a new checkpoint review.
