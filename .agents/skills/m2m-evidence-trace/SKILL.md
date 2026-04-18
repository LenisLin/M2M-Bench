---
name: m2m-evidence-trace
description: Trace M2M claims back to authoritative runs, manifests, tables, and NAS-backed evidence roots.
---

# M2M Evidence Trace

Use this skill when interpreting results, validating figure inputs, checking
manifests, or answering "what evidence supports this claim?"

## Load Order

Read only what the task needs:

1. `docs/governance/repo_conventions.md`
2. `docs/governance/state.md`
3. `docs/governance/local_storage_policy.md`
4. the relevant task contract in `docs/contracts/`
5. plotting docs for figure-facing tasks

## Evidence Ranking

Use this order:

1. audited manifests and stage outputs
2. NAS-backed run directories and manuscript analysis roots
3. contract docs that define output meaning
4. governance docs that define storage and reporting rules

## Core Workflow

1. Name the exact claim, figure panel, table, or path question.
2. Locate the authoritative NAS-backed root first.
3. Find the manifest, audit assertions, CSV, parquet table, or run directory
   that answers the question directly.
4. Report the evidence with explicit file paths and the interpretation boundary.
5. If evidence is missing, state the next check instead of guessing.

## Guardrails

- Do not treat repo-local staging outputs as authoritative evidence.
- Prefer path plus manifest plus table triads over prose-only explanations.
