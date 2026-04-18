---
name: m2m-migration-runbook
description: Apply M2M task-data, path, and run-root rules before touching task-root or cutover workflows.
---

# M2M Migration Runbook

Use this skill when the task touches task data roots, run roots, path policy,
or cutover workflows.

## Load Order

Read these sources first:

1. `docs/governance/runbook.md`
2. `docs/governance/state.md`
3. `docs/governance/local_storage_policy.md`
4. `docs/governance/repo_conventions.md`
5. `docs/data_contracts.md` and `docs/contracts/output-schemas.md` when output
   tables or manifests are affected

## Core Workflow

1. Confirm the active data root and run root.
2. Confirm which stage outputs or manifests the task depends on.
3. Keep staging and cutover work on NAS-backed roots.
4. Keep manifest bundles synchronized with output tables.
5. Verify that the local checkout remains source-only.

## Guardrails

- Do not move audited data into the repo checkout.
- Do not remove tracked source or docs during storage cleanup.
- Do not change benchmark semantics while changing roots or cutover behavior.
