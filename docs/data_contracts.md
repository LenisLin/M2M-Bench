# Data Contracts

This file is an index to the active frozen data contracts. It is not a second
schema source.

## Authoritative contract files

- `docs/contracts/task1_spec.md`: Task1 instance tables, split-half rules, cross alignment, retrieval, and frozen constants.
- `docs/contracts/task2_spec.md`: corrected multisource Task2 snapshot, cohort identity, representation scope, retrieval semantics, and readiness rules.
- `docs/contracts/output-schemas.md`: plot-ready, audit-ready, and S7 project-level output schemas.
- `docs/governance/runbook.md`: stage routing, manifest requirements, and canonical stage outputs.

## Frozen data roots

- `data/task1_snapshot_v1/`: the only Task1 snapshot root.
- `data/task2_snapshot_v1/`: legacy/interim scPerturb-K562 Task2 evidence only.
- `data/task2_snapshot_v2/`: corrected multisource Task2 successor root.

## Evidence-first validation

- Materialized corrected Task2 snapshot facts must be checked against `data/task2_snapshot_v2/snapshot_manifest.json`.
- Stage-level output truth must be checked against `runs/<run_id>/<stage>/run_manifest.json` and `audit_assertions.json`.
- If this index conflicts with the contract files above or audited manifests, the contract files and audited manifests win.
