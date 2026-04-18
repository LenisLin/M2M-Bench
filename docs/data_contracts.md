# Data Contracts

This file indexes the current contract surfaces for task data, result objects,
and manuscript analysis.

## Contract Map

- `docs/redesign_checkpoint.md`: benchmark-wide terminology and figure split
- `docs/contracts/task1_spec.md`: Task1 units, retrieval rules, and metric scope
- `docs/contracts/task2_spec.md`: Task2 units, directions, and cohort rules
- `docs/contracts/output-schemas.md`: stage output tables and manifest payloads
- `docs/manuscript_master.md`: main manuscript structure and figure roles
- `docs/plotting/plotting_preparation_freeze.md`: plot-ready assembly and
  rendering boundary

## Current Naming

- `Task1` unit identity uses `perturbation_gene`.
- `Task2` unit identity uses `anchor_gene`.
- Retrieval rows use `query_instance_id`.
- Pattern tables use `pair_mean_enrichment`.
- Aggregated Task2 tables stay keyed by `anchor_gene`.
- Row-level identity stays in `perturbation_gene`.

## Active Roots

- Task1 data: `/mnt/NAS_21T/ProjectData/M2M/data/task1`
- Task2 data: `/mnt/NAS_21T/ProjectData/M2M/data/task2`
- Stage runs: `/mnt/NAS_21T/ProjectData/M2M/runs`
- Manuscript analysis: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`

## Validation Rule

- Stage tables are validated by the stage manifest, audit assertions, and
  per-table manifest stored under the stage run directory.
- Task data facts are validated against manifests stored inside the active task
  roots.
- If this file disagrees with audited manifests or contract docs, the audited
  manifests and contract docs win.
