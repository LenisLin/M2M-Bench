# Data Contracts

This file is an index to the active frozen data contracts. It is not a second
schema source.

## Authoritative contract files

- `docs/contracts/task1_spec.md`: Task1 instance tables, split-half rules, cross alignment, retrieval, and frozen constants.
- `docs/contracts/task2_spec.md`: corrected multisource Task2 snapshot, cohort identity, representation scope, retrieval semantics, and readiness rules.
- `docs/contracts/output-schemas.md`: stage-level plot-ready, audit-ready, and S7 project-level output schemas. It is not the manuscript-facing canonical object registry.
- `docs/plotting/plotting_preparation_freeze.md`: frozen panel manifest, plot-ready support outputs, and R plotting namespace plan for the plotting-preparation phase.
- `docs/governance/runbook.md`: stage routing, manifest requirements, and canonical stage outputs.
- `docs/manuscript_master.md`: frozen manuscript-facing figure logic, canonical downstream object list, and current support vs historical boundaries.

## Frozen data roots

- `/mnt/NAS_21T/ProjectData/M2M/snapshots/task1_snapshot_v1/`: the only Task1 snapshot root.
- `/mnt/NAS_21T/ProjectData/M2M/snapshots/task2_snapshot_v1/`: legacy/interim scPerturb-K562 Task2 evidence only.
- `/mnt/NAS_21T/ProjectData/M2M/snapshots/task2_snapshot_v2/`: corrected multisource Task2 successor root.

## Evidence-first validation

- Materialized corrected Task2 snapshot facts must be checked against `/mnt/NAS_21T/ProjectData/M2M/snapshots/task2_snapshot_v2/snapshot_manifest.json`.
- Stage-level output truth must be checked against the audited stage manifest and assertion files under the configured NAS-backed run root. `runs/<run_id>/<stage>/...` remains a stage-contract notation and may be absent from a local checkout.
- If this index conflicts with the contract files above or audited manifests, the contract files and audited manifests win.

## Manuscript-facing alignment notes

- Figure 1 is doc-derived in the current frozen phase. There is no dedicated canonical Figure 1 analysis object to hunt for under `analysis/`.
- The frozen Figure 2/Figure 3 canonical downstream object registry lives in `docs/manuscript_master.md` and `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/framework_analysis_manifest.json`.
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/plot_ready/` is the active support-only downstream plotting layer. It is not a new canonical evidence root and should be interpreted through the freeze note plus the `output-schemas.md` plot-ready section.
- `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/` is the retained historical manuscript root. Files there remain discoverable historical evidence, not current canonical outputs.
- The manuscript-facing biological indexing unit is `(dataset, cell_line, target, perturbation_type)`, but the current corrected Task2 code contract still keys cohorts on `mech_key=(dataset, cell_line, target_token)`. Treat that as a documented layering distinction, not as proof that S3-S6 semantics have already changed.
