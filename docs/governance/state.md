# M2M-Bench Project State

Last updated: 2026-04-18

## Objective

Deliver an audit-grade benchmark for:

- `Task1`: modality concordance with `perturbation_type` held fixed
- `Task2`: mechanism concordance between chemical and genetic cohorts inside
  one dataset

## Active Scope

- Task1 internal: `LINCS` and `scPerturb`
- Task1 cross: matched genetic slice between `LINCS` and `scPerturb`
- Task2 core metrics: `LINCS` and `scPerturb` within each `dataset` and
  `cell_line`
- `FM` enters manuscript-facing scope only through the `scPerturb/K562`
  `Figure 3F` local-only panel
- Task1 and Task2 stay separate in manuscript-facing docs

## Active Stage Status

- `S0` to `S2`: Task1 audited stage chain
- `S3` to `S6`: Task2 audited stage chain
- `S7`: project synthesis assembled from audited Task1 and Task2 outputs

## Active Roots

- Task1 data: `/mnt/NAS_21T/ProjectData/M2M/data/task1`
- Task2 data: `/mnt/NAS_21T/ProjectData/M2M/data/task2`
- Stage runs: `/mnt/NAS_21T/ProjectData/M2M/runs`
- Manuscript analysis: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Plot review export: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`

## Active Doc Surface

- `docs/redesign_checkpoint.md`
- `docs/contracts/task1_spec.md`
- `docs/contracts/task2_spec.md`
- `docs/contracts/output-schemas.md`
- `docs/manuscript_master.md`
- `docs/plotting/plotting_preparation_freeze.md`
- `docs/plotting/manuscript_figure_legends.md`

The local checkout holds source and docs only.
