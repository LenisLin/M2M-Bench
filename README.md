# M2M-Bench

M2M-Bench is an audit-grade benchmark for perturbation-response concordance.

## Core Scope

- `Task1`: modality concordance with `perturbation_type` held fixed.
- `Task2`: mechanism concordance between chemical and genetic cohorts inside one
  dataset.
- `Gene` and `Pathway` are the benchmark-wide representation spaces.
- `FM` enters the main manuscript only through the `scPerturb/K562`
  `Figure 3F` local-only panel.
- Figure 1 defines the benchmark.
- Figure 2 carries Task1 main evidence.
- Figure 3 carries Task2 main evidence.
- Task1 and Task2 stay separate in the main text.

## Start Here

- `AGENTS.md`
- `docs/redesign_checkpoint.md`
- `docs/governance/state.md`
- `docs/governance/runbook.md`
- `docs/contracts/task1_spec.md`
- `docs/contracts/task2_spec.md`
- `docs/contracts/output-schemas.md`
- `docs/manuscript_master.md`
- `docs/plotting/plotting_preparation_freeze.md`
- `docs/plotting/manuscript_figure_legends.md`

## Active Roots

- Task1 data: `/mnt/NAS_21T/ProjectData/M2M/data/task1`
- Task2 data: `/mnt/NAS_21T/ProjectData/M2M/data/task2`
- Stage runs: `/mnt/NAS_21T/ProjectData/M2M/runs`
- Manuscript analysis: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Plot review export: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`

## Repo Layout

- `docs/`: benchmark, governance, manuscript, and plotting contracts
- `.agents/skills/`: repo-scoped operating workflows
- `scripts/fm_extractors/`: FM extraction utilities

## Local Checkout

The local checkout is source-only. Audited data, stage outputs, and manuscript
analysis live on NAS-backed roots.
