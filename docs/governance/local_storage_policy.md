# Local Storage Policy

## Goal

Keep the local checkout source-only. Audited data, runs, plot exports, and
manuscript analysis stay on NAS-backed roots.

## Active Roots

- Task1 data: `/mnt/NAS_21T/ProjectData/M2M/data/task1`
- Task2 data: `/mnt/NAS_21T/ProjectData/M2M/data/task2`
- Stage runs: `/mnt/NAS_21T/ProjectData/M2M/runs`
- Manuscript analysis: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Plot review export: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`

## Local Checkout Rules

- Do not store audited task data in the repo checkout.
- Do not store stage run directories in the repo checkout.
- Do not store plot review exports in the repo checkout.
- Use NAS-backed staging roots for temporary build products.

## Verification Rule

After any data-root or run-root update:

- confirm that the local checkout still contains source and docs only
- confirm that the affected NAS-backed root contains manifests and output tables
- record disposable machine-readable inventories only in `/tmp` or a NAS-backed
  staging root
