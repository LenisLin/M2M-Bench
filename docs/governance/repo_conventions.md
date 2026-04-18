# Repository Conventions

## Source Hierarchy

When docs disagree, use this order:

1. Audited manifests and stage outputs
2. `docs/contracts/*.md`
3. `docs/governance/*.md`
4. `docs/manuscript_master.md` and `docs/plotting/*.md`

## Active Roots

- Task1 data: `/mnt/NAS_21T/ProjectData/M2M/data/task1`
- Task2 data: `/mnt/NAS_21T/ProjectData/M2M/data/task2`
- Stage runs: `/mnt/NAS_21T/ProjectData/M2M/runs`
- Manuscript analysis: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Plot review export: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`

## Documentation Rules

- Use the current field names and figure roles from
  `docs/redesign_checkpoint.md`.
- Define every field, metric, stage, and panel when it first appears.
- Remove inactive names instead of carrying multiple names for one object.

## Evidence-First Reporting

- Cite file paths, manifests, tables, or command output.
- If evidence is incomplete, state the missing check directly.

## README Policy

`README.md` is maintained manually as the repo entry document.
