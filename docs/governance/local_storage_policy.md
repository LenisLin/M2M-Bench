# Local Storage Policy

Last updated: 2026-03-12

## Goal

Keep benchmark outputs off local disk when practical while preserving stable local project paths for scripts, audits, and manual inspection.

## Current storage layout

- NAS active-results root:
  - `/mnt/NAS_21T/ProjectData/M2M/runs`
- NAS archive root:
  - `/mnt/NAS_21T/ProjectData/M2M/archive/runs`
- NAS reports archive root:
  - `/mnt/NAS_21T/ProjectData/M2M/archive/reports`

## Current local compatibility strategy

- `runs/` remains a local directory.
- Each migrated authoritative or archive run-id directory is represented locally by a run-id-level symlink under `runs/`.
- `reports` may be represented locally by a root symlink to the NAS archive location.

This policy keeps script-visible paths stable:

- `runs/<run_id>/<stage>/...`
- `reports/...`

## Accepted `runs/` model

- Root-level `runs -> /mnt/NAS_21T/ProjectData/M2M/runs` normalization is not the chosen current strategy.
- The accepted operating model is:
  - local `runs/` remains a real directory
  - each local `runs/<run_id>` entry point is a symlink to either:
    - `/mnt/NAS_21T/ProjectData/M2M/runs/<run_id>` for active authoritative results, or
    - `/mnt/NAS_21T/ProjectData/M2M/archive/runs/<run_id>` for archived historical results
- This mixed active-plus-archive namespace is intentional and stable for current production use.

Any future namespace-unification project would require separate approval and must preserve the existing local `runs/<run_id>` compatibility contract before changing the root model.

## Preservation rule

- “Preserve authoritative artifacts” means:
  - preserve the authoritative data on NAS
  - preserve local path compatibility inside the repository
  - preserve manifests/audit evidence and deterministic verification records

It does **not** require authoritative run bundles to remain on local disk as regular directories.

## Verification rule

After any migration or relinking step:

- verify local-path resolution for authoritative runs
- verify checksum evidence for required manifest/audit/summary/per-query files
- record the result in explicit machine-readable artifacts

## Scope exclusions

- This policy does not authorize deletion or relocation of tracked source files.
- This policy does not authorize cleanup of docs except through explicit documentation updates.
- This policy does not authorize mutation of authoritative data on NAS.
