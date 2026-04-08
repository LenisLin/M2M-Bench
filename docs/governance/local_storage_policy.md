# Local Storage Policy

Last updated: 2026-04-06

## Goal

Keep the local checkout source-only. Code, docs, tests, prompts, and configs may
live in the repository checkout; benchmark snapshots, runs, staging exports, and
figure outputs must live on NAS only.

## Current storage layout

- NAS snapshot root:
  - `/mnt/NAS_21T/ProjectData/M2M/snapshots`
- NAS active-results root:
  - `/mnt/NAS_21T/ProjectData/M2M/runs`
- NAS archive root:
  - `/mnt/NAS_21T/ProjectData/M2M/archive`
- Manuscript analysis root:
  - `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Manuscript support root:
  - `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support`
- Historical manuscript archive root:
  - `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history`

## Local checkout policy

- The repository checkout must not retain authoritative or staging result roots:
  - `data/` snapshots
  - `runs/`
  - `tmp/`
  - `tmp_r_lib/`
  - locally staged plot-ready or figure-export directories
- If a tool needs a staging root, it must stage under NAS, for example
  `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/...`.
- Local symlink compatibility layers are not part of the accepted operating
  model for the current cleanup state.

## Preservation rule

- “Preserve authoritative artifacts” means:
  - preserve the authoritative data on NAS
  - preserve manifests/audit evidence and deterministic verification records

It does **not** require authoritative run bundles or snapshot mirrors to remain
in the local checkout.

## Verification rule

After any migration or relinking step:

- verify that the local checkout contains only source/documentation surfaces
- verify checksum evidence for required manifest/audit/summary/per-query files
- record temporary machine-readable inventories only in a disposable location
  such as `/tmp` or NAS staging, and delete them after verification completes

## Scope exclusions

- This policy does not authorize deletion or relocation of tracked source files.
- This policy does not authorize cleanup of docs except through explicit documentation updates.

## Reference namespaces

- `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready`
- `/mnt/NAS_21T/ProjectData/M2M/Model_Evaluation_Results`

These namespaces remain available for reference and legacy comparison. They are
not the canonical manuscript evidence root and should not replace audited run
or manuscript-analysis discovery paths.
