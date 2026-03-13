## Current Progress Update

Last updated: 2026-03-13

### Authoritative benchmark status

- **Task1 S0 complete and audited**
  - Run: `runs/s0_build_data_inventory_0303/s0_build_data_inventory`
  - Evidence: `run_manifest.json`, `audit_assertions.json`, `task1_data_inventory_long.csv`, `data_source_manifest.csv`

- **Task1 S1 complete and audited**
  - Run: `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics`
  - Status: hard assertions pass; one non-blocking soft audit exception for FM policy presence

- **Task1 S2 complete and audited**
  - Run: `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics`
  - Status: assertions pass; cross-chemical exclusion remains an explicit contract outcome

- **Corrected multisource Task2 S3-S6 complete and audited**
  - S3: `runs/0310_fix_1hae/s3_build_task2_multisource_snapshot`
  - S4: `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource`
  - S5: `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource`
  - S6: `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource`

- **S7 implemented and evidenced**
  - Script: `scripts/s7_project_benchmark_synthesis.py`
  - Implcheck run: `runs/s7_project_benchmark_synthesis_implcheck_20260311_a/s7_project_benchmark_synthesis`

### Storage migration status

- Authoritative active results have been migrated to NAS-backed storage under:
  - `/mnt/NAS_21T/ProjectData/M2M/runs`
- Historical run bundles have been archived under:
  - `/mnt/NAS_21T/ProjectData/M2M/archive/runs`
- Reports have been archived under:
  - `/mnt/NAS_21T/ProjectData/M2M/archive/reports`
- Local path compatibility is currently provided by:
  - a real local `runs/` directory with run-id-level symlinks under `runs/`
  - a root `reports` symlink

### Post-migration verification

- Post-migration verification is complete for both NAS-backed authoritative runs and archive-backed historical runs.
- Operational execution, reconciliation, and verification records are archived under:
  - `/mnt/NAS_21T/ProjectData/M2M/manifests/storage_migration_manifest.json`
  - `/mnt/NAS_21T/ProjectData/M2M/manifests/local_repo_housekeeping/`
- The active repo surface no longer depends on repo-root migration, reconciliation, or verification artifacts.

### Accepted storage strategy

- Root-level `runs -> /mnt/NAS_21T/ProjectData/M2M/runs` normalization is not the chosen current strategy.
- The accepted production model is:
  - keep local `runs/` as a real directory
  - preserve local `runs/<run_id>` compatibility with run-id-level symlinks to NAS active or NAS archive targets as appropriate

### Active interpretation rule

- “Preserve authoritative artifacts” now means:
  - preserve the authoritative data on NAS
  - preserve local project-path usability through symlink-compatible paths
  - do not require authoritative run bundles to remain on local disk as regular directories
