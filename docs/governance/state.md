# M2M-Bench Project State (Onboarding SoT)

Last updated: 2026-03-13 (post-migration root-ops cleanup)

## Project objective
Deliver an audit-grade benchmark for:
- **Task1 (Modality):** internal reliability + cross-modality comparability
- **Task2 (Mechanism):** Chemical↔Genetic consistency within source across corrected multisource Task2 cohorts

FM analyses remain representation-space complements rather than a third primary task.

## Current analysis scope (locked)
- Task1 internal: LINCS + scPerturb eligible cohorts
- Task1 cross: matched LINCS↔scPerturb subset only
- Corrected Task2 core metrics: multisource `LINCS` + `scPerturb`, stratified by `dataset` and `cell_line`
- Task2 FM scope: scPerturb K562 contract subset only
- Legacy `data/task2_snapshot_v1/` and scPerturb-only Task2 runs: historical evidence only, not the corrected authoritative path

## Authoritative benchmark path
- Task1 authoritative audited runs:
  - `runs/s0_build_data_inventory_0303/s0_build_data_inventory`
  - `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics`
  - `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics`
- Corrected Task2 authoritative successor path:
  - `runs/0310_fix_1hae/s3_build_task2_multisource_snapshot`
  - `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource`
  - `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource`
  - `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource`
- S7 authoritative implementation/evidence:
  - `scripts/s7_project_benchmark_synthesis.py`
  - `runs/s7_project_benchmark_synthesis_implcheck_20260311_a/s7_project_benchmark_synthesis`

## Primary active paths
- Execution runbook: `docs/governance/runbook.md`
- Task2 contract: `docs/contracts/task2_spec.md`
- Project positioning: `docs/contracts/project-positioning.md`
- Local storage policy: `docs/governance/local_storage_policy.md`
- Active stage scripts:
  - `scripts/s0_build_data_inventory.py`
  - `scripts/s1_task1_internal_metrics.py`
  - `scripts/s2_task1_cross_metrics.py`
  - `scripts/s3_build_task2_multisource_snapshot.py`
  - `scripts/s4_task2_group_concordance_multisource.py`
  - `scripts/s5_task2_retrieval_multisource.py`
  - `scripts/s6_task2_result_synthesis_multisource.py`
  - `scripts/s7_project_benchmark_synthesis.py`

## Current engineering status
- Corrected multisource Task2 is implemented through S6 and has an implemented S7 project-level synthesis step.
- Authoritative run artifacts are NAS-backed:
  - active results root: `/mnt/NAS_21T/ProjectData/M2M/runs`
  - archive root: `/mnt/NAS_21T/ProjectData/M2M/archive/runs`
  - reports root: `/mnt/NAS_21T/ProjectData/M2M/archive/reports`
- Local path compatibility is preserved through symlinks:
  - local `runs/` remains a real directory
  - `runs/<run_id>` resolves via run-id-level symlinks spanning both NAS active and NAS archive namespaces
  - `reports` resolves via a root symlink
- Root-level `runs -> /mnt/NAS_21T/ProjectData/M2M/runs` normalization is not part of the accepted current operating model.
- Post-migration verification is complete; operational migration, reconciliation, and verification evidence is archived under `/mnt/NAS_21T/ProjectData/M2M/manifests/` and is not part of the active repo-root surface.

## Current stage
- Task1 authoritative through S2
- Corrected Task2 authoritative through S6
- S7 implemented with implcheck evidence
- NAS-backed storage and local path compatibility established

## Next tasks
1. Keep future generated benchmark outputs on the NAS-backed storage path while preserving local `runs/<run_id>` compatibility semantics.
2. Continue later-stage work, if approved, as new stages rather than reviving retired historical later-stage scripts.
3. Treat any future namespace-unification work as an optional separate project, not as pending migration debt.

## Non-goals (do not drift)
- No silent fallback from corrected Task2 successor runs to legacy scPerturb-only runs
- No hidden recomputation of Task1/Task2 core metrics in S7 or later reporting layers
- No deletion or relocation of tracked source/docs files as part of storage cleanup
