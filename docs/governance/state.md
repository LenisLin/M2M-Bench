# M2M-Bench Project State (Onboarding SoT)

Last updated: 2026-04-06 (repo/NAS source-only cleanup pass)

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
- Legacy Task2 snapshot v1 and scPerturb-only Task2 runs: historical evidence only, not the corrected authoritative path

## Authoritative benchmark path
- Task1 authoritative audited runs:
  - `/mnt/NAS_21T/ProjectData/M2M/runs/s0_build_data_inventory_0303/s0_build_data_inventory`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics`
- Corrected Task2 authoritative successor path:
  - `/mnt/NAS_21T/ProjectData/M2M/runs/0310_fix_1hae/s3_build_task2_multisource_snapshot`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource`
- S7 authoritative implementation/evidence:
  - `scripts/s7_project_benchmark_synthesis.py`
  - `/mnt/NAS_21T/ProjectData/M2M/runs/s7_project_benchmark_synthesis_implcheck_20260311_a/s7_project_benchmark_synthesis`

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
- manuscript analysis/support split is now explicit:
  - canonical manuscript analysis root: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
  - active manuscript support root: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support`
  - historical manuscript archive: `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history`
- retained reference namespaces remain available:
  - `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready`
  - `/mnt/NAS_21T/ProjectData/M2M/Model_Evaluation_Results`
- Treat NAS-backed absolute paths as the current source of truth for audited evidence discovery.
- The local checkout is now intended to be source-only: code/docs/tests/configs live locally, while snapshots, runs, staging exports, and figure outputs live on NAS.
- R plotting must run inside conda env `Spatial`; repo-local R library/cache relocation is not part of the accepted current operating model.
- Post-migration verification is complete; operational migration, reconciliation, and verification evidence is archived under `/mnt/NAS_21T/ProjectData/M2M/manifests/` and is not part of the active repo-root surface.

## Current stage
- Task1 authoritative through S2
- Corrected Task2 authoritative through S6
- S7 implemented with implcheck evidence
- NAS-backed storage and local path compatibility established
- F3.6–F3.8 formal modifier result objects (time, dose, target multiplicity) materialized 2026-03-24 via `scripts/manuscript_task2_c2g_modifiers.py`; outputs under `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/`
- Main-figure redesign contract frozen for Figure 1–3 and Extended Figure 1–8
- Plot-data pre-aggregation builder updated in `scripts/manuscript_plot_ready_tables.py`
- R plotting namespace refit to clean on-canvas panels with panel letters only
- Legend prose moved out of figures and into `docs/plotting/manuscript_figure_legends.md`
- NAS-backed plot-ready reductions for the redesigned exemplar panels are the active rendering surface
- The active plotting-document surface is now limited to `docs/plotting/plotting_preparation_freeze.md` and `docs/plotting/manuscript_figure_legends.md`; draft/progress snapshots are not part of the live contract surface
- AVCP-era prompt/template entrypoints have been removed from the active root surface in favor of `prompts/M2M_OPERATING_PROMPT.md` and `docs/governance/repo_conventions.md`

## Next tasks
1. Keep only one latest authoritative result surface on NAS for manuscript support and figure export.
2. Re-check historical docs/scripts that still mention repo-local snapshot/result roots before future reuse.
3. Install any missing R dependencies into conda env `Spatial` rather than reviving repo-local library workarounds.

## Non-goals (do not drift)
- No silent fallback from corrected Task2 successor runs to legacy scPerturb-only runs
- No hidden recomputation of Task1/Task2 core metrics in S7 or later reporting layers
- No deletion or relocation of tracked source/docs files as part of storage cleanup
