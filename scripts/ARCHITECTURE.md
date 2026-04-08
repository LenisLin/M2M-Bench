# Script Architecture

This file is the top-level map for `scripts/s*` and `scripts/manuscript*`.
It distinguishes active canonical pipeline stages, active manuscript canonical
builders, support-only manuscript utilities, and retained historical surfaces.

## 1. Source of Truth and Output Roots

- Authoritative benchmark and manuscript-facing outputs are NAS-backed under `/mnt/NAS_21T/ProjectData/M2M/runs`.
- `runs/<run_id>/<stage>/...` remains the stage-contract notation used by script headers and contracts, but repo-root `runs/` may be absent in a local checkout.
- Figure 1 is doc-derived. There is no dedicated canonical Figure 1 analysis object in the current frozen phase.
- The frozen manuscript-facing registry is `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/framework_analysis_manifest.json`.
- Active support-only manuscript surfaces live under `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/`.
- Historical manuscript-only surfaces live under `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/`.
- `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready` and `/mnt/NAS_21T/ProjectData/M2M/Model_Evaluation_Results` are retained reference namespaces, not the primary manuscript evidence root.

## 2. Active Benchmark Pipeline (S0-S7)

| Script | Status | Manuscript role | Primary outputs | Treat as live canon? |
| --- | --- | --- | --- | --- |
| `scripts/s0_build_data_inventory.py` | active canonical pipeline | upstream inventory only; not direct Figure 2/3 evidence | `task1_data_inventory_long.csv`, `data_source_manifest.csv` | yes |
| `scripts/s1_task1_internal_metrics.py` | active canonical pipeline | primary Task1 internal benchmark evidence feeding Figure 2 | Task1 internal retrieval, leaderboard, attrition, audit bundle | yes |
| `scripts/s2_task1_cross_metrics.py` | active canonical pipeline | primary Task1 cross benchmark evidence feeding Figure 2 | Task1 cross retrieval, leaderboard, alignment proof, attrition, audit bundle | yes |
| `scripts/s3_build_task2_multisource_snapshot.py` | active canonical pipeline | corrected Task2 snapshot builder upstream of Figure 3 | multisource Task2 snapshot, coverage, manifests | yes |
| `scripts/s4_task2_group_concordance_multisource.py` | active canonical pipeline | primary Task2 group benchmark evidence feeding Figure 3 | `task2_group_concordance.csv`, attrition, audit bundle | yes |
| `scripts/s5_task2_retrieval_multisource.py` | active canonical pipeline | primary Task2 direction-specific retrieval evidence feeding Figure 3 | per-query retrieval, summaries, chance checks, attrition | yes |
| `scripts/s6_task2_result_synthesis_multisource.py` | active canonical pipeline | Task2 synthesis layer feeding canonical Figure 3 objects | group/retrieval leaderboards and benchmark summary | yes |
| `scripts/s7_project_benchmark_synthesis.py` | active canonical pipeline | project-level context only; not primary Figure 2/3 evidence | project synthesis tables | yes, but not as primary Figure 2/3 evidence |

## 3. Active Manuscript Canonical Builders

| Script | Status | Manuscript role | Primary outputs | Treat as live canon? |
| --- | --- | --- | --- | --- |
| `scripts/manuscript_framework_analysis_objects.py` | active manuscript canonical builder | materializes the frozen 10-object manuscript backbone and registry | canonical Figure 2/3 objects plus manifest | yes |
| `scripts/manuscript_comparison_statistics.py` | active manuscript canonical builder | downstream comparison-layer statistics over the frozen manuscript objects | `manuscript_comparison_statistics.csv` | yes |

## 4. Support-only Manuscript Builders and Helpers

| Script | Status | Manuscript role | Primary outputs | Treat as live canon? |
| --- | --- | --- | --- | --- |
| `scripts/manuscript_a1_task1_internal_cross.py` | support-only | Figure 2 A1 degradation/support bridge input | `task1_internal_vs_cross_group_bridge*.csv` | no |
| `scripts/manuscript_a3_direction_robustness.py` | support-only | appendix/transparency audit for Task2 directions | `task2_direction_robustness_audit.csv` | no |
| `scripts/manuscript_local_support_tables.py` | support-only | supplementary B5/B8 tables and panel-support ledger refresh | support tables under `docs/` | no |
| `scripts/manuscript_plot_ready_tables.py` | support-only | Figure 1 registry assembly plus live plotting reductions for `2A-2F` and `3B-3F`, alongside `EF4`, `EF6`, and `EF8` support tables and non-main support extracts such as the Task2 scope summary | `plot_ready/figure1/*.csv`, `plot_ready/figure2/*.csv`, `plot_ready/figure3/*.csv`, `plot_ready/extended/*.csv` | no |
| `scripts/manuscript_task1_group_support.py` | helper | shared Task1 group-metric helper for manuscript support builders | no standalone output | no |
| `scripts/manuscript_task2_c2g_modifiers.py` | support-only modifier builder | formal modifier result objects for F3.6–F3.8 (time, dose, target multiplicity); C2G only; non-parametric tests + BH-FDR | `figure3_task2_c2g_time_effect_summary.csv`, `figure3_task2_c2g_time_effect_stats.csv`, `figure3_task2_c2g_dose_effect_summary.csv`, `figure3_task2_c2g_dose_effect_stats.csv`, `figure3_task2_c2g_target_multiplicity_summary.csv`, `figure3_task2_c2g_target_multiplicity_stats.csv`, `figure3_task2_c2g_modifier_comparison_statistics.csv` | no |

## 5. Historical Retained Scripts

| Script | Status | Manuscript role | Primary outputs | Treat as live canon? |
| --- | --- | --- | --- | --- |
| `scripts/s3_build_task2_snapshot.py` | historical retained | legacy scPerturb-K562 snapshot reconstruction | legacy Task2 snapshot v1 artifacts | no |
| `scripts/s4_task2_group_concordance.py` | historical retained | legacy scPerturb-K562 group concordance reconstruction | legacy Task2 group outputs | no |
| `scripts/s5_task2_retrieval.py` | historical retained | legacy scPerturb-K562 retrieval reconstruction | legacy Task2 retrieval outputs | no |
| `scripts/s6_task2_result_synthesis.py` | historical retained | legacy scPerturb-K562 synthesis reconstruction | legacy Task2 synthesis outputs | no |
| `scripts/manuscript_a2_task1_task2_bridge.py` | historical retained | rebuilds the old Task1↔Task2 bridge surface | `task1_task2_group_bridge*.csv` | no |
| `scripts/manuscript_b6_c2g_dose_time_sensitivity.py` | historical retained | rebuilds old B6 C2G dose/time sensitivity support | `task2_c2g_dose_time_sensitivity*.csv` | no |
| `scripts/manuscript_figure_analysis_inputs.py` | historical retained | rebuilds legacy figure-facing overview/backfill files | `figure2_a1_*`, `figure3_a2_*`, `figure3_b6_*` | no |

## 6. Output Directories

- `analysis/`: canonical manuscript-facing objects and retained historical overview files.
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/`: active canonical manuscript analysis root plus current-phase supplementary analysis objects.
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/plot_ready/`: support-only plotting-preparation tables for Figure 1 registries plus main/extended Figure 2/3 panels; downstream of `analysis/` and not a canonical evidence root.
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/group_bridge/`: active support-only Task1 internal-to-cross bridge tables.
- `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/`, `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/group_bridge/`, `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/sensitivity/`, `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/notes/`: retained historical manuscript surfaces.
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/`: appendix/transparency outputs from older support surfaces.
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/b6_dose_time_overlap/task2_dose_time_overlap_audit.csv` is retained appendix evidence with no live `scripts/manuscript*` owner in the current cleanup phase.

R plotting namespace:
- `plotting/R/shared/`: shared IO, theme, labels, palette, orderings, and assembly helpers.
- `plotting/R/figures/`: one script per main or extended figure.
- `plotting/R/config/`: panel labels and palette config.
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/figures/`: active figure export root for the R plotting layer.

## 7. What Is Not Canonical for Figure 2 or Figure 3

- Figure 1 is not backed by a dedicated analysis object.
- Figure 2 does not use Task1↔Task2 bridge language or A2 bridge outputs.
- Figure 2 cell-line/target high-concordance summaries are not canonical objects today.
- Figure 3 keeps `C2G` primary and `G2C` supporting; do not collapse directions.
- Figure 3 uses Task1 internal contextual support only, not Task1 cross.
- `specificity_tier` remains retired as a benchmark-wide modifier surface (no new formal outputs exist).
- Historical B6 dose/time artifacts (`task2_c2g_dose_time_sensitivity*.csv` under `sensitivity/`) and B5 n_targets sensitivity (`docs/task2_c2g_query_n_targets_sensitivity.csv`) remain historical/supplementary-only surfaces.
- F3.6 (time effect), F3.7 (dose effect), and F3.8 (target multiplicity) formal modifier result families have been materialized as supplementary result objects (generated 2026-03-24) by `scripts/manuscript_task2_c2g_modifiers.py`. These are NOT retired placeholders; they are formal manuscript question families, C2G-only and not benchmark-wide.
- Plot-ready support tables and R plotting scripts do not create new benchmark evidence or new canonical Figure 2/Figure 3 analysis objects.
- Historical files may remain on disk; presence on disk does not make them canonical.

## 8. Auxiliary Families

- `scripts/fm_extractors/*` are auxiliary feature-generation utilities outside the frozen Figure 2/3 manuscript backbone.
- They should be documented as upstream utilities, not as manuscript-facing evidence builders.
