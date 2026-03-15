# Manuscript Support Table Specs

## 1. Purpose

This document defines the planned manuscript-only support tables referenced by `docs/manuscript_analysis_expansion_plan.md`. These outputs are downstream manuscript artifacts. They are not new benchmark-stage contracts and they must not replace authoritative S1, S2, S4, S5, or S6 outputs.

## 2. Conventions

- Carry explicit source provenance with `source_stage`, `source_run_id`, or `source_file` where practical.
- Keep `dataset`, `cell_line`, `direction`, and `representation` explicit whenever they are part of the interpretive scope.
- Carry support columns such as `n_pairs`, `n_targets`, `n_queries`, `n_valid`, `N_gallery_*`, and `m_pos_*` whenever they define interpretability.
- Add `local_only_bool` for scPerturb K562-only outputs and `conditional_bool` for overlap-gated outputs.
- For enrichment outputs, require `background_definition`, `background_n`, and a background registry sourced from the matched and evaluable target pool.
- Do not encode benchmark-wide interpretations in table titles or schema labels for local-only outputs.

## 3. Table Specs

| table_name | layer | primary keys | minimum required columns | primary source inputs | notes |
| --- | --- | --- | --- | --- | --- |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a3_direction_robustness/task2_direction_robustness_audit.csv` | high-risk gate | `dataset,cell_line,direction,representation,robustness_slice` | `mean_mrr_corrected`, `n_queries`, `N_gallery_mean`, `N_gallery_max`, `m_pos_mean`, `null_baseline_type`, `reviewer_attack_vector` | S5 summary + per-query, S6 retrieval leaderboard | Keeps corrected retrieval framework and documents whether asymmetry survives support-aware slicing. |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/b6_dose_time_overlap/task2_dose_time_overlap_audit.csv` | high-risk gate | `dataset,cell_line,direction,overlap_slice` | `n_queries`, `n_overlap_queries`, `time_value`, `dose_value`, `overlap_status`, `local_only_bool` | S5 per-query, snapshot `delta_meta.csv` | If overlap is inadequate, this table plus a limitation statement is the stopping output. |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv` | main bridge | `representation,metric_name,paired_slice_id` | `internal_value`, `cross_value`, `delta_value`, `n_pairs`, `pair_unit_type`, `support_scope_note` | S1 summary/per-query, S2 summary/per-query/group cross | No ceiling language; bridge is limited to common-scope genetic rows. |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv` | main bridge | `dataset,cell_line,target_token,representation,metric_family` | `task1_value`, `task2_value`, `delta_value`, `pair_eligible_bool`, `support_counts`, `scope_note` | S1 inputs, S4 target-level group rows, S6 summaries | Shared-slice external-validity bridge only. |
| `docs/manuscript_core_paired_statistics.csv` | statistical packaging | `claim_family,comparison_id,metric_name` | `test_name`, `n_pairs`, `median_delta`, `p_value`, `fdr_value`, `pair_unit_type`, `source_support_table` | `A1-A3` support tables plus S4/S5 detailed rows | This is the manuscript statistics package, not a benchmark summary table. |
| `docs/task2_n_targets_pairwise_stats.csv` | explanation layer | `dataset,cell_line,representation,query_n_targets,metric_name` | `n_queries`, `comparison_group`, `median_delta`, `p_value`, `fdr_value`, `direction` | Existing `task2_c2g_query_n_targets_sensitivity.csv`, S5 per-query | Keep C2G-only and exact `query_n_targets` values. |
| `docs/task2_dataset_cellline_composition_analysis.csv` | explanation layer | `dataset,cell_line,representation` | `n_targets_eligible`, `n_queries`, `n_attrition_rows`, `availability_status`, `composition_note` | coverage, representation registry, S4/S5 attrition, S6 summaries | Composition context should be read before or alongside enrichment. |
| `docs/task2_target_pool_background_registry.csv` | explanation layer | `dataset,cell_line,representation,background_definition` | `background_n`, `target_pool_source`, `inclusion_rule`, `composition_context_file` | S4 target-level rows, S6 summaries, `task2_pairs_coverage.csv` | Required precursor for enrichment; never whole genome. |
| `docs/task2_target_concordance_rankings.csv` | explanation layer | `dataset,cell_line,representation,target_token` | `concordance_score`, `rank_order`, `high_low_set`, `background_definition`, `composition_context_note` | S4 target-level rows, background registry | Serves as the ranked target universe for enrichment. |
| `docs/task2_high_low_concordance_enrichment.csv` | explanation layer | `dataset,cell_line,representation,set_label,annotation_namespace,term_id` | `background_n`, `set_n`, `effect_size`, `p_value`, `fdr_value`, `background_definition` | target rankings, background registry, annotation assets | Must be interpreted with the composition table and matched-background registry. |
| `docs/task2_dose_time_matched_analysis.csv` | explanation layer, conditional | `dataset,cell_line,direction,representation,matched_slice_id` | `n_queries`, `match_rule`, `mean_mrr_corrected`, `delta_vs_unmatched`, `support_ok_bool` | dose-time overlap audit, S5 per-query, S6 retrieval | Only materialize if overlap is real and support survives stratification. |
| `docs/task2_specificity_local_stats.csv` | reviewer-defense, local-only | `representation,specificity_tier,metric_name` | `n_queries`, `mean_mrr_corrected`, `delta_value`, `p_value`, `fdr_value`, `local_only_bool` | existing specificity sensitivity table, S5 per-query, scPerturb metadata | Local-only and not required for the main-text bridge. |
| `docs/task2_fm_k562_local_comparison_summary.csv` | representation layer, local-only | `representation,metric_family,direction` | `mean_value`, `rank_value`, `n_targets_or_queries`, `local_only_bool`, `support_note` | S6 local FM summaries, panel ledger, availability registry | Candidate local main-text comparison within supported scope only. |
| `docs/task2_fm_k562_per_target_robustness.csv` | representation layer, local-only | `representation,target_token,metric_name` | `metric_value`, `paired_reference_representation`, `delta_value`, `evaluable_bool` | S4 target-level rows, S5 local per-query | Per-target robustness layer for FM interpretation. |
| `docs/task2_fm_k562_paired_tests.csv` | representation layer, local-only | `comparison_id,metric_name` | `n_pairs`, `median_delta`, `p_value`, `fdr_value`, `pair_unit_type`, `local_only_bool` | FM per-target robustness table, S5 local rows | Required before any stronger reading of the FM summary table. |
| `docs/manuscript_case_shortlist.csv` | case layer | `case_id,case_class` | `dataset`, `cell_line`, `target_token`, `representation`, `support_tables_used`, `case_rationale`, `local_only_bool` | bridge tables, stats table, S4/S5 rows, snapshot metadata | Case layer starts only after bridge and statistical packaging are complete. |

## 4. Conditional and Local-Only Rules

- `docs/task2_dose_time_matched_analysis.csv` is conditional on `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/b6_dose_time_overlap/task2_dose_time_overlap_audit.csv`.
- `docs/task2_specificity_local_stats.csv`, `docs/task2_fm_k562_local_comparison_summary.csv`, `docs/task2_fm_k562_per_target_robustness.csv`, and `docs/task2_fm_k562_paired_tests.csv` are local-only.
- `docs/manuscript_case_shortlist.csv` is illustrative only and must not substitute for `A1`, `A2`, or `A4`.

## 5. Last-Resort Upstream Note

If existing S1 outputs prove insufficient for `A1`, `A2`, or `A4`, the only allowed upstream touch in this phase is a conditional audited export of already-computed Task1 internal group rows. That export remains last-resort and should not happen by default.
