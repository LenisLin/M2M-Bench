# Manuscript Evidence Extraction Plan

## 1. Objective of the next round

The next round should extract claim-ready evidence pointers and denominator context for the frozen manuscript story without changing the story itself. The goal is to assemble a contract-consistent claim-to-evidence map for Figure 2 and Figure 3 using only authoritative local outputs and frozen snapshot metadata.

Execution order note: the next round should first build a source register, and only then populate the claim-to-evidence map. The source register should enumerate authoritative files, stage ownership, key metric fields, denominator/support fields, scope restrictions, and join dependencies before any claim row is marked extractable.

Follow-on note: any new manuscript-support bridge, robustness, or explanation tables should be staged according to `docs/manuscript_analysis_expansion_plan.md` before they are added to the claim map.

## 2. What will not be done yet

- no new evidence claims or manuscript prose beyond extraction notes
- literature positioning may provide framing context, but authoritative local stage outputs remain the only source of manuscript evidence
- no reruns, recomputation, or contract changes
- no plotting before the claim-to-evidence map is complete
- no use of legacy scPerturb-only Task2 outputs as authoritative corrected Task2 evidence
- no introduction of significance, benchmark-health, or later-stage helper outputs that are not part of the locked Task1/Task2 evidence path
- no over-reading of FM behavior outside the supported scPerturb K562 scope
- no use of S7 to override or replace authoritative Task1 or corrected Task2 stage outputs during claim extraction

## 3. Claim-to-evidence map schema

Use one row per planned claim or figure-panel statement.

The claim-to-evidence map should be populated only after the source register is complete for the relevant Task1 or Task2 slice.

| field | required content |
| --- | --- |
| `claim_id` | Stable identifier such as `C1`, `C2`, `C3`, `C4`, or a scoped subclaim ID. |
| `figure_panel` | Planned manuscript panel target such as `F2a`, `F2b`, `F3c`. |
| `task` | `Task1_internal`, `Task1_cross`, `Task2_group`, `Task2_retrieval`, or `Task2_FM_scope`. |
| `metric` | Contract-consistent metric name from the authoritative source table. |
| `representation` | Representation or representation family carried by the source table. |
| `stratification` | Dataset, direction, perturbation type, cell line, and any approved metadata slice used for interpretation. |
| `source_stage` | Authoritative stage label such as `s1_task1_internal_metrics`, `s2_task1_cross_metrics`, `s4_task2_group_concordance_multisource`, `s5_task2_retrieval_multisource`, `s6_task2_result_synthesis_multisource`, or `task2_snapshot_v2`. |
| `source_file` | Concrete repo path to the authoritative table or metadata file. |
| `denominator_fields` | The exact support columns that must be carried with the claim. |
| `status` | `pending`, `extracted`, `needs_join`, `blocked`, or `rejected_for_main_text`. |
| `notes` | Join instructions, caution codes, scope restrictions, or limitation wording tied to the row. |

Recommended additional columns:

- `source_run_id`
- `source_table_key`
- `caution_codes`
- `support_decision`

## 4. Priority extraction package A: Figure 2 / Task1

Purpose: establish the Task1 reference frame without turning Task1 into a ceiling claim.

Primary claim coverage:

- Claim 1
- Claim 2

Authoritative source files:

- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv`
- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_summary.csv`
- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_attrition.csv`
- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_chance_identity_check.csv`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_leaderboard_long.csv`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_summary.csv`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_attrition.csv`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_alignment_proof.csv`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_chance_identity_check.csv`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_group_cross.parquet` only if a per-group support check is needed
- `docs/contracts/exclusions-and-policies.md` for the frozen cross-chemical exclusion language

Extraction tasks:

- extract Task1 internal group-level and retrieval rows by `dataset_or_direction`, `perturbation_type`, `representation`, and `metric_name`
- extract Task1 cross rows separately by direction and preserve the current genetic-only scope
- record denominator transparency from `n_total`, `n_valid`, `n_excluded_missing_metric_or_mpos0`, `N_gallery_max`, and any metric-specific support counts carried in the source table
- record cross-chemical exclusion as a policy note from the contract, not as a null empirical result
- use chance-identity outputs only as validation support for stable non-random retrieval interpretation

Figure intent to preserve:

- group-level evidence is primary
- retrieval is supporting
- internal versus cross is a calibration comparison, not a formal upper-bound argument

## 5. Priority extraction package B: Figure 3 / Task2

Purpose: extract the main Task2 evidence path for mechanism concordance and its dependence on biological structure.

Primary claim coverage:

- Claim 1
- Claim 3
- Claim 4

Authoritative source files:

- `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv`
- `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_attrition.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary_long.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_chance_identity_check.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet` when query-level stratification is needed
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_concordance_long.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_benchmark_summary_long.csv`
- `data/task2_snapshot_v2/task2_pairs_coverage.csv`
- `data/task2_snapshot_v2/representation_availability_registry.csv`
- `data/task2_snapshot_v2/lincs/derived/delta_meta.csv`
- `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`

Extraction tasks:

- extract group-level Task2 evidence by `dataset`, `cell_line`, `representation`, and metric while preserving metric-specific denominator fields such as `n_targets_total` and `n_targets_metric_valid`
- extract bidirectional retrieval evidence by `dataset`, `cell_line`, `direction`, and `representation`, carrying `gallery_definition_id`, `pos_definition_id`, `n_total`, `n_valid`, `N_gallery_mean`, `N_gallery_max`, and `m_pos` support fields
- use `task2_pairs_coverage.csv` and `task2_group_attrition.csv` to document eligible versus ineligible mechanism cohorts and metric-level support loss
- use `representation_availability_registry.csv` to separate scope policy from attrition
- restrict FM evidence rows to `dataset=scPerturb`, `cell_line=K562`, and approved FM representations only
- join retrieval per-query rows to frozen snapshot metadata when `specificity_tier` or query-level target-complexity stratification is required
- derive `n_targets` from `target_tokens` only where it is not already materialized, and mark that derivation explicitly as metadata extraction rather than new analysis

Specific handling for Figure 3 stratification:

- `dataset` and `cell_line` should come from S4-S6 outputs directly
- `n_targets` can be taken from `query_n_targets` in Task2 per-query retrieval outputs or derived from `target_tokens` in snapshot `delta_meta.csv`
- `specificity_tier` is retained in snapshot `delta_meta.csv` but is not carried as a cohort-level field in S4-S6 summaries; any Figure 3 tier analysis must document whether it is query-level, row-level membership-level, or an approved cohort aggregation

## 6. Decision rules before plotting

- no panel is plot-eligible until its claim row has a concrete `source_file`, explicit denominator fields, and a status other than `pending`
- no claim row is extractable until its source register entry records authoritative file ownership, key metric fields, denominator/support fields, scope restrictions, and join dependencies
- use S6 Task2 synthesis tables for summary extraction, but resolve ambiguities against S4/S5 or frozen snapshot metadata rather than against S7 or ad hoc recomputation
- use S1/S2 Task1 summary tables as the default source; do not recompute Task1 evidence from per-query outputs unless the summary tables are insufficient for denominator reporting
- S7 may be used only as implementation/context evidence and must not override or replace authoritative Task1 or corrected Task2 stage outputs in claim extraction
- keep group-level concordance as the primary evidence layer in both Figure 2 and Figure 3
- keep Task1 cross genetic-only and record the cross-chemical exclusion as contract policy
- keep Task2 retrieval directions separate; do not collapse `C2G` and `G2C`
- do not raw-pool LINCS and scPerturb in Task2 core evidence
- treat high-specificity subsets as sensitivity analyses unless the extraction sheet shows enough support and stability to justify promotion
- do not compare raw similarity magnitudes across representations as if they were on a common effect-size scale
- keep `gallery_definition_id` and `pos_definition_id` with every Task2 retrieval claim so the corrected interpretation cannot drift
- if a claim depends on a stratification key whose aggregation rule is not frozen locally, mark it `blocked` or `needs_join` rather than inferring the missing rule

## 7. Expected outputs of the next round

- a claim-to-evidence map covering all planned Figure 2 and Figure 3 statements
- a source register that lists authoritative files, stage ownership, key metric fields, denominator/support fields, scope restrictions, and join dependencies for every retained claim row
- a short blocker memo for any unresolved join or aggregation issues, especially around `specificity_tier`
- a panel-readiness note that separates main-text-eligible slices from supplementary-only or sensitivity-only slices
- no manuscript draft text and no final plotting artifacts yet
