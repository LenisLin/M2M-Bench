# Manuscript Source Register

This register indexes the frozen local evidence path for the manuscript extraction pass.

Planning note: any new manuscript-only support tables created for bridge, robustness, or explanation work should be registered downstream of the authoritative benchmark outputs and tracked against `docs/manuscript_analysis_expansion_plan.md`.

Included as authoritative for manuscript evidence:
- Task1 S1 internal outputs
- Task1 S2 cross outputs
- corrected Task2 S4, S5, and S6 outputs
- frozen `data/task2_snapshot_v2/` metadata needed for coverage, scope, and join context
- `docs/contracts/exclusions-and-policies.md` only where frozen policy language is required

Explicitly excluded as non-authoritative for Figure 2 and Figure 3 claims:
- S7 project synthesis as a primary evidence source
- legacy `data/task2_snapshot_v1/` and legacy scPerturb-only Task2 stages as corrected Task2 claim evidence

## Task1 S1 Internal Metrics

### `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv`

- Repo path: `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv`
- Stage ownership: Task1 S1 internal metrics; run id `s1_task1_internal_metrics_0303`
- Evidence role: primary Task1 internal summary table for Figure 2 group-level calibration rows; also carries corrected retrieval rows in the same long table
- Main metrics present: `mean_cosine_centroid`; `mean_pcc_centroid`; `mean_edist_biascorr`; `mean_mrr_corrected`; `mean_hit1_corrected`; `mean_hit5_corrected`; `mean_hit10_corrected`
- Denominator/support fields present: `n_total`; `n_valid`; `n_excluded`; `N_gallery_max`; `cross_alignment_contract`
- Scope restrictions: `scope=internal`; `dataset_or_direction` materializes `LINCS` and `scPerturb`; LINCS rows are `Gene` and `Pathway`; scPerturb rows include local `Gene`, `Pathway`, and FM representations, so common-scope scPerturb internal rows exist locally but must be kept separate from the extended FM-only view
- Allowed claim types: Figure 2 internal reference-frame statements; internal group-level calibration; internal retrieval support after `metric_name` filtering
- Required joins, if any: none for summary extraction; use `task1_attrition.csv` for support loss reasons when denominator interpretation needs explanation
- Cautions or interpretation limits: do not treat Task1 as a ceiling for Task2; do not compare raw similarity magnitudes across representations as a common effect-size scale; keep scPerturb `Gene` and `Pathway` rows distinct from the additional FM rows when using this table for common-scope comparison

### `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_summary.csv`

- Repo path: `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_summary.csv`
- Stage ownership: Task1 S1 internal metrics; run id `s1_task1_internal_metrics_0303`
- Evidence role: authoritative internal retrieval summary for Figure 2 support rows
- Main metrics present: `mean_mrr_raw`; `mean_expected_mrr_chance`; `mean_mrr_corrected`; `mean_hit1_raw`; `mean_expected_hit1_chance`; `mean_hit1_corrected`; `mean_hit5_raw`; `mean_expected_hit5_chance`; `mean_hit5_corrected`; `mean_hit10_raw`; `mean_expected_hit10_chance`; `mean_hit10_corrected`
- Denominator/support fields present: `n_total`; `n_valid`; `n_excluded_missing_metric_or_mpos0`; `N_gallery_max`
- Scope restrictions: `scope=internal` only; scPerturb materializes both common-scope `Gene` and `Pathway` rows and additional FM rows, so common-scope and FM-local retrieval slices should be indexed separately
- Allowed claim types: supporting retrieval interpretation for Task1 internal slices; non-random retrieval support after chance correction
- Required joins, if any: none for metric extraction; pair with `task1_chance_identity_check.csv` for corrected-metric validation
- Cautions or interpretation limits: retrieval is supporting evidence only; keep raw, expected, and corrected metrics distinct; do not let the additional scPerturb FM rows replace the local common-scope `Gene` and `Pathway` slice

### `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_attrition.csv`

- Repo path: `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_attrition.csv`
- Stage ownership: Task1 S1 internal metrics; run id `s1_task1_internal_metrics_0303`
- Evidence role: denominator transparency and exclusion-reason support for Task1 internal
- Main metrics present: none; this is an attrition table keyed by `reason`
- Denominator/support fields present: `reason`; `n_dropped`; `n_total_before`; `notes`
- Scope restrictions: internal-only support table; reasons are metric-path specific and should not be conflated across group-level and retrieval uses
- Allowed claim types: denominator disclosure; exclusion-language support; supplementary support notes
- Required joins, if any: join on `scope`; `dataset_or_direction`; `perturbation_type`; `representation`
- Cautions or interpretation limits: absence of a reason row is not itself a metric result; current rows mix LOO retrieval exclusions and split-half group exclusions

### `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_chance_identity_check.csv`

- Repo path: `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_chance_identity_check.csv`
- Stage ownership: Task1 S1 internal metrics; run id `s1_task1_internal_metrics_0303`
- Evidence role: validation-only support that corrected retrieval metrics equal raw minus expected chance within tolerance
- Main metrics present: `delta_mrr`; `abs_delta_mrr`; `delta_hit1`; `abs_delta_hit1`; `delta_hit5`; `abs_delta_hit5`; `delta_hit10`; `abs_delta_hit10`
- Denominator/support fields present: keyed by `scope`; `dataset_or_direction`; `perturbation_type`; `representation`
- Scope restrictions: internal retrieval validation only
- Allowed claim types: validation support for stable non-random retrieval interpretation; not a primary claim table
- Required joins, if any: join to `task1_retrieval_summary.csv` on the shared keys when validation provenance needs to be carried
- Cautions or interpretation limits: use only as a correctness check; do not report these deltas as scientific effect sizes

## Task1 S2 Cross Metrics

### `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_leaderboard_long.csv`

- Repo path: `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_leaderboard_long.csv`
- Stage ownership: Task1 S2 cross metrics; run id `s2_task1_cross_metrics_0303`
- Evidence role: primary Task1 cross summary for Figure 2 cross-calibration rows
- Main metrics present: `mean_cosine_centroid`; `mean_pcc_centroid`; `mean_edist_biascorr`; `mean_mrr_corrected`; `mean_hit1_corrected`; `mean_hit5_corrected`; `mean_hit10_corrected`
- Denominator/support fields present: `n_total`; `n_valid`; `n_excluded`; `N_gallery_max`; `cross_alignment_contract`
- Scope restrictions: `scope=cross`; current materialized evidence is genetic-only; cross group rows are materialized under `dataset_or_direction=LINCS_to_scPerturb`; retrieval rows are directional in both `LINCS_to_scPerturb` and `scPerturb_to_LINCS`
- Allowed claim types: Figure 2 cross-ecosystem genetic reference-frame statements; group-level cross support; directional retrieval support
- Required joins, if any: use `task1_cross_alignment_proof.csv` for eligibility context and `task1_cross_attrition.csv` for metric-support loss notes
- Cautions or interpretation limits: cross chemical is excluded by contract and must not be read as a missing empirical signal; group-level comparison is symmetric and should not be misread as a two-direction leaderboard

### `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_summary.csv`

- Repo path: `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_summary.csv`
- Stage ownership: Task1 S2 cross metrics; run id `s2_task1_cross_metrics_0303`
- Evidence role: authoritative directional cross retrieval summary for Figure 2 support rows
- Main metrics present: `mean_mrr_raw`; `mean_expected_mrr_chance`; `mean_mrr_corrected`; `mean_hit1_raw`; `mean_expected_hit1_chance`; `mean_hit1_corrected`; `mean_hit5_raw`; `mean_expected_hit5_chance`; `mean_hit5_corrected`; `mean_hit10_raw`; `mean_expected_hit10_chance`; `mean_hit10_corrected`
- Denominator/support fields present: `n_total`; `n_valid`; `n_excluded_missing_metric_or_mpos0`; `N_gallery_max`
- Scope restrictions: genetic-only cross retrieval in two directions; direction must remain explicit
- Allowed claim types: supporting Task1 cross retrieval interpretation; cross-direction support statements for Figure 2
- Required joins, if any: pair with `task1_cross_chance_identity_check.csv` for corrected-metric validation
- Cautions or interpretation limits: retrieval is supporting evidence only; do not collapse `LINCS_to_scPerturb` and `scPerturb_to_LINCS`

### `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_attrition.csv`

- Repo path: `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_attrition.csv`
- Stage ownership: Task1 S2 cross metrics; run id `s2_task1_cross_metrics_0303`
- Evidence role: support-loss and exclusion table for cross Task1
- Main metrics present: none; attrition keyed by `reason`
- Denominator/support fields present: `reason`; `n_dropped`; `n_total_before`; `notes`
- Scope restrictions: cross-only support table
- Allowed claim types: denominator disclosure; cross exclusion support notes; energy-distance support caveats
- Required joins, if any: join on `scope`; `dataset_or_direction`; `perturbation_type`; `representation`
- Cautions or interpretation limits: current rows include `cross_matched_keys_lt_5` for chemical exclusion and `edist_insufficient_cells` for cross energy distance; these are support notes, not performance summaries

### `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_alignment_proof.csv`

- Repo path: `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_alignment_proof.csv`
- Stage ownership: Task1 S2 cross metrics; run id `s2_task1_cross_metrics_0303`
- Evidence role: frozen alignment-contract and eligibility proof for cross Task1
- Main metrics present: none; this is a contract-evidence table
- Denominator/support fields present: `cross_alignment_contract`; `perturbation_type`; `n_matched_keys`; `eligible_bool`; `excluded_reason`
- Scope restrictions: cross-only; perturbation-type level eligibility
- Allowed claim types: cross genetic eligibility; cross chemical exclusion support; contract disclosure for alignment key use
- Required joins, if any: none for policy disclosure; combine with `docs/contracts/exclusions-and-policies.md` if claim text needs the frozen policy wording
- Cautions or interpretation limits: current materialized chemical row records `n_matched_keys=0` with `matched_keys_lt5`; do not infer any alternate count or alignment rule

### `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_chance_identity_check.csv`

- Repo path: `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_chance_identity_check.csv`
- Stage ownership: Task1 S2 cross metrics; run id `s2_task1_cross_metrics_0303`
- Evidence role: validation-only support for corrected cross retrieval metrics
- Main metrics present: `delta_mrr`; `abs_delta_mrr`; `delta_hit1`; `abs_delta_hit1`; `delta_hit5`; `abs_delta_hit5`; `delta_hit10`; `abs_delta_hit10`
- Denominator/support fields present: keyed by `scope`; `dataset_or_direction`; `perturbation_type`; `representation`
- Scope restrictions: cross retrieval validation only
- Allowed claim types: validation support for Figure 2 cross retrieval interpretation
- Required joins, if any: join to `task1_cross_retrieval_summary.csv` on the shared keys when carrying corrected-metric validation
- Cautions or interpretation limits: validation-only; not a primary scientific evidence table

### `docs/contracts/exclusions-and-policies.md`

- Repo path: `docs/contracts/exclusions-and-policies.md`
- Stage ownership: frozen project contract policy
- Evidence role: frozen policy language for cross chemical exclusion; alignment-contract wording; scope-policy disclosure
- Main metrics present: none
- Denominator/support fields present: policy statements only
- Scope restrictions: contract language only; not an empirical result table
- Allowed claim types: exclusion-language support; scope-policy disclosure where the manuscript must distinguish policy exclusion from missing evidence
- Required joins, if any: none; cite alongside `task1_cross_alignment_proof.csv` or `task1_cross_attrition.csv` when explaining cross chemical exclusion
- Cautions or interpretation limits: use for policy wording only; current policy prose contains a numeric rationale that is not supported by the materialized S2 support tables, so manuscript-facing planning files should use only the safe wording that cross chemical is excluded by frozen policy/support gating

## Corrected Task2 S4 Group Concordance

### `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv`

- Repo path: `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv`
- Stage ownership: corrected Task2 S4 group concordance multisource; run id `s4_multisource_impl_verify_20260310_c`
- Evidence role: authoritative target-level mechanism-concordance table; upstream source of Task2 group summaries
- Main metrics present: `cosine_centroid`; `pcc_centroid`; `edist_biascorr`; with `*_valid_bool` and `*_na_reason` fields
- Denominator/support fields present: `n_chem_instances_total`; `n_gen_instances_total`; `n_chem_instances_used`; `n_gen_instances_used`; `n_chem_sub`; `n_gen_sub`
- Scope restrictions: keyed by `dataset`; `cell_line`; `target_token`; `representation`; corrected Task2 core evidence remains stratified by `dataset` and `cell_line`; no raw multisource pooling
- Allowed claim types: Figure 3 primary group-level concordance claims; dataset and cell-line stratified mechanism-concordance statements; local scPerturb K562 FM group statements
- Required joins, if any: none for direct target-level evidence; use S6 summaries for representation-level manuscript extraction and resolve ambiguities back to this table
- Cautions or interpretation limits: `specificity_tier` is not carried here; `n_targets` is not materialized here as a cohort field; energy distance is not a cross-representation common scale

### `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_attrition.csv`

- Repo path: `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_attrition.csv`
- Stage ownership: corrected Task2 S4 group concordance multisource; run id `s4_multisource_impl_verify_20260310_c`
- Evidence role: authoritative support-loss table for Task2 group metrics
- Main metrics present: attrition rows keyed by `metric_name` and `reason`
- Denominator/support fields present: `n_chem_before`; `n_chem_after`; `n_gen_before`; `n_gen_after`; `n_chem_removed`; `n_gen_removed`; `notes`
- Scope restrictions: only records rows with attrition events
- Allowed claim types: denominator transparency; UCE valid-mask disclosure in scPerturb K562; support caveats for group-level energy distance
- Required joins, if any: join on `dataset`; `cell_line`; `target_token`; `representation`; `metric_name`
- Cautions or interpretation limits: attrition explanations are support context only; current materialized non-zero FM attrition is UCE-specific in scPerturb K562

## Corrected Task2 S5 Retrieval

### `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv`

- Repo path: `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv`
- Stage ownership: corrected Task2 S5 retrieval multisource; run id `s5_multisource_impl_verify_20260311_a`
- Evidence role: authoritative direction-specific retrieval summary table for Figure 3 support rows
- Main metrics present: `mean_mrr_raw`; `mean_expected_mrr_chance`; `mean_mrr_corrected`; `mean_hit1_raw`; `mean_expected_hit1_chance`; `mean_hit1_corrected`; `mean_hit5_raw`; `mean_expected_hit5_chance`; `mean_hit5_corrected`; `mean_hit10_raw`; `mean_expected_hit10_chance`; `mean_hit10_corrected`
- Denominator/support fields present: `gallery_definition_id`; `pos_definition_id`; `n_total`; `n_valid`; `n_excluded_missing_metric_or_mpos0`; `N_gallery_mean`; `N_gallery_max`; `m_pos_mean`; `m_pos_p50`; `m_pos_p90`
- Scope restrictions: keyed by `dataset`; `cell_line`; `direction`; `representation`; directions must remain separate; no raw dataset pooling
- Allowed claim types: Figure 3 retrieval support; direction-specific mechanism-retrieval statements; local scPerturb K562 FM retrieval support
- Required joins, if any: none for summary extraction; use `task2_retrieval_per_query.parquet` plus snapshot metadata only when query-level stratification is required
- Cautions or interpretation limits: retrieval is supporting evidence only; `gallery_definition_id` and `pos_definition_id` must stay attached to every retrieval claim; do not collapse `C2G` and `G2C`

### `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary_long.csv`

- Repo path: `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary_long.csv`
- Stage ownership: corrected Task2 S5 retrieval multisource; run id `s5_multisource_impl_verify_20260311_a`
- Evidence role: long-format retrieval table for metric-wise extraction and plotting
- Main metrics present: `mean_mrr_corrected`; `mean_hit1_corrected`; `mean_hit5_corrected`; `mean_hit10_corrected`
- Denominator/support fields present: `n_valid`; `N_gallery_mean`; `N_gallery_max`; `m_pos_mean`; `m_pos_p50`; `m_pos_p90`
- Scope restrictions: direction-specific only
- Allowed claim types: metric-wise Figure 3 retrieval support rows; compact extraction into claim map or plotting plans
- Required joins, if any: none for current extraction
- Cautions or interpretation limits: corrected-only long format omits raw and expected retrieval metrics, so resolve correction questions back to `task2_retrieval_summary.csv`

### `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv`

- Repo path: `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv`
- Stage ownership: corrected Task2 S5 retrieval multisource; run id `s5_multisource_impl_verify_20260311_a`
- Evidence role: retrieval support-loss table for query or gallery valid-mask drops
- Main metrics present: none; attrition rows keyed by `reason`
- Denominator/support fields present: `n_query_rows_before`; `n_query_rows_after`; `n_query_rows_removed`; `n_gallery_member_rows_before`; `n_gallery_member_rows_after`; `n_gallery_member_rows_removed`; `n_gallery_items_before`; `n_gallery_items_after`; `n_gallery_items_removed`; `notes`
- Scope restrictions: only non-zero attrition rows are materialized; current table contains scPerturb K562 UCE valid-mask drops
- Allowed claim types: UCE caution support; supplementary attrition disclosure for Figure 3 retrieval
- Required joins, if any: join on `dataset`; `cell_line`; `direction`; `representation`
- Cautions or interpretation limits: this is a support table only; absence of a row should not be over-read as a performance comparison

### `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_chance_identity_check.csv`

- Repo path: `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_chance_identity_check.csv`
- Stage ownership: corrected Task2 S5 retrieval multisource; run id `s5_multisource_impl_verify_20260311_a`
- Evidence role: validation-only support that corrected Task2 retrieval metrics equal raw minus expected chance within tolerance
- Main metrics present: `delta_mrr`; `abs_delta_mrr`; `delta_hit1`; `abs_delta_hit1`; `delta_hit5`; `abs_delta_hit5`; `delta_hit10`; `abs_delta_hit10`
- Denominator/support fields present: keyed by `dataset`; `cell_line`; `direction`; `representation`
- Scope restrictions: retrieval validation only
- Allowed claim types: validation support for corrected retrieval interpretation
- Required joins, if any: join to `task2_retrieval_summary.csv` on the shared keys when validation provenance must be carried
- Cautions or interpretation limits: do not report these deltas as benchmark findings; use only to validate retrieval-correction bookkeeping

### `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet`

- Repo path: `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet`
- Stage ownership: corrected Task2 S5 retrieval multisource; run id `s5_multisource_impl_verify_20260311_a`
- Evidence role: authoritative query-level retrieval evidence for any frozen stratification that is not carried in S5 summaries
- Main metrics present: per-query `mrr_raw`; `hit1_raw`; `hit5_raw`; `hit10_raw`; `expected_mrr_chance`; `expected_hit1_chance`; `expected_hit5_chance`; `expected_hit10_chance`; `mrr_corrected`; `hit1_corrected`; `hit5_corrected`; `hit10_corrected`
- Denominator/support fields present: `query_row_id`; `query_uid`; `gallery_definition_id`; `pos_definition_id`; `N_gallery`; `m_pos`; `rank_true`; `query_n_targets`; `query_target_tokens`; `query_target_token`; `query_time`; `query_dose_value`
- Scope restrictions: query-level table only; `specificity_tier` is not carried directly; direction remains explicit
- Allowed claim types: query-level `n_targets` extraction for retrieval; audit-level retrieval provenance; blocked `specificity_tier` join attempts can resolve row identity from this table
- Required joins, if any: for metadata enrichment join `query_row_id` to snapshot `delta_meta.row_id` within the matching dataset subtree; no frozen aggregation rule exists yet for promoting `specificity_tier` or `n_targets` slices to figure-ready cohort summaries
- Cautions or interpretation limits: do not infer an unapproved aggregation rule from row-level availability; `query_n_targets` is materialized but only query-level

## Corrected Task2 S6 Result Synthesis

### `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_concordance_long.csv`

- Repo path: `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_concordance_long.csv`
- Stage ownership: corrected Task2 S6 result synthesis multisource; run id `s6_multisource_impl_verify_20260311_a`
- Evidence role: long-format group table for manuscript extraction and plotting with explicit metric-validity columns
- Main metrics present: `cosine_centroid`; `pcc_centroid`; `edist_biascorr` through `metric_name` and `metric_value`
- Denominator/support fields present: `metric_valid_bool`; `metric_na_reason`; `n_chem_instances_total`; `n_gen_instances_total`; `n_chem_instances_used`; `n_gen_instances_used`; `n_chem_sub`; `n_gen_sub`
- Scope restrictions: target-level long format; no dataset pooling; no `specificity_tier`
- Allowed claim types: metric-wise target-level extraction; ambiguity resolution for group-level summaries
- Required joins, if any: none for current extraction; resolve manuscript-level representation summaries back to this table if S6 leaderboard interpretation is unclear
- Cautions or interpretation limits: target-level long format is not itself the manuscript ranking table; energy distance remains informational-only for cross-representation reading

### `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv`

- Repo path: `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv`
- Stage ownership: corrected Task2 S6 result synthesis multisource; run id `s6_multisource_impl_verify_20260311_a`
- Evidence role: primary manuscript-facing Task2 group summary by `dataset`; `cell_line`; and `representation`
- Main metrics present: `mean_cosine_centroid`; `mean_pcc_centroid`; `mean_edist_biascorr`; ranking columns for cosine and PCC
- Denominator/support fields present: `n_targets_total`; `n_targets_cosine_valid`; `n_targets_pcc_valid`; `n_targets_edist_valid`; `n_invalid_rows_unique`; `n_invalid_drug_rows_unique`; `n_attrition_target_rows`; `n_attrition_chem_rows_removed_membership`; `n_attrition_gen_rows_removed_membership`
- Scope restrictions: representation ranking is only meaningful within the same `dataset` and `cell_line`; corrected Task2 core interpretation remains dataset-stratified and cell-line-stratified
- Allowed claim types: Figure 3 primary group-level mechanism-concordance rows; local scPerturb K562 FM group panel
- Required joins, if any: use S4 group tables when target-level support must be resolved; use `representation_availability_registry.csv` for scope disclosure
- Cautions or interpretation limits: `mean_edist_biascorr` is informational-only and not a cross-representation common scale; current scPerturb K562 UCE rows carry invalid-row attrition that must be disclosed if used

### `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`

- Repo path: `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`
- Stage ownership: corrected Task2 S6 result synthesis multisource; run id `s6_multisource_impl_verify_20260311_a`
- Evidence role: primary manuscript-facing Task2 retrieval summary by `dataset`; `cell_line`; `direction`; and `representation`
- Main metrics present: `mean_mrr_corrected`; `mean_hit1_corrected`; `mean_hit5_corrected`; `mean_hit10_corrected`; `rank_by_mean_mrr_corrected`
- Denominator/support fields present: `gallery_definition_id`; `pos_definition_id`; `n_total`; `n_valid`; `n_excluded_missing_metric_or_mpos0`; `N_gallery_mean`; `N_gallery_max`; `m_pos_mean`; `m_pos_p50`; `m_pos_p90`
- Scope restrictions: direction-specific; no raw pooling across datasets; ranking only within the same `dataset`; `cell_line`; and `direction`
- Allowed claim types: Figure 3 retrieval support rows; local scPerturb K562 FM retrieval support
- Required joins, if any: resolve any summary ambiguity back to `task2_retrieval_summary.csv`
- Cautions or interpretation limits: keep `C2G` and `G2C` separate; carry `gallery_definition_id` and `pos_definition_id` with the claim

### `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_benchmark_summary_long.csv`

- Repo path: `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_benchmark_summary_long.csv`
- Stage ownership: corrected Task2 S6 result synthesis multisource; run id `s6_multisource_impl_verify_20260311_a`
- Evidence role: synthesis-layer union table that preserves upstream provenance, caution codes, and manuscript-friendly denominator carry-through
- Main metrics present: union of group and retrieval rows through `analysis_family`; `metric_name`; `metric_value`
- Denominator/support fields present: group-side `n_targets_total`; `n_targets_metric_valid`; retrieval-side `n_total`; `n_valid`; `n_excluded_missing_metric_or_mpos0`; `N_gallery_mean`; `N_gallery_max`; `m_pos_mean`; `m_pos_p50`; `m_pos_p90`; plus `source_table`; `source_stage`; `source_run_id`; `caution_codes`
- Scope restrictions: synthesis-only layer; use to inventory panel readiness and caution carry-through rather than to override S4 or S5
- Allowed claim types: manuscript panel-readiness indexing; caution-code extraction; compact figure inventory
- Required joins, if any: none for current extraction; resolve ambiguities against the upstream S4 or S5 authoritative tables
- Cautions or interpretation limits: this is a summary union, not a license to re-aggregate across datasets, cell lines, directions, or representations

## Frozen Task2 Snapshot v2 Metadata

### `data/task2_snapshot_v2/task2_pairs_coverage.csv`

- Repo path: `data/task2_snapshot_v2/task2_pairs_coverage.csv`
- Stage ownership: frozen Task2 snapshot v2 metadata
- Evidence role: authoritative eligibility and coverage fact table for Task2 mechanism cohorts
- Main metrics present: none; coverage facts only
- Denominator/support fields present: `n_chem_instances`; `n_gen_instances`; `is_eligible_bool`; `eligibility_reason`; `source_tag`
- Scope restrictions: keyed by `dataset`; `cell_line`; `target_token`; not a metric table
- Allowed claim types: coverage transparency; eligible-versus-ineligible denominator context; supplementary scope notes
- Required joins, if any: join to S4 or S6 group tables on `dataset`; `cell_line`; `target_token`
- Cautions or interpretation limits: do not use this file as a performance table; eligibility is distinct from representation scope

### `data/task2_snapshot_v2/representation_availability_registry.csv`

- Repo path: `data/task2_snapshot_v2/representation_availability_registry.csv`
- Stage ownership: frozen Task2 snapshot v2 metadata
- Evidence role: authoritative representation-scope registry that separates supported scope from attrition
- Main metrics present: none; scope registry only
- Denominator/support fields present: `availability_status`; `availability_reason`
- Scope restrictions: keyed by `dataset`; `cell_line`; `representation`
- Allowed claim types: FM scope disclosure; LINCS not-applicable-scope disclosure; support notes that scope absence is not attrition
- Required joins, if any: join to S4, S5, or S6 tables on `dataset`; `cell_line`; `representation`
- Cautions or interpretation limits: use only for scope policy; not a performance or denominator-loss table

### `data/task2_snapshot_v2/lincs/derived/delta_meta.csv`

- Repo path: `data/task2_snapshot_v2/lincs/derived/delta_meta.csv`
- Stage ownership: frozen Task2 snapshot v2 metadata; LINCS subtree
- Evidence role: authoritative LINCS row metadata for Task2 query-level joins and complexity context
- Main metrics present: none; metadata only
- Denominator/support fields present: `row_id`; `source_row_index`; `treated_cell_id`; `perturbation_class`; `cell_line`; `target_raw`; `time`; `dose_value`; `specificity_tier`; `dataset_side`; `target_tokens`
- Scope restrictions: LINCS only; `specificity_tier` is currently materialized as `NA`; LINCS supports only `Gene` and `Pathway`
- Allowed claim types: query-level `n_targets` derivation if needed; LINCS metadata joins for retrieval audit context; evidence that LINCS does not support a meaningful `specificity_tier` slice in the current frozen state
- Required joins, if any: join S5 per-query `query_row_id` to `row_id` for LINCS rows
- Cautions or interpretation limits: because `specificity_tier` is `NA` here, any benchmark-wide specificity claim would exceed the frozen local evidence

### `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`

- Repo path: `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`
- Stage ownership: frozen Task2 snapshot v2 metadata; scPerturb K562 subtree
- Evidence role: authoritative scPerturb K562 row metadata for query-level stratification joins
- Main metrics present: none; metadata only
- Denominator/support fields present: `row_id`; `treated_cell_id`; `perturbation_class`; `cell_line`; `target_raw`; `time`; `dose_value`; `specificity_tier`; `dataset_side`; `target_tokens`
- Scope restrictions: scPerturb K562 only; approved FM scope only exists in this subtree
- Allowed claim types: local FM scope context; query-level `n_targets` or `specificity_tier` joins for blocked-or-supplementary analyses; metadata support for K562-only disclosure
- Required joins, if any: join S5 per-query `query_row_id` to `row_id` for scPerturb rows
- Cautions or interpretation limits: `specificity_tier` is row-level only and not carried into S4 through S6 cohort summaries; no frozen aggregation rule currently promotes this metadata to a main-text panel
