# Manuscript Blocker Memo

This memo records unresolved evidence-extraction issues from the frozen local evidence path. It does not change the manuscript story or introduce new aggregation rules.

Sequencing note: blocked versus gated follow-up analyses should now be read together with `docs/manuscript_analysis_expansion_plan.md`, which defines the execution order and stop conditions for the next manuscript-support phase.

## Blocking Issues

### `C3_F3e_specificity_tier`

- Status: blocked
- Why blocked: `specificity_tier` is not carried in S4 through S6 summary tables.
- Available local evidence: row-level `specificity_tier` exists in `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`; LINCS `delta_meta.csv` carries `specificity_tier=NA`.
- Join status: S5 per-query retrieval can be joined to snapshot metadata by `query_row_id -> row_id`.
- Remaining problem: no frozen cohort aggregation rule exists for turning row-level tier labels into a manuscript-ready retrieval or group-level slice; the field is also direction-asymmetric because the intended biological specificity signal lives on chemical-query rows rather than both retrieval directions.
- Extraction decision: keep blocked; do not infer bins; do not promote to main text or supplementary summary claims without an explicit frozen rule.

### `C3_F3f_n_targets_complexity`

- Status: blocked
- Why blocked: `query_n_targets` is already materialized in `task2_retrieval_per_query.parquet`, but only at query level.
- Available local evidence: exact query-level `query_n_targets`; fallback derivation from `target_tokens` also exists in snapshot metadata.
- Remaining problem: no frozen aggregation or binning rule exists for a manuscript-ready complexity panel; the intended chemical-complexity axis is not symmetric across directions, so promoting a generic bidirectional complexity statement would require new interpretation rules.
- Extraction decision: keep blocked; do not create single-target versus multi-target bins or any other pooled complexity summary.

## Ambiguities That Need Caution Language

### `C2_F2g_cross_chemical_policy_exclusion`

- Status: extractable with caution
- Local conflict: the frozen policy prose contains a numeric rationale that is not supported by the current materialized S2 support tables, which record an ineligible chemical cross gate through `n_matched_keys=0` and `cross_matched_keys_lt_5`.
- Safe manuscript handling: state only that cross chemical is excluded by frozen policy/support gating; do not quote a numeric matched-context rationale in manuscript-facing prose.
- Main-text implication: the exclusion statement itself can remain in main text; the numeric rationale cannot.

## Supplementary-Only or Local-Only Rows

### `C2_F2c_internal_group_scperturb_fm`

- Reason: these FM rows remain supplementary-only, and the underlying Task1 internal scPerturb evidence is not FM-only overall: local S1 outputs contain scPerturb `Gene` and `Pathway` rows and that common-scope slice should be indexed separately from the FM-local slice.
- Decision: keep supplementary-only; do not use as the main internal calibration anchor.

### `C2_F2d_internal_retrieval_scperturb_fm`

- Reason: same FM-local status as the group-level scPerturb FM row; retrieval is also support-only by contract. Common-scope scPerturb `Gene` and `Pathway` retrieval rows exist locally and should be handled separately from these FM rows.
- Decision: keep supplementary-only.

### `C3_F3d_coverage_context`

- Reason: `task2_pairs_coverage.csv`, `representation_availability_registry.csv`, and attrition tables support denominator transparency and scope interpretation, but they are not performance tables.
- Decision: keep supplementary-only or figure-note only.

### `C4_F3f_fm_retrieval_k562_local`

- Reason: FM retrieval rows are real and corrected, but Figure 3 keeps group-level concordance as the primary evidence layer and the FM scope is local to scPerturb K562.
- Decision: keep supplementary-only unless a specific local inset panel is explicitly needed.

## Main-Text-Safe With Required Cautions

### `C3_F3a_group_common_scope`

- Safe use: common-scope `Gene` and `Pathway` group-level Task2 statements stratified by `dataset` and `cell_line`
- Required caution: no raw pooling across LINCS and scPerturb; `mean_edist_biascorr` stays informational-only

### `C3_F3b_retrieval_c2g_common_scope` and `C3_F3c_retrieval_g2c_common_scope`

- Safe use: direction-specific retrieval support in common scope
- Required caution: keep `C2G` and `G2C` separate and carry `gallery_definition_id` plus `pos_definition_id`

### `C4_F3e_fm_group_k562_local` and `C4_F3g_fm_scope_policy`

- Safe use: local FM group comparison within scPerturb K562 only
- Required caution: carry K562-only scope language and note UCE attrition where the local FM comparison includes UCE
