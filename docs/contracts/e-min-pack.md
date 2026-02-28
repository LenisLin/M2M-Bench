# E-min Pack (Covariate + Question Layer) — Frozen

Purpose: report auditable facts about ecology/support/bias/uplift using ONLY internal metadata.
No external knowledge bases; no multivariable regression.

## Shared conventions
- All aggregations must be stratified by `metric_family` (no mixing).
- For any "uplift" analysis:
  - baseline = Gene space
  - success_def = BOTH
    - success_def_v1: (amrr_corrected > 0) AND finite
    - success_def_v2: (hit@K_corrected > 0) AND finite  (K is the headline K for retrieval)
  - group types: target, cell_line (NO pair-level by default)

## E1 Ecology Coverage
Output: `E1_ecology_coverage.csv`
Groupby:
- dataset_side in {scPerturb, LINCS, cross}
- task in {S4a,S4b,S5}
- perturbation_type in {Chemical,Genetic,NA}
- space_or_model in {Gene, Pathway, FM:<model>}
Columns (min):
- n_queries_total
- n_targets
- n_cell_lines
- median_q_per_target
- p95_q_per_target
- p99_q_per_target
- top10pct_target_query_share

## E2 Support & Attrition Funnel
Output: `E2_support_and_attrition.csv`
Must include:
- Candidate -> valid_query -> common_support -> headline_gate (e.g., N_min) survival counts
- Explicit policy flags:
  - cross_chemical_excluded=true, excluded_reason="only 2 contexts"
  - cross_alignment_contract="global_idx_lincs + sc_delta_row_idx"

## E3 Preference Ranking
Outputs:
- `E3_preference_rank_summary.csv`
- `E3_preference_rank_details.jsonl`
Definition:
- per-group AMRR macro (group = target or cell_line), with query threshold q_min in {50,100} (or project-configured)
- report top/bottom quantiles (e.g., top10%, bottom10%) within eligible groups
No external enrichment; purely ranking.

## E4 Preference Uplift Audit (formerly "enrichment", re-scoped)
Goal: identify whether success-rate improvement of tested representation is concentrated on some targets/cell_lines.

Outputs:
- `E4_uplift_audit.csv`
- `E4_uplift_audit_details.jsonl`

Inputs: per-query retrieval table with:
- dataset_side, task, perturbation_type, metric_family
- representation identifier (Gene/Pathway/FM:model)
- per-query success outcome for each success_def (see above)
- group keys (target/cell_line)

Procedure:
For each (dataset_side, task, perturbation_type, metric_family, tested_rep) and each group_type in {target, cell_line}:
- baseline_rep = Gene
- compute per-group success_rate_baseline, success_rate_tested
- uplift = tested - baseline
- attach CI (Wilson or bootstrap; choose one and freeze in code; no regression)
- optional hypothesis test (Fisher exact or permutation) + BH-FDR within the family:
  family = (dataset_side, task, perturbation_type, metric_family, tested_rep, group_type, success_def)

E4 Summary columns (min):
- success_def
- group_type
- tested_rep
- n_groups_eligible
- n_queries_used_total
- top_uplift_groups (stringified top-N list)
- bottom_uplift_groups
- fdr_threshold_used (if tests are enabled)

E4 Detail fields (per group):
- group_id, n_queries, success_rate_tested, success_rate_baseline, uplift, ci_low, ci_high, p_value, q_value
