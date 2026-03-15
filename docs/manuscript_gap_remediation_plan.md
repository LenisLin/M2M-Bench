# Manuscript Gap Remediation Plan

## 1. Purpose

This document plans how to close manuscript-supporting evidence gaps without changing frozen Task1 or Task2 contracts. The goal is to separate true blockers from under-extraction, identify the smallest valid downstream remediations, and stop before any action that would redefine benchmark scope, task semantics, or evidence hierarchy.

This is a remediation-planning pass only. It does not promote blocked claims, rerun major stages, or create new manuscript claims.

## 2. Current Gap Inventory

- `G1` Task1 scPerturb internal common-scope under-extraction: the current evidence index treated the scPerturb internal slice as FM-only, but local S1 outputs also contain scPerturb `Gene` and `Pathway` rows.
- `G2` `specificity_tier` as a core Figure 3 axis: the frozen blueprint names it as a primary stratification axis, but authoritative S4-S6 summaries do not carry it.
- `G3` local `specificity_tier` sensitivity: scPerturb K562 C2G query-level rows appear sufficient for a local downstream sensitivity table.
- `G4` `n_targets` as a core Figure 3 axis: per-query values exist, but there is no frozen cohort aggregation rule for turning them into a core summary panel.
- `G5` local `n_targets` retrieval sensitivity: C2G query-level support exists in both datasets, while G2C is structurally non-informative because all queries have `n_targets=1`.
- `G6` cross-chemical exclusion wording conflict: frozen policy prose cites "2 matched contexts" while materialized S2 support tables record `n_matched_keys=0` and generic `matched_keys_lt5` exclusion.
- `G7` FM local-scope panel ambiguity: the local scPerturb K562 FM evidence exists, but group-primary versus retrieval-supporting roles, K562-only scope language, and UCE attrition handling are not yet consolidated into one manuscript-facing support package.
- `G8` denominator/caution packaging for retained main-text slices: support fields and caution codes exist locally, but they are not yet assembled into a panel-facing ledger that can travel cleanly into manuscript prose and figure notes.

## 3. Gap Classification Table

| gap_id | gap_name | current_status | category | affected_claims_or_panels | local_evidence_available | proposed_action | new_script_needed | rerun_needed | risk_level | expected_gain | recommendation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `G1` | Task1 scPerturb internal common-scope under-extraction | Misindexed as FM-only asymmetry | `A` | Figure 2 internal calibration; supplementary Task1 completeness | `task1_leaderboard_long.csv`; `task1_retrieval_summary.csv`; `task1_chance_identity_check.csv`; S1 `run_manifest.json` | Add scPerturb `Gene` and `Pathway` Task1 internal rows to the claim map and source register; keep FM rows separate | No | No | Low | Medium | Pursue immediately |
| `G2` | `specificity_tier` as a core benchmark-wide / group-level Figure 3 axis | Blocked | `D` | Figure 3 main-text specificity axis; Claim 3 precision | `task2_retrieval_per_query.parquet`; `scperturb_k562/derived/delta_meta.csv`; `lincs/derived/delta_meta.csv`; S4-S6 summary tables | Keep blocked for core panels; do not force benchmark-wide or group-level specificity aggregation | No | No | High | Low | Keep as limitation; route recoverable demand to `G3` only |
| `G3` | Local scPerturb K562 C2G `specificity_tier` sensitivity | Blocked but data-rich | `C` | Supplementary Figure 3 sensitivity; discussion support | `task2_retrieval_per_query.parquet`; `scperturb_k562/derived/delta_meta.csv`; S5/S6 retrieval denominators | Build a downstream local-only sensitivity table by exact tier labels for `dataset=scPerturb`, `cell_line=K562`, `direction=C2G` | Yes | No | Medium | Medium-High | Pursue after panel-support ledger exists |
| `G4` | `n_targets` as a core bidirectional/group-level Figure 3 axis | Blocked | `B` | Figure 3 secondary complexity axis in core panels | `task2_retrieval_per_query.parquet`; snapshot `delta_meta.csv`; S4-S6 summaries | Do not treat as frozen core panel yet; first define downstream sensitivity extraction and keep any grouping rule manuscript-local | No | No | Medium | Medium | Do not promote to core panel until downstream sensitivity is reviewed |
| `G5` | Direction-aware `n_targets` retrieval sensitivity | Blocked but query-level values exist | `C` | Figure 3 complexity sensitivity; supplementary or figure-note support | `task2_retrieval_per_query.parquet` with `query_n_targets`; S5/S6 retrieval tables | Build a downstream C2G sensitivity table stratified by exact `query_n_targets`, dataset, representation, and support counts; mark G2C as non-informative | Yes | No | Medium | Medium-High | Pursue in the same pass as `G3` |
| `G6` | Cross-chemical exclusion wording conflict | Caution-heavy credibility issue | `C` | Figure 2 exclusion note; methods and limitations wording | `docs/contracts/exclusions-and-policies.md`; `task1_cross_alignment_proof.csv`; `task1_cross_attrition.csv`; `cross-pairs-alignment-summary.json`; `docs/contracts/task1_spec.md` | Reconcile local wording and remove unsupported numeric language unless a frozen local source supports it | No | No | Low | High | Highest priority before prose |
| `G7` | FM local-scope panel packaging and caution carry-through | Evidence exists, packaging incomplete | `A` | `F3e`; `F3f`; FM scope note | `task2_group_leaderboard.csv`; `task2_retrieval_leaderboard.csv`; `representation_availability_registry.csv`; S4/S5 attrition; S6 caution codes | Create a manuscript-facing local FM panel note or ledger fixing group-primary role, retrieval-supporting role, K562-only scope, and UCE caution handling | No | No | Low | Medium | Pursue after `G6` |
| `G8` | Denominator/caution packaging for retained main-text panels | Under-extracted | `A` | `F3a`; `F3b`; `F3c`; `F3e`; figure notes | `task2_benchmark_summary_long.csv`; `task2_group_leaderboard.csv`; `task2_retrieval_leaderboard.csv`; `task2_group_attrition.csv`; `task2_retrieval_attrition.csv`; `task2_pairs_coverage.csv`; `representation_availability_registry.csv` | Create a panel-support ledger carrying denominator fields, caution codes, and scope-vs-attrition distinctions for retained slices | Yes | No | Low | High | Pursue early; it de-risks all later prose |

## 4. Remediation Actions

### `G1` Task1 scPerturb internal common-scope under-extraction

Why category `A`:

- This is not a missing-data problem. The S1 run already materialized scPerturb `Gene` and `Pathway` rows in both `task1_leaderboard_long.csv` and `task1_retrieval_summary.csv`.
- The gap arose in the evidence-index layer, not in upstream computation.

Local files that prove recoverability:

- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv`
- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_summary.csv`
- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_chance_identity_check.csv`
- `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/run_manifest.json`

Smallest valid remediation path:

- Update `docs/manuscript_claim_to_evidence_map.csv` to add scPerturb internal common-scope `Gene` and `Pathway` rows.
- Update `docs/manuscript_source_register.md` notes so the scPerturb internal common-scope slice is explicitly represented.
- Keep FM-specific Task1 internal rows supplementary-only rather than using them as a stand-in for all scPerturb internal evidence.

Target output file(s):

- `docs/manuscript_claim_to_evidence_map.csv`
- `docs/manuscript_source_register.md`
- `docs/manuscript_blocker_memo.md` only if wording still implies FM-only scPerturb internal evidence

Expected manuscript impact:

- Strengthens Figure 2 completeness and removes a misleading asymmetry in the current evidence index.
- Mostly affects supplement and backup figure logic; it is not required to change the frozen Figure 2 main story.

### `G2` `specificity_tier` as a core benchmark-wide / group-level Figure 3 axis

Why category `D`:

- The field is absent from authoritative S4-S6 summaries.
- LINCS snapshot metadata carries `specificity_tier=NA`, so benchmark-wide specificity reporting cannot be supported from both datasets.
- A group-level specificity analysis would require a new membership-to-cohort aggregation rule that is not frozen locally.

Local files that support non-recoverability:

- `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv`
- `data/task2_snapshot_v2/lincs/derived/delta_meta.csv`
- `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`

Smallest valid remediation path:

- None for a core benchmark-wide or group-level panel.
- Keep this as a limitation and avoid promising a benchmark-wide specificity axis in results prose unless the manuscript explicitly narrows the scope to the local sensitivity path in `G3`.

Target output file(s):

- No new result file recommended
- Limitation wording only

Expected manuscript impact:

- Prevents an invalid scope expansion.
- Protects Claim 3 from over-promising a benchmark-wide specificity result that the current local evidence cannot support.

### `G3` Local scPerturb K562 C2G `specificity_tier` sensitivity

Why category `C`:

- The necessary data already exist locally at query level, but the relevant summary must be created downstream.
- The tier labels are already frozen categories in metadata, so the local-only sensitivity can avoid inventing new bins.

Local files that suggest recoverability:

- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet`
- `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`

Minimal downstream analysis:

- Join `query_row_id` in Task2 per-query retrieval to scPerturb K562 `delta_meta.row_id`.
- Restrict to `dataset=scPerturb`, `cell_line=K562`, `direction=C2G`.
- Aggregate corrected retrieval metrics by exact `specificity_tier`, `representation`, and denominator fields.
- Report this only as a local downstream sensitivity table, not as benchmark-wide evidence.

Target output file(s):

- `docs/manuscript_support_tables/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv`
- Optional short readme or note: `docs/manuscript_support_tables/README.md`

New code / rerun:

- A small downstream extractor script is recommended, for example `scripts/manuscript_extract_support_tables.py`.
- No upstream rerun is needed.

Expected manuscript impact:

- Can strengthen supplementary biological interpretation and discussion discipline.
- May support a local scPerturb specificity sensitivity note, but should not replace the current benchmark-wide limitation.

### `G4` `n_targets` as a core bidirectional/group-level Figure 3 axis

Why category `B`:

- The data exist locally in per-query form, but there is no frozen summary rule that turns them into a core panel.
- `G2C` is structurally non-informative because all queries have `query_n_targets=1`.
- Group-level complexity interpretation would require a new aggregation rule beyond the current frozen summaries.

Local files that show both the opportunity and the constraint:

- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`
- `data/task2_snapshot_v2/lincs/derived/delta_meta.csv`
- `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`

Smallest valid remediation path:

- Do not force `n_targets` into a core bidirectional or group-level panel.
- Route complexity support through the downstream sensitivity path in `G5`.
- Keep any display grouping explicitly manuscript-local if later needed.

Target output file(s):

- No core-result file recommended at this stage
- Decision note can remain in `docs/manuscript_gap_remediation_plan.md` and later in the blocker memo if needed

Expected manuscript impact:

- Prevents accidental contract drift while preserving a valid path for local sensitivity analysis.

### `G5` Direction-aware `n_targets` retrieval sensitivity

Why category `C`:

- `query_n_targets` is already materialized in the authoritative Task2 per-query retrieval output.
- A C2G-only, dataset-stratified downstream table can be built without touching upstream stages.

Local files that suggest recoverability:

- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`

Minimal downstream analysis:

- Restrict to `direction=C2G`.
- Aggregate corrected retrieval metrics by exact `query_n_targets`, `dataset`, `cell_line`, and `representation`.
- Keep LINCS and scPerturb separate.
- Treat `G2C` as a documented non-informative branch rather than forcing a symmetric panel.
- If a later plotted display needs coarse bins, define them in a manuscript-only sensitivity note after inspecting the exact-value table, not before.

Target output file(s):

- `docs/manuscript_support_tables/task2_c2g_query_n_targets_sensitivity.csv`
- Optional counts-only companion: `docs/manuscript_support_tables/task2_c2g_query_n_targets_counts.csv`

New code / rerun:

- Use the same small downstream extractor script as `G3`.
- No upstream rerun is needed.

Expected manuscript impact:

- Strongest recoverable route for the target-complexity story.
- Best suited to supplement, figure note, or a bounded sensitivity panel rather than an immediate core benchmark-wide claim.

### `G6` Cross-chemical exclusion wording conflict

Why category `C`:

- This is a local reconciliation problem, not an upstream analysis problem.
- The credibility risk is real because the current frozen policy prose and the current materialized support tables disagree on the numeric rationale.

Local files that matter:

- `docs/contracts/exclusions-and-policies.md`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_alignment_proof.csv`
- `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_attrition.csv`
- `data/task1_snapshot_v1/cross_contract/cross-pairs-alignment-summary.json`
- `docs/contracts/task1_spec.md`

Smallest valid remediation path:

- First, check whether any frozen local artifact supports the numeric "2 matched contexts" wording.
- If no such artifact exists, patch the policy language to numeric-free wording aligned with the materialized contract outcome: chemical cross is excluded because it fails the current eligibility/support gate.
- Update the source register and claim map notes to remove any temptation to quote the unsupported number.

Target output file(s):

- `docs/contracts/exclusions-and-policies.md`
- `docs/manuscript_source_register.md`
- `docs/manuscript_claim_to_evidence_map.csv`
- Optional reconciliation note: `docs/manuscript_cross_exclusion_reconciliation_note.md`

Expected manuscript impact:

- High credibility gain with minimal analytical risk.
- Should be resolved before any results are written into prose.

### `G7` FM local-scope panel packaging and caution carry-through

Why category `A`:

- All required local evidence already exists.
- The gap is that group evidence, retrieval evidence, scope policy, and UCE attrition are not yet assembled into a single manuscript-facing package.

Local files that support immediate extraction:

- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`
- `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_attrition.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv`
- `data/task2_snapshot_v2/representation_availability_registry.csv`

Smallest valid remediation path:

- Produce a local FM scope ledger that records: `dataset=scPerturb`; `cell_line=K562`; group-level is primary; retrieval is supporting; UCE carries attrition notes; LINCS FM absence is scope policy rather than attrition.
- This can be done as documentation or a lightweight CSV without new benchmark computation.

Target output file(s):

- `docs/manuscript_fm_scope_panel_note.md`
- Optional CSV: `docs/manuscript_support_tables/task2_fm_k562_scope_ledger.csv`

Expected manuscript impact:

- Makes FM reporting safer and more concise.
- Preserves the current local-only FM interpretation while reducing ambiguity during results writing.

### `G8` Denominator/caution packaging for retained main-text panels

Why category `A`:

- The denominator fields and caution codes already exist in current S4-S6 and snapshot metadata.
- The missing piece is a panel-facing ledger that makes those fields easy to carry into prose, captions, and supplementary methods.

Local files that support immediate packaging:

- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_benchmark_summary_long.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv`
- `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv`
- `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_attrition.csv`
- `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv`
- `data/task2_snapshot_v2/task2_pairs_coverage.csv`
- `data/task2_snapshot_v2/representation_availability_registry.csv`

Smallest valid remediation path:

- Create one manuscript-only ledger keyed by claim row or panel slice.
- Carry forward denominator fields, caution codes, availability-scope notes, and attrition notes for every retained main-text slice.
- Use that ledger to keep `GROUP_EDIST_NOT_CROSS_REP`, `GROUP_EDIST_PARTIAL_TARGET_VALIDITY`, and `UCE_ATTRITION_PRESENT` from being dropped in later prose.

Target output file(s):

- `docs/manuscript_panel_support_ledger.csv`
- Optional note: `docs/manuscript_panel_support_ledger.md`

New code / rerun:

- A small downstream extractor script is useful for reproducibility.
- No upstream rerun is needed.

Expected manuscript impact:

- High. This lowers the risk of denominator drift and keeps caution language synchronized with the actual retained slices.

## 5. Priority Ranking

1. Resolve `G6` cross-chemical exclusion wording conflict.
Reason: it is a credibility risk in current manuscript-facing text and does not require any new analysis.

2. Build `G8` panel-support ledger.
Reason: it strengthens every retained Figure 3 slice and reduces the chance of denominator/caution drift when results prose is drafted.

3. Fix `G1` Task1 scPerturb internal common-scope under-extraction.
Reason: low-risk, immediate correction of the evidence index; improves Figure 2 completeness without changing the story.

4. Build `G5` direction-aware `n_targets` retrieval sensitivity.
Reason: it is the strongest recoverable route for the target-complexity interpretation and directly supports Figure 3 biological reading.

5. Build `G3` local scPerturb K562 C2G `specificity_tier` sensitivity.
Reason: it can strengthen a local sensitivity story, but it should remain local and should not outrank the broader denominator/complexity work.

6. Package `G7` FM local-scope panel support.
Reason: worthwhile, but lower priority than the main Figure 3 common-scope and wording-stability remediations.

7. Do not pursue `G2` as a core benchmark-wide panel and do not promote `G4` into a frozen core axis before the downstream sensitivity tables are reviewed.
Reason: both would risk contract drift or overclaiming.

## 6. Stop Rules

- Stop if the remediation would require rerunning major upstream benchmark stages rather than using existing local outputs.
- Stop if the remediation would redefine Task1 or Task2 scope, task semantics, or benchmark contracts.
- Stop if the remediation requires raw pooling across LINCS and scPerturb for core Task2 claims.
- Stop if the remediation would generalize FM evidence beyond the supported scPerturb K562 scope.
- Stop if the remediation depends on inventing a benchmark-wide aggregation rule that is not frozen locally.
- Stop if the remediation only improves cosmetic symmetry and does not materially strengthen Figure 2, Figure 3, denominator transparency, or manuscript credibility.
- Stop if the remediation still cannot produce exact denominator/support fields for the proposed manuscript slice.
- Stop if the only way to resolve a wording conflict is to quote a number that lacks frozen local support.

## 7. Recommended Next Execution Order

1. Reconcile the cross-chemical exclusion wording and remove unsupported numeric language unless a frozen local source supports it.
2. Correct the evidence index for Task1 scPerturb internal common-scope `Gene` and `Pathway` rows.
3. Create one shared downstream support extractor script that reads only existing S1/S4/S5/S6 outputs and snapshot metadata.
4. Use that script first to emit `docs/manuscript_panel_support_ledger.csv`.
5. Use the same script to emit `task2_c2g_query_n_targets_sensitivity.csv`.
6. Use the same script to emit `task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv`.
7. Package the local FM scope ledger or note so K562-only scope and UCE cautions are carried forward cleanly.
8. Reassess whether the downstream sensitivity tables are strong enough for supplementary-only promotion or figure-note use; if not, keep them as limitations and stop there.
