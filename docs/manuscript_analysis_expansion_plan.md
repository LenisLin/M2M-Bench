# Manuscript Analysis Expansion Plan

## 1. Purpose

Create `docs/manuscript_analysis_expansion_plan.md` as the planning source of truth for the next manuscript-support phase. This phase plans document, code, and output changes only; it does not run new analyses, invent results, rewrite benchmark contracts, or broaden the benchmark beyond what is needed to support manuscript argumentation.

- Highest-priority analyses: `A3` direction robustness and `B6` dose/time overlap gating, then the bridge analyses `A1-A2`, then statistical packaging in `A4`.
- Highest-priority document changes: `docs/manuscript_blueprint.md`, `docs/manuscript_evidence_extraction_plan.md`, `docs/manuscript_source_register.md`, `docs/manuscript_claim_to_evidence_map.csv`, `docs/manuscript_blocker_memo.md`, `docs/manuscript_panel_support_ledger.csv`, plus a new manuscript-support table spec doc.
- Highest-priority code changes: new downstream `scripts/manuscript_analysis_support.py`, new downstream `scripts/manuscript_enrichment_support.py`, and only a conditional last-resort S1 export if existing S1 outputs are proven insufficient for `A1/A2/A4`.

## 2. Guiding Principles

- Preserve frozen benchmark contracts and keep `Task1` and `Task2` core definitions unchanged.
- Prioritize manuscript-critical analyses before reviewer-defense and optional exploratory work.
- Separate main-bridge analyses, explanation-layer analyses, representation-layer analyses, and case-layer illustration in both wording and execution order.
- Prefer downstream manuscript-support analyses over upstream stage changes; no upstream touch should happen by default.
- Treat `docs/manuscript_blueprint.md` as the canonical revised-proposal replacement; `docs/manuscript_proposal.md` is absent locally and is not a planned update in this phase.
- Keep new manuscript-support outputs in `docs/` root for consistency with the existing `docs/task2_*` support tables.
- Use matched and evaluable target pools as enrichment backgrounds; never use the whole genome.
- Keep S7 as context only, not as the primary evidence source for new manuscript claims.
- This phase exists to strengthen manuscript logic, not to launch a second benchmark-expansion project.

## 3. Analysis Inventory

`Deps`: `OA` overlap audit, `NAR` new downstream aggregation rule, `EST` existing summary table, `QJ` query-level table join, `MST` new manuscript-support table.

| ID | Analysis | Status | Need / support / reviewer question | Inputs | Outputs | Gates | Class | Risk | reviewer_attack_vector | Deps |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `A1` | Task1 internal vs cross common-scope comparison | Partial | Modality-level measurable comparability bridge for Result 2 / Figure 2; answers "is supported cross comparability weaker than internal comparability on the same common-scope genetic rows without invoking ceiling language?" | S1 leaderboard + per-query; S2 cross leaderboard + `task1_group_cross.parquet` + per-query | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv` | Gene and Pathway only; genetic-only matched rows; use downstream support tables first, and consider S1 internal-group export only if current S1 outputs prove insufficient | Required, main bridge | High | `virtual_cell_scope_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `A2` | Task1 vs Task2 shared-slice bridge comparison | Missing | Shared-slice external-validity bridge for Results 2-3; answers "when the same interpretable and common-scope slice is available, how much modality-level comparability carries into mechanism concordance?" | Existing S1, S4, and S6 outputs; S1 internal-group export only if needed; S4 target-level Task2; S6 group and retrieval summaries | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv` | Match only on shared `dataset`, `cell_line`, `target_token`, and `representation`; no generic "Task1 vs Task2 worse" framing; no ceiling language | Required, main bridge | High | `virtual_cell_scope_challenge` | OA:N NAR:Y EST:Y QJ:N MST:Y |
| `A3` | C2G vs G2C direction-specific robustness audit | Partial | Reviewer-defense for Figure 3 retrieval; answers "is direction asymmetry only a gallery or null artifact?" | S5 summary + per-query + chance check; S6 retrieval leaderboard | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a3_direction_robustness/task2_direction_robustness_audit.csv` | Keep corrected retrieval framework; stratify by `N_gallery`, `m_pos`, corrected null baseline, dataset, and cell line | Reviewer-defense critical | High | `asymmetry_artifact_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `A4` | Paired/statistical framework for core claims | Missing | Statistical packaging for the bridge analyses; answers "are main manuscript bridge claims backed by paired tests rather than raw means?" | `A1-A3` support tables; S4 target-level rows; S5 per-query; last-resort S1 internal-group export only if required | `docs/manuscript_core_paired_statistics.csv` | Execute after `A1-A3`; pair only on legitimate matched units; Wilcoxon primary, sign-test fallback, FDR by claim family | Required, statistical packaging | High | `denominator_fragility_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `B5` | Target-complexity analysis (`n_targets`) | Partial | Explanation layer for Result 3; answers "does C2G concordance vary with target complexity?" | Existing `docs/task2_c2g_query_n_targets_sensitivity.csv`; S5 per-query; S6 retrieval leaderboard | updated `docs/task2_c2g_query_n_targets_sensitivity.csv`; `docs/task2_n_targets_pairwise_stats.csv` | C2G only; exact `query_n_targets`; G2C reported as non-informative | Required, explanation layer | Medium | `denominator_fragility_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `B6` | Dose/time covariate-overlap audit | Missing | High-risk gate for Result 3; answers "is a controlled dose/time comparison actually supportable?" | S5 per-query; snapshot `delta_meta.csv` | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/b6_dose_time_overlap/task2_dose_time_overlap_audit.csv` | Audit exact `dataset`, `cell_line`, `time`, and `dose_value` overlap in evaluable query pools before any controlled comparison | Required, high-risk gate | High | `dose_time_confounding_challenge` | OA:self NAR:N EST:N QJ:Y MST:Y |
| `B7` | Matched-subset / controlled dose-time analysis | Blocked | Conditional explanation follow-up; answers "do dose/time-balanced subsets change the interpretation once overlap is proven?" | `B6` audit output; S5 per-query; S6 retrieval summaries | `docs/task2_dose_time_matched_analysis.csv` or limitation-only note | Run only if `B6` finds exact shared strata and support remains adequate after stratification | Required, explanation layer, conditional | High | `dose_time_confounding_challenge` | OA:Y NAR:Y EST:Y QJ:Y MST:Y |
| `B8` | Specificity local sensitivity | Partial | Local-only reviewer-defense note; answers "does scPerturb K562 C2G behavior shift by specificity tier?" | Existing `docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv`; S5 per-query; scPerturb `delta_meta.csv` | updated specificity table; `docs/task2_specificity_local_stats.csv` | scPerturb K562, C2G only; explicitly local-only; not required for the core main-text argument if `B5`, `B6`, `B10`, and `B9` already provide sufficient explanation | Reviewer-defense, local-only | Medium | `denominator_fragility_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `B10` | Dataset / cell-line composition analysis | Partial | Explanation-layer context; answers "are high and low concordance slices structurally biased by coverage or composition?" | `task2_pairs_coverage.csv`; representation registry; S4 and S5 attrition; S6 summaries | `docs/task2_dataset_cellline_composition_analysis.csv` | Treat LINCS as multi-cell-line and scPerturb as K562-local only; no raw dataset pooling; interpret `B9` only after or alongside this output | Required, explanation layer | Medium | `survivor_bias_enrichment_challenge` | OA:N NAR:Y EST:Y QJ:N MST:Y |
| `B9` | Enrichment for high-/low-concordance target sets | Missing | Explanation layer for Result 3; answers "are high and low concordance targets biologically structured after controlling for support and composition?" | S4 target-level Task2; S6 summaries; `B10` composition output; matched background registry; target annotations | `docs/task2_target_concordance_rankings.csv`; `docs/task2_high_low_concordance_enrichment.csv` | Build matched and evaluable background first; read only after or alongside `B10`; stop if coverage or composition context is inadequate | Required, explanation layer | High | `survivor_bias_enrichment_challenge` | OA:N NAR:Y EST:Y QJ:N MST:Y |
| `C11` | Local FM comparison within supported scPerturb K562 scope | Partial | Candidate local main-text comparison for Result 4; answers "within supported scope, do FM choices change the local Task2 reading?" | S6 group and retrieval leaderboards; panel ledger; availability registry; attrition tables | `docs/task2_fm_k562_local_comparison_summary.csv` | scPerturb K562 only; Gene and Pathway anchors retained; no benchmark-wide FM generalization | Candidate main-text, local-only | Medium | `local_fm_outlier_challenge` | OA:N NAR:Y EST:Y QJ:N MST:Y |
| `C12` | Per-target FM robustness / paired nonparametric testing | Missing | Robustness gate for any stronger interpretation of `C11`; answers "are FM differences target-stable or just aggregate shifts?" | S4 target-level rows; S5 local per-query rows; `C11` summary | `docs/task2_fm_k562_per_target_robustness.csv`; `docs/task2_fm_k562_paired_tests.csv` | Common evaluable local target and query pool only; paired tests only; stronger `C11` interpretation should not proceed without this layer | Required before stronger `C11` interpretation | High | `local_fm_outlier_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `D15` | Modality-good but mechanism-bad cases | Missing | Case-layer illustration of the external-validity gap; answers "where is Task1 bridge support present but Task2 support weak on the matched slice?" | `A1-A2` bridge tables; `A4` stats; S4 and S5 rows; Task1 support tables | `docs/manuscript_case_shortlist.csv` with `case_class=modality_good_mechanism_bad` | Run only after `A1/A2`, `A4`, and denominator constraints are settled; illustration only, not claim replacement | Case layer, explanatory illustration | High | `case_selection_bias_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `D13` | Representative aligned cases | Missing | Case-layer illustration; answers "what does a clean aligned case look like once the quantitative bridge is established?" | `A1-A2` bridge tables; `A4` stats; S4 and S5 rows; snapshot metadata | same shortlist, `case_class=aligned` | Same start conditions as `D15`; shortlist only from audited common-scope rows with denominator fields | Case layer, explanatory illustration | Medium | `case_selection_bias_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `D14` | Representative mechanism-mismatch cases | Missing | Case-layer illustration; answers "what does a mechanism mismatch look like when support is adequate?" | `A1-A2` bridge tables; `A4` stats; S4 and S5 rows; snapshot metadata | same shortlist, `case_class=mechanism_mismatch` | Same start conditions as `D15`; illustration only after quantitative bridge tables exist | Case layer, explanatory illustration | Medium | `case_selection_bias_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |
| `D16` | Representation-switch cases if present | Conditional | Optional FM case-layer illustration; answers "are there local cases where FM reverses the representation reading?" | `C11-C12` outputs; S4 and S5 local rows | same shortlist, `case_class=representation_switch` | Run only if paired FM tests show real local reversals after support filtering and after the other case-layer prerequisites are met | Optional exploratory | Medium | `case_selection_bias_challenge` | OA:N NAR:Y EST:Y QJ:Y MST:Y |

## 4. Document Change Plan

- Create `docs/manuscript_analysis_expansion_plan.md` as the execution-planning document for this phase.
- Create optional `docs/manuscript_analysis_file_change_inventory.csv` as an analysis-to-file matrix.
- No update planned for `docs/manuscript_proposal.md` because it is absent locally; note in the plan that `docs/manuscript_blueprint.md` is the canonical revised-proposal replacement.
- Update `docs/manuscript_blueprint.md` to clarify the bridge layer, explanation layer, local FM layer, and case layer.
- Update `docs/manuscript_evidence_extraction_plan.md` with extraction packages in the same phase order used in Section 7.
- Update `docs/manuscript_source_register.md` to register new support tables, composition-before-enrichment logic, and last-resort S1 export provenance if ever used.
- Update `docs/manuscript_claim_to_evidence_map.csv` with bridge-table rows, local-only cautions, and stats-table references.
- Update `docs/manuscript_blocker_memo.md` so blocked items become gated analyses with explicit start and stop conditions.
- Update `docs/manuscript_panel_support_ledger.csv` to carry direction-robustness notes, composition-before-enrichment notes, and local-only labels.
- Create `docs/manuscript_support_table_specs.md` to define keys and schemas for every new manuscript-only support table.

## 5. Code Change Plan

- Inspect `scripts/manuscript_support_remediation.py`: reuse helpers where possible, but keep downstream manuscript-support workflows as the default path.
- Inspect `scripts/s2_task1_cross_metrics.py`: no semantic change planned; use existing `task1_group_cross.parquet` as the cross-side anchor.
- Inspect `scripts/s5_task2_retrieval_multisource.py`: no stage edit planned because current per-query fields already support the needed downstream joins.
- Create `scripts/manuscript_analysis_support.py`: build `A1-A4`, `B5-B8`, `B10`, `C11-C12`, and `D13-D16` support tables from audited S1, S2, S4, S5, and S6 outputs only.
- Create `scripts/manuscript_enrichment_support.py`: build the matched background registry, per-target concordance ranking table, and enrichment outputs for `B9`, with `B10` composition context read in by default.
- Conditional last-resort option: touch `scripts/s1_task1_internal_metrics.py` only if existing S1 outputs are proven insufficient for `A1/A2/A4`; if needed, add `task1_group_internal.parquet` as an audited export of already-computed rows rather than a new benchmark stage.
- Upstream reruns: none by default; rerun S1 only if the last-resort export is approved and still necessary after exhausting downstream options.
- Validation: enforce unique pair keys, denominator conservation against source tables, explicit `n_pairs`, `median_delta`, `p`, and `FDR` in stats tables, matched-background equality checks for enrichment, and automatic local-only labels on specificity and FM outputs.

## 6. Output / Table Plan

- High-risk gates: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a3_direction_robustness/task2_direction_robustness_audit.csv`, `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/b6_dose_time_overlap/task2_dose_time_overlap_audit.csv`.
- Main bridge analyses: `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a1_task1_internal_cross/task1_internal_vs_cross_common_scope_comparison.csv`, `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a2_task1_task2_bridge/task1_task2_shared_scope_bridge_table.csv`.
- Statistical packaging: `docs/manuscript_core_paired_statistics.csv`.
- Biological / explanation layer: updated `docs/task2_c2g_query_n_targets_sensitivity.csv`, `docs/task2_n_targets_pairwise_stats.csv`, `docs/task2_dataset_cellline_composition_analysis.csv`, `docs/task2_target_pool_background_registry.csv`, `docs/task2_target_concordance_rankings.csv`, `docs/task2_high_low_concordance_enrichment.csv`, conditional `docs/task2_dose_time_matched_analysis.csv`, updated `docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv`, `docs/task2_specificity_local_stats.csv`.
- Local representation-effect layer: `docs/task2_fm_k562_local_comparison_summary.csv`, `docs/task2_fm_k562_per_target_robustness.csv`, `docs/task2_fm_k562_paired_tests.csv`.
- Case-study layer: `docs/manuscript_case_shortlist.csv` with `case_class` values `modality_good_mechanism_bad`, `aligned`, `mechanism_mismatch`, and `representation_switch` when present.

## 7. Priority Ranking

Phase 1: high-risk gates

- `A3` C2G vs G2C direction-specific robustness audit.
- `B6` dose/time covariate-overlap audit.

Phase 2: main bridge analyses

- `A1` Task1 internal vs cross common-scope comparison.
- `A2` Task1 vs Task2 shared-slice bridge comparison.

Phase 3: statistical packaging

- `A4` paired/statistical framework for the bridge analyses.

Phase 4: biological / explanatory layer

- `B5` target-complexity analysis.
- `B10` dataset / cell-line composition analysis.
- `B9` enrichment for high-/low-concordance target sets.
- `B7` matched-subset / controlled dose-time analysis only if `B6` passes.
- `B8` specificity local sensitivity.

Phase 5: local representation-effect layer

- `C11` local FM comparison within supported scPerturb K562 scope.
- `C12` per-target FM robustness / paired testing.

Phase 6: case-study layer

- `D15` modality-good but mechanism-bad cases.
- `D13` representative aligned cases.
- `D14` representative mechanism-mismatch cases.
- `D16` representation-switch cases only if supported.

## 8. Stop Rules

- Do not run `B7` if `B6` finds no exact shared `dataset`, `cell_line`, `time`, and `dose_value` strata in the evaluable query pool; emit only the overlap audit table plus a limitation statement.
- Do not interpret `B9` without a matched and evaluable background registry and `B10` composition context.
- Do not promote any local-only analysis (`B8`, `C11`, `C12`, `D16`) to benchmark-wide language.
- Do not touch S1 upstream by default; only consider the audited internal-group export after proving that existing S1 outputs cannot support `A1/A2/A4` downstream.
- Do not keep a main bridge or statistical packaging claim (`A1`, `A2`, `A4`) if the required matched-unit support is absent and the fallback would reduce the claim to raw summary means.
- Do not run enrichment (`B9`) without a matched and evaluable target-pool background registry and acceptable annotation coverage; if unavailable, stop at background and composition reporting only.
- Do not start the case layer (`D13-D16`) until `A1/A2` are complete, `A4` is available, and major denominator and support constraints have already been resolved.
- Do not keep a stratified analysis when denominator support collapses after filtering; downgrade it to reviewer-defense, supplementary, or limitation-only.
- Stop if a proposed step only improves cosmetic symmetry or broadens scope without strengthening manuscript logic; this phase is manuscript support work, not a second benchmark-expansion project.
- Do not use the whole genome as an enrichment background, do not force a dose/time analysis without overlap, and do not infer benchmark-wide conclusions from K562-local FM or specificity results.
