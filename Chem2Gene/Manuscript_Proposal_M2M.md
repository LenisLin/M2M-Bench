# M2M-Bench Manuscript Proposal (Draft for Co-author Review)

**Version:** v0.4 (2026-02-14)  
**Project:** M2M-Bench (*From Modality Concordance to Mechanism Fidelity*)  
**Goal of this document:** provide a concrete writing/figure proposal with real values from current outputs.

---

## 1) Proposed Results Structure (8 sections, 4 main figures)

### R0. Motivation (no dedicated main figure; introduced with Fig1)
- Core problem: current perturbation studies often mix **modality gap** (bulk vs scRNA-seq) with **mechanism gap** (drug vs gene).
- Our benchmark logic: first establish modality validity (Task1), then quantify mechanism fidelity (Task2), then stress-test representation robustness (Task3).

### R0b. Claim hierarchy + analysis units (Fig1)
- **Primary contribution:** reusable audit-and-evaluation standard (not a new predictor).
- **Primary finding:** cross-modality concordance is limited and context-dependent.
- **Analysis units explicitly defined:** `row`, `context`, `matched pair`, `query`, `Task2 context`, `class`.

### R1. Study design + data composition (**Fig1**)
- Global data composition and attrition from task0→task1.
- Why available matched contexts are sparse and biased.

### R2. Task1 — Modality concordance (**Fig2**)
- Group-wise gap, instance-level retrieval, and confounder decomposition.
- Set-level centroid analysis as robustness check.

### R3. Task2 — Mechanism fidelity in same-domain only (**Fig3**)
- Context-level robustness classes and protocol sensitivity.
- Evidence-driven case studies for robust vs fragile mechanism transfer.

### R4. Task3 — Foundation model stress test (**Fig4**)
- Does representation quality align with biological fidelity?
- Direction/view-dependent model behavior.

### R5. Impact and implications (discussion-only, **no figure**)
- Practical recommendations for virtual cell / perturbation prediction ecosystem.
- Reporting standards and interpretation boundaries.

### R6. Current empirical snapshot (2026-02-14 rerun)
- Task1 matched pairs: full `670`, strict protocol subset `313` (`46.7%`), strict subset entirely `Genetic / exact_cond`.
- Task1 global cosine: full gene/path `0.0238 / 0.0599`; strict gene/path `0.0060 / 0.0527`.
- Task1 retrieval after balanced evaluation (fixed gallery=256, true=1) drops markedly for chemical groups:
  - Chemical mean MRR (across tracks/directions): `0.598 -> 0.104`
  - Chemical mean Top1 (across tracks/directions): `0.489 -> 0.062`
- Task2 class composition remains: `Robust_High=567`, `Intermediate=1020`, `Robust_Low=1020`, `Protocol_Sensitive=454`.
- Task3 remains direction/view-dependent with no universal best model.

### R7. Reviewer-critical robustness add-ons (new)
- Task1 retrieval now reported as **raw + balanced + random baseline**.
- Task1 now includes **cross-vs-within effect-size calibration**.
- Task1 now includes **continuous protocol sensitivity** and **leave-HEPG2-out** sensitivity.
- Task2 now includes **target-confidence / polypharmacology stratification**.
- Task3 now includes **model meta-analysis** (`size_class`, `perturbation_trained`).

---

## 2) Main Figures (Fig1–Fig4) with panel-level plan

## Fig1. Study design + data composition

**Scientific objective:** establish why a modality-first benchmark is necessary and what data support is actually available.

### Panel plan
- **Fig1A (design schematic):** M2M task flow: Task1 (modality) → Task2 (mechanism) → Task3 (FM stress test).
- **Fig1A2 (analysis-unit box):** define `row/context/matched pair/query/class` and task-specific context keys.
- **Fig1B (attrition funnel):** `2,175,012 rows / 108,713 contexts` → `17,078 rows / 84 contexts` → `670 matched pairs`.
- **Fig1C (overlap composition):**
  - Chemical overlap contexts: `2` (all A549: `HSP90AA1`, `SIRT1`)
  - Genetic overlap contexts: `82` (HEPG2 `73/82 = 89.0%`, JURKAT `7`, K562 `2`)
- **Fig1D (context drop reasons):** dominant losses are `cell_and_target_absent_in_other_source` and `target_not_shared_in_other_source`.
- **Fig1E (matching protocol gaps):**
  - Chemical matched by nearest dose/time: dose log-diff mean `1.014`, time abs-diff mean `124.5h`
  - Genetic matched under exact condition (`0` dose/time gap)
  - Report strict protocol subset (dose_logdiff ≤ 0.5, time_absdiff ≤ 24h) in supplement.
- **Fig1F (benchmarkability map):**
  - `Z0_NoCrossSourceOverlap`: Chemical `16,761`, Genetic `91,868` contexts
  - `Z1_ChemicalNearestOnly`: `2` contexts
  - `Z3_GeneticExactComparable`: `82` contexts

### Required data files
- `outputs/task1_audit/analysis/stage_summary.csv`
- `outputs/task1/qc/tables/qc_context_overlap_summary.csv`
- `outputs/task1/qc/tables/qc_overlap_contexts_by_cell.csv`
- `outputs/task1_audit/analysis/context_drop_reasons_summary.csv`
- `outputs/task1_audit/analysis/protocol_gap_summary.csv`
- `outputs/task1_audit/analysis/match_type_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/benchmarkability_zone_summary.csv`

---

## Fig2. Task1 modality concordance

**Scientific objective:** quantify cross-modality agreement and identify where apparent alignment is real vs structural artifact.

### Panel plan
- **Fig2A (group-wise gap summary):**
  - Global mean cosine: gene `0.0238`, pathway `0.0599`
  - FDR-significant vs null for cosine metrics, but absolute concordance remains low.
- **Fig2A2 (strict protocol subset contrast):**
  - Strict subset (`n=313`): global mean cosine gene/path `0.0060 / 0.0527`
  - Show full vs strict contrast to quantify protocol-mismatch effect.
- **Fig2B (retrieval performance heatmap):**
  - Chemical `LINCS→sc` gene: MRR `0.9009`, Top1 `0.8599`
  - Genetic `LINCS→sc` gene: MRR `0.0941`, Top1 `0.0479`
  - Show both directions and gene/path tracks.
- **Fig2C (balanced retrieval control):**
  - Report balanced-candidate metrics (fixed gallery size + fixed positives per query).
  - Use `balanced_gallery_size=256` and `balanced_true_per_query=1` (default), with repeated subsampling.
  - Highlight raw vs balanced drop (chemical mean MRR `0.598 -> 0.104`).
  - Add explicit random baseline:
    - theoretical `MRR_random=H256/256=0.0239`, `Top1_random=1/256=0.003906`
    - balanced chemical still above random (`MRR 0.104`, `4.36×` random), but far below raw.
- **Fig2D (context contribution / sample-structure explanation):**
  - Chemical average true-ratio in gallery: `0.2618` vs Genetic `0.0123`
  - Example contrast:
    - `A549||Chemical||HSP90AA1`: true-ratio `0.9188`, sc→LINCS MRR `1.000`
    - `A549||Chemical||SIRT1`: true-ratio `0.0812`, sc→LINCS MRR `0.278`
- **Fig2E (effect-size calibration):**
  - Cross vs within replicate consistency:
    - Chemical gene: within-avg `0.125` vs cross-matched `0.039` (ratio `0.314`)
    - Genetic gene: within-avg `0.055` vs cross-matched `0.006` (ratio `0.109`)
  - Interpretation: statistical significance should be interpreted against within-modality consistency ceilings.
- **Fig2F (continuous protocol sensitivity):**
  - Chemical Spearman trends:
    - `time_absdiff` vs `cosine_gene`: rho `-0.255`, p `<1e-5`
    - `time_absdiff` vs `cosine_path`: rho `-0.250`, p `<1e-5`
    - `dose_logdiff` vs `cosine_gene`: rho `+0.169`, p `0.0013`
  - Show line curves (`dose/time` mismatch bin centers vs mean cosine), not cutoff-only comparison.
- **Fig2G (composition-bias sensitivity):**
  - Pairwise ALL vs leave-HEPG2:
    - ALL gene/path `0.0238/0.0599`
    - leave-HEPG2 gene/path `0.0377/0.0630`
  - Genetic-only flips after leave-HEPG2 (`0.0060/0.0527 -> -0.0087/-0.0275`).
- **Fig2H (confounder effects):**
  - For gene-gap: dose/time/match score/target/cell effects significant in permutation tests.
  - Variance decomposition indicates target dominates explained variance.
- **Fig2I (set-level centroid reliability):**
  - Global set-level cosine (unweighted): gene `0.0457`, pathway `0.1555`
  - Case examples:
    - High: `HEPG2-Genetic-HSPA5` gene `0.401`, path `0.748`
    - Low: `HEPG2-Genetic-STIL` gene `-0.111`, path `-0.292`

### Required data files
- `outputs/task1/analysis/modality_gap_summary.csv`
- `outputs/task1/analysis/strict/modality_gap_summary_strict.csv`
- `outputs/task1/analysis/significance_results.csv`
- `outputs/task1/retrieval/analysis/retrieval_summary.csv`
- `outputs/task1/retrieval/analysis/retrieval_per_query.csv`
- `outputs/task1/retrieval/analysis/retrieval_null_summary.csv`
- `outputs/task1/qc/tables/qc_strict_subset_counts.csv`
- `outputs/task1/confounder/analysis/factor_permutation_tests.csv`
- `outputs/task1/confounder/analysis/variance_contribution_gap_gene.csv`
- `outputs/task1/set_level_allctx/analysis/set_level_summary.csv`
- `outputs/task1/set_level_allctx/analysis/set_level_context_metrics.csv`
- `outputs/task1_reviewer_fixes/analysis/retrieval_dual_report.csv`
- `outputs/task1_reviewer_fixes/analysis/effect_calibration_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/protocol_continuous_spearman.csv`
- `outputs/task1_reviewer_fixes/analysis/protocol_continuous_curves.csv`
- `outputs/task1_reviewer_fixes/analysis/leave_one_cell_out_pairwise.csv`

---

## Fig3. Task2 mechanism fidelity (same-domain only)

**Scientific objective:** identify when drug↔gene mechanism transfer is robust, weak, or protocol-sensitive under controlled domain scope.

### Panel plan
- **Fig3A (global statistical tests):**
  - Same-domain-only primary analysis.
  - Systema vs Standard not significant at global level (gene/path tracks).
  - Pathway vs gene differences statistically detectable but small effect.
- **Fig3B (context class composition):**
  - `Robust_High = 567`
  - `Intermediate = 1020`
  - `Robust_Low = 1020`
  - `Protocol_Sensitive = 454`
- **Fig3C (representative context cases):**
  - Strong examples (N≥20): `YAPC-PSMB1 0.9986 (n=636)`, `HT29-PSMB1 0.9979 (n=723)`, `HT29-MAP2K1 0.9760`
  - Weak examples (N≥20): `A375-PYGM 0.0960`, `HT29-PYGM 0.1102`
- **Fig3D (protocol sensitivity and enrichment):**
  - Protocol-sensitive contexts tested: `429`
  - Significant dose associations: `36`; significant time associations: `5`
  - Target enrichments:
    - Robust_High: `MTOR`, `PSMB1`, `MAP2K1`
    - Protocol_Sensitive: `NR3C1`, `HSP90AA1`
- **Fig3E (target confidence / polypharmacology tiers):**
  - `Tier_A_HighConfidence`: `84` contexts (`6` targets), Robust_High fraction `52.4%`
  - `Tier_C_AmbiguousOrPolypharm`: `31` contexts (`3` targets), Protocol_Sensitive fraction `48.4%`
  - Enrichment (FDR<0.05):
    - Robust_High enriched in Tier_A (log2OR `+2.37`)
    - Protocol_Sensitive enriched in Tier_C (log2OR `+2.47`)

### Required data files
- `outputs/task2_nodomain/run_manifest_task2_nodomain.json`
- `outputs/task2_nodomain/analysis/Step1_Tests_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step2_Context_Labels_NoDomain.csv`
- `outputs/task2_nodomain/analysis/task2_top_context_examples_n20.csv`
- `outputs/task2_nodomain/analysis/Step5_Protocol_Correlations_ProtocolSensitive_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step5_Enrichment_Targets_RobustHigh_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step5_Enrichment_Targets_ProtocolSensitive_NoDomain.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_summary.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_enrichment.csv`

---

## Fig4. Task3 foundation model stress test

**Scientific objective:** test whether model representation quality is stable and biologically faithful across directions/views.

### Panel plan
- **Fig4A (retrieval scoreboard by direction/view):**
  - No universal winner; top model depends on direction/view.
  - Standard CRISPR→Drug top: `PCA200 (MRR 0.2022)`
  - Standard Drug→CRISPR top: `scBERT (MRR 0.2358)`
  - Systema Drug→CRISPR top: `UCE (MRR 0.2266)`
- **Fig4B (pairwise centroid cosine scoreboard):**
  - Standard top: `UCE (0.1462)`
  - Systema top: `scGPT (0.1878)` and `UCE (0.1866)` close second
- **Fig4C (lollipop: Systema − Standard, retrieval CRISPR→Drug):**
  - Significant positive gains in multiple tracks (e.g., UCE, PCA100/200/50, TahoeX1_3b).
  - Non-significant/weak gains for some models (e.g., scBERT, Geneformer, scFoundation).
- **Fig4D (quality control / consistency):**
  - Task3 audit checks passed: `18/18`.
  - Ensure figure claims are reproducible and not artifact of missing/filtered rows.
- **Fig4E (model meta-analysis):**
  - Add model-level aggregate (`mean_scaled_best`) with metadata (`size_class`, `perturbation_trained`).
  - Current top overall: `scGPT (0.684)`, followed by `PCA200 (0.647)`.
  - No robust universal perturbation-trained advantage across slices (best p≈`0.058`, non-significant).
  - Size effects are metric-dependent (e.g., pairwise `Mean_NegEDist` shows strong positive trend; retrieval trends mixed).

### Required data files
- `/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Analysis/Task3_Retrieval_Summary.csv`
- `/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Analysis/Task3_Pairwise_Summary.csv`
- `/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified/viz_fig2_lollipop.csv`
- `outputs/task3_audit/analysis/task3_consistency_checks.csv`
- `outputs/task3_audit/analysis/task3_key_summary.csv`
- `outputs/task3_meta/analysis/task3_model_scoreboard_meta.csv`
- `outputs/task3_meta/analysis/task3_meta_perturbation_training_tests.csv`
- `outputs/task3_meta/analysis/task3_meta_size_trend_tests.csv`

---

## 3) Supplementary Figures (8 figures total: S1–S8)

### S1. Detailed attrition accounting
- Stage-wise row/context loss tables.
- Source-specific kept vs dropped at stage3.
- Benchmarkability zone map table and counts.

### S2. Context overlap atlas
- Cell × modality overlap heatmaps.
- Full target lists for chemical/genetic overlaps.

### S3. Matching quality diagnostics
- Match score distribution.
- Dose/time mismatch distributions and scatter.
- Strict protocol subset counts and summary statistics.

### S4. Task1 inference robustness
- Null distribution summaries.
- Bootstrap CIs and permutation details by strata.
- Cross-vs-within effect-size calibration tables.
- Continuous protocol sensitivity curves and 2D dose-time grids.
- Leave-one-cell-out sensitivity tables.

### S5. Task1 retrieval deep-dive
- Per-query and per-context retrieval distributions.
- Candidate-space snapshots and 2D embedding plots.
- Balanced retrieval evaluation (fixed gallery size / fixed positives).

### S6. Task1 set-level analyses
- Context-level centroid metrics full table.
- Set-level vs pairwise comparison scatter/hist.

### S7. Task2 deep diagnostics
- Step4 tracer full outputs.
- Enrichment tables (targets/cells) and protocol-correlation details.
- Target confidence / polypharmacology tier mapping and enrichment.

### S8. Task3 full scoreboard and audit
- Full track/view/direction scoreboards.
- Supplementary lollipop comparisons.
- File-level manifest and audit report.
- Model metadata, size-class trends, and perturbation-training tests.

---

## 4) R5 Discussion Proposal (no figure)

### 4.1 Key claims we can defend with current evidence
- Cross-modality comparability is limited and uneven; modality effects cannot be ignored.
- Cross-modality cosine should be interpreted against within-modality ceilings; cross is a small fraction of within replicate consistency.
- Mechanism transferability is context-dependent, not globally stable.
- Protocol effects (especially dose) are a major source of instability.
- FM representation quality is task-direction dependent; no single model dominates all objectives.

### 4.2 Implications for virtual cell / perturbation prediction
- Claims of mechanism generalization should be gated by modality-valid contexts.
- Reporting should be context-level (not only global averages).
- Benchmark outputs should include coverage, protocol, and confounder diagnostics by default.

### 4.3 Boundaries and caution
- Current overlap is sparse (especially chemical contexts), so conclusions should explicitly state data-coverage limits.
- High retrieval in small/imbalanced contexts should not be overinterpreted as broad biological concordance.
- Target-tier analysis currently uses a compact heuristic map (can be replaced by curated pharmacology annotations in final submission).

---

## 5) Risk Register and Mitigations (summary for reviewers)

1. **Sparse chemical overlap (2 contexts).**  
   Mitigation: constrain claims to diagnostic conclusions; emphasize per-context reporting and avoid global chemical generalization.
2. **Composition bias (HEPG2 dominates genetic overlap).**  
   Mitigation: implemented cell-stratified and leave-one-cell-out sensitivity tables (including leave-HEPG2-out).
3. **Protocol mismatch (dose/time gaps).**  
   Mitigation: implemented strict subset + continuous dose/time sensitivity curves.
4. **Metric susceptibility to systematic variation.**  
   Mitigation: implemented balanced retrieval + null-calibrated + theoretical-random baselines.
5. **Candidate-space imbalance inflating retrieval.**  
   Mitigation: fixed gallery size + fixed positives per query (balanced retrieval).
6. **Drug–gene MoA ambiguity.**  
   Mitigation: implemented target-confidence / polypharmacology tier stratification and enrichment reporting.
7. **Readout mismatch (bulk vs sc).**  
   Mitigation: report platform-specific baselines and avoid over-claiming modality gaps as biology.
8. **Same-domain restriction.**  
   Mitigation: justify as first-order identifiable setting; cross-domain remains exploratory.

---

## 6) Immediate next actions (for co-author review round)

1. Confirm Fig1 unit-definition box and Fig1 benchmarkability map design.
2. Freeze Fig2 reviewer-critical panels (balanced random baseline, calibration, continuous sensitivity, LOO).
3. Replace heuristic target-tier map with curated annotation file if available.
4. Draft Result section text directly against this updated panel map.
