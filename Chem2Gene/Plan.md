# Project Proposal: M2M-Bench (Modality-to-Mechanism Benchmark) v2.4.0

**Former name:** Chem2Gen-Bench  
**Full subtitle:** *From Modality Concordance to Mechanism Fidelity: A Unified Perturbation Benchmark*

---

# Part 1. Motivation and Scientific Questions

This project addresses a core translational problem in perturbation biology:

> When we observe transcriptomic differences, are we measuring true mechanism differences, or assay/modality artifacts?

We focus on two linked scientific questions:

1. **Modality Concordance Question (Q1):**  
   For the same cell line under the same perturbation condition, how different are profiles between **LINCS bulk** and **scPerturb single-cell**?
2. **Mechanism Fidelity Question (Q2):**  
   After controlling key confounders, how faithfully does a **drug perturbation** phenocopy its matched **genetic perturbation**?

Why this order matters:

- If Q1 is unstable, Q2 can be over- or under-estimated.
- Therefore Q1 is the validity layer; Q2 is the biological mechanism layer.

---

# Part 2. Unified Framework (Already Established in Chem2Gene)

The project keeps and extends your existing framework with three analysis engines:

1. **Group-wise comparison** (distribution-level effect size and significance)
2. **Instance-level retrieval** (query-level recoverability of true identity/target)
3. **Confounder analysis** (cell, source, dose/time, tier, protocol sensitivity)

This same framework is applied consistently to Q1 and Q2.

## 2.1 Claim hierarchy (locked for manuscript)

- **Primary contribution (what this paper contributes):**
  - a reusable **evaluation-and-audit standard** for modality→mechanism benchmarking,
  - including coverage audit, protocol-mismatch diagnostics, candidate-bias controlled retrieval, and confounder-aware interpretation.
- **Primary empirical finding (what the benchmark finds):**
  - cross-modality concordance is low and strongly context-dependent under current overlap.
- **Secondary findings:**
  - mechanism fidelity is context/tier dependent;
  - FM performance is direction/view dependent without a universal winner.

## 2.2 Analysis-unit dictionary (for Fig1 + Methods)

- `row`:
  - one perturbation profile in unified metadata (`unified_meta`), tied to one tensor pointer.
- `context` (Task1):
  - one tuple `(cell_std, modality, target_std)`.
- `candidate row` (Task1):
  - one `row` inside an overlap context where both LINCS and scPerturb exist.
- `matched pair` (Task1):
  - one LINCS row ↔ one scPerturb row matched within a context by `exact_cond` (genetic) or nearest dose/time (chemical).
- `query` (Task1 retrieval):
  - one candidate row evaluated against opposite-source gallery.
- `instance` (Task2 Step1):
  - one retrieval/pairwise observation keyed by `(Track, Scenario, View, Direction, Cell, Target, Source_Chem, Source_Gene, Dose, Time, CondID)`.
- `context` (Task2 Step1 aggregated):
  - one tuple `(Track, View, Cell, Target, Source_Chem, Source_Gene)`.
- `labeled_context` (Task2 Step2):
  - one tuple `(Cell, Target)` with `Mean_Success`, `Peak_Success`, and `Performance_Class`.
- `class` (Task2):
  - one context label in `{Robust_High, Intermediate, Robust_Low, Protocol_Sensitive}`.

## 2.3 Primary endpoints + non-claims (submission lock)

- **Task1 primary endpoints:**
  - pairwise `mean cosine_gene` with cross-vs-within calibration,
  - balanced retrieval `MRR` (gallery=256, true=1) with random baseline.
- **Task2 primary endpoints:**
  - class composition on `labeled_context` table,
  - protocol-sensitive association counts.
- **Task3 primary endpoints:**
  - retrieval `Mean_MRR`,
  - pairwise `Mean_CentroidCosine`.
- **Non-claims to state explicitly:**
  - no broad chemical generalization beyond sparse overlap,
  - no global drug↔gene mechanism equivalence,
  - no universal best FM.

---

# Part 3. Data Scope and Representation

## 3.1 Data assets

- **LINCS**: bulk perturbation transcriptomics
- **scPerturb**: single-cell perturbation transcriptomics
- **Matched K562 gold set** (Task 3 model stress-testing)

## 3.2 Harmonized feature space

- Shared curated gene universe (~2.5k genes)
- Standardized metadata fields (e.g., `cell_std`, `target_std`, `modality`, `source_db`, dose/time fields)
- Canonical field semantics:
  - `target_std`: **target gene label only** (standardized string), not including treatment type
  - `modality` / treatment type: perturbation class (`Chemical` or `Genetic`)
  - Recommended composite key when needed: (`cell_std`, `modality`, `target_std`)

## 3.3 Subset design

- **L1:** same cell + same target (mechanism gold pairs)
- **L2:** different cell + same target (target generalization pool)
- **L3:** same cell + different target (background pool for Systema-style debiasing)
- **M1 (new):** same cell + same perturbation identity + matched conditions across LINCS/scPerturb (modality benchmark core set)

## 3.4 Signature views

1. **Standard view**: raw differential signatures
2. **Systema view**: background-subtracted signatures for stress/systematic component control

## 3.5 Canonical data structure (Task 0 outputs)

Task 0 (`scripts/task0_curate_data.py`) writes a stable structure:

- `outputs/task0_curated/metadata/`
  - `unified_meta.parquet`
  - `unified_meta.csv` (parquet-free fallback)
  - `level1_pairs_meta.parquet`
  - `level1_pairs_unique.parquet`
  - `level1_pair_stats.parquet`
  - `level2_chem_pool.parquet`, `level2_gene_pool.parquet`, `level2_pair_stats.parquet`
  - `level3_chem_pool.parquet`, `level3_gene_pool.parquet`, `level3_pair_stats.parquet`
- `outputs/task0_curated/tensors/Level1_Tensors_Split/`
  - `L1_chem_x.pt`, `L1_chem_y_gene.pt`, `L1_chem_y_path.pt`
  - `L1_gene_x.pt`, `L1_gene_y_gene.pt`, `L1_gene_y_path.pt`
  - `manifest_level1_tensors.json`
- `outputs/task0_curated/bundle/m2m_task0_bundle.pt` (optional unified object)
- `outputs/task0_curated/run_manifest_task0.json` (run config + summary + output pointers)

Core metadata fields used across tasks:

- identity/context: `cell_std`, `target_std`, `modality`, `source_db`, `pair_id`
- technical pointers: `uid`, `global_idx`, `chunk_file`, `chunk_idx`
- protocol/context covariates: `dose_val`, `time_val`, `tissue`, `disease`, `cond_id`

Important identity note:

- `pair_id = cell_std:target_std` is gene-centered and does **not** include treatment type.
- Treatment type should be carried/stratified explicitly via `modality` (or `modality_chem` / `modality_gene` in paired tables).

Tensor semantics:

- `x`: baseline representation
- `y_gene`: gene-level perturbation delta
- `y_path`: pathway-level perturbation delta

---

# Part 4. Three-Task Project Structure

## 4.0 Task-wise data usage (what data to use)

### Task 1 — Modality Concordance (Q1)

Use:

- `unified_meta.parquet` as the source table for cross-source matching.
- A Task 1 matched table (**M1**) derived from `unified_meta` with keys:
  - same `cell_std`,
  - same `modality` (treatment type),
  - same `target_std` (gene label),
  - source contrast fixed as LINCS vs scPerturb,
  - best-available protocol matching (`dose_val`, `time_val`, `cond_id`).
- Raw deltas from source tensors/chunks for matched instances.

Do not use as primary input:

- L2/L3 pools (these are not the primary modality concordance unit).

### Task 2 — Mechanism Fidelity (Q2)

Use:

- `level1_pairs_meta.parquet` and `level1_pairs_unique.parquet` for gold pair logic.
- L1 split tensors (`L1_chem_*`, `L1_gene_*`) for pairwise/retrieval core metrics.
- `level2_*` pools for target generalization diagnostics.
- `level3_*` pools for Systema-style background estimates.
- `level*_pair_stats.parquet` for QC and denominator checks.

Evaluation unit:

- For each (`cell_std`, `target_std`), compare **Chemical vs Genetic** effects (treatment type separated by design).

### Task 3 — Foundation Model Stress Test

Use:

- Matched K562 evaluation-set data (Task 3 specific, independent from Task 0 pools).
- FM embedding outputs from `8_Task3_*.py`.
- Cell-level delta artifacts from `9_Task3_PrepareAnalysis.py`.
- Analysis/report tables from `10_Task3_AnalyzeOutputs.py` and `11_Task3_OutputCollect.py`.

Target filtering rule:

- Keep targets that appear in **both** modalities on K562 (at least one Drug perturbation and one Genetic perturbation for the same target gene), with controls retained.

Interpretation rule:

- Task 3 model conclusions must be interpreted against Task 2 stability classes and Task 1 modality reliability context.

## Task 1 — Modality Concordance Benchmark (Q1)

**Goal:** quantify LINCS↔scPerturb agreement for the same perturbation context.

### 1A. Group-wise modality gap
- Compare matched pairs in M1 using `cosine`, `cosine_z`, `l2`, and `edist`.
- Build stratified null distributions and confidence intervals.
- Run strict protocol subset as a required sensitivity analysis:
  - Chemical: keep only pairs with controlled `dose_logdiff` and `time_absdiff`.
  - Genetic: keep `exact_cond` pairs.

### 1B. Cross-source retrieval
- LINCS query → scPerturb gallery and reverse.
- Evaluate Top-k, MRR, and failure strata.
- Report both raw retrieval and balanced retrieval (fixed gallery size + fixed positives).
- Always include random baselines:
  - empirical null from permutation,
  - theoretical balanced baseline (`Top1 = 1/N`, `MRR = H_N/N` when one positive).

### 1C. Confounder decomposition
- Quantify effect of cell context, target class, dose/time mismatch, source, and perturbation class.
- Label contexts as **Robust / Conditional / Unstable** for downstream interpretation.

### 1D. Reviewer-critical robustness extensions
- **Effect-size calibration:** compare cross-modality cosine against within-modality replicate cosine (LINCS-within, sc-within).
- **Continuous protocol sensitivity:** model cosine trends along continuous `dose_logdiff` and `time_absdiff` (not only strict cutoff).
- **Composition-bias sensitivity:** report cell-stratified and leave-one-cell-out summaries (including `leave-HEPG2-out`).
- **Benchmarkability map:** classify contexts into comparability zones (`NoOverlap`, `ChemicalNearest`, `GeneticExact`).

**Primary output:** modality reliability table per context (used as interpretation guardrail in Task 2).

---

## Task 2 — Chemical–Genetic Mechanism Fidelity (Q2; original Chem2Gene core)

**Goal:** measure and explain drug↔gene transferability after confounder-aware control.

**Primary analysis contract (updated):**

- Task 2 main hypothesis testing **does not use SameDomain vs CrossDomain as a primary contrast**.
- Domain/source fields are retained for diagnostics, provenance, and sensitivity analysis only.
- Main reporting axis is mechanism fidelity under matched chemical↔genetic contexts, not platform-domain stratification.

### 2A. Pairwise gap quantification
- Retain existing pairwise metrics pipeline (`3_Task1_Pairwise_Test.py`).
- Compare Standard vs Systema behavior.

### 2B. Multi-scenario retrieval
- Retain existing retrieval pipeline (`4_Task1_Retrieval_Test.py`).
- Scenarios: target retrieval, cell dependence, global specificity.

### 2C. Attribution and protocol sensitivity
- Retain/extend `5_Task2_Analysis_Attribution.py`.
- Produce Consistently_Good / Consistently_Bad / Conditionally_Good classes.
- Evaluate dose/time/context sensitivity and target-tier behavior.

### 2D. Target-confidence / polypharmacology stratification
- Add explicit target-tier annotation (`HighConfidence`, `ModerateOrUnknown`, `AmbiguousOrPolypharm`).
- Report class composition and enrichment under tier stratification.

**Interpretation rule:** report mechanism conclusions jointly with Task 1 modality reliability.

---

## Task 3 — Foundation Model Stress Test

**Goal:** test whether FM latent spaces preserve biological fidelity instead of hallucinated alignment.

### 3A. Embedding extraction
- Existing scripts: `8_Task3_scGPT.py`, `8_Task3_geneformer.py`, `8_Task3_scBERT.py`, `8_Task3_UCE.py`, `8_Task3_STATE.py`, `8_Task3_Tahoe-x1.py`, `8_Task3_scFoundation.py`

### 3B. Delta construction and evaluation
- `9_Task3_PrepareAnalysis.py` for cell-level deltas
- `10_Task3_AnalyzeOutputs.py` for retrieval + pairwise metrics
- `11_Task3_OutputCollect.py` for unified reporting

### 3C. Fidelity-to-ground-truth check
- Correlate FM implied alignment with Task 2 ground-truth stability labels.

### 3D. Model meta-analysis
- Relate Task3 score slices to model attributes:
  - `size_class` (baseline / fm_standard / fm_large),
  - `perturbation_trained` flag.
- Report trend tests and group contrasts at the score-slice level.

---

# Part 5. Detailed Step-by-Step Execution Plan

## Phase A — Project initialization

1. Freeze metadata schema and naming dictionary (`cell_std`, `target_std`, dose/time harmonization rules).
2. Freeze gene universe and preprocessing defaults.
3. Freeze null-model and bootstrap settings for statistical consistency.

## Phase B — Task 1 implementation (Modality Concordance)

4. Build M1 matched index (same cell + same perturbation + best condition match).
5. Run group-wise modality metrics on full matched pairs.
6. Run strict protocol subset metrics and compare with full results.
7. Run bidirectional retrieval with raw and balanced candidate-space settings.
8. Run confounder decomposition and produce modality reliability labels.
9. Compute effect-size calibration (cross vs within replicates).
10. Compute continuous dose/time sensitivity curves and trend tests.
11. Run cell-stratified + leave-one-cell-out sensitivity (including HEPG2 leave-out).
12. Build benchmarkability zone map and coverage summary.
13. Publish Task 1 summary tables and QC diagnostics.

## Phase C — Task 2 execution (Mechanism Fidelity)

14. Run pairwise metrics on L1/L2/L3 with existing hardened pipeline.
15. Run retrieval scenarios A/B/C with existing multi-scenario pipeline.
16. Run attribution/protocol sensitivity module and generate context classes.
17. Add target-confidence/polypharmacology stratification outputs.
18. Integrate Task 1 reliability labels into Task 2 interpretation tables.

## Phase D — Task 3 execution (FM Stress Test)

19. Generate FM embeddings on matched evaluation set.
20. Compute cell-level deltas and Standard/Systema variants.
21. Run retrieval/pairwise FM scoreboard and subgroup analyses.
22. Run model meta-analysis (`size_class`, `perturbation_trained`).
23. Link FM fidelity scores to Task 2 ground-truth classes.

## Phase E — Release packaging

24. Produce unified result bundle (tables, metadata, manifests, logs).
25. Update version + change log in this `Plan.md`.
26. Tag release in GitHub using semantic versioning.

---

# Part 6. Coding and Collaboration Rules

## 6.1 Code quality rules (required)

- Code must be **readable**: clear function boundaries and meaningful names.
- Code must be **logical**: each script has explicit input, processing, output sections.
- Code must be **well-annotated**: short docstrings/comments explaining intent, assumptions, and non-obvious steps.
- Avoid silent behavior: key filters/drops must be logged to console and/or output files.
- Keep reproducibility: fixed random seeds for null/permutation/bootstrapping workflows.
- For analysis-heavy tasks, prefer **single, end-to-end data-analysis scripts** with explicit step blocks (load → transform → test → export), so users can read/modify workflow in one file.
- If helper modules are used, each script must still expose clear stage outputs (CSV/MD/manifest) without requiring deep code tracing.

## 6.2 Output and traceability rules

- Every major script writes:
  - parameter snapshot,
  - data manifest,
  - summary QC table.
- No analysis result is considered final without provenance fields (`source_db`, `view`, `track`, `task`, version).

## 6.3 Semantic versioning policy (GitHub)

Use `MAJOR.MINOR.PATCH`:

- **MAJOR**: breaking changes in schema, task logic, or interpretation contract.
- **MINOR**: backward-compatible new task modules/metrics/outputs.
- **PATCH**: bug fixes, refactors, documentation-only clarifications without schema breaks.

Release actions for every merged update:

1. bump version,
2. update `Plan.md` change log,
3. create Git tag/release note.

---

# Part 7. Project Change Log (maintained in Plan.md)

| Date | Version | Type | Summary | Breaking? | Owner |
| --- | --- | --- | --- | --- | --- |
| 2026-02-14 | v2.4.0 | MINOR | Added anti-drift reviewer-robustness layer: Task1 context-definition audit outputs (`task1_context_definition_comparison.csv`, `task1_metacontext_map.csv`), expanded protocol deconfounding specs in partial Spearman (`cell/target` + optional mismatch/match-score controls), new retrieval sensitivity script across balanced gallery sizes/seeds (`scripts/task1_retrieval_sensitivity.py`), and project-level alignment audit (`scripts/project_alignment_audit.py`) with PASS/FAIL checks against manuscript-facing claims. Updated workflow docs accordingly. | No | Team |
| 2026-02-14 | v2.3.0 | MINOR | Added submission-critical robustness upgrades: Task1 deconfounded protocol sensitivity via partial Spearman (`protocol_continuous_partial_spearman.csv`), strict-subset composition table (`strict_subset_composition.csv`), and full LINCS chemical internal-consistency analysis (`lincs_internal_consistency_*.csv`); marked Task2 target-tier layer as explicit sensitivity analysis with mapping QC output (`task2_target_tier_mapping_qc.csv`); clarified Task3 aggregate score definition with explicit export (`task3_mean_scaled_best_definition.csv`); added reproducibility-pack builder (`scripts/build_reproducibility_pack.py`) generating manifest/seed/environment/one-command artifacts under `outputs/reproducibility/analysis`. | No | Team |
| 2026-02-14 | v2.2.0 | MINOR | Added reviewer-critical robustness modules and outputs: Task1 add-on analyses (`scripts/task1_reviewer_addons.py`) for retrieval raw+balanced+random baseline reporting, cross-vs-within effect-size calibration, continuous dose/time sensitivity, leave-one-cell-out (including HEPG2) sensitivity, and benchmarkability-zone mapping; Task2 target-confidence/polypharmacology stratification (`scripts/task2_target_tier_analysis.py`); Task3 FM meta-analysis (`scripts/task3_fm_meta_analysis.py`). Updated plan contracts to lock claim hierarchy and analysis-unit definitions. | No | Team |
| 2026-02-14 | v2.1.1 | PATCH | Refined Task1 documentation contract with mandatory full-vs-strict protocol comparison and raw-vs-balanced retrieval reporting; appended current empirical snapshot from the 2026-02-14 rerun for co-author discussion. | No | Team |
| 2026-02-14 | v2.1.0 | MINOR | Added strict protocol subset analysis for Task1 (dose/time thresholds) and balanced retrieval evaluation to control candidate-space imbalance; added CSV export option for `unified_meta` and parquet-safe fallbacks; exposed new flags via Task1 pipeline and updated manuscript proposal to reflect mitigations. | No | Team |
| 2026-02-13 | v2.0.2 | PATCH | Updated Task2 analysis script (`scripts/task2_mechanism_analysis.py`) to support explicit domain scope and default to `same_only`, so mechanism analysis can be constrained to SameDomain data when required. Added Task3 engineering cleanup: `scripts/task3_results_audit.py` (objective consistency checks) and `scripts/task3_pipeline.py` (readable wrapper for Task3 scripts 9/10/11 + audit). Updated `docs/analysis_workflow.md` with Task3 audit/pipeline commands. | No | Team |
| 2026-02-13 | v2.0.1 | PATCH | Refactored Task0/Task1 script layer for readability and maintainability: added unified script bootstrap (`scripts/_script_bootstrap.py`), upgraded Task0/Task1 pipeline entrypoints (`scripts/task0_pipeline.py`, `scripts/task1_pipeline.py`) with stronger path fallback and manifest conventions, and aligned legacy Task0/Task1 wrappers to readable/compatible entry behavior. Updated `docs/analysis_workflow.md` with recommended Task0/Task1 pipeline commands. | No | Team |
| 2026-02-13 | v2.0.0 | MAJOR | Updated Task 2 interpretation contract: removed SameDomain vs CrossDomain as a primary comparison; re-centered Task 2 on mechanism-fidelity analysis with source/domain used only for diagnostics. Added readable analysis-style scripts: `scripts/task01_attrition_audit.py`, `scripts/task2_mechanism_analysis.py`, and `scripts/task2_results_audit.py`. Added Task0→Task1 attrition audit outputs and Task2 result consistency audit outputs. | Yes | Team |
| 2026-02-13 | v1.1.6 | PATCH | Added Task 1 retrieval parquet→CSV fallback loader and wrote `m1_candidates.csv` during group-wise pipeline for environment portability. | No | Team |
| 2026-02-13 | v1.1.5 | PATCH | Added Task 1 instance-level retrieval pipeline (`src/m2m_bench/task1/retrieval_instance.py`) with per-query and summary CSV outputs, null baseline summary, and example retrieval cases. Added `docs/task1_retrieval_spec.md`. | No | Team |
| 2026-02-12 | v1.1.4 | PATCH | Implemented Task 1 standalone pipeline (`src/m2m_bench/task1/groupwise_gap.py`) with independent data extraction, M1 matching, QC tables/plots, and CSV-first group-wise modality gap outputs. Added `docs/task1_data_spec.md`. | No | Team |
| 2026-02-12 | v1.1.3 | PATCH | Clarified canonical field semantics: `target_std` is gene-only, treatment type is `modality`; formalized task-wise matching keys for Task1/Task2/Task3. | No | Team |
| 2026-02-12 | v1.1.2 | PATCH | Added canonical Task 0 data structure documentation and explicit task-wise data usage guide (Task1/Task2/Task3) to improve project navigation. | No | Team |
| 2026-02-12 | v1.1.1 | PATCH | Added configurable Task 0 curation code under `src/m2m_bench/task0` with CLI entrypoint and local-output-friendly path defaults. | No | Team |
| 2026-02-12 | v1.1.0 | MINOR | Renamed project to M2M-Bench; updated title/subtitle; added project naming dictionary for consistent manuscript/code usage. | No | Team |
| 2026-02-12 | v1.0.0 | MAJOR | Renamed project to PerFi-Bench; unified two-question architecture (modality + mechanism); formalized 3-task structure and governance rules. | Yes | Team |

**Change log rule:** every meaningful update must append one row here.

---

# Part 8. Name Rationale and Alternatives

## Recommended current name

**M2M-Bench (Modality-to-Mechanism Benchmark)**

Why:

- directly encodes Task 1 (Modality) and Task 2 (Mechanism),
- compact and memorable acronym with clear conceptual symmetry,
- aligns with the project's ordered logic: concordance first, fidelity second.

## Alternatives (if journal branding preference changes)

- **Perturbation Fidelity Benchmark (PerFi-Bench)**
- **ChemGen Fidelity Benchmark**
- **Cross-Modality Perturbation Fidelity Benchmark**
- **Drug–Gene Transfer Fidelity Benchmark**

---

# Part 9. Project Naming Dictionary (canonical usage)

- `project_name_short`: **M2M-Bench**
- `project_name_long`: **Modality-to-Mechanism Benchmark**
- `project_subtitle`: **From Modality Concordance to Mechanism Fidelity: A Unified Perturbation Benchmark**
- `legacy_name`: **Chem2Gen-Bench**
- `module_name_task2_legacy`: **Chem2Gene (mechanism fidelity module/history)**

Usage rule:

- In repository docs and manuscript drafts, use: **M2M-Bench**.
- At first mention in formal writing, use: **M2M-Bench (formerly Chem2Gen-Bench)**.

---

# Part 10. Current Empirical Snapshot (2026-02-14 rerun + reviewer add-ons)

Task1 baseline:

- matched pairs (full): `670`
- strict protocol subset: `313` (`46.7%`)
- strict subset composition: `Genetic exact_cond only` (chemical strict pairs: `0`)
- global mean cosine (full): gene/path `0.0238 / 0.0599`
- global mean cosine (strict): gene/path `0.0060 / 0.0527`

Task1 retrieval dual report (`outputs/task1_reviewer_fixes/analysis/retrieval_dual_report.csv`):

- balanced setting: gallery `256`, positives `1`
- theoretical balanced random baselines: `Top1=1/256=0.003906`, `MRR=H256/256=0.023923`
- chemical mean MRR: raw `0.598` → balanced `0.104` (still `4.36×` random)
- chemical mean Top1: raw `0.489` → balanced `0.062`
- genetic mean MRR: raw `0.073` → balanced `0.032` (`1.36×` random)

Task1 effect-size calibration (`outputs/task1_reviewer_fixes/analysis/effect_calibration_summary.csv`):

- chemical gene track:
  - within-average `0.125`, cross-matched `0.039`, ratio `0.314`
- genetic gene track:
  - within-average `0.055`, cross-matched `0.006`, ratio `0.109`
- interpretation: cross-modality agreement is far below within-modality replicate consistency.

Task1 protocol continuous sensitivity (`outputs/task1_reviewer_fixes/analysis/protocol_continuous_spearman.csv`):

- `dose_logdiff` vs `cosine_gene`: rho `+0.169`, p `0.0013`
- `time_absdiff` vs `cosine_gene`: rho `-0.255`, p `<1e-5`
- `time_absdiff` vs `cosine_path`: rho `-0.250`, p `<1e-5`
- note: strict cutoff-only reporting is insufficient; continuous mismatch trends are non-trivial.

Task1 composition-bias sensitivity:

- leave-HEPG2-out (pairwise):
  - ALL mean cosine gene/path `0.0238 / 0.0599`
  - leave-HEPG2 mean cosine gene/path `0.0377 / 0.0630`
  - genetic-only drops from `0.0060 / 0.0527` to `-0.0087 / -0.0275`
- leave-HEPG2-out (retrieval, mean over 8 groups):
  - MRR `0.3358 → 0.3100`
  - balanced MRR `0.0684 → 0.0641`

Task1 benchmarkability map (`outputs/task1_reviewer_fixes/analysis/benchmarkability_zone_summary.csv`):

- `Z0_NoCrossSourceOverlap`: Chemical `16,761` contexts; Genetic `91,868` contexts
- `Z1_ChemicalNearestOnly`: `2` contexts (A549 `HSP90AA1`, `SIRT1`)
- `Z3_GeneticExactComparable`: `82` contexts

Task2 target-tier stratification (`outputs/task2_nodomain/target_tier/analysis/`):

- target tiers:
  - `Tier_A_HighConfidence`: `84` contexts (`6` targets)
  - `Tier_B_ModerateOrUnknown`: `2,946` contexts (`407` targets)
  - `Tier_C_AmbiguousOrPolypharm`: `31` contexts (`3` targets)
- significant enrichments (FDR<0.05):
  - `Robust_High` enriched in `Tier_A` (log2OR `+2.37`)
  - `Protocol_Sensitive` enriched in `Tier_C` (log2OR `+2.47`)

Task3 FM meta-analysis (`outputs/task3_meta/analysis/`):

- top overall mean scaled score: `scGPT (0.684)`, followed by `PCA200 (0.647)`
- no robust universal gain for perturbation-trained models across slices
  (best p-value `0.058`, non-significant after multiple comparisons).
- size-class trends are metric-dependent:
  - strong positive trend only for pairwise `Mean_NegEDist` (rho `0.699`, p `3.4e-08`)
  - retrieval metrics show weak/mixed size effects.
