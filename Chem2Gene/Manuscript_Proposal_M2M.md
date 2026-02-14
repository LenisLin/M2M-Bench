# M2M-Bench Manuscript Proposal (Draft for Co-author Review)

**Version:** v0.5 (2026-02-14)  
**Project:** M2M-Bench (*From Modality Concordance to Mechanism Fidelity*)  
**Goal:** finalize a reviewer-robust Results/Methods figure map with locked units, endpoints, and reproducibility outputs.

---

## 1) Claim hierarchy, units, and endpoints (locked)

### 1.1 Claim hierarchy
- **Primary contribution:** a reusable **audit + evaluation standard** for perturbation benchmarking.
- **Primary empirical finding:** cross-modality concordance is low and context-dependent under current overlap.
- **Secondary findings:** mechanism fidelity is context-dependent; FM behavior is direction/view-dependent.

### 1.2 Analysis units and denominators
- **Task1 `row`:** one perturbation profile in unified metadata.
- **Task1 `context`:** `(cell_std, modality, target_std)`.
- **Task1 `matched_pair`:** one LINCS row matched to one scPerturb row within a context.
- **Task1 `query`:** one retrieval query row against opposite-source gallery.
- **Task2 `instance`:** one row in `Step1_L1_Instance_Tidy_NoDomain.csv` keyed by `(Track, Scenario, View, Direction, Cell, Target, Source_Chem, Source_Gene, Dose, Time, CondID)`.
- **Task2 `context` (Step1):** `(Track, View, Cell, Target, Source_Chem, Source_Gene)` aggregated row.
- **Task2 `labeled_context` (Step2):** `(Cell, Target)` with class label.
- **Task2 `class`:** `{Robust_High, Intermediate, Robust_Low, Protocol_Sensitive}`.
- **Current denominators:** Task1 overlap contexts `84`, matched pairs `670`; Task2 Step1 contexts `9852`; Task2 labeled contexts `3061`.

### 1.3 Primary endpoints (submission lock)
- **Task1 primary endpoints:**
  1) pairwise `mean cosine_gene` (with cross-vs-within calibration),  
  2) retrieval `MRR_balanced` (gallery `256`, true `1`) with random baseline.
- **Task2 primary endpoints:**
  1) context-level `Mean_Success`/class composition on labeled contexts,  
  2) protocol-sensitive fraction and protocol-association counts.
- **Task3 primary endpoints:**
  1) retrieval `Mean_MRR` by direction/view,  
  2) pairwise `Mean_CentroidCosine` by view.

---

## 2) Proposed Results structure (submission-facing)

### R0. Motivation
- Separate **modality gap** from **mechanism gap** before biological claims.

### R1. Study design + benchmarkability (Fig1)
- Data coverage, attrition, overlap scarcity, and comparability zones.

### R2. Task1 modality concordance (Fig2; 4 panels only)
- Main story: low cross-modality agreement, retrieval inflation under imbalance, protocol sensitivity, and composition bias.

### R3. Task2 mechanism fidelity (Fig3)
- Same-domain mechanism fidelity with context classing and protocol-sensitive strata.

### R4. Task3 FM stress test (Fig4)
- Direction/view-dependent model behavior plus metadata-based interpretation.

### R5. Discussion boundaries
- Explicit scope limits and non-claims.

### R6. Reproducibility package
- Seeds/manifests/environment/one-command paths as first-class outputs.

---

## 3) Main figures (Fig1–Fig4)

## Fig1. Study design + benchmarkability + unit box

### Panel plan
- **Fig1A:** Task flow (Task1 → Task2 → Task3).
- **Fig1A2:** analysis-unit/denominator box (Task1 vs Task2 objects explicitly separated).
- **Fig1B:** attrition funnel (`2,175,012 rows / 108,713 contexts` → `17,078 candidate rows / 84 overlap contexts` → `670 matched pairs`).
- **Fig1C:** overlap composition (Chemical `2` contexts; Genetic `82`, HEPG2-dominant).
- **Fig1D:** drop reasons (`cell_and_target_absent_in_other_source`, `target_not_shared_in_other_source`).
- **Fig1E:** protocol mismatch diagnostics (chemical nearest matching vs genetic exact matching).
- **Fig1F:** benchmarkability map (`Z0`, `Z1`, `Z3` zones).

### Required files
- `outputs/task1_audit/analysis/stage_summary.csv`
- `outputs/task1_audit/analysis/context_drop_reasons_summary.csv`
- `outputs/task1_audit/analysis/protocol_gap_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/benchmarkability_zone_summary.csv`

---

## Fig2. Task1 modality concordance (main figure reduced to 4 panels)

### Main panels only
- **Fig2A (Gap + calibration):**
  - Pairwise cross-modality cosine (gene/path) + cross-vs-within calibration.
  - Cross is a minority fraction of within-modality consistency.
- **Fig2B (Retrieval fairness):**
  - raw vs balanced vs random baseline (`Top1_random=1/256`, `MRR_random≈0.0239`).
  - Chemical mean MRR `0.598 -> 0.104` after balancing.
- **Fig2C (Protocol sensitivity, deconfounded):**
  - continuous mismatch curves + partial Spearman controlling for `cell_std` and `target_std`.
  - Example (gene): `time_absdiff` partial rho `-0.312` (p `1.7e-9`), `dose_logdiff` partial rho `+0.151` (p `0.0043`).
- **Fig2D (Composition bias):**
  - ALL vs leave-HEPG2 sensitivity:
    - ALL gene/path `0.0238/0.0599`
    - leave-HEPG2 gene/path `0.0377/0.0630`
    - genetic-only flips to negative (`-0.0087/-0.0275`).

### Supplement-only Task1 panels
- strict subset comparison (**explicitly mostly/entirely genetic**),
- context true-ratio examples,
- confounder decomposition,
- set-level centroid analyses,
- LINCS internal consistency tables.

### Required files
- `outputs/task1/analysis/modality_gap_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/retrieval_dual_report.csv`
- `outputs/task1_reviewer_fixes/analysis/effect_calibration_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/protocol_continuous_spearman.csv`
- `outputs/task1_reviewer_fixes/analysis/protocol_continuous_partial_spearman.csv`
- `outputs/task1_reviewer_fixes/analysis/protocol_continuous_curves.csv`
- `outputs/task1_reviewer_fixes/analysis/leave_one_cell_out_pairwise.csv`
- `outputs/task1_reviewer_fixes/analysis/strict_subset_composition.csv`
- `outputs/task1_reviewer_fixes/analysis/lincs_internal_consistency_summary.csv`

---

## Fig3. Task2 mechanism fidelity (same-domain main analysis)

### Panel plan
- **Fig3A:** global tests (Systema vs Standard, gene/path contrasts).
- **Fig3B:** class composition on labeled contexts (`N=3061`).
- **Fig3C:** representative robust/fragile contexts.
- **Fig3D:** protocol-sensitive context enrichment (dose/time associations).
- **Fig3E:** target-tier/polyrisk stratification (**sensitivity analysis layer**, not primary endpoint).

### Current tier snapshot
- `Tier_A_HighConfidence`: `84` contexts (`6` targets), Robust_High fraction `52.4%`.
- `Tier_C_AmbiguousOrPolypharm`: `31` contexts (`3` targets), Protocol_Sensitive fraction `48.4%`.
- Mapping currently heuristic; external curated map remains preferred for submission final.

### Required files
- `outputs/task2_nodomain/analysis/Step1_L1_Instance_Tidy_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step1_L1_Context_Aggregated_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step2_Context_Labels_NoDomain.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_summary.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_enrichment.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_mapping_qc.csv`

---

## Fig4. Task3 foundation model stress test

### Panel plan
- **Fig4A:** retrieval scoreboard by direction/view.
- **Fig4B:** pairwise centroid cosine scoreboard.
- **Fig4C:** Systema − Standard lollipop (CRISPR→Drug).
- **Fig4D:** audit consistency checks.
- **Fig4E:** metadata-aware meta-analysis (`size_class`, `perturbation_trained`).

### Definition lock
- `mean_scaled_best` = mean(`Scaled_Best`) across metric slices for one model; equal weight per slice.
- Use as cross-slice trend score, not as a single-task endpoint.

### Required files
- `outputs/task3_audit/analysis/task3_consistency_checks.csv`
- `outputs/task3_meta/analysis/task3_model_scoreboard_meta.csv`
- `outputs/task3_meta/analysis/task3_mean_scaled_best_definition.csv`
- `outputs/task3_meta/analysis/task3_meta_perturbation_training_tests.csv`
- `outputs/task3_meta/analysis/task3_meta_size_trend_tests.csv`

---

## 4) Supplementary (focused additions)

### S1–S3
- Attrition accounting, overlap atlas, protocol mismatch diagnostics.

### S4 (Task1 robustness)
- strict subset caveat table (`strict_subset_composition.csv`),
- continuous + partial Spearman protocol analysis,
- leave-one-cell-out results,
- LINCS internal consistency:
  - same `(cell,target)` weighted cosine gene/path: `0.124 / 0.157`,
  - same `(cell,target,dose,time)` weighted cosine gene/path: `0.354 / 0.436`.

### S5–S8
- retrieval deep-dive, set-level analyses, Task2 deep diagnostics/tier mapping, Task3 full audit/meta tables.

---

## 5) Current empirical snapshot (2026-02-14 rerun)

- **Task1:** matched pairs `670`; strict `313` and strict composition is `100% Genetic`.
- **Task1 pairwise cosine:** full gene/path `0.0238 / 0.0599`.
- **Task1 retrieval (chemical):** MRR `0.598 -> 0.104` (raw→balanced), Top1 `0.489 -> 0.062`.
- **Task1 protocol deconfounding (chemical):**
  - partial rho(`time_absdiff`, `cosine_gene`) `-0.312` (p `1.7e-9`),
  - partial rho(`dose_logdiff`, `cosine_gene`) `+0.151` (p `0.0043`).
- **Task2 class counts (N=3061):** `Robust_High=567`, `Intermediate=1020`, `Robust_Low=1020`, `Protocol_Sensitive=454`.
- **Task3 meta:** top `mean_scaled_best` model is `scGPT (0.684)`, then `PCA200 (0.647)`; no robust universal perturbation-trained advantage (best p≈`0.058`).

---

## 6) Claims we do NOT make (must state explicitly)

- We do **not** claim broad chemical cross-modality generalization beyond current sparse overlap.
- We do **not** claim global drug↔gene mechanism equivalence across contexts.
- We do **not** claim a universal best FM across all directions/views/objectives.
- We do **not** treat target-tier heuristic mapping as definitive pharmacology annotation.

---

## 7) Reproducibility package (implemented)

Outputs from `scripts/build_reproducibility_pack.py`:
- `outputs/reproducibility/analysis/repro_manifest_inventory.csv`
- `outputs/reproducibility/analysis/repro_seed_registry.csv`
- `outputs/reproducibility/analysis/repro_environment.json`
- `outputs/reproducibility/analysis/repro_one_command_paths.md`

Current status:
- manifest coverage: `8/8` found,
- seed entries extracted: `4`.

---

## 8) Immediate next actions

1. Freeze Fig2 main 4-panel design and move remaining Task1 diagnostics to supplement.
2. Replace heuristic Task2 target-tier map with curated external annotations if available.
3. Draft Results text with locked units/denominators and explicit non-claims.
4. Prepare co-author review with reproducibility pack attached.
