# Project Proposal: M2M-Bench (Modality-to-Mechanism Benchmark) v1.1.0

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

---

# Part 3. Data Scope and Representation

## 3.1 Data assets

- **LINCS**: bulk perturbation transcriptomics
- **scPerturb**: single-cell perturbation transcriptomics
- **Matched K562 gold set** (Task 3 model stress-testing)

## 3.2 Harmonized feature space

- Shared curated gene universe (~2.5k genes)
- Standardized metadata fields (e.g., `cell_std`, `target_std`, `modality`, `source_db`, dose/time fields)

## 3.3 Subset design

- **L1:** same cell + same target (mechanism gold pairs)
- **L2:** different cell + same target (target generalization pool)
- **L3:** same cell + different target (background pool for Systema-style debiasing)
- **M1 (new):** same cell + same perturbation identity + matched conditions across LINCS/scPerturb (modality benchmark core set)

## 3.4 Signature views

1. **Standard view**: raw differential signatures
2. **Systema view**: background-subtracted signatures for stress/systematic component control

---

# Part 4. Three-Task Project Structure

## Task 1 — Modality Concordance Benchmark (Q1)

**Goal:** quantify LINCS↔scPerturb agreement for the same perturbation context.

### 1A. Group-wise modality gap
- Compare matched pairs in M1 using `cosine`, `cosine_z`, and `edist` (Standard/Systema).
- Build stratified null distributions and confidence intervals.

### 1B. Cross-source retrieval
- LINCS query → scPerturb gallery and reverse.
- Evaluate Top-k, MRR, and failure strata.

### 1C. Confounder decomposition
- Quantify effect of cell context, target class, dose/time mismatch, source, and perturbation class.
- Label contexts as **Robust / Conditional / Unstable** for downstream interpretation.

**Primary output:** modality reliability table per context (used as interpretation guardrail in Task 2).

---

## Task 2 — Chemical–Genetic Mechanism Fidelity (Q2; original Chem2Gene core)

**Goal:** measure and explain drug↔gene transferability after confounder-aware control.

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

---

# Part 5. Detailed Step-by-Step Execution Plan

## Phase A — Project initialization

1. Freeze metadata schema and naming dictionary (`cell_std`, `target_std`, dose/time harmonization rules).
2. Freeze gene universe and preprocessing defaults.
3. Freeze null-model and bootstrap settings for statistical consistency.

## Phase B — Task 1 implementation (Modality Concordance)

4. Build M1 matched index (same cell + same perturbation + best condition match).
5. Run group-wise modality metrics (Standard/Systema).
6. Run bidirectional cross-source retrieval.
7. Run confounder decomposition and produce modality reliability labels.
8. Publish Task 1 summary tables and QC diagnostics.

## Phase C — Task 2 execution (Mechanism Fidelity)

9. Run pairwise metrics on L1/L2/L3 with existing hardened pipeline.
10. Run retrieval scenarios A/B/C with existing multi-scenario pipeline.
11. Run attribution/protocol sensitivity module and generate context classes.
12. Integrate Task 1 reliability labels into Task 2 interpretation tables.

## Phase D — Task 3 execution (FM Stress Test)

13. Generate FM embeddings on matched evaluation set.
14. Compute cell-level deltas and Standard/Systema variants.
15. Run retrieval/pairwise FM scoreboard and subgroup analyses.
16. Link FM fidelity scores to Task 2 ground-truth classes.

## Phase E — Release packaging

17. Produce unified result bundle (tables, metadata, manifests, logs).
18. Update version + change log in this `Plan.md`.
19. Tag release in GitHub using semantic versioning.

---

# Part 6. Coding and Collaboration Rules

## 6.1 Code quality rules (required)

- Code must be **readable**: clear function boundaries and meaningful names.
- Code must be **logical**: each script has explicit input, processing, output sections.
- Code must be **well-annotated**: short docstrings/comments explaining intent, assumptions, and non-obvious steps.
- Avoid silent behavior: key filters/drops must be logged to console and/or output files.
- Keep reproducibility: fixed random seeds for null/permutation/bootstrapping workflows.

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
