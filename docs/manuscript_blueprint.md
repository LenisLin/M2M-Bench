# Manuscript Blueprint

## 1. Purpose

This document fixes the manuscript-level positioning, result roles, figure roles, claim boundaries, and evidence requirements for the current local M2M-Bench repository state. It preserves the frozen Task1/Task2 story, incorporates the accepted internal-review revisions, and does not add new results, rerun analyses, or redefine audited contracts.

Execution note: the next manuscript-support analysis sequencing, gating, and output planning now lives in `docs/manuscript_analysis_expansion_plan.md`. This blueprint remains the result-structure source of truth.

## 2. Paper type and scope

The target manuscript is a Brief Communication / Analysis-style short paper. The main paper should carry 2-3 figures and a compact claims structure. Supplementary material should carry denominator detail, sensitivity detail, audit context, and local-scope comparisons that are necessary for review but not central to the main scientific arc.

The paper reads from the current local repository state as source of truth. In the current manuscript contract:

- Task1 is authoritative through S2.
- Corrected multisource Task2 is authoritative through S6.
- S7 exists as project-level synthesis implementation evidence only; it is not the primary evidence source for Figure 2 or Figure 3 claims.
- If summary ambiguity exists, authoritative evidence resolves back to Task1 S1/S2 outputs and corrected Task2 S4/S5/S6 outputs rather than to S7.
- Legacy `data/task2_snapshot_v1/` and legacy scPerturb-only Task2 stages remain historical evidence only.

## 3. Primary positioning

M2M-Bench is an evaluation and boundary-definition framework for transcriptome-centric virtual cell modeling. It is not a full virtual cell model, but it is also not merely a representation paper. Its scientific role is to define and test the external-validity gap between within-modality predictive success and mechanism-level biological portability.

In manuscript terms, M2M-Bench evaluates what perturbation-response signal is measurable under modality-preserving settings and what signal remains portable under mechanism-shifting settings, using audited task contracts, denominator-aware reporting, and explicit scope restrictions rather than a new predictor.

The central conceptual contribution is the decomposition of perturbation-response concordance into two distinct evaluation layers:

- modality concordance
- mechanism concordance

Core framing sentence:

M2M-Bench addresses the external-validity gap between within-modality predictive success and mechanism-level biological portability by defining an audited evaluation framework in which modality-preserving concordance and mechanism-level concordance are measured separately and interpreted under explicit workflow, scope, and evidence contracts.

## 4. Results architecture

### Result 1

Benchmark design defines the problem, workflow, audit path, and evidence contract.

Result 1 is not a generic methods section. The benchmark design itself is part of the scientific contribution: it specifies the problem definition, task decomposition, workflow, scope boundaries, audit trail, denominator discipline, representation-scope rules, and frozen evidence path needed to distinguish measurable modality-preserving signal from mechanism-portable signal.

This result must establish:

- Task1 as the modality-concordance reference frame
- Task2 as the mechanism-concordance and external-validity layer
- the audited workflow from staged outputs to manuscript evidence
- the scope rules for common-scope `Gene`/`Pathway` comparisons and local FM use
- the evidence contract separating authoritative corrected outputs from historical or implementation-only material

### Result 2

Quantitative comparison of modality concordance and mechanism concordance defines the main benchmark contrast.

Result 2 must be written as a quantitative result, not as an abstract contrast. It requires explicit evidence for:

- Task1 internal comparisons
- Task1 cross comparisons in the supported genetic setting
- Task2 group-level comparisons
- Task2 direction-specific retrieval comparisons
- dataset- and cell-line-stratified reporting
- separate `C2G` and `G2C` reporting under the corrected retrieval contract
- representative aligned cases and representative misaligned cases tied to audited slices

The interpretive role of Result 2 is to show that modality-preserving signal can be measurable and stable while mechanism-level portability is more constrained. Where shared slices permit it, the manuscript should treat the Task1-to-Task2 contrast as a matched external-validity comparison rather than as a loose side-by-side summary.

### Result 3

Biological and experimental covariates explain heterogeneity in mechanism concordance.

Result 3 interprets Task2 heterogeneity through biological covariates, experimental covariates, and enrichment-style readouts rather than through a single pooled mechanism-concordance number. It should include:

- target complexity analysis
- dose/time analysis
- enrichment analysis
- local specificity analysis where the frozen scope supports it

Result 3 must keep the current local support boundaries explicit:

- dose/time can support matched or controlled claims only after a covariate-overlap audit shows that the compared slices share evaluable support
- if overlap is not demonstrated, dose/time remains descriptive context rather than controlled evidence
- local specificity analysis remains local where the authoritative summaries do not support a benchmark-wide cohort definition
- enrichment interpretation must use the matched/evaluable target pool as background, not the full genome or an unrestricted universe

### Result 4

Conditional representation effects are tested as a local benchmark comparison inside the supported FM scope.

Result 4 remains local to the supported scPerturb K562 Task2 subset, but it should be written as a concrete comparison rather than as a disclaimer section. The result asks:

- whether FM changes group-level conclusions
- whether FM changes retrieval conclusions
- whether FM changes target-level ranking or only local numerical summaries

Group-level evidence remains primary and retrieval remains supporting, but the manuscript should directly test whether FM alters the local interpretation of Task2 rather than treating restricted scope as the main message.

## 5. Main figures

### Figure 1

Result 1: benchmark design as scientific contribution.

Required content:

- problem definition: Task1 as modality concordance; Task2 as mechanism concordance and external-validity layer
- workflow diagram from audited stages to manuscript evidence
- scope diagram covering Task1 internal, Task1 cross, corrected multisource Task2, and restricted FM scope
- evidence-contract diagram separating authoritative corrected outputs, implementation-only synthesis, and historical materials
- denominator and audit principles, including representation-scope rules and explicit exclusions

Interpretive role:

- show that the benchmark design itself defines the scientific question and the admissible evidence
- make workflow, scope, audit, and evidence-contract logic part of the contribution rather than background methods

### Figure 2

Result 2 reference frame: modality concordance establishes the comparison baseline.

Required content:

- Task1 internal versus cross comparisons
- genetic cross as the supported cross-ecosystem result
- group-level metrics as the primary evidence layer
- retrieval as supporting evidence only
- denominator and support transparency for every reported slice

Interpretive role:

- show stable non-random signal in supported modality-preserving settings
- establish the calibrated reference frame for reading Task2
- avoid ceiling language for Task2 while making the Task1-to-Task2 contrast quantitatively usable

### Figure 3

Results 2-4: mechanism concordance, covariate interpretation, and local FM comparison.

Required content:

- Task2 group-level concordance in common scope
- Task2 direction-specific retrieval in `C2G` and `G2C`
- dataset and cell-line stratification
- representative aligned and misaligned cases
- target-complexity interpretation and other supported biological/external covariates
- dose/time comparisons only where a covariate-overlap audit supports matched reading
- enrichment-style interpretation with matched/evaluable target-pool backgrounds
- restricted-scope K562/scPerturb FM comparison

Interpretive role:

- make Task2 the main external-validity and biological-interpretation layer
- explain heterogeneity through covariates and case structure rather than through a pooled benchmark average
- test whether local FM scope changes conclusions or only local summaries

## 6. Core claims

### Claim 1

M2M-Bench contributes a benchmark design that defines the problem boundary, workflow, scope, audit path, and evidence contract for transcriptome-centric virtual cell evaluation.

### Claim 2

M2M-Bench exposes an external-validity gap: within-modality predictive success does not automatically transfer to mechanism-level biological portability, and this gap becomes visible only when Task1 and Task2 are evaluated separately and quantitatively.

### Claim 3

Task2 heterogeneity is biologically and experimentally structured. Its interpretation depends on dataset, cell line, target complexity, supported specificity context, and covariate-audited experimental context rather than on a single mechanism-concordance summary.

### Claim 4

Foundation-model effects are local and conditional: within the supported scPerturb K562 scope, the key question is whether FM changes group conclusions, retrieval conclusions, or only local target-level summaries.

## 7. Analysis and evidence requirements

The manuscript requires paired or statistical comparisons where valid. Summary-number reporting alone is not sufficient for the retained result path.

Minimum requirements:

- Task1 versus Task2 contrasts must use paired, matched, or otherwise explicit statistical comparisons wherever shared slices make that valid.
- Task1 internal, Task1 cross, Task2 group, and Task2 retrieval comparisons must all carry denominator and support fields rather than mean values alone.
- Task2 retrieval must keep `C2G` and `G2C` separate and include robustness checks against candidate-pool and gallery-size artifacts under the current corrected retrieval logic.
- Any `C2G` versus `G2C` asymmetry is a robustness-analysis question, not a reason to discard the corrected retrieval framework.
- Local FM comparisons require per-target robustness checks so the manuscript can tell whether FM changes target ranking or only local aggregate summaries.
- Robustness checks must be reported where mathematically relevant, including retrieval correction validation, denominator attrition, valid-mask effects, and other support-sensitive comparisons.
- Dose/time claims require a covariate-overlap audit before any matched or controlled comparison is stated.
- Enrichment analyses must use the matched/evaluable target pool as the background universe.

## 8. Claim guardrails

- Figure 1 must present benchmark design as part of the result structure, not as a generic methods appendix.
- Figure 2 must show stable non-random signal; weak-support slices stay out of the main claim path.
- Task1 is a calibrated reference frame, not a formal upper bound for Task2.
- Task1 cross remains genetic-only in the current contract; cross-chemical exclusion is a contract outcome, not a missing empirical result.
- Task2 core evidence remains stratified by `dataset` and `cell_line`; raw pooling across LINCS and scPerturb is not a core-result interpretation.
- Retrieval must follow the corrected contract-consistent interpretation, including explicit `gallery_definition_id`, `pos_definition_id`, and separate direction reporting.
- `C2G`/`G2C` asymmetry should be tested for robustness against candidate-pool structure, but it must not be assumed to be purely mathematical artifact.
- Group-level concordance remains the primary evidence layer; retrieval supports interpretation but does not replace group-level conclusions.
- Target-complexity claims must respect the frozen local aggregation scope and must not imply unsupported benchmark-wide symmetry from a local or direction-specific sensitivity.
- Local specificity analysis must stay local unless an authoritative cohort definition exists in the corrected summaries.
- Raw similarity magnitudes across representations must not be over-interpreted as a common effect-size scale.
- S7 provides project-level implementation context only and must not replace Task1 S1/S2 or corrected Task2 S4/S5/S6 as the evidence source for Figure 2 or Figure 3.
- FM results remain restricted to the currently supported scPerturb K562 scope and must not be generalized to LINCS or benchmark-wide behavior.

## 9. Discussion and limitation boundaries

The discussion should stay inside benchmark-evaluation claims. It should not drift into predictor novelty, mechanism causality, or platform-agnostic transfer claims that exceed the current audited contracts and evidence.

In this framing, M2M-Bench does not claim to solve transcriptome-centric virtual cell modeling. It defines a benchmark boundary for what current perturbation-response modeling can support under audited, scope-limited evaluation, especially when moving from modality-preserving success to mechanism-level portability.

Required limitation boundaries:

- The manuscript does not become a generative virtual-cell benchmark in this revision.
- The manuscript does not depend on zero-shot or OOD generation tasks for validity; those tasks are outside the current audited scope.
- The current corrected retrieval framework remains in force; asymmetry analyses are robustness checks, not a metric reset.
- Corrected Task2 is multisource, but its core results remain stratified by `dataset` and `cell_line`.
- FM evidence is limited to the scPerturb K562 subset.
- Local specificity, target-complexity, and dose/time interpretations must stay inside their supported audited scopes.
- Enrichment claims remain constrained to matched/evaluable target pools.
- Legacy scPerturb-only Task2 outputs are historical evidence only and must not be substituted for corrected Task2 claims.
