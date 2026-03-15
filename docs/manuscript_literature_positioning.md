# Manuscript Literature Positioning

## 1. Purpose

This document positions M2M-Bench relative to current virtual-cell and perturbation-prediction literature without making evidence claims beyond the local repository state. It is a framing note for manuscript planning, not a substitute for local Task1 or corrected Task2 evidence.

## 2. Recommended top-level positioning

Preferred framing:

M2M-Bench is an evaluation and boundary-definition framework for transcriptome-centric virtual cell modeling.

Operational meaning for the manuscript:

- M2M-Bench evaluates perturbation-response concordance in transcriptomic readouts.
- It asks which signals remain measurable under modality-preserving settings and which remain portable under mechanism-shifting settings.
- It is not a full virtual cell model, not a complete cell simulator, and not a paper about building a new predictor.

## 3. What M2M-Bench Adds Conceptually

- It separates modality concordance from mechanism concordance instead of treating perturbation-response agreement as a single construct.
- It distinguishes within-modality comparability from cross-mechanism portability.
- It focuses on concordance structure and evaluation boundaries, not only prediction accuracy.
- It provides auditable benchmark contracts, denominator-aware reporting, and scope-limited interpretation rules.
- It adds an external-validity layer by asking not only whether responses can be predicted, but whether those signals remain portable across ecosystem and mechanism boundaries.

## 4. Continuity with Existing Literature

### Arc / Virtual Cell Framing

What that work contributes:

- Arc frames the virtual cell problem around predicting transcriptomic state shifts from a starting transcriptome plus perturbation.
- Arc's public Virtual Cell Challenge emphasizes benchmark datasets and rigorous standards for evaluating virtual cell models.
- Arc's first public challenge also states that current models are not yet consistently outperforming naive baselines across all metrics.

What question it leaves open:

- How evaluation should separate easier modality-preserving agreement from harder portability across ecosystem and mechanism boundaries.

How M2M-Bench complements rather than duplicates it:

- M2M-Bench does not propose a virtual cell model.
- It contributes an evaluation layer that asks which transcriptomic perturbation signals are stable enough to serve as a reference frame and which remain constrained under mechanism-shifting comparisons.

### Systema

What that work contributes:

- Systema argues that common evaluation metrics can be inflated by systematic variation.
- It also emphasizes that predicting unseen genetic perturbations is harder than standard metrics can suggest.

What question it leaves open:

- Whether apparently strong signal in transcriptomic perturbation evaluation remains stable once modality-preserving structure is separated from mechanism-level portability and cross-ecosystem constraints.

How M2M-Bench complements rather than duplicates it:

- M2M-Bench turns that concern into an explicit benchmark design question by separating Task1 reference-frame evidence from Task2 mechanism-concordance evidence.
- Its role is not to replace predictor benchmarking, but to add an auditable boundary on how far standard concordance interpretations can safely travel.

### Wei et al.

What that work contributes:

- Wei et al. benchmark single-cell perturbation prediction methods across cellular-context generalization and perturbation-generalization settings.
- Their results highlight that no single method works best across all datasets and that generalizability remains limited.

What question it leaves open:

- How much of the remaining performance heterogeneity reflects modality-specific comparability versus true mechanism portability across matched perturbation classes.

How M2M-Bench complements rather than duplicates it:

- M2M-Bench is model-agnostic and asks a different evaluation question.
- Rather than ranking prediction methods, it separates calibration-style modality concordance from the harder mechanism-concordance layer that constrains transfer claims.

### scDrugMap

What that work contributes:

- scDrugMap benchmarks foundation models for single-cell drug response prediction.
- It distinguishes pooled-data from cross-data evaluation and shows that evaluation design materially changes the conclusions.

What question it leaves open:

- Which conclusions remain stable when dataset scope, representation scope, and mechanism interpretation are all made explicit and audited.

How M2M-Bench complements rather than duplicates it:

- M2M-Bench makes dataset-stratified interpretation, scope-limited FM reading, and mechanism-versus-modality separation central to the benchmark design.
- It therefore serves as a complement to FM benchmarking, not as a benchmark-wide verdict on foundation models in general.

## 5. Safe Manuscript Language

### Intro-Style Positioning Sentences

- M2M-Bench is designed as an evaluation framework for transcriptome-centric virtual cell modeling rather than as a virtual cell model itself.
- The benchmark asks not only whether perturbation-response signal can be recovered in transcriptomic space, but also whether that signal remains portable across ecosystem and mechanism boundaries.
- In this framing, Task1 provides a modality-preserving reference frame, whereas Task2 tests the more constrained mechanism-concordance layer.

### Discussion-Style Significance Sentences

- M2M-Bench provides a boundary-definition layer for interpreting transcriptomic perturbation modeling claims under audited and scope-limited evaluation.
- The benchmark suggests that modality-preserving concordance and mechanism-level concordance should not be treated as interchangeable evidence for virtual cell progress.
- Read this way, Task2 functions as an external-validity check on perturbation-response modeling rather than as a leaderboard of representations.

### Limitation-Safe Sentences

- These results do not establish a full virtual cell model and should be interpreted only within the audited benchmark scopes reported here.
- Stronger concordance in Task1 should not be read as guaranteeing mechanism concordance in Task2.
- Any interpretation of foundation-model behavior remains local to the supported scPerturb K562 subset and should not be generalized benchmark-wide.

## 6. Overclaim Guardrails

Do not write this:

- `M2M-Bench is a virtual cell model`
- `M2M-Bench proves current virtual cell models fail because training data are flawed`
- `FM results generalize benchmark-wide`
- `modality concordance guarantees mechanism concordance`
- `M2M-Bench shows that current virtual cell methods do not work`
- `Task1 provides an upper bound for Task2`

Write this instead:

- `M2M-Bench is an evaluation and boundary-definition framework for transcriptome-centric virtual cell modeling.`
- `M2M-Bench defines a benchmark boundary for what current perturbation-response modeling can support under the audited local scopes; it does not assign a single cause for current limits.`
- `Foundation-model findings are interpretable only within the supported scPerturb K562 scope.`
- `Modality concordance and mechanism concordance should be evaluated separately.`
- `The benchmark identifies where current evaluation signal is stable enough to support interpretation and where portability remains more constrained.`
- `Task1 serves as a reference frame for Task2, not as a formal ceiling.`
