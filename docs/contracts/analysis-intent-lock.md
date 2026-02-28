# M2M v2 Analysis Intent Lock (Canonical)

This document is the canonical analysis contract for all subsequent coding and runs.
If any script behavior conflicts with this file, this file takes precedence.

---

## 1) Scientific Scope (Locked)

M2M v2 has two core tasks:

- **Task1 (Modality):** quantify modality consistency.
  - internal consistency first (within data source),
  - cross-modality consistency second (matched subset only).
- **Task2 (Mechanism):** quantify Chemical↔Genetic mechanism consistency within source.

No standalone “Task3” main axis in v2. FM analyses are integrated into Task1/Task2 as representation-space complements.

---

## 2) Data Scope and Cohort Rules (Locked)

### 2.1 Task1

- **Task1 internal (primary):**
  - Use **full eligible internal cohorts** for LINCS and scPerturb (not m1 matched subset).
  - Analyze LINCS Chemical, LINCS Genetic, scPerturb Chemical, scPerturb Genetic.
  - Representation spaces:
    - delta gene,
    - delta pathway,
    - FM-delta (scPerturb only).

- **Task1 cross-modality (secondary but required):**
  - Use **matched LINCS↔scPerturb subset only** (m1/matched-pair scope).
  - Purpose: quantify cross-modality comparability under explicit overlap constraints.
  - Chemical cross-modality is retained as diagnostic evidence even if sparse.
  - Direction policy:
    - Primary endpoint direction: `LINCS -> scPerturb`.
    - Sensitivity direction: `scPerturb -> LINCS` (always computed and reported).

### 2.2 Task2

- **Task2 delta-space:**
  - Use full eligible LINCS and scPerturb mechanism contexts (within-source Chemical↔Genetic).
  - Multi-target chemical handling uses explicit target-expansion matching (`explode_chemical` policy family).
  - Reliability gate:
    - context is valid only if both modalities satisfy `n_rows >= min_rows_per_modality` under selected context scope.
    - invalid contexts are excluded from primary summaries but retained in diagnostics tables.

- **Task2 FM-space:**
  - FM branch runs on **K562 contract subset** only (as defined by benchmark construction constraints).
  - Metrics and view definitions must stay isomorphic to Task2 delta-space.

### 2.3 scPerturb FM dataset construction

- scPerturb FM branch uses control/treated cell sets built from raw counts + processed metadata.
- FM embedding delta is defined as:
  - `delta_z = z_treated - mean(z_control_pool)`.
- Control-pool construction and treated/control linkage must be explicit, auditable, and deterministic under seed.

### 2.4 Context validity and denominator policy

- All context-level tables must carry:
  - `n_rows_chemical`,
  - `n_rows_genetic`,
  - `context_valid` flag.
- Primary context-level summaries include only `context_valid == True`.
- Supplemental diagnostics must include full pre-gate candidate context universe.

---

## 3) Analysis Units (Locked)

- `row`: one perturbation profile.
- `context`:
  - Task1 internal: source/modality/cell/target granularity.
  - Task1 cross: matched LINCS-scPerturb comparison key.
  - Task2: source-specific mechanism context key for Chemical↔Genetic matching.
- `query`: one retrieval query instance.
- `gallery`: retrieval candidate set for one query.

All summary tables must expose denominators (`n_rows`, `n_contexts`, `n_queries`, `n_gallery_mean` or equivalent).

---

## 4) Shared Metric Contract (Task1 & Task2)

## 4.1 Pairwise branch (retain, do not remove)

- Centroid-level cosine similarity.
- Spearman correlation.
- HVG-Top50 Spearman:
  - within context, rank genes by `|mean_chem - mean_gen|`,
  - select top-50 genes,
  - compute Spearman on top-50.
- Legacy-compatible E-distance; report harmonized `NegEDist = -EDist` where needed.

Task2 primary interpretation is centroid-level; instance-level pairwise remains supplementary.

## 4.2 Retrieval branch (retain, do not remove)

- Retrieval mode for mainline analysis:
  - **instance → centroid**.
- Report:
  - `MRR`,
  - `Top1`, `Top5`, `Top10`,
  - per-query outputs for audit.
- Balanced retrieval:
  - fixed gallery size,
  - fixed positives per query,
  - repeated subsampling with seed control.
- Random baseline under balanced setting:
  - `Topk_random = min(1, k / G)`,
  - `MRR_random = H_G / G`.

`two_stage` retrieval is removed from the v2 mainline implementation.

## 4.3 Space-uplift test family (locked)

- Space-uplift comparisons must be evaluated in a predefined family:
  - `gene_delta -> pathway_delta`,
  - `gene_delta -> fm_delta` (when fm branch exists for the same cohort),
  - `pathway_delta -> fm_delta` (optional, when both are available).
- Statistical protocol:
  - paired Wilcoxon signed-rank (primary),
  - BH-FDR correction within each task × endpoint family.
- No post-hoc metric family expansion without updating this contract first.

---

## 5) View Policy: Standard vs Systema

- `Standard` is primary view.
- `Systema` is secondary robustness view.
- Both views are computed for the same cohorts/metrics.
- View comparison uses paired statistical tests + FDR correction.

---

## 6) Covariate and Statistics Package (Task1 & Task2 Unified)

The covariate module must be symmetric across Task1 and Task2:

- Univariate association:
  - Spearman (`metric` vs covariate).
- Multivariable deconfounding:
  - joint model / permutation-compatible effect attribution,
  - report independent contribution (coefficient/partial effect/ΔR² style output).
- Enrichment analysis:
  - high- vs low-consistency group enrichment on cell/target/target-family/multi-target properties,
  - effect size + p-value + BH/FDR.
- Stability sensitivity:
  - leave-one-cell-out / leave-one-target-out summaries where applicable.

Minimum interpretation coverage:
- coverage bias (cell/target/sample-size),
- protocol dependence (dose/time where defined),
- representation effect (space uplift and view uplift),
- context enrichment (high/low consistency strata).

---

## 7) Output Contract (Minimum Required)

Each module run must emit:

- `run_manifest_*.json` (full config + seed),
- `*_key_summary.csv` (core denominators + run-mode tags),
- pairwise detail + summary tables,
- retrieval per-query + summary tables,
- balanced/random baseline tables (if retrieval enabled),
- covariate statistics tables (when covariate module is invoked).

---

## 8) Mandatory Intent-Gate Checklist (Run Before/After Any Major Change)

Before editing or running:

1. Is this Task1 internal run using full internal cohort, not matched subset?
2. Is this Task1 cross run explicitly limited to matched subset?
3. Is Task2 delta using full eligible mechanism cohorts?
4. Is Task2 FM restricted to K562 contract subset?
5. Are both pairwise and retrieval branches enabled?
6. Is retrieval mode for mainline set to `instance_to_centroid`?
7. Are Standard and Systema both produced?
8. Are seeds/manifests captured?

After run:

1. Do key summaries report expected cohort-mode tags and denominators?
2. Do pairwise and retrieval outputs both exist?
3. Do balanced/random baseline outputs exist for retrieval?
4. Are covariate/enrichment outputs present for the covariate phase?
5. Any mismatch between declared analysis scope and actual input cohort? (must be zero tolerance)

---

## 9) Engineering Runtime Policy (Hard Requirements)

These rules are mandatory for all new/modified scripts in M2M v2.

### 9.1 Terminal observability (must-have)

- Every long-running script must print:
  - run start banner (task/mode/input/output/seed/jobs),
  - progress heartbeat (at least per major split unit, e.g. dataset/context block),
  - completion summary (elapsed time + key denominators).
- Silent loops are disallowed for large runs.
- Progress output must be deterministic and human-readable (no ambiguous status text).

### 9.2 Memory safety (must-have)

- Large arrays/tables must be processed in split units (dataset/chunk/shard), not held globally when avoidable.
- Intermediate objects must be released after each split unit (`del` + `gc.collect()` for heavy objects).
- For very large jobs, write split-level intermediates first, then merge; avoid monolithic in-memory accumulation.
- Any change that increases peak memory must document expected RAM range in run docs.

### 9.3 Large-job split + parallel policy (must-have)

- Any heavy analysis module must expose explicit split/parallel knobs (e.g., `--dataset-jobs`, `--chunk-size`, or shard args).
- Split key must follow analysis semantics (Task1/Task2 typically split by `dataset_std` first).
- Parallel execution must preserve deterministic merge behavior and output schema.
- If environment blocks process-based parallelism, thread-based parallel fallback is required and must be logged.

### 9.4 Readability and annotation policy (must-have)

- Each script must include module-level docstring with:
  - analysis intent,
  - primary inputs/outputs,
  - denominator semantics.
- Non-trivial functions must include concise comments on algorithm intent (not line-by-line restatement).
- Public data structures/configs must be type-annotated.
- Hidden magic constants are disallowed; expose them as named config parameters.

### 9.5 Execution checklist (enforcement)

Before accepting any pipeline/script change:

1. Does runtime log progress clearly enough for remote monitoring?
2. Is peak memory controlled via split/chunk processing?
3. Does module support parallel/split execution for large cohorts?
4. Are comments/docstrings sufficient for handoff-level readability?
5. Are output schemas unchanged or explicitly versioned?

If any answer is “No”, the change is incomplete.

---

## 10) Change-Control Rule

Any future code change must reference this file and explicitly state:

- which section(s) it implements or modifies,
- whether the change affects Task1 internal, Task1 cross, Task2 delta, or Task2 FM scope,
- whether denominator semantics changed.
