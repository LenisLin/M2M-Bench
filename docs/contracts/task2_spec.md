# Task2 Contract

**Mechanism Concordance under Explode-Matching + R-INSTANCE Retrieval**

## Status Note

The currently materialized `data/task2_snapshot_v1/k562` bundle and the current
`scripts/s3_build_task2_snapshot.py`, `scripts/s4_task2_group_concordance.py`,
`scripts/s5_task2_retrieval.py`, and `scripts/s6_task2_result_synthesis.py`
pipeline are preserved as **legacy/interim scPerturb-K562 evidence**. They are
valid as bounded historical evidence, but they are **not** the final corrected
Task2 implementation.

The corrected Task2 contract is **multisource** from core metric computation
onward. It must include both `LINCS` and `scPerturb`. The frozen corrected
successor snapshot root is `data/task2_snapshot_v2/`. Within that root, the
scPerturb contribution is materialized as `scperturb_k562/` by reusing the
audited legacy K562 bundle from `data/task2_snapshot_v1/k562/`, while the LINCS
contribution is derived from `data/task1_snapshot_v1/`.

---

## 0. Scope & Goal

Task2 evaluates **mechanism concordance** between **Chemical** and **Genetic**
perturbations within source dataset under matched biological context. The
authoritative Task2 analysis/cohort key is:

* `mech_key = (dataset, cell_line, target_token)`

Corrected Task2 v1 studies both first-class source datasets:

* `dataset = LINCS`
* `dataset = scPerturb`

Task2 produces:

1. Group-level concordance between Chemical and Genetic cohorts per `mech_key`
2. Instance-level retrieval in both directions with multi-positive chance correction
3. Dataset-stratified and cell-line-stratified benchmark summaries for reporting

Corrected Task2 v1 core metrics are computed **stratified by `dataset` and
`cell_line`**. LINCS and scPerturb are **not** raw-pooled in S4 or S5. Any
Task2-overall aggregate is a later synthesis-layer product and must be labeled
explicitly as a derived aggregate.

---

## 1. Identity Model: Row Identity vs Analysis/Cohort Identity

Task2 must distinguish **row identity** from **analysis/cohort identity**.

### 1.1 Snapshot row identity

Snapshot row identity is dataset-specific and representation-specific evidence
identity, for example:

* `row_id`
* `treated_cell_id`
* `query_uid`
* `global_idx_lincs`
* dataset-native row indices carried through `delta_meta` or per-query outputs

Row identity is used for audit evidence, delta alignment, FM validity, and
per-query retrieval provenance.

### 1.2 Analysis/cohort identity

The Task2 analysis/cohort key is:

* `mech_key = (dataset, cell_line, target_token)`

This key defines which Chemical and Genetic instances belong to the same
mechanism-level comparison cohort.

### 1.3 Explode-membership implication

Chemical instances may map to multiple `target_token` memberships through
`target_tokens`, but this does **not** create new snapshot rows. Membership is
relation-based:

* a Chemical row keeps one row identity
* that same row may contribute to multiple `mech_key` cohorts

This distinction must be preserved in snapshot, core metrics, and synthesis
outputs.

---

## 2. Representation-Scope Policy

Representation availability is a **scope policy**, not an attrition outcome.

### 2.1 LINCS

Supported representations:

* `Gene`
* `Pathway`

Unsupported on LINCS:

* `scgpt`
* `geneformer`
* `scbert`
* `scfoundation`
* `uce`
* `state`
* `tahoe-x1`

Absence of these rows on LINCS is **not-applicable scope**, not missing data and
not attrition.

### 2.2 scPerturb K562 subset

Supported representations:

* `Gene`
* `Pathway`
* `scgpt`
* `geneformer`
* `scbert`
* `scfoundation`
* `uce`
* `state`
* `tahoe-x1`

This frozen registry is the only FM scope approved for corrected Task2 v1.

---

## 3. Inputs (Snapshot) & Fail-fast Readiness

### 3.1 Corrected successor snapshot root

The corrected Task2 contract reads from a versioned successor root:

* `data/task2_snapshot_v2/`

The legacy/interim root that exists today is:

* `data/task2_snapshot_v1/`

The legacy root is preserved for historical evidence but is insufficient on its
own for corrected final Task2 because it does not yet include the required
multisource Task2 snapshot contract.

### 3.2 Required corrected inputs

**LINCS (derived from Task1-internal LINCS full set):**

* `data/task2_snapshot_v2/lincs/task2_lincs_pairs.csv`
* LINCS Gene and Pathway delta/meta assets reused from `data/task1_snapshot_v1/lincs/`
  according to the frozen Task1 contracts

**scPerturb K562 (materialized by reusing the audited legacy K562 bundle):**

* `data/task2_snapshot_v2/scperturb_k562/CRISPR_counts.pt`
* `data/task2_snapshot_v2/scperturb_k562/CRISPR_meta.csv`
* `data/task2_snapshot_v2/scperturb_k562/Drug_counts.pt`
* `data/task2_snapshot_v2/scperturb_k562/Drug_meta.csv`
* `data/task2_snapshot_v2/scperturb_k562/Common_Targets_K562.csv`
* `data/task2_snapshot_v2/scperturb_k562/shared_var_names.csv`

**Derived scPerturb K562 artifacts:**

* `data/task2_snapshot_v2/scperturb_k562/derived/pair_list.parquet`
* `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv`
* `data/task2_snapshot_v2/scperturb_k562/derived/gene_delta.npy`
* `data/task2_snapshot_v2/scperturb_k562/derived/pathway_delta.npy`

**scPerturb K562 FM artifacts (frozen registry only):**

* `data/task2_snapshot_v2/scperturb_k562/fm/<model>/fm_delta.npy`
* `data/task2_snapshot_v2/scperturb_k562/fm/<model>/fm_delta_meta.csv`
* `data/task2_snapshot_v2/scperturb_k562/fm/<model>/delta_operator_policy.json`

### 3.3 Fail-fast conditions

* Any required corrected successor path missing -> `TASK2_SNAPSHOT_NOT_READY`
* LINCS gating inventory missing -> `TASK2_LINCS_NOT_BUILT`
* scPerturb K562 derived artifacts missing -> `TASK2_K562_NOT_BUILT`
* Any required scPerturb K562 FM model missing outputs -> `TASK2_K562_FM_INCOMPLETE`

Explicitly documented non-blocking provenance notes do not fail readiness by
themselves unless they also violate one of the required contracts above.

---

## 4. Frozen Constants & Global Randomness

```text
GLOBAL_SEED = 619
EDIST_MAX_N = 256
K_CONTROLS_MAX = 50
MIN_CONTROLS_REQUIRED = 1
HIT_K_LIST = {1, 5, 10}
chance_identity_tol = 1e-12
rank_tie_rule = strict  # rank = 1 + #{score > true_score}
```

Randomness is initialized once per script using `GLOBAL_SEED`. No derived seeds.

---

## 5. Target Tokenization & Explode-Matching

### 5.1 Chemical tokenization

For each Chemical instance:

* `target_raw = clean_target_mapped`
* `target_tokens = unique(sorted(split(target_raw, ';', strip=True)))`
* `n_targets = len(target_tokens)`
* retain any available metadata such as `specificity_tier`, `time`, and `dose_value`

### 5.2 Genetic tokenization

For each Genetic instance:

* `target_token = clean_target_mapped`
* the genetic row contributes to exactly one `target_token`

### 5.3 Membership rule

A Chemical instance `i` contributes to `mech_key=(dataset, cell_line, t)` iff:

* `t in target_tokens(i)`

This same membership rule is reused in:

* S4 group-level cohort construction
* S5 retrieval positives definition
* any later stratified helper analysis

---

## 6. Task2 Data Construction

### 6.1 LINCS (no delta recomputation)

LINCS Task2 input is derived from Task1-internal LINCS full coverage, not from
Task1-cross subsets.

LINCS inclusion requires all of the following for the same source dataset:

* same `cell_line`
* same `target_token`
* at least one Chemical instance present
* at least one Genetic instance present

`task2_lincs_pairs.csv` is the gating list of eligible LINCS analysis/cohort
keys. LINCS retains `time` and `dose_value` as instance metadata where present,
but those fields do not enter the Task2 analysis/cohort key.

### 6.2 scPerturb K562 build

The scPerturb contribution to corrected Task2 v1 is the K562 subset. It remains
the only scope where FM representations are defined.

#### 6.2.1 Treated and control definitions

From `CRISPR_meta.csv` and `Drug_meta.csv`:

* treated rows: `benchmark_group == "Test"`
* control rows: `benchmark_group == "Control"`

#### 6.2.2 Control pool policy

* CRISPR: pool = all CRISPR Control rows
* Drug: pool = Drug Control rows with `time == treated.time`
* controls do not match target; controls are untargeted

#### 6.2.3 Pairing evidence

Pairing must remain auditable and deterministic:

* no replacement sampling
* if `pool_size == 0`, soft-fail into attrition
* else use all controls when `pool_size <= K_CONTROLS_MAX`
* else sample `K_CONTROLS_MAX` controls without replacement using `GLOBAL_SEED`
* persist explicit pairing evidence via `pair_list.parquet`

#### 6.2.4 Derived outputs

The scPerturb K562 build must persist:

* `delta_meta.csv`
* `gene_delta.npy`
* `pathway_delta.npy`
* FM delta bundles for the frozen FM registry

---

## 7. S4 Group-level Mechanism Concordance

For each `mech_key=(dataset, cell_line, target_token)` and each supported
representation in that dataset scope:

### 7.1 Cohort definition

* Chemical cohort:
  * `X = { delta_i | target_token in target_tokens(i), delta_valid_bool=True }`
* Genetic cohort:
  * `Y = { delta_j | genetic target_token == target_token, delta_valid_bool=True }`

If either side is empty, emit attrition with an explicit reason and do not
compute the metric row.

### 7.2 Metrics

* `cosine = cosine(centroid(X), centroid(Y))`
* `pcc = PearsonCorr(centroid(X), centroid(Y))`
* `edist = bias_corrected_energy_distance(X_sub, Y_sub)`

`edist` remains a valid within-row descriptive metric but is **not**
cross-representation comparable and must never be used to create a
cross-representation leaderboard rank.

### 7.3 Stratification policy

Corrected Task2 v1 group metrics are computed:

* stratified by `dataset`
* stratified by `cell_line`
* keyed by `target_token`

Raw pooling across LINCS and scPerturb is not allowed in S4.

---

## 8. S5 Retrieval

Retrieval remains two-direction and chance-corrected:

* `C2G` = Chemical -> Genetic
* `G2C` = Genetic -> Chemical
* Corrected Task2 v1 S5 is target-level retrieval for both `LINCS` and `scPerturb`.

Directions must remain separate in all summary and leaderboard outputs.

### 8.1 Common scoring

* `score(q, c) = -||q - c||^2`
* strict rank:
  * `rank_true = 1 + count(score_all > best_positive_score)`

### 8.2 C2G

* query = Chemical instance delta
* gallery = Genetic target centroids within the same `(dataset, cell_line)` scope
* positives = `target_tokens(query) intersect gallery_targets`

### 8.3 G2C

* query = Genetic instance delta
* gallery = Chemical target centroids within the same `(dataset, cell_line)` scope
* positives = the chemical gallery target matching the query target token
* `N_gallery = number of chemical target centroids within the same dataset/cell_line slice after validity filtering`
* this supersedes the prior chemical-context-count interpretation
* corrected Task2 S5 v1 does not depend on reconstructing chemical-context galleries from `Drug_meta.csv`, `perturbation_raw`, or any perturbagen-context key
* `time` and `dose_value` remain instance metadata / descriptive fields only, not gallery-defining keys

### 8.4 Chance correction

For each query, compute exact expectations from `(N_gallery, m_pos)` and emit:

* raw metrics
* expected chance metrics
* corrected metrics

### 8.5 Stratification policy

Corrected Task2 v1 retrieval metrics are computed:

* stratified by `dataset`
* stratified by `cell_line`
* stratified by `direction`

Raw pooling across LINCS and scPerturb is not allowed in S5.

---

## 9. Attrition & Scope Semantics

Task2 must distinguish:

* **scope policy**: not-applicable representation or dataset scope
* **attrition**: data that should have existed in-scope but failed eligibility or validity

Examples:

* no LINCS `scgpt` rows -> scope policy, not attrition
* empty Chemical or Genetic cohort for an eligible in-scope `mech_key` -> attrition
* FM `valid_mask=False` within scPerturb K562 -> attrition

---

## 10. Corrected Task2 v1 Output Checklist

### 10.1 Core metric outputs

1. `task2_pairs_coverage.csv`
2. `task2_group_concordance.csv`
3. `task2_group_attrition.csv`
4. `task2_retrieval_per_query.parquet`
5. `task2_retrieval_summary.csv`
6. `task2_retrieval_summary_long.csv`
7. `task2_retrieval_attrition.csv`
8. `task2_chance_identity_check.csv`

### 10.2 Task-level synthesis outputs

1. `task2_group_concordance_long.csv`
2. `task2_group_leaderboard.csv`
3. `task2_retrieval_leaderboard.csv`
4. `task2_benchmark_summary_long.csv`

### 10.3 Optional/helper audit inventory

* `task2_post_build_inventory.csv`

This file is helper/audit inventory only in corrected Task2 v1. It is not a
required core metric output unless separately re-approved.

---

## 11. Deferred or Approval-needed Outputs

The following are **not** part of corrected Task2 v1 core or task-level
synthesis contract unless separately approved:

* Task2 significance outputs
* benchmark-health outputs
* project-level synthesis outputs
* additional stratified helper tables such as lift/time-dose summaries beyond
  the locked corrected Task2 v1 tables above

If any of these return to the immediate Task2 stage map, that must be approved
explicitly before script implementation begins.
