# Task2 Proposal

**Mechanism Concordance under Explode-Matching + R-INSTANCE Retrieval**

## 0. Scope & Goal

Task2 evaluates **mechanism concordance** between **Chemical** and **Genetic** perturbations under matched contexts. The ecological unit is the **mechanism key**:

* `mech_key = (dataset, cell_line, target_token)`

Task2 produces:

1. **Group-level concordance** between Chemical vs Genetic distributions/centroids per `mech_key`
2. **Instance-level retrieval** in both directions with **multi-positive** chance correction (`m_pos > 1`)
3. **Stratified interpretability** by `specificity_tier` and `n_targets` (multi-target degree)
4. **Dose/time reporting** (distribution tables + C2G per-query metadata; G2C context gallery explicitly uses dose/time)

---

## 1. Inputs (Snapshot) & Fail-fast Readiness

### 1.1 Required inputs

**LINCS (Task2-LINCS derived from Task1-internal LINCS full set):**

* `data/task2_snapshot_v1/lincs/task2_lincs_pairs.csv`

**K562 (scPerturb_K562, from Evaluation_Set_K562; strict row alignment already manually verified):**

* `data/task2_snapshot_v1/k562/CRISPR_counts.pt`
* `data/task2_snapshot_v1/k562/CRISPR_meta.csv` (has `cell_barcode`)
* `data/task2_snapshot_v1/k562/Drug_counts.pt`
* `data/task2_snapshot_v1/k562/Drug_meta.csv` (rowname is cell_id)
* `data/task2_snapshot_v1/k562/Common_Targets_K562.csv`
* `data/task2_snapshot_v1/k562/shared_var_names.csv`

**Derived artifacts (must be created by Task2 build step):**

* `data/task2_snapshot_v1/k562/derived/pair_list.parquet`
* `data/task2_snapshot_v1/k562/derived/delta_meta.csv`
* `data/task2_snapshot_v1/k562/derived/gene_delta.npy`
* `data/task2_snapshot_v1/k562/derived/pathway_delta.npy`

**FM models (mandatory; list frozen from task1_snapshot_v1/fm_delta):**
`geneformer, pca50, pca100, pca200, scbert, scfoundation, scgpt, state, tahoex1_3b, uce`

Per model must exist:

* `data/task2_snapshot_v1/k562/fm/<model>/fm_delta.npy`
* `data/task2_snapshot_v1/k562/fm/<model>/fm_delta_meta.csv` (must include `row_id`, `treated_cell_id`, `valid_mask`, `n_controls_used`)
* `data/task2_snapshot_v1/k562/fm/<model>/delta_operator_policy.json`

### 1.2 Fail-fast conditions

* Any required path missing → `TASK2_SNAPSHOT_NOT_READY`
* K562 derived artifacts missing → `TASK2_K562_NOT_BUILT`
* Any required FM model missing outputs → `TASK2_K562_FM_INCOMPLETE`

---

## 2. Frozen Constants & Global Randomness

```text
GLOBAL_SEED = 619
EDIST_MAX_N = 256
K_CONTROLS_MAX = 50
MIN_CONTROLS_REQUIRED = 1
HIT_K_LIST = {1, 5, 10}
chance_identity_tol = 1e-12
rank_tie_rule = strict  # rank = 1 + #{score > true_score}
```

Randomness is initialized **once** per script using `GLOBAL_SEED`. No derived seeds.

---

## 3. Target Tokenization & Explode-Matching (No Instance Duplication)

### 3.1 Tokenization

For each Chemical instance:

* `target_raw = clean_target_mapped`
* `target_tokens = unique(sorted(split(target_raw, ';', strip=True)))`
* `n_targets = len(target_tokens)`
* retain `specificity_tier` (missing → `"NA"`)

For each Genetic instance:

* `target_token = clean_target_mapped` (must be single token; Control rows excluded)

### 3.2 Explode-Matching membership rule

Chemical instances are **not duplicated** as multiple rows. Instead, membership is relation-based:

A Chemical instance `i` contributes to `mech_key=(dataset, cell_line, t)` iff:

* `t ∈ target_tokens(i)`

This rule is used consistently in:

* Group-level cohort construction (`X` set)
* Retrieval positives definition (`m_pos` computation)
* Tier/n_targets stratification joins

---

## 4. Task2 Data Construction

## 4.1 LINCS (no delta recomputation)

* Uses precomputed LINCS deltas (Gene; Pathway via W projection) and meta.
* Chemical side instances retain `time` and `dose_value` as metadata (do not define groups by them in group-level concordance; see retrieval gallery definition).

`task2_lincs_pairs.csv` is the gating list of `(cell_line, target_token)` that have both Chemical and Genetic availability in LINCS.

## 4.2 K562 (scPerturb_K562) build (pair-consistent Gene/Pathway/FM)

### 4.2.1 Treated & control definitions

From `CRISPR_meta.csv` / `Drug_meta.csv`:

* Treated rows: `benchmark_group == "Test"`
* Control rows: `benchmark_group == "Control"`

### 4.2.2 Control pool (frozen)

* CRISPR: pool = all CRISPR Control
* Drug: pool = Drug Control with `time == treated.time`
* Control does **not** match target (controls are untargeted)

### 4.2.3 Pairing policy (v1.3 key fix)

No replacement sampling.

For each treated:

* If `pool_size == 0`: **soft-fail**, drop this treated, log attrition
* Else:

  * If `pool_size <= K_CONTROLS_MAX`: take **all** controls (sorted by control_cell_id) for determinism
  * If `pool_size > K_CONTROLS_MAX`: sample `K_CONTROLS_MAX` controls **without replacement** using RNG(seed=619)
  * Define `n_controls_used = min(pool_size, K_CONTROLS_MAX)`

Output `pair_list.parquet` schema (each row is a treated-control pair):

* `treated_cell_id`
* `control_cell_id`
* `control_rank` (1..n_controls_used)
* `n_controls_used` (repeated per treated for audit)
* `dataset_side` ∈ {CRISPR, DRUG}
* `perturbation_class` ∈ {Genetic, Chemical}
* `cell_line="K562"`
* `target_raw`, `target_tokens`
* `time`, `dose_value` (Drug; CRISPR may be NA)
* `specificity_tier` (Drug; CRISPR may be NA)
* `seed=619`

### 4.2.4 Gene/Pathway delta (pair-consistent)

For each treated:

* `ctrl_mean = mean(counts[controls(t)])` using `n_controls_used(t)`
* `gene_delta = counts[t] - ctrl_mean`
* `pathway_delta = gene_delta @ W`

Output:

* `k562/derived/gene_delta.npy`, `k562/derived/pathway_delta.npy`
* `k562/derived/delta_meta.csv` includes:

  * `row_id`, `treated_cell_id`, `perturbation_class`, `cell_line`, `target_raw`, `time`, `dose_value`, `specificity_tier`, `n_controls_used`

### 4.2.5 FM embedding + FM delta (pair-consistent; dynamic control_counts)

For each FM model:

* Compute embeddings for all cells involved in `pair_list` (treated ∪ controls) via frozen embedding pipeline.
* For each treated:

  * `control_sum = Σ emb(control_cells)`
  * `control_counts = n_controls_used(treated)`  **(must be dynamic)**
  * Apply `delta_operator_interface(...)` with operator policy (`operator_type`, `eps`) to obtain:

    * `fm_delta` and `valid_mask`

Output per model:

* `k562/fm/<model>/fm_delta.npy`
* `k562/fm/<model>/fm_delta_meta.csv` includes:

  * `row_id`, `treated_cell_id`, `valid_mask`, `n_controls_used`

### 4.2.6 K562 overlap evidence (mandatory)

`k562_overlap_report.csv` reports overlap of treated universes across:

* Gene / Pathway
* each FM model after `valid_mask`

---

## 5. Group-level Mechanism Concordance (Main Analysis)

For each `mech_key=(dataset, cell_line, target_token)` and each `representation`:

### 5.1 Cohort definition

* Chemical side:

  * `X = { delta_i | chemical instance i, target_token ∈ target_tokens(i), delta_valid_bool=True }`
* Genetic side:

  * `Y = { delta_j | genetic instance j, target_token == target_token, delta_valid_bool=True }`

If either side empty → `mech_attrition.csv` with reason `missing_one_side`.

### 5.2 Metrics (frozen)

* `cosine = cosine(centroid(X), centroid(Y))`
* `pcc = PearsonCorr(centroid(X), centroid(Y))`

  * Always computed; may be NA for degenerate vectors; record `pcc_na_reason` if needed.
* `edist = bias_corrected_energy_distance(X_sub, Y_sub)`

  * Each side subsample ≤ 256 using single global RNG seed table (precomputed indices)
  * Uses squared Euclidean distance per spec.

Output: `task2_group_concordance.csv`
Must include denominators:

* `n_chem_instances_used`, `n_gen_instances_used`, `n_chem_sub`, `n_gen_sub`

---

## 6. Retrieval (R-INSTANCE; Two-direction; m_pos>1)

### 6.1 Common scoring and tie-handling

* score: `score(q,c) = -||q - c||^2`
* strict rank:

  * `rank = 1 + count(score_all > true_score)` (ties do not count as better)
* raw metrics:

  * `MRR_raw = 1/rank_true`
  * `Hit@K_raw = 1[rank_true <= K]` for K ∈ {1,5,10}

### 6.2 C2G (Chemical → Genetic)

**Query**: chemical treated instance delta (not duplicated)
**Gallery**: genetic target centroids within same `(dataset, cell_line)` scope

Positives:

* `Pos(i) = target_tokens(i) ∩ gallery_targets`
* `m_pos = |Pos(i)|` (can be >1)
* `rank_true = min_{t ∈ Pos(i)} rank(t)`

**Per-query must record** (explicit reviewer-facing metadata):

* `query_time`, `query_dose_value` (NA if missing)
* `query_target_raw`, `query_n_targets`
* `m_pos`, `N_gallery`

### 6.3 G2C (Genetic → Chemical; context-centroid gallery)

**Query**: genetic instance delta with target `t`
**Gallery items**: chemical **context centroids**, context key:

* `chem_context_id = (cell_line, perturbation_raw, time, dose_value)`

Each context centroid must carry:

* `context_target_tokens = union(target_tokens of chemical instances in this context)`
* `n_targets_context`

Positives:

* `Pos(t) = { context ∈ gallery : t ∈ context_target_tokens }`
* `m_pos = |Pos(t)|` (can be >1)
* `rank_true = min_{context ∈ Pos(t)} rank(context)`

### 6.4 Multi-positive chance correction (mandatory)

For each query:

* Use exact expectations with `(N_gallery, m_pos)`:

  * `expected_mrr = E[1/R | N_gallery, m_pos]`
  * `expected_hitK = P(R<=K | N_gallery, m_pos)`
* corrected:

  * `MRR_corr = MRR_raw - expected_mrr`
  * `Hit@K_corr = Hit@K_raw - expected_hitK`

### 6.5 Identity check (mandatory; tol=1e-12)

Output `chance_identity_check.csv` verifying:

* `MRR_corr == MRR_raw - expected_mrr` within tol
* `Hit@K_corr == Hit@K_raw - expected_hitK` within tol
  for all summary rows.

### 6.6 Retrieval outputs

* `task2_retrieval_per_query.parquet` (must include `m_pos`, `N_gallery`, and required metadata)
* `task2_retrieval_summary.csv` (denominators + means)
* `chance_identity_check.csv`

---

## 7. Dose/Time Reporting (Pre-binning Phase)

### 7.1 Distribution table (mandatory)

`task2_lincs_time_dose_distribution.csv` includes:

* unique `time`, `dose_value`, and `(time,dose_value)` counts
* at least global; optionally per cell_line

No binning is frozen at v1.3; bins may be set in later revision.

---

## 8. Significance Testing & Ecological Unit (Target-level)

### 8.1 Ecological unit declaration (frozen)

* Statistical testing is performed at **target-level (mech_key rows)**.
* Due to explode-matching, mech_keys may share chemical source instances; p/q are **not** interpreted as independent-instance inference.

### 8.2 Gene baseline (mandatory)

For each dataset × metric × direction:

* Compare `rep != Gene` vs `Gene` using paired tests:

  * group-level: paired by `mech_key`
  * retrieval: paired by query_id within same direction and cell_line scope (or by mech_key if query is mech_key-defined)

Multiplicity correction:

* BH-FDR within family `(dataset, test_scope, metric_name, direction)`.

Output: `task2_significance.csv`

---

## 9. specificity_tier & n_targets Stratification (Mandatory)

Output: `task2_tier_stratified_summary.csv`:

* Stratify chemical instances by:

  * `specificity_tier` (include `"NA"`)
  * `n_targets` bins: {1, 2–3, ≥4}
* Report both group-level and retrieval summaries with denominators:

  * `n_mech_keys`
  * `n_unique_chemical_instances`
  * `dist_m_pos` (for retrieval)

---

## 10. Denominator Transparency Summary (Mandatory)

Output: `task2_denominator_summary.csv` including:

* `n_mech_keys`
* `n_unique_chemical_instances`
* `n_unique_genetic_instances`
* `dist_n_targets` (chemical)
* `dist_m_pos` (C2G and G2C separately)
* `dist_N_gallery` (C2G and G2C separately)

---

## 11. Attrition Policy (Soft-fail)

Must output:

* `mech_attrition.csv` (missing_one_side, insufficient_subsample_for_edist, etc.)
* `k562_attrition.csv` (no_controls_available, FM_valid_mask_drop, etc.)

No run termination unless a fail-fast readiness condition is violated.

---

## 12. Final Output Checklist (v1.3)

**Core:**

1. `task2_group_concordance.csv`
2. `task2_retrieval_per_query.parquet`
3. `task2_retrieval_summary.csv`
4. `chance_identity_check.csv`
5. `task2_significance.csv`
6. `task2_tier_stratified_summary.csv`
7. `task2_denominator_summary.csv`
8. `task2_lincs_time_dose_distribution.csv`

**K562 evidence:**
9) `k562/derived/pair_list.parquet`
10) `k562/derived/delta_meta.csv`
11) `k562/derived/gene_delta.npy`
12) `k562/derived/pathway_delta.npy`
13) `k562/fm/<model>/*` for all required models
14) `k562_overlap_report.csv`
15) `k562_attrition.csv`
16) `mech_attrition.csv`
