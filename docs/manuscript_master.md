# Manuscript Master

Single canonical source-of-truth for manuscript structure, frozen figure logic, active analysis plan, implementation specs, source register, and writing guardrails.

Structured registry files kept separately: `manuscript_claim_to_evidence_map.csv`, `manuscript_panel_support_ledger.csv`.

---

## § 1 — Paper Type and Scope

- Target: Brief Communication / Analysis-style short paper. 2–3 main figures, compact claims structure. Supplementary carries denominator detail, sensitivity detail, audit context, and local-scope comparisons.
- Evidence chain: Task1 is authoritative through S2; corrected multisource Task2 is authoritative through S6; S7 is project-level synthesis implementation context only (not primary evidence for Figure 2 or Figure 3).
- Excluded as non-authoritative: S7 as primary evidence; legacy `data/task2_snapshot_v1/`; legacy scPerturb-only Task2 stages.

---

## § 2 — Primary Positioning

M2M-Bench is an evaluation and boundary-definition framework for transcriptome-centric virtual cell modeling. It is not a virtual cell model and not a predictor paper.

Core framing: M2M-Bench defines and tests the external-validity gap between within-modality predictive success and mechanism-level biological portability, by measuring modality concordance and mechanism concordance separately under explicit workflow, scope, and evidence contracts.

Central decomposition: **modality concordance** (Task1) vs. **mechanism concordance** (Task2).

---

## § 3 — Figure 2 Structure (Frozen)

Figure 2 = Task1 only. No Task2 bridge content in Figure 2. A1 bridge belongs here as the Figure 2 / Result 2 reference frame.

**Sub-panel scope (frozen):**

| Sub-panel | Content |
|---|---|
| F2-1 | Task1 internal evaluable scope |
| F2-2 | Task1 internal overall performance |
| F2-3 | Interpretation of Task1 internal performance |
| F2-4 | Task1 internal representation comparison: LINCS = Gene/Pathway only; scPerturb = Gene/Pathway (common scope, main text) + all available FM embeddings (supplementary panels only) |
| F2-5 | Task1 internal cell line / target preference |
| F2-6 | Task1 cross evaluable scope (genetic-only) |
| F2-7 | Task1 cross overall performance |
| F2-8 | Interpretation of Task1 cross performance |
| F2-9 | Task1 cross representation comparison: Gene/Pathway only; **NO FM in Task1 cross** |
| F2-10 | Task1 cross cell line / target preference |
| F2-11 | Internal→cross degradation summarized by dataset / cell_line deltas (A1 bridge output) |

Figure 2 interpretive role: show stable non-random signal in modality-preserving settings; establish the calibrated reference frame for reading Task2; do not use ceiling language for Task2.

---

## § 4 — Figure 3 Structure (Frozen)

Figure 3 = Task2 + Task1→Task2 bridge. Figure 3 covers Task2 itself and then Task1-to-Task2 bridge / mechanism transfer.

**Sub-panel scope (frozen):**

| Sub-panel | Content |
|---|---|
| F3-1 | Task2 evaluable scope |
| F3-2 | Overall performance for G2C and C2G — split in figures, unified in prose |
| F3-3 | Task2 representation comparison: LINCS = Gene/Pathway only; scPerturb = Gene/Pathway (common scope) + all available FM embeddings |
| F3-4 | Interpretation of Task2 performance |
| F3-5 | Chemical-side dose/time covariates and their effect on C2G retrieval against matched genetic target centroids (B6 output; C2G only; not symmetric matching; not benchmark-wide) |
| F3-6 | Task2 cell line / target preference |
| F3-7 | Task1-to-Task2 bridge and mechanism transfer (A2 output — belongs in Figure 3, not Figure 2) |

Figure 3 interpretive role: make Task2 the main external-validity and biological-interpretation layer; explain heterogeneity through covariates and case structure; test whether local FM scope changes conclusions or only local summaries.

---

## § 5 — Results Architecture

**Result 1** — Benchmark design as scientific contribution. Task1 = modality-concordance reference frame; Task2 = mechanism-concordance and external-validity layer. Establishes audited workflow, scope rules, evidence contract.

**Result 2** — Quantitative comparison of modality concordance and mechanism concordance. Requires A1 (Task1 internal-vs-cross) and A2 (Task1-vs-Task2) group-level bridges on shared `(dataset, cell_line, target)` groups.

**Result 3** — Biological and experimental covariates explain Task2 heterogeneity. Includes target complexity, C2G-only dose/time downstream sensitivity (B6), enrichment-style readouts, local specificity where frozen scope supports it.

**Result 4** — Conditional representation effects within supported scPerturb K562 FM scope. Tests whether FM changes group conclusions, retrieval conclusions, or only local target-level summaries.

---

## § 6 — Core Claims

**C1** — M2M-Bench contributes a benchmark design that defines the problem boundary, workflow, scope, audit path, and evidence contract for transcriptome-centric virtual cell evaluation.

**C2** — M2M-Bench exposes an external-validity gap: within-modality predictive success does not automatically transfer to mechanism-level biological portability, and this gap becomes visible only when Task1 and Task2 are evaluated separately and quantitatively.

**C3** — Task2 heterogeneity is biologically and experimentally structured. Interpretation depends on dataset, cell line, target complexity, supported specificity context, and covariate-audited experimental context rather than a single mechanism-concordance summary.

**C4** — Foundation-model effects are local and conditional: within the supported scPerturb K562 scope, the key question is whether FM changes group conclusions, retrieval conclusions, or only local target-level summaries.

---

## § 7 — Active Analysis Inventory

| ID | Analysis | Status | Figure role | Output path | Active boundary |
|---|---|---|---|---|---|
| A1 | Task1 internal-vs-cross group bridge | **Materialized** | Figure 2 / Result 2 reference frame | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/group_bridge/task1_internal_vs_cross_group_bridge.csv` | Shared `(dataset, cell_line, target)` groups; Gene/Pathway only; genetic-only cross |
| A2 | Task1-vs-Task2 group bridge | **Materialized** | **Figure 3 / Task1→Task2 bridge** | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/group_bridge/task1_task2_group_bridge.csv` | Shared `(dataset, cell_line, target)` groups |
| B6 | C2G dose/time downstream sensitivity | **Materialized** | Figure 3 / Result 3 covariate interpretation | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/sensitivity/task2_c2g_dose_time_sensitivity.csv` | C2G only; downstream sensitivity only; not benchmark-wide |
| A3 | Direction robustness appendix | Accepted | Appendix / transparency support | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a3_direction_robustness/task2_direction_robustness_audit.csv` | Appendix / transparency scope only; do not widen |
| B5 | n_targets sensitivity | Partial | Result 3 downstream sensitivity | `docs/task2_c2g_query_n_targets_sensitivity.csv` | C2G only; downstream sensitivity only; not benchmark-wide |
| B8 | Specificity-tier sensitivity | Partial | Local appendix / reviewer-defense | `docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv` | scPerturb K562 C2G only; local-only; not benchmark-wide |
| R1 | Claim/evidence registration refresh | Planned | Register support tables before drafting | `manuscript_claim_to_evidence_map.csv`; `manuscript_source_register` section | Update only after B5/B8 tables exist locally |

Execution order: A1 → A2 → B6 → A3 → R1 → drafting.

---

## § 8 — Implementation Specs (A1 / A2 / B6)

### A1

- Goal: group-level Task1 internal-versus-cross bridge on shared `(dataset, cell_line, target)` groups.
- Inputs: Task1 internal and cross outputs (S1/S2); last-resort audited S1 group export if current outputs are insufficient.
- Common-scope: Gene and Pathway only; genetic-only cross.
- Output schema: `dataset`, `cell_line`, `target`, agreed metric fields, internal values, cross values, deltas, denominator/support fields, scope notes.
- Validation: required-column assertions; lawful common-scope filtering; denominator checks against frozen source tables.

### A2

- Goal: group-level Task1-vs-Task2 bridge on shared `(dataset, cell_line, target)` groups.
- Inputs: Task1 (S1/S2) and Task2 (S4/S5/S6) outputs on shared groups.
- Output schema: `dataset`, `cell_line`, `target`, agreed metric fields, Task1 values, Task2 values, deltas, denominator/support fields, scope notes.
- Validation: required-column assertions; lawful shared-group validation; denominator checks.

### B6

- Goal: downstream chemical-side dose/time sensitivity for Result 3 covariate interpretation.
- Required scope: `direction=C2G` only.
- Inputs: S5 per-query retrieval; LINCS `delta_meta.csv`; scPerturb `delta_meta.csv`.
- Output schema: `query_row_id`, `target`, `dataset`, `cell_line`, `direction`, `representation`, `time`, `dose_value`, `covariate_slice_id`, `n_queries_in_covariate_slice`, corrected retrieval metric, `sensitivity_only_bool`.
- Validation: metadata join uniqueness; explicit accounting for missing covariates; denominator checks.

### Global coding principles

- Downstream manuscript-support only; no upstream stage mutation.
- Explicit local-only and sensitivity-only labeling.
- No silent coercion, no silent join loss, no silent fillna.
- Status vocabulary: `successful` / `partial` / `unavailable` / `blocked`.

---

## § 9 — Source Register (Condensed)

Authoritative evidence path for Figure 2 and Figure 3:

### Task1 S1 Internal — run `s1_task1_internal_metrics_0303`

| File | Role | Key caution |
|---|---|---|
| `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_leaderboard_long.csv` | Primary Task1 internal group summary; internal side of A1 | Keep scPerturb Gene/Pathway rows distinct from FM rows |
| `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_retrieval_summary.csv` | Task1 internal retrieval summary | Common-scope and FM-local retrieval slices must be indexed separately |
| `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_attrition.csv` | Denominator transparency / exclusion support | Support notes only |
| `runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics/task1_chance_identity_check.csv` | Corrected-metric validation | Validation only; not a primary claim table |

### Task1 S2 Cross — run `s2_task1_cross_metrics_0303`

| File | Role | Key caution |
|---|---|---|
| `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_leaderboard_long.csv` | Primary Task1 cross summary; cross side of A1 | Genetic-only; cross chemical excluded by contract |
| `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_retrieval_summary.csv` | Directional cross retrieval | Keep LINCS_to_scPerturb and scPerturb_to_LINCS separate |
| `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_attrition.csv` | Cross support-loss table | Chemical exclusion: `cross_matched_keys_lt_5`; do not quote numeric rationale |
| `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_alignment_proof.csv` | Alignment-contract eligibility proof | Chemical row: `n_matched_keys=0`; use only to state exclusion by policy/support gating |
| `runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics/task1_cross_chance_identity_check.csv` | Corrected cross retrieval validation | Validation only |
| `docs/contracts/exclusions-and-policies.md` | Frozen policy language for cross chemical exclusion | Policy wording only; do not quote numeric rationale from policy prose |

### Task2 S4 Group Concordance — run `s4_multisource_impl_verify_20260310_c`

| File | Role | Key caution |
|---|---|---|
| `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_concordance.csv` | Authoritative target-level mechanism-concordance; Task2 side of A2 | No `specificity_tier`; no `n_targets` cohort field |
| `runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource/task2_group_attrition.csv` | Task2 group support-loss | UCE attrition in scPerturb K562 must be disclosed |

### Task2 S5 Retrieval — run `s5_multisource_impl_verify_20260311_a`

| File | Role | Key caution |
|---|---|---|
| `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary.csv` | Authoritative direction-specific retrieval summary; downstream sensitivity inputs | Keep C2G and G2C separate; carry `gallery_definition_id` and `pos_definition_id` |
| `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_summary_long.csv` | Long-format retrieval for metric-wise extraction | Corrected only; resolve correction questions back to summary |
| `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_attrition.csv` | Retrieval support-loss | UCE valid-mask drops in scPerturb K562 |
| `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_chance_identity_check.csv` | Corrected Task2 retrieval validation | Validation only |
| `runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource/task2_retrieval_per_query.parquet` | Query-level retrieval; inputs for B5/B6/B8 | Any dose/time or specificity use stays downstream sensitivity-only |

### Task2 S6 Result Synthesis — run `s6_multisource_impl_verify_20260311_a`

| File | Role | Key caution |
|---|---|---|
| `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_leaderboard.csv` | Primary Figure 3 group-level summary; local scPerturb K562 FM group panel | Ranking within same dataset/cell_line only; UCE attrition must be disclosed |
| `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_retrieval_leaderboard.csv` | Primary Figure 3 retrieval summary; local FM retrieval support | Keep C2G and G2C separate |
| `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_group_concordance_long.csv` | Long-format group table for metric-wise extraction | Target-level only; not the manuscript ranking table |
| `runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource/task2_benchmark_summary_long.csv` | Synthesis union table for panel-readiness indexing | Summary only; not a license to re-aggregate |

### Task2 Snapshot v2 Metadata

| File | Role | Key caution |
|---|---|---|
| `data/task2_snapshot_v2/task2_pairs_coverage.csv` | Eligibility and coverage fact table | Coverage transparency only; not a metric table |
| `data/task2_snapshot_v2/representation_availability_registry.csv` | Representation-scope registry | LINCS FM absence = not_applicable_scope (not attrition) |
| `data/task2_snapshot_v2/lincs/derived/delta_meta.csv` | LINCS query-level metadata | `specificity_tier=NA` for LINCS; benchmark-wide specificity claim not supported |
| `data/task2_snapshot_v2/scperturb_k562/derived/delta_meta.csv` | scPerturb K562 query-level metadata | FM scope exists only here; any specificity use stays local K562 C2G |

---

## § 10 — Claim Guardrails, Retired Axes, and Stop Rules

### Core guardrails

- **FM scope**: FM evidence is local to scPerturb K562 only. Do not generalize FM results to LINCS or benchmark-wide behavior.
- **LINCS FM absence**: LINCS does not support FM comparison. LINCS FM absence is a scope policy (`not_applicable_scope`), not attrition.
- **Task1 cross**: genetic-only in current contract. Cross-chemical exclusion is a contract outcome, not a missing empirical result. State only that it is excluded by frozen policy/support gating; do not quote a numeric matched-context rationale.
- **Task1 not a ceiling**: Task1 is a calibrated reference frame for Task2, not a formal upper bound.
- **No raw dataset pooling**: Task2 core evidence remains stratified by dataset and cell_line; raw pooling across LINCS and scPerturb is not a core-result interpretation.
- **Retrieval contract**: keep C2G and G2C separate; carry `gallery_definition_id` and `pos_definition_id` with every retrieval claim.
- **scPerturb FM vs. common-scope**: scPerturb internal Gene/Pathway rows must remain distinct from the additional FM rows in S1 outputs.
- **S7**: implementation context only; must not replace Task1 S1/S2 or corrected Task2 S4/S5/S6 as evidence.

### Retired axes (do not promote to benchmark-wide core panels)

- **`specificity_tier`** is retired as a benchmark-wide or group-level core Figure 3 axis. Reason: absent from authoritative S4–S6 summaries; LINCS `specificity_tier=NA` means benchmark-wide specificity cannot be supported from both datasets. Route only to B8: local scPerturb K562 C2G sensitivity.
- **`n_targets`** is retired as a benchmark-wide bidirectional core panel. Reason: G2C is structurally non-informative (all queries have `query_n_targets=1`); group-level complexity aggregation rule is not frozen. Route only to B5: downstream C2G sensitivity.
- **dose/time**: downstream C2G-only sensitivity readout for covariate interpretation. The question is whether chemical-side dose/time covariates affect C2G retrieval success against matched genetic target centroids. This is **not** symmetric chemical-genetic matching. Not benchmark-wide.
- **B5/B8**: remain local/sensitivity-only; do not promote to benchmark-wide claims.

### Stop rules

Stop any remediation or analysis expansion if it would:
1. Require rerunning major upstream benchmark stages (S1–S6).
2. Redefine Task1 or Task2 scope, task semantics, or benchmark contracts.
3. Raw-pool LINCS and scPerturb for core Task2 claims.
4. Generalize FM evidence beyond the supported scPerturb K562 scope.
5. Invent a benchmark-wide aggregation rule not frozen locally.
6. Quote a numeric rationale that lacks frozen local source support.

---

## § 11 — Support Table Schemas (Condensed)

Active and pending downstream manuscript-support tables. These are downstream artifacts only; they do not replace S1–S6 authoritative outputs.

| Table | Status | Scope | Key columns | Source inputs |
|---|---|---|---|---|
| `/mnt/.../group_bridge/task1_internal_vs_cross_group_bridge.csv` | **Materialized** | A1 bridge; Gene/Pathway; genetic cross only | `dataset`, `cell_line`, `target`, `metric_name`, `representation`, `internal_value`, `cross_value`, `cross_minus_internal`, `internal_n_queries`, `cross_n_queries`, `provenance_note` | Task1 S1/S2 shared groups |
| `/mnt/.../group_bridge/task1_task2_group_bridge.csv` | **Materialized** | A2 bridge | `dataset`, `cell_line`, `target`, `metric_family`, `task1_value`, `task2_value`, `task2_minus_task1`, `task1_support_n`, `task2_support_n`, `provenance_note` | Task1 S1/S2; Task2 S4/S6 shared groups |
| `/mnt/.../sensitivity/task2_c2g_dose_time_sensitivity.csv` | **Materialized** | B6; C2G only; not benchmark-wide | `direction`, `representation`, `target`, `time`, `dose_value`, `covariate_slice_id`, `n_queries_in_covariate_slice`, `retrieval_value`, `sensitivity_only_bool` | S5 per-query; snapshot delta_meta |
| `/mnt/.../a3_direction_robustness/task2_direction_robustness_audit.csv` | Accepted | A3; appendix/transparency | `dataset`, `cell_line`, `direction`, `representation`, `mean_mrr_corrected`, `n_queries`, `N_gallery_mean`, `support_note` | Existing A3 output |
| `docs/task2_c2g_query_n_targets_sensitivity.csv` | Partial | B5; C2G only | `dataset`, `cell_line`, `representation`, `query_n_targets`, `direction`, `n_queries`, `mean_mrr_corrected`, `sensitivity_only_bool` | S5 per-query; target metadata |
| `docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv` | Partial | B8; scPerturb K562 C2G only | `representation`, `specificity_tier`, `dataset`, `cell_line`, `direction`, `n_queries`, `mean_mrr_corrected`, `local_only_bool`, `sensitivity_only_bool` | S5 per-query; scPerturb K562 metadata |

Additional planned tables (pending B5/B8 completion): `task2_n_targets_pairwise_stats.csv`, `task2_specificity_local_stats.csv`, FM local comparison and robustness tables. See original support_table_specs for full schemas if needed.

---

## § 12 — Literature Positioning and Overclaim Guardrails

### Safe positioning

- M2M-Bench is an evaluation and boundary-definition framework for transcriptome-centric virtual cell modeling, not a virtual cell model.
- Task1 provides a modality-preserving reference frame; Task2 tests the mechanism-concordance layer.
- M2M-Bench defines a benchmark boundary for what current perturbation-response modeling can support under audited, scope-limited evaluation.

### How M2M-Bench complements existing work

- **Arc / Virtual Cell Challenge**: M2M-Bench contributes an evaluation layer separating easier modality-preserving agreement from harder portability across ecosystem and mechanism boundaries.
- **Systema**: M2M-Bench turns concern about metric inflation into an explicit benchmark design question by separating Task1 reference-frame evidence from Task2 mechanism-concordance evidence.
- **Wei et al.**: M2M-Bench is model-agnostic; rather than ranking prediction methods, it separates calibration-style modality concordance from the mechanism-concordance layer constraining transfer claims.
- **scDrugMap**: M2M-Bench makes dataset-stratified interpretation, scope-limited FM reading, and mechanism-versus-modality separation central to the benchmark design.

### Do NOT write

- `M2M-Bench is a virtual cell model`
- `FM results generalize benchmark-wide`
- `modality concordance guarantees mechanism concordance`
- `M2M-Bench shows that current virtual cell methods do not work`
- `Task1 provides an upper bound for Task2`

### Write instead

- `M2M-Bench is an evaluation and boundary-definition framework for transcriptome-centric virtual cell modeling.`
- `Foundation-model findings are interpretable only within the supported scPerturb K562 scope.`
- `Modality concordance and mechanism concordance should be evaluated separately.`
- `Task1 serves as a reference frame for Task2, not as a formal ceiling.`

---

## § 13 — Discussion and Limitation Boundaries

The discussion must stay inside benchmark-evaluation claims. Do not drift into predictor novelty, mechanism causality, or platform-agnostic transfer claims that exceed the current audited contracts and evidence.

Required limitation statements:
- The manuscript does not become a generative virtual-cell benchmark in this revision.
- The manuscript does not depend on zero-shot or OOD generation tasks for validity; those tasks are outside the current audited scope.
- The current corrected retrieval framework remains in force; asymmetry analyses are robustness checks, not a metric reset.
- Corrected Task2 is multisource, but its core results remain stratified by dataset and cell_line.
- FM evidence is limited to the scPerturb K562 subset.
- Local specificity, target-complexity, and dose/time interpretations must stay inside their supported audited scopes.
- Enrichment claims remain constrained to matched/evaluable target pools.
- Legacy scPerturb-only Task2 outputs are historical evidence only.
