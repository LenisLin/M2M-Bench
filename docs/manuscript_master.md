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

## § 3 — Current Canonical Downstream Object Contract

Authoritative current-phase registry:
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/framework_analysis_manifest.json`

Canonical downstream object root:
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/`

| Figure | Canonical object | Current role | Notes |
|---|---|---|---|
| Figure 2 | `figure2_task1_scope_summary.csv` | Task1 evaluable-scope summary | Covers Task1 internal and Task1 cross frozen scope boundaries |
| Figure 2 | `figure2_task1_performance_structure.csv` | Task1 performance-structure summary | Carries Task1 group and retrieval structure within frozen scope |
| Figure 2 | `figure2_task1_internal_to_cross_degradation_summary.csv` | Task1 internal-to-cross degradation summary | Figure 2 only; downstream A1 degradation summary |
| Figure 3 | `figure3_task2_scope_summary.csv` | Task2 evaluable-scope summary | Dataset/cell_line/representation scope with FM policy boundary |
| Figure 3 | `figure3_task2_direction_support_summary.csv` | Task2 direction-support summary | Canonical direction-specific retrieval support object |
| Figure 3 | `figure3_task2_performance_structure.csv` | Task2 performance-structure summary | Canonical benchmark performance structure object |
| Figure 3 | `figure3_task2_cell_line_pattern_summary.csv` | Task2 cell_line pattern summary | Canonical Task2 pattern object |
| Figure 3 | `figure3_task2_target_pattern_summary.csv` | Task2 target pattern summary | Canonical Task2 pattern object |
| Figure 3 | `figure3_task1_internal_contextual_support_summary.csv` | Task1 internal contextual support summary | Figure 3 uses Task1 internal support only; no Task1 cross and no A2 pairwise bridge |

Figure 1 source note:
- There is no dedicated canonical Figure 1 analysis object in the frozen current phase.
- Figure 1 is doc-derived from this file plus `docs/contracts/project-positioning.md`, `docs/contracts/task1_spec.md`, `docs/contracts/task2_spec.md`, and governance prose.
- Missing Figure 1 CSV/JSON outputs are therefore not repo debt in this cleanup phase.
- The plotting-preparation freeze for Figure 1 registries and the full main/extended panel manifest is tracked in `docs/plotting/plotting_preparation_freeze.md`.
- The support-only Figure 1 registry builder is `scripts/manuscript_plot_ready_tables.py`, which writes to `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/plot_ready/figure1/`.

Figure 1 final layout:
- Main: `1B` hierarchical benchmark composition panel
- Deferred: workflow schematic prompt and data-filtering schematic prompt
- Supplementary / methods support: `1C` lawful scope matrix and `1D`
  representation/modifier availability
- Deleted from the main figure: `1E`

Retained but non-canonical historical surfaces include:
- A2 Task1↔Task2 pairwise bridge tables and A2 overview files
- standalone C2G-oriented Figure 3 covariate/support files
- B6 overview files and C2G sensitivity detail/summary files
- the merged Task2 descriptive/preference-support placeholder file
- legacy manuscript-history manifest and figure-input note files

---

## § 4 — Figure 2 Structure (Frozen)

Figure 2 = Task1 only. Task1 cross remains canonical here only. No Figure 3 bridge wording belongs in Figure 2.

| Figure 2 panel role | Canonical object | Frozen scope note |
|---|---|---|
| Task1 evaluable scope | `figure2_task1_scope_summary.csv` | Task1 internal keeps FM where lawful; Task1 cross remains Gene/Pathway only |
| Task1 performance structure | `figure2_task1_performance_structure.csv` | Internal scPerturb FM remains local to Task1 internal; any LINCS-involving slice excludes FM |
| Internal→cross degradation | `figure2_task1_internal_to_cross_degradation_summary.csv` | Downstream A1 summary only; not a Figure 3 bridge |

Figure 2 contract notes:
- Task1 cross is canonical for Figure 2 only.
- `2A` is a stacked scope bar panel, not a tile matrix.
- `2B` is the compact internal scoreboard and should retain both LINCS and
  scPerturb Gene/Pathway rows; scPerturb FM rows remain local-scope only.
- Any scPerturb-only comparison includes FM.
- Any comparison involving LINCS excludes FM.
- Figure 2 does not use Task1↔Task2 pairwise bridge language.
- Figure 2 final panel list is `2A` stacked Task1 scope, `2B` scoreboard-style
  internal benchmark matrix, `2C` unified matched-unit Gene-vs-Pathway comparison,
  `2D` internal-to-cross trend curves, `2E` paired cell-line enrichment exemplars,
  and `2F` paired target-enrichment exemplars.
- `2E` and `2F` consume support-only `plot_ready/` reductions built in Python;
  underpowered filtering is mandatory before R plotting, the main panels are
  internal-only, and each dataset contributes `3` high + `3` low exemplars
  with shared LINCS/scPerturb entities highlighted on the shared axis.

---

## § 5 — Figure 3 Structure (Frozen)

Figure 3 = Task2 + Task1 internal contextual support only.

| Figure 3 panel role | Canonical object | Frozen scope note |
|---|---|---|
| Task2 evaluable scope | `figure3_task2_scope_summary.csv` | Current lawful Task2 scope by dataset, cell_line, and representation |
| Task2 direction support | `figure3_task2_direction_support_summary.csv` | Direction-specific retrieval support; canonical and separate from performance structure |
| Task2 performance structure | `figure3_task2_performance_structure.csv` | Canonical benchmark performance structure object |
| Task2 cell_line pattern summary | `figure3_task2_cell_line_pattern_summary.csv` | Canonical Task2 pattern surface |
| Task2 target pattern summary | `figure3_task2_target_pattern_summary.csv` | Canonical Task2 pattern surface; retrieval rows are C2G direction aggregated by canonical target; G2C target-anchored retrieval absent |
| Task1 internal contextual support | `figure3_task1_internal_contextual_support_summary.csv` | Internal-only contextual support; excludes Task1 cross and any Task1↔Task2 pairwise bridge |

Figure 3 contract notes:
- C2G is the primary direction axis for Figure 3 main-text reading in the frozen current phase.
- G2C remains a supporting direction-specific axis and should not be framed as a co-headline replacement for C2G.
- No Task1↔Task2 pairwise bridge object is canonical in this phase.
- No standalone C2G-side experiment object is canonical in this phase.
- Historical C2G sensitivity artifacts may remain on disk but are non-canonical.
- `specificity_tier` remains retired as a modifier surface (no formal new outputs). Historical B6 dose/time artifacts and B5 n_targets sensitivity tables remain non-canonical. The F3.6–F3.8 formal modifier result families (time, dose, target multiplicity) have been materialized as supplementary result objects; they are not retired placeholders (see §10 and §11).
- Any scPerturb-only comparison includes FM.
- Any comparison involving LINCS excludes FM.
- Figure 3 final main panel list is `3B` scoreboard-style common-scope benchmark matrix,
  `3C` LINCS cell-line enrichment exemplars, `3D` target enrichment exemplars, `3E` paired Gene-vs-Pathway boxplots,
  and `3F` FM-local metric trade-off scatter panels.
- The Task2 scope summary remains support/reference data only; it is no longer rendered as a main `3A` panel.
- `3E` is backed by `figure3_panel_3e_gene_vs_pathway_paired.csv`,
  and `3F` is backed by `figure3_panel_3f_fm_local_tradeoff.csv`.
- `3F` remains `scPerturb/K562` local only.
- `3J` is background/reference only and is moved out of the main figure to
  supplementary/reference support.
- `3C` and `3D` are rendered as paired enrichment exemplar panels derived from
  the authoritative pattern summaries rather than as percentile heatmaps.
- `3C` main-text scope is LINCS-only; scPerturb K562 remains represented in the
  `3F` local panel rather than being forced into the `3C` exemplar ranking.

---

## § 6 — Results Architecture

**Result 1** — Benchmark design as scientific contribution. Task1 = modality-concordance reference frame; Task2 = mechanism-concordance and external-validity layer. Establishes audited workflow, scope rules, evidence contract.

**Result 2** — Figure 2 Task1 reference frame using the current canonical Task1 scope, performance-structure, and internal-to-cross degradation objects.

**Result 3** — Figure 3 Task2 scope, performance-structure, and pattern
summaries, with paired common-scope comparison and FM-local trade-off reading
kept separate from the supplementary Task1 contextual reference surface.

**Result 4** — Conditional representation effects within supported scPerturb K562 FM scope. FM wording remains local-scope only and does not promote historical C2G sensitivity artifacts to canonical status.

---

## § 7 — Core Claims

**C1** — M2M-Bench contributes a benchmark design that defines the problem boundary, workflow, scope, audit path, and evidence contract for transcriptome-centric virtual cell evaluation.

**C2** — M2M-Bench exposes an external-validity gap: within-modality predictive success does not automatically transfer to mechanism-level biological portability, and this gap becomes visible only when Task1 and Task2 are evaluated separately and quantitatively.

**C3** — Task2 heterogeneity is biologically and experimentally structured. Interpretation depends on dataset, cell line, target complexity, supported specificity context, and covariate-audited experimental context rather than a single mechanism-concordance summary.

**C4** — Foundation-model effects are local and conditional: within the supported scPerturb K562 scope, the key question is whether FM changes group conclusions, retrieval conclusions, or only local target-level summaries.

---

## § 8 — Manuscript-facing Script Surface

| Script surface | Current role | Status | Notes |
|---|---|---|---|
| `scripts/manuscript_framework_analysis_objects.py` | Authoritative current-phase canonical object builder and registry writer | Live canonical | Declares the canonical object set and historical/non-canonical inventory |
| `scripts/manuscript_a1_task1_internal_cross.py` | Figure 2 A1 support builder for internal→cross degradation inputs | Live support only | Supports Figure 2 only |
| `scripts/manuscript_task1_group_support.py` | Shared Task1 helper module for A1 and retained historical bridge builders | Live helper | No canonical output role by itself |
| `scripts/manuscript_a3_direction_robustness.py` | Appendix/transparency audit builder | Live appendix support | A3 remains appendix/transparency only |
| `scripts/manuscript_local_support_tables.py` | Supplementary B5/B8 tables and panel-support ledger refresh | Live supplementary support | B5/B8 are non-canonical supplementary support |
| `scripts/manuscript_a2_task1_task2_bridge.py` | Historical A2 pairwise bridge rebuilder | Historical/non-canonical | Retained only so historical A2 artifacts can be reconstructed if needed |
| `scripts/manuscript_b6_c2g_dose_time_sensitivity.py` | Historical B6 C2G sensitivity rebuilder | Historical/non-canonical | Retained only for historical support reconstruction |
| `scripts/manuscript_figure_analysis_inputs.py` | Historical figure-input backfill helper | Historical/non-canonical | Produces retained overview files only; not the canonical builder |

---

## § 9 — Source Register (Condensed)

Authoritative manuscript-facing evidence roots use NAS-backed paths only.

| Source surface | Current role | Path |
|---|---|---|
| Canonical manuscript registry | Frozen current-phase canonical/historical inventory | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/framework_analysis_manifest.json` |
| Task1 S1 internal root | Authoritative Task1 internal evidence | `/mnt/NAS_21T/ProjectData/M2M/runs/s1_task1_internal_metrics_0303/s1_task1_internal_metrics` |
| Task1 S2 cross root | Authoritative Task1 cross evidence | `/mnt/NAS_21T/ProjectData/M2M/runs/s2_task1_cross_metrics_0303/s2_task1_cross_metrics` |
| Task2 S4 group root | Authoritative Task2 group evidence | `/mnt/NAS_21T/ProjectData/M2M/runs/s4_multisource_impl_verify_20260310_c/s4_task2_group_concordance_multisource` |
| Task2 S5 retrieval root | Authoritative Task2 retrieval evidence | `/mnt/NAS_21T/ProjectData/M2M/runs/s5_multisource_impl_verify_20260311_a/s5_task2_retrieval_multisource` |
| Task2 S6 synthesis root | Authoritative Task2 synthesis evidence | `/mnt/NAS_21T/ProjectData/M2M/runs/s6_multisource_impl_verify_20260311_a/s6_task2_result_synthesis_multisource` |
| Task2 snapshot registry | Frozen Task2 scope and FM policy registry | `data/task2_snapshot_v2/` |
| Manuscript phase1 appendix root | Appendix/transparency audit outputs | `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1` |

No current manuscript-facing source or registry surface should depend on repo-root `runs/`.

---

## § 10 — Claim Guardrails, Retired Axes, and Stop Rules

### Core guardrails

- **Figure 1 source model**: Figure 1 is doc-derived from contracts/governance/manuscript prose. Do not invent a fake Figure 1 analysis object during cleanup.
- **FM scope**: FM evidence is local to scPerturb K562 only. Do not generalize FM results to LINCS or benchmark-wide behavior.
- **LINCS FM absence**: LINCS does not support FM comparison. LINCS FM absence is a scope policy (`not_applicable_scope`), not attrition.
- **Task1 cross**: genetic-only in current contract. Cross-chemical exclusion is a contract outcome, not a missing empirical result. State only that it is excluded by frozen policy/support gating; do not quote a numeric matched-context rationale.
- **Figure 3 Task1 support**: Figure 3 uses Task1 internal contextual support only. Do not reintroduce Task1 cross or Task1↔Task2 pairwise bridge language into the live Figure 3 surface.
- **Task1 not a ceiling**: Task1 is a calibrated reference frame for Task2, not a formal upper bound.
- **Biological indexing vs code keying**: manuscript-facing biological indexing is `(dataset, cell_line, target, perturbation_type)`, while current Task2 S3-S6 cohort code remains keyed on `mech_key=(dataset, cell_line, target_token)`. Cleanup documents this distinction and does not silently reinterpret S3-S6 semantics.
- **No raw dataset pooling**: Task2 core evidence remains stratified by dataset and cell_line; raw pooling across LINCS and scPerturb is not a core-result interpretation.
- **Retrieval contract**: keep C2G and G2C separate; carry `gallery_definition_id` and `pos_definition_id` with every retrieval claim.
- **scPerturb FM vs. common-scope**: scPerturb internal Gene/Pathway rows must remain distinct from the additional FM rows in S1 outputs.
- **S7**: implementation context only; must not replace Task1 S1/S2 or corrected Task2 S4/S5/S6 as evidence.

### Retired axes (do not promote to benchmark-wide core panels)

- **`specificity_tier`** is retired as a benchmark-wide or group-level core Figure 3 axis. Reason: absent from authoritative S4–S6 summaries; LINCS `specificity_tier=NA` means benchmark-wide specificity cannot be supported from both datasets. Route only to B8: local scPerturb K562 C2G sensitivity.
- **`n_targets`** is retired as a benchmark-wide bidirectional core panel. Reason: G2C is structurally non-informative (all queries have `query_n_targets=1`); group-level complexity aggregation rule is not frozen. Route only to B5: downstream C2G sensitivity.
- **dose/time (B6 historical)**: the old B6 historical dose/time sensitivity artifacts (`sensitivity/task2_c2g_dose_time_sensitivity*.csv`) remain non-canonical. They are **not** current canonical Figure 3 objects, not symmetric chemical-genetic matching, and not benchmark-wide.
- **F3.6–F3.8 formal modifier result families** (distinct from B6 historical artifacts): F3.6 (time effect), F3.7 (dose effect), and F3.8 (target multiplicity) have been materialized as formal supplementary result objects (generated 2026-03-24) by `scripts/manuscript_task2_c2g_modifiers.py`. These are **not** retired placeholders. They are formal manuscript question families, C2G-only and not benchmark-wide. See §11 for paths.
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

## § 11 — Support and Historical Table Surface (Condensed)

Current canonical manuscript-facing downstream objects:

| Table | Current status | Current role |
|---|---|---|
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure2_task1_scope_summary.csv` | Canonical | Figure 2 Task1 scope summary |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure2_task1_performance_structure.csv` | Canonical | Figure 2 Task1 performance structure |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure2_task1_internal_to_cross_degradation_summary.csv` | Canonical | Figure 2 internal→cross degradation summary |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_scope_summary.csv` | Canonical | Figure 3 Task2 scope summary |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_direction_support_summary.csv` | Canonical | Figure 3 Task2 direction support |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_performance_structure.csv` | Canonical | Figure 3 Task2 performance structure |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_cell_line_pattern_summary.csv` | Canonical | Figure 3 Task2 cell_line pattern summary |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_target_pattern_summary.csv` | Canonical | Figure 3 Task2 target pattern summary |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task1_internal_contextual_support_summary.csv` | Canonical | Figure 3 Task1 internal contextual support |

Retained supplementary or historical surfaces:

Downstream plotting support layer:
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/plot_ready/` is the support-only plotting-preparation root.
- These plot-ready tables are reductions/registries over frozen canonical or support objects and do not create new canonical evidence.
- The declared plot-ready contract is frozen in `docs/plotting/plotting_preparation_freeze.md` and `docs/contracts/output-schemas.md`.
- Repo-local visual staging under `tmp/manuscript_visual_revision_v2/` may host support-only diagnostic or exemplar panels that clarify biology-structured heterogeneity without changing the canonical evidence roots.
- `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/` is the retained historical manuscript root for backfills, bridge leftovers, and retired notes.

| Table | Current status | Current role |
|---|---|---|
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_c2g_time_effect_summary.csv` | Supplementary/formal modifier result | F3.6 time-effect modifier result family — C2G only; LINCS 6h/24h primary; scPerturb 24h-only (not testable) |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_c2g_time_effect_stats.csv` | Supplementary/formal modifier result | F3.6 time-effect statistical tests (Kruskal-Wallis + Spearman trend; BH-FDR) |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_c2g_dose_effect_summary.csv` | Supplementary/formal modifier result | F3.7 dose-effect modifier result family — C2G only; 5 log10 bins; scPerturb 2 bins only; dose unit assumed µM |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_c2g_dose_effect_stats.csv` | Supplementary/formal modifier result | F3.7 dose-effect statistical tests (Kruskal-Wallis + Spearman trend; BH-FDR) |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_c2g_target_multiplicity_summary.csv` | Supplementary/formal modifier result | F3.8 target-multiplicity modifier result family — C2G only; bins {1, 2, 3+}; G2C excluded (structurally single-target) |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_c2g_target_multiplicity_stats.csv` | Supplementary/formal modifier result | F3.8 target-multiplicity statistical tests (Mann-Whitney U + Kruskal-Wallis; BH-FDR) |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis/figure3_task2_c2g_modifier_comparison_statistics.csv` | Supplementary/formal modifier result | Consolidated BH-FDR comparison statistics across F3.6/F3.7/F3.8 modifier families |
| `docs/task2_c2g_query_n_targets_sensitivity.csv` | Supplementary/non-canonical | B5 local C2G-only sensitivity support |
| `docs/task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv` | Supplementary/non-canonical | B8 local scPerturb K562 C2G-only sensitivity support |
| `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_phase1/a3_direction_robustness/task2_direction_robustness_audit.csv` | Appendix/transparency | A3 direction robustness audit |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure2_a1_metric_overview.csv` | Historical/non-canonical | Legacy Figure 2 A1 overview |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure2_a1_dataset_cellline_overview.csv` | Historical/non-canonical | Legacy Figure 2 A1 overview |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/group_bridge/task1_task2_group_bridge.csv` | Historical/non-canonical | Retained A2 pairwise bridge detail |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/group_bridge/task1_task2_group_bridge_summary.csv` | Historical/non-canonical | Retained A2 pairwise bridge summary |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_a2_bridge_support_matrix.csv` | Historical/non-canonical | Retained A2 bridge support matrix |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/sensitivity/task2_c2g_dose_time_sensitivity.csv` | Historical/non-canonical | Retained B6 sensitivity detail |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/sensitivity/task2_c2g_dose_time_sensitivity_summary.csv` | Historical/non-canonical | Retained B6 sensitivity summary |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_a2_metricpair_overview.csv` | Historical/non-canonical | Retained A2 overview |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_a2_dataset_cellline_overview.csv` | Historical/non-canonical | Retained A2 overview |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_task2_c2g_covariate_summary.csv` | Historical/non-canonical | Retained standalone C2G summary |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_b6_metric_dose_time_overview.csv` | Historical/non-canonical | Retained B6 overview |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_b6_dataset_cellline_overview.csv` | Historical/non-canonical | Retained B6 overview |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_b6_supported_slices.csv` | Historical/non-canonical | Retained B6 overview |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure3_task2_preference_support_summary.csv` | Historical/non-canonical | Retired merged descriptive/preference-support placeholder |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/manifest.json` | Historical/non-canonical | Legacy manifest superseded by `framework_analysis_manifest.json` |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/analysis/figure_analysis_update.md` | Historical/non-canonical | Legacy note for old figure-input backfills |
| `/mnt/NAS_21T/ProjectData/M2M/archive/manuscript_history/notes/manuscript_active_experiment_update.md` | Historical/non-canonical | Legacy note containing A2/B6-era wording |

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
- Support-only diagnostic or exemplar panels may sharpen interpretation, but they do not upgrade the paper to mechanism-causality evidence.
- Legacy scPerturb-only Task2 outputs are historical evidence only.
