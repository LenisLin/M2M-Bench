# M2M-Bench Analysis Workflow (Readable Script Style)

This file records the recommended runnable order for the current refactored analysis scripts.

---

## 1) Task0 curation (recommended entrypoint)

Purpose: generate canonical Task0 outputs used by all downstream tasks.

```bash
python scripts/task0_pipeline.py \
  --run-curation \
  --no-run-attrition-audit \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/task0_curated
```

If parquet engines are unavailable, run:

```bash
python scripts/task0_curate_data.py \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/task0_curated \
  --no-save-parquet \
  --save-unified-csv
```

Main outputs:
- `outputs/task0_curated/metadata/unified_meta.parquet`
- `outputs/task0_curated/metadata/unified_meta.csv` (parquet-free fallback)
- `outputs/task0_curated/metadata/level1_pairs_meta.parquet`
- `outputs/task0_curated/tensors/Level1_Tensors_Split/*.pt`
- `outputs/task0_curated/run_manifest_task0.json`

---

## 2) Task1 full pipeline (recommended entrypoint)

Purpose: one command for group-wise gap + retrieval + confounder + set-level outputs.

```bash
python scripts/task1_pipeline.py \
  --task0-unified-meta-path ./outputs/task0_curated/metadata/unified_meta.parquet \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --task1-output-dir ./outputs/task1 \
  --retrieval-output-dir ./outputs/task1/retrieval \
  --confounder-output-dir ./outputs/task1/confounder \
  --set-level-output-dir ./outputs/task1/retrieval_set_level \
  --strict-enable \
  --retrieval-run-balanced-eval
```

Main outputs:
- `outputs/task1/analysis/modality_gap_per_pair.csv`
- `outputs/task1/retrieval/analysis/retrieval_summary.csv`
- `outputs/task1/confounder/analysis/factor_permutation_tests.csv`
- `outputs/task1/confounder/analysis/context_reliability_labels.csv`
- `outputs/task1/retrieval_set_level/analysis/set_level_context_metrics.csv`

---

## 3) Task0 → Task1 attrition audit

Purpose: quantify data drop from Task0 unified pool to Task1 overlap/matching.

```bash
python scripts/task01_attrition_audit.py \
  --bundle-path ./outputs/task0_curated/bundle/m2m_task0_bundle.pt \
  --matched-pairs-path ./outputs/task1/data/m1_matched_pairs.csv \
  --output-dir ./outputs/task1_audit
```

Main outputs:
- `outputs/task1_audit/analysis/stage_summary.csv`
- `outputs/task1_audit/analysis/context_drop_reasons_summary.csv`
- `outputs/task1_audit/analysis/stage3_drop_reason_summary.csv`
- `outputs/task1_audit/analysis/protocol_gap_summary.csv`

---

## 4) Task1 set-level retrieval sanity (already generated in this repo)

Purpose: avoid only query-level interpretation and check context-level retrieval.

Main outputs:
- `outputs/task1/retrieval_set_level/centroid_retrieval_summary.csv`
- `outputs/task1/retrieval_set_level/set_level_summary_from_query_aggregation.csv`

---

## 5) Task2 mechanism analysis (No DomainType main contrast)

Purpose: run mechanism analysis with readable script flow and remove Same/Cross domain primary comparison.
Current recommendation is `same_only` domain scope, so analysis focuses on mechanism difference under same-domain data.

```bash
python scripts/task2_mechanism_analysis.py \
  --pairwise-wide-path /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task1_Metrics/Task1_Pairwise_Metrics_Wide.csv \
  --retrieval-per-query-path /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task1_Retrieval_MultiScenario/Task1_Retrieval_MultiScenario_PerQuery.csv \
  --domain-scope same_only \
  --output-dir ./outputs/task2_nodomain
```

Main outputs:
- `outputs/task2_nodomain/analysis/Step1_L1_Instance_Tidy_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step1_L1_Context_Aggregated_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step1_Tests_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step2_Context_Labels_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step4_CaseStudy_Tracer_NoDomain.csv`
- `outputs/task2_nodomain/analysis/Step5_*.csv`

---

## 6) Task2 legacy-output audit (quality check)

Purpose: verify existing old Task2 outputs are internally consistent.

```bash
python scripts/task2_results_audit.py \
  --task2-dir /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task2_Unified \
  --output-dir ./outputs/task2_audit
```

Main outputs:
- `outputs/task2_audit/analysis/task2_consistency_checks.csv`
- `outputs/task2_audit/analysis/task2_primary_slice_summary.csv`
- `outputs/task2_audit/analysis/task2_audit_report.md`

---

## 7) Task3 results audit (recommended)

Purpose: verify completed Task3 outputs are internally consistent before interpretation.

```bash
python scripts/task3_results_audit.py \
  --deltas-dir /mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Deltas_K562_CellLevel \
  --analysis-dir /mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Analysis \
  --unified-dir /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified \
  --output-dir ./outputs/task3_audit
```

Main outputs:
- `outputs/task3_audit/analysis/task3_consistency_checks.csv`
- `outputs/task3_audit/analysis/task3_key_summary.csv`
- `outputs/task3_audit/analysis/task3_track_view_direction_summary.csv`
- `outputs/task3_audit/analysis/task3_audit_report.md`

---

## 8) Task3 full pipeline wrapper (optional rerun)

Purpose: rerun Script 9/10/11 with one readable entrypoint and optional audit.

```bash
python scripts/task3_pipeline.py \
  --run-prepare-deltas \
  --run-analysis \
  --run-collect-unified \
  --run-audit
```

For dry-run preview:
```bash
python scripts/task3_pipeline.py --dry-run
```

---

## 9) Task1 reviewer add-ons (must-run for manuscript robustness)

Purpose: generate reviewer-critical robustness outputs without rerunning heavy Task1 core pipeline.

```bash
python scripts/task1_reviewer_addons.py \
  --m1-candidates-path ./outputs/task1/data/m1_candidates.csv \
  --per-pair-path ./outputs/task1/analysis/modality_gap_per_pair.csv \
  --retrieval-summary-path ./outputs/task1/retrieval/analysis/retrieval_summary.csv \
  --retrieval-null-path ./outputs/task1/retrieval/analysis/retrieval_null_summary.csv \
  --retrieval-per-query-path ./outputs/task1/retrieval/analysis/retrieval_per_query.csv \
  --context-overlap-path ./outputs/task1_audit/analysis/context_overlap_counts.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/task1_reviewer_fixes
```

Main outputs:
- `outputs/task1_reviewer_fixes/analysis/retrieval_dual_report.csv`
- `outputs/task1_reviewer_fixes/analysis/effect_calibration_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/protocol_continuous_spearman.csv`
- `outputs/task1_reviewer_fixes/analysis/protocol_continuous_partial_spearman.csv`
- `outputs/task1_reviewer_fixes/analysis/leave_one_cell_out_pairwise.csv`
- `outputs/task1_reviewer_fixes/analysis/benchmarkability_zone_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/strict_subset_composition.csv`
- `outputs/task1_reviewer_fixes/analysis/lincs_internal_consistency_summary.csv`
- `outputs/task1_reviewer_fixes/analysis/task1_context_definition_comparison.csv`
- `outputs/task1_reviewer_fixes/analysis/task1_metacontext_map.csv`

---

## 10) Task2 target-tier stratification (recommended)

Purpose: add target-confidence / polypharmacology tier analysis on top of Task2 context labels.

```bash
python scripts/task2_target_tier_analysis.py \
  --context-labels-path ./outputs/task2_nodomain/analysis/Step2_Context_Labels_NoDomain.csv \
  --output-dir ./outputs/task2_nodomain/target_tier
```

Main outputs:
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_summary.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_enrichment.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_map_used.csv`
- `outputs/task2_nodomain/target_tier/analysis/task2_target_tier_mapping_qc.csv`

---

## 11) Task3 model meta-analysis (recommended)

Purpose: relate Task3 performance to model attributes (`size_class`, `perturbation_trained`).

```bash
python scripts/task3_fm_meta_analysis.py \
  --scoreboard-long-path /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified/viz_scoreboard_long.csv \
  --output-dir ./outputs/task3_meta
```

Main outputs:
- `outputs/task3_meta/analysis/task3_model_scoreboard_meta.csv`
- `outputs/task3_meta/analysis/task3_meta_perturbation_training_tests.csv`
- `outputs/task3_meta/analysis/task3_meta_size_trend_tests.csv`
- `outputs/task3_meta/analysis/task3_mean_scaled_best_definition.csv`

---

## 12) Reproducibility pack (recommended before manuscript freeze)

Purpose: gather seeds, run manifests, environment, and one-command rerun paths.

```bash
python scripts/build_reproducibility_pack.py \
  --output-dir ./outputs/reproducibility
```

Main outputs:
- `outputs/reproducibility/analysis/repro_manifest_inventory.csv`
- `outputs/reproducibility/analysis/repro_seed_registry.csv`
- `outputs/reproducibility/analysis/repro_environment.json`
- `outputs/reproducibility/analysis/repro_one_command_paths.md`

---

## 13) Task1 retrieval sensitivity (recommended for reviewer robustness)

Purpose: test balanced retrieval stability across gallery sizes and seeds.

```bash
python scripts/task1_retrieval_sensitivity.py \
  --m1-candidates-path ./outputs/task1/data/m1_candidates.csv \
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \
  --output-dir ./outputs/task1_reviewer_fixes/retrieval_sensitivity \
  --gallery-sizes 64,128,256,512 \
  --seeds 42,43
```

Main outputs:
- `outputs/task1_reviewer_fixes/retrieval_sensitivity/analysis/retrieval_sensitivity_long.csv`
- `outputs/task1_reviewer_fixes/retrieval_sensitivity/analysis/retrieval_sensitivity_summary.csv`
- `outputs/task1_reviewer_fixes/retrieval_sensitivity/analysis/retrieval_sensitivity_modality_compact.csv`

---

## 14) Project alignment audit (recommended after major iteration)

Purpose: ensure manuscript-facing core claims remain aligned with generated outputs.

```bash
python scripts/project_alignment_audit.py \
  --output-dir ./outputs/project_alignment_audit
```

Main outputs:
- `outputs/project_alignment_audit/analysis/project_alignment_checks.csv`
- `outputs/project_alignment_audit/analysis/project_alignment_summary.csv`
- `outputs/project_alignment_audit/analysis/project_alignment_report.md`

---

## Notes

- Refactored scripts use explicit stage blocks and CSV-first outputs for easier debugging.
- Domain/source fields are retained for diagnostics, but Task2 main hypothesis testing is domain-agnostic.
- Legacy step-only entrypoints are still available in `scripts/task0_curate_data.py` and `scripts/task1_*.py`.
