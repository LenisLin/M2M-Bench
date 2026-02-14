# Task 1 Data Specification (M2M-Bench)

## Goal

Task 1 measures modality difference between LINCS (bulk) and scPerturb (scRNA-seq) under matched context.

## Canonical matching key

- `cell_std`
- `modality` (`Chemical` or `Genetic`, derived from `pert_type`)
- `target_std` (gene label only)
- `source_db` is contrasted as `LINCS` vs `scPerturb`

Protocol fields for matching quality:

- `dose_val`
- `time_val`
- `cond_id`

## Inputs

- `outputs/task0_curated/metadata/unified_meta.parquet`
- raw tensors from processed source:
  - `processed/LINCS_Processed/LINCS_Engine1_TrainData.pt`
  - `processed/scPerturb_Processed/scPerturb_Engine1_TrainData_part*.pt`

## Task 1 generated artifacts

- `outputs/task1/data/task1_index.parquet`
- `outputs/task1/data/m1_candidates.parquet`
- `outputs/task1/data/m1_candidates.csv`
- `outputs/task1/data/m1_matched_pairs.parquet`
- `outputs/task1/data/m1_matched_pairs.csv`
- `outputs/task1/data/m1_tensors/*.npy`

QC tables and plots:

- `outputs/task1/qc/tables/*.csv`
- `outputs/task1/qc/figures/*.png`

Analysis tables:

- `outputs/task1/analysis/modality_gap_per_pair.csv`
- `outputs/task1/analysis/modality_gap_summary.csv`
- `outputs/task1/analysis/modality_gap_by_strata.csv`
- `outputs/task1/analysis/null_distribution_summary.csv`
- `outputs/task1/analysis/significance_results.csv`
- `outputs/task1/analysis/effect_size_table.csv`

Run manifest:

- `outputs/task1/run_manifest_task1.json`
