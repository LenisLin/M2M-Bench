# Task 1 Retrieval Specification (M2M-Bench)

## Goal

Instance-level cross-source retrieval to quantify modality concordance:

- LINCS → scPerturb
- scPerturb → LINCS

within each treatment modality (`Chemical`, `Genetic`).

## Retrieval label key

- `label_key = cell_std || modality || target_std`

For each query, success is defined by retrieving any gallery sample with the same `label_key`.

## Inputs

- `outputs/task1/data/m1_candidates.parquet`
  - fallback supported: `outputs/task1/data/m1_candidates.csv`
- source tensors:
  - LINCS: `processed/LINCS_Processed/LINCS_Engine1_TrainData.pt`
  - scPerturb: `processed/scPerturb_Processed/scPerturb_Engine1_TrainData_part*.pt`

## Tracks

- `gene`: `y_delta_gene`
- `path`: `y_delta_pathway`

## Metrics

- Per query:
  - `true_rank`
  - `mrr`
  - `success_topk` (k configurable; default 1,5,10)
  - `top1_is_true`
- Group summary (`track`, `direction`, `modality`):
  - mean/median rank
  - MRR
  - Top-k rates
- Null baseline:
  - permutation-style random-rank simulation based on class counts
  - p-values vs null for MRR and Top-k

## Outputs

- `outputs/task1/retrieval/analysis/retrieval_per_query.csv`
- `outputs/task1/retrieval/analysis/retrieval_summary.csv`
- `outputs/task1/retrieval/analysis/retrieval_null_summary.csv`
- `outputs/task1/retrieval/analysis/retrieval_examples.csv`
- `outputs/task1/retrieval/qc/retrieval_input_counts.csv`
- `outputs/task1/retrieval/run_manifest_task1_retrieval.json`
