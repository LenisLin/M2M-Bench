# Task1 Snapshot v1 (Canonical Inputs)

This file enumerates the canonical inputs to freeze for Task1 rebuild.

## scPerturb gene delta (2477) + meta
- CRISPR_delta_gene_recomputed.npy
- CRISPR_delta_meta_recomputed.csv   (alignment key: delta_row_idx)
- Drug_delta_gene_recomputed.npy
- Drug_delta_meta_recomputed.csv     (alignment key: delta_row_idx)

## scPerturb pathway delta (50) + meta
- CRISPR_delta_pathway_recomputed.npy + CRISPR_delta_meta_recomputed.csv (delta_row_idx)
- Drug_delta_pathway_recomputed.npy   + Drug_delta_meta_recomputed.csv   (delta_row_idx)

## scPerturb FM deltas (multiple models; Chemical+Genetic)
Root: data/task1_snapshot_v1/fm_delta/<fm_model>/
For each fm_model in the frozen list:
- Genetic_delta_cell.npy + Genetic_delta_meta.csv  (row_in_modality_matrix)
- Chemical_delta_cell.npy + Chemical_delta_meta.csv (row_in_modality_matrix)

Frozen FM_MODELS:
- geneformer, pca50, pca100, pca200, scbert, scfoundation, scgpt, state, tahoex1_3b, uce

## LINCS bulk gene delta + meta
- LINCS_Engine1_TrainData.pt (720231 x 2477)
- LINCS_Engine1_Meta.csv

## Task1-cross pairing contract (Genetic only)
- cross_pairs_genetic_contract.(csv|parquet)
Required fields:
- global_idx_lincs (indexes LINCS meta and pt rows)
- sc_m2m_row_idx
- sc_delta_row_idx (indexes CRISPR delta arrays)
- expected_cell_line, expected_target
- lincs_cell_line_from_meta, lincs_target_from_meta
- sc_cell_line_from_meta, sc_target_from_meta

## LINCS pathway projection policy (project_on_load)
- data/task1_snapshot_v1/pathway/hallmark-w-2477x50.npy + sha256
- data/task1_snapshot_v1/pathway/lincs-pathway-policy.json
- LINCS_gene_alignment.csv + sha256 (gene order contract)
