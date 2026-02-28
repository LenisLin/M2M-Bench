# Task Definitions

## Representations (delta space)
All comparisons are performed on delta representations:
- delta = treated - control (definition is fixed by the preprocessing contract).

Spaces:
1) Gene space: 2477-dim delta vectors (gene-aligned contract).
2) Pathway space: 50-dim Hallmark projection, computed as:
   pathway_delta = gene_delta @ W   (W: 2477x50, fixed hash + file).
3) Foundation model (FM) delta: embedding deltas (model-specific dims), stored as arrays + meta.

## Task1: Cross-modality / Cross-ecosystem consistency
### Task1-internal (within dataset side)
- scPerturb internal:
  - Genetic (CRISPR) and Chemical (Drug) must both exist.
  - Evaluate in Gene/Pathway spaces and multiple FM spaces.
- LINCS internal:
  - bulk only, no FM.
  - Evaluate in Gene/Pathway spaces.

### Task1-cross (S5; cross dataset)
- Focus: Genetic only (Chemical cross excluded by design because matched contexts are too sparse).
- Cross pairing is defined by an explicit pairing contract file (see docs/snapshots/task1-snapshot-v1.md).
- Task1-cross (S5): Genetic only. Chemical cross is excluded by policy (docs/25).
- Cross alignment contract is single and frozen: global_idx_lincs + sc_delta_row_idx (docs/60).

## Task2: Chemical vs Genetic mechanism consistency (same cell_line+target)
- Not in scope until Task2 inputs are complete.
