# Task1 Contract

## Scope

Task1 evaluates modality concordance while holding `perturbation_type` fixed.

### Active Task1 Slices

- `LINCS` internal chemical
- `LINCS` internal genetic
- `scPerturb` internal chemical
- `scPerturb` internal genetic
- matched `LINCS` to `scPerturb` genetic cross slice

### Active Representation Sets

- `LINCS` internal: `Gene`, `Pathway`
- `scPerturb` internal: `Gene`, `Pathway`, `scgpt`, `geneformer`, `scbert`,
  `scfoundation`, `uce`, `state`, `tahoe-x1`
- Task1 cross: `Gene`, `Pathway`

## Unit And Identity

- Task1 internal unit:
  `(dataset, cell_line, perturbation_type, perturbation_gene)`
- Task1 cross unit:
  `(cell_line, perturbation_type, perturbation_gene)`
- `perturbation_gene` is the Task1 unit identity field.
- Chemical `time` and `dose` stay as metadata and do not enter the Task1 unit.
- A multi-target chemical uses the canonical `perturbation_gene` target-set
  string from `docs/redesign_checkpoint.md`.

## Group Analysis

- Internal group analysis uses deterministic split-half cohorts.
- Cross group analysis uses matched units shared by `LINCS` and `scPerturb`.
- Group metrics are `cosine similarity`, `PCC`, and bias-corrected
  `e_distance`.
- `e_distance` uses cell-wise squared Euclidean distances.
- Tables mark `underpowered_for_e_distance=true` when a comparison side has
  fewer than `2` instances.

## Retrieval

- Retrieval queries are single instances.
- Retrieval galleries are centroids built from lawful Task1 units.
- Task1 retrieval is single-positive.
- Task1 internal retrieval uses leave-one-out true centroids.
- `n_positive_keys = 1` for every lawful Task1 query.
- Task1 cross alignment uses `global_idx_lincs` on the `LINCS` side and
  `sc_delta_row_idx` on the `scPerturb` side.
- Primary retrieval metrics are corrected `Hit@1`, `Hit@3`, `Hit@5`, and
  corrected `MRR`.

## Active Output Surface

- `task1_group_concordance_long.csv`
- `task1_retrieval_per_query.parquet`
- `task1_retrieval_summary.csv`
- `task1_leaderboard_long.csv`
- `task1_cross_alignment_proof.csv`

See `docs/contracts/output-schemas.md` for field-level table contracts.
