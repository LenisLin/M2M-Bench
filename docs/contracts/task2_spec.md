# Task2 Contract

## Scope

Task2 evaluates mechanism concordance between chemical and genetic cohorts
inside one dataset.

### Active Datasets

- `LINCS`
- `scPerturb`

### Active Representation Sets

- `LINCS`: `Gene`, `Pathway`
- `scPerturb/K562`: `Gene`, `Pathway`, `scgpt`, `geneformer`, `scbert`,
  `scfoundation`, `uce`, `state`, `tahoe-x1`

## Unit And Membership

- Task2 unit:
  `(dataset, cell_line, anchor_gene)`
- `anchor_gene` is the Task2 unit field.
- A lawful Task2 unit contains at least one chemical member and at least one
  genetic member.
- A chemical instance joins a Task2 unit when `anchor_gene` belongs to the
  chemical target set encoded in `perturbation_gene`.
- `perturbation_gene` remains the row identity field.
- `query_instance_id` is the retrieval query identifier.

## Data Surface

- Active Task2 data root: `/mnt/NAS_21T/ProjectData/M2M/data/task2`
- Active Task2 cohort and membership outputs:
  - `task2_pairs_coverage.csv`
  - `task2_row_membership.parquet`
  - `delta_meta.csv`
  - `representation_availability_registry.csv`
  - `snapshot_manifest.json`

## Group Analysis

- Task2 group analysis compares chemical and genetic cohorts within each
  `(dataset, cell_line, anchor_gene)` unit.
- Core group metrics are `cosine similarity`, `PCC`, and bias-corrected
  `e_distance`.
- Task2 core group metrics are computed within each `dataset` and `cell_line`.

## Retrieval

- `C2G`: query = one chemical instance keyed by `query_instance_id`; positives =
  lawful genetic centroids keyed by `anchor_gene`
- `G2C`: query = one genetic instance keyed by `query_instance_id`; positive =
  the chemical centroid for the same `(dataset, cell_line, anchor_gene)` unit
- `C2G` and `G2C` stay as separate directions in every Task2 output table.
- Primary retrieval metrics are corrected `Hit@1`, `Hit@3`, `Hit@5`, and
  corrected `MRR`.

## Synthesis Rule

Task2 synthesis combines audited group and retrieval summaries without changing
their unit definitions or denominator fields.

## Active Output Surface

- `task2_pairs_coverage.csv`
- `task2_group_concordance_long.csv`
- `task2_group_leaderboard.csv`
- `task2_retrieval_per_query.parquet`
- `task2_retrieval_summary_long.csv`
- `task2_retrieval_leaderboard.csv`
- `task2_benchmark_summary_long.csv`

See `docs/contracts/output-schemas.md` for field-level table contracts.
