# M2M-Bench Redesign Checkpoint

This file defines the current benchmark system and the terminology shared across
repo docs.

## Benchmark Position

- M2M-Bench is a benchmark/evaluation paper.
- Figure 1 defines the benchmark.
- Figure 2 carries Task1 main evidence.
- Figure 3 carries Task2 main evidence.
- `Gene` and `Pathway` are the benchmark-wide representation spaces.
- `FM` enters the main manuscript only through the `scPerturb/K562`
  `Figure 3F` local-only panel.
- Task1 and Task2 stay separate in the main text.

## Shared Terminology

- `Task1`: modality concordance with `perturbation_type` held fixed.
- `Task2`: mechanism concordance between chemical and genetic cohorts inside one
  dataset.
- `anchor_gene`: the Task2 unit gene in `(dataset, cell_line, anchor_gene)`.
- `perturbation_gene`: the benchmark instance identity field. For genetic rows
  it is one perturbed gene. For chemical rows it is one canonical target-set
  string.
- `query_instance_id`: the per-query identifier carried by retrieval tables.
- `C2G`: chemical-to-genetic retrieval in Task2.
- `G2C`: genetic-to-chemical retrieval in Task2.
- `pair_mean_enrichment`: the signed summary used to rank entities in the
  Figure `2E`, `2F`, `3C`, and `3D` pattern panels.
- `active data roots`: `/mnt/NAS_21T/ProjectData/M2M/data/task1` and
  `/mnt/NAS_21T/ProjectData/M2M/data/task2`.
- `active run roots`: `/mnt/NAS_21T/ProjectData/M2M/runs`,
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`, and
  `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`.

## Instance Model

- A benchmark instance is a delta-space object, not a raw row.
- Every perturbation row maps to one aggregated delta instance.
- The instance fields are `dataset`, `cell_line`, `perturbation_type`,
  `perturbation_gene`, `time`, and `dose`.
- `perturbation_type` takes `chemical` or `genetic`.
- For chemical instances, `perturbation_gene` uses uppercase tokens, duplicate
  removal, stable ordering, and `|` as the delimiter.
- Multi-target chemicals stay as one instance row. Membership expansion for
  Task2 happens at the cohort-construction layer.

## Analysis Families And Metrics

- `Group` analysis is directionless.
- `Instance retrieval` is directional.
- Group metrics are `PCC`, `cosine similarity`, and bias-corrected
  `e_distance`.
- `e_distance` uses cell-wise squared Euclidean distances.
- Retrieval metrics are corrected `Hit@1`, `Hit@3`, `Hit@5`, and corrected
  `MRR`.
- Raw retrieval metrics may be stored with the same tables.
- Chance correction uses uniform random ranking over the lawful gallery
  conditioned on `gallery_size` and `n_positive_keys`.

## Task1 Contract

- Task1 internal unit:
  `(dataset, cell_line, perturbation_type, perturbation_gene)`.
- Task1 cross unit:
  `(cell_line, perturbation_type, perturbation_gene)` with both datasets
  present.
- Chemical `time` and `dose` stay as metadata and do not enter the Task1 unit.
- Multi-target chemicals use exact target-set equality in Task1.
- Task1 retrieval is single-positive.
- Task1 internal retrieval uses leave-one-out true centroids.
- `n_positive_keys` is fixed to `1` for every lawful Task1 query.

## Task2 Contract

- Task2 unit:
  `(dataset, cell_line, anchor_gene)`.
- A lawful Task2 unit contains at least one chemical member and at least one
  genetic member.
- A chemical row joins a Task2 unit when `anchor_gene` is a member of the
  chemical target set.
- `perturbation_gene` remains the row identity field in Task2 result tables.
- Aggregated Task2 tables are keyed by `anchor_gene`.
- Task2 directions are `C2G` and `G2C`.
- `G2C` uses one positive chemical centroid per lawful Task2 unit.
- `C2G` uses chemical instance queries keyed by `query_instance_id`.
- Task2 core metrics are computed within each `dataset` and `cell_line`.

## Stage Map

- `S0`: data inventory
- `S1`: Task1 internal metrics
- `S2`: Task1 cross metrics
- `S3`: Task2 multisource cohort build
- `S4`: Task2 group concordance
- `S5`: Task2 retrieval
- `S6`: Task2 synthesis
- `S7`: project synthesis assembled from audited Task1 and Task2 outputs

## Figure Contract

### Figure 2

- `2A`: Task1 lawful-scope composition
- `2B`: Task1 shared matched-unit scoreboard
- `2C`: Task1 internal-to-cross degradation
- `2D`: Task1 paired `Gene` versus `Pathway` comparison
- `2E`: Task1 cell-line pattern ranked by `pair_mean_enrichment`
- `2F`: Task1 `perturbation_gene` pattern ranked by `pair_mean_enrichment`

### Figure 3

- `3A`: Task2 lawful-scope composition
- `3B`: Task2 performance backbone
- `3C`: Task2 cell-line pattern ranked by `pair_mean_enrichment`
- `3D`: Task2 `anchor_gene` pattern ranked by `pair_mean_enrichment`
- `3E`: Task2 `C2G` paired `Gene` versus `Pathway` comparison
- `3F`: `scPerturb/K562` `FM` local-only absolute-performance panel

## Rendering Thresholds

- Figure `2D` and Figure `3E` render surfaces with `n_pairs >= 3`.
- Figure `2E`, `2F`, `3C`, and `3D` render entities with `support_n >= 3`.
- Figure `3F` uses the approved `scPerturb/K562` scope only.
