# Manuscript Master

## Paper Shape

- M2M-Bench is a benchmark/evaluation paper.
- Figure 1 defines the benchmark.
- Figure 2 presents Task1 main evidence.
- Figure 3 presents Task2 main evidence.
- `FM` enters the main manuscript only through the `scPerturb/K562`
  `Figure 3F` local-only panel.
- Task1 and Task2 stay separate in the main text.

## Core Claims

1. M2M-Bench defines an audited benchmark for perturbation-response concordance.
2. Task1 measures modality concordance with `perturbation_type` held fixed
   across internal and cross settings.
3. Task2 measures mechanism concordance between chemical and genetic cohorts
   inside one dataset, with core metrics reported within each `dataset` and
   `cell_line`.

## Figure 1

Figure 1 defines the benchmark question, workflow, and lawful comparison scope.

## Figure 2

- `2A`: Task1 lawful-scope composition
- `2B`: Task1 shared matched-unit scoreboard
- `2C`: Task1 internal-to-cross degradation
- `2D`: Task1 paired `Gene` versus `Pathway` comparison
- `2E`: Task1 cell-line pattern ranked by `pair_mean_enrichment`
- `2F`: Task1 `perturbation_gene` pattern ranked by `pair_mean_enrichment`

Figure 2 uses:

- `n_pairs >= 3` for `2D`
- `support_n >= 3` for `2E` and `2F`

## Figure 3

- `3A`: Task2 lawful-scope composition
- `3B`: Task2 performance backbone
- `3C`: Task2 cell-line pattern ranked by `pair_mean_enrichment`
- `3D`: Task2 `anchor_gene` pattern ranked by `pair_mean_enrichment`
- `3E`: Task2 `C2G` paired `Gene` versus `Pathway` comparison
- `3F`: `scPerturb/K562` `FM` local-only absolute-performance panel

Figure 3 uses:

- `C2G` as the first Task2 retrieval direction in the main text
- `G2C` as the second Task2 retrieval direction
- `n_pairs >= 3` for `3E`
- `support_n >= 3` for `3C` and `3D`
- the approved `scPerturb/K562` scope only for `3F`

## Active Manuscript Analysis Root

`/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`

## Wording Rules

- use `anchor_gene` for Task2 unit identity
- use `perturbation_gene` for row identity
- use `query_instance_id` for retrieval query rows
- use `pair_mean_enrichment` for ranked pattern summaries
