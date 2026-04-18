# Project Positioning

## Paper Type

M2M-Bench is a benchmark/evaluation paper for transcriptome-centered
perturbation-response concordance.

## Central Question

The benchmark asks how consistently perturbation signals align across:

- modality concordance with `perturbation_type` held fixed in `Task1`
- mechanism concordance between chemical and genetic cohorts inside one dataset
  in `Task2`

## Task Split

- `Task1` measures modality concordance with `perturbation_type` held fixed.
- `Task2` measures mechanism concordance between chemical and genetic cohorts
  within one dataset.
- The manuscript keeps these two task families separate.

## Representation Split

- `Gene` and `Pathway` define the benchmark-wide representation comparison.
- `FM` appears in the main manuscript only through the `scPerturb/K562`
  `Figure 3F` local-only panel.

## Figure Split

- Figure 1: benchmark definition
- Figure 2: Task1 main evidence
- Figure 3: Task2 main evidence

## Language Rules

- Use `anchor_gene` for the Task2 unit field.
- Use `perturbation_gene` for row identity.
- Use `query_instance_id` for retrieval query rows.
- Use `pair_mean_enrichment` for ranked pattern summaries.
- Keep `C2G` as the main Task2 retrieval direction in Figure 3 and keep `G2C`
  as the second Task2 direction.
