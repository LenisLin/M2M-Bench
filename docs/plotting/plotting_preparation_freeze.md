# Plotting Preparation Freeze

## Scope

- Python assembles plot-ready tables.
- R renders the final panels.
- Plotting does not change benchmark semantics.

## Active Roots

- Manuscript analysis:
  `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_active/analysis`
- Plot review export:
  `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current`

## Main Panel Set

### Figure 1

- benchmark question
- workflow
- lawful comparison scope

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

## Assembly Rules

- Python owns pairing, thresholding, ranking, significance joins, and panel
  membership.
- R owns factor ordering, labels, themes, panel composition, and vector export.
- R does not recompute panel-defining summaries.

## Threshold Rules

- `2D` and `3E` render surfaces with `n_pairs >= 3`.
- `2E`, `2F`, `3C`, and `3D` render entities with `support_n >= 3`.
- `3F` uses the approved `scPerturb/K562` scope only.
