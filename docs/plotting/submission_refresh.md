# Submission Refresh Renderer

This renderer creates a publication-oriented replacement set for the legacy PDFs
under `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready/Figures_Submission`.

## Source Data

- Root: `/mnt/NAS_21T/ProjectData/M2M/R_Vis_Ready`
- Figure 1 refresh reads the data-description tables plus Task 1 pairwise and
  retrieval summaries.
- Figure 2 refresh reads the frozen effect summary, Task 2 tracer summaries,
  Task 2 target-enrichment tables, and protocol-correlation tables.
- Figure 3 refresh reads the Task 3 target-gain, scoreboard, and comparison
  tables.

## Rendered Outputs

Running the refresh pipeline writes vector PDFs to:

- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/submission_refresh/figure1_submission_refresh.pdf`
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/submission_refresh/figure2_submission_refresh.pdf`
- `/mnt/NAS_21T/ProjectData/M2M/runs/manuscript_support/submission_refresh/figure3_submission_refresh.pdf`

## Design Changes

- The radial sunburst is replaced with ordered bar-based structure views.
- The bubble summary is replaced by a dot-heatmap that separates effect size
  from significance more cleanly.
- Dense tracer and protocol figures are reframed as ribbons, dual-tail
  enrichment lollipops, and ranked signed point plots.
- Task 3 summary figures are rebuilt as heatmap and dot-matrix views with
  stable semantic colors and no on-canvas `Figure X` title blocks.

## Run

From the repo root:

```bash
conda run -n Spatial Rscript plotting/R/render_submission_refresh.R
```

Or from any directory:

```bash
conda run -n Spatial Rscript /home/lenislin/Experiment/projects/M2M/plotting/R/render_submission_refresh.R
```

Optional environment variables:

- `M2M_SUBMISSION_REFRESH_DATA_ROOT`
- `M2M_SUBMISSION_REFRESH_OUTPUT_ROOT`
