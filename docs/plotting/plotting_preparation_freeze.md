# Plotting Preparation and Rendering Freeze

This file is the implementation-facing freeze note for manuscript plotting and
final R rendering. It operationalizes the accepted figure plan without creating
new canonical benchmark evidence.

## Scope

- Figure 1 is doc/snapshot-derived and consumes plot-ready registries built by
  `scripts/manuscript_plot_ready_tables.py`; the current main-text figure is a
  reduced `1B` composition panel only.
- Figure 1 workflow and data-filtering schematics are deferred as separate
  markdown-described assets and are not forced into the current main render.
- Figure 2 and Figure 3 main panels consume panel-specific `plot_ready/`
  tables derived from the frozen canonical `analysis/` objects.
- R plotting code lives under `plotting/R/`; main-figure scripts read only
  `plot_ready/` plus config.
- Python owns any new pairing, binning, display-set selection, and
  significance-join logic required by the redesigned panels.
- The active review/export contract for this revision writes only to
  `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current/{analysis,plot_ready,figures}`;
  historical `/mnt/.../manuscript_support` outputs remain untouched during review export.
- This layer does not change benchmark logic, regenerate analysis outputs, or
  rewrite manuscript prose.

## Final figure layout

### Figure 1

- Main: `1B` compact standalone sunburst panels (`1B_LINCS`, `1B_scPerturb`)
  with hierarchy `Task -> Perturbation type -> Cell-line / target context`
- Deferred: `1A` workflow schematic prompt and data-filtering schematic prompt
- Supplementary / methods support: `1C` lawful scope matrix, `1D`
  representation/modifier availability
- Deleted from main: `1E`

### Figure 2

- `2A`: stacked Task1 scope bars
- `2B`: scoreboard-style internal benchmark matrix
- `2C`: unified matched-unit Gene-vs-Pathway comparison across internal and cross contexts
- `2D`: family-separated (`Group`, `Retrieval`) internal-to-cross paired summary comparison
- `2E`: four-block internal-only cell-line enrichment suitability panel
- `2F`: four-block internal-only target enrichment suitability panel

### Figure 3

- `3B`: compact paired Task2 benchmark comparison with dataset blocks and cell-line rows
- `3C`: LINCS-only ranked cell-line suitability panel
- `3D`: dataset-faceted ranked target suitability panel (`LINCS` and `scPerturb`)
- `3E`: performance-only matched-target Gene-vs-Pathway comparison, rendered as dataset x metric small panels
- `3F`: FM-local absolute-performance panels (backed by `figure3_panel_3f_fm_local_tradeoff.csv`)
- Removed from main: `3A` Task2 scope bars, retained only as support/reference data
- Removed from main: `3J` Task1 contextual support, retained only as
  supplementary/reference support

### Extended figures

- `EF1`: Figure 1 contract/support registry
- `EF2`: Task1 internal full decomposition
- `EF3`: Task1 cross plus degradation full decomposition
- `EF4`: Task1 dataset x perturbation-type enrichment decomposition
- `EF5`: Task2 full performance decomposition
- `EF6`: ranked Task2 suitability decomposition
- `EF7`: Task2 modifier analyses
- `EF8`: corrected Task2 representation plus contextual support decomposition
- `EF9`: support-vs-suitability diagnostics plus target exemplar archetypes

## Plot-ready outputs

`scripts/run_spatial_figure_revision.sh` materializes the review-facing
downstream plotting layer under
`/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current/{analysis,plot_ready,figures}`.

### Figure 1 registries

- `figure1/figure1_panel_1b_dataset_context_coverage.csv`
- `figure1/figure1_panel_1c_lawful_scope_matrix.csv`
- `figure1/figure1_panel_1d_representation_modifier_availability.csv`
- `figure1/figure1_panel_1e_result_object_map.csv`
- `extended/extended_figure1_support_registry.csv`

### Figure 2 / Figure 3 panel-ready reductions

- `figure2/figure2_panel_2a_task1_scope.csv`
- `figure2/figure2_panel_2b_internal_performance_overview.csv`
- `figure2/figure2_panel_2c_gene_vs_pathway_matched_units.csv`
- `figure2/figure2_panel_2d_internal_to_cross_degradation.csv`
- `figure2/figure2_panel_2e_cell_line_high_concordance_summary.csv`
- `figure2/figure2_panel_2f_target_high_concordance_summary.csv`
- `extended/extended_figure4_cell_line_high_concordance_full.csv`
- `extended/extended_figure4_target_high_concordance_full.csv`
- `figure3/figure3_panel_3b_c2g_performance_overview.csv`
- `figure3/figure3_panel_3c_cell_line_pattern.csv`
- `figure3/figure3_panel_3d_target_pattern_summary.csv`
- `figure3/figure3_panel_3e_gene_vs_pathway_paired.csv` (backs main `3E`)
- `extended/extended_figure6_target_pattern_full.csv`
- `figure3/figure3_panel_3f_fm_local_tradeoff.csv` (backs main `3F`)
- `figure3/figure3_panel_3j_task1_contextual_support_reference.csv`
- `extended/extended_figure8_contextual_support_full.csv`
- `extended/extended_figure9_support_vs_suitability.csv`
- `extended/extended_figure9_target_exemplars.csv`

## Frozen implementation rules

- `2A` must preserve an explicit `ALL` cross-eligibility bar rather than
  coercing those rows into `FM`.
- `2B` and `3B` are scoreboard matrices: Python defines the cell membership
  and row membership, but `3B` now renders paired Gene/Pathway comparison
  marks rather than fill-score tiles.
- `2C` must be rebuilt from Task1 bridge-summary long rows plus
  `manuscript_comparison_statistics.csv`, not from the old cross-overview
  shortcut, and it replaces the old `2G` split.
- `1B` must use fixed dataset-level sunburst staging with deterministic
  ordering, small-slice folding to `Others`, fixed hierarchy
  `Task -> Perturbation type -> Context`, and direct on-arc labels instead of
  floating repel labels.
- `2D` must use shared-unit Task1 bridge detail as a direct Internal-vs-Cross
  paired summary surface, separate `Group` and `Retrieval` into distinct family
  blocks, and render any eDist-family metric from the logged transform.
  Significance annotations come from the frozen comparison-statistics layer
  rather than on-canvas recomputation.
- `2E` and `2F` must stage internal-only enrichment by
  `dataset x perturbation_type x entity x representation`, then apply dual-tail
  selection (`3` high + `3` low) independently within each dataset x
  perturbation-type block. Plot-ready tables must expose `row_block`,
  `col_block`, `support_n`, `shared_flag`, and `selection_tail`.
- `3B` must use dataset blocks, one row per cell line, and one paired
  Gene/Pathway comparison unit per metric; it no longer uses the legacy
  scoreboard-lane grammar.
- `3C` and `3D` are staged as ranked suitability panels derived from the
  pattern summaries. `3C` is LINCS-only in the main figure; `3D` retains both
  datasets and highlights shared targets across LINCS and scPerturb.
- `3C` and `3D` must carry `row_rank`, `suitability_score`, and right-side
  support labels in plot-ready staging; they no longer use dual-tail exemplar
  selection in the main figure.
- `3E` must use the Task2 target-level performance surface only, with paired
  units keyed by `(dataset, cell_line, target, metric_name, direction)`.
  Both LINCS and scPerturb therefore compare Gene vs Pathway on matched
  `(cell_line, target)` units rather than cell-line-only aggregates, and the
  frozen significance layer comes from
  `figure3_task2_performance_target_level_representation`.
- The Task2 scope summary remains available, but `3A` is removed from the
  live main-panel export contract.
- `3F` remains strictly `scPerturb/K562` local and uses metric-pair trade-off
  scatter panels, with builder-side scope assertions on `dataset` and
  `cell_line`.
- `3J` is supplementary/reference only and does not appear in the main figure.
- `EF9A` is diagnostic-only and must not be described as a new inferential family.
- `EF9B` is exemplar-only and must preserve a deterministic selection rule over the frozen target-pattern surface.
- `plot_ready/` remains a support-only downstream layer and not a new
  canonical evidence layer.

## Final rendering contract

- `plotting/R/render_spatial_revision_panels.R` is the sole review-export entry
  point for this revision cycle.
- Review export is panel-only: every main and supplementary panel is emitted as
  its own editable vector PDF.
- The export writes a panel manifest to
  `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current/figures/panel_export_manifest.tsv`, and
  QC/re-render state is tracked there rather than in handwritten notes.
- Final manuscript export is editable vector PDF only.
- Panel widths may exceed the legacy `178 mm` default for standalone review
  panels, but each panel still has an explicit width/height contract in the R
  panel-dimension registry.
- No panel titles, subtitles, panel footers, or figure-level explanatory text
  blocks are rendered inside the figure canvas.
- Panel letters, axes, legends, and data marks remain on-canvas.
- Figure legends are maintained in
  `docs/plotting/manuscript_figure_legends.md`, not in the figure canvas.
- Main Figure 2 and Figure 3 scripts consume `plot_ready/` objects and must
  not perform panel-defining large-table summarization inside ggplot code.
- Final figure exports target downstream manual refinement in Adobe; raster
  layers are disallowed on the final PDF path unless separately justified.
- R rendering must run in conda env `Spatial`; missing packages must be
  installed into that environment rather than routed through repo-local
  libraries or cache relocations.

## Plot-data / R boundary

Python plot-ready builders own:

- main-panel scope normalization and panel-specific extract assembly
- underpowered filtering
- pairing, confidence-interval, and significance joins
- binning and display-set selection
- rank and percentile calculations that define panel membership or order
- FM-local hybrid target/query joins
- repeated `perturbation_type` collapse for contextual support
- Figure 1 doc-derived registries

R rendering owns:

- CSV loading from `plot_ready/`
- factor ordering and label formatting
- scales, themes, legends, and panel composition
- editable vector-PDF export
- lightweight column selection / renaming / formatting adjustments only

Forbidden in R:

- on-the-fly summarization of raw large support tables
- reapplication of underpowered filtering
- recalculation of paired statistics or confidence intervals
- recomputation of plot-defining hybrid joins or rank membership
- repo-local library injection or repo-local cache redirection

## R plotting namespace

- Shared helpers: `plotting/R/shared/`
- Figure scripts: `plotting/R/figures/`
- Config files: `plotting/R/config/`
- Export helper: `plotting/R/shared/export.R`
- Figure dimension policy: `plotting/R/config/figure_dimensions.csv`
- Panel-letter registry: `plotting/R/config/panel_labels.csv`
- Legend prose: `docs/plotting/manuscript_figure_legends.md`
- Frozen palette config:
  - `plotting/R/config/representation_palette.csv`
  - `plotting/R/config/palette_values.csv`

## Palette freeze

### Section-band / container palette

- task-design blue fill: `#D9ECF7`
- task-design blue border: `#6EA6D7`
- datasets lavender fill: `#ECE8F3`
- datasets lavender border: `#9C7AC7`
- methods gray fill: `#E6E7E8`
- methods gray border: `#8F8F8F`
- evaluation rose fill: `#F4E4E8`
- evaluation rose border: `#D47F92`

### Auxiliary count palette

- context-count blue: `#5E88C1`
- perturbation-count magenta: `#C76898`
- light perturbation pink: `#F1A6B5`

### Representation palette

- Gene: `#5E88C1`
- Pathway: `#C76898`
- FM: `#9C7AC7`

### Direction palette

- C2G: `#C76898`
- G2C: `#8F8F8F`

### Perturbation palette

- drug / chemical: `#C76898`
- genetic / CRISPR: `#5E88C1`

### Status palette

- significant: `#111111`
- non-significant: `#BDBDBD`
- underpowered: `#D9D9D9`
- excluded: `#F0F0F0`

### Heatmap gradient

- low: `#F7F7F7`
- mid: `#C9D9EE`
- high: `#5E88C1`

### FM disambiguation

- FM keeps the frozen lavender-family hue.
- FM is distinguished from Pathway by hollow points, dashed lines, and
  outline-only fills where both appear together.
- FM must not receive a new hue to solve local overlay conflicts.

## Typography and export defaults

- Shared manuscript figure font: `Arial`
- Shared theme default: `theme_m2m(base_family = "Arial")`
- Shared export default: `export_editable_pdf(..., font_family = "Arial")`
- Final manuscript export should fail clearly if `Arial` cannot be resolved as
  an installed system font on the local system.

## Figure-dimension policy

The frozen default dimension contract is stored in
`plotting/R/config/figure_dimensions.csv`.

- Default width for every figure: `178 mm`
- Width hard limit: `180 mm`
- Heights are figure-specific defaults and may be overridden explicitly during
  render calls when layout pressure requires it.
- No figure is split before rendering in this freeze.

## Export helper contract

Shared helper:

`plotting/R/shared/export.R::export_editable_pdf(plot, output_path, width_mm = 178, height_mm, device = "pdf", transparent_bg = FALSE, font_family = "Arial")`

Required behavior:

- PDF-only export path
- mm-to-inch conversion internally
- hard stop for width above `180 mm`
- no rasterization by default
- Adobe-friendly vector PDF options via `grDevices::cairo_pdf()`
- resolve `Arial` through `systemfonts::match_font()` before final export
