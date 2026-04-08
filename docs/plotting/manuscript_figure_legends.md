# Manuscript Figure Legends

## Figure 1

### Panel 1B
Compact sunburst panel, rendered separately for LINCS and scPerturb, showing how the benchmark decomposes from Task to perturbation type and then to ranked cell-line or target-context slices. Segment area is proportional to eligible support, small slices are explicitly folded into Others, center labels report the dataset-level support summary, the outer target ring now keeps the cell-line hue instead of washing out to near-white, and the color legend is exported as a standalone companion PDF.

### Deferred schematic note
The workflow schematic and the data-filtering schematic are described in markdown only and are intentionally not embedded in the current Figure 1 main render.

## Figure 2

### Panel 2A
Stacked bars summarize the Task1 slices that are actually evaluable under the frozen contract. Internal rows retain lawful FM-local scope where permitted, whereas cross-chemical eligibility is shown explicitly as an excluded bar rather than being reframed as attrition.

### Panel 2B
Scoreboard-style benchmarking matrix for Task1 internal performance. Rows are methods, columns are dataset-metric benchmark groups, bar height tracks performance magnitude, and grayscale intensity tracks within-group rank.

### Panel 2C
Unified matched-unit Gene-versus-Pathway comparison for Task1 across Internal LINCS, Internal scPerturb, LINCS -> scPerturb, and scPerturb -> LINCS. Facets preserve explicit metrics while the plotted points and intervals retain the paired matched-unit structure used for the test.

### Panel 2D
Internal-only matched-unit summary panel for Task1. Group-concordance and Retrieval families are displayed in separate facet blocks with independent value ranges; Gene and Pathway are shown as compact internal summaries over shared `(dataset, cell_line, target)` units; significance labels report only the internal paired test result (`q`, `n`, and signed median delta), and any eDist-derived row is plotted on the logged transform.

### Panel 2E
Internal-only paired enrichment suitability panel for Task1 cell lines. Rows are split by perturbation type (`Chemical`, `Genetic`), columns are split by dataset (`LINCS`, `scPerturb`), and each dataset x perturbation-type block contributes `3` high and `3` low cell-line examples selected independently. Gene and Pathway are plotted as paired points joined by a dumbbell segment on the enrichment axis, right-side labels report support `n`, and shared exemplars are flagged directly in-row rather than separated into a secondary heatmap.

### Panel 2F
Internal-only paired enrichment suitability panel for Task1 targets. The panel mirrors Panel 2E, but the displayed entities are target-level rows rather than cell lines; dual-tail selection is still performed independently inside each dataset x perturbation-type block, and the plotted targets remain atomic target labels rather than pooled target-family aliases.

## Figure 3

### Panel 3B
Two-scoreboard Panel 3B for Task2. The upper board keeps the LINCS common-scope benchmark view with rows as cell lines and paired Gene/Pathway lanes across group and C2G retrieval metrics, whereas the lower board switches scPerturb to a K562 local scoreboard with rows as representations so the panel no longer collapses to a single misleading point.

### Panel 3C
LINCS-only ranked suitability panel for Task2 cell lines. Cell lines are ordered by a paired suitability summary derived from the pattern surface, Gene and Pathway are shown as a dumbbell pair on the same axis, and the right margin reports support `n` for each ranked row. This panel no longer uses dual-tail exemplar selection.

### Panel 3D
Dataset-faceted ranked suitability panel for Task2 targets. Targets are ordered within each dataset by the paired suitability summary, shared LINCS/scPerturb targets are flagged explicitly, and the panel uses the same shared-axis Gene-versus-Pathway dumbbell encoding as Panels 2E, 2F, and 3C.

### Panel 3E
Performance-only paired Gene-versus-Pathway comparison for Task2, arranged as dataset x metric small panels. Raw points and connecting segments represent matched `(cell_line, target)` units, compact intervals summarize the distribution within each representation, and the top annotation reports only the frozen paired test result (`q`, `n`, `delta`). Both LINCS and scPerturb therefore use target-level matched units rather than cell-line-only aggregation, and the panel is restricted to the performance layer rather than the pattern layer.

### Panel 3F
FM-local absolute-performance panels for the supported scPerturb K562 subset only. Each point-and-interval summarizes one method on the absolute metric scale, dashed vertical guides mark the Gene reference value for that metric, and right-side annotations retain the Gene-reference significance layer without turning the panel into a benchmark-wide delta-versus-Gene effect plot.

## Extended Figures

### Extended Figure 1
Support registry for the Figure 1 contract, covering the composition table, lawful-scope matrix, representation/modifier availability table, and retained object-map bookkeeping that were removed from the main figure.

### Extended Figure 4
Shared-only Task1 enrichment heatmaps aligned to the revised Figure 2E and Figure 2F reading. `EF4A` summarizes shared cell-line rows and `EF4B` summarizes shared target rows, while the full decomposition remains available as support tables rather than crowded PDF panels.

### Extended Figure 6
Matrix-style Task2 suitability decomposition. `EF6A` expands the LINCS cell-line surface behind Panel 3C with one row per cell line, and `EF6B` expands the dataset-faceted target surface behind Panel 3D with one row per target, replacing the older repeated-row dumbbell grammar.

### Extended Figure 7
Modifier-family supplementary figure covering the formal Task2 C2G time, dose, and target-multiplicity result families plus their unified comparison statistics. Descriptive panels now use matched `(cell_line, target)` target-level units rather than cell-line-only summaries, and the `EF7D` summary switches to a dose/time volcano-style view.

### Extended Figure 8
Supplementary reference block combining the corrected target-level Gene-versus-Pathway comparison statistics (`EF8A`), FM-local detail synchronized to Panel 3F (`EF8B`), and the moved Task1 contextual reference surface (`EF8C`). `EF8A` follows the same matched `(cell_line, target)` significance contract as Panel 3E rather than the older cell-line-only aggregation.

### Extended Figure 9
Support-only interpretation block for the biology-facing readout of Figure 3. `EF9A` diagnoses whether suitability summaries are merely support-size artifacts by plotting support versus suitability for both the cell-line and target surfaces, while `EF9B` collects representative target archetypes spanning shared anchors, cross-dataset gaps, and dataset-specific high-suitability examples. This figure remains interpretive support and does not promote the manuscript to a mechanism-causality claim.
