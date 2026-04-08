# Manuscript Panel Revision Note V2

This note records the revision-facing contract for the standalone panel export
cycle under `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current/`. It is intended as the working
bridge between plot implementation and later legend writing.

- Review root: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current/`
- Panel manifest: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current/figures/panel_export_manifest.tsv`
- Semantic audit: `/mnt/NAS_21T/ProjectData/M2M/runs/_staging/manuscript_visual_revision_current/figures/spatial_revision_semantics.md`
- Export mode: standalone panel PDF only

Panels not called out below remain on their prior scientific contract and are
mainly tracked through the manifest audit layer (`EF1A-E`, `EF2A/C`,
`EF5A-C`, `EF7A-D`).

## Markdown-only notes

### Figure 1 data filtering note

- The Figure 1 data-filtering schematic is no longer exported as a PDF panel in
  the current revision cycle.
- Keep the filtering logic in markdown only for manual redrawing:
  Task1 support comes from `task1_data_inventory_long.csv`;
  Task2 support comes from `task2_pairs_coverage.csv`;
  the revised `1B` collapses task identity and stages the hierarchy as
  `Perturbation type -> Cell line -> Target gene`.
- Cell-line colors are assigned dataset-locally to the top-10 supported cell
  lines, with all remaining cell lines folded into `Others`.

### Figure 3B target-tier audit note

- The current canonical `figure3_panel_3b_c2g_performance_overview.csv` does
  not expose a benchmark-wide target-tier field suitable for trustworthy
  labeling such as `The Cleanest Hits`.
- Those labels exist only in the retired submission-refresh artifact
  (`plotting/R/submission_refresh/figure3Refresh.R`) and are therefore treated
  as non-canonical for the live manuscript panel.

### EF5A direction support note

- `EF5A` is treated as markdown-only support documentation in this revision.
- The plotted supplementary surface should focus on `C2G` full breakdown only:
  `EF5B = LINCS C2G full performance breakdown`,
  `EF5C = scPerturb/K562 C2G full performance breakdown`.

## Panel 1B

- 科学问题: 数据集内部的 benchmark 组成如何沿 `Perturbation type -> Cell line -> Target gene` 分解。
- 统计单元: dataset 内部的 combined Task1/Task2 support slice；task 来源保留在汇总统计里，但不再单独占据内圈。
- 图型: 紧凑 sunburst，每个 dataset 一个 standalone panel。
- Legend 提示: 强调中心是 dataset-level support summary，内圈是 perturbation，中圈是 dataset-local top-10 cell lines，外圈是 target staging；颜色 legend 单独导出。
- 注意事项: cell-line 颜色按 dataset-local top-10 固定，`Others` 用中性灰；target 只对高支持切片上弧标签，不再追求满标注。

## Panel 2D

- 科学问题: 在共享 Task1 internal bridge units 上，Gene 与 Pathway 的内部表现如何概括展示。
- 统计单元: shared dataset x cell_line x target internal triplets。
- 图型: internal-only summary panel，`Group` 与 `Retrieval` 分开分面。
- Legend 提示: 显著性标签只报告 internal matched-unit paired result (`q`, `n`, `delta`)。
- 注意事项: `eDist` 进入可视化与统计前统一走 log transform；不再保留 internal-vs-cross 连线语义。

## Panel 2E

- 科学问题: 哪些 cell lines 更适合 Task1 cell perturbation enrichment 分析。
- 统计单元: `dataset x perturbation_type x cell_line x representation` 汇总后，按 block 独立做 `3 high + 3 low` 双尾选择。
- 图型: 四宫格 paired dumbbell / lollipop；rows=`Chemical/Genetic`，cols=`LINCS/scPerturb`。
- Legend 提示: 右侧 `n=` 为 support size；同一行的 Gene 与 Pathway 是 paired enrichment score，不是两个独立榜单。
- 注意事项: `shared_flag` 只用于标识跨数据集共享 exemplar；选择语义必须保留 `selection_tail`。

## Panel 2F

- 科学问题: 哪些 targets 更适合 Task1 enrichment 分析。
- 统计单元: `dataset x perturbation_type x target x representation` 汇总后，按 block 独立做双尾选择。
- 图型: 与 2E 同构的四宫格 paired dumbbell / lollipop。
- Legend 提示: target 是原子 target label，不再把多个 target family 混成一个展示单元。
- 注意事项: shared target 标记沿用主图语法；support 列与 selection tail 语义必须保留。

## Panel 3B

- 科学问题: 在 Task2 lawful benchmark 上，LINCS common scope 与 scPerturb/K562 local scope 应分别如何展示。
- 统计单元: LINCS 部分为 dataset x cell_line x metric 的 common benchmark unit；scPerturb 部分为 K562 local representation x metric unit。
- 图型: 双 scoreboard 组合；上半部分是 LINCS common benchmark，下半部分是 scPerturb local scoreboard。
- Legend 提示: 这是 benchmark/local scoreboard 组合，不是统一 y 轴的棒棒糖图。
- 注意事项: scPerturb 不再被 cell-line inventory 强行拉成只有一个点的伪 benchmark 视图。

## Panel 3C

- 科学问题: 在 LINCS 中，哪些 cell lines 更适合 Task2 pattern-level suitability 分析。
- 统计单元: LINCS cell-line level suitability summary，由 pattern surface 汇总到 `row_rank` 与 `suitability_score`。
- 图型: ranked paired dumbbell / lollipop。
- Legend 提示: 排序依据是 suitability summary，而不是 dual-tail exemplar selection。
- 注意事项: 右侧 support 列要保留；主图只画 LINCS，不把 scPerturb 混入同一排名。

## Panel 3D

- 科学问题: 在 LINCS 与 scPerturb 中，哪些 targets 更适合 Task2 pattern-level suitability 分析。
- 统计单元: dataset-faceted target suitability summary，按 dataset 内部独立排序。
- 图型: dataset-faceted ranked paired dumbbell / lollipop。
- Legend 提示: shared target 需显式标注，排序依据固定为 suitability summary。
- 注意事项: 这不是 exemplar top/bottom 面板；`row_rank` 与 `shared_flag` 都是 panel-defining 字段。

## Panel 3E

- 科学问题: 在 Task2 performance layer 上，Gene 与 Pathway 是否存在系统性差异。
- 统计单元: matched `(cell_line, target)` units；键为 `(dataset, cell_line, target, metric_name, direction)`。
- 图型: `dataset x metric` 小面板 distribution comparison，叠加 raw paired points 与紧凑 summary interval。
- Legend 提示: 顶部标注只解释 performance-layer paired test (`q`, `n`, `delta`)；不要混入 pattern-layer 叙述。
- 注意事项: `scPerturb` 的 `n` 必须是 K562 上合法 matched targets 数，不能再写成 cell-line-only `1`。

## Panel 3F

- 科学问题: 在受支持的本地 FM scope 中，FM 方法与 Gene baseline 的绝对表现如何。
- 统计单元: `scPerturb/K562` local-only representation x analysis_family x metric rows。
- 图型: absolute-performance point-range panel，配 Gene reference 垂直虚线。
- Legend 提示: 明确标注 `local-only`，避免读者把它误解成 benchmark-wide comparison。
- 注意事项: 右侧显著性注释只是对 Gene reference 的补充，不应把坐标轴改写成 delta-vs-Gene。

## Panel EF4A

- 科学问题: 共享 cell lines 在 Task1 enrichment 中的扩展表面是什么样。
- 统计单元: shared-only `dataset x perturbation_type x cell_line x representation` rows。
- 图型: shared-only heatmap。
- Legend 提示: 说明 full surface 结果保留在表格中，热图只承担 shared overview。
- 注意事项: 颜色读的是 enrichment score，不再复用主图的 paired dumbbell 语法。

## Panel EF4B

- 科学问题: 共享 targets 在 Task1 enrichment 中的扩展表面是什么样。
- 统计单元: shared-only `dataset x perturbation_type x target x representation` rows。
- 图型: shared-only heatmap。
- Legend 提示: 这是 shared overview；完整 target surface 仍以表格保留。
- 注意事项: target 保持 atomic display，但不再把上千个 shared targets 全部挤进图中。

## Panel EF6A

- 科学问题: LINCS cell-line suitability 的扩展矩阵表面是什么样。
- 统计单元: LINCS cell-line suitability rows。
- 图型: entity-by-metric heatmap，每个 cell line 只占一行。
- Legend 提示: 颜色表示列内 percentile，不再重复绘制同一 entity 的多条 y 轴文本。
- 注意事项: 继续按整体 suitability 排序，但把重复 row 压缩成矩阵。

## Panel EF6B

- 科学问题: Task2 target suitability 的扩展矩阵表面是什么样。
- 统计单元: dataset-faceted target suitability rows。
- 图型: dataset-faceted heatmap，每个 target 只占一行。
- Legend 提示: 颜色表示列内 percentile，面板主要承担 dense detail 阅读而不是 exemplar 叙事。
- 注意事项: 不恢复旧版 top/bottom exemplar grammar，也不再让同一 target 在 y 轴重复出现多次。

## Panel EF8A

- 科学问题: 修正后的 Task2 performance-layer Gene-vs-Pathway 配对统计在补充材料里如何以更可读的 signed-delta 形式展开。
- 统计单元: 与 3E 相同的 matched `(cell_line, target)` target-level performance units。
- 图型: supplementary signed-delta comparison panel，沿 `dataset | metric` 展开。
- Legend 提示: 明确它直接消费修正后的 target-level stats，而不是旧版 cell-line summary 配对；颜色表示方向，标签表示 paired-test 结果。
- 注意事项: `n_units`、`test_status`、`effect_size`、`median_delta`、`bh_q` 口径要与 3E 完全一致。

## Panel EF8B

- 科学问题: FM-local reference detail 如何与主图 3F 保持一致。
- 统计单元: `scPerturb/K562` local FM support rows。
- 图型: 与 3F 对齐的 supplementary local-detail panel。
- Legend 提示: 延续 `local-only` framing，不把它升级为 broader benchmark claim。
- 注意事项: 风格与注释体系需和 3F 同步。

## Panel EF8C

- 科学问题: 被移出主图的 Task1 contextual support reference 如何以逐表示法比较的方式补充呈现。
- 统计单元: Task1 internal contextual support reference rows，按 dataset x metric x representation 聚合到 distribution surface。
- 图型: multi-panel boxplot comparison surface。
- Legend 提示: 说明它是 reference surface，不参与主图 3E 的 performance-layer 结论；这里读的是 representation-wise support distribution。
- 注意事项: contextual support 的统计对象与 3E/EF8A 不同，不能混写；LINCS 与 scPerturb 的 lawful representation inventory 不同，因此 panel 允许不同 x 轴稀疏度。

## Panel EF9A

- 科学问题: Task2 的 suitability 结构是否只是 support size 的投影。
- 统计单元: Task2 cell-line 与 target full-surface rows，经同一 suitability summary 汇总后保留 `support_n` 与 shared flag。
- 图型: support-vs-suitability diagnostic scatter，按 surface 分面。
- Legend 提示: 这是诊断面，不是新的显著性家族，也不是机制因果图。
- 注意事项: 只用于排查 small-denominator 或 support-dominated 解释风险；不得把趋势线写成生物学规律。

## Panel EF9B

- 科学问题: 哪些 target archetypes 最能代表 shared anchor、cross-dataset gap 与 dataset-specific high-suitability 三类模式。
- 统计单元: 从 Task2 target pattern full surface 选择的 target exemplars；每个 exemplar 在 dataset 内保留 Gene/Pathway paired enrichment 与 support。
- 图型: dataset-faceted paired dumbbell exemplar panel。
- Legend 提示: exemplar 是 support-only illustration，用于帮助读者理解 pattern surfaces，不是机制 proof。
- 注意事项: selection rule 必须固定且可重放；不得把 exemplar 面板写成 causal mechanism validation。
