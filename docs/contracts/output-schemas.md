## Task1 Proposal — 追加小节 A：Outputs & Field Contracts (Plot-ready + Audit-ready)

### A.1 Plot-ready（R 友好长表；用于 Figure 1/2）

**(T1-P1) `task1_leaderboard_long.csv`（主 leaderboard/scoreboard）**

* **主键**：`scope, dataset_or_direction, perturbation_type, representation, metric_name`

  * `scope ∈ {internal, cross}`
  * `dataset_or_direction`：internal 取 `{LINCS, scPerturb}`；cross 取 `{LINCS_to_scPerturb, scPerturb_to_LINCS}`
* **必需字段**：

  * `metric_value`
  * 分母：`n_total, n_valid, n_excluded, N_gallery_max`
  * `cross_alignment_contract`（cross 行必填）
* **metric_name（冻结枚举）**：

  * retrieval corrected：`mean_mrr_corrected, mean_hit1_corrected, mean_hit5_corrected, mean_hit10_corrected`
  * group-level：`mean_cosine_centroid, mean_pcc_centroid, mean_edist_biascorr`

**(T1-P2) `task1_internal_vs_cross_boxplot.parquet`（箱线图数据源，轻量）**

* **每行**：一个 query 的一个指标
* **必需字段**：`scope, dataset_or_direction, perturbation_type, representation, metric_name, metric_value`
* **建议字段**：`m_pos, N_gallery`

**(T1-P3) `task1_data_inventory_long.csv`（Task1 数据生态分布）**

* **主键**：`dataset, perturbation_type, cell_line, target_token`
* **字段**：`n_instances, source_db`

---

### A.2 Audit-ready（bounded audit 最小证据链）

**(T1-A1) `task1_retrieval_per_query.parquet`（最小可证伪证据）**

* **主键**：`scope, dataset_or_direction, perturbation_type, representation, query_uid`
* **必需字段**：

  * `cell_line, target_token`
  * `N_gallery, N_gallery_max, gallery_group_ids_hash`
  * `m_pos, rank_true`（strict tie）
  * raw：`mrr_raw, hit1_raw, hit5_raw, hit10_raw`
  * expected：`expected_mrr_chance, expected_hit1_chance, expected_hit5_chance, expected_hit10_chance`
  * corrected：`mrr_corrected, hit1_corrected, hit5_corrected, hit10_corrected`
  * policy：`loo_policy`（internal=`strict_recompute`；cross=`disjoint_no_leakage`）
  * cross：`cross_alignment_contract`（cross 行必填）

**(T1-A2) `task1_retrieval_summary.csv`**

* **键**：`scope, dataset_or_direction, perturbation_type, representation`
* **必需字段**：

  * 分母：`n_total, n_valid, n_excluded_missing_metric_or_mpos0, N_gallery_max`
  * 均值：`mean_mrr_raw, mean_expected_mrr_chance, mean_mrr_corrected`
  * hit@k（k=1/5/10）同理：raw/expected/corrected

**(T1-A3) `task1_chance_identity_check.csv`（tol=1e-12）**

* **键**：同 summary
* 字段：`delta_mrr, abs_delta_mrr, delta_hit1/5/10, abs_delta_hit1/5/10`

**(T1-A4) `task1_attrition.csv`（soft-fail 归因）**

* 键：`scope, dataset_or_direction, perturbation_type, representation, reason`
* 字段：`n_dropped, n_total_before, notes`
* reason 枚举建议：`missing_metric, m_pos0_after_intersection, cross_matched_keys_lt_5, group_too_small_for_split_half, edist_insufficient_cells, ...`

**(T1-A5) `task1_cross_alignment_proof.csv`**

* 必需字段：

  * `cross_alignment_contract="global_idx_lincs + sc_delta_row_idx"`
  * `n_matched_keys`（按 perturbation_type）
  * `eligible_bool`（`n_matched_keys >= 5`）
  * `excluded_reason`（如 `<5`）

---

## Task1 Proposal — 追加小节 B：Assertions / Fail-fast / Evidence Bundle

### B.1 Frozen assertions（必要且可证伪）

1. Cross alignment：`cross_alignment_contract` 恒等正确（cross 全行）
2. Cross eligibility gate：对每个 perturbation_type，若 `n_matched_keys<5` → 必须写入 attrition 且不入指标计算（本数据下 chemical 将满足该条件并被排除，这是事实结果）
3. LOO policy：internal=`strict_recompute`；cross=`disjoint_no_leakage`
4. Chance identity：`max(abs_delta_*) <= 1e-12`（mrr 与 hit@1/5/10）
5. Denominator 完整：`n_total == n_valid + n_excluded_missing_metric_or_mpos0`；`n_valid==0` 时均值必须 NA
6. E-distance：`EDIST_MAX_N=256`，subsample 由 seed=619 单次生成并复用；多线程只负责计算
7. Split-half（internal）：group size <4 → soft-fail NA + attrition（不终止全局）

### B.2 Evidence bundle（每个 run 最少输出）

* `run_manifest.json`（含 seed、inputs 列表、script_path、git_head、config）
* `audit_assertions.json`（断言 pass/fail + ≤5 counterexamples）
* `manifest.json`（本 step 文件清单）
* 本节 A.1/A.2 定义的 canonical outputs 全套

---

---

## Task2 Contract — Outputs & Field Contracts (Corrected Multisource v1)

### A.0 Schema preamble

Task2 must distinguish:

* **row identity fields** for audit, snapshot, and per-query evidence, such as
  `row_id`, `treated_cell_id`, `global_idx_lincs`, and `query_uid`
* **analysis/cohort keys** for benchmark metrics, centered on
  `dataset, cell_line, target_token`

Chemical rows may belong to multiple `target_token` memberships through
explode-matching, but that does not create new snapshot row identities.

All Task2 analysis-facing tables must carry `dataset`.

Representation-scope policy is explicit:

* LINCS supports `Gene` and `Pathway` only
* scPerturb K562 supports `Gene`, `Pathway`, `scgpt`, `geneformer`, `scbert`,
  `scfoundation`, `uce`, `state`, and `tahoe-x1`

Unsupported representation rows are **not-applicable scope**, not silent missing
data and not attrition.

Corrected Task2 v1 S4 and S5 metrics are stratified by `dataset` and
`cell_line`. Raw cross-dataset pooling is not part of the v1 core metric
contract.

### A.1 Plot-ready + synthesis-layer canonical tables

**(T2-P1) `task2_pairs_coverage.csv`（配对覆盖事实表）**

* 每行：`dataset, cell_line, target_token`
* 字段：`n_chem_instances, n_gen_instances, is_eligible_bool, source_tag`

**(T2-P2) `task2_post_build_inventory.csv`（helper / audit inventory，可选）**

* 每行：`dataset, perturbation_type, cell_line, representation, stage`
* 字段：`n_total_instances, n_valid_instances, n_dropped, drop_reason_breakdown`
* 口径：辅助盘点表，不是 corrected Task2 v1 的必需 core output

**(T2-P3) `task2_group_concordance_long.csv`（target-level group 主表）**

* 主键：`dataset, cell_line, target_token, representation, metric_name`
* 字段：`metric_value`
* 分母：`n_chem_instances_used, n_gen_instances_used, n_chem_sub, n_gen_sub`
* `metric_name`：`cosine_centroid, pcc_centroid, edist_biascorr`

**(T2-P4) `task2_group_leaderboard.csv`（representation-level group summary）**

* 主键：`dataset, cell_line, representation`
* 字段：`mean_cosine_centroid, mean_pcc_centroid, mean_edist_biascorr`
* 约束：`mean_edist_biascorr` 可保留，但必须 informational-only，不得 rank

**(T2-P5) `task2_retrieval_summary_long.csv`（direction-specific retrieval 主表）**

* 主键：`dataset, cell_line, direction, representation, metric_name`
* 字段：`metric_value`
* 分母：`n_valid, N_gallery_mean, N_gallery_max, m_pos_mean, m_pos_p50, m_pos_p90`

**(T2-P6) `task2_retrieval_leaderboard.csv`（direction-specific leaderboard）**

* 主键：`dataset, cell_line, direction, representation`
* 字段：`mean_mrr_corrected, mean_hit1_corrected, mean_hit5_corrected, mean_hit10_corrected`
* 排名：仅允许在同一 `dataset, cell_line, direction` 内排名

**(T2-P7) `task2_benchmark_summary_long.csv`（synthesis-layer union）**

* 主键：`analysis_family, dataset, cell_line, direction, representation, metric_name`
* group rows 来源：`task2_group_leaderboard.csv`
* retrieval rows 来源：`task2_retrieval_leaderboard.csv`
* 约束：不得从 target-level `task2_group_concordance_long.csv` 直接构造 group leaderboard rows

### A.2 Deferred / approval-needed Task2 tables

下列输出不属于 corrected Task2 v1 immediate canonical tables，除非另行审批：

* `task2_group_lift_long.csv`
* `task2_time_dose_distribution.csv`
* `rep_comparison_effects_long.csv`
* `benchmark_health_contribution.csv`
* `benchmark_health_summary.csv`

---

## Project Benchmark Synthesis — Outputs & Field Contracts (S7)

### S7.1 `project_input_registry.csv`

* 每行：一个 direct-ingest stage 或 transitive-provenance stage
* 主键：`provenance_role, task, stage_name, run_id`
* 必需字段：

  * `stage_dir`
  * `reporting_tables_json`
  * `validation_tables_json`
  * `manifest_path`
  * `audit_assertions_path`
  * `stage_manifest_path`
  * `audit_status`
  * `non_pass_assertions_json`
  * `notes`

### S7.2 `project_benchmark_summary_long.csv`

* 主键：
  `task, task_scope, analysis_family, dataset, cross_direction, cell_line, perturbation_type, direction, representation_raw, metric_name, row_origin`
* 必需字段：

  * `representation_canonical`
  * `metric_value`
  * `rank_value`
  * `rank_basis_metric_name`
  * `leaderboard_eligible_bool`
  * `cross_representation_comparable_bool`
  * `ingest_table, ingest_stage, ingest_run_id`
  * `source_table, source_stage, source_run_id`
  * `caution_codes`
* 分母/支持集字段按上游 reporting 表保留；若上游不存在，对应列可为空

### S7.3 `project_axis_score_inputs_long.csv`

* 每行：一个 scorecard axis family 输入或 contract-exclusion row
* 主键：`task_axis, family_id, representation_raw, row_origin`
* 必需字段：

  * `task, task_scope, dataset, cross_direction, cell_line, perturbation_type, direction`
  * `representation_canonical`
  * `metric_name, metric_value`
  * `coverage_status`
  * `aggregation_weight`
  * `cross_alignment_contract`
  * `ingest_table, ingest_stage, ingest_run_id`
  * `source_table, source_stage, source_run_id`
  * `caution_codes`

### S7.4 `project_representation_scorecard.csv`

* 主键：`representation_canonical`
* 必需字段：

  * `task1_axis_value_raw`
  * `task1_axis_observed_families`
  * `task1_axis_contract_exclusion_families`
  * `task1_axis_expected_families`
  * `task2_axis_value_raw`
  * `task2_axis_observed_families`
  * `task2_axis_contract_exclusion_families`
  * `task2_axis_expected_families`
  * `x_percentile`
  * `y_percentile`
  * `scorecard_eligible_bool`
  * `caution_codes`
* 约束：

  * percentile universe 只包含 `scorecard_eligible_bool=True` 的表示
  * `scorecard_eligible_bool=False` 的行必须保留观测到的 Task1 家族计数，且不得伪造 Task2 轴分数

### A.3 Audit-ready（bounded audit 最小证据链）

**(T2-A1) `pair_list.parquet`（scPerturb K562 pairing 证据）**

* 行身份字段：`treated_cell_id, control_cell_id, control_rank`
* 必需字段：`n_controls_used, dataset_side, perturbation_class, cell_line, target_raw, target_tokens, time, dose_value, specificity_tier, seed`
* 规则：无放回；`n_controls_used=min(pool,50)`；pool=0 -> attrition

**(T2-A2) `delta_meta.csv` + `gene_delta.npy` + `pathway_delta.npy`**

* `delta_meta.csv` 必含：`row_id, treated_cell_id, perturbation_class, cell_line, target_raw, time, dose_value, specificity_tier, n_controls_used`

**(T2-A3) `fm/<model>/fm_delta_meta.csv`（每模型）**

* 必含：`row_id, treated_cell_id, valid_mask, n_controls_used`

**(T2-A4) `task2_group_concordance.csv`**

* 主键：`dataset, cell_line, target_token, representation`
* 必需字段：`cosine_centroid, pcc_centroid, edist_biascorr`
* 分母：`n_chem_instances_used, n_gen_instances_used, n_chem_sub, n_gen_sub`

**(T2-A5) `task2_group_attrition.csv`**

* 主键：`dataset, cell_line, target_token, representation, reason`
* 字段：`n_dropped, n_total_before, notes`

**(T2-A6) `task2_retrieval_per_query.parquet`（multi-positive 证据链）**

* 主键：`dataset, cell_line, direction, representation, query_uid`
* 必需字段：
  * `query_uid`
  * `cell_line`
  * C2G：`query_target_raw, query_n_targets, query_time, query_dose_value`
  * G2C：`query_target_token`
  * `N_gallery, m_pos, rank_true`
  * raw / expected / corrected（mrr 与 hit@1/5/10 全套）
  * `pos_definition_id`
* 语义冻结：
  * G2C gallery items are chemical target centroids within the same `dataset, cell_line` slice
  * G2C does not depend on `Drug_meta.csv`, `perturbation_raw`, or any perturbagen-context key

**(T2-A7) `task2_retrieval_summary.csv`**

* 主键：`dataset, cell_line, direction, representation`
* 必需字段：`n_total, n_valid, n_excluded_missing_metric_or_mpos0, N_gallery_max`
* 均值字段：raw / expected / corrected 全套
* G2C 口径：`N_gallery = number of chemical target centroids within the same dataset/cell_line slice after validity filtering`
* 上述 G2C `N_gallery` 定义取代旧的 chemical-context-count 解释

**(T2-A8) `task2_retrieval_attrition.csv`**

* 主键：`dataset, cell_line, direction, representation, reason`
* 字段：`n_dropped, n_total_before, notes`

**(T2-A9) `task2_chance_identity_check.csv`（tol=1e-12）**

* 主键：`dataset, cell_line, direction, representation`
* 字段：`delta_mrr, abs_delta_mrr, delta_hit1/5/10, abs_delta_hit1/5/10`

## Task2 Contract — Multi-positive Retrieval Assertions

### B.1 Frozen definitions

* `score(q,c) = -||q-c||^2`
* `rank_true = 1 + #{ item : score(item) > best_positive_score }`
* `m_pos = |Pos|`
* chance correction uses exact expectation from `(N_gallery, m_pos)`
* these retrieval definitions apply uniformly to `LINCS` and `scPerturb`

### B.2 C2G（Chemical -> Genetic）

```text
Pos := target_tokens(query) intersect gallery_targets
m_pos := |Pos|
best_positive_score := max_{p in Pos} score(q, centroid[p])
rank_true := 1 + #{t in gallery_targets : score(q, centroid[t]) > best_positive_score}
raw := {MRR, Hit@1/5/10}; expected := E[.|N_gallery,m_pos]; corrected := raw-expected
```

### B.3 G2C（Genetic -> Chemical）

```text
gallery items are chemical target centroids within the same dataset and cell_line
N_gallery := number of chemical target centroids within the same dataset/cell_line slice after validity filtering

Pos := {t : t == query_target_token}
m_pos := |Pos|
best_positive_score := max_{p in Pos} score(q, centroid[p])
rank_true := 1 + #{t in gallery : score(q, centroid[t]) > best_positive_score}
raw/expected/corrected as above
```

Corrected Task2 S5 v1 does not reconstruct chemical-context galleries from
`Drug_meta.csv`, `perturbation_raw`, or any perturbagen-context key. `time` and
`dose_value` remain instance metadata / descriptive fields only, not
gallery-defining keys.

### B.4 Frozen assertions

1. `dataset` must exist on every Task2 analysis-facing row.
2. Target-level group rows must carry `dataset, cell_line, target_token`.
3. Retrieval summary and leaderboard rows must carry `dataset, cell_line, direction`.
4. `pos_definition_id` must be fixed and versioned.
5. Denominator rows entering retrieval must satisfy `1 <= m_pos <= N_gallery`.
6. chance identity must satisfy `max(abs_delta_*) <= 1e-12`.
7. Corrected Task2 v1 must not raw-pool LINCS and scPerturb in S4 or S5.
