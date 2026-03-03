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

## Task2 Proposal — 追加小节 A：Outputs & Field Contracts (Plot-ready + Audit-ready)

### A.1 Plot-ready（R 友好长表；用于 Figure 3/4）

**(T2-P1) `task2_pairs_coverage.csv`（配对覆盖事实表）**

* 每行：`dataset, cell_line, target_token`
* 字段：`n_chem_instances, n_gen_instances, is_eligible_bool, source_tag`（LINCS 派生 / K562 common targets）

**(T2-P2) `task2_post_build_inventory.csv`（build 后真实有效样本量盘点）**

* 每行：`dataset, perturbation_type, cell_line, representation, stage`
* 字段：`n_total_instances, n_valid_instances, n_dropped, drop_reason_breakdown`
* 目的：避免 S0 阶段对不存在的 Task2 snapshot 做统计

**(T2-P3) `task2_group_concordance_long.csv`（跨机制 group-level 主表）**

* 主键：`dataset, cell_line, target_token, representation, metric_name`
* 字段：`metric_value`
* 分母：`n_chem_instances_used, n_gen_instances_used, edist_n_chem_sub, edist_n_gen_sub`
* metric_name：`cosine_centroid, pcc_centroid, edist_biascorr`（PCC 允许但需注释解释性有限）

**(T2-P4) `task2_retrieval_summary_long.csv`（C2G/G2C leaderboard）**

* 主键：`direction, dataset, representation, metric_name`（direction ∈ {C2G,G2C}）
* 字段：`metric_value`
* 分母：`n_valid, N_gallery_mean, N_gallery_max, m_pos_mean, m_pos_p50, m_pos_p90`

**(T2-P5) `task2_group_lift_long.csv`（cell_line/target 偏好）**

* 主键：`dataset, direction, representation, group_type, group_value`
* 字段：`n_support, success_rate, global_success_rate, lift_vs_global, success_def`
* success_def（冻结）：`success_bool_both := (mrr_corrected>0) AND (hit10_corrected>0)`

**(T2-P6) `task2_time_dose_distribution.csv`（分布表；bins 后置）**

* 每行：`dataset, cell_line, time_value, dose_value, n_instances, n_unique_contexts`

**(T2-P7) `rep_comparison_effects_long.csv`（Gene baseline effect + p/q）**

* 主键：`task, scope_or_direction, dataset, metric_name, representation`
* 字段：`effect_vs_gene, n_overlap_pairs, test_method, p_value, q_value`

**(T2-P8) Benchmark health（见 Story 合同）**

* `benchmark_health_contribution.csv`, `benchmark_health_summary.csv`（含 gini）

---

### A.2 Audit-ready（bounded audit 最小证据链）

**(T2-A1) `pair_list.parquet`（K562 pairing 证据）**

* 必需字段：`treated_cell_id, control_cell_id, control_rank, n_controls_used, dataset_side, perturbation_class, cell_line, target_raw, target_tokens, time, dose_value, specificity_tier, seed`
* 规则：无放回；`n_controls_used=min(pool,50)`；pool=0→attrition

**(T2-A2) `delta_meta.csv` + `gene_delta.npy` + `pathway_delta.npy`**

* `delta_meta.csv` 必含：`row_id, treated_cell_id, perturbation_class, cell_line, target_raw, time, dose_value, specificity_tier, n_controls_used`

**(T2-A3) `fm/<model>/fm_delta_meta.csv`（每模型）**

* 必含：`row_id, treated_cell_id, valid_mask, n_controls_used`

**(T2-A4) `task2_retrieval_per_query.parquet`（multi-positive 证据链）**

* 主键：`direction, dataset, representation, query_uid`
* 必需字段：

  * `cell_line`
  * C2G：`query_target_raw, query_n_targets, query_time, query_dose_value`
  * G2C：`query_target_token`（genetic token）
  * `N_gallery, m_pos, rank_true`（strict tie）
  * raw/expected/corrected（mrr 与 hit@1/5/10 全套）
  * `pos_definition_id`（固定字符串，防漂移）

**(T2-A5) `task2_chance_identity_check.csv`（tol=1e-12）**

* delta/abs_delta 全套（mrr 与 hit@1/5/10）

**(T2-A6) `task2_attrition.csv`**

* reason 枚举建议：`no_controls_available, fm_valid_mask_drop, m_pos0_after_intersection, missing_one_side, ...`

---

## Task2 Proposal — 追加小节 B：Multi-positive Retrieval Contract（伪代码 + 断言）

### B.1 冻结定义

* score：`score(q,c) = -||q-c||^2`（Gene/Pathway/FM delta 均用欧式距离；FM delta 已由 operator 得到）
* strict rank：`rank_true = 1 + #{ item : score(item) > best_positive_score }`
* multi-positive：`m_pos = |Pos|`；若 `m_pos==0` → attrition
* chance correction：精确期望基于 `(N_gallery, m_pos)`；identity check tol=1e-12

### B.2 C2G（Chemical → Genetic）

```text
Pos := target_tokens(query) ∩ gallery_targets
m_pos := |Pos|
best_positive_score := max_{p in Pos} score(q, centroid[p])
rank_true := 1 + #{t in gallery_targets : score(q, centroid[t]) > best_positive_score}
raw := {MRR, Hit@1/5/10}; expected := E[·|N_gallery,m_pos]; corrected := raw-expected
```

### B.3 G2C（Genetic → Chemical；R-INSTANCE context centroid）

```text
gallery items are chemical contexts:
  context_id = (cell_line, perturbation_raw, time, dose_value)
  context_target_tokens = union(target_tokens of member chemical instances)

Pos := {ctx : query_target_token ∈ context_target_tokens(ctx)}
m_pos := |Pos|
best_positive_score := max_{p in Pos} score(q, centroid[p])
rank_true := 1 + #{ctx in gallery : score(q, centroid[ctx]) > best_positive_score}
raw/expected/corrected as above
```

### B.4 断言（冻结）

1. `pos_definition_id` 必须为固定字符串（例如 `C2G_tokens_intersect_gallery_targets_v1`、`G2C_target_in_context_tokens_v1`）
2. 进入分母的行必须满足 `1 <= m_pos <= N_gallery`
3. chance identity：`max(abs_delta_*) <= 1e-12`
4. **显著性检验交集（冻结）**：paired test 必须在 `valid_uids_gene ∩ valid_uids_rep` 上进行；输出 `n_overlap_pairs`；若不足阈值 → soft-fail（p/q NA + reason）