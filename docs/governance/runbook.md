## 1. Pipeline Step Map（脚本数量与步骤严格对齐）

**原则**：每个 step 只有一个入口脚本；新增 step 才允许新增脚本；逻辑复用必须进入 `src/m2mbench/`。

| Step | Script (唯一入口)                                   | Purpose                                     | Canonical outputs (must exist)                                                                                                                             |
| ---- | ----------------------------------------------- | ------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| S0   | `scripts/s0_build_data_inventory.py`            | 盘点原始输入生态（不依赖 Task2 snapshot）                | `result_data_inventory_long.csv`, `data_source_manifest.csv`                                                                                               |
| S1   | `scripts/s1_task1_internal_metrics.py`          | Task1 internal 指标                           | `task1_retrieval_per_query.parquet`, `task1_retrieval_summary.csv`, `task1_chance_identity_check.csv`, `task1_leaderboard_long.csv`, `task1_attrition.csv` |
| S2   | `scripts/s2_task1_cross_metrics.py`             | Task1 cross matched pairs（eligibility gate） | `task1_cross_*` 全套 + `task1_cross_alignment_proof.csv` + attrition                                                                                         |
| S3   | `scripts/s3_build_task2_snapshot.py`            | 派生 Task2 pairs + K562 pairing/delta/FM      | `task2_lincs_pairs.csv`, `pair_list.parquet`, `delta_meta.csv`, `fm_delta_meta.csv`, `task2_post_build_inventory.csv`, `task2_pairs_coverage.csv`          |
| S4   | `scripts/s4_task2_group_concordance.py`         | Task2 group concordance                     | `task2_group_concordance.csv`, `task2_group_concordance_long.csv`, `task2_group_lift_long.csv`                                                             |
| S5   | `scripts/s5_task2_retrieval.py`                 | Task2 retrieval（C2G/G2C multi-positive）     | `task2_retrieval_per_query.parquet`, `task2_retrieval_summary.csv`, `task2_chance_identity_check.csv`, `task2_retrieval_summary_long.csv`                  |
| S6   | `scripts/s6_rep_comparison_and_significance.py` | Gene baseline 显著性（valid UID 交集）             | `rep_comparison_effects_long.csv`（含 `n_overlap_pairs`）                                                                                                     |
| S7   | `scripts/s7_benchmark_health_and_story_pack.py` | Benchmark health + story pack 汇总            | `benchmark_health_contribution.csv`, `benchmark_health_summary.csv`, `story_pack/`                                                                         |

> 注：S0 不包含 `task2_post_build_inventory.csv`（由 S3 生成），以消除“Task2 snapshot 尚未存在”的冲突。

---

## 2. 目录结构（冻结）

* `data/task1_snapshot_v1/`：只读；Task1 输入快照
* `data/task2_snapshot_v1/`：只读（由 S3 build 后冻结）；包含 K562 derived 与 FM outputs
* `src/m2mbench/`：复用库代码（metrics、chance correction、io schema、assertions）
* `scripts/`：仅入口脚本（S0–S7）；**禁止复制大段逻辑**
* `runs/<run_id>/<Sx_name>/`：每步独立输出目录
* `runs/_registry/`：运行登记（见清理策略）
* `docs/`：contracts、pipeline、governance、story 等文档

---

## 3. 脚本代码形式（冻结）

每个 `scripts/Sx_*.py` 必须满足：

1. 单一入口：`if __name__ == "__main__": main()`
2. 顶部 docstring 必写：Inputs / Outputs / Frozen constants / Attrition rules
3. 统一 CLI：至少包含 `--project-root`, `--run-id`, `--seed`（默认 619）
4. 输出三件套：

   * `run_manifest.json`（必须含 `run_id, stage, script_path, git_head, config(seed等), inputs(paths)`）
   * `audit_assertions.json`（断言结果 + counterexamples ≤5）
   * `manifest.json`（本 step 文件清单 + size）
5. 所有输出只能写入 `runs/<run_id>/<Sx_name>/`（不得写到其他目录）
6. 多线程可用于计算加速，但不得影响 subsample/seed/排序决定（随机性仅由 seed=619 决定）

---

## 4. 结果表的“canonical”原则（冻结）

* 所有 plot-ready 表必须是 **long-format**（便于 R 画图）
* 所有 summary 表必须携带分母字段（n_total/n_valid/n_excluded/m_pos/N_gallery）
* per-query 表用于审计，不强制 long，但字段必须完整
* 任何“便捷别名文件”不得替代 canonical outputs

---

## 5. 旧结果清理策略（冻结）

**原则**：脚本内部不得自动删除 runs；由显式 cleanup 脚本执行。

* 维护 registry：`runs/_registry/story_runs.csv`，字段：

  * `run_id, stage, status ∈ {experimental,frozen,superseded}, created_at, notes`
* 新 run 默认 `experimental`
* 仅当人工确认后可标记为 `frozen`
* `scripts/cleanup_runs.py --status experimental --keep-last N`

  * 只删除 registry 中 status=experimental 的 runs
  * 不允许删除 frozen runs
  * 删除前必须打印待删除清单（dry-run 默认）

---

## 6. Step 级关键断言摘要（冻结）

* S0：输入文件与关键列存在
* S1：LOO strict_recompute；chance identity <=1e-12；分母守恒
* S2：cross alignment contract；eligibility gate（matched_keys<5 → attrition）；chance identity <=1e-12
* S3：无放回 pairing；动态 `n_controls_used`；FM delta 依赖实际 counts；输出 post-build inventory
* S6：paired test 在 valid UID 交集；输出 n_overlap_pairs
* S7：dominance 只用 corrected；连续 corrected 指标 Gini 前非负化；含 gini