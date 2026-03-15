## 1. Pipeline Step Map

**原则**：每个 step 只有一个入口脚本；新增 step 才允许新增脚本；逻辑复用必须进入 `src/m2mbench/`。

### 1.1 当前保留的 legacy / interim Task2 stages

这些 stage 对应当前已审计的 scPerturb-K562 Task2 历史证据，必须保留，但不再代表 corrected Task2 v1 的最终 stage map。

| Legacy step | Script | Status | Purpose |
| ---- | ---- | ---- | ---- |
| S3-legacy | `scripts/s3_build_task2_snapshot.py` | preserved | scPerturb-K562 snapshot build only |
| S4-legacy | `scripts/s4_task2_group_concordance.py` | preserved | scPerturb-K562 group concordance only |
| S5-legacy | `scripts/s5_task2_retrieval.py` | preserved | scPerturb-K562 retrieval only |
| S6-legacy | `scripts/s6_task2_result_synthesis.py` | preserved | scPerturb-K562 synthesis only |

### 1.2 Corrected successor stage map for Task2

下表是 corrected multisource Task2 的 authoritative stage map。当前脚本文件名与 canonical outputs 以本地已物化实现为准。

| Step | Script | Purpose | Canonical outputs (must exist) |
| ---- | ---- | ---- | ---- |
| S0 | `scripts/s0_build_data_inventory.py` | 盘点原始输入生态（不依赖 Task2 snapshot） | `task1_data_inventory_long.csv`, `data_source_manifest.csv` |
| S1 | `scripts/s1_task1_internal_metrics.py` | Task1 internal 指标 | `task1_retrieval_per_query.parquet`, `task1_retrieval_summary.csv`, `task1_chance_identity_check.csv`, `task1_leaderboard_long.csv`, `task1_attrition.csv` |
| S2 | `scripts/s2_task1_cross_metrics.py` | Task1 cross matched pairs（eligibility gate） | `task1_cross_*` 全套 + `task1_cross_alignment_proof.csv` + attrition |
| S3 | `scripts/s3_build_task2_multisource_snapshot.py` | corrected Task2 multisource snapshot build | `task2_pairs_coverage.csv`, optional `task2_post_build_inventory.csv`, snapshot manifests |
| S4 | `scripts/s4_task2_group_concordance_multisource.py` | corrected Task2 group concordance | `task2_group_concordance.csv`, `task2_group_attrition.csv` |
| S5 | `scripts/s5_task2_retrieval_multisource.py` | corrected Task2 target-level retrieval（C2G/G2C multi-positive） | `task2_retrieval_per_query.parquet`, `task2_retrieval_summary.csv`, `task2_retrieval_summary_long.csv`, `task2_retrieval_attrition.csv`, `task2_chance_identity_check.csv` |
| S6 | `scripts/s6_task2_result_synthesis_multisource.py` | corrected Task2 synthesis / reporting | `task2_group_concordance_long.csv`, `task2_group_leaderboard.csv`, `task2_retrieval_leaderboard.csv`, `task2_benchmark_summary_long.csv` |
| S7 | `scripts/s7_project_benchmark_synthesis.py` | project-level benchmark synthesis | `project_input_registry.csv`, `project_benchmark_summary_long.csv`, `project_axis_score_inputs_long.csv`, `project_representation_scorecard.csv` |

### 1.3 Retired later-stage entries

早期探索性的 later-stage 产物不再属于 corrected Task2 v1 的
authoritative stage map。若未来要恢复这些分析，必须作为单独审批后的
later-stage work 重新引入，而不是恢复旧入口脚本。

---

## 2. Directory Structure

* `data/task1_snapshot_v1/`：只读；Task1 输入快照
* `data/task2_snapshot_v1/`：只读；当前 legacy/interim scPerturb-K562 Task2 证据
* `data/task2_snapshot_v2/`：corrected multisource Task2 successor snapshot root（当前已冻结并物化）
* `src/m2mbench/`：复用库代码（metrics、chance correction、io schema、assertions）
* `scripts/`：仅入口脚本；禁止复制大段逻辑
* `runs/<run_id>/<stage>/`：项目内的本地兼容路径；当前允许通过 symlink 解析到 NAS-backed active/archive storage
* `runs/_registry/`：运行登记
* `docs/`：contracts、pipeline、governance、story 等文档

---

## 2.1 Storage Backing

**原则**：authoritative benchmark outputs 可以迁移出本地磁盘，但本地项目路径兼容性必须保持。

* active authoritative results:
  * `/mnt/NAS_21T/ProjectData/M2M/runs`
* historical/archive results:
  * `/mnt/NAS_21T/ProjectData/M2M/archive/runs`
* archived reports:
  * `/mnt/NAS_21T/ProjectData/M2M/archive/reports`
* local compatibility:
  * `runs/` 当前保持为本地目录
  * `runs/<run_id>` 按 run-id 提供 symlink，可分别指向 active 或 archive namespace
  * `reports` 当前可作为指向 NAS archive 的 symlink
* root-level `runs` symlink 不是当前运行策略；若未来要统一 namespace，必须作为单独项目审批，并先保证现有本地 `runs/<run_id>` 兼容路径不被破坏

参见：`docs/governance/local_storage_policy.md`

---

## 3. Script Contract

每个 `scripts/Sx_*.py` 必须满足：

1. 单一入口：`if __name__ == "__main__": main()`
2. 顶部 docstring 必写：Inputs / Outputs / Frozen constants / Attrition rules
3. 统一 CLI：至少包含 `--project-root`, `--run-id`, `--seed`（默认 619）
4. 输出三件套：
   * `run_manifest.json`
   * `audit_assertions.json`
   * `manifest.json`
5. 所有输出只能写入 `runs/<run_id>/<stage>/`
6. 多线程只用于计算加速，不得改变 seed、排序、subsample 或分母口径

---

## 4. Canonical Output Principles

* 所有 plot-ready / synthesis-ready 表优先使用 long-format
* 所有 summary 表必须携带分母字段
* per-query 表用于审计证据，不强制 long，但字段必须完整
* representation not applicable 必须作为 scope policy 明示，不能伪装成 attrition
* corrected Task2 v1 的 core metrics 必须按 `dataset` 与 `cell_line` 分层，禁止 raw-pool LINCS 与 scPerturb

---

## 5. Run Preservation

**原则**：脚本内部不得自动删除 runs；由显式 cleanup 脚本执行。

* registry：`runs/_registry/story_runs.csv`
* 字段：`run_id, stage, status, created_at, notes`
* legacy/interim Task2 runs 必须保留并可标记为 `superseded`
* corrected successor runs 另行登记，不覆盖 legacy 历史记录
* 若 authoritative runs 已迁移到 NAS，只要 `runs/<run_id>` 本地兼容路径仍成立，就视为 preservation contract 满足

---

## 6. Step-level Assertion Summary

* S0：输入文件与关键列存在
* S1：LOO strict_recompute；chance identity <=1e-12；分母守恒
* S2：cross alignment contract；eligibility gate；chance identity <=1e-12
* S3：LINCS + scPerturb 双源 snapshot readiness；representation-scope policy 明确；pairing 证据可审计
* S4：按 `dataset, cell_line, target_token` 计算；禁止 raw-pool；`edist` 不得作为 cross-representation rank
* S5：按 `dataset, cell_line, direction` 计算；C2G 使用 genetic target centroids，G2C 使用 chemical target centroids；G2C `N_gallery` 计数的是 validity filtering 后的 chemical target centroids；不再依赖 chemical-context reconstruction；C2G 与 G2C 不得合并；chance identity <=1e-12
* S6：仅消费审计通过的 S4/S5 输出；`task2_benchmark_summary_long.csv` 只能从 leaderboard 层构造
* S7：只做跨任务 / 项目级 synthesis；不回写 Task2 core metrics
