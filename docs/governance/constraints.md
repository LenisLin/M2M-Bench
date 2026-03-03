# Data Constraints（Task1 已冻结可用；Task2 需从 Task1-internal + K562 raw 重建）

## 0. 状态与范围（冻结口径）

* **Task1 输入数据已冻结且可用**：`data/task1_snapshot_v1/` 为唯一来源（本轮不修改其内容）。
* **Task2 输入数据尚未准备好**：本文件仅规定 Task2 应如何从：

  1. **Task1-internal 的 LINCS 全量数据**派生出 Task2-LINCS 子集；
  2. **scPerturb 的 K562 evaluation set（外部路径）**重建 Task2-K562（含 gene/pathway/FM）。
* Task2 的最终数据应迁移到：`data/task2_snapshot_v1/`（外部路径只作为临时 source，不允许长期依赖）。

---

## 1. Task1 snapshot v1（完整定义，不得省略）

### 1.1 Snapshot layout（必须存在）

根目录：`data/task1_snapshot_v1/`

* `lincs/`
* `scperturb_delta/`
* `fm_delta/`
* `pathway/`
* `cross_contract/`

### 1.2 统一标准化（Normalization）

所有表中涉及的字符串字段统一规则（加载时执行）：

* 缺失值：`NaN -> "NA"`
* `cell_line/target`：`strip()`；建议统一大小写（例如 upper），但必须在**整个项目一致**
* `perturbation_type`：只允许 `"Genetic"` 或 `"Chemical"`
* `dataset`：只允许 `"LINCS"` 或 `"scPerturb"`

---

### 1.3 LINCS（dataset=LINCS）

#### 1.3.1 LINCS meta（CSV）

路径（示例）：`data/task1_snapshot_v1/lincs/lincs-engine1-meta.csv`

**必需列**

* `pert_type`（用于映射 perturbation_type）
* `target`
* `cell_line`
* `dose_val`（允许缺失）
* `time_val`（允许缺失）

**canonical 行键**

* `global_idx_lincs = row_number (0..N-1)`
* `canonical_query_uid_str = str(global_idx_lincs)`

#### 1.3.2 LINCS gene delta（PT）

路径（示例）：`data/task1_snapshot_v1/lincs/lincs-engine1-gene-delta.pt`

**必需结构**

* torch load 为 dict，包含 `y_delta_gene`
* `y_delta_gene` 为 2D 矩阵 `[N_rows, 2477]`

**一致性约束（fail-fast）**

* `N_rows == len(lincs-engine1-meta.csv)`

**表示范围约束（数据层）**

* LINCS 不存在 FM 表示（任何 LINCS 表/文件出现 `FM:` 视为输入异常）。

---

### 1.4 scPerturb（dataset=scPerturb）

#### 1.4.1 scPerturb meta（两份）

* Genetic：`scperturb-crispr-delta-meta.csv`
* Chemical：`scperturb-drug-delta-meta.csv`

**必需列**

* `delta_row_idx`（int；canonical 行键）
* `target_std`
* `cell_std`
* `dose_val`（允许缺失）
* `time_val`（允许缺失）

**一致性约束（fail-fast）**

* `delta_row_idx` 唯一
* `delta_row_idx` 必须在对应 delta 矩阵行范围内

**canonical 行键**

* `canonical_query_uid_str = str(delta_row_idx)`

#### 1.4.2 scPerturb gene delta（NPY）

* Genetic gene：`scperturb-crispr-gene-delta.npy`
* Chemical gene：`scperturb-drug-gene-delta.npy`

**约束**

* 2D `[N_rows, 2477]`
* 所有 `delta_row_idx ∈ [0, N_rows)`

#### 1.4.3 scPerturb pathway delta（NPY）

* Genetic pathway：`scperturb-crispr-pathway-delta.npy`
* Chemical pathway：`scperturb-drug-pathway-delta.npy`

**约束**

* 2D `[N_rows, D_pathway]`
* 与 gene delta 同一 `delta_row_idx` 行语义一致（同一实例在两种表示下行键一致）

---

### 1.5 Pathway 投影资产（W + policy）

* `data/task1_snapshot_v1/pathway/hallmark-w-2477x50.npy`
* `data/task1_snapshot_v1/pathway/lincs-pathway-policy.json`

**W 约束（fail-fast）**

* `W.shape == (2477, 50)`

**policy 约束（fail-fast）**

* `mode == "project_on_load"`
* policy 中记录的 `W_sha256`（或兼容字段名）必须等于 W 文件 sha256

---

### 1.6 FM delta（Task1；scPerturb-side）

目录：`data/task1_snapshot_v1/fm_delta/<model_name>/`

每个 model 至少包含：

* `<PerturbationType>_delta_cell.npy`（2D）
* `<PerturbationType>_delta_meta.csv`
* （强制）`delta_operator_policy.json`

`*_delta_meta.csv` 必需列：

* `row_in_modality_matrix`（int）
* `target_std`
* `cell_std`

**关键一致性约束（fail-fast）**

* `row_in_modality_matrix` 必须是对应 scPerturb meta 的 `delta_row_idx` 子集（同 perturbation_type 下）
* `row_in_modality_matrix` 唯一
* `delta_operator_policy.json` 必须包含：

  * `operator_type ∈ {"euclidean","spherical_normalized"}`
  * `eps > 0`
  * 与 Task1 的 `delta_operator_interface` 语义一致（阈值若可配置必须写入 policy）

**canonical 行键（FM）**

* `canonical_query_uid_str = str(row_in_modality_matrix)`

---

### 1.7 Cross contract（Task1-cross 输入）

路径（示例）：`data/task1_snapshot_v1/cross_contract/...`

**必需列**

* `global_idx_lincs`
* `sc_delta_row_idx`
* `expected_target`
* `expected_cell_line`

**一致性约束（fail-fast）**

* expected fields 与两侧 meta 完全一致
* 索引不越界

---

## 2. Task2 数据（当前未就绪；必须从以下来源重建）

> 重要：Task2 的 LINCS 子集必须从 **Task1-internal 的 LINCS 全量**派生，而不是从 cross 子集派生（否则会丢实例，影响机制一致性覆盖）。

### 2.1 Task2-LINCS（从 Task1-internal LINCS 派生）

#### 2.1.1 Task2-LINCS 准入规则（机制匹配）

在 LINCS 内，定义机制匹配键：

* `task2_key = (cell_line, target)`

准入条件（必须同时满足）：

* 对同一 `task2_key`，LINCS 中 **存在至少 1 个 Chemical 实例**，且 **存在至少 1 个 Genetic 实例**。

> dose/time 不进入匹配键；仅作为 instance 元数据保留。

#### 2.1.2 Task2-LINCS 派生清单（Task2 snapshot 的入口表）

目标路径：`data/task2_snapshot_v1/lincs/task2_lincs_pairs.csv`

建议字段（最小可审计）：

* `cell_line`, `target`
* `n_chemical`, `n_genetic`
* `chemical_global_idx_list`（JSON string 或外联展开表）
* `genetic_global_idx_list`

一致性约束（fail-fast）：

* 所有 idx 在 `[0, N_lincs)` 内
* `cell_line/target` 不允许为 `"NA"`（除非你明确允许；建议默认禁止）

覆盖统计（必须记录到 meta，不 fail-fast）：

* `n_task2_keys_total`
* `n_chemical`/`n_genetic` 的分布摘要（min/median/max）

---

### 2.2 Task2-scPerturb（仅 K562 子集；当前外部路径）

#### 2.2.1 当前可用文件（外部 evaluation set）

路径：`/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562`

你已列出该目录包含：

* `CRISPR_counts.pt`
* `CRISPR_meta.csv`
* `Drug_counts.pt`
* `Drug_meta.csv`
* `Common_Targets_K562.csv`
* `shared_var_names.csv`

#### 2.2.2 上游 raw 数据源（用于重算 delta 与 FM）

raw/processed 源路径：

* `/mnt/NAS_21T/ProjectData/OSMOSIS/processed/scPerturb_Processed`

**原因（你已指出）**：现有处理好的 delta/聚合版本缺少 (treated, control) pairing 证据；若仅基于其输出重算 FM embedding/delta，可能与 gene/pathway delta 的 pairing 不一致，导致 Task2 机制一致性结论不可审计。

**因此，Task2-K562 采用策略：从 raw/processed 源重算 gene/pathway/FM 的 embedding 与 delta，并在 snapshot 中持久化 pairing 证据（见 2.3）。**

---

### 2.3 Task2-K562：delta 构造与 pairing 证据（必须可审计）

你给的 delta 构造逻辑：对每个 treated cell，选择若干 control cells（同 dataset、同 cell_line、同 target），并计算 delta（gene/pathway/FM 都需一致）。

**数据约束要求：Task2 snapshot 必须持久化以下至少之一（二选一）**：

A) **显式 pairing 列表（最强可审计）**

* `treat_cell_id -> [control_cell_id_1..k]`

或 B) **可复算的 control 汇总证据（更省空间）**

* `control_sums`（按 treated cell 对齐）
* `control_counts`（按 treated cell 对齐）
* 以及 `control_selection_policy.json`（记录 control 选择过滤规则 + GLOBAL_SEED=619 + k/采样策略）

> FM delta 计算在 Task2 仍沿用 Task1 的 `delta_operator_interface`，因此 FM 的 `delta_operator_policy.json` 也必须存在（operator_type/eps 等）。

---

### 2.4 Task2 snapshot 目标目录结构（同意迁移）

目标根：`data/task2_snapshot_v1/`

建议结构：

* `lincs/`

  * `task2_lincs_pairs.csv`
  * （可选但推荐）`task2_lincs_instances.parquet`（展开到 instance 表，便于复用 Task1 管线）
* `k562/`

  * `CRISPR_counts.pt`, `CRISPR_meta.csv`
  * `Drug_counts.pt`, `Drug_meta.csv`
  * `Common_Targets_K562.csv`, `shared_var_names.csv`
  * `raw_source_manifest.json`（记录上游 raw 路径与版本/文件清单）
  * `derived/`

    * `gene_delta.npy`, `pathway_delta.npy`
    * `delta_meta.csv`（统一 schema：cell_id、perturbation_type、cell_line、target、dose/time、row_id 等）
  * `fm/<model_name>/`

    * `treated_embeddings.npy`（或等价）
    * `control_sums.npy`, `control_counts.npy`
    * `fm_delta.npy`
    * `fm_delta_meta.csv`
    * `delta_operator_policy.json`
    * `control_selection_policy.json`（若采用方案 B）

---

## 3. 数据使用小结（给 CODEX 的“防迷路”矩阵）

| Task / Scope                   | dataset           | 使用的输入数据来源                                                   | representation                     |
| ------------------------------ | ----------------- | ----------------------------------------------------------- | ---------------------------------- |
| Task1-internal                 | LINCS             | `data/task1_snapshot_v1/lincs/*` + `pathway/W`              | Gene, Pathway                      |
| Task1-internal                 | scPerturb         | `data/task1_snapshot_v1/scperturb_delta/*` + `fm_delta/*`   | Gene, Pathway, FM:*                |
| Task1-cross                    | LINCS ↔ scPerturb | `task1_snapshot_v1/cross_contract/*` + 两侧 meta+delta        | Gene, Pathway（FM 不进入 cross，除非另行冻结） |
| Task2 (mechanism)              | LINCS             | **从 Task1-internal LINCS 全量派生**：`task2_lincs_pairs.csv`     | Gene, Pathway                      |
| Task2 (mechanism, K562 subset) | scPerturb-K562    | 外部源 → 迁移到 `data/task2_snapshot_v1/k562/`，并从 raw 重算 delta/FM | Gene, Pathway, FM:*（K562 子集）       |

---

## 4. Task2 就绪判定（用于实现层 fail-fast）

Task2 代码入口应先检查：

* `data/task2_snapshot_v1/` 存在且包含：

  * `lincs/task2_lincs_pairs.csv`
  * `k562/` 的基础文件（counts/meta/common targets/var names）
  * `k562/derived/`（至少 gene/pathway delta 与 meta）**或**能够明确触发“重算管线”（由 scripts 驱动）
  * 如要做 FM：`k562/fm/<model>/` 必须包含 operator policy +（pairing 证据 A 或 B）

若缺失：直接报 `TASK2_SNAPSHOT_NOT_READY`，禁止隐式回退到 Task1 snapshot。
