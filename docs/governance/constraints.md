# Data Constraints（Task1 已冻结可用；Task2 corrected multisource successor 已物化）

## 0. 状态与范围（冻结口径）

* **Task1 输入数据已冻结且可用**：`data/task1_snapshot_v1/` 为唯一来源（本轮不修改其内容）。
* **Task2 当前已物化的本地状态仅为 legacy/interim scPerturb-K562 子集**：`data/task2_snapshot_v1/k562/` 保留为历史证据，不可被重写为 corrected final Task2。
* **Corrected Task2 contract 是 multisource**：必须同时纳入：
  1. **Task1-internal 的 LINCS 全量数据**派生出的 Task2-LINCS 子集；
  2. **复用已审计 legacy K562 证据并重定址到 corrected successor root 的** Task2-scPerturb-K562 子集（含 gene/pathway/FM）。
* Corrected Task2 successor snapshot 必须使用独立版本根目录，而不是静默覆写 `data/task2_snapshot_v1/`。当前冻结根目录为 `data/task2_snapshot_v2/`。

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

## 2. Task2 数据（corrected multisource successor contract）

> 重要：Task2 的 LINCS 子集必须从 **Task1-internal 的 LINCS 全量**派生，而不是从 cross 子集派生（否则会丢实例，影响机制一致性覆盖）。

> 重要：`dataset, cell_line, target_token` 是 Task2 的 **analysis/cohort key**，不是 snapshot 行身份键。诸如 `global_idx_lincs`、`row_id`、`treated_cell_id`、`query_uid` 等属于 row identity / audit identity。

### 2.1 Task2-LINCS（从 Task1-internal LINCS 派生）

#### 2.1.1 Task2-LINCS 准入规则（机制匹配）

在 LINCS 内，定义 Task2 analysis/cohort key：

* `mech_key = (dataset="LINCS", cell_line, target_token)`

准入条件（必须同时满足）：

* 对同一 `mech_key`，LINCS 中 **存在至少 1 个 Chemical 实例**；
* 对同一 `mech_key`，LINCS 中 **存在至少 1 个 Genetic 实例**。

> dose/time 不进入匹配键；仅作为 instance 元数据保留。

**row identity（LINCS）**

* `global_idx_lincs` 仍然是实例层 canonical 行键
* 一个 Chemical 行可通过 `target_tokens` 贡献给多个 `target_token` cohort，但不复制行身份

#### 2.1.2 Task2-LINCS 派生清单（Task2 snapshot 的入口表）

目标路径：`data/task2_snapshot_v2/lincs/task2_lincs_pairs.csv`

建议字段（最小可审计）：

* `dataset="LINCS"`, `cell_line`, `target_token`
* `n_chemical`, `n_genetic`
* `chemical_global_idx_list`（JSON string 或外联展开表）
* `genetic_global_idx_list`

一致性约束（fail-fast）：

* 所有 idx 在 `[0, N_lincs)` 内
* `cell_line/target_token` 不允许为 `"NA"`（除非你明确允许；建议默认禁止）

覆盖统计（必须记录到 meta，不 fail-fast）：

* `n_task2_keys_total`
* `n_chemical`/`n_genetic` 的分布摘要（min/median/max）

---

### 2.2 Task2-scPerturb（corrected v1 仅 K562 子集）

#### 2.2.1 当前可用文件（外部 evaluation set）

路径：`/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562`

你已列出该目录包含：

* `CRISPR_counts.pt`
* `CRISPR_meta.csv`
* `Drug_counts.pt`
* `Drug_meta.csv`
* `Common_Targets_K562.csv`
* `shared_var_names.csv`

#### 2.2.2 当前物化策略（冻结）

当前 corrected multisource Task2 的 scPerturb K562 子树**不是**从 raw 重新计算；
它通过复用已审计的 legacy K562 证据并重定址到 corrected successor root
来物化：

* source: `data/task2_snapshot_v1/k562/`
* destination: `data/task2_snapshot_v2/scperturb_k562/`

复用内容必须保留：

* 基础输入文件（counts/meta/common targets/shared var names）
* `derived/pair_list.parquet`
* `derived/delta_meta.csv`
* `derived/gene_delta.npy`
* `derived/pathway_delta.npy`
* `fm/<model>/fm_delta.npy`
* `fm/<model>/fm_delta_meta.csv`
* `fm/<model>/delta_operator_policy.json`

authoritative evidence:

* `runs/0310_fix_1hae/s3_build_task2_multisource_snapshot/run_manifest.json::scperturb_reuse_summary`
* `data/task2_snapshot_v2/snapshot_manifest.json::datasets.scPerturb.subtree`

**analysis/cohort key（scPerturb K562）**

* `mech_key = (dataset="scPerturb", cell_line="K562", target_token)`

**row identity（scPerturb K562）**

* `row_id`, `treated_cell_id`, `query_uid` 等保持实例层身份
* Chemical 行可通过 `target_tokens` 参与多个 `target_token` cohort，但不复制 snapshot 行身份

---

### 2.3 Task2-K562：delta 构造与 pairing 证据（必须可审计）

Task2-K562 的 delta 构造逻辑：对每个 treated cell，选择若干 control cells，并计算 delta（gene/pathway/FM 都需一致）。

控制组冻结规则：

* CRISPR：all CRISPR controls
* Drug：Drug controls with matched `time`
* controls 不按 target 匹配；controls 是 untargeted

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

目标根：`data/task2_snapshot_v2/`

建议结构：

* `lincs/`

  * `task2_lincs_pairs.csv`
  * （可选但推荐）`task2_lincs_instances.parquet`（展开到 instance 表，便于复用 Task1 管线）
* `scperturb_k562/`

  * `CRISPR_counts.pt`, `CRISPR_meta.csv`
  * `Drug_counts.pt`, `Drug_meta.csv`
  * `Common_Targets_K562.csv`, `shared_var_names.csv`
  * `derived/`

    * `pair_list.parquet`
    * `gene_delta.npy`, `pathway_delta.npy`
    * `delta_meta.csv`（统一 schema：cell_id、perturbation_type、cell_line、target、dose/time、row_id 等）
  * `fm/<model_name>/`

    * `fm_delta.npy`
    * `fm_delta_meta.csv`
    * `delta_operator_policy.json`

---

### 2.5 Representation scope policy（冻结）

* LINCS：只允许 `Gene`, `Pathway`
* scPerturb K562：允许 `Gene`, `Pathway`, `scgpt`, `geneformer`, `scbert`, `scfoundation`, `uce`, `state`, `tahoe-x1`
* 不支持的表示必须视为 scope policy，而不是 attrition

---

## 3. 数据使用小结（给 CODEX 的“防迷路”矩阵）

| Task / Scope                   | dataset           | 使用的输入数据来源                                                   | representation                     |
| ------------------------------ | ----------------- | ----------------------------------------------------------- | ---------------------------------- |
| Task1-internal                 | LINCS             | `data/task1_snapshot_v1/lincs/*` + `pathway/W`              | Gene, Pathway                      |
| Task1-internal                 | scPerturb         | `data/task1_snapshot_v1/scperturb_delta/*` + `fm_delta/*`   | Gene, Pathway, FM:*                |
| Task1-cross                    | LINCS ↔ scPerturb | `task1_snapshot_v1/cross_contract/*` + 两侧 meta+delta        | Gene, Pathway（FM 不进入 cross，除非另行冻结） |
| Task2 (mechanism, corrected)   | LINCS             | **从 Task1-internal LINCS 全量派生**：`task2_lincs_pairs.csv`     | Gene, Pathway                      |
| Task2 (mechanism, corrected)   | scPerturb-K562    | 复用已审计 legacy K562 证据并物化到 `data/task2_snapshot_v2/scperturb_k562/` | Gene, Pathway, scgpt, geneformer, scbert, scfoundation, uce, state, tahoe-x1 |
| Task2 (legacy/interim)         | scPerturb-K562    | `data/task2_snapshot_v1/k562/`                                    | current audited historical subset  |

---

## 4. Task2 就绪判定（用于实现层 fail-fast）

Corrected Task2 代码入口应先检查：

* `data/task2_snapshot_v2/` 存在且包含：
  * `lincs/task2_lincs_pairs.csv`
  * `scperturb_k562/` 的基础文件（counts/meta/common targets/var names）
  * `scperturb_k562/derived/`（至少 pair_list、gene/pathway delta 与 meta）
  * 如要做 FM：`scperturb_k562/fm/<model>/` 必须包含 operator policy + pairing 证据

额外约束：

* 所有 Task2 analysis-facing 数据必须可还原到 `dataset, cell_line, target_token`
* LINCS 不允许出现 FM 表示
* corrected Task2 v1 core metrics 必须按 `dataset` 与 `cell_line` 分层，禁止 raw-pool LINCS 与 scPerturb

当前 `data/task2_snapshot_v1/` 只满足 legacy/interim scPerturb-K562 历史证据口径，不满足 corrected final Task2 readiness。

若 corrected successor 缺失：直接报 `TASK2_SNAPSHOT_NOT_READY`，禁止隐式回退到 Task1 snapshot 或 legacy `task2_snapshot_v1/` 充当 final Task2。
