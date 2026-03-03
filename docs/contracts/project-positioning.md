# M2M-Bench 项目总览

## 1. 项目目标与总体框架

M2M-Bench 的目标是建立一个**可审计、可复现、可比较**的 benchmark 框架，用于评估扰动响应（perturbation response）在不同数据源与表示空间下的**一致性（concordance）**与**可比性（comparability）**。项目输出以“事实性比较结论”为主：在明确协议与分母透明的前提下，不同表示在 Task1（模态一致性）与 Task2（机制一致性）下的 concordance 证据。

框架由两个任务组成：

* **Task1：模态一致性（Modality Concordance）**
  比较时保持 perturbation_type 一致（Chemical↔Chemical、Genetic↔Genetic），包含两个子项：

  1. **Internal**：同一数据集内部一致性（LINCS 内；scPerturb 内）
  2. **Cross**：跨数据集一致性（LINCS ↔ scPerturb 的 matched 子集），但仍保持 perturbation_type 一致

* **Task2：机制一致性（Mechanism Concordance）**
  在同一数据集内部比较 **Genetic vs Chemical** 的一致性，评测协议与 Task1 复用；此外 Task2 还包含一个受限子集：**K562 subset**，在该子集上允许使用 FM 表示进行机制一致性评估（FM 主要面向单细胞，因此只在 scPerturb 侧定义；具体约束写入 Task2 细则）。

---

## 2. 数据与对象抽象（保持）

### 2.1 Instance（单样本）最小字段

每条扰动样本记录：

* `dataset ∈ {LINCS, scPerturb}`
* `perturbation_type ∈ {Genetic, Chemical}`
* `target`
* `cell_line`
* `delta_vector`（某表示空间中的扰动响应向量）
* Chemical 可选协变量：`dose_val, time_val`（仅 instance 记录，不进入 centroid）

### 2.2 Group（分组键）

* `group_key = (dataset, perturbation_type, cell_line, target)`

---

## 3. 任务定义与预期结论产物（FM 与 matched 说明更新）

### 3.1 表示空间（Representation）

* `Gene`：基因空间 delta（线性向量空间，可直接做差/聚合）
* `Pathway`：`Gene @ W`（线性投影空间）
* `FM:*`：foundation model embedding delta

  * **Task1：仅在 scPerturb internal 使用**（这些 FM 面向单细胞数据开发）
  * **Task2：在 K562 subset（受限子集）允许使用 FM**（用于 Genetic vs Chemical 的机制一致性评估）

### 3.2 Matched 的语义（Cross 与 Task2）

* **Task1-cross matched**：在跨数据集（LINCS ↔ scPerturb）中，选择在 `(perturbation_type, cell_line, target)` 上可对应的实例集合，并通过冻结对齐契约映射到具体行索引（细则在 Task1 伪代码 spec 中写死）。
* **Task2 matched**：在同一数据集内按 `(dataset, cell_line, target)` 匹配 Chemical 与 Genetic 两个子集，形成机制对比对。

---

## 4. 评测协议（强调：Task1-internal / Task1-cross / Task2 复用同一套方法）

**关键点：**

* 三个 scope（Task1-internal、Task1-cross、Task2）使用**同一套分析“模块”**：

  1. group-level concordance（主分析）
  2. instance-level retrieval（辅助视角）
  3. 显著性与偏好富集（解释性必产物）
* 差异只来自：数据域（internal/cross）、匹配集合（cross contract 或机制匹配）、以及允许的 representation 集合（FM 仅限特定 scope）。

### 4.1 Group-level Concordance（主分析，所有 metric 在此层完成）

在每个表示空间的 delta 上，构建 group-level 对象并计算一致性指标：

1. **group 对象构建**

* 对每个 `group_key` 收集该组所有 instances 的 `delta_vector`
* 构建：

  * `centroid`（聚合向量；聚合函数冻结于任务细则）
  * `distribution`（该 group 的 instance delta 集合）

2. **Task1-internal（同数据集内部一致性：采用 split-half）**

* 对每个 group，将其 instance 集合随机（确定性 seed）分为两半：A/B
* 分别计算 `centroid_A` 与 `centroid_B`，并在 group-level 上计算 concordance：

  * `cosine(centroid_A, centroid_B)`
  * `PCC(centroid_A, centroid_B)`
  * `E-distance(distribution_A, distribution_B)`（energy distance，bias-corrected 默认）
* 在 cohort 层汇总这些 group-level concordance 的分布（均值/分位数/支持集）。

3. **Task1-cross（跨数据集 matched group 一致性）**

* 对 matched group 对（LINCS 与 scPerturb）计算：

  * `cosine(centroid_L, centroid_S)`
  * `PCC(centroid_L, centroid_S)`
  * `E-distance(distribution_L, distribution_S)`（若定义为 distribution-level；否则仅以 centroid-level 输出并在细则中说明）
* 输出每个 matched group 对的 concordance 以及 cohort-level 汇总。

4. **Task2（同数据集内 Chemical vs Genetic 机制一致性）**

* 对 `(dataset, cell_line, target)` 匹配的 Chemical/Genetic 两组，计算：

  * centroid-level cosine / PCC
  * distribution-level E-distance
* 在 K562 subset（受限）上，允许在 FM 表示空间重复同样计算（具体 scope 在 Task2 细则中冻结）。

### 4.2 Instance-level Retrieval（辅助视角：仍在 delta 空间，但不用于 cosine/PCC）

Instance-level retrieval 用于提供“可解释的实例层证据”，但**不承载 cosine/PCC/E-distance 的主要结论**：

* query 为 instance delta
* gallery 为 group centroid 集合
* 输出 `rank_true` 与检索指标：MRR、Hit@1/5/10（raw 与 chance-corrected），并记录客观的 `N_gallery`
* retrieval 的 score 规则与 metric_family 绑定（细则中定义），但 **cosine/PCC 不在 instance-level 做**（避免不必要的计算负担）。

> 这里的 retrieval metric_family 可以选择一个“轻量主力”（例如欧氏类/内积类之一），并在 Task1/Task2 细则中冻结；总览不在此处锁定实现细节。

---

## 5. 指标族与必选集合（按你要求：cosine/PCC/E-distance 只用于 group-level）

### 5.1 Cosine similarity（仅 group-level）

用于 centroid–centroid 一致性（internal split-half、cross matched、task2 chemical-vs-genetic）。

### 5.2 PCC（仅 group-level）

用于 centroid–centroid 一致性（与 cosine 并列，提供不同视角）。

### 5.3 E-distance（Energy distance，bias-corrected 默认；仅 group-level）

用于 distribution–distribution 一致性（internal split-half、cross matched、task2 chemical-vs-genetic）。
其定义遵循你给出的 energy distance 公式与 bias-corrected 版本；距离计算采用平方欧氏距离（细则中写死：标准/偏差修正版本的使用范围与极小样本处理）。

---

## 6. 解释性产物（必须产物；显著性以 Gene 为 baseline）

### 6.1 表示比较的显著性检验（必须）

* 所有显著性对比统一以 **Gene 作为 baseline**：

  * internal：Pathway vs Gene；FM vs Gene（仅适用 scope：scPerturb internal；Task2 K562 subset）
  * cross：Pathway vs Gene
  * task2：Chemical vs Genetic 的 concordance 差异（并可在不同表示空间比较，但仍以 Gene baseline 给出主结论）

* 输出必须包含：effect size、p-value、q-value（BH-FDR）、以及 N（支持集透明）。

### 6.2 cell_line / target 偏好富集（必须）

* 以 group-level 或 instance-level 的 success/high-concordance 标记定义 positives
* 对 cell_line 与 target 做富集检验，输出 effect size + p/q + 分母透明。

---

## 7. 目录结构（保持 v0.3）

* `docs/`：设计图与合同（总览、Task1/Task2 伪代码 spec、指标定义、角色与工程约束）
* `data/`：冻结输入数据（snapshot）
* `src/`：可复用实现（IO、delta 加载、metrics、group-level、retrieval、显著性、富集、audit）
* `scripts/`：与分析 step 严格对齐的薄入口（只做参数解析与调用 src；禁止无限新建）
* `runs/`：运行产物（可删可复现）

# Story Proposal

## 0. 一句话摘要

M2M-Bench 将“虚拟细胞（virtual cell）”评估中的 **模态鲁棒性（Modality Concordance）** 与 **机制可迁移性（Mechanism Concordance）** 强制解耦，并在真实药物生态（多靶点、剂量/时间）下，用 **multi-positive 结构归一化的检索指标** 与 **基准健康度（dominance/coverage）诊断**，给出可行动的模型选择与使用边界结论；可审计证据链仅作为支撑性基础设施。

---

## 1. 背景与问题：为什么现有 virtual cell 评测被混杂

virtual cell 模型通常声称具有广义泛化能力，但多数评测把两个本质不同的因素混在一起：

1. **测量/平台差异（modality/measurement shift）**：同机制扰动在不同平台、不同处理流程中是否一致？
2. **机制差异（mechanism shift）**：药物扰动与基因扰动是否可互相替代（同靶点、同细胞系下的可迁移性）？

单一 leaderboard 分数无法区分“模型只是拟合了平台/流程”还是“确有机制层面的可迁移性”，导致研究结论不可解释、不可对比、不可用于转化场景。

---

## 2. 核心观点：两轴诊断范式

我们沿两条正交轴评估 virtual cell：

* **Task1：Modality Concordance**
  *同机制类型*（chemical-chemical 或 genetic-genetic），评估同平台内部一致性（internal）与跨平台一致性（cross）。
* **Task2：Mechanism Concordance**
  *同 context*（dataset、cell_line、target_token），评估 chemical ↔ genetic 的一致性。

这提供一个诊断视角：模型在 Task1 强而 Task2 弱意味着“平台鲁棒但机制不可迁移”；反之亦然。

---

## 3. 贡献（面向生物/virtual cell 社区）

### Contribution 1 — 2D 诊断图谱（从 1D 排行榜到可诊断地图）

我们用 Task1（模态）与 Task2（机制）形成二维诊断平面，并用同一套指标/分母报告框架，使模型失败模式可被明确定位。

#### 2D Map 坐标定义（冻结）

为避免 Task1/Task2 方差结构不同导致的视觉误导，二维坐标采用 **percentile normalization**（0–1）：

* 先定义每轴的主标量（采用 **chance-corrected retrieval** 的宏平均作为统一量纲）：

  * `S1 = macro_avg(mean_mrr_corrected over Task1 families)`
  * `S2 = macro_avg(mean_mrr_corrected over Task2 families)`
    其中“family”的定义由各自 summary 表的键（scope×dataset_or_direction×perturbation_type×representation 等）给出，采用 **macro-average**（每个 family 同权），并在表中显式输出分母。
* 再映射为二维坐标：

  * `X = percentile_rank(S1 among compared reps/models)`
  * `Y = percentile_rank(S2 among compared reps/models)`
* 原始 `S1,S2` 同步在 leaderboard 表中保留，便于复核与解释。

### Contribution 2 — 真实药物生态：保留多靶点而不引入分母幻觉

LINCS 中多靶点药物广泛存在。简单删除含分号样本会造成覆盖崩溃与严重选择偏差。我们采用 **explode-matching（关系展开，不复制实例）**：

* chemical instance 具有 `target_tokens`（由 `clean_target_mapped` 分号拆分去重排序）
* 在 Task2 的 `mech_key=(dataset, cell_line, target_token)` 构建中，只要 `target_token ∈ target_tokens(instance)`，就视为该 instance 对该 mech_key 有效匹配
* **不复制** chemical instance 行（避免虚假样本膨胀）
* 保留 `specificity_tier` 与 `n_targets` 作为分层元数据，用于解释性分析

**边界声明（冻结）**：`clean_target_mapped` 与 `specificity_tier` 视为外部注释；我们不做机制因果推断，仅报告在该注释体系下的评测事实，并通过 `tier/n_targets` 分层做敏感性描述。

### Contribution 3 — multi-positive chance correction：让检索指标可比

explode-matching 与 context gallery 会自然产生 `m_pos>1`。若不处理，MRR/Hit@K 在不同 `m_pos`、不同 `N_gallery` 下不可比。我们强制：

* per-query 输出 `m_pos` 与 `N_gallery`
* 使用精确期望 `E[MRR|N,m]`、`P(hit@K|N,m)` 做校正
* 输出 identity check（容差 `1e-12`）确保实现闭环
* 空假设明确为标准 retrieval null：**随机均匀排名**

### Contribution 4 — Benchmark Health：dominance/coverage 诊断

基准结论可能被少数高频 target/药物 context 主导。我们新增健康度输出：

* topK 贡献占比曲线
* **Gini coefficient**（0–1）作为 headline dominance
* effective_n（推荐）与 entropy（可选）

**硬约束（冻结）**：dominance 权重只允许来自 **chance-corrected** 指标；连续 corrected 指标计算 Gini 前必须非负化（例如 `max(0, mrr_corrected)`）；推荐 headline 使用天然非负的 `success_bool_both_sum`。

### Contribution 5 — 可审计证据链（支撑性基础设施）

我们提供最小证据链（per-query、attrition、chance identity、pairing evidence、alignment proof）以支持上述贡献可复核，不将“可审计”作为主贡献主张。

---

## 4. 核心研究问题与对应图表（不做机制推断）

1. 数据生态：纳入了哪些数据、cell line/target 分布、polypharmacology 分布（旭日图/饼图）
2. Task1 internal：LINCS/scPerturb 内部一致性（leaderboard）
3. Task1 cross：跨模态一致性及其相对 internal 的差距（leaderboard + boxplot）
4. Task2：跨机制一致性；哪些 cell line/target 更稳定；dose/time 分布（leaderboard + lift maps + 分布表）
5. Gene/Pathway/FM：相对 Gene baseline 的 effect + 显著性 + 分层（p/q）
6. Benchmark health：是否被少数高频 target/context 主导（dominance curves + Gini/effective_n）

---

## 5. 统计与显著性（冻结、有限）

* 不做机制因果推断；仅输出评测事实与分层描述
* 以 Gene 为基线做 representation 对比
* **Paired testing（冻结）**：

  * Task1：以 `canonical_query_uid`（或等价稳定 per-query uid）做严格成对
  * Task2 group：以 `mech_key=(dataset, cell_line, target_token)` 成对
  * Task2 retrieval：同 direction 内按 query set 成对（显式输出 `n_overlap_pairs`）
* **有效样本交集（冻结）**：所有 paired test 必须在 `Gene_valid_uids ∩ Rep_valid_uids` 上进行（尤其 FM valid_mask_drop 情况）
* 多重检验：BH-FDR within 预注册 family（family key 写入结果表）

---

## 6. 输出契约（Story 层级最低要求）

* Plot-ready 长表（R 直接作图）+ Audit-ready 最小证据链（bounded audit）
* 每个 story claim 必须能在这些表上重现；其他中间文件允许存在但不替代这些 canonical outputs

---

## 7. 关键限制（显式）

* 外部注释噪声（tier/target mapping）：通过分层报告呈现敏感性，不做生物推断
* dose/time：本版先输出分布与在 G2C gallery 中显式建模；bins/regression 后置到正式实验阶段再冻结
