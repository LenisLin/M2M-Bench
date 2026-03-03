# Task1 计划

## 0) 全局常量（冻结）

```pseudo
GLOBAL_SEED = 619
EDIST_MAX_N = 256
HIT_K_LIST = [1, 5, 10]
EDIST_WORKERS = <config/arg>
```

初始化（每个脚本入口第一件事）：

```pseudo
random.seed(GLOBAL_SEED)
np.random.seed(GLOBAL_SEED)
rng = np.random.default_rng(GLOBAL_SEED)
```

---

## 1) 统一键与基础表（唯一真源）

### 1.1 Instance 表：`task1_instances.parquet`

每行 = 一个 instance 在一个 representation 下的 delta 记录（只记录索引与 meta，不强制把向量塞进表）：

**必备列**

* `dataset ∈ {LINCS, scPerturb}`
* `perturbation_type ∈ {Genetic, Chemical}`
* `cell_line`, `target`
* `representation ∈ {Gene, Pathway, FM:*}`
* `canonical_query_uid`（str）

  * LINCS = `global_idx_lincs`
  * scPerturb Gene/Pathway = `sc_delta_row_idx`
  * FM = `row_in_modality_matrix`
* `delta_row_index`（用于 load delta 向量的行号）
* `delta_valid_bool`（FM delta_operator valid_mask；Gene/Pathway 恒 True）
* `dose_val`, `time_val`（可选；只记录）

**分组键**

```pseudo
group_key = (dataset, perturbation_type, cell_line, target)
```

### 1.2 Cross 配对表：`task1_cross_pairs_validated.csv`

从 frozen cross contract 构建并校验：

**输入 contract 至少含**

* `global_idx_lincs`, `sc_delta_row_idx`, `expected_target`, `expected_cell_line`

**校验（fail-fast）**

```pseudo
for row in cross_contract:
  assert lincs_meta[target,cell_line] == expected_*
  assert sc_meta[target,cell_line]    == expected_*
```

**输出列**

* `perturbation_type`（来自 meta；两侧必须一致，否则 fail-fast）
* `cell_line`, `target`
* `global_idx_lincs`, `sc_delta_row_idx`

### 1.3 Cross 规则化筛选：`cross_attrition.csv`

```pseudo
for pert_type in {Chemical, Genetic}:
  n_match_key = count_unique( (cell_line,target,pert_type) in validated_pairs )
  if n_match_key < 5:
     mark SKIP (pert_type)
     write cross_attrition row: excluded_reason="matched_keys_lt5"
```

> 规则：SKIP 的模态在 cross 下**不计算任何结果表行**（或输出 NA 行但必须显式标注 excluded）。

---

## 2) 预生成索引（单线程；保证多线程可复现）

所有涉及抽样/切分的模块，统一先生成索引表（**只调用一次 rng**，且遍历顺序固定）：

### 2.1 固定遍历顺序

```pseudo
groups = sorted(unique(group_key))  # 稳定排序
```

### 2.2 split-half 索引：`split_half_index.parquet`

对每个 `(group_key, representation)`：

```pseudo
idx = instance_indices_of(group_key, representation) where delta_valid_bool=True
n = len(idx)

if n < 4:
  record split_half_status = "NA_n_lt_4"
  continue

# 固定：先对 idx 按 canonical_query_uid 排序，避免输入顺序漂移
idx_sorted = sort_by(canonical_query_uid, idx)

perm = rng.permutation(n)          # 仅在这里消耗 rng
half = n // 2
A = idx_sorted[perm[:half]]
B = idx_sorted[perm[half:]]
record A_indices, B_indices
```

### 2.3 E-distance 抽样索引：`edist_subsample_index.parquet`

对每个需要 energy distance 的集合（internal 的 A/B；cross 的 L/S）：

```pseudo
subsample(X_indices):
  if len(X_indices) <= EDIST_MAX_N: return X_indices
  # 固定：先排序，再抽样
  X_sorted = sort_by(canonical_query_uid, X_indices)
  perm = rng.permutation(len(X_sorted))
  return X_sorted[perm[:EDIST_MAX_N]]
```

---

## 3) Group-level 指标计算（cosine / PCC / Energy distance）

### 3.1 Internal（split-half）

输出：`task1_group_internal.parquet`

对每个 `(group_key, representation)`：

```pseudo
if split_half_status == "NA_n_lt_4":
  output metrics = NA
  write group_attrition row: excluded_reason="split_half_requires_n>=4"
  continue

A_idx, B_idx = load split_half_index
A_vec = load_delta_vectors(A_idx)
B_vec = load_delta_vectors(B_idx)

centroid_A = mean(A_vec)
centroid_B = mean(B_vec)

cosine = dot(centroid_A, centroid_B)/(||A||*||B||)
pcc    = pearson_corr(centroid_A, centroid_B)  # always compute; note in report

A_sub = subsample(A_idx)  # from edist_subsample_index
B_sub = subsample(B_idx)
if len(A_sub) < 2 or len(B_sub) < 2:
   edist = NA
   write group_attrition row: excluded_reason="edist_requires_n>=2_per_split"
else:
   edist = energy_distance_biascorr(load_vec(A_sub), load_vec(B_sub))  #见 3.3

write output row with: n_total, n_A, n_B, n_A_sub, n_B_sub, cosine, pcc, edist
```

### 3.2 Cross（matched）

输出：`task1_group_cross.parquet`

对每个 `perturbation_type`：

```pseudo
if SKIP by cross_attrition: continue

# 构建 matched 子集分布
LINCS_instances  = instances where dataset=LINCS and canonical_query_uid in validated_pairs.global_idx_lincs
scP_instances    = instances where dataset=scPerturb and canonical_query_uid in validated_pairs.sc_delta_row_idx

# 按 match_key=(pert_type,cell_line,target) 聚合（dataset 不同但 match_key 相同）
for each match_key:
   L_idx = indices of LINCS instances for this match_key, delta_valid_bool=True
   S_idx = indices of scP instances   for this match_key, delta_valid_bool=True

   if len(L_idx)==0 or len(S_idx)==0: continue  # 或输出 NA 并记录 reason

   centroid_L = mean(load_vec(L_idx))
   centroid_S = mean(load_vec(S_idx))

   cosine = cosine(centroid_L, centroid_S)
   pcc    = pearson_corr(centroid_L, centroid_S)

   L_sub = subsample(L_idx)  # EDIST_MAX_N=256
   S_sub = subsample(S_idx)
   if len(L_sub)<2 or len(S_sub)<2:
       edist = NA; record reason
   else:
       edist = energy_distance_biascorr(load_vec(L_sub), load_vec(S_sub))

   write output row with: n_L,n_S,n_L_sub,n_S_sub, cosine,pcc,edist
```

### 3.3 Energy distance（bias-corrected；squared Euclidean 距离）

输入 `X ∈ R^{N×d}`, `Y ∈ R^{M×d}`（已 subsample，N,M ≤ 256）：

```pseudo
sqdist(a,b) = sum_k (a_k - b_k)^2

deltaXY = (1/(N*M)) * sum_i sum_j sqdist(x_i, y_j)

sigmaX  = (1/(N*(N-1))) * sum_{i!=j} sqdist(x_i, x_j)
sigmaY  = (1/(M*(M-1))) * sum_{i!=j} sqdist(y_i, y_j)

E = 2*deltaXY - sigmaX - sigmaY
return E
```

并行要求：对不同 group 的 edist 计算可多线程/多进程；但索引已预生成，因此结果与调度无关。

---

## 4) Retrieval（euclid + LOO + chance correction；输出 Hit@1/5/10 + MRR）

输出：

* `task1_retrieval_per_query.parquet`
* `task1_retrieval_summary.csv`
* `chance_identity_check.csv`

### 4.1 Cohort 定义

```pseudo
cohort_key = (data_scope, dataset, perturbation_type, direction, representation)
```

internal：

* `(internal, LINCS, Chemical, internal, Gene/Pathway)`
* `(internal, LINCS, Genetic, internal, Gene/Pathway)`
* `(internal, scPerturb, Chemical, internal, Gene/Pathway/FM:*)`
* `(internal, scPerturb, Genetic, internal, Gene/Pathway/FM:*)`

cross：对每个未 SKIP 的 `perturbation_type` 与方向（如有）分别作为 cohort。

### 4.2 internal centroid + LOO

```pseudo
for cohort:
  instances = filter(task1_instances by cohort_key and delta_valid_bool=True)
  build group sums/counts:
     sum_g, count_g for each group_key

  centroids = sum_g / count_g

  for each query q in instances:
     g = q.group_key
     if count_g <= 1:
        exclude query with reason="loo_requires_count>1"
        continue

     true_centroid = (sum_g - q.vec)/(count_g - 1)   # strict LOO
     candidate_centroids = centroids (optionally capped by N_gallery_max deterministic rule)
     ensure true group centroid is included (if cap removes it, replace last with true)
```

### 4.3 Score / rank / metrics（euclid）

统一 score：

```pseudo
score(q,c) = -||q - c||^2
```

rank 与 raw 指标：

```pseudo
rank_true = 1 + count{ score(q,c_i) > score(q,true_centroid) }
MRR_raw = 1/rank_true
HitK_raw(K) = 1[rank_true <= K]  for K in {1,5,10}
N_gallery = number of candidate centroids
```

### 4.4 chance correction（exact；m_pos=1）

```pseudo
expected_mrr = E[1/R] under random ranking with 1 positive among N_gallery
expected_hitK(K) = P(R<=K) under random ranking with 1 positive among N_gallery

MRR_corr = MRR_raw - expected_mrr
HitK_corr = HitK_raw - expected_hitK
```

**summary identity check（容差）**

```pseudo
delta = mean(MRR_corr) - (mean(MRR_raw) - mean(expected_mrr))
assert |delta| <= 1e-12
```

---

## 5) 显著性（必须产物；Gene baseline）

输出：`task1_significance.csv`

规则：所有比较以 Gene 为 baseline（同 scope、同 dataset/pert_type、同 metric）。

### 5.1 group-level 显著性

对 `metric ∈ {cosine, pcc, edist}`：

```pseudo
paired_support = groups with metric_Gene and metric_rep both finite
effect_size = mean(metric_rep - metric_Gene) on paired_support
p_value = wilcoxon_signed_rank(metric_rep, metric_Gene)  # paired
q_value = BH_FDR within family
```

### 5.2 retrieval 显著性

对 `metric ∈ {mrr_corr, hit1_corr, hit5_corr, hit10_corr}`：

```pseudo
paired_support = intersection(canonical_query_uid) between Gene and rep
effect_size = mean(metric_rep - metric_Gene)
p_value = wilcoxon_signed_rank(...)
q_value = BH_FDR within family
```

family（冻结建议，写死即可）：

```pseudo
family_key = (scope, dataset_or_pair, perturbation_type, metric_name)
```

---

## 6) 富集偏好（必须产物：cell_line 与 target）

输出：

* `task1_enrichment_cell_line.csv`
* `task1_enrichment_target.csv`

positives 定义（冻结）：

```pseudo
success_bool = (MRR_corr > 0) AND (Hit10_corr > 0)
```

对每个 cohort（scope×dataset×pert_type×representation）与 entity_type：

```pseudo
for entity in unique(entity_value):
   2x2 table:
     a = #pos in entity
     b = #neg in entity
     c = #pos not in entity
     d = #neg not in entity
   p_value = fisher_exact(a,b,c,d)
   effect_size = odds_ratio(a,b,c,d)  # or risk_diff; choose one and freeze
q_value = BH_FDR within family=(cohort_key, entity_type)
```

---

## 7) Attrition / NA 规则（统一）

输出：`group_attrition.csv`, `cross_attrition.csv`, `retrieval_attrition.csv`

* internal group-level：`n<4` → NA + `split_half_requires_n>=4`
* edist：任一侧 subsample <2 → NA + `edist_requires_n>=2`
* cross：`n_match_key<5` → SKIP + `matched_keys_lt5`
* retrieval：`count_g<=1` → query excluded + `loo_requires_count>1`
* FM：`delta_valid_bool=False` → excluded + `invalid_delta`
