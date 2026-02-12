import os
import scanpy as sc
import pandas as pd
import numpy as np
import torch
import gc

# ================= 配置 =================
# 原始数据路径
DATAPATH = "/mnt/NAS_21T/ProjectData/OSMOSIS/raw/scPerturb_Processed"
# 输出根目录
BASE_OUTPUT_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen"
# 数据集名称
CRISPR_Name = "ReplogleWeissman2022_K562_essential"
drug_Name = "SrivatsanTrapnell2020_K562"

# 创建专属子目录，保持整洁
OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "Benchmark_Datasets", "Evaluation_Set_K562")
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_adata(name_):
    print(f">>> Loading {name_} ...")
    adata_ = sc.read_h5ad(os.path.join(DATAPATH, "Cleaned_" + name_ + ".h5ad"))
    obs_path = os.path.join(DATAPATH, "Cleaned_" + name_ + "_obs.csv")
    
    if os.path.exists(obs_path):
        obs_ = pd.read_csv(obs_path, index_col=0)
        # 索引对齐检查
        if len(obs_) == len(adata_):
            adata_.obs = obs_
        else:
            print(f"[WARN] Obs length ({len(obs_)}) mismatch with Adata ({len(adata_)}). Keeping original obs.")
    return adata_

# ================= 1. 加载数据 =================
CRISPR_adata = load_adata(CRISPR_Name)
Drug_adata = load_adata(drug_Name)

# 格式化 Target 列
CRISPR_adata.obs["target"] = CRISPR_adata.obs["gene"].astype(str)
Drug_adata.obs["target"] = Drug_adata.obs["target"].astype(str)

# ================= 2. 提取与清洗靶点 =================

# --- A. 清洗 CRISPR 靶点 ---
exclude_crispr = ["non-targeting", "nan", "None"]
crispr_genes = set(CRISPR_adata.obs["target"].unique())
crispr_genes = {g for g in crispr_genes if g not in exclude_crispr}
print(f"\n[CRISPR] Unique Targets found: {len(crispr_genes)}")

# --- B. 清洗 Drug 靶点 ---
exclude_drug = ["control", "nan", "None"]
drug_targets_raw = set(Drug_adata.obs["target"].unique())

drug_genes_expanded = set()
drug_target_map = {} 

for t_str in drug_targets_raw:
    if t_str in exclude_drug:
        continue
    
    split_targets = t_str.split('_')
    for gene in split_targets:
        gene = gene.strip()
        drug_genes_expanded.add(gene)
        
        if gene not in drug_target_map:
            drug_target_map[gene] = []
        drug_target_map[gene].append(t_str)

print(f"\n[Drug] Unique Gene Targets (expanded): {len(drug_genes_expanded)}")

# ================= 3. 计算交集 =================

common_genes = crispr_genes.intersection(drug_genes_expanded)
sorted_common = sorted(list(common_genes))

print("\n" + "="*40)
print(f"🎉 FOUND OVERLAP: {len(common_genes)} GENES")
print("="*40)

if len(common_genes) == 0:
    raise ValueError("No overlapping targets found! Check gene symbol naming conventions.")

# 保存交集列表
pd.DataFrame(sorted_common, columns=["Common_Target"]).to_csv(os.path.join(OUTPUT_DIR, "Common_Targets_K562.csv"), index=False)

# ================= 4. 基因特征对齐 (Feature Alignment) =================
print("\n>>> Aligning Features (Genes)...")
common_vars = np.intersect1d(CRISPR_adata.var_names, Drug_adata.var_names)
print(f"   CRISPR Vars: {len(CRISPR_adata.var_names)}")
print(f"   Drug Vars:   {len(Drug_adata.var_names)}")
print(f"   Intersected: {len(common_vars)}")

# 【关键检查】如果交集过小，说明可能一个是Symbol一个是ID
if len(common_vars) < 1000:
    raise ValueError("Feature intersection is too small! One dataset might be using Ensembl IDs while the other uses Symbols.")

# 保存基因列表
pd.Series(common_vars).to_csv(os.path.join(OUTPUT_DIR, "shared_var_names.csv"), header=["gene_symbol"], index=False)

# ================= 5. 构建 CRISPR 子集 =================
print("\n[Processing CRISPR Data]...")

META_COLS_CRISPR = ['target', 'perturbation_type', 'perturbation_raw']

# 筛选
target_mask_crispr = CRISPR_adata.obs["target"].isin(set(common_genes) | {"non-targeting"})
subset_crispr = CRISPR_adata[target_mask_crispr, common_vars].copy()

# Metadata
meta_crispr = subset_crispr.obs[META_COLS_CRISPR].copy()
meta_crispr["benchmark_group"] = meta_crispr["target"].apply(lambda x: "Control" if x == "non-targeting" else "Test")
meta_crispr.to_csv(os.path.join(OUTPUT_DIR, "CRISPR_meta.csv"))

# Tensor (Save dense to save loading time later)
X_crispr = subset_crispr.X
if hasattr(X_crispr, "toarray"):
    X_crispr = X_crispr.toarray()
    
# 检查是否含 NaN
if np.isnan(X_crispr).any():
    print("[WARN] CRISPR data contains NaNs! Zero-filling...")
    X_crispr = np.nan_to_num(X_crispr)

torch.save(torch.tensor(X_crispr, dtype=torch.float32), os.path.join(OUTPUT_DIR, "CRISPR_counts.pt"))
print(f"  - Saved CRISPR subset: {X_crispr.shape}")

# 释放内存
del subset_crispr, X_crispr, CRISPR_adata
gc.collect()

# ================= 6. 构建 Drug 子集 =================
print("\n[Processing Drug Data]...")

META_COLS_DRUG = ['target', 'perturbation_type', 'perturbation_raw', 'time', 'dose_value']

# 辅助函数
def classify_specificity(raw_label):
    if raw_label == "control": return "Control"
    targets = raw_label.split('_')
    n = len(targets)
    if n == 1: return "The Cleanest Hits"
    if n > 4: return "The Promiscuous Hits"
    prefixes = [t[:3] for t in targets]
    if len(set(prefixes)) == 1: return "The Family Hits"
    return "The Promiscuous Hits"

def map_back_to_clean(dirty_label):
    if dirty_label == "control": return "control"
    parts = set(dirty_label.split('_'))
    hits = [g for g in sorted_common if g in parts]
    return ";".join(hits) if hits else "Other"

# 筛选
valid_drug_labels = set()
for gene in common_genes:
    valid_drug_labels.update(drug_target_map.get(gene, []))
valid_drug_labels.add("control")

target_mask_drug = Drug_adata.obs["target"].isin(valid_drug_labels)
subset_drug = Drug_adata[target_mask_drug, common_vars].copy()

# Metadata
meta_drug = subset_drug.obs[META_COLS_DRUG].copy()
meta_drug["benchmark_group"] = meta_drug["target"].apply(lambda x: "Control" if x == "control" else "Test")
meta_drug["clean_target_mapped"] = meta_drug["target"].apply(map_back_to_clean)
meta_drug["specificity_tier"] = meta_drug["target"].apply(classify_specificity)

print("\n[Specificity Tier Stats]:")
print(meta_drug[meta_drug["benchmark_group"]=="Test"]["specificity_tier"].value_counts())

meta_drug.to_csv(os.path.join(OUTPUT_DIR, "Drug_meta.csv"))

# Tensor
X_drug = subset_drug.X
if hasattr(X_drug, "toarray"):
    X_drug = X_drug.toarray()

if np.isnan(X_drug).any():
    print("[WARN] Drug data contains NaNs! Zero-filling...")
    X_drug = np.nan_to_num(X_drug)

torch.save(torch.tensor(X_drug, dtype=torch.float32), os.path.join(OUTPUT_DIR, "Drug_counts.pt"))
print(f"  - Saved Drug subset: {X_drug.shape}")

print(f"\n>>> ✅ All Benchmark Data Successfully Saved to: {OUTPUT_DIR}")