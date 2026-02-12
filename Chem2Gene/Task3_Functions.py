import os
import torch
import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity
import gc
import sys

# ================= 配置 =================
# 数据集目录 (存放 meta.csv)
EVAL_DATA_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562"

# 结果输出目录
OUTPUT_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

class Chem2GenEvalDataLoader:
    """
    负责加载我们定义好的 Chem2Gen-FM-Eval 数据集
    """
    def __init__(self, data_dir):
        print(f">>> Loading Benchmark Data from {data_dir}...")
        
        # 1. Load Gene Features
        self.var_names = pd.read_csv(os.path.join(data_dir, "shared_var_names.csv"))["gene_symbol"].tolist()
        
        # 2. Load CRISPR Data
        self.g_counts = torch.load(os.path.join(data_dir, "CRISPR_counts.pt"))
        self.g_meta = pd.read_csv(os.path.join(data_dir, "CRISPR_meta.csv"), index_col=0)
        
        # 3. Load Drug Data
        self.c_counts = torch.load(os.path.join(data_dir, "Drug_counts.pt"))
        self.c_meta = pd.read_csv(os.path.join(data_dir, "Drug_meta.csv"), index_col=0)
        
        # 4. Load Common Targets List (The Evaluation Task List)
        self.tasks = pd.read_csv(os.path.join(data_dir, "Common_Targets_K562.csv"))["Common_Target"].tolist()
        
        print(f"    Loaded {len(self.var_names)} genes.")
        print(f"    CRISPR: {self.g_counts.shape[0]} cells (Control + Test)")
        print(f"    Drug:   {self.c_counts.shape[0]} cells (Control + Test)")
        print(f"    Tasks:  {len(self.tasks)} common targets to evaluate.")

    def get_control_indices(self):
        """返回所有 Control 组的索引"""
        g_ctrl_idx = np.where(self.g_meta["benchmark_group"] == "Control")[0]
        c_ctrl_idx = np.where(self.c_meta["benchmark_group"] == "Control")[0]
        return g_ctrl_idx, c_ctrl_idx

    def get_target_indices(self, target_name):
        """返回指定 Target 的索引 (Genetic & Chemical)"""
        # Genetic: 直接匹配 target 列
        g_idx = np.where(self.g_meta["target"] == target_name)[0]
        
        # Chemical: 需要匹配 clean_target_mapped 包含该 target 的行
        # (处理多靶点映射: "AURKA;AURKB" 应该被 AURKA 选中)
        # 这里的逻辑是：只要这个药物被映射的目标里包含当前 target，就算。
        c_mask = self.c_meta["clean_target_mapped"].apply(lambda x: target_name in x.split(";"))
        c_idx = np.where(c_mask)[0]
        
        # 同时提取该 Drug 的 metadata (用于后续分析 Tier)
        # 如果有多个不同的 drug label 命中 (比如 dirty drug)，我们取出现次数最多的 tier 作为代表，或者保留 tier 信息
        tier = "Unknown"
        if len(c_idx) > 0:
            # 简单策略：取这些样本中大多数的 tier
            tier = self.c_meta.iloc[c_idx]["specificity_tier"].mode()[0]
            
        return g_idx, c_idx, tier

class LatentSpaceEvaluator:
    """
    轻量级评测器：只负责加载 Embedding 和 Metadata，执行潜空间算术。
    """
    def __init__(self, dataset_dir):
        print(f">>> Initializing Evaluator based on metadata from: {dataset_dir}")
        
        # 1. 加载 Metadata (这是地图，告诉我们要取哪几行)
        self.g_meta = pd.read_csv(os.path.join(dataset_dir, "CRISPR_meta.csv"), index_col=0)
        self.c_meta = pd.read_csv(os.path.join(dataset_dir, "Drug_meta.csv"), index_col=0)
        
        # 2. 加载任务列表 (共有的靶点)
        self.tasks = pd.read_csv(os.path.join(dataset_dir, "Common_Targets_K562.csv"))["Common_Target"].tolist()
        
        # 3. 预计算 Control 组索引 (加速后续计算)
        self.g_ctrl_idx = np.where(self.g_meta["benchmark_group"] == "Control")[0]
        self.c_ctrl_idx = np.where(self.c_meta["benchmark_group"] == "Control")[0]
        
        print(f"    Loaded Metadata: {len(self.g_meta)} Genetic cells, {len(self.c_meta)} Chemical cells.")
        print(f"    Evaluation Tasks: {len(self.tasks)} targets.")

    def load_embeddings(self, model_name, g_emb_path, c_emb_path):
        """
        加载特定模型生成的 Embedding 文件
        """
        print(f"\n>>> [{model_name}] Loading Embeddings...")
        
        # 支持 .npy 或 .pt
        if g_emb_path.endswith('.npy'):
            self.g_emb = np.load(g_emb_path)
            self.c_emb = np.load(c_emb_path)
        else:
            self.g_emb = torch.load(g_emb_path).cpu().numpy()
            self.c_emb = torch.load(c_emb_path).cpu().numpy()
            
        self.current_model_name = model_name
        
        # Sanity Check: 行数必须对齐
        assert self.g_emb.shape[0] == len(self.g_meta), \
            f"Shape Mismatch! Genetic Meta has {len(self.g_meta)} rows, but Embedding has {self.g_emb.shape[0]}"
        assert self.c_emb.shape[0] == len(self.c_meta), \
            f"Shape Mismatch! Drug Meta has {len(self.c_meta)} rows, but Embedding has {self.c_emb.shape[0]}"
            
        print(f"    ✅ Shape check passed. Latent Dim: {self.g_emb.shape[1]}")

    def compute_latent_scores(self):
        """
        核心逻辑：计算 Latent Delta Cosine
        """
        print(f">>> [{self.current_model_name}] Computing Latent Fidelity Scores...")
        
        # 1. 计算全局基线质心 (Global Baseline Centroids)
        # 这一步消除了模型本身的 bias (比如模型倾向于把所有 embedding 放在某个象限)
        Z_g_ctrl = np.mean(self.g_emb[self.g_ctrl_idx], axis=0)
        Z_c_ctrl = np.mean(self.c_emb[self.c_ctrl_idx], axis=0)
        
        results = []
        
        for target in tqdm(self.tasks, desc="Scoring"):
            # --- A. 获取 Genetic 组索引 ---
            g_idx = np.where(self.g_meta["target"] == target)[0]
            
            # --- B. 获取 Chemical 组索引 (处理多靶点映射) ---
            # 逻辑：clean_target_mapped 列包含目标基因 (例如 "AURKA;AURKB" 包含 "AURKA")
            c_mask = self.c_meta["clean_target_mapped"].apply(lambda x: target in x.split(";"))
            c_idx = np.where(c_mask)[0]
            
            if len(g_idx) == 0 or len(c_idx) == 0:
                continue
            
            # --- C. 获取该药物的特异性分级 (用于后续分析) ---
            # 取众数 (Mode)
            tier = self.c_meta.iloc[c_idx]["specificity_tier"].mode()[0]
            
            # --- D. 潜空间算术 (Latent Arithmetic) ---
            # Mean Pooling
            Z_g_trt = np.mean(self.g_emb[g_idx], axis=0)
            Z_c_trt = np.mean(self.c_emb[c_idx], axis=0)
            
            # Subtract Baseline (Delta)
            Delta_g = Z_g_trt - Z_g_ctrl
            Delta_c = Z_c_trt - Z_c_ctrl
            
            # Cosine Similarity
            sim = cosine_similarity(Delta_g.reshape(1, -1), Delta_c.reshape(1, -1))[0][0]
            
            results.append({
                "Model": self.current_model_name,
                "Target": target,
                "Specificity_Tier": tier,
                "Latent_Cosine": sim,
                "N_Cells_G": len(g_idx),
                "N_Cells_C": len(c_idx)
            })
            
        return pd.DataFrame(results)

# ================= 使用示例 =================
# if __name__ == "__main__":
#     # 1. 初始化评测器
#     evaluator = LatentSpaceEvaluator(EVAL_DATA_DIR)
    
#     all_scores = []
    
#     # 假设你已经跑完了 scGPT，文件存放在 results/scgpt/ 下
#     # MODEL 1: scGPT
#     evaluator.load_embeddings(
#         "scGPT", 
#         "/path/to/results/scgpt_CRISPR_emb.npy", 
#         "/path/to/results/scgpt_Drug_emb.npy"
#     )
#     scores_scgpt = evaluator.compute_latent_scores()
#     all_scores.append(scores_scgpt)
    
#     # MODEL 2: Geneformer
#     evaluator.load_embeddings(
#         "Geneformer", 
#         "/path/to/results/gf_CRISPR_emb.npy", 
#         "/path/to/results/gf_Drug_emb.npy"
#     )
#     scores_gf = evaluator.compute_latent_scores()
#     all_scores.append(scores_gf)
    
#     # 保存最终大表
#     final_df = pd.concat(all_scores, ignore_index=True)
#     final_df.to_csv(os.path.join(OUTPUT_DIR, "Task3_AllModels_Benchmark.csv"), index=False)