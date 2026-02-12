#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import gc
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

import scanpy as sc  # 用于读取 reference h5ad

# =========================
# User Config
# =========================
EVAL_DATA_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562"
OUTPUT_DIR    = "/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# scBERT local repo (你需要改成你本机 clone 的路径)
SCBERT_REPO_DIR = "/home/lenislin/Experiment/projects/Chem2Gen/benchmark/scBERT"

# scBERT checkpoint + gene2vec + demo reference
MODEL_CKPT  = "/home/lenislin/Experiment/projects/Chem2Gen/benchmark/scBERT/panglao_pretrain.pth"
GENE2VEC_NPY = "/home/lenislin/Experiment/projects/Chem2Gen/benchmark/scBERT/data/gene2vec_16906.npy"
REF_H5AD     = "/home/lenislin/Experiment/projects/Chem2Gen/benchmark/scBERT/data/Zheng68K.h5ad"

# token/bin config
BIN_NUM = 5          # 和官方 pretrain.py/predict.py 默认一致（bin_num=5）:contentReference[oaicite:11]{index=11}
GENE_NUM = 16906     # 你 demo 验证过 Zheng68K var_names length=16906
SEQ_LEN  = GENE_NUM + 1
CLASS    = BIN_NUM + 2

BATCH_SIZE = 24       # 16907 序列较长，建议从 1-4 试，显存够再调大
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# =========================
# Import scBERT
# =========================
sys.path.insert(0, SCBERT_REPO_DIR)
os.chdir(SCBERT_REPO_DIR)

# 尽量保证 gene2vec 文件在 repo 根目录可见（很多实现会用相对路径加载）
gene2vec_link = os.path.join(SCBERT_REPO_DIR, "gene2vec_16906.npy")
if not os.path.exists(gene2vec_link):
    try:
        os.symlink(GENE2VEC_NPY, gene2vec_link)
        print(f"[Info] Symlinked gene2vec -> {gene2vec_link}")
    except Exception as e:
        print(f"[Warn] Failed to symlink gene2vec: {e}")
        print("       If scBERT fails to load gene2vec positional embedding, "
              "please copy gene2vec_16906.npy into the scBERT repo root manually.")

from performer_pytorch import PerformerLM  # noqa: E402


# =========================
# Load Chem2Gen Eval Data
# =========================
def load_eval_data(eval_dir: str):
    var_names = pd.read_csv(os.path.join(eval_dir, "shared_var_names.csv"))["gene_symbol"].tolist()

    g_counts = torch.load(os.path.join(eval_dir, "CRISPR_counts.pt"))
    g_meta   = pd.read_csv(os.path.join(eval_dir, "CRISPR_meta.csv"), index_col=0)

    c_counts = torch.load(os.path.join(eval_dir, "Drug_counts.pt"))
    c_meta   = pd.read_csv(os.path.join(eval_dir, "Drug_meta.csv"), index_col=0)

    return var_names, g_counts, g_meta, c_counts, c_meta


# =========================
# Preprocess (normalize + log1p base2) then discretize
# =========================
def normalize_log1p_base2_dense(X: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """
    X: (B, G) raw counts, dense float
    returns: (B, G) float after normalize_total + log1p(base=2)
    """
    totals = X.sum(axis=1, keepdims=True)
    scale = np.ones_like(totals, dtype=np.float32)
    nonzero = totals.squeeze() > 0
    scale[nonzero] = (target_sum / totals[nonzero]).astype(np.float32)

    Xn = X * scale
    # log1p base 2: log2(1 + x)
    return np.log2(1.0 + Xn, dtype=np.float32)


def discretize_to_tokens(X_log: np.ndarray, bin_num: int) -> np.ndarray:
    """
    将 float expression 映射到 int token:
      token = floor(X_log), clip to [0..bin_num]
    这与官方 forward 前会再次 clip > bin_num 的行为保持一致。:contentReference[oaicite:12]{index=12}
    """
    T = np.floor(X_log).astype(np.int64)
    T[T < 0] = 0
    T[T > bin_num] = bin_num
    return T


# =========================
# Map genes to reference order (16906)
# =========================
def build_gene_mapper(shared_genes, ref_genes):
    ref_index = {g: i for i, g in enumerate(ref_genes)}
    src_cols = []
    ref_cols = []
    for j, g in enumerate(shared_genes):
        i = ref_index.get(g, None)
        if i is not None:
            src_cols.append(j)
            ref_cols.append(i)
    return np.array(src_cols, dtype=np.int64), np.array(ref_cols, dtype=np.int64)


def map_counts_batch_to_ref(counts_batch: torch.Tensor,
                            src_cols: np.ndarray,
                            ref_cols: np.ndarray,
                            ref_gene_num: int) -> np.ndarray:
    """
    counts_batch: (B, G_shared) torch tensor
    returns: (B, ref_gene_num) dense float32
    """
    X = counts_batch.cpu().numpy().astype(np.float32)
    out = np.zeros((X.shape[0], ref_gene_num), dtype=np.float32)
    out[:, ref_cols] = X[:, src_cols]
    return out


# =========================
# Load scBERT model
# =========================
def load_scbert_pretrained(ckpt_path: str, pos_embed_g2v: bool = True) -> PerformerLM:
    model = PerformerLM(
        num_tokens=CLASS,
        dim=200,
        depth=6,
        max_seq_len=SEQ_LEN,
        heads=10,
        local_attn_heads=0,
        g2v_position_emb=pos_embed_g2v,
    )

    ckpt = torch.load(ckpt_path, map_location="cpu")
    state = ckpt["model_state_dict"] if isinstance(ckpt, dict) and "model_state_dict" in ckpt else ckpt

    # 去掉 DDP 前缀
    new_state = {}
    for k, v in state.items():
        nk = k.replace("module.", "")
        new_state[nk] = v

    missing, unexpected = model.load_state_dict(new_state, strict=False)
    print(f"[scBERT] Loaded ckpt. missing={len(missing)}, unexpected={len(unexpected)}")

    model = model.to(DEVICE)
    model.eval()
    return model


# =========================
# Extract embeddings
# =========================
@torch.no_grad()
def extract_embeddings_for_counts(model: PerformerLM,
                                  counts: torch.Tensor,
                                  src_cols: np.ndarray,
                                  ref_cols: np.ndarray,
                                  out_path: str):
    """
    输出: (N_cells, 200) 默认取最后 appended token 的 hidden 作为 cell embedding
    """
    # 临时把 to_out 替换为 Identity，让 forward 返回 (B, L, 200)
    # predict.py 的分类 head 正是接收 (B, L, 200) 作为输入:contentReference[oaicite:13]{index=13}
    orig_to_out = model.to_out
    model.to_out = torch.nn.Identity()

    N = counts.shape[0]
    embs = []

    for start in tqdm(range(0, N, BATCH_SIZE), desc=f"Embedding -> {os.path.basename(out_path)}"):
        end = min(start + BATCH_SIZE, N)
        batch_counts = counts[start:end]

        # 1) gene order mapping -> (B, 16906)
        X_ref = map_counts_batch_to_ref(batch_counts, src_cols, ref_cols, GENE_NUM)

        # 2) normalize + log1p(base2)
        X_log = normalize_log1p_base2_dense(X_ref)

        # 3) discretize to tokens
        X_tok = discretize_to_tokens(X_log, BIN_NUM)  # (B, 16906)

        # 4) append 0 at end (len=16907)，与官方 pretrain/predict 一致:contentReference[oaicite:14]{index=14}
        cls_col = np.zeros((X_tok.shape[0], 1), dtype=np.int64)
        seq = np.concatenate([X_tok, cls_col], axis=1)

        seq_t = torch.from_numpy(seq).to(DEVICE, dtype=torch.long)

        # 5) forward -> (B, L, 200)
        hidden = model(seq_t)

        # 6) take last token embedding as cell embedding
        cell_emb = hidden[:, -1, :]  # (B, 200)
        embs.append(cell_emb.detach().cpu().numpy())

        # 控内存
        del X_ref, X_log, X_tok, seq, seq_t, hidden, cell_emb
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        gc.collect()

    model.to_out = orig_to_out

    emb_all = np.vstack(embs)
    np.save(out_path, emb_all)
    print(f"[Done] Saved embeddings: {out_path} | shape={emb_all.shape}")


def main():
    print(">>> Loading Chem2Gen Eval data ...")
    shared_genes, g_counts, g_meta, c_counts, c_meta = load_eval_data(EVAL_DATA_DIR)
    print(f"    shared genes: {len(shared_genes)}")
    print(f"    CRISPR cells: {g_counts.shape[0]} | genes={g_counts.shape[1]}")
    print(f"    Drug   cells: {c_counts.shape[0]} | genes={c_counts.shape[1]}")

    print(">>> Loading reference gene order from Zheng68K.h5ad ...")
    ref = sc.read_h5ad(REF_H5AD)
    ref_genes = ref.var_names.tolist()
    assert len(ref_genes) == GENE_NUM, f"REF_H5AD gene_num mismatch: {len(ref_genes)} != {GENE_NUM}"
    print(f"    ref genes: {len(ref_genes)}")

    src_cols, ref_cols = build_gene_mapper(shared_genes, ref_genes)
    print(f">>> Gene intersection (shared_var_names ∩ ref_genes): {len(src_cols)}")

    print(">>> Loading scBERT pretrained model ...")
    os.chdir("/home/lenislin/Experiment/projects/Chem2Gen/benchmark/scBERT/performer_pytorch")
    model = load_scbert_pretrained(MODEL_CKPT, pos_embed_g2v=True)

    # Extract CRISPR
    out_g = os.path.join(OUTPUT_DIR, "scBERT_CRISPR_emb.npy")
    extract_embeddings_for_counts(model, g_counts, src_cols, ref_cols, out_g)

    # Extract Drug
    out_c = os.path.join(OUTPUT_DIR, "scBERT_Drug_emb.npy")
    extract_embeddings_for_counts(model, c_counts, src_cols, ref_cols, out_c)


if __name__ == "__main__":
    main()
