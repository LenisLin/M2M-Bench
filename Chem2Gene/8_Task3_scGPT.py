import os
import json
import inspect
import numpy as np
import pandas as pd
import torch
import anndata as ad

import scgpt as scg
from scgpt.utils import set_seed  # 若没有可删掉

from Task3_Functions import OUTPUT_DIR, EVAL_DATA_DIR, Chem2GenEvalDataLoader

# ================= 配置 =================
MODEL_DIR = "/home/lenislin/Experiment/projects/Chem2Gen/benchmark/scGPT/checkpoints"
ARGS_FILE = os.path.join(MODEL_DIR, "args.json")

BATCH_SIZE = 96
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# 为了可复现：注意 scGPT 的 embedding pipeline 在长度超过 max_length 时会 sampling :contentReference[oaicite:3]{index=3}
set_seed(0)
torch.manual_seed(0)
np.random.seed(0)

def _build_adata(counts_t: torch.Tensor, meta_df: pd.DataFrame, var_names: list[str]) -> ad.AnnData:
    """
    把你的 (N_cells, N_genes) counts + meta + shared_var_names 转成 scGPT 的 AnnData 输入。
    embed_data 会按 gene_col 去做 vocab 映射与过滤。:contentReference[oaicite:4]{index=4}
    """
    X = counts_t.detach().cpu().numpy().astype(np.float32)
    obs = meta_df.copy()

    # embed_data 要求 cell_type_key 存在（默认 cell_type）:contentReference[oaicite:5]{index=5}
    if "cell_type" not in obs.columns:
        obs["cell_type"] = "NA"

    adata = ad.AnnData(X=X, obs=obs)
    adata.var["gene_symbol"] = var_names  # 你的 shared_var_names.csv 是 gene_symbol
    return adata

def _call_embed_data(adata: ad.AnnData, max_length: int, batch_size: int):
    """
    兼容不同 scGPT 版本：
    - 新版 embed_data 有 use_fast_transformer 参数（文档里有）:contentReference[oaicite:6]{index=6}
    - 老版可能没有该参数（但核心流程一致）:contentReference[oaicite:7]{index=7}
    """
    # 不同版本 embed_data 入口位置可能不同：scg.tasks.embed_data vs scgpt.tasks.cell_emb.embed_data
    if hasattr(scg, "tasks") and hasattr(scg.tasks, "embed_data"):
        embed_fn = scg.tasks.embed_data
    else:
        from scgpt.tasks.cell_emb import embed_data as embed_fn  # v0.1.x 文档路径 :contentReference[oaicite:8]{index=8}

    sig = inspect.signature(embed_fn)
    kwargs = dict(
        adata_or_file=adata,
        model_dir=MODEL_DIR,
        gene_col="gene_symbol",
        max_length=max_length,
        batch_size=batch_size,
        device=DEVICE,
        return_new_adata=True,
    )
    if "use_fast_transformer" in sig.parameters:
        # 没装 flash-attn / 不想用就 False
        kwargs["use_fast_transformer"] = False  # 或 True（如果你的环境支持）:contentReference[oaicite:9]{index=9}

    return embed_fn(**kwargs)

def run_inference():
    base_loader = Chem2GenEvalDataLoader(EVAL_DATA_DIR)

    # 读取 max_length：优先用 args.json 里的 max_seq_len（若不存在则用 1200）
    with open(ARGS_FILE, "r") as f:
        cfg = json.load(f)
    max_length = int(cfg.get("max_seq_len", cfg.get("max_length", 1200)))

    print(f">>> scGPT embed config: max_length={max_length}, batch_size={BATCH_SIZE}, device={DEVICE}")

    # ========== CRISPR ==========
    adata_g = _build_adata(base_loader.g_counts, base_loader.g_meta, base_loader.var_names)
    g_emb_adata = _call_embed_data(adata=adata_g, max_length=max_length, batch_size=BATCH_SIZE)
    g_emb = np.asarray(g_emb_adata.X, dtype=np.float32)

    save_g = os.path.join(OUTPUT_DIR, "scGPT_CRISPR_emb.npy")
    np.save(save_g, g_emb)
    print(f"✅ Saved CRISPR embeddings: {save_g} | shape={g_emb.shape}")

    # ========== Drug ==========
    adata_c = _build_adata(base_loader.c_counts, base_loader.c_meta, base_loader.var_names)
    c_emb_adata = _call_embed_data(adata_c, max_length=max_length, batch_size=BATCH_SIZE)
    c_emb = np.asarray(c_emb_adata.X, dtype=np.float32)

    save_c = os.path.join(OUTPUT_DIR, "scGPT_Drug_emb.npy")
    np.save(save_c, c_emb)
    print(f"✅ Saved Drug embeddings: {save_c} | shape={c_emb.shape}")

    # Sanity checks (与你的 evaluator 对齐要求)
    assert g_emb.shape[0] == len(base_loader.g_meta)
    assert c_emb.shape[0] == len(base_loader.c_meta)

if __name__ == "__main__":
    run_inference()
