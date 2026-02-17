"""
Fixed Geneformer GRN Validation Script
=======================================
Bug fix: The original geneformer_grn_validation.py used gc104M dictionaries and
treated adata.var_names (Ensembl IDs) as gene symbols, causing complete namespace
mismatch → AUROC=0.50, P/R/F1=0.00.

This script uses the SAME dictionary set (gc30M) and mapping logic as the working
geneformer_grn_pipeline.py, adapted into the validation framework structure.
"""

import os, sys, pickle, json, time, warnings
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from collections import defaultdict
import scipy.sparse as sp
import torch
from transformers import BertForMaskedLM
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score, recall_score, f1_score

warnings.filterwarnings("ignore")
np.random.seed(42)
torch.manual_seed(42)

# ── Paths (same as working pipeline) ──
MODEL_DIR = Path(r"D:\openclaw\intelligence-augmentation\models\Geneformer\Geneformer-V1-10M")
DICT_DIR = Path(r"D:\openclaw\intelligence-augmentation\models\Geneformer\geneformer\gene_dictionaries_30m")
DATA_PATH = Path(r"D:\openclaw\intelligence-augmentation\data\brain_scrna\DLPFC_11k.h5ad")
TRRUST_PATH = Path(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\external\networks\trrust_human.tsv")
DOROTHEA_PATH = Path(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\external\networks\dorothea_human.tsv")
OUT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\multi_model_proper")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Device: {DEVICE}", flush=True)

# ── Step 1: Load dictionaries (gc30M — same as working pipeline) ──
print("\n=== Step 1: Loading Geneformer gc30M dictionaries ===", flush=True)

with open(DICT_DIR / "token_dictionary_gc30M.pkl", "rb") as f:
    token_dict = pickle.load(f)  # ensembl_id -> token_id
with open(DICT_DIR / "gene_name_id_dict_gc30M.pkl", "rb") as f:
    gene_name_id_dict = pickle.load(f)  # gene_symbol -> ensembl_id
with open(DICT_DIR / "gene_median_dictionary_gc30M.pkl", "rb") as f:
    gene_median_dict = pickle.load(f)  # ensembl_id -> median

print(f"  Tokens: {len(token_dict)}, Gene names: {len(gene_name_id_dict)}, Medians: {len(gene_median_dict)}", flush=True)

# Inverse maps
id_to_token = {v: k for k, v in token_dict.items()}
id_to_gene = {v: k for k, v in gene_name_id_dict.items()}

def token_to_gene(tid):
    ens = id_to_token.get(tid, None)
    if ens is None: return None
    return id_to_gene.get(ens, ens)

# ── Step 2: Load data ──
print("\n=== Step 2: Loading DLPFC data ===", flush=True)
adata = sc.read_h5ad(DATA_PATH)
print(f"  Shape: {adata.shape}", flush=True)

# ── Step 3: Tokenize (same logic as working pipeline) ──
print("\n=== Step 3: Tokenizing cells ===", flush=True)

gene_to_token = {}
gene_to_median = {}
for i, g in enumerate(adata.var_names):
    ens = g.split('.')[0]
    if ens in token_dict and ens in gene_median_dict:
        med = gene_median_dict[ens]
        if med > 0:
            gene_to_token[i] = token_dict[ens]
            gene_to_median[i] = med

valid_gene_idx = sorted(gene_to_token.keys())
print(f"  {len(valid_gene_idx)}/{adata.n_vars} genes have token mapping", flush=True)

X = adata.X
if not sp.issparse(X):
    X = sp.csr_matrix(X)
X = X.tocsc()

all_tokens = []
MAX_LEN = 2048
MIN_GENES = 50

valid_idx_arr = np.array(valid_gene_idx)
token_arr = np.array([gene_to_token[i] for i in valid_gene_idx])
median_arr = np.array([gene_to_median[i] for i in valid_gene_idx])

X_valid = X[:, valid_idx_arr].tocsr()
print(f"  Subsetted matrix: {X_valid.shape}", flush=True)

for ci in range(adata.n_obs):
    row = X_valid[ci]
    if sp.issparse(row):
        row = row.toarray().flatten()
    else:
        row = np.asarray(row).flatten()
    nz = np.nonzero(row)[0]
    nz_vals = row[nz]
    
    if len(nz) < MIN_GENES:
        continue
    
    norm_vals = nz_vals / median_arr[nz]
    order = np.argsort(-norm_vals)[:MAX_LEN]
    toks = token_arr[nz[order]].tolist()
    all_tokens.append(toks)

print(f"  Total valid cells: {len(all_tokens)}", flush=True)

# ── Step 4: Load model ──
print("\n=== Step 4: Loading Geneformer V1-10M ===", flush=True)
model = BertForMaskedLM.from_pretrained(str(MODEL_DIR), output_attentions=True)
model = model.to(DEVICE)
model.eval()
print(f"  Model loaded: {model.config.num_hidden_layers}L x {model.config.num_attention_heads}H", flush=True)

# ── Step 5: Load reference GRNs ──
print("\n=== Step 5: Loading reference GRNs ===", flush=True)
trrust = pd.read_csv(TRRUST_PATH, sep='\t', header=None, names=['source','target','type','pmid'])
trrust_edges = set(zip(trrust['source'].str.upper(), trrust['target'].str.upper()))
print(f"  TRRUST: {len(trrust_edges)} unique edges", flush=True)

dorothea = pd.read_csv(DOROTHEA_PATH, sep='\t')
dorothea_edges = set(zip(dorothea.iloc[:,0].str.upper(), dorothea.iloc[:,1].str.upper()))
print(f"  DoRothEA: {len(dorothea_edges)} unique edges", flush=True)

# ── Step 6: Attention extraction (same as working pipeline) ──
print("\n=== Step 6: Attention-based GRN inference ===", flush=True)

def extract_attention_grn(cell_tokens_list, model, device, batch_size=4, max_seq=512):
    edge_scores = defaultdict(float)
    edge_counts = defaultdict(int)
    n = len(cell_tokens_list)
    
    for b_start in range(0, n, batch_size):
        b_end = min(b_start + batch_size, n)
        batch = [t[:max_seq] for t in cell_tokens_list[b_start:b_end]]
        
        max_len = max(len(t) for t in batch)
        input_ids = torch.zeros(len(batch), max_len, dtype=torch.long)
        attn_mask = torch.zeros(len(batch), max_len, dtype=torch.long)
        
        for i, toks in enumerate(batch):
            input_ids[i, :len(toks)] = torch.tensor(toks)
            attn_mask[i, :len(toks)] = 1
        
        input_ids = input_ids.to(device)
        attn_mask = attn_mask.to(device)
        
        with torch.no_grad():
            out = model(input_ids=input_ids, attention_mask=attn_mask, output_attentions=True)
        
        attn = torch.stack(out.attentions).mean(dim=(0, 2)).cpu().numpy()
        
        for i, toks in enumerate(batch):
            seq_len = len(toks)
            a = attn[i, :seq_len, :seq_len]
            
            genes = []
            for tid in toks:
                g = token_to_gene(tid)
                genes.append(g.upper() if g else None)
            
            TOP_K_PER_ROW = 20
            for r in range(seq_len):
                if genes[r] is None:
                    continue
                row = a[r].copy()
                row[r] = -1
                top_idx = np.argpartition(row, -TOP_K_PER_ROW)[-TOP_K_PER_ROW:]
                for c in top_idx:
                    if genes[c] is None:
                        continue
                    edge = (genes[r], genes[c])
                    edge_scores[edge] += row[c]
                    edge_counts[edge] += 1
        
        del out, attn, input_ids, attn_mask
        if device == "cuda":
            torch.cuda.empty_cache()
        
        if (b_end) % 50 == 0 or b_end == n:
            print(f"    {b_end}/{n} cells, {len(edge_scores)} edges", flush=True)
    
    return {e: edge_scores[e] / edge_counts[e] for e in edge_scores}

def evaluate_grn(edge_scores, ref_edges, name=""):
    ref_genes = set()
    for s, t in ref_edges:
        ref_genes.add(s); ref_genes.add(t)
    
    filtered = {e: s for e, s in edge_scores.items() if e[0] in ref_genes and e[1] in ref_genes}
    if len(filtered) == 0:
        print(f"    {name}: No overlap!", flush=True)
        return {}
    
    sorted_edges = sorted(filtered.items(), key=lambda x: x[1], reverse=True)
    labels = [1 if e in ref_edges else 0 for e, _ in sorted_edges]
    scores = [s for _, s in sorted_edges]
    
    n_pos = sum(labels)
    results = {'n_filtered': len(filtered), 'n_positive': n_pos}
    print(f"    {name}: {len(filtered)} edges, {n_pos} positives ({n_pos/len(filtered)*100:.2f}%)", flush=True)
    
    if 0 < n_pos < len(labels):
        results['auroc'] = roc_auc_score(labels, scores)
        results['auprc'] = average_precision_score(labels, scores)
        print(f"    AUROC={results['auroc']:.4f}, AUPRC={results['auprc']:.4f}", flush=True)
    
    for k in [500, 1000, 5000]:
        if k > len(sorted_edges): continue
        top_set = set(e for e, _ in sorted_edges[:k])
        tp = len(top_set & ref_edges)
        p = tp / k; r = tp / len(ref_edges) if ref_edges else 0
        f1 = 2*p*r/(p+r) if (p+r) > 0 else 0
        results[f'top{k}_p'] = p; results[f'top{k}_r'] = r; results[f'top{k}_f1'] = f1
        print(f"    Top-{k}: P={p:.4f} R={r:.4f} F1={f1:.4f}", flush=True)
    
    return results

# ── Step 7: Run validation at multiple cell counts ──
print("\n=== Step 7: Validation experiments ===", flush=True)

cell_counts = [200, 500, 1000]
all_results = {}

for nc in cell_counts:
    print(f"\n--- {nc} cells ---", flush=True)
    t0 = time.time()
    
    actual_n = min(nc, len(all_tokens))
    idx = np.random.choice(len(all_tokens), actual_n, replace=False)
    subset = [all_tokens[i] for i in idx]
    
    edge_scores = extract_attention_grn(subset, model, DEVICE, batch_size=4, max_seq=512)
    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s, Edges: {len(edge_scores)}", flush=True)
    
    r = {'n_cells': actual_n, 'n_edges': len(edge_scores), 'time_s': elapsed}
    for ref_name, ref_set in [('trrust', trrust_edges), ('dorothea', dorothea_edges)]:
        r[ref_name] = evaluate_grn(edge_scores, ref_set, ref_name)
    all_results[nc] = r

# ── Save results ──
out_path = OUT_DIR / "geneformer_grn_validation_fixed_results.json"
with open(out_path, "w") as f:
    json.dump(all_results, f, indent=2, default=str)

print("\n" + "="*80)
print("FIXED VALIDATION RESULTS")
print("="*80)
for nc in cell_counts:
    r = all_results[nc]
    ta = r.get('trrust', {}).get('auroc', 'N/A')
    da = r.get('dorothea', {}).get('auroc', 'N/A')
    print(f"  {nc} cells: TRRUST AUROC={ta}, DoRothEA AUROC={da}")
print(f"\nResults saved to: {out_path}")
print("Bug fix: Now uses gc30M dictionaries with correct Ensembl ID mapping")
