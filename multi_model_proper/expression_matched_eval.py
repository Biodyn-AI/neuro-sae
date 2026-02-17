"""
Expression-Matched Negative Evaluation for Attention-Based GRN Recovery.

For each TRRUST positive edge (TF→target), we sample matched negatives where:
- The negative TF has similar marginal attention (±20%) to the positive TF
- The negative target has similar marginal attention (±20%) to the positive target

We use attention marginals as the matching statistic because:
1. In scGPT, attention marginals correlate strongly with expression level
2. Matching on attention marginals is MORE conservative than matching on expression
   (it controls for exactly the confound we're worried about)

If attention AUROC > 0.50 on matched negatives: GENUINE pairwise signal
If attention AUROC ≈ 0.50: no pairwise structure beyond marginals
"""

import numpy as np
import json
from sklearn.metrics import roc_auc_score
from collections import defaultdict
import os

np.random.seed(42)

# --- Load data ---
base = r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\outputs\head_layer_checkpoint_replication"
genes = json.load(open(os.path.join(base, "gene_list_1200.json")))
gene2idx = {g: i for i, g in enumerate(genes)}
n_genes = len(genes)

# Attention scores: (12 layers, 8 heads, 1200, 1200)
scores = np.load(os.path.join(base, "brain", "attention_scores_head_layer.npy"))
print(f"Attention shape: {scores.shape}, dtype: {scores.dtype}")

# Aggregate attention: mean over all layers and heads
attn_mean = scores.mean(axis=(0, 1)).astype(np.float32)  # (1200, 1200)
print(f"Mean attention matrix: {attn_mean.shape}")

# Also try best head (layer 13 head 14 equivalent, or find best)
# For now use the mean across all heads as the primary score

# --- Load TRRUST ---
trrust_path = os.path.join(r"D:\openclaw\biodyn-nmi-paper\multi_model_proper", "trrust_rawdata.human.tsv")
gene_set = set(genes)
positive_edges = set()
with open(trrust_path) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            tf, target = parts[0], parts[1]
            if tf in gene_set and target in gene_set and tf != target:
                positive_edges.add((tf, target))

positive_edges = list(positive_edges)
print(f"TRRUST positive edges: {len(positive_edges)}")

# --- Compute marginal statistics ---
# Row mean = mean attention GIVEN by gene i (as source/TF)
row_mean = attn_mean.mean(axis=1)  # (1200,)
# Col mean = mean attention RECEIVED by gene j (as target)
col_mean = attn_mean.mean(axis=0)  # (1200,)

# Also compute variance of attention per gene (for variance baseline)
attn_var_row = attn_mean.var(axis=1)
attn_var_col = attn_mean.var(axis=0)

print(f"Row mean range: [{row_mean.min():.6f}, {row_mean.max():.6f}]")
print(f"Col mean range: [{col_mean.min():.6f}, {col_mean.max():.6f}]")

# --- Expression-matched sampling ---
def find_matched_genes(gene_idx, marginals, tolerance=0.20):
    """Find genes with similar marginal attention (±tolerance fraction)."""
    val = marginals[gene_idx]
    low = val * (1 - tolerance)
    high = val * (1 + tolerance)
    # Handle near-zero values
    if val < 1e-8:
        high = 1e-7
    matched = np.where((marginals >= low) & (marginals <= high))[0]
    return matched

N_MATCHED = 50  # negatives per positive
results_per_edge = []
all_labels = []
all_attn_scores = []
all_marginal_product = []
all_var_product = []

positive_set = set((gene2idx[tf], gene2idx[tgt]) for tf, tgt in positive_edges)

for tf_name, tgt_name in positive_edges:
    tf_idx = gene2idx[tf_name]
    tgt_idx = gene2idx[tgt_name]
    
    # Find matched TFs and targets
    matched_tfs = find_matched_genes(tf_idx, row_mean)
    matched_tgts = find_matched_genes(tgt_idx, col_mean)
    
    # Remove self and actual positives
    neg_candidates = []
    for mt in matched_tfs:
        for mg in matched_tgts:
            if mt != mg and (mt, mg) not in positive_set:
                neg_candidates.append((mt, mg))
    
    if len(neg_candidates) < 5:
        print(f"  WARNING: Only {len(neg_candidates)} matched negatives for {tf_name}->{tgt_name}, relaxing tolerance")
        matched_tfs = find_matched_genes(tf_idx, row_mean, tolerance=0.40)
        matched_tgts = find_matched_genes(tgt_idx, col_mean, tolerance=0.40)
        neg_candidates = [(mt, mg) for mt in matched_tfs for mg in matched_tgts 
                         if mt != mg and (mt, mg) not in positive_set]
    
    n_sample = min(N_MATCHED, len(neg_candidates))
    if n_sample == 0:
        print(f"  SKIP: No matched negatives for {tf_name}->{tgt_name}")
        continue
    
    chosen = [neg_candidates[i] for i in np.random.choice(len(neg_candidates), n_sample, replace=False)]
    
    # Positive
    all_labels.append(1)
    all_attn_scores.append(float(attn_mean[tf_idx, tgt_idx]))
    all_marginal_product.append(float(row_mean[tf_idx] * col_mean[tgt_idx]))
    all_var_product.append(float(attn_var_row[tf_idx] * attn_var_col[tgt_idx]))
    
    # Negatives
    for ni, nj in chosen:
        all_labels.append(0)
        all_attn_scores.append(float(attn_mean[ni, nj]))
        all_marginal_product.append(float(row_mean[ni] * col_mean[nj]))
        all_var_product.append(float(attn_var_row[ni] * attn_var_col[nj]))

all_labels = np.array(all_labels)
all_attn_scores = np.array(all_attn_scores)
all_marginal_product = np.array(all_marginal_product)
all_var_product = np.array(all_var_product)

print(f"\nTotal samples: {len(all_labels)} ({all_labels.sum()} pos, {(1-all_labels).sum()} neg)")

# --- Compute AUROCs ---
auroc_attention = roc_auc_score(all_labels, all_attn_scores)
auroc_marginal = roc_auc_score(all_labels, all_marginal_product)
auroc_variance = roc_auc_score(all_labels, all_var_product)

print(f"\n=== EXPRESSION-MATCHED AUROC ===")
print(f"Attention (mean all heads):   {auroc_attention:.4f}")
print(f"Marginal product baseline:    {auroc_marginal:.4f}")
print(f"Variance product baseline:    {auroc_variance:.4f}")

# --- Also try per-head best ---
best_head_auroc = 0
best_head = None
for layer in range(scores.shape[0]):
    for head in range(scores.shape[1]):
        head_attn = scores[layer, head].astype(np.float32)
        head_scores = []
        head_labels = []
        pos_idx = 0
        neg_idx = 0
        for i, lab in enumerate(all_labels):
            # Reconstruct which gene pairs we used
            pass  # Need to store indices
        
# Redo with stored indices for per-head analysis
all_pairs = []  # (tf_idx, tgt_idx)
all_labels2 = []
np.random.seed(42)  # Reset for reproducibility

for tf_name, tgt_name in positive_edges:
    tf_idx = gene2idx[tf_name]
    tgt_idx = gene2idx[tgt_name]
    
    matched_tfs = find_matched_genes(tf_idx, row_mean)
    matched_tgts = find_matched_genes(tgt_idx, col_mean)
    
    neg_candidates = [(mt, mg) for mt in matched_tfs for mg in matched_tgts 
                     if mt != mg and (mt, mg) not in positive_set]
    
    if len(neg_candidates) < 5:
        matched_tfs = find_matched_genes(tf_idx, row_mean, tolerance=0.40)
        matched_tgts = find_matched_genes(tgt_idx, col_mean, tolerance=0.40)
        neg_candidates = [(mt, mg) for mt in matched_tfs for mg in matched_tgts 
                         if mt != mg and (mt, mg) not in positive_set]
    
    n_sample = min(N_MATCHED, len(neg_candidates))
    if n_sample == 0:
        continue
    
    chosen = [neg_candidates[i] for i in np.random.choice(len(neg_candidates), n_sample, replace=False)]
    
    all_pairs.append((tf_idx, tgt_idx))
    all_labels2.append(1)
    for ni, nj in chosen:
        all_pairs.append((ni, nj))
        all_labels2.append(0)

all_labels2 = np.array(all_labels2)
all_pairs = np.array(all_pairs)

# Per-head AUROC
print(f"\n=== PER-HEAD AUROC (expression-matched) ===")
head_aurocs = {}
for layer in range(scores.shape[0]):
    for head in range(scores.shape[1]):
        head_attn = scores[layer, head].astype(np.float32)
        head_vals = np.array([head_attn[i, j] for i, j in all_pairs])
        try:
            auc = roc_auc_score(all_labels2, head_vals)
        except:
            auc = 0.5
        head_aurocs[f"L{layer}_H{head}"] = auc

# Sort and show top 10
sorted_heads = sorted(head_aurocs.items(), key=lambda x: x[1], reverse=True)
for name, auc in sorted_heads[:10]:
    print(f"  {name}: {auc:.4f}")

best_head_name, best_head_auroc = sorted_heads[0]
print(f"\nBest head: {best_head_name} = {best_head_auroc:.4f}")
print(f"Median head AUROC: {np.median(list(head_aurocs.values())):.4f}")

# --- Also try with actual expression data from the immune h5ad ---
# Use row/col sums of attention as expression proxy (already done via marginals)

# --- Save results ---
output = {
    "experiment": "expression_matched_negative_evaluation",
    "description": "AUROC computed against negatives matched on marginal attention (±20%), controlling for expression-driven confounds",
    "n_positive_edges": int(all_labels2.sum()),
    "n_negative_edges": int((1 - all_labels2).sum()),
    "n_total": len(all_labels2),
    "negatives_per_positive": N_MATCHED,
    "matching_tolerance": 0.20,
    "matching_statistic": "attention_marginals (row/col means, proxy for expression)",
    "auroc_attention_mean_all_heads": round(float(auroc_attention), 4),
    "auroc_marginal_product_baseline": round(float(auroc_marginal), 4),
    "auroc_variance_product_baseline": round(float(auroc_variance), 4),
    "auroc_best_single_head": round(float(best_head_auroc), 4),
    "best_head": best_head_name,
    "auroc_median_head": round(float(np.median(list(head_aurocs.values()))), 4),
    "top10_heads": {k: round(v, 4) for k, v in sorted_heads[:10]},
    "interpretation": (
        "Marginal-product baseline is ~0.50 by construction (matched negatives). "
        "If attention AUROC > 0.50, this indicates genuine pairwise structure beyond marginal biases."
    )
}

out_path = os.path.join(r"D:\openclaw\biodyn-nmi-paper\multi_model_proper", "expression_matched_evaluation.json")
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)

print(f"\nResults saved to {out_path}")
print(json.dumps(output, indent=2))
