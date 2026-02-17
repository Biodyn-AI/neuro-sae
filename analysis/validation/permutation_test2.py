"""Permutation test - TF-restricted evaluation matching original analysis."""
import json, numpy as np, pandas as pd
from scipy import stats
np.random.seed(42)
print("Start", flush=True)

with open(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\outputs\head_layer_checkpoint_replication\gene_list_1200.json") as f:
    gene_list = json.load(f)
gene_to_idx = {g: i for i, g in enumerate(gene_list)}
n_genes = len(gene_list)

trrust = pd.read_csv(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\external\networks\trrust_human.tsv", 
                      sep='\t', header=None, names=['tf','target','mode','pmid'])

# Get all TRRUST TFs that have edges mapping to our gene list
trrust_pos = []  # (tf_idx, target_idx)
trrust_tfs = set()
for _, row in trrust.iterrows():
    if row['tf'] in gene_to_idx and row['target'] in gene_to_idx:
        i, j = gene_to_idx[row['tf']], gene_to_idx[row['target']]
        if i != j:
            trrust_pos.append((i, j))
            trrust_tfs.add(i)
trrust_pos = list(set(trrust_pos))
n_pos = len(trrust_pos)
tf_list = sorted(trrust_tfs)
print(f"TRRUST edges: {n_pos}, TFs: {len(tf_list)}", flush=True)
print(f"TF genes: {[gene_list[t] for t in tf_list]}", flush=True)

# Evaluation universe: all TF -> any gene (excluding self)
# This matches standard GRN AUROC evaluation
eval_pairs = []  # (tf_idx, target_idx)
eval_labels = []
pos_set = set(trrust_pos)
for tf in tf_list:
    for g in range(n_genes):
        if g != tf:
            eval_pairs.append((tf, g))
            eval_labels.append(1 if (tf, g) in pos_set else 0)
eval_labels = np.array(eval_labels)
print(f"Eval universe: {len(eval_pairs)} pairs, {eval_labels.sum()} positive", flush=True)

# Load attention
print("Memory-mapping attention...", flush=True)
attn = np.lib.format.open_memmap(
    r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp\outputs\head_layer_checkpoint_replication\brain\attention_scores_head_layer.npy",
    mode='r')

# Per-layer AUROC
print("Computing per-layer AUROCs...", flush=True)
for L in range(12):
    layer_mean = np.zeros((n_genes, n_genes), dtype=np.float64)
    for H in range(8):
        layer_mean += attn[L, H].astype(np.float64)
    layer_mean /= 8
    
    scores = np.array([layer_mean[i, j] for i, j in eval_pairs])
    
    from sklearn.metrics import roc_auc_score
    auroc = roc_auc_score(eval_labels, scores)
    print(f"  L{L}: AUROC = {auroc:.6f}", flush=True)

# All-layer mean
print("\nAll-layer mean...", flush=True)
mean_attn = np.zeros((n_genes, n_genes), dtype=np.float64)
for L in range(12):
    for H in range(8):
        mean_attn += attn[L, H].astype(np.float64)
mean_attn /= 96

scores_all = np.array([mean_attn[i, j] for i, j in eval_pairs])
obs_auroc = roc_auc_score(eval_labels, scores_all)
print(f"All-layer mean AUROC: {obs_auroc:.6f}", flush=True)

# Permutation test on all-layer mean
# Randomly assign n_pos positives among the eval universe
n_eval = len(eval_pairs)
ranks = scores_all.argsort().argsort().astype(np.float64) + 1
n_neg_eval = n_eval - n_pos

# Observed
pos_indices = np.where(eval_labels == 1)[0]
obs_rank_sum = ranks[pos_indices].sum()
obs_U = obs_rank_sum - n_pos * (n_pos + 1) / 2

# Analytical
sigma = np.sqrt((n_eval + 1) / (12 * n_pos * n_neg_eval))
z = (obs_auroc - 0.5) / sigma
p_anal = 1 - stats.norm.cdf(z)
print(f"\nAnalytical: z={z:.2f}, p={p_anal:.2e}", flush=True)

# Empirical permutation
print("Running 10000 permutations...", flush=True)
n_perms = 10000
null_aurocs = np.zeros(n_perms)
for p in range(n_perms):
    rand_idx = np.random.choice(n_eval, size=n_pos, replace=False)
    rs = ranks[rand_idx].sum()
    null_aurocs[p] = (rs - n_pos*(n_pos+1)/2) / (n_pos * n_neg_eval)
    if (p+1) % 2000 == 0:
        print(f"  {p+1}/{n_perms}", flush=True)

emp_p = (null_aurocs >= obs_auroc).sum() / n_perms
z_emp = (obs_auroc - null_aurocs.mean()) / null_aurocs.std()

print(f"\n=== RESULTS ===")
print(f"Observed AUROC: {obs_auroc:.4f}")
print(f"Null: {null_aurocs.mean():.4f} +/- {null_aurocs.std():.4f}")
print(f"Null range: [{null_aurocs.min():.4f}, {null_aurocs.max():.4f}]")
print(f"Z-score: {z_emp:.2f}")
print(f"Analytical p: {p_anal:.2e}")
print(f"Empirical p: {emp_p}")
if emp_p == 0:
    print(f"P < {1/n_perms}")
print("DONE", flush=True)
