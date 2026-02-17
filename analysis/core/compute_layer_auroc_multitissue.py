"""
Compute per-layer AUROC against TRRUST and DoRothEA for brain, kidney, whole_human tissues.
Uses pre-computed attention_scores_head_layer.npy from head_layer_checkpoint_replication.
"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score

BASE = Path(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp")
ATLAS_DIR = BASE / "outputs" / "head_layer_checkpoint_replication"
GENE_LIST = ATLAS_DIR / "gene_list_1200.json"

# Load gene list
with open(GENE_LIST) as f:
    genes = json.load(f)
n_genes = len(genes)
gene_to_idx = {g: i for i, g in enumerate(genes)}
print(f"Loaded {n_genes} genes")

# Load ground truth
def load_trrust():
    """Load TRRUST v2 from the project's data."""
    trrust_path = BASE / "data" / "ground_truth" / "trrust_rawdata.human.tsv"
    if not trrust_path.exists():
        # Try alternative paths
        for p in [BASE / "data" / "trrust_rawdata.human.tsv",
                  BASE / "external" / "trrust_rawdata.human.tsv"]:
            if p.exists():
                trrust_path = p
                break
    if not trrust_path.exists():
        print(f"TRRUST not found, searching...")
        import subprocess
        result = subprocess.run(
            ["powershell", "-Command", 
             f'Get-ChildItem -Path "{BASE}" -Recurse -File -Include "*trrust*" | Select-Object FullName'],
            capture_output=True, text=True
        )
        print(result.stdout)
        return None
    
    df = pd.read_csv(trrust_path, sep='\t', header=None, 
                     names=['source', 'target', 'mode', 'pmid'])
    edges = set()
    for _, row in df.iterrows():
        s, t = row['source'], row['target']
        if s in gene_to_idx and t in gene_to_idx:
            edges.add((gene_to_idx[s], gene_to_idx[t]))
    print(f"TRRUST: {len(edges)} edges in vocabulary")
    return edges

def load_dorothea():
    """Load DoRothEA from the project."""
    # Try to find dorothea file
    import subprocess
    result = subprocess.run(
        ["powershell", "-Command",
         f'Get-ChildItem -Path "{BASE}" -Recurse -File -Include "*dorothea*","*DoRothEA*" | Select-Object FullName'],
        capture_output=True, text=True
    )
    files = [l.strip() for l in result.stdout.strip().split('\n') if l.strip() and 'FullName' not in l and '---' not in l]
    
    for fpath in files:
        fpath = fpath.strip()
        if not fpath or not Path(fpath).exists():
            continue
        print(f"Trying DoRothEA file: {fpath}")
        try:
            if fpath.endswith('.tsv'):
                df = pd.read_csv(fpath, sep='\t')
            elif fpath.endswith('.csv'):
                df = pd.read_csv(fpath)
            else:
                continue
            # Try to extract source-target pairs
            if 'source' in df.columns and 'target' in df.columns:
                edges = set()
                for _, row in df.iterrows():
                    s, t = str(row['source']), str(row['target'])
                    if s in gene_to_idx and t in gene_to_idx:
                        edges.add((gene_to_idx[s], gene_to_idx[t]))
                if edges:
                    print(f"DoRothEA: {len(edges)} edges in vocabulary from {fpath}")
                    return edges
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    # Try loading from the src module
    try:
        sys.path.insert(0, str(BASE))
        from src.eval.dorothea import load_dorothea as _load_dorothea
        df = _load_dorothea()
        edges = set()
        for _, row in df.iterrows():
            s, t = str(row['source']), str(row['target'])
            if s in gene_to_idx and t in gene_to_idx:
                edges.add((gene_to_idx[s], gene_to_idx[t]))
        print(f"DoRothEA (via module): {len(edges)} edges in vocabulary")
        return edges
    except Exception as e:
        print(f"Could not load DoRothEA via module: {e}")
    
    return None


def compute_layer_auroc(attn_scores, positive_edges, n_genes, max_eval_pairs=50000, seed=42):
    """
    Compute per-layer AUROC.
    attn_scores: shape (n_layers, n_heads, n_genes, n_genes) or (n_layers*n_heads, n_genes, n_genes)
    """
    if attn_scores.ndim == 3:
        # Reshape: assume it's (n_layers*n_heads, n_genes, n_genes)
        total = attn_scores.shape[0]
        # scGPT has 12 layers, ~n heads
        # Try to infer
        for n_layers in [12, 6, 24]:
            if total % n_layers == 0:
                n_heads = total // n_layers
                attn_scores = attn_scores.reshape(n_layers, n_heads, n_genes, n_genes)
                print(f"  Reshaped to ({n_layers}, {n_heads}, {n_genes}, {n_genes})")
                break
        else:
            print(f"  Cannot reshape {total} into layers*heads")
            return None
    
    n_layers = attn_scores.shape[0]
    n_heads = attn_scores.shape[1]
    
    # Build positive set
    pos_list = list(positive_edges)
    pos_set = positive_edges
    
    if len(pos_list) < 5:
        print(f"  Too few positive edges ({len(pos_list)}), skipping")
        return None
    
    # Sample negative edges
    rng = np.random.default_rng(seed)
    n_neg = min(len(pos_list) * 10, max_eval_pairs - len(pos_list))
    neg_pairs = set()
    attempts = 0
    while len(neg_pairs) < n_neg and attempts < n_neg * 100:
        i = rng.integers(0, n_genes)
        j = rng.integers(0, n_genes)
        if i != j and (i, j) not in pos_set and (i, j) not in neg_pairs:
            neg_pairs.add((i, j))
        attempts += 1
    neg_list = list(neg_pairs)
    
    all_pairs = pos_list + neg_list
    labels = np.array([1]*len(pos_list) + [0]*len(neg_list))
    src_idx = np.array([p[0] for p in all_pairs])
    tgt_idx = np.array([p[1] for p in all_pairs])
    
    print(f"  Evaluating {len(pos_list)} pos + {len(neg_list)} neg = {len(all_pairs)} pairs")
    
    results = {}
    
    # Per-layer AUROC (aggregating heads by mean)
    for layer in range(n_layers):
        # Mean across heads for this layer
        layer_scores = attn_scores[layer].mean(axis=0)  # (n_genes, n_genes)
        scores = layer_scores[src_idx, tgt_idx]
        
        if np.std(scores) == 0:
            results[f"L{layer}"] = 0.5
            continue
        
        try:
            auroc = roc_auc_score(labels, scores)
            results[f"L{layer}"] = auroc
        except:
            results[f"L{layer}"] = 0.5
    
    # Per-head AUROC for best layer
    best_layer = max(range(n_layers), key=lambda l: results.get(f"L{l}", 0.5))
    for head in range(n_heads):
        head_scores = attn_scores[best_layer, head]
        scores = head_scores[src_idx, tgt_idx]
        if np.std(scores) == 0:
            results[f"L{best_layer}_H{head}"] = 0.5
            continue
        try:
            auroc = roc_auc_score(labels, scores)
            results[f"L{best_layer}_H{head}"] = auroc
        except:
            results[f"L{best_layer}_H{head}"] = 0.5
    
    # All-layer mean
    all_layer_scores = attn_scores.mean(axis=(0, 1))  # (n_genes, n_genes)
    scores = all_layer_scores[src_idx, tgt_idx]
    try:
        results["all_layers_mean"] = roc_auc_score(labels, scores)
    except:
        results["all_layers_mean"] = 0.5
    
    return results


# Load ground truth
print("\n=== Loading ground truth ===")
trrust_edges = load_trrust()
dorothea_edges = load_dorothea()

# Process each tissue
tissues = ["brain", "kidney", "whole_human"]
all_results = {}

for tissue in tissues:
    tissue_dir = ATLAS_DIR / tissue
    attn_path = tissue_dir / "attention_scores_head_layer.npy"
    
    if not attn_path.exists():
        print(f"\n=== {tissue}: attention scores not found ===")
        continue
    
    print(f"\n=== Processing {tissue} ===")
    print(f"  Loading {attn_path} ({attn_path.stat().st_size / 1e6:.1f} MB)...")
    attn = np.load(attn_path)
    print(f"  Shape: {attn.shape}, dtype: {attn.dtype}")
    
    tissue_results = {}
    
    if trrust_edges:
        print(f"\n  --- TRRUST ---")
        res = compute_layer_auroc(attn, trrust_edges, n_genes)
        if res:
            tissue_results["trrust"] = res
            # Print layer results
            layer_keys = sorted([k for k in res if k.startswith("L") and "_H" not in k], 
                              key=lambda x: int(x[1:]))
            print(f"  Per-layer AUROC (TRRUST):")
            for k in layer_keys:
                marker = " ***" if res[k] == max(res[lk] for lk in layer_keys) else ""
                print(f"    {k}: {res[k]:.4f}{marker}")
            print(f"    All-layers mean: {res.get('all_layers_mean', 'N/A')}")
    
    if dorothea_edges:
        print(f"\n  --- DoRothEA ---")
        res = compute_layer_auroc(attn, dorothea_edges, n_genes)
        if res:
            tissue_results["dorothea"] = res
            layer_keys = sorted([k for k in res if k.startswith("L") and "_H" not in k],
                              key=lambda x: int(x[1:]))
            print(f"  Per-layer AUROC (DoRothEA):")
            for k in layer_keys:
                marker = " ***" if res[k] == max(res[lk] for lk in layer_keys) else ""
                print(f"    {k}: {res[k]:.4f}{marker}")
            print(f"    All-layers mean: {res.get('all_layers_mean', 'N/A')}")
    
    all_results[tissue] = tissue_results
    del attn  # Free memory

# Save results
output_path = Path(r"D:\openclaw\biodyn-nmi-paper\multi_tissue_layer_auroc.json")
# Convert to serializable
serializable = {}
for tissue, tres in all_results.items():
    serializable[tissue] = {}
    for ref, res in tres.items():
        serializable[tissue][ref] = {k: float(v) for k, v in res.items()}

with open(output_path, 'w') as f:
    json.dump(serializable, f, indent=2)
print(f"\nResults saved to {output_path}")

# Summary table
print("\n\n=== SUMMARY TABLE ===")
print(f"{'Tissue':<15} {'Reference':<12} {'Best Layer':<12} {'Best AUROC':<12} {'All-Mean':<12}")
print("-" * 63)
for tissue, tres in all_results.items():
    for ref, res in tres.items():
        layer_keys = sorted([k for k in res if k.startswith("L") and "_H" not in k],
                          key=lambda x: int(x[1:]))
        if layer_keys:
            best_k = max(layer_keys, key=lambda k: res[k])
            best_v = res[best_k]
            all_mean = res.get("all_layers_mean", float("nan"))
            print(f"{tissue:<15} {ref:<12} {best_k:<12} {best_v:<12.4f} {all_mean:<12.4f}")
