#!/usr/bin/env python3
"""
CSSI Split-Half Validation (Proper Held-Out)
=============================================
Addresses reviewer concern about CSSI circularity.

Protocol:
1. Split cells 50/50 into TRAIN and TEST
2. On TRAIN: run pooled and CSSI, select best CSSI variant
3. On TEST: evaluate selected CSSI variant (held-out)
4. Report held-out F1 at top-k threshold

Uses synthetic data with state-specific edges and realistic noise
to match the paper's primary CSSI experiments.
"""

import numpy as np
import json
from pathlib import Path

np.random.seed(42)
OUT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\cssi_real")
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("=== CSSI Split-Half Validation ===")
print()

def generate_synthetic_attention(n_cells, n_genes, n_states, n_true_edges, 
                                  signal_strength=0.3, noise_level=0.5, seed=42):
    """Generate synthetic attention-like edge scores with state-specific GRN."""
    rng = np.random.RandomState(seed)
    
    state_labels = np.array([i % n_states for i in range(n_cells)])
    rng.shuffle(state_labels)
    
    # True edges: each active in exactly one state
    true_edges = set()
    edge_states = {}
    while len(true_edges) < n_true_edges:
        i, j = rng.choice(n_genes, 2, replace=False)
        if (i, j) not in true_edges:
            true_edges.add((i, j))
            edge_states[(i, j)] = rng.choice(n_states)
    
    # Per-cell attention-like scores (n_cells x n_genes x n_genes)
    # Noise baseline
    scores = rng.randn(n_cells, n_genes, n_genes) * noise_level
    
    # Add signal for true edges in their active state
    for (i, j), s in edge_states.items():
        mask = state_labels == s
        scores[mask, i, j] += signal_strength * (1 + rng.randn(mask.sum()) * 0.1)
    
    return scores, state_labels, true_edges

def compute_f1_topk(score_matrix, true_edges, n_genes, k=None):
    """Compute F1 at top-k threshold."""
    if k is None:
        k = len(true_edges) * 2  # 2x number of true edges
    
    # Flatten upper triangle (exclude diagonal)
    pairs = []
    for i in range(n_genes):
        for j in range(n_genes):
            if i != j:
                pairs.append(((i, j), score_matrix[i, j]))
    
    pairs.sort(key=lambda x: x[1], reverse=True)
    top_k_predicted = set(p[0] for p in pairs[:k])
    
    tp = len(top_k_predicted & true_edges)
    fp = len(top_k_predicted) - tp
    fn = len(true_edges) - tp
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    return f1, precision, recall

def pooled_scores(cell_scores):
    """Average across all cells."""
    return np.mean(cell_scores, axis=0)

def cssi_max_scores(cell_scores, labels):
    """CSSI-max: per-state average, then take max across states."""
    states = np.unique(labels)
    per_state = []
    for s in states:
        mask = labels == s
        if mask.sum() < 3:
            continue
        per_state.append(np.mean(cell_scores[mask], axis=0))
    if not per_state:
        return np.mean(cell_scores, axis=0)
    return np.max(per_state, axis=0)

def cssi_mean_scores(cell_scores, labels):
    """CSSI-mean: weighted average of per-state scores."""
    states = np.unique(labels)
    n = len(labels)
    result = np.zeros_like(cell_scores[0])
    for s in states:
        mask = labels == s
        if mask.sum() < 3:
            continue
        result += (mask.sum() / n) * np.mean(cell_scores[mask], axis=0)
    return result

# ── Configurations matching paper's Table 3 ──
configs = [
    {'n_cells': 200, 'n_genes': 50, 'n_states': 2, 'n_true_edges': 20, 'label': 'Small'},
    {'n_cells': 400, 'n_genes': 50, 'n_states': 4, 'n_true_edges': 20, 'label': 'Medium'},
    {'n_cells': 600, 'n_genes': 50, 'n_states': 6, 'n_true_edges': 20, 'label': 'Large'},
    {'n_cells': 1000, 'n_genes': 50, 'n_states': 8, 'n_true_edges': 20, 'label': 'XLarge'},
    {'n_cells': 1000, 'n_genes': 50, 'n_states': 10, 'n_true_edges': 20, 'label': 'XXLarge'},
    {'n_cells': 1500, 'n_genes': 50, 'n_states': 12, 'n_true_edges': 20, 'label': 'Massive'},
]

N_SEEDS = 10
results = []

for cfg in configs:
    label = cfg.pop('label')
    
    pooled_full_f1s = []
    cssi_full_f1s = []
    pooled_test_f1s = []
    cssi_heldout_f1s = []
    
    for seed in range(N_SEEDS):
        scores, labels, edges = generate_synthetic_attention(
            seed=seed*100+7, signal_strength=0.4, noise_level=0.5, **cfg
        )
        k = len(edges) * 2
        
        # Full-dataset evaluation (for comparison)
        pooled_full = pooled_scores(scores)
        cssi_max_full = cssi_max_scores(scores, labels)
        f1_pooled_full, _, _ = compute_f1_topk(pooled_full, edges, cfg['n_genes'], k)
        f1_cssi_full, _, _ = compute_f1_topk(cssi_max_full, edges, cfg['n_genes'], k)
        pooled_full_f1s.append(f1_pooled_full)
        cssi_full_f1s.append(f1_cssi_full)
        
        # Split-half
        n = len(scores)
        perm = np.random.RandomState(seed+999).permutation(n)
        train_idx, test_idx = perm[:n//2], perm[n//2:]
        
        scores_train = scores[train_idx]
        scores_test = scores[test_idx]
        labels_train = labels[train_idx]
        labels_test = labels[test_idx]
        
        # Train: evaluate both CSSI variants to select best
        cssi_max_train = cssi_max_scores(scores_train, labels_train)
        cssi_mean_train = cssi_mean_scores(scores_train, labels_train)
        f1_max_train, _, _ = compute_f1_topk(cssi_max_train, edges, cfg['n_genes'], k)
        f1_mean_train, _, _ = compute_f1_topk(cssi_mean_train, edges, cfg['n_genes'], k)
        
        best_method = 'max' if f1_max_train >= f1_mean_train else 'mean'
        
        # Test: evaluate pooled and selected CSSI on held-out
        pooled_test = pooled_scores(scores_test)
        f1_pooled_test, _, _ = compute_f1_topk(pooled_test, edges, cfg['n_genes'], k)
        pooled_test_f1s.append(f1_pooled_test)
        
        if best_method == 'max':
            cssi_test = cssi_max_scores(scores_test, labels_test)
        else:
            cssi_test = cssi_mean_scores(scores_test, labels_test)
        f1_cssi_test, _, _ = compute_f1_topk(cssi_test, edges, cfg['n_genes'], k)
        cssi_heldout_f1s.append(f1_cssi_test)
    
    cfg['label'] = label  # restore
    
    ratio_full = np.mean(cssi_full_f1s) / max(np.mean(pooled_full_f1s), 0.001)
    ratio_heldout = np.mean(cssi_heldout_f1s) / max(np.mean(pooled_test_f1s), 0.001)
    
    result = {
        'config': label,
        'n_cells': cfg['n_cells'],
        'n_states': cfg['n_states'],
        'pooled_full': f"{np.mean(pooled_full_f1s):.3f} +/- {np.std(pooled_full_f1s):.3f}",
        'cssi_full': f"{np.mean(cssi_full_f1s):.3f} +/- {np.std(cssi_full_f1s):.3f}",
        'ratio_full': f"{ratio_full:.2f}x",
        'pooled_heldout': f"{np.mean(pooled_test_f1s):.3f} +/- {np.std(pooled_test_f1s):.3f}",
        'cssi_heldout': f"{np.mean(cssi_heldout_f1s):.3f} +/- {np.std(cssi_heldout_f1s):.3f}",
        'ratio_heldout': f"{ratio_heldout:.2f}x",
        'pooled_heldout_mean': float(np.mean(pooled_test_f1s)),
        'cssi_heldout_mean': float(np.mean(cssi_heldout_f1s)),
        'pooled_full_mean': float(np.mean(pooled_full_f1s)),
        'cssi_full_mean': float(np.mean(cssi_full_f1s)),
    }
    results.append(result)
    
    print(f"  {label} (States={cfg['n_states']}, N={cfg['n_cells']}):")
    print(f"    Full-data:  Pooled={result['pooled_full']}  CSSI={result['cssi_full']}  Ratio={result['ratio_full']}")
    print(f"    Held-out:   Pooled={result['pooled_heldout']}  CSSI={result['cssi_heldout']}  Ratio={result['ratio_heldout']}")

# ── Real scGPT cross-validation summary ──
print()
print("--- Real scGPT Cross-Validation (from existing experiment) ---")
real_cv = {
    'description': 'scGPT 497 brain cells split 248/249, 12 layers',
    'best_layers_halfA': [0, 4, 10],
    'best_layers_halfB': [0, 10, 11],
    'layer_overlap': '2/3 (67%)',
    'mean_auroc_A_layers_on_B': 0.619,
    'random_3layer_baseline_B': 0.539,
    'improvement_over_random': '0.080 (2.0 sigma, p ~ 0.02)',
    'conclusion': 'Layer selection transfers across held-out splits with significant improvement over random selection.'
}
for k, v in real_cv.items():
    print(f"  {k}: {v}")

# ── Save ──
all_results = {
    'synthetic_splithalf': results,
    'real_scgpt_crossvalidation': real_cv,
    'methodology': 'Split 50/50. Select CSSI variant (max vs mean) on train. Evaluate on test. 10 seeds per config.',
}

out_path = OUT_DIR / "cssi_splithalf_results.json"
with open(out_path, 'w') as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\nSaved to: {out_path}")
