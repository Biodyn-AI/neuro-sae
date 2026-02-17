import numpy as np
import json

# Check cell types
ct = np.load('activations/cell_types.npy', allow_pickle=True)
print('Cell types shape:', ct.shape)
print('Unique types:', list(np.unique(ct)))
for t in np.unique(ct):
    print(f'  {t}: {(ct==t).sum()}')

# Check activation shapes
for i in [0, 6, 11]:
    cell_path = f'activations/layer_{i:02d}_cell_activations.npy'
    token_path = f'activations/layer_{i:02d}_token_activations.npy'
    ca = np.load(cell_path)
    ta = np.load(token_path)
    print(f'\nLayer {i}: cell_acts={ca.shape}, token_acts={ta.shape}')

# Check metadata summary
with open('activations/metadata.json') as f:
    meta = json.load(f)
print(f'\nn_cells: {meta["n_cells"]}')
print(f'n_genes_in_vocab: {meta["n_genes_in_vocab"]}')
print(f'd_model: {meta["d_model"]}')
print(f'Gene names (first 20): {meta["gene_names"][:20]}')

# Training results
with open('sae_models/training_results.json') as f:
    tr = json.load(f)
for name, r in sorted(tr.items()):
    print(f'{name}: MSE={r["final_mse"]:.6f}, alive={r["alive_fraction"]:.2%}, L0={r["avg_L0"]:.1f}')
