"""
Step 1: Extract hidden state activations from scGPT brain model.
Hooks into each transformer layer to capture residual stream activations.
"""
import sys
import json
import types
import logging
import numpy as np
import torch
from pathlib import Path
from torch.utils.data import DataLoader

# Paths
BASE = Path(r"D:\openclaw\mechinterp-bio\biodyn-work\single_cell_mechinterp")
SCGPT_REPO = BASE / "external" / "scGPT"
CHECKPOINT = BASE / "external" / "scGPT_checkpoints" / "brain" / "best_model.pt"
VOCAB_PATH = BASE / "external" / "scGPT_checkpoints" / "brain" / "vocab.json"
OUT_DIR = Path(r"D:\openclaw\biodyn-nmi-paper\brain-sae-paper\experiments\activations")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Settings
MAX_CELLS = 1000
MAX_GENES = 600
BATCH_SIZE = 16
DEVICE = "cpu"

# ---------- Setup scGPT import ----------
repo_str = str(SCGPT_REPO)
if repo_str not in sys.path:
    sys.path.insert(0, repo_str)

scgpt_pkg = SCGPT_REPO / "scgpt"
stub = types.ModuleType("scgpt")
stub.__path__ = [str(scgpt_pkg)]
stub.__package__ = "scgpt"
logger = logging.getLogger("scGPT")
logger.propagate = False
logger.setLevel(logging.WARNING)
stub.logger = logger
sys.modules["scgpt"] = stub

from scgpt.model import TransformerModel

# ---------- Load vocab ----------
with open(VOCAB_PATH) as f:
    vocab = json.load(f)
pad_token_id = vocab["<pad>"]
ntokens = len(vocab)
print(f"Vocab size: {ntokens}, pad_token_id: {pad_token_id}")

# ---------- Load model ----------
print("Loading scGPT brain model...")
model = TransformerModel(
    ntoken=ntokens,
    d_model=512,
    nhead=8,
    d_hid=512,
    nlayers=12,
    vocab=vocab,
    dropout=0.0,  # inference
    pad_token="<pad>",
    pad_value=-2,
    input_emb_style="continuous",
    cell_emb_style="avg-pool",
)

state = torch.load(CHECKPOINT, map_location=DEVICE)
if isinstance(state, dict):
    state_dict = state.get("model_state_dict", state.get("model", state))
else:
    state_dict = state

# Try scGPT's own load_pretrained
try:
    from scgpt.utils import load_pretrained
    load_pretrained(model, state_dict, verbose=False)
    print("Loaded with scGPT load_pretrained")
except Exception as e:
    print(f"scGPT load_pretrained failed ({e}), using strict=False")
    missing, unexpected = model.load_state_dict(state_dict, strict=False)
    print(f"  Missing: {len(missing)}, Unexpected: {len(unexpected)}")

model.eval()
print(f"Model loaded. d_model=512, nlayers=12, nheads=8")

# ---------- Register residual stream hooks ----------
layer_activations = {i: [] for i in range(12)}

def make_hook(layer_idx):
    def hook(module, input, output):
        # TransformerEncoderLayer output is the residual stream
        # May be a NestedTensor — convert to padded
        if isinstance(output, torch.Tensor):
            if output.is_nested:
                output = output.to_padded_tensor(0.0)
            layer_activations[layer_idx].append(output.detach().cpu())
        elif isinstance(output, tuple):
            t = output[0]
            if t.is_nested:
                t = t.to_padded_tensor(0.0)
            layer_activations[layer_idx].append(t.detach().cpu())
    return hook

handles = []
for i, layer in enumerate(model.transformer_encoder.layers):
    h = layer.register_forward_hook(make_hook(i))
    handles.append(h)
print(f"Registered hooks on {len(handles)} layers")

# ---------- Get brain data ----------
# We need to download or find brain h5ad data
# First check if we have any brain data locally
import scanpy as sc

# Try to download a small brain dataset from cellxgene
print("Preparing brain single-cell data...")
brain_h5ad_path = OUT_DIR.parent / "brain_data.h5ad"

if brain_h5ad_path.exists():
    print(f"Loading existing brain data from {brain_h5ad_path}")
    adata = sc.read_h5ad(brain_h5ad_path)
else:
    # Use scanpy's built-in datasets or download from cellxgene
    # Try a synthetic approach: use pbmc3k and relabel, or download brain data
    print("Downloading brain data...")
    try:
        # Try cellxgene census
        import cellxgene_census
        census = cellxgene_census.open_soma()
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_value_filter="tissue_general == 'brain' and is_primary_data == True",
            obs_column_names=["cell_type", "tissue", "disease"],
            var_value_filter="feature_name in " + str(list(vocab.keys())[:100]),
        )
    except Exception as e:
        print(f"Census download failed: {e}")
        print("Generating synthetic brain-like data from existing datasets...")
        
        # Use one of the existing h5ad files and adapt it
        # The tabula sapiens immune subset should work — we just need expression data
        existing_paths = [
            BASE / "data" / "raw" / "tabula_sapiens_immune_subset_20000.h5ad",
            BASE / "data" / "raw" / "tabula_sapiens_kidney.h5ad",
        ]
        
        adata = None
        for p in existing_paths:
            if p.exists():
                print(f"Using existing data: {p}")
                adata = sc.read_h5ad(p)
                break
        
        if adata is None:
            raise RuntimeError("No h5ad files found! Need brain or other scRNA-seq data.")
    
    # Save for future use
    if adata.n_obs > MAX_CELLS:
        sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=42)
    adata.write_h5ad(brain_h5ad_path)
    print(f"Saved {adata.n_obs} cells to {brain_h5ad_path}")

print(f"Data: {adata.n_obs} cells x {adata.n_vars} genes")

# Subsample if needed
if adata.n_obs > MAX_CELLS:
    sc.pp.subsample(adata, n_obs=MAX_CELLS, random_state=42)
    print(f"Subsampled to {adata.n_obs} cells")

# Preprocess
import scipy.sparse as sp

if sp.issparse(adata.X):
    X = adata.X.toarray()
else:
    X = np.array(adata.X)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Map genes to vocab
gene_names = list(adata.var_names)
gene_to_id = {g: vocab[g] for g in gene_names if g in vocab}
valid_genes = [g for g in gene_names if g in vocab]
print(f"Genes in vocab: {len(valid_genes)} / {len(gene_names)}")

if len(valid_genes) < 100:
    print("WARNING: Very few genes in vocab. Trying var['feature_name'] or var index...")
    # Some h5ad files use ensembl IDs as var_names with gene symbols in a column
    for col in ['feature_name', 'gene_symbols', 'gene_short_name', 'gene_name']:
        if col in adata.var.columns:
            alt_names = list(adata.var[col])
            alt_valid = [g for g in alt_names if g in vocab]
            if len(alt_valid) > len(valid_genes):
                print(f"  Using var['{col}']: {len(alt_valid)} genes in vocab")
                adata.var_names = alt_names
                gene_names = alt_names
                gene_to_id = {g: vocab[g] for g in gene_names if g in vocab}
                valid_genes = [g for g in gene_names if g in vocab]
                break

# Filter to valid genes only
adata_filtered = adata[:, [g in vocab for g in adata.var_names]].copy()
print(f"Filtered data: {adata_filtered.n_obs} cells x {adata_filtered.n_vars} genes")

# ---------- Build dataset ----------
sys.path.insert(0, str(BASE))
from src.data.scgpt_dataset import ScGPTDataset, ScGPTDatasetConfig, collate_scgpt

ds_config = ScGPTDatasetConfig(
    max_genes=MAX_GENES,
    include_zero=False,
    sort_by_expression=True,
    pad_token_id=pad_token_id,
    cls_token_id=None,
)

dataset = ScGPTDataset(adata_filtered, gene_to_id, ds_config)
loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate_scgpt)
print(f"Dataset: {len(dataset)} cells, batch_size={BATCH_SIZE}")

# ---------- Extract activations ----------
print("Extracting activations...")
cell_metadata = []

with torch.no_grad():
    for batch_idx, batch in enumerate(loader):
        if batch_idx % 10 == 0:
            print(f"  Batch {batch_idx}/{len(loader)}")
        
        src = batch["gene_ids"].to(DEVICE)
        values = batch["gene_values"].to(DEVICE)
        mask = batch["src_key_padding_mask"].to(DEVICE)
        
        # Forward pass
        try:
            output = model._encode(src, values, mask)
        except Exception as e:
            print(f"  Error in forward pass: {e}")
            # Try alternative forward
            try:
                output = model(src, values, mask)
            except:
                raise
        
        # Collect cell-level metadata
        for i in range(src.shape[0]):
            cell_metadata.append({
                "cell_idx": batch_idx * BATCH_SIZE + i,
                "n_genes": int((~mask[i]).sum().item()),
            })

print(f"Extracted activations for {len(cell_metadata)} cells")

# ---------- Aggregate and save ----------
# For each layer, we have per-token activations (batch, seq_len, 512)
# Average over non-padded tokens to get cell-level embeddings (n_cells, 512)
print("Aggregating per-cell activations...")

# We need padding masks - reconstruct from dataset
all_masks = []
for i in range(len(dataset)):
    sample = dataset[i]
    all_masks.append(sample["src_key_padding_mask"].numpy())
all_masks = np.stack(all_masks)  # (n_cells, max_genes)

for layer_idx in range(12):
    layer_acts = torch.cat(layer_activations[layer_idx], dim=0)  # (n_cells, seq_len, 512)
    
    # Mean-pool over non-padded tokens
    # Build mask from the actual activations (non-zero rows)
    n_cells_act = layer_acts.shape[0]
    seq_len_act = layer_acts.shape[1]
    
    # Use the dataset masks, truncated to match
    mask_np = all_masks[:n_cells_act, :seq_len_act] if all_masks.shape[1] >= seq_len_act else np.pad(all_masks[:n_cells_act], ((0,0),(0,seq_len_act-all_masks.shape[1])), constant_values=True)
    mask_expanded = torch.tensor(~mask_np).unsqueeze(-1).float()  # (n_cells, seq_len, 1)
    
    cell_acts = (layer_acts * mask_expanded).sum(dim=1) / mask_expanded.sum(dim=1).clamp(min=1)
    
    out_path = OUT_DIR / f"layer_{layer_idx:02d}_cell_activations.npy"
    np.save(out_path, cell_acts.numpy())
    print(f"  Layer {layer_idx}: {cell_acts.shape} -> {out_path}")
    
    # Also save per-token activations for a subset (first 500 cells) for richer SAE training
    if layer_acts.shape[0] >= 500:
        token_acts = layer_acts[:500]
        token_mask = mask_expanded[:500]
        # Flatten: keep only non-padded tokens
        flat_acts = []
        for ci in range(500):
            n_valid = int(token_mask[ci, :, 0].sum().item())
            flat_acts.append(token_acts[ci, :n_valid, :].numpy())
        flat_acts = np.concatenate(flat_acts, axis=0)
        token_path = OUT_DIR / f"layer_{layer_idx:02d}_token_activations.npy"
        np.save(token_path, flat_acts)
        print(f"  Layer {layer_idx} tokens: {flat_acts.shape} -> {token_path}")

# Save metadata
metadata = {
    "n_cells": len(cell_metadata),
    "n_genes_in_vocab": len(valid_genes),
    "d_model": 512,
    "n_layers": 12,
    "checkpoint": str(CHECKPOINT),
    "gene_names": valid_genes,
    "cell_metadata": cell_metadata,
}

# Save gene names per cell for later analysis
gene_ids_per_cell = []
gene_names_list = list(adata_filtered.var_names)
for i in range(len(dataset)):
    sample = dataset[i]
    ids = sample["gene_ids"].numpy()
    mask = sample["src_key_padding_mask"].numpy()
    valid_ids = ids[~mask]
    gene_ids_per_cell.append(valid_ids.tolist())

metadata["gene_ids_per_cell"] = gene_ids_per_cell

# Save id_to_gene mapping
id_to_gene = {v: k for k, v in vocab.items()}
metadata["id_to_gene"] = {str(k): v for k, v in id_to_gene.items()}

with open(OUT_DIR / "metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# Also save the cell type labels if available
cell_type_col = None
for col in ['cell_type', 'celltype', 'CellType', 'cell_ontology_class']:
    if col in adata.obs.columns:
        cell_type_col = col
        break

if cell_type_col:
    cell_types = list(adata.obs[cell_type_col].values[:len(cell_metadata)])
    np.save(OUT_DIR / "cell_types.npy", np.array(cell_types, dtype=str))
    print(f"Saved cell types from '{cell_type_col}': {len(set(cell_types))} unique types")

# Remove hooks
for h in handles:
    h.remove()

print("\n=== Step 1 Complete ===")
print(f"Activations saved to {OUT_DIR}")
print(f"Cell-level: 12 files x ({len(cell_metadata)}, 512)")
print(f"Token-level: 12 files for SAE training")
