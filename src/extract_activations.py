"""
Step 1: Extract residual-stream activations from scGPT.

Processes single-cell data through a pre-trained scGPT model and saves
per-layer activations for subsequent SAE training.

Usage:
    python -m src.extract_activations \
        --data data/tabula_sapiens_immune_1000.h5ad \
        --checkpoint external/scGPT_checkpoints/brain/best_model.pt \
        --vocab external/scGPT_checkpoints/brain/vocab.json \
        --scgpt-repo external/scGPT \
        --output-dir outputs/activations \
        --max-cells 1000 \
        --max-genes 200
"""
import argparse
import json
import logging
import sys
import types
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import torch


def setup_scgpt(scgpt_repo: Path):
    """Import scGPT's TransformerModel from the cloned repo."""
    repo_str = str(scgpt_repo)
    if repo_str not in sys.path:
        sys.path.insert(0, repo_str)
    pkg = scgpt_repo / "scgpt"
    stub = types.ModuleType("scgpt")
    stub.__path__ = [str(pkg)]
    stub.__package__ = "scgpt"
    logger = logging.getLogger("scGPT")
    logger.propagate = False
    logger.setLevel(logging.WARNING)
    stub.logger = logger
    sys.modules["scgpt"] = stub
    from scgpt.model import TransformerModel

    return TransformerModel


def main(args):
    import scanpy as sc

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load vocab
    with open(args.vocab) as f:
        vocab = json.load(f)
    pad_token_id = vocab["<pad>"]
    print(f"Vocab: {len(vocab)} tokens")

    # Load model
    TransformerModel = setup_scgpt(Path(args.scgpt_repo))
    model = TransformerModel(
        ntoken=len(vocab),
        d_model=512,
        nhead=8,
        d_hid=512,
        nlayers=12,
        vocab=vocab,
        dropout=0.0,
        pad_token="<pad>",
        pad_value=-2,
        input_emb_style="continuous",
        cell_emb_style="avg-pool",
    )
    state = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    state_dict = (
        state.get("model_state_dict", state.get("model", state))
        if isinstance(state, dict)
        else state
    )
    try:
        from scgpt.utils import load_pretrained

        load_pretrained(model, state_dict, verbose=False)
    except Exception:
        model.load_state_dict(state_dict, strict=False)
    model.eval()
    print("Model loaded (d=512, 12 layers)")

    # Register hooks on all transformer layers
    layer_activations = {i: [] for i in range(12)}

    def make_hook(idx):
        def hook(module, input, output):
            t = output[0] if isinstance(output, tuple) else output
            if hasattr(t, "is_nested") and t.is_nested:
                t = t.to_padded_tensor(0.0)
            layer_activations[idx].append(t.detach().cpu())

        return hook

    handles = [
        layer.register_forward_hook(make_hook(i))
        for i, layer in enumerate(model.transformer_encoder.layers)
    ]

    # Load and preprocess data
    adata = sc.read_h5ad(args.data)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Map genes to vocab (handle Ensembl IDs via feature_name column)
    gene_names = list(adata.var_names)
    valid_genes = [g for g in gene_names if g in vocab]
    if len(valid_genes) < 100:
        for col in ["feature_name", "gene_symbols", "gene_short_name"]:
            if col in adata.var.columns:
                import pandas as pd

                alt = list(adata.var[col].astype(str))
                if sum(1 for g in alt if g in vocab) > len(valid_genes):
                    adata.var_names = pd.Index(alt)
                    adata.var_names_make_unique()
                    valid_genes = [g for g in adata.var_names if g in vocab]
                    break

    gene_to_id = {g: vocab[g] for g in valid_genes}
    adata = adata[:, [g in vocab for g in adata.var_names]].copy()
    gene_name_list = [g for g in adata.var_names if g in gene_to_id]
    gene_id_list = [gene_to_id[g] for g in gene_name_list]
    print(f"Data: {adata.n_obs} cells, {len(valid_genes)} genes in vocab")

    X = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    n_cells = min(adata.n_obs, args.max_cells)

    # Extract cell type labels
    cell_types = []
    for col in ["cell_type_mapped", "cell_type", "celltype"]:
        if col in adata.obs.columns:
            cell_types = list(adata.obs[col].values[:n_cells])
            break

    # Forward pass in batches
    print("Extracting activations...")
    batch_size = 16
    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        batch_genes, batch_values = [], []
        for ci in range(start, end):
            expr = X[ci]
            nonzero = np.where(expr > 0)[0]
            if len(nonzero) > args.max_genes:
                top_idx = nonzero[np.argsort(expr[nonzero])[-args.max_genes :]]
            else:
                top_idx = nonzero
            gids = [gene_id_list[j] for j in top_idx if j < len(gene_id_list)]
            vals = [float(expr[j]) for j in top_idx if j < len(gene_id_list)]
            batch_genes.append(gids)
            batch_values.append(vals)

        max_len = max(len(g) for g in batch_genes)
        if max_len == 0:
            continue
        src = torch.full((len(batch_genes), max_len), pad_token_id, dtype=torch.long)
        values = torch.zeros(len(batch_genes), max_len)
        mask = torch.ones(len(batch_genes), max_len, dtype=torch.bool)
        for i, (gids, vals) in enumerate(zip(batch_genes, batch_values)):
            src[i, : len(gids)] = torch.tensor(gids)
            values[i, : len(gids)] = torch.tensor(vals)
            mask[i, : len(gids)] = False

        with torch.no_grad():
            try:
                model(src, values, mask)
            except Exception:
                model._encode(src, values, mask)

        if start % 200 == 0:
            print(f"  {start}/{n_cells} cells")

    # Save per-layer activations
    print("Saving activations...")
    for layer_idx in range(12):
        if not layer_activations[layer_idx]:
            continue
        all_acts = torch.cat(layer_activations[layer_idx], dim=0)
        cell_acts = all_acts.mean(dim=1).numpy()[:n_cells]
        np.save(output_dir / f"layer_{layer_idx}_cell_acts.npy", cell_acts)

        # Token-level activations for SAE training (subsample to 50k)
        token_acts = all_acts[: min(500, all_acts.shape[0])].reshape(-1, 512).numpy()
        if len(token_acts) > 50000:
            idx = np.random.choice(len(token_acts), 50000, replace=False)
            token_acts = token_acts[idx]
        np.save(output_dir / f"layer_{layer_idx}_token_acts.npy", token_acts)
        print(f"  Layer {layer_idx}: cells={cell_acts.shape}, tokens={token_acts.shape}")

    np.save(output_dir / "cell_types.npy", np.array(cell_types, dtype=str))
    with open(output_dir / "gene_info.json", "w") as f:
        json.dump(
            {
                "gene_names": gene_name_list,
                "gene_to_id": {g: int(v) for g, v in gene_to_id.items()},
            },
            f,
        )

    for h in handles:
        h.remove()
    print(f"Done. Activations saved to {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract scGPT activations")
    parser.add_argument("--data", required=True, help="Path to .h5ad file")
    parser.add_argument("--checkpoint", required=True, help="Path to scGPT best_model.pt")
    parser.add_argument("--vocab", required=True, help="Path to scGPT vocab.json")
    parser.add_argument("--scgpt-repo", required=True, help="Path to cloned scGPT repo")
    parser.add_argument("--output-dir", default="outputs/activations")
    parser.add_argument("--max-cells", type=int, default=1000)
    parser.add_argument("--max-genes", type=int, default=200)
    main(parser.parse_args())
