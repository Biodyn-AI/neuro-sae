# Sparse Autoencoders for Single-Cell Foundation Models

Code and data for the paper:

> **Sparse Autoencoders Reveal Interpretable Cell-Type Programs in Single-Cell Foundation Model Representations**
> Ihor Kendiukhov, Department of Computer Science, University of Tubingen
> *Journal of Biomedical Informatics* (under review)

We apply sparse autoencoders (SAEs) — a mechanistic interpretability technique from AI safety — to decompose the internal representations of scGPT into biologically interpretable features aligned with immune cell-type programs.

## Repository Structure

```
neuro-sae/
├── src/                          # Core source code
│   ├── model.py                  #   SAE model definition
│   ├── extract_activations.py    #   Step 1: extract scGPT activations
│   ├── train_saes.py             #   Step 2: train SAEs
│   └── analyze_features.py       #   Step 3: analyze features
├── scripts/                      # Helper scripts
│   ├── download_data.sh          #   Download data & checkpoints
│   ├── prepare_immune_subset.py  #   Build 1000-cell immune subset
│   └── run_pipeline.sh           #   Run full pipeline end-to-end
├── results/                      # Pre-computed analysis results (JSON)
├── figures/                      # Paper figures (PDF)
├── data/                         # Data directory (populated by download script)
├── requirements.txt              # Python dependencies
└── README.md
```

## Quick Start

### 1. Clone and install dependencies

```bash
git clone https://github.com/Biodyn-AI/neuro-sae.git
cd neuro-sae
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Download data and model checkpoints

```bash
bash scripts/download_data.sh
```

This downloads:
- **scGPT source code** from [bowang-lab/scGPT](https://github.com/bowang-lab/scGPT) → `external/scGPT/`
- **scGPT brain checkpoint** (196 MB) from Google Drive → `external/scGPT_checkpoints/brain/`
- **Tabula Sapiens** blood and bone marrow datasets (~3.4 GB) from [CZ CELLxGENE](https://cellxgene.cziscience.com/) → `data/`
- **1,000-cell immune subset** matching the paper's Table 2 composition → `data/tabula_sapiens_immune_1000.h5ad`

### 3. Run the full pipeline

```bash
bash scripts/run_pipeline.sh
```

This runs all three steps sequentially:

1. **Extract activations** — processes 1,000 immune cells through scGPT, saving per-layer residual-stream activations
2. **Train SAEs** — trains sparse autoencoders at layers 0/6/11, λ ∈ {1, 3, 10}, plus expansion factor variants
3. **Analyze features** — computes all paper results: cell-type AUROC, bootstrap CIs, feature ablation, PCA comparison, activation statistics

Results are saved as JSON files in `results/`.

### 4. Run individual steps

Each step can also be run independently with custom arguments:

```bash
# Step 1: Extract activations
python -m src.extract_activations \
    --data data/tabula_sapiens_immune_1000.h5ad \
    --checkpoint external/scGPT_checkpoints/brain/best_model.pt \
    --vocab external/scGPT_checkpoints/brain/vocab.json \
    --scgpt-repo external/scGPT \
    --output-dir outputs/activations

# Step 2: Train SAEs (custom config)
python -m src.train_saes \
    --activations-dir outputs/activations \
    --output-dir outputs/sae_models \
    --layers 11 \
    --dict-sizes 2048 4096 \
    --lambdas 1.0 3.0 10.0

# Step 3: Analyze features
python -m src.analyze_features \
    --activations-dir outputs/activations \
    --sae-dir outputs/sae_models \
    --output-dir results
```

## Pre-computed Results

The `results/` directory contains pre-computed analysis outputs from the paper:

| File | Description |
|------|-------------|
| `pca_analysis.json` | PCA comparison metrics (Table 6 in paper) |
| `bootstrap_ci.json` | Bootstrap 95% CIs per cell type (Table 5) |
| `ablation_results.json` | Feature ablation selectivity (Section 4.5) |
| `hyperparameter_sensitivity.json` | Expansion factor sensitivity (Table 7) |
| `layer_characterization.json` | Per-layer resolved/unresolved cell types (Section 4.2) |
| `activation_statistics.json` | scGPT vs GPT-2 activation distributions (Table 8) |

## Requirements

- Python 3.9+
- PyTorch 2.0+
- ~4 GB disk space for data and checkpoints
- CPU is sufficient (no GPU required)

Full dependency list in `requirements.txt`.

## Citation

```bibtex
@article{kendiukhov2026sparse,
  title={Sparse Autoencoders Reveal Interpretable Cell-Type Programs in
         Single-Cell Foundation Model Representations},
  author={Kendiukhov, Ihor},
  journal={Journal of Biomedical Informatics},
  year={2026}
}
```

## License

This project is for research purposes. The scGPT model and Tabula Sapiens data are subject to their respective licenses.
