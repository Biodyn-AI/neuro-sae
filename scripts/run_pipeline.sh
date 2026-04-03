#!/usr/bin/env bash
# Run the full analysis pipeline: extract activations, train SAEs, analyze.
#
# Prerequisites:
#   1. Run scripts/download_data.sh first.
#   2. Install dependencies: pip install -r requirements.txt
#
# Usage:
#   bash scripts/run_pipeline.sh

set -euo pipefail
cd "$(dirname "$0")/.."

echo "=== Step 1: Extract scGPT activations ==="
python -m src.extract_activations \
    --data data/tabula_sapiens_immune_1000.h5ad \
    --checkpoint external/scGPT_checkpoints/brain/best_model.pt \
    --vocab external/scGPT_checkpoints/brain/vocab.json \
    --scgpt-repo external/scGPT \
    --output-dir outputs/activations

echo ""
echo "=== Step 2: Train SAEs (main configs) ==="
python -m src.train_saes \
    --activations-dir outputs/activations \
    --output-dir outputs/sae_models \
    --layers 0 6 11 \
    --dict-sizes 2048 \
    --lambdas 1.0 3.0 10.0

echo ""
echo "=== Step 2b: Train SAEs (expansion factor sensitivity) ==="
python -m src.train_saes \
    --activations-dir outputs/activations \
    --output-dir outputs/sae_models \
    --layers 11 \
    --dict-sizes 1024 4096 \
    --lambdas 3.0

echo ""
echo "=== Step 3: Analyze features ==="
python -m src.analyze_features \
    --activations-dir outputs/activations \
    --sae-dir outputs/sae_models \
    --output-dir results

echo ""
echo "=== Pipeline complete ==="
echo "Results in results/"
ls results/*.json
