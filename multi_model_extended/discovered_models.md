# Discovered Single-Cell Foundation Models

## Available Models for Multi-Model Validation

### 1. ✅ scVI - Variational Autoencoder Baseline
- **Status**: Installing (scvi-tools)
- **Type**: Variational autoencoder (non-attention baseline)
- **Strength**: Widely used, well-documented, good for comparison
- **Architecture**: VAE-based, no attention mechanism
- **Use**: Baseline to compare against attention-based models

### 2. ✅ UCE (Universal Cell Embeddings) - HuggingFace Available  
- **Status**: Found on HuggingFace
- **Models Available**: 
  - `minwoosun/uce-100m` (100M parameters)
  - `minwoosun/uce-650m` (650M parameters)  
- **Type**: Transformer-based foundation model
- **Source**: Rosen et al. 2024 implementation
- **Strength**: Multiple sizes, likely has attention patterns

### 3. ✅ C2S-Pythia - Multi-Cell Task Model
- **Status**: Found on HuggingFace
- **Model**: `vandijklab/C2S-Pythia-410m-diverse-single-and-multi-cell-tasks`
- **Type**: Pythia-based transformer (410M parameters)
- **Strength**: Specifically trained on diverse single-cell tasks
- **Architecture**: GPT-style transformer with attention

### 4. ❓ scBERT - Need to investigate further
- **Status**: No direct pip package found
- **Alternative**: May need to implement from scratch or find GitHub repo
- **Type**: BERT-based single-cell model

## Implementation Priority

### Phase 1: Quick Wins (High Success Probability)
1. **scVI** - Install and run baseline comparison  
2. **UCE-100m** - Smaller model, easier to run on 6GB VRAM

### Phase 2: Transformer Comparison  
3. **C2S-Pythia** - Alternative attention-based model
4. **UCE-650m** - Larger model if memory allows

### Phase 3: Research (If Time)
5. **scBERT** - Manual implementation from papers/GitHub

## Technical Constraints
- GPU: RTX 2060 6GB VRAM
- Environment: conda bioinfo with Python 3.10.19
- Memory: Need efficient batching strategies
- Models: Focus on smaller variants first

## Experimental Plan
For each working model:
1. **Basic Functionality Test**: Load model and run inference
2. **Scaling Analysis**: 200 vs 500-1000 cells  
3. **Pattern Extraction**: Attention/representation analysis
4. **GRN Inference**: Compare gene regulatory network quality
5. **Architecture Documentation**: Interpretability differences