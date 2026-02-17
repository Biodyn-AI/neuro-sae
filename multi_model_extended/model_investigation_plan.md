# Multi-Model Extension Plan

## Target Models for NMI Paper Extension

### 1. scBERT - BERT-based single-cell model
- **Status**: Investigating
- **Installation**: Check pip availability 
- **Input**: Binned gene expression
- **Architecture**: BERT transformer variant

### 2. UCE (Universal Cell Embeddings)
- **Status**: Investigating
- **Source**: Rosen et al. 2024
- **Installation**: Check HuggingFace availability
- **Architecture**: Transformer-based

### 3. scFoundation 
- **Status**: Investigating
- **Source**: Hao et al. 2024
- **Installation**: Check availability
- **Architecture**: Large-scale foundation model

### 4. scVI
- **Status**: Investigating 
- **Architecture**: Variational autoencoder (non-attention baseline)
- **Installation**: Should be available via scvi-tools

## Experiment Plan for Each Working Model

1. **Scaling Behavior Test**
   - 200 vs 500-1000 cells
   - Monitor performance degradation/stability
   - Compare with Geneformer results

2. **Attention/Representation Pattern Extraction**
   - Extract attention weights (for transformer models)
   - Analyze representation patterns (for VAE models)
   - GRN inference quality assessment

3. **Architecture Documentation**
   - Document key differences affecting interpretability
   - Compare tokenization strategies
   - Assess mechanistic interpretability potential

## System Configuration
- Conda env: C:\Users\Agent\miniconda3\envs\bioinfo\python.exe
- GPU: RTX 2060 6GB VRAM
- Memory constraints: Need efficient batching

## Success Criteria
- At least 1-2 additional models running and analyzed
- Complete documentation of all attempts (including failures)
- LaTeX sections and comprehensive report