# Adversarial Review Log — brain-sae-paper

## Review 1 — 2026-02-17T16:45 (Initial)

### SCORE: 2/10 (Not publishable — fundamental issues)

### Blocker List

#### BLOCKER-1: DATA IS NOT BRAIN DATA
- Paper claims Allen Human Brain Atlas (~500K cells) and SEA-AD (~1.2M cells)
- **Reality**: Data is Tabula Sapiens immune subset (1000 cells: macrophages, T cells, B cells, NK cells, neutrophils, etc.)
- Only 1 microglial cell out of 1000 is a brain cell type
- The step1_extract_activations.py fallback path loaded immune data instead of brain data
- **This is the single most critical issue**: all biological claims about brain cell types are unfounded
- **Fix**: Must either obtain actual brain data via cellxgene census or reframe paper entirely

#### BLOCKER-2: SAEs ARE NOT SPARSE (L0 = 1200-2400)
- Paper claims target sparsity L0 ≈ 10-50 active features
- **Reality**: Best L0 = 1258 (layer00 d2048 λ=0.1), worst L0 = 2414 (layer06 d4096 λ=0.01)
- For 2048-feature SAEs, 60-94% of features are simultaneously active — this is NOT sparse
- For 4096-feature SAEs, 49-59% active — still not sparse
- The λ values used (0.01, 0.03, 0.1) are far too low
- **Fix**: Retrain with λ ∈ {1.0, 3.0, 10.0, 30.0} to achieve L0 ≈ 10-50

#### BLOCKER-3: ENTIRE RESULTS SECTION IS PLACEHOLDERS
- Every quantitative result uses \placeholder{} macro
- All 5 figures are \fbox{} text boxes, not actual figures
- Table 1 values are explicitly labeled as placeholders
- The paper claims "expected outcomes based on the analysis framework" — no actual analysis
- **Fix**: Run analysis, generate real numbers, create real figures

#### MAJOR-1: ANALYSIS PIPELINE NEVER COMPLETED
- step3_analyze_features.py has been launched ~15 times but never completed
- No all_results.json produced
- The compute_feature_gene_associations function correlates each feature with each of 38,607 genes — extremely slow for 2048+ features
- **Fix**: Optimize analysis code, use vectorized operations, limit to meaningful gene subsets

#### MAJOR-2: GENEFORMER NOT INCLUDED
- Paper discusses Geneformer throughout Methods and Results
- Only scGPT was actually run
- Claims about cross-model consistency are fabricated
- **Fix**: Either include Geneformer or remove all Geneformer claims

#### MAJOR-3: SAMPLE SIZE (N=1000) 
- Paper claims ~500,000 AHBA + ~1.2M SEA-AD cells
- Only 1000 cells extracted
- Token activations: 295,895 (from 1000 cells × ~296 genes average)
- SAE training used 50,000 subsampled tokens
- Too small for credible foundation model analysis
- **Fix**: Increase to at least 10,000-50,000 cells if possible

#### MAJOR-4: d4096 SAEs WORSE THAN d2048
- Higher-dimensional SAEs have higher MSE than lower-dimensional ones
- layer00 d2048 λ=0.01: MSE=0.0016 vs layer00 d4096 λ=0.01: MSE=0.0086
- Suggests d4096 models are undertrained or need different hyperparameters
- **Fix**: More epochs, lower LR, or longer warmup for larger models

#### MAJOR-5: NO STATISTICAL VALIDATION
- Paper describes permutation tests, FDR correction, AUROC analyses
- None have been implemented or run
- No null distributions computed
- **Fix**: Implement actual statistical framework

#### MINOR-1: Companion paper citation is "In preparation" (kendiukhov2024nmi)
#### MINOR-2: SEA-AD citation year mismatch (2023 in text vs 2021 in bib)
#### MINOR-3: Model architecture details have placeholders embedded
#### MINOR-4: No code/data availability statement is concrete
#### MINOR-5: Abstract claims results that don't exist
#### MINOR-6: Discussion speculates without evidence base

---

---

## Review 2 — 2026-02-17T17:30 (Post Iteration 1)

### SCORE: 5/10 (Major improvement; honest paper with real results)

### Changes Made:
1. **BLOCKER-2 FIXED**: Retrained 9 SAEs with λ∈{1,3,10} across layers {0,6,11}
   - λ=10 achieves L0≈48-54 (genuine sparsity)
   - λ=3 achieves L0≈219-251 (moderate sparsity with R²>0.92)
   - λ=1 achieves L0≈553-592 (mild sparsity with R²>0.97)
   
2. **BLOCKER-3 PARTIALLY FIXED**: Replaced all placeholders with real computed values
   - Table 1 has real training stats
   - Cell-type AUROC results populated
   - Gene set enrichment counts are real
   - Figures still TBD (no plotting infrastructure)
   
3. **BLOCKER-1 FIXED**: Completely reframed paper
   - Title changed: no longer claims brain-specific
   - Data honestly described as Tabula Sapiens immune cells
   - Claims matched to evidence: immune cell types, not brain types
   - Removed Geneformer claims (MAJOR-2 fixed)
   - Removed all brain-specific GWAS, synaptic, cortical layer analysis
   
4. **MAJOR-1 FIXED**: Analysis pipeline completes in ~2 min
   
5. **MAJOR-5 PARTIALLY FIXED**: Fisher's exact with BH FDR correction implemented

### Remaining Issues:
- MAJOR: R²_cell is negative for layer 0 (token-trained SAE doesn't transfer well to cell-level)  
- MAJOR: No figures (paper has no visual elements)
- MAJOR: PCA comparison shows PCA is competitive/slightly better than SAE features
- MINOR: Discussion could better discuss the sparsity-specificity tradeoff
- MINOR: More λ values (e.g. 5, 7) could find the optimal point
- MINOR: Need to report token-level R² clearly distinct from cell-level

### Score delta: 2→5 (+3)

## Fix Priority (Impact × Feasibility)

1. **BLOCKER-2** — Fix SAE sparsity (retrain, ~1h CPU)  
2. **BLOCKER-3** — Run analysis, populate real results  
3. **MAJOR-1** — Fix/optimize analysis pipeline  
4. **BLOCKER-1** — Address data issue (reframe or fix)  
5. **MAJOR-2** — Remove Geneformer claims  
6. **MAJOR-5** — Implement statistics  
7. **MAJOR-3** — Address sample size  
8. **MAJOR-4** — Fix d4096 training  
9. **MINOR items** — Cleanup  
