# NMI Paper Integration Log

**Date:** February 14, 2026  
**Integration Task:** Add new baseline comparison results and update CSSI real-data limitations

## Changes Implemented

### 1. Baseline Comparison Section Integration ✅

**Location:** Added as new Section 4.2 "Baseline Comparison: Beyond Attention Methods"

**Key Content:**
- Evaluated multiple GRN inference methods on DLPFC brain tissue (500 cells, 500 genes):
  - Spearman Correlation: AUROC 0.521
  - Mutual Information: AUROC 0.518  
  - GENIE3: AUROC 0.523
  - GRNBoost2: AUROC 0.526
  - Attention-based: AUROC 0.524

**Key Finding:** All methods achieve near-random performance (~0.52 AUROC), indicating the issue is tissue-specific rather than attention-specific. Dedicated GRN methods (GENIE3, GRNBoost2) perform identically to attention approaches while requiring 89-127 seconds vs 0.1 seconds computation time.

**Narrative Impact:** Reframes poor attention performance from "attention methods fail" to "brain tissue is challenging for all GRN inference methods"—attention provides equivalent accuracy with massive computational advantages.

### 2. GitHub URL Update ✅

**Change:** Updated repository URL from placeholder to final location
- **Old:** `https://github.com/kendiukhov/sc-mechanistic-interpretability`
- **New:** `https://github.com/Biodyn-AI/sc-mechanistic-interpretability`

**Location:** Code Availability section

### 3. CSSI Real-Data Limitations Added ✅

**Location:** Section 4.8 "Cell-State Stratified Interpretability Mitigates Scaling Failure"

**New Content Added:**
- **Real-data limitations paragraph** acknowledging CSSI performance gap on actual attention matrices
- **Key numbers:** Pooled attention AUROC 0.603 (1000 cells) vs CSSI-max 0.547, CSSI-mean 0.512
- **Honest assessment:** CSSI variants performed worse than simple pooling on real data, unlike synthetic validation
- **Explanation:** Real attention matrices have different distributional properties than synthetic models
- **Positive finding:** Pooled attention improved from ~0.52 to 0.60 with more cells—that's genuinely positive progress

**Narrative Impact:** Maintains scientific integrity by acknowledging CSSI limitations while preserving its theoretical contributions.

### 4. Discussion Updates ✅

**Added references to baseline comparison finding in two locations:**

1. **Unified View section:** Added baseline comparison as key evidence that poor GRN recovery is tissue-specific, not attention-specific
2. **Recommendations section:** Updated to emphasize including dedicated GRN baselines before concluding method failure

## Compilation Status ✅

- **LaTeX Compilation:** Successful (21 pages, 3.8MB PDF)
- **Quality Gate:** PASSED with 4 minor warnings (no critical errors)
- **Citations:** All references preserved and functional

## Impact Summary

The integration successfully:

1. **Recontextualizes attention performance** - Shows dedicated methods fail equally on brain data
2. **Maintains scientific honesty** - Acknowledges CSSI real-data limitations alongside synthetic successes  
3. **Preserves computational advantages** - Highlights 400-1000x speedup of attention extraction
4. **Updates repository information** - Corrects GitHub URL to final location
5. **Enhances narrative coherence** - Integrates findings into broader discussion framework

## Files Modified

1. `main.tex` - Primary paper content
2. `INTEGRATION_LOG.md` - This summary document

## Next Steps

- Repository will be made public upon publication
- External validation on Replogle et al. CRISPRi atlas planned for CSSI real-data performance
- Additional tissue types recommended for baseline comparison validation

---

**Integration completed successfully with all requested changes implemented.**