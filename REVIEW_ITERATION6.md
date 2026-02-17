# Adversarial Review — Iteration 6
**Paper:** A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models  
**Reviewer model:** Claude Opus 4  
**Date:** 2026-02-15  
**Previous ratings:** 5.5 → 6.5 → 7.0 → 7.5 → 7.8  
**Current Rating: 8.1 / 10** (+0.3)

---

## Executive Summary

The three-round-old stratified baseline edge-set mismatch has been resolved — not by computing matched-edge AUROCs, but by the more intellectually honest route of withdrawing the claim entirely. The paper now explicitly states: "We therefore do not claim that attention outperforms stratified correlation; rather, both approaches recover above-chance regulatory signal when properly stratified, with the comparison requiring matched edge sets for definitive conclusions." This is the correct response. The paper crosses 8.0 for the first time.

---

## Assessment of R5 Fixes

### Fix 1 (from R3→R4→R5): Stratified baseline edge-set mismatch — ✅ RESOLVED
The offending claim is withdrawn. The text now includes an explicit bold caveat: **"Crucially, these two AUROCs are not directly comparable"** with clear explanation of why (35 HVG-restricted edges vs. 8,330 full-vocabulary edges, biased subset enrichment). This is the right call — withdrawing an unsupported claim is better than forcing an apples-to-apples comparison that might not be feasible.

### Fix 2: Restructure around core contributions — ❌ NOT DONE
The paper retains its 12-analysis structure. Still reads as a comprehensive audit rather than a focused contribution paper.

### Fix 3: Bootstrap CIs on Table 4 — ❌ NOT DONE
The suspiciously tight 0.010 AUROC range across all layers/tissues in Table 4 still lacks bootstrap CIs.

### Fix 4: Unsupervised CSSI on real data — ❌ NOT DONE
Still only oracle cell-type labels tested on real data. Leiden-based CSSI remains a suggestion without evidence.

### Fix 5: Abstract trimmed — ❌ NOT DONE
Still ~550 words.

---

## Remaining Issues

### 1. **Paper Structure / Identity Crisis** ⚠️ MEDIUM (carried from R5)

The paper is 25+ pages with 12 co-equal analyses. At a top venue, reviewers will ask: "What is the single main contribution?" The three candidates are:
- (a) Attention encodes co-expression, not regulation (mechanistic insight)
- (b) Layer-stratified analysis recovers regulatory signal where naive aggregation fails (practical method)
- (c) Cross-tissue/cross-architecture replication establishing boundary conditions (empirical finding)

Each is publishable. Together, they compete for space and attention. The remaining 9 analyses (pseudotime, batch leakage, calibration, ortholog transfer, etc.) are valuable quality-control checks but dilute the narrative. This is the primary barrier to 8.5+.

### 2. **Table 4 Variance Still Unexplained** ⚠️ LOW-MEDIUM

The 0.010 total AUROC range across 36 cells (12 layers × 3 tissues) for TRRUST, and similar tightness for DoRothEA, remains suspicious without bootstrap CIs. Could be a ceiling effect of the 27-edge evaluation set, or it could be genuine — but the reader can't tell.

### 3. **No Unsupervised CSSI Validation** ⚠️ LOW-MEDIUM

The practical use case for CSSI is when cell-type labels are unavailable. All real-data results use oracle annotations. Leiden clustering is mentioned as "a practical alternative" but never tested. This limits practical impact claims.

### 4. **Abstract Too Long** ⚠️ LOW

550 words. A 300-350 word version focusing on the three core findings would be sharper.

### 5. **Scaling Analysis Model Mismatch** ⚠️ LOW

The detailed 9-point scaling curve is on Geneformer; the positive CSSI/layer results are on scGPT. Table 5 provides coarse scGPT scaling but not at the same resolution. Minor but worth noting.

---

## What Has Improved Since R5

1. **Stratified baseline claim withdrawn** — The oldest unresolved issue is now properly handled
2. **Explicit non-comparability caveat** — Bold text, clear reasoning, no overreach
3. **Honest framing throughout** — The paper consistently under-claims rather than over-claims

---

## Path to 8.5: Prioritized Fixes

| Rank | Fix | Effort | Impact |
|------|-----|--------|--------|
| 1 | Restructure: promote 3 core contributions, demote 9 supporting analyses to "Supplementary Audit" | Medium | High |
| 2 | Add bootstrap CIs to Table 4 | Low | Medium |
| 3 | Run unsupervised CSSI (Leiden) on scGPT-18L brain data | Medium | Medium |
| 4 | Trim abstract to 350 words | Low | Low |

---

## Rating Justification: 8.1/10

**Up from 7.8:** The resolution of the stratified baseline mismatch — by withdrawing the unsupported claim rather than forcing a questionable comparison — demonstrates scientific maturity. Every claim in the paper now has proportionate evidence. The co-expression≠regulation finding is well-supported across architectures. The layer-stratified recovery (AUROC 0.694-0.706 in scGPT-18L) is properly caveated as architecture-specific. The cross-tissue replication (AUROC ~0.72 TRRUST, ~0.84 DoRothEA) is the strongest positive result and is appropriately presented.

**Still below 8.5 because:** The paper's diffuse 12-analysis structure prevents any single contribution from achieving the depth and clarity expected at a top venue. The lack of bootstrap CIs on Table 4 and the absence of unsupervised CSSI validation are addressable gaps. The abstract remains too long.

**For 8.5:** Restructure around core contributions. The current paper reads like a comprehensive technical report; a top-venue paper needs a sharper narrative arc.

**For 9.0+:** Would require either (a) a striking positive result (e.g., CSSI with Leiden labels matching oracle performance on real data), or (b) external validation on Replogle et al. perturbation atlas, or (c) a significantly more focused paper that goes deep on 2-3 findings instead of broad on 12.

**Bottom line:** The paper is now scientifically sound throughout — no claim exceeds its evidence, limitations are honestly stated, and the oldest methodological concern has been properly resolved. The trajectory from 5.5 to 8.1 across six rounds reflects genuine, consistent improvement. The remaining gap to top-venue acceptance is structural (paper organization) rather than scientific (claim validity).
