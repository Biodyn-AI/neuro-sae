# Adversarial Review — Iteration 2
**Paper:** A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models  
**Reviewer model:** Claude Opus 4  
**Date:** 2026-02-15  
**Previous rating:** 5.5/10 (Iteration 1)  
**Current Rating: 6.5 / 10** (+1.0)

---

## Executive Summary

The iteration 1 fixes are largely well-executed. CSSI is now honestly framed as diagnostic, synthetic results are labeled as internal consistency checks, claims-evidence language has been downgraded throughout, and the co-expression vs. regulation finding has been expanded into a proper section (3.12). These changes resolve the most egregious framing problems. However, several substantive issues remain unaddressed, and the fixes themselves introduce a few new problems. The paper is now more honest but still fundamentally limited by single-dataset validation of its strongest claim.

---

## Assessment of Iteration 1 Fixes

### Fix 1: CSSI reframed as diagnostic — ✅ WELL DONE
CSSI is now consistently called a "diagnostic and layer-selection framework" throughout the abstract, methods, results, and discussion. The 1.85× figure is now properly contextualized with "synthetic experiments with oracle cell-state labels." The key sentence "CSSI's primary contribution is diagnostic—identifying these signal-rich layers—rather than improving their already-strong performance" is clear and honest. The abstract still leads with "improving GRN recovery by up to 1.85×" before the diagnostic framing, which slightly undermines the fix, but the caveat follows immediately.

### Fix 2: Synthetic validation relabeled — ✅ WELL DONE
Section 3.11 is now titled "Synthetic Ground-Truth Internal Consistency Checks" and the text explicitly states "these experiments confirm the internal logic of our framework rather than providing independent external validation, since the synthetic generator encodes assumptions..." This is exactly the right framing.

### Fix 3: Claims-evidence language downgraded — ✅ MOSTLY DONE
Section 3.3 is now "Preliminary Evidence of Non-Additivity" (good). Pseudotime section correctly reports null finding after FDR. Perturbation section says "weak, exploratory evidence." However, some residual overclaiming persists (see below).

### Fix 4: Co-expression vs. regulation expanded — ✅ WELL DONE
New Section 3.12 is the strongest addition. It connects co-expression finding to scaling, CSSI, and multi-model results in a coherent mechanistic narrative. The layer-dependent transition from co-expression to regulation encoding is a genuinely interesting insight.

---

## Remaining Critical Issues

### 1. **Single-Dataset Validation of the Core Finding Remains Unaddressed** ⚠️ CRITICAL
**Impact: Very High | Previously: Issue #1, #3**

The L13-L14 AUROC 0.694-0.706 finding still rests entirely on 497 brain cells from one dataset against TRRUST. This was the #1 fix recommendation in iteration 1 ("Validate L13-L14 finding on ≥2 additional tissues with stratified correlation baselines") and it has NOT been done. The cross-validation (248/249 split) is helpful but insufficient—it's still the same 497 cells. The instability noted in iteration 1 persists: the primary analysis finds L13-L14 optimal, the training split finds L17, L16, L13. The paper now acknowledges this ("specific layer rankings show sensitivity to data composition") but acknowledgment ≠ resolution.

**Why this matters:** The entire constructive contribution of the paper (CSSI as diagnostic → later layers contain regulatory signal) depends on this single experiment. One dataset, one tissue, 497 cells.

**Fix:** Run the layer-stratified analysis on at least 2 additional tissues from Tabula Sapiens (immune, kidney—data already in hand). Even if results are weaker, it would establish whether the depth-dependent pattern generalizes.

### 2. **Reference Database Saturation Analysis Remains Circular** ⚠️ HIGH
**Impact: High | Previously: Issue #5**

This was flagged in iteration 1 and NOT fixed. The noise model `σ = 0.05 + 0.002 × cell_count` still linearly increases noise with cell count by construction. The conclusion that "scaling failure is a genuine phenomenon" and "complete references show greater degradation" is entirely an artifact of the noise model. The paper still claims these results "definitively demonstrate" genuine scaling failure—this is the strongest language in the paper applied to its weakest evidence.

**Fix:** Either (a) remove the word "definitively" and reframe as "consistent with" scaling failure under the assumed noise model, or (b) justify the linear noise-cell count relationship from empirical data, or (c) remove the section entirely.

### 3. **Missing Stratified Correlation Baseline Still Not Added** ⚠️ HIGH  
**Impact: High | Previously: Issue #11**

The paper never tests whether simple Spearman correlation *within cell types* achieves comparable AUROC to the layer-stratified attention approach. This is critical because if stratified correlation matches L13-L14 attention AUROC (0.694), then the regulatory signal comes from stratification alone, not from anything attention-specific. The baseline comparison (Section 3.2) uses pooled correlations on brain tissue. A fair test would compute correlations within each of the 7 cell types in the 497-cell brain dataset and compare to Table 5.

**Fix:** Compute stratified Spearman correlation AUROC on the same 497 brain cells, same 7 cell types, same TRRUST edges. One afternoon of work; high impact on interpretation.

### 4. **FDR Correction Still Mixes Synthetic and Real Tests** ⚠️ MEDIUM
**Impact: Medium | Previously: Issue #7**

The Wilcoxon test on synthetic CSSI (p = 2.5 × 10⁻¹¹) is still in the same 47-test FDR family as real-data tests. This inflates the number of "surviving" tests and makes the real-data correction slightly less stringent. The paper acknowledges framework-level correction but doesn't separate families.

**Fix:** Report synthetic and real-data tests in separate FDR families with a brief justification.

---

## New Issues Introduced by Fixes

### 5. **Cross-Validation Section Is Self-Contradictory** ⚠️ MEDIUM-HIGH
**Impact: Medium-High**

The scGPT cross-validation paragraph states: "The top 3 layers from Half A were [0, 4, 10], while Half B yielded [0, 10, 11]." But earlier the paper claims L13-L14 are the best layers. Layers 0, 4, 10, 11 ≠ layers 13-14. The paper tries to reconcile this by saying "later layers (13+) consistently outperform earlier layers" but the cross-validation *literally selected* layers 0, 4, and 10 as top performers—all early/middle layers. This directly contradicts the core claim.

The paper also claims "2/3 layer overlap" between splits as evidence of robustness, but layers [0,4,10] vs [0,10,11] is cherry-picked framing. The overlap with the primary analysis (L13-L14) is 0/3 for both splits.

**Fix:** Either (a) explain why the cross-validation selects different layers than the primary analysis (different evaluation metric? different cell composition? different head aggregation?), or (b) honestly acknowledge that layer selection is unstable and the L13-L14 finding may be dataset-specific.

### 6. **Repetition and Length** ⚠️ MEDIUM
**Impact: Medium**

The fixes added substantial defensive text, making the paper significantly longer. Key findings are now stated 3-4 times each (abstract, results, discussion, conclusions). The CSSI section alone (3.10) is ~2,500 words. The abstract is ~450 words—extremely long. The paper would benefit from a 20-30% cut focusing on removing repetitive caveats and consolidating.

### 7. **"Real-Data-Structured Validation" Framing Is Misleading** ⚠️ MEDIUM
**Impact: Medium**

Section 3.10 includes "Real-data-structured validation" using Tabula Sapiens cell-type proportions but with *synthetic* edge scores. Calling this "real-data-structured" is technically accurate but misleading—a reviewer expecting actual real data will be disappointed. The subsequent "Real attention matrix validation" paragraph is the actual real-data result.

**Fix:** Rename to "Biologically-structured synthetic validation" to avoid confusion with the actual real-data results that follow.

---

## Moderate Issues (Carried Forward)

### 8. **Geneformer Evaluated on Different Tissue Than scGPT** 
**Previously: Issue #9 | Still present**

Table 8 compares scGPT (Tabula Sapiens) with Geneformer (DLPFC brain) side-by-side. The caveat is mentioned but the table presentation still implies direct comparability.

### 9. **Post-Hoc Power Analysis Retained**
**Previously: Issue #12 | Still present**

The perturbation section still includes post-hoc power analysis despite its well-documented limitations.

### 10. **Scope vs. Depth Trade-off**
**Previously: Issue #6 | Partially addressed**

The co-expression section (3.12) helps unify findings, but the paper still has 3-4 analyses yielding null results after correction. The structure hasn't changed.

---

## What Has Improved

1. **Intellectual honesty is now excellent.** The paper is one of the most self-critical I've reviewed. Every major limitation is flagged, often preemptively.
2. **The co-expression finding (Section 3.12)** is now properly developed and serves as a unifying mechanistic narrative. This is the paper's strongest contribution.
3. **CSSI framing** is now appropriate and defensible.
4. **The synthetic validation** is properly caveated and no longer presented as independent evidence.
5. **Language calibration** is much improved throughout—"preliminary evidence," "exploratory," "consistent with" replace prior overclaiming.

---

## Actionable Fixes Ranked by Impact

| Rank | Fix | Effort | Impact |
|------|-----|--------|--------|
| 1 | Validate layer-depth pattern on ≥2 additional tissues | High | Critical |
| 2 | Add stratified-correlation baseline on same 497 cells | Low | High |
| 3 | Fix or remove saturation analysis (circular noise model) | Low | High |
| 4 | Resolve cross-validation contradiction (L0,4,10 vs L13-14) | Low | High |
| 5 | Separate FDR into real vs. synthetic families | Low | Medium |
| 6 | Rename "real-data-structured" → "biologically-structured synthetic" | Low | Medium |
| 7 | Cut 20-30% of repetitive text, especially in abstract and CSSI section | Medium | Medium |
| 8 | Match Geneformer/scGPT tissue for Table 8, or add explicit non-comparability warning | Medium | Medium |
| 9 | Remove post-hoc power analysis | Low | Low |

---

## Rating Justification: 6.5/10

**Up from 5.5: Framing is now honest, the co-expression finding is well-developed, and the paper no longer overclaims CSSI's real-data performance.**

**Still below 7 because:** (1) The core positive finding (L13-L14 regulatory signal) still rests on a single dataset with unstable cross-validation, (2) the saturation analysis remains circular, (3) the cross-validation results actually contradict the primary layer selection, and (4) no stratified-correlation baseline exists to determine whether attention adds value beyond simple stratification.

**For a top venue:** Fixes 1-4 would bring this to ~7.5-8.0. Fix 1 (multi-tissue validation) is the make-or-break item. If L13-L14 depth-dependence replicates across tissues, this is a solid contribution. If it doesn't, the paper becomes a well-executed negative result paper about attention-based GRN inference—still publishable but with different framing.

**Bottom line:** The paper has improved from "overclaimed" to "honestly uncertain." The next step is reducing the uncertainty with additional data, not additional caveats.
