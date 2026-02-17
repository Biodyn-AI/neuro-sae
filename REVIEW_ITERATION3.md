# Adversarial Review — Iteration 3
**Paper:** A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models  
**Reviewer model:** Claude Opus 4  
**Date:** 2026-02-15  
**Previous ratings:** 5.5 → 6.5  
**Current Rating: 7.0 / 10** (+0.5)

---

## Executive Summary

The stratified correlation baseline is the most important addition since iteration 2. The 0.09 AUROC gap (stratified correlation 0.604 vs. attention L13-L14 0.694) is the first real evidence that attention captures something beyond simple co-expression stratified by cell type. The saturation language fix and FDR family separation are adequate. However, several structural problems remain that prevent this from reaching 8/10, and the cross-validation contradiction from iteration 2 has been *explained* but not *resolved*.

---

## Assessment of Iteration 2 Fixes

### Fix 1: Stratified correlation baseline — ✅ CRITICAL FIX, WELL EXECUTED
The comparison (pooled 0.582, stratified CSSI-max 0.595, stratified CSSI-mean 0.604, attention L13-L14 0.694) is exactly what was needed. The 0.09 gap is meaningful and survives the obvious "stratification explains everything" objection. **However**, three caveats weaken the interpretation:

- The stratified baseline uses only 35 mappable TRRUST edges across 13 TFs (stated in text), while the attention analysis uses 8,330 edges. This >200× difference in evaluation set size is never reconciled. Are these the same edges? Different subsets? The comparison is only valid if both methods are evaluated on identical edge sets.
- The stratified baseline uses "absolute Spearman correlations" while attention uses raw attention weights. Different scoring functions on different edge sets ≠ fair comparison.
- N=497 cells, 7 cell types, 35 edges. The confidence interval on a 0.09 AUROC difference with 35 positive edges is enormous. No bootstrap CI is reported for this specific comparison.

**Verdict:** Directionally correct and important, but the comparison needs tightening before it can bear the interpretive weight placed on it.

### Fix 2: Saturation analysis language — ✅ ADEQUATE
The caveat about the linear noise-cell count assumption is now present. The word "definitively" appears to have been softened to "consistent with." This resolves the circularity objection at the language level, though the analysis itself remains constructed.

### Fix 3: FDR families — ✅ ADEQUATE
The paper now acknowledges that synthetic and real tests are conceptually distinct families and notes that separating them doesn't change which real-data results survive. This is honest and sufficient.

### Fix 4: Cross-validation layer instability — PARTIALLY ADDRESSED ⚠️
A new paragraph ("Reconciling cross-validation and primary analysis layer rankings") offers two explanations: (a) different evaluation criteria (transfer vs. within-sample AUROC), and (b) reduced sample size per fold (~248 cells, below the 300-cell sensitivity region). Both are plausible. But:

- The cross-validation top layers are [0, 4, 10] and [0, 10, 11]. The primary analysis optimum is L13-L14. The paper claims "later layers (13+) consistently outperform earlier layers" but the cross-validation *selected* L0 as a top layer in both folds. Layer 0 is the earliest possible layer.
- The reconciliation paragraph implicitly concedes that the L13-L14 finding may be sample-size-dependent. If 248 cells yield different optimal layers than 497 cells, how robust is the 497-cell result? Would 1,000 cells yield yet another optimum?

---

## Remaining Critical Issues

### 1. **Single-Dataset, Single-Tissue Core Finding** ⚠️ CRITICAL (carried from R1, R2)
**Impact: Very High**

This remains the #1 problem and has not been addressed across three review rounds. The L13-L14 AUROC 0.694-0.706 finding comes from:
- 1 dataset (Tabula Sapiens)
- 1 tissue (brain/DLPFC)
- 497 cells
- 7 cell types
- Against TRRUST (which the paper repeatedly notes has poor brain-specific coverage)

The paper has immune and kidney data already in hand. Running the layer-stratified analysis on these tissues would either (a) replicate the depth-dependent pattern and dramatically strengthen the paper, or (b) fail to replicate, which is equally important to know. The persistent avoidance of this obvious validation is the single biggest weakness.

**Required fix:** Layer-stratified AUROC analysis on ≥1 additional tissue. This is non-negotiable for a strong venue.

### 2. **Stratified Baseline Comparison Is Not Apples-to-Apples** ⚠️ HIGH (new)
**Impact: High**

As noted above, the stratified correlation baseline uses 35 TRRUST edges while the attention analysis uses 8,330. The paper states this explicitly ("35 mappable TRRUST edges across 13 TFs") but never addresses the implications:

- With 35 edges, AUROC variance is ~0.05-0.08 (standard for small binary classification). The 0.09 gap may not be statistically significant.
- The 8,330 edges used for attention analysis likely include many edges with zero expression in these 497 brain cells, deflating the baseline but not attention (which can assign non-zero weights regardless of expression).
- **Fix:** Report both methods on the identical edge set with bootstrap CIs on the AUROC difference.

### 3. **Cross-Validation Contradicts Primary Finding More Than Acknowledged** ⚠️ HIGH (carried, partially addressed)
**Impact: High**

The reconciliation paragraph is reasonable but incomplete. The fundamental problem: if 248 cells select L0 as a top layer and 497 cells select L13, the layer selection is a function of sample size, not a stable architectural property. The paper frames depth-dependence as a "consistent architectural principle" but the cross-validation evidence directly contradicts this.

The proposed explanation (scaling dynamics favor different layers at different N) is actually a *stronger* objection than the original instability concern—it means the L13-L14 recommendation is specific to N≈500 and may not generalize to other dataset sizes practitioners might use.

**Fix:** Either (a) run cross-validation with larger folds (e.g., 80/20 instead of 50/50 to keep N>350 per fold), or (b) honestly state that layer selection is N-dependent and provide guidelines for different sample sizes.

### 4. **The Paper Is Too Long and Repetitive** ⚠️ MEDIUM (carried from R2)
**Impact: Medium**

The abstract is ~500 words. Key findings are stated 4+ times (abstract, results, discussion, conclusions). The CSSI section alone is ~3,000 words. The paper has 12 analyses, many yielding null or weak results. A reviewer at a top venue will note the low signal-to-noise ratio in the paper itself.

Conservative estimate: 30% could be cut without losing substance. The ortholog transfer, pseudotime, and batch leakage sections are thorough but could be condensed to 1-2 paragraphs each in a "boundary condition analyses" omnibus section.

### 5. **Geneformer and scGPT Evaluated on Different Tissues** ⚠️ MEDIUM (carried)
**Impact: Medium**

Table 8 compares scGPT (Tabula Sapiens) with Geneformer (DLPFC brain). The caveat is mentioned but the side-by-side table presentation implies comparability. The "convergent failure" claim is weakened when the two models are tested on different data.

---

## Minor Issues

### 6. **"Real-data-structured validation" naming**
Still misleading. These are synthetic data with realistic cell-type proportions. Call them "biologically-structured synthetic" experiments.

### 7. **Post-hoc power analysis retained**
The perturbation section still includes post-hoc power analysis, which is methodologically questionable and adds length.

### 8. **Custom synthetic generator vs. SERGIO**
The justification for not using SERGIO is reasonable but runs ~200 words. The core argument ("SERGIO doesn't model attention matrices") could be stated in one sentence.

### 9. **Missing positive controls**
The limitations section honestly acknowledges this, which is good. But for a paper making primarily negative claims, the absence of positive controls means every null finding is ambiguous: is the method failing, or is the evaluation framework failing?

---

## What Has Improved Since R2

1. **Stratified correlation baseline** — The single most important addition. Even with caveats, it provides the first evidence that attention layers add value beyond stratified correlations.
2. **Saturation language** — No longer overclaims circularity analysis.
3. **FDR honesty** — Acknowledges the synthetic/real family issue without pretending it's fully resolved.
4. **Cross-validation explanation** — Plausible, even if incomplete.

---

## Path to 8/10: Prioritized Fixes

| Rank | Fix | Effort | Impact | Notes |
|------|-----|--------|--------|-------|
| 1 | Layer-stratified analysis on ≥1 additional tissue | Medium | Critical | Data in hand (immune, kidney). Make-or-break. |
| 2 | Apples-to-apples stratified baseline (same edges, bootstrap CI) | Low | High | Tighten the 0.09 gap claim. |
| 3 | Larger cross-validation folds (80/20) or N-dependent guidance | Low | High | Resolve the L0-vs-L13 contradiction. |
| 4 | Cut 25-30% of text, especially repetition and null-result sections | Medium | Medium | Respect reviewer time. |
| 5 | Match tissues for multi-model comparison | Medium | Medium | Run Geneformer on Tabula Sapiens or scGPT on DLPFC. |

---

## Rating Justification: 7.0/10

**Up from 6.5:** The stratified correlation baseline (0.09 AUROC gap) provides genuinely new evidence that attention captures regulatory structure beyond co-expression. The saturation and FDR fixes close minor holes. The paper is now intellectually honest throughout.

**Still below 8 because:** (1) The core positive finding still rests on a single tissue/dataset despite data for additional tissues being available—this is now a three-round-old unaddressed critical issue. (2) The stratified baseline comparison has a >200× edge set size mismatch that undermines the 0.09 gap. (3) Cross-validation selects fundamentally different layers (L0 vs L13), and the reconciliation is a concession that layer selection is N-dependent rather than a resolution. (4) The paper is 30% too long for its signal content.

**For 8/10:** Fix #1 (multi-tissue validation) is essential. If the depth-dependent pattern replicates on immune tissue—even partially (e.g., later layers outperform earlier by 0.05+ AUROC)—the paper becomes a genuinely strong contribution. Fix #2 (matched edge sets + CI) would solidify the attention-vs-correlation comparison. These two fixes together would bring the rating to 7.5-8.0.

**For 8.5+:** Additionally cut length by 25%, resolve the cross-validation contradiction with larger folds, and match tissues across models.

**Bottom line:** The paper has progressed from "overclaimed" (5.5) through "honestly uncertain" (6.5) to "one experiment away from strong" (7.0). The stratified baseline was the right move. Now validate the core finding on a second tissue.
