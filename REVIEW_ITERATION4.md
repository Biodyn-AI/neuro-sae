# Adversarial Review — Iteration 4
**Paper:** A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models  
**Reviewer model:** Claude Opus 4  
**Date:** 2026-02-15  
**Previous ratings:** 5.5 → 6.5 → 7.0  
**Current Rating: 7.5 / 10** (+0.5)

---

## Executive Summary

The cross-tissue replication (Section 3.13, Table 4) is the most consequential addition since iteration 3. AUROC ~0.70 across brain, kidney, and whole-human tissues directly addresses the #1 critical issue carried across three prior reviews ("single-dataset, single-tissue core finding"). The honest reporting that the depth-dependent pattern did NOT replicate in the 12-layer model is commendable and scientifically valuable. However, this non-replication introduces a new tension with the paper's framing, and several issues remain that prevent 8/10.

---

## Assessment of Iteration 3 Fixes

### Fix 1: Cross-tissue replication — ✅ CRITICAL FIX, WELL EXECUTED
**Impact: Resolves the #1 concern from R1–R3.**

AUROC ~0.70 across brain, kidney, and whole-human tissue is exactly the experiment I've been requesting for three rounds. Key strengths:
- Three distinct tissue contexts with consistent performance
- Bootstrap CIs confirm cross-tissue differences are non-significant (Δ = 0.001, 95% CI [−0.003, 0.005])
- Honest reporting that the depth-dependent pattern doesn't replicate

**However**, the replication comes with important caveats that weaken it:

1. **Different architecture than the primary analysis.** The replication uses a 12-layer, 8-head model with 1,200 HVGs and 27 TRRUST edges. The primary CSSI analysis uses an 18-layer model with 1,000 HVGs and 8,330 TRRUST edges. This is a >300× difference in evaluation edge count. The 0.70 AUROC in the 12-layer model and the 0.694 AUROC in the 18-layer model are from different architectures evaluated on different edge sets—their numerical similarity may be coincidental.

2. **27 TRRUST edges is very small.** With 27 positive edges, the standard error on AUROC is ~0.05–0.08. The consistency across tissues (all ~0.70) could reflect the specific properties of these 27 edges rather than a general phenomenon. Are these 27 edges a representative sample of regulatory relationships, or are they biased toward well-characterized, highly conserved interactions that any reasonable method would recover?

3. **The depth-dependent pattern—the paper's most architecturally interesting finding—did NOT replicate.** Brain showed a weak *inverse* trend (ρ = −0.78, p = 0.003), kidney showed no trend, whole-human was flat. This is not a partial replication—it's a non-replication of the depth-dependent pattern. The paper frames this diplomatically as "architecture-specific," but the implication is significant: the L13–L14 recommendation from Section 3.10 may be specific to one architecture and one tissue configuration.

### Fix 2: Abstract and discussion updates — ✅ ADEQUATE
The abstract now includes the cross-tissue replication and the depth-dependent non-replication caveat. The discussion in Section 4.5 (External Validation and Generalization) appropriately references the cross-tissue results as mitigating tissue-specificity concerns.

---

## Remaining Critical Issues

### 1. **The Depth-Dependent Pattern Non-Replication Undermines Core CSSI Narrative** ⚠️ HIGH (new)
**Impact: High**

The paper's most interesting mechanistic claim is that regulatory signal concentrates in later transformer layers (L13–L14), with a monotonic depth-AUROC relationship (Spearman ρ = 0.85 for L6–L17). Section 3.12 (co-expression→regulation transition) builds an elegant narrative around this finding.

But Table 4 shows the 12-layer model has essentially NO depth dependence: all layers achieve ~0.70 ± 0.003. The paper explains this as "architecture-specific" but doesn't adequately grapple with the implications:

- If regulatory signal is uniformly distributed in a 12-layer model but concentrated in L13–L14 in an 18-layer model, the "progressive co-expression→regulation transition" narrative (Section 3.12) is architecture-dependent, not a general principle about how transformers learn regulatory structure.
- The practical recommendation "focus on later layers" (Section 3.10 Decision implication) only applies to the specific 18-layer scGPT architecture on brain tissue. With the 12-layer model, any layer works equally well.
- The CSSI diagnostic framework loses much of its value when all layers are equivalent—there's nothing to diagnose.

**Required fix:** The paper needs to honestly downgrade the depth-dependent claims from "architectural principle" to "architecture-specific observation in one model configuration." The Section 3.12 narrative about progressive co-expression→regulation transition should be explicitly qualified as applying only to the 18-layer model. The practical recommendations should be updated to reflect that layer selection may not matter for all architectures.

### 2. **27-Edge Evaluation Set Is Underpowered for Cross-Tissue Claims** ⚠️ HIGH (new)
**Impact: High**

The cross-tissue replication evaluates against only 27 TRRUST-mapped edges. Compare this to the primary analysis's 8,330 edges—a 308× difference. With 27 edges:

- AUROC has high variance (SE ~0.05–0.08)
- A few consistently high-scoring edges could drive the entire AUROC
- The consistency across tissues could reflect these 27 edges being "easy" (well-conserved, highly expressed TF–target pairs that any method would recover)

The paper should:
- Report which 27 edges these are and whether they're biased toward housekeeping/constitutive regulation
- Report per-edge scores to check if a handful of edges dominate
- Explain why only 27 edges map when the primary analysis achieves 8,330 (is it the 1,200 HVG restriction? If so, can a larger HVG set recover more edges?)

### 3. **Stratified Baseline Comparison Still Not Apples-to-Apples** ⚠️ MEDIUM (carried from R3)
**Impact: Medium**

The 35-edge vs. 8,330-edge mismatch in the stratified correlation baseline comparison was flagged in R3 and remains unaddressed. The 0.09 AUROC gap (stratified correlation 0.604 vs. attention L13–L14 0.694) is evaluated on different edge sets. No bootstrap CI on the AUROC difference is reported.

### 4. **Cross-Validation Contradiction Partially Addressed But Still Problematic** ⚠️ MEDIUM (carried from R3, partially addressed)
**Impact: Medium**

The reconciliation paragraph offers plausible explanations, but the cross-tissue replication actually *deepens* this problem. In the 12-layer model, brain's best layer is L0 (AUROC 0.705), while kidney's best is L6 (0.704). In the 18-layer model, the best layers are L13–L14. And in cross-validation of the 18-layer model, the best layers are L0, L4, L10. There is no stable "regulatory layer"—it changes with architecture, tissue, and sample size. The paper should state this directly rather than presenting the L13–L14 finding as the headline result.

### 5. **Paper Length and Repetition** ⚠️ MEDIUM (carried from R2, R3)
**Impact: Medium**

The paper remains excessively long. The abstract is ~550 words. Key claims are stated 4–5 times. The new Section 3.13 adds ~400 words but could be condensed to a table + 2 paragraphs. The co-expression vs. regulation section (3.12) repeats findings from 3.10 and 3.11 at length.

---

## Minor Issues

### 6. **Two Different scGPT Models Without Clear Labeling**
The paper uses an 18-layer scGPT model (primary CSSI analysis) and a 12-layer scGPT model (cross-tissue replication) without consistently labeling which results come from which model. Table 3 (per-layer AUROC) is from the 18-layer model; Table 4 is from the 12-layer model. A reader could easily conflate the two. Add clear model identifiers (e.g., "scGPT-18L" vs. "scGPT-12L") throughout.

### 7. **Missing Statistical Comparison Between 18-Layer and 12-Layer Models**
The paper presents the 0.694 AUROC (18-layer, 8,330 edges) and 0.70 AUROC (12-layer, 27 edges) as if they're comparable or reinforcing. They're from different models on different edge sets. A formal comparison would require evaluating both models on the same edge set in the same tissue.

### 8. **"Architecture-Specific" Is Under-Theorized**
The paper notes the depth-dependent pattern is "architecture-specific" but offers no mechanistic hypothesis for *why* a 12-layer model distributes regulatory signal uniformly while an 18-layer model concentrates it. Some speculation would be valuable—is it model depth? Training data? Number of heads (8 vs. unknown for the 18-layer model)?

### 9. **Cross-Tissue Causal Analysis Is Underpowered (p = 0.24)**
The cross-tissue causal AUPR analysis (16 edges, 5 positives, permutation p = 0.24) doesn't reach significance. This is reported honestly but adds length without adding insight. Consider condensing or moving to supplement.

---

## What Has Improved Since R3

1. **Cross-tissue replication** — The #1 requested experiment, now executed. AUROC ~0.70 across three tissues is a genuine positive result that substantially strengthens the paper.
2. **Honest depth-dependent non-replication** — The paper could have cherry-picked or downplayed the flat layer profile. Instead it reports it prominently with correlation statistics. This is exemplary scientific practice.
3. **Architecture-dependent framing** — The paper acknowledges the depth pattern is architecture-specific rather than universal, though it could go further (see Issue #1).
4. **Updated abstract and discussion** — The cross-tissue replication and its caveats are integrated throughout.

---

## Path to 8/10: Prioritized Fixes

| Rank | Fix | Effort | Impact | Notes |
|------|-----|--------|--------|-------|
| 1 | Downgrade depth-dependent claims to architecture-specific observation; update Section 3.12 narrative | Low | High | The co-expression→regulation transition story must be qualified |
| 2 | Characterize the 27 evaluation edges (which TFs? constitutive vs. specific?) | Low | High | Critical for interpreting Table 4 |
| 3 | Consistent model labeling (scGPT-18L vs. scGPT-12L) throughout | Low | Medium | Prevents reader confusion |
| 4 | Apples-to-apples stratified baseline (same edges, bootstrap CI on Δ) | Low | Medium | Carried from R3 |
| 5 | Cut 25% of text, especially repetition between 3.10, 3.11, 3.12 | Medium | Medium | Respect reviewer time |

---

## Rating Justification: 7.5/10

**Up from 7.0:** The cross-tissue replication resolves the most critical issue from three prior review rounds. AUROC ~0.70 across brain, kidney, and whole-human tissue is genuine evidence that attention-based regulatory signal is not a tissue-specific artifact. The honest reporting of the depth-dependent non-replication is exemplary. The paper now has a real positive finding (attention captures regulatory structure above chance, consistently across tissues) backed by multi-tissue evidence.

**Still below 8 because:** (1) The depth-dependent non-replication undermines the paper's most interesting mechanistic narrative, and this contradiction isn't fully confronted—Sections 3.10 and 3.12 still read as if depth-dependence is a general architectural principle rather than a model-specific observation. (2) The cross-tissue replication uses only 27 edges, making the AUROC estimates noisy and potentially dominated by a few "easy" edges. (3) The two scGPT models (12-layer vs. 18-layer) are not clearly distinguished, risking conflation of results. (4) The stratified baseline comparison (35 vs. 8,330 edges) remains unresolved. (5) The paper is still 25–30% too long.

**For 8/10:** Fix #1 (honest downgrade of depth-dependent claims) and Fix #2 (characterize the 27 edges) are essential. If the 27 edges include diverse TF–target relationships (not just housekeeping genes), and the depth-dependent claims are properly scoped to the 18-layer model, the paper becomes a strong, honest contribution to the field. Add clear model labels (Fix #3) and the paper is ready for a strong venue.

**For 8.5+:** Additionally resolve the stratified baseline edge-set mismatch with bootstrap CIs, cut 25% of length, and provide some mechanistic hypothesis for why depth-dependence varies across architectures.

**Bottom line:** The paper has progressed from "overclaimed" (5.5) → "honestly uncertain" (6.5) → "one experiment away" (7.0) → "multi-tissue validated but narrative needs updating" (7.5). The cross-tissue replication was the right experiment and the result is positive. The remaining gap to 8/10 is primarily about honest framing: the depth-dependent pattern didn't replicate, and the paper needs to fully absorb that implication rather than treating it as a minor caveat.
