# Adversarial Review — Iteration 5
**Paper:** A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models  
**Reviewer model:** Claude Opus 4  
**Date:** 2026-02-15  
**Previous ratings:** 5.5 → 6.5 → 7.0 → 7.5  
**Current Rating: 7.8 / 10** (+0.3)

---

## Executive Summary

The paper has addressed four of five priority fixes from R4. Depth-dependent claims are now properly scoped as "architecture-specific." The 27 TRRUST edges are characterized (14 TFs including STAT3, TP53, MYC, SP1, etc.). Model labeling is consistent (scGPT-18L/scGPT-12L throughout). The paper is trimmed to ~25 pages. The addition of DoRothEA cross-tissue results (~0.84 AUROC, 123 edges) substantially strengthens the cross-tissue replication by providing a larger evaluation set. The remaining gap to 8.0+ is no longer about missing experiments—it's about structural clarity, one unresolved methodological issue, and the paper's identity crisis between "comprehensive negative results" and "CSSI as constructive contribution."

---

## Assessment of R4 Fixes

### Fix 1: Depth-dependent claims downgraded — ✅ DONE WELL
Throughout the paper, the depth-dependent L13–L14 finding is now consistently described as "architecture-specific" with explicit caveats that it applies only to scGPT-18L. Section 3.12 (co-expression→regulation transition) is properly qualified. The decision implications in Section 3.10 now recommend "layer-stratified evaluation for each new architecture" rather than "focus on later layers." This is exactly what was needed.

### Fix 2: 27-edge characterization — ✅ DONE WELL
The 14 TFs are now named (STAT3, TP53, NFKB1, SP1, MYC, RELA, JUN, EGR1, CREB1, HIF1A, FOS, E2F1, ETS1, GATA2) with explicit commentary that they span ubiquitous regulators and context-specific factors. This substantially mitigates the concern about "easy" housekeeping edges dominating.

### Fix 3: Consistent model labeling — ✅ DONE
scGPT-18L and scGPT-12L used consistently throughout. No more ambiguity.

### Fix 4: Paper trimmed — ✅ PARTIALLY DONE
Down from ~29 to ~25 pages. Still long for what it delivers (see Issue #3 below), but the reduction is real.

### Fix 5 (from R3): Apples-to-apples stratified baseline — ❌ NOT ADDRESSED
The 35-edge vs. 8,330-edge mismatch in the stratified correlation comparison remains. The text still says "this comparison involves different edge sets (35 vs. 8,330 edges)" without resolving it. No bootstrap CI on the AUROC difference. This was flagged in R3 and R4.

### NEW: DoRothEA cross-tissue results — ✅ EXCELLENT ADDITION
The DoRothEA evaluation (123 edges, ~0.84 AUROC across all tissues) is the single most impactful change in this iteration. It:
- Provides 4.5× more evaluation edges than the 27-edge TRRUST set
- Achieves higher AUROC with tighter cross-tissue agreement
- Substantially strengthens the claim that regulatory signal recovery is tissue-invariant
- ChIP-seq-derived edges are more robust ground truth than literature-curated TRRUST

This directly addresses the R4 concern about the 27-edge set being underpowered.

---

## Remaining Issues Preventing 8.0

### 1. **The Stratified Baseline Edge-Set Mismatch Is Now the Oldest Unresolved Issue** ⚠️ MEDIUM-HIGH (carried R3→R4→R5)

The claim that "deep attention layers capture regulatory structure beyond what stratified correlations recover" (the ~0.09 AUROC gap) rests on comparing 35-edge stratified correlation AUROC against 8,330-edge attention AUROC. These are different evaluation sets—the gap could be entirely driven by edge-set composition. After three rounds of flagging, this needs to be either:
- (a) Resolved by evaluating both methods on the same edge set, or
- (b) Explicitly acknowledged as an unresolvable limitation and the claim withdrawn

This matters because the stratified correlation baseline is the key comparison that justifies attention's value over simple methods. If the ~0.09 gap is an artifact of edge-set mismatch, attention offers no advantage over cell-type-stratified Spearman correlation—which would be a significant reframing.

### 2. **The Paper's Identity Crisis: Comprehensive Audit vs. CSSI Method Paper** ⚠️ MEDIUM

The paper tries to be two things simultaneously:
1. A comprehensive negative-results audit of attention-based GRN inference (12 analyses, mostly showing limitations)
2. A methods paper introducing CSSI as a constructive contribution

These two identities compete. The 12-analysis audit is thorough but diffuse—no single finding gets the depth it deserves. CSSI is the constructive contribution but its real-data impact is modest (diagnostic only; +0.033 AUROC at best on real attention matrices). The cross-tissue replication (Section 3.13) is arguably the strongest positive result but it's buried as analysis #13 of 13.

A top venue reviewer would likely say: "This paper reports many things but advances no single thing deeply enough." The mediation bias analysis has N=16. The pseudotime analysis doesn't survive FDR. The perturbation validation doesn't survive FDR. The cross-tissue consistency has limited pairs. Each analysis is interesting but underpowered individually.

**Suggestion:** The paper would be stronger if restructured around its three genuinely novel contributions: (1) attention encodes co-expression not regulation, (2) layer-stratified analysis recovers regulatory signal where naive aggregation fails, (3) cross-tissue/cross-architecture replication establishing boundary conditions. The other 8-9 analyses could be condensed into a "supplementary audit" section.

### 3. **CSSI's Real-Data Value Remains Primarily Diagnostic** ⚠️ MEDIUM

On real scGPT-18L attention matrices:
- Best CSSI improvement: +0.033 at L10 (modest)
- At the top layers (L13–L14): CSSI provides zero improvement (Δ ≤ 0.001)
- The paper honestly states "CSSI's primary value on real data is diagnostic"

On scGPT-12L:
- All layers are equivalent → nothing to diagnose
- CSSI adds no value

The synthetic results (1.85× improvement) are compelling but use oracle cell-state labels and a custom generator designed to validate the framework. The gap between synthetic promise and real-data impact is significant. A reviewer could reasonably argue: "CSSI is a stratification heuristic whose main contribution is telling you which layers to look at—a useful observation, but not a major methodological advance."

### 4. **Scaling Analysis Uses Geneformer, Not scGPT** ⚠️ LOW-MEDIUM

The detailed 9-point scaling analysis (Section 3.1) is conducted on Geneformer, while the core positive results (CSSI, layer-stratified AUROC) are on scGPT. The paper doesn't present equivalent scaling curves for scGPT. Table 5 shows scGPT at 200/500/1000 cells but with coarse resolution. Since the paper's main positive contribution involves scGPT, the scaling behavior of scGPT specifically is important and under-characterized.

### 5. **Cross-Tissue Table 4 (scGPT-12L): The Variance Is Suspiciously Low** ⚠️ LOW-MEDIUM

All 36 cells in Table 4 (12 layers × 3 tissues, TRRUST) fall within 0.711–0.721. This 0.010 total range across all layers and all tissues is remarkably tight—tighter than typical AUROC measurement noise with 27 edges. This could indicate:
- (a) The 27-edge set is so specific that all methods score them similarly (ceiling effect)
- (b) The balanced negative sampling (10:1 ratio, seed-controlled) creates consistent evaluation conditions
- (c) The attention matrices genuinely encode these 27 edges with similar strength everywhere

The DoRothEA results (123 edges, similarly tight: 0.839–0.846) somewhat mitigate this by replicating the pattern with a larger edge set. But the tightness is still notable—it suggests that what's being measured may be properties of these specific TF–target pairs rather than general regulatory signal quality. Bootstrap CIs on these values would help clarify whether the apparent uniformity is within noise.

---

## Minor Issues

### 6. **The Custom Synthetic Generator Justification Is Too Defensive**
The 200-word justification for not using SERGIO (Section 2.11) reads as anticipating criticism rather than making a clear methodological choice. The reasons are valid—SERGIO doesn't model attention matrices—but the defensive tone weakens the presentation. A single sentence would suffice: "We use a custom generator because testing attention-matrix-level predictions requires explicit control over the attention–regulation mapping, which SERGIO does not provide."

### 7. **Abstract Length**
Still ~550 words. For a 25-page paper, this is reasonable, but a 350-word version focusing on the three core findings (co-expression≠regulation, layer-stratified recovery, cross-tissue replication) would be stronger.

### 8. **The "12 Complementary Analyses" Framing**
Listing 12 analyses in the introduction with numbered descriptions reads like a laundry list. Several of these (pseudotime, batch leakage, calibration) are standard quality-control checks that could be grouped. The paper would benefit from a hierarchy: 3 core contributions + 9 supporting analyses.

### 9. **Biological Cell-Type Annotations vs. Unsupervised Clustering**
CSSI is evaluated with oracle cell-type labels (supervised) and briefly mentions "Leiden clustering provides a practical alternative." But no unsupervised CSSI results are shown on real data. Practitioners without ground-truth cell-type labels—the realistic use case—have no evidence that CSSI works for them.

---

## What Has Improved Since R4

1. **Depth-dependent claims properly scoped** — No more overreach about "general architectural principle"
2. **27-edge characterization** — Named TFs, diversity acknowledged
3. **DoRothEA cross-tissue results** — 123 edges with ~0.84 AUROC is strong, tissue-invariant evidence
4. **Consistent model labeling** — scGPT-18L / scGPT-12L throughout
5. **Paper length reduced** — ~25 pages, still long but improved

---

## Path to 8.5: Prioritized Fixes

| Rank | Fix | Effort | Impact |
|------|-----|--------|--------|
| 1 | Resolve stratified baseline edge-set mismatch (evaluate on same edges or withdraw claim) | Low | High |
| 2 | Restructure around 3 core contributions; condense supporting analyses | Medium | High |
| 3 | Add bootstrap CIs to Table 4 AUROC values | Low | Medium |
| 4 | Show unsupervised CSSI (Leiden) on real data | Medium | Medium |
| 5 | Trim abstract to 350 words | Low | Low |

---

## Rating Justification: 7.8/10

**Up from 7.5:** The DoRothEA cross-tissue results (123 edges, ~0.84 AUROC, tissue-invariant) are a meaningful strengthening of the paper's positive claims. The 27-edge characterization and proper scoping of depth-dependent claims resolve two of the highest-priority issues from R4. The paper is now scientifically honest throughout—no claim exceeds its evidence.

**Still below 8.0 because:** (1) The stratified baseline edge-set mismatch—the oldest unresolved issue—still prevents a clean claim that attention outperforms stratified correlation. (2) The paper's diffuse structure (12 analyses of roughly equal weight) prevents any single contribution from landing with sufficient depth. (3) CSSI's real-data value is honestly described as diagnostic, which is useful but modest for a central contribution. (4) The suspiciously low variance in Table 4 needs bootstrap CIs to confirm it's not a ceiling/floor artifact of the small evaluation set.

**For 8.0:** Fix #1 (stratified baseline) is essential—it's been three rounds. Either show the comparison on matched edges or explicitly state the claim is unsupported.

**For 8.5+:** Restructure the paper around its strongest contributions. The co-expression≠regulation insight, layer-stratified recovery, and cross-tissue replication are each publishable findings. The remaining 9 analyses are valuable but belong in a supporting role.

**Bottom line:** The paper has progressed from "overclaimed" (5.5) → "honestly uncertain" (6.5) → "one experiment away" (7.0) → "multi-tissue validated" (7.5) → "properly scoped with strong cross-tissue evidence" (7.8). The trajectory is consistently positive. The gap to 8.0 is a single unresolved methodological issue (stratified baseline). The gap to 8.5+ is structural—the paper needs to decide what it is and commit to that identity.
