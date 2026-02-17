# Review Iteration 7 — Adversarial Review

**Paper:** A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models  
**Rating: 8.4 / 10** (previous: 8.1)  
**Trajectory:** 5.5 → 6.5 → 7.0 → 7.5 → 7.8 → 8.1 → **8.4**

---

## Summary of Changes

Full structural restructuring around 3 core contributions, with supporting analyses explicitly deemphasized. Paper trimmed and reorganized.

---

## Does the 3-Contribution Structure Work?

**Yes — this is the single biggest improvement the paper has made.** The three pillars are:

1. Attention encodes co-expression, not regulation (the diagnostic)
2. Layer-stratified CSSI recovers signal (the method)
3. Cross-tissue replication confirms it's real (the validation)

This creates a clean narrative arc: problem → solution → proof. The reader no longer drowns in nine parallel analyses of equal apparent weight. The "Supporting Analyses" subsection header with the explicit framing paragraph ("The following analyses provide quality control...") is well executed.

**Minor structural issue:** The supporting analyses still consume ~60% of the Results section by word count. The deemphasis is structural (clear labeling, explicit framing) but not proportional. A reader scanning the paper could still feel overwhelmed. Consider whether some supporting analyses could move to supplementary material in the final submission.

---

## Abstract Assessment

**Crisp and effective.** The abstract now opens with "organized around three core findings" and delivers each with a bold-face label. Key numbers are present. The final sentence about supporting analyses is appropriately brief.

**One issue:** The abstract is ~200 words, which is fine for NMI, but the sentence about supporting analyses is a laundry list ("scaling behavior, mediation bias, detectability theory, perturbation validation, cross-species transfer, pseudotime directionality, batch leakage, and uncertainty calibration"). This reads like a table of contents rather than a summary. Consider replacing with something like "...supported by comprehensive quality-control analyses spanning scaling behavior, cross-species transfer, and technical artifact assessment."

---

## What Improved

1. **Narrative clarity.** The introduction now explicitly states the three contributions with section references. The discussion mirrors this structure. A reader can extract the paper's story in 60 seconds.

2. **Honest scoping of CSSI.** The paper now correctly positions CSSI as primarily a diagnostic tool ("CSSI's primary value on real data is diagnostic") rather than overselling it as a universal performance booster. The circularity caveat for held-out validation is explicitly flagged. This is mature scientific writing.

3. **Cross-tissue replication as Contribution 3.** Elevating this from a supporting analysis to a core finding was the right call. The scGPT-12L results (AUROC ~0.72 TRRUST, ~0.84 DoRothEA across brain/kidney/whole-human) are the paper's strongest evidence that something real is being recovered.

4. **Architecture-specificity acknowledged.** The paper now clearly states that the depth-dependent pattern (L13-14 in scGPT-18L) does NOT replicate in scGPT-12L, and that this is informative rather than problematic. Good science.

5. **Framework-level multiple testing correction.** The 47-test BH-FDR correction is a significant methodological strength. The honest reporting that pseudotime doesn't survive correction (p=0.124) builds credibility.

---

## Remaining Issues

### Major

1. **The scGPT-12L TRRUST evaluation set is tiny (27-28 edges).** This is acknowledged in the table caption but under-discussed. With only 27-28 edges and 14 TFs, the cross-tissue "consistency" could partly reflect the dominance of a few well-known TFs (STAT3, TP53, NFKB1, SP1, MYC) that are constitutively expressed across tissues. The AUROC of ~0.72 on 28 edges has wide confidence intervals (±0.03-0.05 as noted), making cross-tissue differences of <0.002 essentially noise. **The paper should more prominently discuss whether this replication is meaningful or trivially expected given the evaluation set composition.** This is the paper's most important validation — it needs to be bulletproof.

2. **Tension between Core Finding 1 and Core Finding 3.** Finding 1 says attention encodes co-expression, not regulation (ρ≈0 with TRRUST). Finding 3 says attention recovers regulatory signal (AUROC ~0.72 TRRUST). These use different models (scGPT-18L pooled vs. scGPT-12L per-layer), different gene sets, and different evaluation edge counts (8,330 vs. 27-28). The paper explains this (layer selection, architecture differences), but the surface-level contradiction deserves a more explicit reconciliation paragraph, ideally in the Discussion. A skeptical reviewer will flag this.

3. **Title remains too long.** "A Systematic Framework for Mechanistic Interpretability of Single-Cell Foundation Models: Attention-Based Gene Regulatory Network Inference and Cell-State Stratified Analysis" is 23 words. NMI typically favors shorter titles. Suggestion: "Attention in Single-Cell Foundation Models Encodes Co-Expression, Not Regulation: A Diagnostic Framework" or similar.

### Minor

4. **The baseline comparison (Section 3.4) uses brain tissue with 500 cells and 500 genes.** This is a restricted setting. GENIE3/GRNBoost2 were designed for larger gene sets. The claim that "all methods show similar poor performance" could be an artifact of the restricted gene set. Acknowledge this.

5. **Synthetic validation (Section 3.12) is buried.** It's a one-paragraph section after all supporting analyses. Given that it validates three theoretical predictions, it deserves slightly more visibility — perhaps as a subsection under Core Finding 2, since it directly supports CSSI.

6. **The "Competing Effects" biological interpretation of scaling** (below 200 cells = insufficient power, 200-750 = sweet spot, >750 = dilution) is stated as fact but is one of several possible explanations. The alternative explanation paragraph (overfitting) is good, but the "competing effects" framing in the main text should be softened to "consistent with" rather than "reflects."

7. **DoRothEA AUROC (~0.84) vs TRRUST (~0.72) difference** in scGPT-12L is noted but under-explained. The paper attributes it to "larger evaluation set and ChIP-seq-derived nature." But ChIP-seq edges are physical binding, not functional regulation — this actually makes the higher DoRothEA score *less* impressive for regulatory inference, not more. Discuss.

8. **Missing: computational cost section.** The paper mentions attention extraction takes 0.1s vs. 89-127s for GENIE3/GRNBoost2, but there's no systematic discussion of CSSI's computational overhead (stratification + per-stratum scoring + aggregation).

---

## New Issues from Restructuring

- **Cross-tissue consistency (Section 3.7) vs. Cross-tissue replication (Section 3.3).** These are different analyses (mediation-based cross-tissue consistency vs. layer-stratified AUROC replication) but the naming is confusingly similar. The mediation-based Section 3.7 shows *limited* cross-tissue consistency, while Core Finding 3 shows *strong* cross-tissue replication. A reader may see contradiction. Add a brief note distinguishing what's being compared in each.

- **Section numbering flow.** The "Supporting Analyses" header uses `\subsection*` (unnumbered) followed by numbered subsections. This creates a visual break that works but may confuse automated reference systems.

---

## Score Justification

| Criterion | Score | Notes |
|-----------|-------|-------|
| Novelty | 8/10 | Co-expression insight is genuinely useful; CSSI is practical |
| Rigor | 8.5/10 | Framework-level FDR correction, honest null results, caveats |
| Clarity | 8.5/10 | Major improvement from restructuring; 3-pillar narrative works |
| Significance | 8/10 | Important for the growing foundation model interpretability field |
| Completeness | 9/10 | Comprehensive supporting analyses; limitations well-documented |
| Presentation | 8/10 | Title too long; some redundancy in supporting analyses |

**Overall: 8.4/10**

The restructuring was the right move and executed well. The paper now tells a coherent story. The remaining gap to 9+ is: (1) the fragility of the scGPT-12L replication on 27-28 edges, (2) the unexplained tension between Finding 1 and Finding 3, and (3) presentation polish (title, proportional deemphasis of supporting material).

---

## Priority Actions for Iteration 8

1. **Add explicit reconciliation** of Finding 1 (ρ≈0) vs Finding 3 (AUROC 0.72) in the Discussion
2. **Discuss the 27-28 edge evaluation set limitation** more prominently in Section 3.3
3. **Shorten the title**
4. **Trim the abstract's supporting-analysis list** to a summary phrase
5. **Distinguish Section 3.7 from Section 3.3** with a clarifying note
