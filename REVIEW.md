# Adversarial Review: "A Comprehensive Framework for Mechanistic Interpretability of Single-Cell Foundation Models"

## Overall Recommendation: **REJECT**

This paper attempts to provide a "comprehensive framework" for mechanistic interpretability of single-cell foundation models but suffers from fundamental flaws in experimental design, statistical rigor, and novelty claims that render it unsuitable for publication at a top-tier venue.

---

## MAJOR ISSUES (Rejection-Worthy)

### 1. **Catastrophic Multiple Testing Problem**
The paper conducts **dozens of statistical tests** across 12 different analyses without any framework-level multiple testing correction. The authors acknowledge this explicitly: "we do not apply a framework-level correction across all analyses" claiming the analyses are "distinct." This is statistically irresponsible. Key findings like pseudotime directionality (p=0.068) and perturbation validation would likely disappear under proper Bonferroni or FDR correction across the entire study.

**Specific example**: The cross-tissue consistency analysis reports "only two of six pair-granularity comparisons surviving FDR control at α = 0.05" - meaning 67% of their tests failed even within-analysis correction. With framework-level correction, likely zero findings would survive.

### 2. **Deeply Flawed Synthetic Validation Design**
The authors acknowledge creating a "custom generator" instead of using the established SERGIO simulator, claiming they need "explicit control over the mapping from ground-truth regulatory weights to synthetic attention matrices." This is circular reasoning - they've engineered their synthetic data to confirm their hypotheses about attention weights. SERGIO is the community standard for exactly this type of validation, and avoiding it appears designed to circumvent rigorous testing.

### 3. **Cherry-Picked Baseline Comparisons**
Section 5.2 reveals the authors tested attention methods on brain tissue where **ALL methods fail equally** (AUROC 0.521-0.526 for Spearman, MI, GENIE3, GRNBoost2, and attention). This undermines their entire narrative - if state-of-the-art GRN methods achieve identical performance to attention, then the "failure" is tissue-specific, not method-specific. The authors bury this crucial finding and continue to criticize attention methods throughout the paper.

**Critical flaw**: Why not test on tissues where GRN methods actually work? This appears designed to make attention methods look bad by choosing the worst possible evaluation context.

### 4. **Trivial Methodological Contribution (CSSI)**
Cell-State Stratified Interpretability (CSSI) is simply **cell-type stratification before analysis** - an obvious preprocessing step that any competent practitioner would consider. The authors present this as a novel "constructive solution," but it's basic good practice. The improvement factors (1.85x) are meaningless when baseline performance is near-random (AUROC ~0.5).

### 5. **Severe Sample Size Issues**
- Pseudotime analysis: Only **56 TF-target pairs** across 3 lineages. This is grossly underpowered.
- Cross-species analysis: Limited to lung tissue from 2 species with different cell compositions, confounding species effects with dataset effects.
- Scaling analysis: Only **3 cell count points** (200, 1000, 3000) - insufficient to characterize scaling behavior.

### 6. **Contradictory Core Narrative**
The paper simultaneously claims:
- Attention-based GRN inference "fails across architectures" (Section 6.12)  
- CSSI improves GRN recovery using attention matrices (Section 6.10)
- Real scGPT attention achieves "AUROC 0.683-0.694" which "substantially exceeds" baselines

These claims are logically inconsistent. If attention truly fails fundamentally, CSSI shouldn't work.

### 7. **Misleading Statistical Presentation**
The authors repeatedly present AUROC values around 0.5 as meaningful when they're statistically indistinguishable from random chance. Bootstrap confidence intervals frequently include 0.5, but this is buried in dense text. The entire scaling failure analysis could be noise fluctuation around random performance.

### 8. **Unfalsifiable Ground Truth Problem**
The authors acknowledge that TRRUST/DoRothEA reference databases are incomplete and context-specific, then use failure against these databases to conclude attention methods don't work. This creates an unfalsifiable situation - any method failing against incomplete ground truth can claim the ground truth is wrong rather than the method failing.

---

## MINOR ISSUES

### Writing and Presentation
- **Excessive density**: 725 lines, 12 analyses - impossible to digest meaningfully
- **Inconsistent methodology**: Different tissues, cell counts, and evaluation metrics across sections
- **Overstated claims**: "Comprehensive framework" for what amounts to mostly negative results
- **Poor figure quality**: Many plots are difficult to interpret with crowded legends and small text

### Technical Issues
- Geneformer evaluation limited to brain tissue where all methods fail
- Batch leakage analysis uses correlation-based edges, not attention-derived edges
- "Frozen archives" mentioned for mediation analysis - not reproducible
- Calibration analysis uses 45.8% positive rate that may not reflect true regulatory sparsity

### Missing Related Work
- Insufficient coverage of causal inference methodology
- Limited connection to mechanistic interpretability beyond attention weights  
- Missing recent work on GRN inference benchmarking
- No discussion of why attention should encode regulation in the first place

---

## SPECIFIC LINE-LEVEL SUGGESTIONS

**Line 45-47**: The claim about "millions of cells" is vague - specify exact training corpus sizes.

**Line 234**: "Systematic scaling analysis" with only 3 points is not systematic.

**Line 445**: The non-additivity finding (62.5% rate) lacks confidence intervals.

**Line 567**: Cross-tissue Spearman correlations ranging from -0.44 to 0.71 suggests these measurements are noise.

**Line 892**: The pseudotime "failure rate" of 79% is misleading - these are well-established regulatory relationships, so the failure is in the validation method, not the biology.

**Tables 4-5**: AUROC confidence intervals should be prominently displayed, not buried in text.

**Figure 1**: The scaling failure plot needs error bars and significance testing between adjacent points.

---

## SUGGESTED EXPERIMENTS TO STRENGTHEN THE PAPER

1. **Proper positive controls**: Test on tissues/contexts where established GRN methods actually work, not just brain where everything fails.

2. **SERGIO validation**: Re-run synthetic experiments using the community-standard SERGIO simulator instead of custom synthetic data.

3. **Framework-level statistical correction**: Apply Bonferroni or FDR correction across all statistical tests in the framework.

4. **Dense scaling analysis**: Test 8-10 cell count points to properly characterize scaling behavior, not just 3 points.

5. **Cross-architecture CSSI**: Test CSSI on Geneformer and other architectures, not just scGPT.

6. **Perturbation-aware baselines**: Compare against methods that actually use perturbation data (CellOracle, etc.) rather than just correlation methods.

7. **Independent external validation**: Test on completely held-out datasets not used in method development.

---

## VERDICT

While the authors tackle an important problem, the execution is fundamentally flawed. The synthetic validation appears engineered to support their conclusions, the statistical analysis ignores multiple testing, and the core methodological contribution (CSSI) is trivial. The paper reads more like an extended negative result disguised as a "comprehensive framework."

Most critically, the cherry-picked baseline comparison (Section 5.2) reveals that the entire premise may be wrong - if dedicated GRN methods fail equally on the test data, then criticizing attention methods is unfair. This should have been the starting point for a more nuanced analysis, not buried as a minor finding.

The paper would require complete redesign with proper statistical controls, established synthetic benchmarks, and testing on datasets where methods actually work to be suitable for publication at a top venue.

**Recommendation: Reject**