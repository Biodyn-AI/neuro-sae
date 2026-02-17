# REVIEW ROUND 2: Second Adversarial Review 
## "A Comprehensive Framework for Mechanistic Interpretability of Single-Cell Foundation Models"

---

## EXECUTIVE SUMMARY

This second review examines the paper after the authors addressed issues from the first round. While they have made substantial improvements—particularly implementing framework-level FDR correction and better integrating their narrative—**significant new issues emerge that were missed in the first review, and some fixes are incomplete or introduce new problems**.

**VERDICT: MAJOR REVISION REQUIRED**

The paper is closer to publication quality but requires addressing several critical issues before acceptance.

---

## VERIFICATION OF FIRST REVIEW FIXES

### ✅ **PROPERLY ADDRESSED:**

**1. Multiple Testing Correction**: The authors now implement proper Benjamini-Hochberg FDR correction across all 47 statistical tests framework-wide (Section 2.3). They honestly report which findings survive correction and which don't. This is exemplary statistical practice.

**2. Narrative Integration**: The contradictory claims about attention failure vs. CSSI success are now well-integrated. The authors clearly explain that attention contains regulatory signal but requires cell-state stratification to extract it properly.

**3. Baseline Comparison Framing**: Section 5.2 now prominently presents the finding that all methods fail equally on brain tissue and correctly interprets this as tissue-specific rather than method-specific limitations.

### ⚠️ **PARTIALLY ADDRESSED:**

**4. CSSI Novelty Claims**: The authors provide much more theoretical justification (Section 2.11), but the core remains basic cell-type stratification. The "theoretical grounding" is mostly mathematical formalization of an obvious biological insight.

**5. Synthetic Validation**: They provide detailed justification for avoiding SERGIO but their custom generator is still specifically designed to validate their hypotheses about attention matrices.

### ❌ **UNRESOLVED:**

**6. Sample Size Issues**: Still only 56 TF-target pairs for pseudotime, 3 cell counts for scaling analysis, etc. The fundamental statistical power problems remain.

---

## NEW CRITICAL ISSUES (Not Identified in First Review)

### **1. SEVERE METHODOLOGICAL INCONSISTENCY: Mixed Edge Score Types**

**Critical flaw**: The paper claims to evaluate "attention-based GRN inference" but several key analyses use **correlation-based edge scores**, not attention-derived ones:

- **Cross-species transfer (Section 6.7)**: Uses Spearman correlations, not attention
- **Pseudotime analysis (Section 6.8)**: Uses expression correlations, not attention  
- **Batch leakage (Section 6.9)**: Uses Pearson correlations, not attention

This creates a **fundamental interpretability problem**: Are the limitations due to attention methods specifically, or edge-based GRN inference generally? The paper conflates these throughout.

**Evidence**: In Section 6.7: "edge scores were computed as Spearman correlation-based edge scores computed independently in each species" - this is NOT testing attention-based inference.

**Impact**: The reader cannot determine which findings apply to foundation model interpretability vs. general GRN inference limitations.

### **2. STATISTICALLY QUESTIONABLE FDR CORRECTIONS**

While the authors implemented framework-level FDR correction (good), there are **serious problems with the correction strategy**:

**Problem A: Selective Reporting Post-Correction**
- Original pseudotime finding: p=0.068, becomes non-significant after correction (adjusted p=0.124)
- Authors still discuss this as a "trend" and include it in their conclusions
- This is statistically inappropriate - corrected non-significant results should be reported as null findings

**Problem B: Effect Size Interpretation Without Statistical Significance**
- They report Cohen's d=1.58 for pseudotime directionality despite p>0.05 after correction
- Large effect sizes without statistical significance in small samples are often spurious
- 56 pairs is grossly underpowered for reliable effect size estimation

**Problem C: Multiple Correction Levels**
- They apply both within-analysis FDR (e.g., "No individual pair reached significance after FDR correction" in pseudotime)
- AND framework-level FDR across all analyses
- This is **double-penalization** - statistically inappropriate

### **3. CIRCULAR REASONING IN CSSI VALIDATION**

**Critical flaw in real attention validation**: The authors claim CSSI works on "real scGPT attention matrices from 497 human brain cells" achieving "AUROC 0.694-0.706" (Section 6.10). However:

**The validation is circular**: 
1. They identify layers 13-14 as having the highest regulatory signal
2. They apply CSSI to these pre-selected layers  
3. They conclude CSSI works because it performs well on these layers

**Missing control**: What happens when you apply CSSI to the layers that DON'T work (layers 0-6 with AUROC ~0.51-0.61)? If CSSI is truly solving the heterogeneity problem, it should improve even bad layers.

**Selection bias**: The "best individual heads achieving AUROC 0.706" is likely p-hacking across 18 layers × multiple heads. No correction for this multiple testing is applied.

### **4. UNADDRESSED CAUSAL INFERENCE GAPS**

**Fundamental conceptual flaw**: The paper assumes that recovering edges from reference databases (TRRUST, DoRothEA) validates causal regulatory inference. This is wrong for several reasons:

**Problem A: Reference databases contain correlational relationships**
- Many TRRUST edges are discovered via co-expression studies, not causal perturbation
- Using co-expression to validate co-expression is circular

**Problem B: Missing temporal validation**
- Even when they attempt temporal validation (pseudotime), it fails for 79% of edges
- They dismiss this as a "limitation of pseudotime" rather than their methods

**Problem C: Confounding experimental validation**
- Their perturbation validation (Section 6.5) shows weak/inconsistent results
- Only 1 of 4 datasets shows significant consistency, and it disappears after confound adjustment in most cases

### **5. OVERLOOKED COMPUTATIONAL VALIDATION PROBLEMS**

**Attention weight extraction issues not addressed**:

**Problem A: Layer-dependent findings are not validated across models**
- scGPT shows regulatory signal in layers 13-14
- Geneformer analysis doesn't examine layer-specific patterns  
- Cross-model generalization of layer effects is untested

**Problem B: Head-level analysis lacks statistical rigor**
- Authors report individual head AUROC up to 0.706 but provide no confidence intervals
- No correction for testing multiple heads per layer
- Could be random fluctuation rather than genuine regulatory signal

**Problem C: Missing negative controls**
- No analysis of attention weights on random/permuted data
- No validation that high-performing heads are consistently the same across datasets
- No control for whether observed patterns could arise from random attention matrices

### **6. BIOLOGICAL IMPLAUSIBILITY OF KEY FINDINGS**

**Cross-species conservation claims are problematic**:

**Problem A: Species comparison confounded by dataset differences**
- Human lung: 65,847 cells from Tabula Sapiens (10X)
- Mouse lung: 9,409 cells from Krasnow (Smart-seq2)  
- Different protocols, cell compositions, and sample sizes confound species effects

**Problem B: Conservation interpretation is backward**
- Authors claim "lineage-specifying TFs transfer well; signaling-responsive TFs do not"
- But this pattern could reflect dataset artifacts (batch effects, protocol differences) rather than genuine biology

**Problem C: Missing proper controls**
- No analysis of within-species, cross-protocol conservation as control
- No validation that "conserved" edges replicate in independent same-species datasets

### **7. INADEQUATE TREATMENT OF NEGATIVE RESULTS**

**Problem**: The paper presents mostly negative results (attention fails, scaling is inverse, validation methods don't work) as a "comprehensive framework" but lacks proper statistical power analysis for null findings.

**Specific issues**:
- Geneformer AUROC confidence intervals include 0.5 (random), but authors don't formally test non-inferiority to random
- "Near-random performance" is presented as definitive failure without power analysis
- Multiple null results without considering that true effect sizes might be smaller than detectable with current sample sizes

---

## BIOLOGICAL/DOMAIN-SPECIFIC CONCERNS

### **1. QUESTIONABLE GROUND TRUTH ASSUMPTIONS**

The authors acknowledge that TRRUST/DoRothEA are incomplete but still use failure against them as definitive evidence. This creates several problems:

**Problem A**: True regulatory networks are likely much sparser than current methods can detect
**Problem B**: Context-specificity means tissue-specific regulatory programs may not appear in general databases
**Problem C**: Temporal dynamics of regulation are poorly captured in static reference databases

### **2. CELL-TYPE COMPOSITION CONFOUNDS NOT FULLY ADDRESSED**

While the authors identify batch/donor leakage, they underexamine how cell-type composition differences drive their key findings:

- Cross-tissue differences may reflect cell-type composition, not regulatory differences
- Cross-species differences are partially driven by differential cell proportions
- Scaling effects could be due to changing cell-type ratios with sample size

---

## STATISTICAL/METHODOLOGICAL CONCERNS

### **1. BOOTSTRAP CONFIDENCE INTERVALS WITHOUT FINITE SAMPLE CORRECTION**

The authors use 10,000 bootstrap resamples throughout but don't address finite sample bias:
- Bootstrap CIs can be biased for small samples (n<100)
- Multiple analyses have small sample sizes where bootstrap assumptions may be violated
- No bias-corrected and accelerated (BCa) bootstrap intervals provided

### **2. MISSING EFFECT SIZE REPORTING AND INTERPRETATION**

While some analyses report effect sizes (Cohen's d), many don't:
- Cross-tissue correlations lack standardized effect sizes
- CSSI improvements reported as ratios (1.85×) without confidence intervals
- Clinical significance vs. statistical significance not distinguished

---

## PRESENTATION AND CLARITY ISSUES

### **1. MISLEADING FIGURE PRESENTATIONS**

**Figure quality issues**:
- Many confidence intervals are mentioned in text but not shown in figures
- Statistical significance markers (*, **) are inconsistently applied
- Error bars are sometimes standard error, sometimes confidence intervals (inconsistent)

### **2. INCONSISTENT TERMINOLOGY**

- "Edge scores" sometimes means attention weights, sometimes correlations
- "GRN inference" conflates causal discovery with edge ranking
- "Validation" includes both cross-reference comparison and experimental perturbation

---

## RECOMMENDATIONS FOR REVISION

### **ESSENTIAL CHANGES (Required for Acceptance):**

1. **Methodological consistency**: Clearly separate attention-based analyses from correlation-based analyses. The paper should either focus on attention methods OR general edge-based GRN inference, not conflate them.

2. **Fix FDR correction strategy**: Choose either within-analysis OR framework-level correction, not both. Report corrected non-significant results as null findings.

3. **Add proper negative controls**: Test CSSI on known-bad layers, analyze random attention matrices, validate head consistency across datasets.

4. **Expand sample sizes**: At minimum, double the pseudotime analysis pairs and add intermediate scaling points.

5. **Cross-model layer analysis**: Validate that layer-dependent regulatory signals generalize to Geneformer and other models.

### **STRONGLY RECOMMENDED CHANGES:**

6. **Power analysis for null findings**: Formally test whether "near-random" performance is statistically indistinguishable from random chance.

7. **Independent external validation**: Test key findings (CSSI, layer effects, scaling failure) on completely held-out datasets.

8. **Biological plausibility analysis**: Add proper controls for cross-species comparison, validate within-species consistency.

### **MINOR IMPROVEMENTS:**

9. **Consistent figure presentation**: Show confidence intervals in figures, not just text
10. **Effect size reporting**: Add standardized effect sizes throughout
11. **Terminology consistency**: Define and consistently use technical terms

---

## STRENGTHS OF THE REVISION

The authors deserve credit for:
- Implementing proper framework-level statistical correction
- Providing detailed methodological justifications  
- Being more honest about limitations and negative results
- Better integrating their complex narrative
- Adding baseline comparisons that strengthen rather than weaken their arguments

---

## CONCLUSION

This paper addresses an important problem and makes genuine contributions to foundation model interpretability. The first revision addressed several critical issues successfully. However, **new methodological problems have emerged that are equally serious**:

1. **Methodological inconsistency** (mixing attention and correlation analyses)
2. **Circular CSSI validation** (pre-selecting good layers)
3. **Inappropriate FDR correction strategy** (double-penalization)  
4. **Inadequate negative controls** (missing power analysis)

The paper is much improved but requires another round of careful revision to achieve publication quality at a top venue.

**RECOMMENDATION: MAJOR REVISION REQUIRED**

---

## LINE-SPECIFIC ISSUES

**Section 2.11 (CSSI Methods)**: The claim of "rigorous theoretical grounding" oversells basic cell-type stratification.

**Section 6.7**: Should clarify early that this uses correlation, not attention-based edges.

**Section 6.10**: The layer pre-selection creates circular validation that undermines the CSSI conclusions.

**Table 6**: AUROC values should include confidence intervals, not just be mentioned in text.

**Discussion**: Should acknowledge that most findings are about general GRN inference limitations, not attention-specific issues.