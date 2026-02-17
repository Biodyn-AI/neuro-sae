# Unified Nature Machine Intelligence Paper

## Paper Selection and Rationale

This unified paper combines the strongest results from 5 key Biodyn mechanistic interpretability papers to create a comprehensive framework for interpreting single-cell foundation models. The selection was based on:

1. **Novel, counterintuitive findings** that challenge field assumptions
2. **Strong theoretical foundations** with practical applications  
3. **Coherent narrative** covering the full scope of mechanistic interpretability
4. **Complementary methodological approaches** that strengthen each other
5. **High impact potential** for Nature Machine Intelligence readership

## Selected Papers

### 1. Subproject 06 - Scaling Failure Analysis
**Original target**: Neural Networks  
**Key contribution**: Demonstrates that increasing cell counts from 200→1000→3000 can systematically degrade rather than improve GRN recovery across all model tiers (small/medium/large). This challenges the fundamental assumption that "more data = better interpretability."

**Why included**: This is a striking, counterintuitive finding that will grab NMI reviewers' attention and has major implications for how the field approaches data collection and analysis.

### 2. Subproject 22 - Bias Bounds for Patching-Based Mediation  
**Original target**: PLOS Computational Biology  
**Key contribution**: Provides theoretical framework for quantifying bias in single-component mediation analysis when interaction effects are non-negligible. Shows that standard ranking approaches can be systematically misleading.

**Why included**: Strong theoretical foundation that addresses a fundamental methodological issue affecting most mechanistic interpretability studies. Provides practical tools (bias bounds, ranking certificates) that researchers can immediately apply.

### 3. Subproject 21 - Detectability Phase Diagram
**Original target**: Annals of Applied Statistics  
**Key contribution**: Establishes closed-form statistical foundations for determining when mechanistic signals are recoverable under realistic experimental constraints. Shows regime-dependent performance of attention vs intervention signals.

**Why included**: Provides the statistical rigor that NMI expects. Moves the field from ad-hoc evaluation to principled statistical analysis. Complements the bias analysis by addressing when signals are detectable at all.

### 4. Subproject 20 - Invariant Causal Edges Across Tissues
**Original target**: PLOS Computational Biology  
**Key contribution**: Reveals substantial heterogeneity in mechanistic insights across tissue types (Spearman -0.44 to 0.71), challenging assumptions about transferability of regulatory relationships.

**Why included**: Addresses biological validity and generalization—key concerns for NMI. Shows that mechanistic claims may not transfer as broadly as assumed, requiring context-specific validation.

### 5. Subproject 11 - Counterfactual Perturbation Consistency
**Original target**: Genome Biology  
**Key contribution**: Develops framework for validating mechanistic interpretations against perturbation experiments. Shows improved consistency when accounting for interaction effects.

**Why included**: Provides validation methodology that connects computational predictions to experimental reality. Demonstrates that the framework's interaction-aware approaches yield more reliable biological predictions.

## Unified Narrative

These papers combine into a coherent story:

1. **Foundation**: Statistical frameworks for when mechanistic signals are detectable and how to quantify bias (Papers 2, 3)
2. **Core finding**: Scaling behavior is counterintuitive—more data can hurt rather than help (Paper 1)  
3. **Biological validity**: Mechanistic insights don't transfer reliably across contexts (Paper 4)
4. **Validation**: Interaction-aware approaches improve consistency with experimental validation (Paper 5)

## Why This Combination for NMI

**Scope**: Covers methodology, theory, empirical findings, and biological validation—exactly what NMI wants to see.

**Impact**: Challenges fundamental assumptions in the field while providing practical solutions.

**Novelty**: Multiple counterintuitive findings supported by rigorous analysis.

**Utility**: Provides immediately applicable tools and guidelines for the community.

**Biological insight**: Demonstrates that proper statistical analysis reveals important limitations in how we interpret foundation models for biology.

## Papers NOT Included

Several strong papers were not included to maintain focus and coherence:

- **Subproject 25 - Evaluation Invariance Metrics**: Excellent methodological work but overlaps with detectability analysis
- **Subproject 28 - Cell State Conditional Circuits**: Interesting but more specialized application
- **Subproject 30+ papers**: Many focus on specific applications rather than fundamental methodology

## Target Impact

This unified paper should:
1. **Change practice** in how researchers interpret foundation models
2. **Establish standards** for statistical rigor in mechanistic interpretability  
3. **Provide tools** that researchers can immediately apply
4. **Challenge assumptions** that currently guide the field
5. **Bridge** computational methodology and biological insight

The combination creates a comprehensive framework that is greater than the sum of its parts while maintaining the rigorous standards expected by Nature Machine Intelligence.