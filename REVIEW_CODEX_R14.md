# Codex Review R14 (adversarial)

**Paper:** *Attention in Single-Cell Foundation Models Reflects Co-Expression, Not Regulation: A Diagnostic Framework with Null-Model Benchmarking*  
**File reviewed:** `D:\openclaw\biodyn-nmi-paper\main.tex`  
**Round:** 4 (prior ratings: 5.5 → 6.5 → 7.0)

## Overall score (1–10)
**7.6 / 10**

This is a meaningful improvement over the previous version primarily because the paper now **actually runs** the expression-/marginal-matched negative evaluation and uses it to resolve the central ambiguity. However, the new result is still **statistically and methodologically under-specified** (small positive set, unclear uncertainty, matching definition not fully transparent), so it does not yet fully “close” the interpretability claim.

## Has the expression-matched evaluation addressed the #1 blocker?
**Mostly yes (≈80% addressed).**

The previous #1 blocker was: *unmatched AUROC could be entirely due to expression/marginal confounds; without a matched-negative evaluation, you can’t claim pairwise signal.*

The updated manuscript now includes an explicit matched-negative experiment:
- **Matched negatives:** “50 per positive, ±20% tolerance” (matched on **marginal attention** for TF and target)
- **Result:** attention AUROC **0.623** vs marginal-product baseline **0.504** (chance)
- **Interpretation:** provides **direct evidence of residual pairwise structure** beyond marginal effects.

That is the right diagnostic and it materially changes credibility: the paper no longer rests on purely conceptual confound arguments.

## Major strengths (what is now working)
1. **Central narrative coherence improved.** The paper now cleanly reconciles: (i) co-expression dominance and (ii) above-chance AUROC on curated benchmarks, via gene-level null models *plus* the matched-negative test that isolates pairwise signal.
2. **Correct posture about magnitude.** You explicitly call the remaining signal “modest” and avoid overclaiming.
3. **Null-model benchmarking is now sharper.** Showing simple gene-level products ≥ attention on unmatched AUROC sets up the need for matched negatives very well.

## Remaining major issues (still adversarial)
### 1) Matching criterion is not fully convincing / could be circular
You describe this as **“expression-matched negative sampling”**, but the operational definition is **matching on marginal attention** (within ±20%).
- If marginal attention is itself already a function of expression prevalence/detection/etc., that may be fine *as a control for the specific confound in your pipeline*, but it is not the same as matching on expression statistics.
- Calling it “expression-matched” risks pushback. Reviewers may ask: *why not match directly on detection rate / mean expression / variance (or a multivariate propensity score) rather than on attention marginals?*

**Actionable fix:** rename to **“marginal-attention-matched negatives”** (or do both: match on expression metrics and show similar residual AUROC).

### 2) Very small positive set (n = 27)
The matched-negative claim hinges on only **27 TRRUST positives** in the scGPT-12L cross-tissue setup.
- With n=27, AUROC variance is large; a 0.623 estimate could swing materially with small perturbations.
- You should report **confidence intervals** and/or a **paired test** (e.g., DeLong CI for AUROC; permutation over labels within the matched sets) for the 0.623 vs 0.504 gap.

**Actionable fix:** add CI/p-values and ideally repeat with a larger positive set (e.g., DoRothEA, expanded TRRUST mapping, or relax constraints to include more TFs/targets).

### 3) Baseline choice in the matched setting is too narrow
You compare attention to a “marginal-product baseline” that is constructed to be chance under the matching. That’s expected.

**Actionable fix:** include at least one baseline that is *not* trivially neutralized by the matching procedure, e.g.
- a simple pairwise score from expression (correlation/MI) computed **within the same matched candidate pool**,
- or a conditional model that uses gene-level covariates but can still rank within the matched group.

This would clarify whether the residual 0.623 is “attention-specific” or just “any weak pairwise statistic survives.”

### 4) Insufficient methodological detail for reproduction
The matched-negative procedure is described in prose but not at the level where another lab can reproduce without code.

**Actionable fix:** specify explicitly:
- how “marginal attention” is computed (row/column mean across cells? across heads? normalization?),
- whether matching is per-head or pooled,
- whether negatives are sampled with/without replacement,
- how AUROC is aggregated across the 50 negatives per positive (flattened vs per-positive ranking then averaged),
- random seed / repeats.

## Minor comments
- Cross-tissue TRRUST evaluation with 27–28 edges is framed appropriately as limited, but it is still surprisingly small; consider moving it from “evidence” tone toward “sanity check” tone.
- The paper would benefit from a compact “Decision tree” figure: **(1) unmatched AUROC + null models → (2) matched-negative test → (3) only then interpret pairwise signal**.

## Bottom line
The #1 blocker (no matched-negative evaluation) is **no longer blocking** because the paper now contains a concrete, on-the-record experiment showing **above-chance pairwise signal after controlling marginals**.

But the current implementation is **too small and under-specified** to carry the full argumentative weight you put on it. With (i) clearer naming, (ii) uncertainty statistics, and (iii) a larger/replicated matched-positive set, this could jump into the **8+** range.
