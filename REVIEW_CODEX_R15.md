# REVIEW_CODEX_R15 (adversarial)

**Paper:** *Attention in Single-Cell Foundation Models Reflects Co-Expression, Not Regulation: A Diagnostic Framework with Null-Model Benchmarking*

## Overall score (1–10)
**7.8 / 10**

This draft is materially stronger than prior rounds: it squarely addresses the biggest prior circularity concern (expression-matched evaluation uses *actual expression statistics*, not attention marginals), adds bootstrap CIs, and now explicitly caveats the small 27-edge evaluation regime. The narrative is also more coherent: (1) pooled attention ≈ co-expression, (2) CSSI localizes benchmark-discriminative signal, (3) cross-tissue invariance + gene-level nulls show benchmark degeneracy, with expression-matching providing the key “pairwise-signal survives” counterpoint.

For *top-venue* acceptance, the remaining blockers are mostly (i) **evidentiary strength** around the key “pairwise signal” claim (n=27 positives), (ii) **internal consistency / numerical mismatches**, and (iii) **scope/positioning**: it reads like a bundle of 12 analyses with uneven maturity rather than a single crisp, validated contribution.

---

## Summary of contributions (as I interpret them)
1. **Main negative result:** pooled attention edges correlate with co-expression (ρ≈0.31–0.42) and not with TRRUST labels (ρ≈0), explaining AUROC≈0.5 pooled GRN recovery.
2. **CSSI diagnostic:** layer- and state-stratification can localize where benchmark-discriminative signal sits (notably in scGPT-18L late layers), with strong gains in synthetic settings.
3. **Null-model benchmarking & invariance:** cross-tissue + cross-layer AUROC invariance (scGPT-12L) is consistent with a benchmark confound; simple gene-level products match/exceed attention on standard negatives.
4. **Key “rescue” test:** expression-matched negatives (mean expression + detection rate matching) reduce gene-level baselines to near-chance while attention retains AUROC ~0.65 with wide CI.

---

## Major strengths
- **The central critique is important and timely.** The paper targets a real failure mode: treating attention-as-regulation in single-cell transformers without strong controls.
- **Null-model framing is a strong conceptual contribution.** Demonstrating that benchmark AUROC is recoverable from trivial gene-level statistics is a valuable diagnostic pattern others can reuse.
- **Expression-matched evaluation is the right move.** Matching on *expression statistics* (not attention-derived marginals) is the correct way to test whether pairwise structure exists beyond gene-level prominence.
- **Honest limitation language has improved.** Multiple sections explicitly call out underpower, circularity risks, and benchmark enrichment effects.

---

## Top-venue blockers (what still prevents an 8.5–9+)

### 1) Internal numerical inconsistency around the key matched-negative result (serious)
The paper currently reports **multiple different numbers** for the same conceptual result:
- In **Core Finding 3**: attention AUROC **0.646** (CI [0.539, 0.747]) vs expression baseline **0.522**.
- In the **Conclusions**: attention retains AUROC **0.623** vs baseline **0.504**.
- In the **Introduction/Core contrib 3 preview**: attention retains AUROC **0.623** vs baseline **0.504**.

This is not a cosmetic issue: the *entire* “pairwise signal exists” claim hangs on this experiment. Top venues will not tolerate key-number drift across abstract/introduction/results/conclusion.

**Fix:** pick one canonical matched-negative protocol + model/layer aggregation; report it consistently everywhere. If there are multiple variants (expression-matched vs attention-matched; different tissues; different layer pooling), present them as separate rows in a table with explicit labels.

### 2) Key positive claim remains underpowered (n=27 positives) and looks fragile
Even with corrected matching, the matched-negative analysis is based on **27 TRRUST edges**. Wide CI is acknowledged, but for a top venue this still reads like: “we found modest signal, maybe.”

What would convince a skeptical reviewer:
- A **larger positive edge set** in the same experimental regime, *without reintroducing circularity*. Options:
  - Use **larger curated resources** with functional/perturbation grounding (even if noisier) and show robustness of the matched-negative gap.
  - Construct a **held-out perturbation-derived edge set** (Perturb-seq) to define positives in a way orthogonal to expression.
  - If limited to TRRUST, expand across tissues / settings and show the matched-negative effect replicates (same matching protocol) in ≥2 contexts.
- A **per-TF stratified analysis**: if only 14 TFs exist, show the effect is not driven by 1–2 TFs. A leave-one-TF-out sensitivity or mixed-effects model would help.

### 3) Benchmark degeneracy argument is plausible but could be made tighter/cleaner
The claim “standard AUROC reflects gene prominence” is believable, and the null-model table is strong. But for a top venue, I’d expect:
- A clearer separation of **(a) benchmark bias** vs **(b) model failure** vs **(c) biology (context-specific regulation)**.
- More explicit discussion of **AUROC pathologies** under heavy class imbalance / edge sampling, and why AUROC is the right summary.
- A stronger justification for the negative sampling strategy in the unmatched evaluation (how are negatives sampled? from what universe? are self-edges excluded? are TF→target directionality constraints enforced?).

### 4) CSSI remains partially circular / exploratory in real data
You do flag this (“same dataset for layer identification and reporting”), but the real-data CSSI claim still risks being read as “post hoc layer picking.”

**Fix:** tighten the presentation:
- Make split-half (or cross-validation) the *primary* real-data CSSI result, and put the exploratory full-data layer curves as descriptive.
- State explicitly whether **hyperparameters/variants** (max/mean/range/deviation) were selected on held-out data.

### 5) Scope is extremely broad; uneven maturity across the 9 supporting analyses
The paper currently tries to be:
- a critique of attention-as-regulation,
- a new diagnostic method (CSSI),
- a benchmark audit (null models + invariance),
- a mediation/patching bias paper,
- a cross-species conservation paper,
- a pseudotime audit,
- a batch leakage audit,
- and a calibration/conformal paper.

Individually, several analyses are interesting, but in aggregate this reads like an “all my experiments” compendium. Top venues often reward focus.

**Suggestion:** decide what the *one-sentence thesis* is, then demote the rest to appendix/supplement, or explicitly frame as a “diagnostic framework” paper with a minimal set of *validated* core components.

### 6) The multiple-testing story may confuse reviewers
A framework-wide BH correction across **47 tests** is conservative, but it’s also unusual to apply a single correction across heterogeneous, partially exploratory analyses. Reviewers may suspect either:
- cherry-picking (if some are corrected and others not), or
- over-penalizing and rendering the inferential claims moot.

You partially address this, but I’d recommend:
- Distinguish **confirmatory** vs **exploratory** analyses.
- Apply FDR within families (e.g., perturbation validation tests as one family), while being transparent about total testing.

### 7) Reproducibility / “submission ready” issues
Not scientific blockers, but top venues do care:
- Zenodo DOI is “to be assigned”; repo is “reviewer access.” That’s fine for a draft, but for submission you need stable artifacts.
- The paper claims very specific runtimes and datasets; ensure the repo actually reproduces every figure/table deterministically.

---

## Specific technical questions a reviewer will ask (preempt them)
1. **Definition of co-occurrence/co-expression metric** used to correlate with attention: is it Pearson/Spearman across cells, dropout co-detection, or something else? How does it relate to the attention edge score definition?
2. **Edge universe + negative sampling**: From what set are negatives drawn? Are TFs restricted to known TF lists? Is directionality respected? Are genes outside HVG set excluded? These choices can drive AUROC.
3. **Attention edge scoring**: exact pooling across heads/layers; are rows/cols normalized? Are diagonals removed? Is attention symmetrized or directional?
4. **Expression matching**: which statistics exactly (mean of rank-value encoded? raw counts? log1p normalized?) and is matching done within tissue/cell-state? What happens when matches fail (replacement, fewer negatives)?
5. **Bootstrap**: resampling over edges assumes edge independence, which is false (shared TFs/targets). For n=27, dependence can distort CIs.

---

## Concrete improvements that would most increase the score
1. **Unify the matched-negative experiment numbers everywhere**; add a single table: protocol → AUROC (CI) for attention, expression baselines, attention-marginal baseline.
2. **Strengthen the matched-negative evidence**:
   - replicate across ≥2 tissues *with the same matching protocol*;
   - add leave-one-TF-out stability; and/or
   - move to a larger positive set (even if noisier) and show consistent gap.
3. **Make CSSI validation cleanly held-out** on real data (primary result), with exploratory layer scans clearly labeled.
4. **Trim or reframe supporting analyses** so the paper doesn’t look like 12 partially-connected mini-papers.
5. **Clarify evaluation design** (edge universe, negative sampling, directionality, independence).

---

## Bottom line
The paper is now in the “strong workshop / borderline conference” range for an ML+bio venue: the critique is meaningful, the controls are improving, and the matched-negative design is the right direction.

What still blocks top-venue acceptance is **not** the high-level thesis, but the **strength and cleanliness** of the pivotal evidence (small n, dependence, and inconsistent reporting), plus the **scope/positioning** that makes it harder for reviewers to see a single, sharply validated contribution.
