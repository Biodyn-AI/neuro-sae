# Codex Adversarial Review R13 (Round 3)

**Paper:** *Attention in Single-Cell Foundation Models Reflects Co-Expression, Not Regulation: A Diagnostic Framework with Null-Model Benchmarking*

**Overall rating (1–10): 7.0/10**

This revision is notably clearer and the internal contradiction around “regulatory recovery” is largely resolved. The reframing toward **benchmark-discriminative signal** + **degenerate evaluation** is directionally correct and makes the paper more honest. That said, for a *top-venue* acceptance the work still needs to (i) tighten causal/empirical claims, (ii) fix evaluation-design weaknesses that are currently presented as “future work,” and (iii) improve methodological specificity so results are auditably reproducible from the text.

Below are the key remaining blockers and concrete fixes.

---

## 1) Major acceptance blocker: the core “degenerate benchmark” claim needs a *controlled* evaluation, not just a post-hoc baseline

You correctly show that gene-level nulls (detection/mean/variance products) match or beat attention AUROC on curated benchmarks. However, you then state that **expression-matched negative sampling is needed** and stop there.

For a top venue, reviewers will likely say: *“Then do it.”* As written, the paper’s strongest conclusion (“AUROC on curated benchmarks is invalid evidence of regulation”) still rests on an **argument** plus a **set of un-matched baselines**, rather than a controlled test.

**What to add (minimum):**
- Implement **expression-matched negative sampling** (match TF and target marginal statistics separately; not just the pair product). At minimum match on detection rate and mean; ideally include variance too.
- Report AUROC **and** AUPRC under matched negatives for:
  - attention pooled,
  - best-layer / best-head,
  - correlation baseline,
  - your dumb gene-level nulls.
- Show whether attention retains *any* advantage when negatives are matched.

**Why this matters:** Without this, a reviewer can argue that your null models are “too dumb” or that product baselines overfit the benchmark idiosyncrasies. Matching removes that debate.

---

## 2) scGPT-12L cross-tissue section is statistically too small to bear the narrative weight

You do acknowledge the TRRUST evaluation set is **27–28 edges** (14 TFs). Good. But the section still uses that result as a central pillar of the “cross-tissue invariance → degeneracy” argument.

With 27 positives, AUROC is extremely noisy and also vulnerable to per-TF composition effects (a few TFs can dominate). The paper should not lean on this as strongly without strengthening the evaluation.

**Fixes:**
- Expand beyond the 27-edge TRRUST slice. If mapping is the issue, explicitly say why only 27 edges remain after filtering and what the universe is.
- Provide **per-TF AUROC** (or leave-one-TF-out AUROC) to show invariance is not an artifact of a few TFs.
- Consider adding additional references (e.g., OmniPath, RegNetwork, ChEA) *with the same matching control*.

---

## 3) The paper still reads like it sometimes conflates three distinct objects: co-expression, attention edges, and “regulatory signal”

You now say “benchmark-discriminative signal,” which helps. But some sentences still imply that later layers “encode regulation,” then later retract to “gene-level prominence.” Example: the CSSI results call it “regulatory signal,” then immediately caveat that it is not regulatory structure.

**Fix:**
- Use one consistent hierarchy:
  1) *attention correlates with co-expression* (pairwise correlational structure),
  2) *benchmarks are enriched for prominent genes* (marginal structure),
  3) *AUROC can be high for non-regulatory reasons*.
- Reserve the word “regulatory” only for claims validated against perturbations, time-course, or properly controlled benchmarks.

This is mostly editorial, but top-venue reviewers will be very sensitive to overstated mechanistic language.

---

## 4) Methodological clarity gaps: a skeptical reviewer cannot reconstruct your edge scoring and correlation analyses from the text

Several core computations are still underspecified:

**Attention edge score definition:**
- How exactly do you aggregate attention?
  - across heads (mean? max? weighted?),
  - across layers (when pooled),
  - across cells (mean of per-cell matrices? median? rank aggregation?),
  - normalization (row-stochastic attention already sums to 1; do you symmetrize? remove diagonal?).
- How do you map attention tokens to genes when using HVG subsets? What happens to genes not in vocabulary?

**Co-expression / co-occurrence definition:**
- You state correlation with “gene co-occurrence.” Is that:
  - Pearson/Spearman across cells of expression vectors,
  - binary detection co-occurrence,
  - or something else?
- If it’s detection-based, do you control for library size / sequencing depth?

**Ground truth mapping:**
- What is the negative edge universe for AUROC?
  - all TF–gene pairs within HVGs?
  - only TFs present in TRRUST/DoRothEA?
  - are self-edges excluded?

These choices can swing results dramatically and will be asked.

**Fix:** add an explicit “Edge scoring + evaluation protocol” subsection with pseudocode-level clarity.

---

## 5) Multiple-testing correction across “47 framework-level tests” is unusual and may invite criticism

A global BH correction across heterogeneous analyses (scaling, pseudotime, calibration, cross-species) may be seen as:
- either *unnecessarily conservative* (masking true positives),
- or *methodologically incoherent* (different families of hypotheses).

Right now it’s asserted, but not justified as a principled choice vs. per-section families.

**Fixes:**
- Report both: per-analysis correction (within each family) **and** global correction in a supplement.
- Clarify what counts as a “test” and how you prevent post-hoc cherry-picking of which p-values enter the family.

---

## 6) The “baseline comparison” conclusion is currently too strong / possibly misleading

You write that GENIE3/GRNBoost2 near-random AUROC suggests “tissue-specific characteristics or inappropriate benchmarking.” That may be true, but this section can be attacked because:
- GENIE3/GRNBoost2 are not designed to recover TRRUST edges specifically in a heterogeneous adult brain snapshot;
- hyperparameters / candidate TF lists / pre-filtering can make them look arbitrarily bad;
- evaluation set incompleteness can flatten differences.

**Fix:**
- Reframe as: “Under this benchmark design, multiple methods collapse to ~0.5, consistent with benchmark limitations and/or mismatch to biological context.”
- Add at least one condition where baselines *should* work better (e.g., a perturbation dataset, or a simpler homogeneous cell line subset). A positive control is essential.

---

## 7) CSSI: strong idea, but needs sharper positioning and stronger held-out validation on real data

You correctly flag circularity and call the real-data layer selection exploratory, but then still report layer peaks prominently.

**Fixes:**
- Promote the held-out protocol from “split-half” to something more robust:
  - stratified splits by donor/cell type, or
  - evaluation on a second brain cohort / independent tissue.
- Clarify what CSSI variants are (“mean/range/deviation”)—these appear without formal definition. Right now it reads like unpublished heuristics.

---

## 8) Cross-model section: “stable attention patterns are stably uninformative” is catchy but needs care

The Geneformer stability claim is interesting, but a reviewer may ask whether stability is:
- stability of attention matrices across cells,
- across bootstrap samples,
- across tissues,
- or across runs.

Also, cosine similarity 0.979±0.001 is so high that it suggests either a heavy averaging procedure or a restricted subset. That’s fine, but must be explained.

**Fix:** define the stability metric precisely and relate it to the edge scoring pipeline.

---

## 9) Presentation / narrative tightening

A top-venue reviewer will likely praise the negative-results framing, but still ask for:
- A clearer “what should practitioners do tomorrow?” box that is *not* just “be cautious.”
- A concise definition of the paper’s main contribution: is it (a) debunking attention-as-regulation, (b) proposing CSSI, or (c) diagnosing benchmark pathology? Right now it’s all three.

Consider explicitly positioning this as a **benchmark/evaluation paper** with CSSI as a diagnostic tool, rather than a mechanistic-interpretability paper that sometimes sounds like it is finding real regulatory circuits.

---

## Summary: what still needs fixing for top-venue acceptance

1. **Do expression-matched negative evaluation** (not just propose it). This is the single biggest upgrade.
2. Strengthen or de-emphasize conclusions relying on **27-edge** TRRUST slices; add per-TF robustness.
3. Add a crisp, auditable description of **edge scoring + evaluation protocol**.
4. Clarify language so “signal localization” never reads as “regulatory recovery.”
5. Rework multiple testing correction rationale and reporting.
6. Add at least one **positive control** setting where regulation recovery should be possible / meaningfully different.

With these, I think the paper can move from “interesting and mostly correct critique” to “decisive evaluation framework paper” suitable for a strong ML-for-biology venue.
