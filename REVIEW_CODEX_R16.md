# Codex Adversarial Review — Round 16 (updated `main.tex`)

**Overall rating (1–10): 8.1 / 10**

This revision is substantially stronger than prior rounds: the headline results are now **numerically consistent**, the **evaluation methodology is mostly specified end-to-end**, and the paper is more honest about **post-hoc selection** and **non-independence**. I can largely reconstruct how the key AUROCs were obtained from the text alone (especially the expression-matched evaluation), which was the central reproducibility gap previously.

The main remaining issues are (i) a couple of **internal/architectural inconsistencies** (notably Geneformer layer/head counts in the scaling-methods paragraph), and (ii) a few places where **“what universe of candidate edges”** or **exact scoring definitions** are still implicit rather than explicitly written in one canonical “Evaluation Protocol” block.

---

## What improved (relative to the prior 7.8 version)

1. **Number consistency fixed.** The expression-matched result is consistently reported as **attention AUROC 0.646** vs **baseline 0.522**, with **CI [0.539, 0.747]** repeated consistently (Introduction/Core Contribution 3, Cross-tissue section, Discussion, Conclusions).

2. **Evaluation methodology is now mostly reproducible.** The paper now states:
   - **Negative universe:** “all non-TRRUST gene pairs among the 1,200 HVGs.”
   - **Matching criteria:** for each positive edge, sample up to 50 negatives with both TF and target matched on **mean expression** and **detection rate** within **±20%**.
   - **Uncertainty:** **1,000 bootstrap resamples** of the combined positive+negative set.
   - **Caveat:** positives not fully independent (shared TFs), potentially inflating effective sample size.

3. **CSSI held-out split-half is foregrounded** with explicit “post-hoc selection concern” language and a concrete split-half procedure.

Net effect: key claims are more falsifiable and less “benchmark AUROC = regulation” naïve.

---

## Internal consistency check (key numbers & cross-references)

### Core triad
- **Finding 1 (pooled attention dominated by co-expression):** consistent framing; multiple sections repeat ρ(co-occurrence)≈0.31–0.42 and ρ(ground truth)≈0, with pooled AUROC≈0.5.
- **Finding 2 (CSSI/layer localization):** scGPT-18L later layers L13–L14 AUROC **0.694–0.706** repeatedly consistent (abstract, core finding section, discussion, conclusions).
- **Finding 3 (cross-tissue invariance + null models):** scGPT-12L cross-tissue AUROC ≈0.72 TRRUST / ≈0.84 DoRothEA, null models match/exceed (Table `dumb_baselines`).

### Expression-matched evaluation
- All occurrences now agree on **0.646 vs 0.522** and CI **[0.539, 0.747]**.

### Remaining tension (mostly rhetorical)
- The text sometimes says “across scGPT … ρ≈0 with TRRUST” while later emphasizing deep-layer AUROC≈0.70. You do qualify this as “**pooled attention**” vs “layer-stratified,” but I’d still suggest tightening one sentence to avoid an easy reviewer nitpick:
  - e.g., say explicitly “**pooled across layers/heads** yields ρ≈0 / AUROC≈0.5; layer-specific scores can be above chance.”

---

## Major issues (still blocking a 9+)

### 1) Apparent Geneformer architecture inconsistency in Methods (reproducibility + credibility)
In **Methods → Scaling Behavior Analysis**, you write:

> “Attention-derived edge scores were extracted from **Layer 13** (pooled across all **18 heads**) … using Geneformer V1-10M …”

This reads like a **copy-over from scGPT-18L** (18 layers/heads) rather than Geneformer’s actual configuration (commonly 12-layer BERT-like). If this is wrong, it’s a serious “small detail” that undermines trust in the scaling experiment spec.

**Action:** verify Geneformer’s layer/head count and correct the sentence (or explicitly explain that you used a particular internal layer index/variant). This is the clearest remaining internal inconsistency I noticed.

### 2) Evaluation protocol is distributed across sections (harder to replicate than necessary)
You now give the critical details, but they’re spread between **Cross-tissue** and **Discussion/Intro**. For reproducibility, it would help to have a single compact block:

- candidate edge universe (TF×target? all directed pairs? how do you treat direction?)
- what counts as a “negative” in unmatched evaluation (all non-benchmark pairs? sampled? ratio?)
- whether AUROC is computed on *directed* edges or undirected gene pairs
- exact attention score definition (average over heads? over cells? how do you aggregate per layer?)

Right now, a careful reader can infer much of this, but a *reproducer* still has to hunt.

### 3) Small-edge-set cross-tissue TRRUST evaluation remains fragile
You do acknowledge the limitation (27–28 edges, 14 TFs) and report wide CIs. That’s good. But the paper still leans on cross-tissue invariance rhetorically.

This is now framed as “diagnostic of degenerate evaluation,” which is reasonable, but reviewers may still object that invariance could be an artifact of the tiny set. The DoRothEA 123-edge set helps, but consider strengthening with:

- a sensitivity analysis showing invariance persists under random subsampling of DoRothEA to 27 edges (or conversely, show how noisy the TRRUST estimate is)
- explicitly reporting **the exact number of positives per tissue** (27 vs 28) and why it differs.

---

## Minor issues / polish

- **Terminology drift:** sometimes “co-expression correlation” vs “co-occurrence” are used; ensure you define exactly what the co-expression statistic is (Pearson/Spearman? on what transformed data? pooled across cells?)
- **Multiple-testing correction:** you state 47 tests and some corrected p-values that “survive.” The main text doesn’t show the raw p-values for those particular claims (likely in appendix). That’s okay, but ensure appendix table is easy to audit.
- **CSSI split-half (real data):** you report AUROC 0.619 vs 0.539 and p≈0.02. Specify the exact test (permutation on edges? bootstrap?) and the sampling unit (edges) to avoid confusion.

---

## Reproducibility verdict

**Mostly reproducible** for the key headline results (especially the expression-matched analysis), *provided the code repository exists as claimed.* The manuscript alone is now close to “sufficient,” but a couple of details (Geneformer architecture line; centralized evaluation protocol) should be fixed for high-confidence independent reproduction.

---

## Bottom line

This is now a coherent, internally consistent narrative with appropriate caveats, and the methodological additions materially improve auditability. Fix the Geneformer layer/head inconsistency and consolidate the evaluation protocol, and the paper will likely clear a 9/10 on “internal consistency + reproducibility.”
