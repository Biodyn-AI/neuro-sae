# Codex Adversarial Review — Round 8 (R18)

**Manuscript:** `main.tex` (updated)

## Overall rating (1–10)
**8.5 / 10** (up from 8.3).

You’re in the zone where improvements are increasingly about tightening claims, reducing self-contradictions, and pre-empting reviewer attack lines rather than adding new analyses. I do think you’re approaching **diminishing returns**: another full iteration might net +0.1–0.3 unless you address a couple of high-leverage narrative/rigor issues called out below.

## What improved since last round (noticed in the text)
1. **Title softening** (“Primarily Reflects Co-Expression”) lowers the “overclaim / strawman” risk and better matches the nuance later (esp. expression-matched AUROC > baseline).
2. **Small-n caveat** before cross-tissue results helps; you now explicitly flag the 27-edge TRRUST set and CI width.
3. **[ATTN]/[CORR] tags** are genuinely helpful and reduce a likely reviewer complaint (“you call it attention work but half the paper is correlation”).
4. **p vs q convention** is clearer; the framework-level correction is now explicitly described.
5. **Cell threshold softened to 150–300 range** matches the scaling section’s quantitative story and avoids a brittle “magic number”.

## Main strengths (as a reviewer)
- **Co-expression dominance finding** is clear, multi-model, and framed as a mechanistic/evaluation pitfall rather than “models are bad”. This is the best contribution.
- **Null-model benchmarking** (gene-level baselines matching AUROC) is the right “kill shot” for naive benchmark interpretations.
- **Expression-matched negative sampling** is an important “save” that prevents the paper from becoming purely negative: it supports the nuanced conclusion “pairwise structure exists but modest”.
- The paper is unusually explicit about **boundary conditions** and **what the analyses do/do not test**.

## High-leverage vulnerabilities / adversarial critique
These are the points I expect a tough reviewer to press hardest.

### 1) Internal tension: “pooled attention AUROC≈0.5” vs “per-layer AUROC≈0.72” vs “null-models ≥0.75”
You do try to reconcile this (Discussion: “Reconciling…”), but it still reads *confusing on first pass* because the headline metrics shift across sections/models/setups.

**Adversarial reading:** the paper may appear to “move the goalposts” (switching models, edge sets, vocabularies, tissues), so the reader can’t tell whether attention is (i) useless, (ii) moderately useful, or (iii) just a proxy for expression.

**Fix (high ROI, mostly writing):** add a small **one-paragraph “metric map”** early (end of Intro or start of Results) that explicitly enumerates:
- which claim uses which model (scGPT-18L vs scGPT-12L vs Geneformer),
- which evaluation edge counts (8,330 vs 27–28 vs 123 vs 483),
- which score pooling (pooled across layers/heads vs per-layer),
- and which negative sampling (unmatched vs expression-matched).

This can be done without adding data; it reduces the “shell game” impression.

### 2) Cross-tissue section: the *small evaluation set* is so small it risks being dismissed
You already caveat this, but the narrative still leans on the invariance observation.

**Adversarial reading:** 27 edges (14 TFs) is small enough that invariance across layers/tissues could be an artifact of (a) a few TFs dominating, (b) correlated edges per TF, and (c) bootstrap on non-independent edges. A reviewer could argue: “this is not a cross-tissue replication; it’s a toy sanity check.”

**Mitigation (writing + one robustness check if available):**
- You already mention non-independence of positives; consider moving that caveat *earlier* and making it more prominent.
- If you have it (maybe already computed): report **TF-block bootstrap** or **leave-one-TF-out** sensitivity for the expression-matched AUROC gap. Even one sentence like “the attention>baseline gap persists under leave-one-TF-out” would blunt this.

### 3) Multiple testing correction story may invite criticism
You claim “Key findings survive correction: scaling failure (q=0.011), mediation non-additivity (q=0.003), cross-species conservation (q=0.001)”. In the main text, some p-values are astronomically small; others are modest and then declared non-significant after BH.

**Adversarial reading:** A reviewer may ask for a **table** listing the 47 tests with p/q, because otherwise the correction feels like a rhetorical shield rather than an auditable procedure.

**Fix:** You already `\input{appendix_statistical_tests}`; ensure it **contains the full test ledger** and that the main text points to it when making any “survives correction” statement.

### 4) Claims around “regulatory” vs “pairwise” signal need careful wording
You now say expression-matched evaluation provides “direct evidence of genuine pairwise structure beyond gene-level confounds” (good), but the manuscript sometimes drifts into “regulatory information” language.

**Adversarial reading:** “Pairwise structure” is not the same as “regulatory”. It could be chromatin co-accessibility proxies, cell-type mixture, shared pathway activity, etc.

**Fix:** replace any remaining “pairwise regulatory information” phrasing with **“pairwise benchmark-discriminative structure”** unless you can justify causality.

### 5) Some supporting analyses are (by design) correlation-based; you defend this, but it still dilutes focus
The [CORR] tagging helps a lot. Still, a reviewer might ask: “why include ortholog/pseudotime/batch/calibration if they’re not attention?”

**Fix:** In Supporting Analyses preamble, add 1–2 sentences making the inclusion principle explicit:
> These analyses are not claims about attention; they quantify *the baseline difficulty and confound structure* any edge-scoring method (including attention) must overcome.

You already say something like this; making it even more direct will pre-empt the complaint.

## Minor issues / polish opportunities
- **Methods: scGPT-12L cross-tissue TRRUST edges “27–28”**: explain why it varies (mapping differences per tissue?) so it doesn’t look sloppy.
- **AUROC reporting precision**: some tables show 0.72 everywhere; consider reporting 3 decimals to avoid “copy-paste” skepticism.
- **“dedicated GRN methods … achieve identical near-random performance on brain tissue”**: make explicit that this is for *that dataset + that preprocessing*; otherwise experts in GRN inference may push back.

## Are we at diminishing returns?
**Yes, mostly.** The paper is now coherent and defensible. The next improvements with the highest marginal payoff are:
1) a “metric map” paragraph/box to remove the perceived inconsistency across sections,
2) one robustness line for the tiny 27-edge cross-tissue evaluation (TF-block or LOOTO),
3) ensuring the appendix provides a full auditable p/q ledger.

If you do those, I could see the rating moving to **8.7–9.0**. Without them, further tinkering will likely be cosmetic.
