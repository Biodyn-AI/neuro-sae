# JMLR Review (Codex)

## Overall recommendation: Reject (encourage resubmission after major restructuring and new experiments)

The manuscript is ambitious and contains substantial work, but in its current form it still has multiple rejection-level weaknesses: internal logical contradictions, over-claimed contributions relative to evidence, and missing experiments needed to support the central thesis for JMLR.

## Major issues

1. Core narrative is internally contradictory.
- The paper claims architecture-independent near-random failure of attention-based GRN inference (main.tex:36, main.tex:65, main.tex:617, main.tex:645).
- Baseline section says poor performance is tissue-specific rather than attention-specific (main.tex:242-main.tex:246).
- CSSI section later claims strong recoverable signal from real attention matrices (AUROC 0.694–0.706) (main.tex:535, main.tex:565, main.tex:567).
- Discussion reframes again to "unstratified methods fail" (main.tex:676).
These are materially different claims and currently coexist without reconciliation.

2. Claimed constructive contribution (CSSI) is over-stated relative to real-data evidence.
- Largest improvements are synthetic/correlation-proxy (main.tex:483-main.tex:515).
- On real attention, strongest layers show ~zero or slightly negative CSSI delta (main.tex:553-main.tex:555, main.tex:563).
- Conclusion still states CSSI "eliminates scaling failure" broadly (main.tex:722).

3. Multiple-testing correction is not fully auditable.
- Methods claim framework-wide BH correction over 47 tests (main.tex:109).
- Results mix raw and adjusted reporting.
- No complete test registry table (hypothesis, raw p, adjusted p, effect size).

4. Several central analyses are not attention-based despite paper framing.
- Cross-species, pseudotime, and batch sections explicitly use correlation-based edges.
- These are useful context but do not directly validate/critique mechanistic interpretability of FM attention.

5. Potential evaluation leakage in layer/head selection.
- Best layers/heads are identified and emphasized on the same benchmark used for reporting.
- No nested selection/held-out validation to show generalization of selected layers/heads.

6. "Fundamental limitation across architectures" claim is too strong.
- Deep evidence exists mainly for scGPT + Geneformer.
- scVI/C2S are explicitly preliminary and limited.
- This does not support broad "fundamental" language.

7. Internal bookkeeping inconsistency.
- Methods specify scaling points 200/500/1000/3000 (main.tex:79).
- Limitations later states only 3 points (200/1000/3000) (main.tex:715).

8. Perturbation validation is weak for a headline contribution.
- Effects are small/marginal and none survive BH.
- Yet perturbation consistency is positioned as a core pillar in title/abstract framing.

9. Writing clarity still harms credibility.
- Repetitive "important methodological note" blocks and repeated thesis restatements.
- Frequent rhetorical overreach ("first comprehensive," "fundamental," "eliminates") beyond demonstrated scope.
- Minor consistency/polish issues (e.g., TRRUST/TRRust spelling around main.tex:533).

## Missing experiments needed to avoid rejection

1. Strict held-out validation of layer/head/CSSI strategy (nested selection on one dataset, evaluation on separate dataset/tissue).
2. Attention-based versions of cross-species/pseudotime/batch analyses, or explicit de-scoping of these claims.
3. Full multiplicity appendix/table for all 47 tests.
4. Independent external validation beyond the current core datasets.
5. Fair stratified baseline study: apply comparable stratification to dedicated GRN methods to test whether CSSI benefit is attention-specific.

## Final verdict

Reject in current form.

Rationale: The manuscript has potentially publishable components, but current claim scope exceeds demonstrated evidence, and unresolved internal contradictions reduce trust in conclusions. A focused resubmission could be competitive after claim narrowing, stricter out-of-sample validation, and complete statistical auditability.
