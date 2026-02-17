# Codex Review Round 10

## Overall rating: 3/10

Strong negative audit paper, weak positive method paper. As written, it overstates causal/mechanistic recovery and under-delivers on independent validation.

## (1) Underpowered results
1. Scaling uses only 2 repeats per cell-count condition, admits overlapping 95% CIs. Yet asserts "critical threshold."
2. Mediation: only 16 run-pairs — pilot, not prevalence estimate.
3. Perturbation validation: zero findings surviving framework-level correction.
4. Real attention validation on 497 brain cells while arguing tissue heterogeneity is a dominant confound.
5. Disconnect between statistical outcomes and rhetorical certainty.

## (2) Circular validation
1. Correlation-like edge scores against TRRUST/DoRothEA — acknowledged but undermines headline inference.
2. Layer-selection/CSSI: discovery and reporting on same dataset.
3. "Anti-circular" arguments still internally selected: same tissue, same benchmark family.

## (3) Alternative explanations
1. Reference mismatch (brain-specific vs non-neural curation)
2. Technical/compositional confounding (donor/batch/cell-type leakage)
3. Objective mismatch (attention tracks co-expression, not causality)
4. Selection effects ("best layer/head" with weak generalization)
5. Synthetic confirmation bias

## (4) Weakest result
Real-data CSSI as "constructive resolution" — gains at top layers essentially zero. Practical rescue claim carried by synthetic/oracle-label settings.

## (5) One experiment to fix
Fully external, preregistered, perturbation-grounded validation with frozen pipeline on independent dataset using perturbation-defined causal edges (not TRRUST/DoRothEA).
