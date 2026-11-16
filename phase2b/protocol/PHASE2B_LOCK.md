# Phase 2B analysis lock

Locked on 2026-11-14 before any model of partial R² on structural category.

The question is whether genes with a predicted principal-versus-alternative structural difference have a higher gene-level partial R² than genes whose predicted protein sequences differ without a supported structural difference. The association is not a causal effect.

## Categories

Categories come from the Phase 2A pair table. A pair is usable when `coordinate_note` is `version_matched`. A usable pair has a supported structural difference when its Pfam class is `principal_only`, `alternative_only`, or `mixed`, or when the SignalP presence values differ, or when the TMHMM span counts differ.

- S2: at least one usable pair with different sequences and a supported structural difference.
- S1: at least one usable pair with different sequences and no such structural difference, and no S2 pair.
- S0: usable pairs exist, every usable pair has identical sequences, and the gene has no unresolved differing pair.
- ANNOTATION_UNRESOLVED: every remaining candidate gene, including genes with no measured alternative protein sequence.

S2 outranks S1, and S1 outranks S0. Missing annotation is not placed in S1 or S0.

The historical Phase 2A counts are checks, not targets to be forced: 10,795 candidate genes, 39,490 measured protein-coding transcripts, 10,742 genes with a principal SwissProt or TrEMBL accession, 9,608 genes with a different predicted sequence, 4,829 genes with a supported structural difference, and 10,228 genes with a Pfam accession on the principal protein.

## Primary sample

Literal A and Literal B are fit separately. A row is one successful Phase 1A fit (`fit_status = OK`) for an S1 or S2 gene. Partial R² is the outcome copied from the frozen fit file. M0 and M1 are not refit.

The primary parameter is the coefficient on the S2 indicator in

partial R² ~ S2 + cancer indicators + prespecified covariates,

with BRCA as the reference cancer. Uncertainty is cluster-robust at the gene identifier, using the Liang–Zeger meat and the finite-sample factor G/(G−1) × (N−1)/(N−K). A positive coefficient means a higher average partial R² in S2 inside this model.

## Covariates

The candidate covariates, in descending priority, are:

1. Principal protein length.
2. Number of Pfam accessions on the principal protein.
3. Fraction of the principal sequence covered by merged Pfam intervals.
4. Number of measured alternative protein-coding transcripts.
5. Number of measured isoform rows.
6. Complete-case count N for that gene and cancer.
7. Protein-observation fraction N / candidate cases.
8. Median log2 reconstructed gene expression in that cancer, using the same literal inverse as the outcome.

SignalP and TMHMM values are not covariates. They are part of the S2 definition.

Before the outcome regression, covariates missing on more than 5% of the estimation rows are removed. Variance-inflation factors are then computed on the remaining covariates, without partial R². While any VIF exceeds 10, the lowest-priority covariate among those above 10 is removed and the factors are recomputed. The S2 indicator is never removed. The surviving set is the adjustment set for every later sensitivity of that transformation.

## Practical rule

The practical magnitude is 0.005 partial R². It was set from the Phase 1B description of an approximately 0.01 gap between the lowest and highest cancer means. Half of that gap is the smallest difference treated as comparable to the heterogeneity already described. A confidence interval that excludes zero is not sufficient.

## Common support

The comparison is not qualified when either class has fewer than 500 genes with at least one successful fit, or when the absolute standardized difference in principal length exceeds 1 and the overlapping central window contains fewer than 300 genes in either class. The window runs from the larger of the two class-specific 10th percentiles to the smaller of the two 90th percentiles, applied to length and to Pfam count together.

## What is reported around the primary coefficient

- Unadjusted stacked difference from a cancer-indicator model without covariates.
- Gene-equal-weight difference of each gene’s unweighted mean partial R².
- Cluster-bootstrap percentile interval, 1,000 draws, seed 20261114.
- Leave-one-cancer-out refits.
- Length-quintile unadjusted differences, then the unweighted mean of strata that contain at least 30 genes in each class.
- The overlapping-window refit.
- A refit after removing genes with only one measured alternative protein-coding transcript.
- A one-step cluster influence screen. Genes with an absolute one-step change in the S2 coefficient above 2/sqrt(G) are removed in one sensitivity.

Within each cancer, the adjusted S2 coefficient uses the same covariates except the cancer indicators, with HC1 intervals. Benjamini–Hochberg adjustment is applied to those ten cancer tests within each transformation. A rank comparison of gene-mean partial R² is secondary.

Domain truncation and domain-boundary change are not estimated. Phase 2A qualified exact sequence equality and did not qualify a residue alignment. Shared Pfam accessions are counted as gained or absent. A change in the raw interval length on two unaligned sequences is not called a truncation.

## Decision labels

Apply the first matching label.

1. TECHNICAL_NO_RESULT when the historical counts do not reproduce or the primary design cannot be estimated.
2. PRIMARY_COMPARISON_NOT_QUALIFIED when the common-support rule fails.
3. STRUCTURAL_ASSOCIATION_SUPPORTED when, for both inverses, the adjusted estimate is at least 0.005, its cluster-robust interval lies entirely above zero, every leave-one-cancer-out estimate stays positive, and the overlapping-window estimate keeps the same sign and stays above zero.
4. STRUCTURAL_ASSOCIATION_TRANSFORMATION_DEPENDENT when the two inverses disagree in sign, or only one inverse meets the supported rule.
5. STRUCTURAL_ASSOCIATION_CANCER_SPECIFIC when both inverses show at least one cancer with an adjusted estimate of at least 0.005 and an interval above zero, and at least one other cancer with an adjusted estimate of at most zero and an interval whose upper bound is below 0.005.
6. NO_PRACTICAL_STRUCTURAL_ASSOCIATION otherwise.

Phase 2C is not authorized under any label.
