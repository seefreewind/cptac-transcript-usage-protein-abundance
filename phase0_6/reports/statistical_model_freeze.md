# Statistical model freeze

Status: frozen on 2026-11-14. Decision: **PARTIAL**. Phase 1 has not been started. This document does not authorize it.

The gene-level RNA–protein comparison is feasible. The primary transcript-usage representation is defined. Its number is still a pair, because the isoform inverse is unresolved. That is the qualification that keeps this phase at PARTIAL. A GO decision would require one inverse that the processing record identifies. A NO-GO decision would require that a usage proportion cannot be defined. Both literal inverses define one.

## What is carried forward

The hypothesis, outcome, covariates, 71-case gate, cancer-specific fits, random-effects summary, and negative-result rule are the ones written in `phase0_5/reports/statistical_analysis_plan.md` on 2026-11-13, before any test of transcript usage against protein. This freeze adds the literal-inverse pair required by Phase 0.6 and the representation choice made from the form of the estimand.

Within a cancer, among case-matched tumors, the principal-transcript proportion from the isoform matrix is tested for additional explanation of gene-level protein abundance beyond log2 of the isoform-derived gene sum, WES purity, and gene-level copy number. The result is a statement about gene-level protein abundance. It is not a measurement of protein function, and it is not an isoform-specific protein abundance.

## Names of the inverses

Phase 0.6 literal A is X = 2^Y, and a stored 0 becomes 1. Phase 0.6 literal B is X = 2^Y − 1, and a stored 0 stays 0. Literal B is the same function as Phase 0.5 Candidate B.

Phase 0.5 Candidate C keeps a stored 0 at 0 and sets positive values to 2^Y. That clipped reading is a labeled sensitivity. It is not a third primary, and it is not literal A. The primary unresolved pair is literal A together with literal B. Details and the evidence are in `transformation_resolution_report.md`.

## Primary representation

The primary explanatory variable is the raw proportion of the GENCODE v34 principal transcript:

P = X_principal / G, where G is the sum of X over the isoform rows of that gene.

The principal transcript is the unique MANE Select transcript when the gene has one, and otherwise the unique APPRIS principal-1 transcript. MANE Select is kept when the tags disagree. The 1,001 genes with no unique principal stay out. Thirty-seven principal transcripts are absent from the isoform matrix, so those genes have no contrast in these files. The contrast is defined for the 10,795 principal genes that have at least two isoform rows.

Centered log-ratio transformation of the gene’s transcripts is a sensitivity. Its pseudocount is half the smallest positive linear value in that gene, so the pseudocount is computed inside whichever inverse is being fit. The isometric log-ratio is not a primary model. The reason is the estimand: one prespecified degree of freedom, interpretable as the share of the principal transcript. No association with protein was used to make this choice.

## Models

For each cancer, each eligible gene, and each literal inverse, using the same cases in the two fits of that inverse:

- M0: protein ~ log2(G) + WES_purity + CNV
- M1: protein ~ log2(G) + P + WES_purity + CNV

The outcome stays on the published protein log2 scale. CNV is the gene-level log2 copy-number ratio, matched on the versioned ENSG identifier. A case absent from the CNV matrix is dropped. WES purity is the primary purity covariate. A case with WES purity missing is dropped. The combined ESTIMATE score is the sensitivity replacement for purity. No phenotype table contains a batch or plex column, so none is added.

The primary gene covariate is log2(G). The published RSEM value replaces it only in a sensitivity. RSEM is a different estimator, scaled to an upper quartile of 1,500 before log2.

The primary estimand, until one inverse is identified, is the pair of partial R² values of P. For each inverse, ΔR² = R²(M1) − R²(M0), with the adjusted R² of each model and the nested comparison. The reported result is the within-cancer distribution of that partial R², then a random-effects summary across cancers. Cancers are not stacked into one regression of 1,027 cases.

## Complete cases

Inside one inverse, M0 and M1 use the same cases.

The case lists of the two inverses differ where every isoform row is stored 0. Literal A gives each of those rows the value 1, so G equals the number of rows and P equals 1 over that number. Literal B gives G = 0, and the case is out of M0 and M1. The count of those gene-samples is reported with the fit. A sensitivity of the A fit drops them, which puts the A case list on the same gene-samples as B.

G = 0 is not imputed. Protein NA is not imputed. The four OV cases missing both CNV and WES purity stay out: `01OV007`, `01OV018`, `11OV002`, `17OV013`. The other triple cases missing purity or CNV stay listed in the Phase 0.5 cohort report and stay out of the fits that require the missing covariate.

## Inclusion

Cases are the exact intersection of tumor gene RNA, tumor isoform RNA, and tumor protein. The Phase 0 counts are the reference: BRCA 121, CCRCC 103, COAD 96, GBM 99, HNSCC 108, LSCC 108, LUAD 110, OV 82, PDAC 105, UCEC 95. The match is case-level.

A gene enters the primary distribution when at least 71 cases have protein, log2(G), P, WES purity, and CNV all observed, and the gene has at least two isoform rows. The same minimum applies in every cancer, including OV. OV is summarized on its own. In the Phase 0.5 zero-preserving count, OV has 4,250 genes at or above 71, with median eligible N 78 and mean eligible N 61.4. Genes below 71 are left out of the primary distribution. If a cancer contributed no gene at that gate, the report would say so. OV does contribute under the published count.

Those published counts treat a gene-sample as eligible when G > 0 under the zero-preserving sum. Every gene in that table remains eligible under literal A. Literal A can add all-zero gene-samples and can therefore move additional genes across 71. That recount is part of the fit, and it has not been run. The planning anchor is unchanged: for one added predictor, λ = (1.959964 + 0.841621)² = 7.848880, and partial R² = f² / (1 + f²) with f² = λ / N. At the published median eligible N, that detectable partial R² is 0.067–0.091. The historical anchor of 0.10 is a planning number.

## Inference and sensitivities

The primary claim is the within-cancer mean partial R² and its cross-cancer random-effects summary. A gene-level screen, if shown, uses Benjamini–Hochberg within that cancer and is labeled as a screen. A case-resampled interval is computed within the cancer. Any predictive check resamples cases within a cancer, and the gene filter is fit inside the training cases.

A within-cancer mean partial R² whose interval covers values below 0.02, or a cross-cancer summary that does so, means this prespecified contrast adds no practically meaningful increment beyond gene-level RNA, purity, and CNV in these data. That result is about this contrast.

Prespecified sensitivities, each labeled as such:

- the same partial R² with the RSEM value in place of log2(G)
- the A fit after dropping all-zero gene-samples
- Phase 0.5 Candidate C, the clipped reading
- centered log-ratio in place of P
- ESTIMATE combined score in place of WES purity
- a refit that drops CNV only when that gene’s variance-inflation factor for CNV exceeds 10
- a refit that drops a case whose absolute Pearson residual in M0 exceeds 4; the primary fit keeps the case
- the same partial R² when the outcome is the mean of the observed phosphosites of that gene. Phosphosite matrices are supplementary. They are validated for all ten cancers. They do not identify an isoform.

A gene whose least-squares fit fails, or whose design matrix is rank-deficient, is counted as not estimated and is left out of the mean. It is not assigned a partial R² of 0.

## What stays outside the fit

Differential expression, differential transcript usage, isoform-switch lists, gene ranking, pathway enrichment, survival models, and machine-learning fits are outside this plan. TCGA, GTEx, confirmation cohorts, and DIA cohorts are not added. The six open SEPEP matrices are patient-level quantification groups. They contain no peptide sequence that can be checked against one isoform. CCRCC, COAD, GBM, and HNSCC SEPEP files returned HTTP 403 and were not retrieved by any other route. No isoform-specific protein abundance is claimed.

## Decision

PARTIAL. The matrices, the 1,027-case reconstruction, the principal-transcript rule, and the gene-level model are in place. The isoform inverse remains a pair, and the pair changes the principal proportion by more than 0.05 in 0.614–0.639 of the gene-samples where both proportions exist. Phase 1 starts only if the project owner authorizes it.
