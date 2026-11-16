# Statistical analysis plan

Status: prespecified on 2026-11-13, before any test of transcript usage against protein abundance. Phase 1 has not started. This document does not authorize it.

The historical Phase 0 anchor (residual partial R² = 0.10, two-sided α = 0.05, power 0.80, about 71 informative pairs for one added predictor) is a planning assumption. It is not an observed effect.

## 1. Primary hypothesis

Within a cancer, among case-matched tumors, a transcript-usage feature derived from the isoform matrix explains additional variation in gene-level protein abundance beyond a gene-level RNA covariate derived from that same isoform matrix, tumor purity, and gene-level copy number.

This is a statement about incremental explanation of gene-level protein abundance. It is not evidence that protein function was lost, and it is not a measurement of isoform-specific protein abundance.

## 2. Primary outcome

For gene g in cancer c, the outcome Y is the gene-level protein value in `*_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` for a case in the triple-matched set. The unit stays on the published log2 scale. Zeros are not recoded, and NA stays missing.

## 3. Primary explanatory variable

One prespecified contrast, not the full vector of transcript proportions.

Let the isoform matrix be the v1.2 file `*_RNAseq_isoform_FPKM_log2_Tumor.txt`. Do not sum those log2 values. Do not divide a log2 transcript value by a log2 gene value.

The linear abundance of transcript t is recovered by a documented inverse, applied to the stored value y:

- Candidate A, used when the full-matrix census finds negative values: stored y is treated as log2 of a positive FPKM, with 0 retained as a separate zero code. Then X = 0 when y = 0, and X = 2^y when y > 0. A negative y would also invert as 2^y. This candidate is unavailable if the census has no negatives and has values in (0, 1), because an unclipped log2 of an FPKM below 1 would be negative.
- Candidate B: X = 0 when y = 0, and X = 2^y − 1 otherwise. This is the inverse of log2(FPKM + 1).

The filename says `FPKM_log2`. A different object in the same bucket is explicitly named `log2(x+1)` and is not this file. The full-matrix census is now in `missingness_report.csv`. In all ten tumor isoform matrices the numeric values have no negatives, the minimum positive value is 0.01, every numeric value lies on a 0.01 grid, and values in (0, 1) are present. An unclipped log2(FPKM) cannot produce that pattern. Candidate B and Candidate C both remain arithmetically possible. They are fit as a prespecified pair. Neither is chosen after looking at its association with protein.

For a case and a gene, G is the sum of X across isoform-matrix rows of that gene. P_principal = X_principal / G. The principal transcript is chosen from GENCODE v34 before any protein model is fit. A gene with exactly one MANE Select transcript uses that transcript. If MANE Select and APPRIS principal-1 point at different transcripts, MANE Select is used. A gene with no MANE Select transcript uses APPRIS principal-1 only when exactly one transcript carries that tag. In this v34 file, 3,311 genes have more than one transcript tagged appris_principal_1. A gene in that set that also has exactly one MANE Select transcript uses the MANE transcript. A gene in that set with no MANE Select transcript is not assigned a transcript by file order. That unresolved set contains 1,001 genes. A gene with no unique principal transcript, a gene with fewer than two isoform-matrix rows, and a case with G = 0 are not in the primary regression for that gene. The count of G = 0 cases is reported and is not imputed.

The RSEM gene matrix is not the primary gene covariate. The primary gene covariate is log2(G) when G > 0. The published RSEM value is a sensitivity covariate only.

## 4. Primary model

For each cancer and each eligible gene, using the same cases in both fits:

- M0: protein ~ log2(G) + WES_purity + CNV
- M1: protein ~ log2(G) + P_principal + WES_purity + CNV

CNV is the gene-level log2 copy-number ratio from `*_WES_CNV_gene_ratio_log2.txt`, matched on the versioned ENSG identifier. If that case is absent from the CNV matrix, the case is dropped for models that include CNV. It is not filled from a neighboring column.

The primary estimand is the partial R² of P_principal, with ΔR² = R²(M1) − R²(M0), the adjusted R² of each model, and the nested comparison of M1 against M0. The reported result is the within-cancer distribution of that partial R² (mean, median, and a case-resampled interval), then a random-effects summary across cancers. Cancers are not stacked into one 1,027-case regression.

## 5. Secondary outcomes

- The same partial R² when the outcome is a gene-level summary of the phosphosite matrix, defined only as the mean of observed sites for that gene in that case. Sites are not tests of isoform identity.
- ΔR² when the RSEM gene value replaces log2(G) as the gene covariate.
- A centered-log-ratio sensitivity on the transcripts of a gene, using a fixed pseudocount of half the smallest positive linear value in that gene, reported only as a sensitivity and not as the primary contrast. The isometric log-ratio is not the primary model: the principal-transcript fraction is the interpretable one-degree-of-freedom contrast.

No secondary outcome is a discovery ranking, a pathway test, or a survival test.

## 6. Covariates

WES_purity is a column in every phenotype table. Across the ten phenotype tables, the fraction of rows with WES_purity equal to NA is 0 for CCRCC and PDAC, 0.009–0.018 for HNSCC, LSCC, LUAD, and GBM, 0.045 for COAD, 0.060 for OV, 0.095 for UCEC, and 0.098 for BRCA. These fractions are a property of the covariate table. They were counted before any protein-versus-usage model. A case whose WES_purity is NA is dropped from M0 and M1. ESTIMATE stromal, immune, and combined scores are present in the same tables and are almost complete; the combined ESTIMATE score is the sensitivity replacement, not the primary purity covariate. The BRCA phenotype header has no column whose name contains “batch” or “plex,” so TMT batch is not available as a covariate from that table. It is not invented. The same screen will be repeated per cancer when the cohort tables are built, and a batch column is included only if one is actually present.

Purity and CNV are retained even though they can correlate with RNA. Variance-inflation factors are reported. A sensitivity fit drops CNV only when a gene’s VIF for CNV exceeds 10. Age, sex, and stage are not in the primary model; stage is entirely empty for GBM, and grade is entirely empty for BRCA, COAD, and GBM.

## 7. Inclusion and exclusion

- Cases: the exact intersection of case IDs in tumor gene RNA, tumor isoform RNA, and tumor protein. Phase 0 targets are the audit reference (BRCA 121, CCRCC 103, COAD 96, GBM 99, HNSCC 108, LSCC 108, LUAD 110, OV 82, PDAC 105, UCEC 95). A discrepancy is reported. Identifiers are not put back to force the total to 1,027.
- Gene key: the exact versioned ENSG string is the primary identifier. A version-stripped match is reported beside it. Version-stripped identifiers are used for a model only when the versioned intersection is smaller and the cohort report shows that the loss is the version suffix. That choice is an identifier rule. It is not a search over genes.
- The match is case-level. It is not described as the same aliquot.
- A gene enters M0/M1 only when at least 71 cases have protein, log2(G), P_principal, WES_purity, and CNV all observed, and the gene has at least two isoform rows.
- OV is not forced through that same gate. Its Phase 0 average-gene effective N was 61.4. OV is analyzed separately. Its genes enter only when the full-matrix count of complete cases reaches 71. If too few OV genes meet that count, OV contributes no primary estimate, and that fact is reported.
- Unmatched IDs from Phase 0 stay unmatched: BRCA `03BR011`, COAD `20CO006`, OV `01OV029`, the ten COAD RNA-only cases, and the three COAD case-list-only IDs.
- Confirmation, DIA, retrospective TCGA proteomes, GTEx, and independent GEO series are not added.

## 8. Missing data

Protein NA, phosphosite NA, and purity NA are missing. They are not imputed in the primary fit. A case with a missing covariate is left out of that gene’s model. Transcript and gene-RNA exact zeros stay zeros. Infinite or non-numeric tokens, if the census finds any, exclude that cell.

## 9. Compositional data

Transcript fractions in a gene sum to 1 whenever G > 0, so they are compositional. The primary model uses one contrast, P_principal, and does not enter every fraction as its own predictor. Zeros, low linear abundance, and genes with many transcripts are handled by the eligibility rule and by reporting the number of transcripts per gene. A CLR fit is sensitivity only.

## 10. Multiple testing

The primary claim is the within-cancer mean partial R² and its cross-cancer summary, not a significant gene list. If a gene-level screen is shown, it uses Benjamini–Hochberg on the nested-model p values within that cancer, and it is labeled as a screen. Bonferroni control at 5 × 10⁻⁶ is not claimed to be powered for a residual partial R² of 0.10. Phase 0’s approximation put that genome-wide target near a residual partial R² of 0.23 at N = 100.

## 11. Power

For one added predictor, the large-sample non-centrality at two-sided α = 0.05 and power 0.80 is λ = (1.959964 + 0.841621)² = 7.848880. Then f² = λ / N and partial R² = f² / (1 + f²). N here is the number of complete cases for that gene, after protein, purity, and CNV missingness, not the column count of the matrix.

M1 has an intercept, log2(G), P_principal, purity, and CNV, so five coefficients. The one-coefficient approximation is the planning number. The cohort report replaces it with the complete-case N from the full matrices once those matrices are on disk. A gene below 71 complete cases is not in the primary distribution.

## 12. Cross-validation

Any predictive check resamples cases within a cancer. Normalization, the choice of which genes meet the 71-case rule, and any gene filter are fit inside the training cases only. The held-out quantity is the squared error of M1 versus M0 on the held-out cases, using coefficients from the training cases. There is no pooling of cancers inside a fold. No independent paired transcript-plus-protein cohort was verified in Phase 0, so an external validation fold is not claimed.

## 13. Cancer-specific analysis

Every model is fit inside one cancer. Sample size, missingness, and the partial-R² distribution are reported per cancer.

## 14. Meta-analysis

The cancer-level mean partial R² values are combined with a random-effects summary. Weights use the within-cancer sampling variance of that mean. A fixed-effect pooled regression of all cases is not the primary analysis.

## 15. Sensitivity analyses

- RSEM gene value in place of log2(G).
- Candidate A versus Candidate B when both inverses remain possible.
- Models that drop CNV when VIF exceeds 10.
- ESTIMATE combined score in place of WES_purity.
- The 483 cases with a single GDC sample labeled Primary Tumor, reported separately and not substituted for the case-level set. That subset leaves UCEC, LUAD, and HNSCC too small for the 71-case rule and is not the primary cohort.

## 16. Negative-result interpretation

A within-cancer mean partial R² whose interval covers values below 0.02, or a cross-cancer summary that does so, means this prespecified contrast does not add a practically meaningful increment beyond gene-level RNA, purity, and CNV in these data. That result does not prove that no transcript feature is informative, and it does not prove that protein isoforms are unchanged. It is the result the primary design can return.

Outlier rule: a case whose protein residual in M0 has absolute Pearson residual above 4 is refit without that case as a sensitivity. The primary fit keeps the case.

Convergence: a gene whose least-squares fit fails, or whose design matrix is rank-deficient, is counted as not estimated and is left out of the mean. It is not assigned a zero partial R².

## What this plan does not do

No differential expression, no differential transcript usage, no isoform-switch ranking, no pathway enrichment, no survival model, and no machine-learning training. Phase 1 starts only if the project owner authorizes it after the Phase 0.5 decision.
