# Phase 0.5 decision

**PARTIAL**

The 1,027-case cohort rebuilds exactly, and the gene-level model can be specified. All 70 required matrices are validated. The isoform values admit two inverses, and no public table supplies a peptide that can be checked against one isoform. Phase 1 has not been started.

## 1. Executive summary

The ten LinkedOmicsKB Release 1 tumor matrices for gene RNA, isoform RNA, and gene-level protein were downloaded from the original May 2023 S3 freeze, checked against Content-Length, and read in full. Their case identifiers intersect at the Phase 0 counts: BRCA 121, CCRCC 103, COAD 96, GBM 99, HNSCC 108, LSCC 108, LUAD 110, OV 82, PDAC 105, UCEC 95. The total is 1,027. No identifier was added to force that total, and the 16 identifiers outside the intersection stay outside it.

Every protein row uses a versioned ENSG identifier that is also a row in the gene-RNA matrix. The isoform matrix has 108,777 rows in every cancer, and every row maps to GENCODE v34. The stored isoform values are non-negative, sit on a 0.01 grid, and include numbers between 0 and 1. An unclipped log2(FPKM) cannot produce that pattern. log2(FPKM+1) and a clipped or rounded log2(FPKM) both can. For the prespecified principal-transcript fraction, those two inverses differ by more than 0.05 in 34–44% of gene-cases where the fraction is defined. The statistical plan therefore fits both. It does not pick one after seeing the protein values.

Every cancer has at least 3,971 genes with 71 or more cases that also have gene RNA, protein, a positive linear transcript sum, a unique principal transcript, WES purity, and CNV. OV has 4,250 such genes. Its median eligible N is 78. Its mean is 61.4, because missingness is uneven across genes. OV is not given a lower gate.

## 2. Acquisition

Source: `https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/`. No newer freeze was substituted. Files were written to `.partial` and renamed only after the byte count matched Content-Length, the header parsed, and a whole-file scan found no NUL bytes. URLs, HTTP status, Content-Length, ETag, and SHA256 are in `phase0_5/manifests/`.

| Layer | Validated tumor matrices |
| --- | ---: |
| Clinical meta | 10 / 10 |
| Phenotype | 10 / 10 |
| Gene RNA, RSEM coding UQ log2 | 10 / 10 |
| Isoform RNA, FPKM log2 | 10 / 10 |
| Global proteome | 10 / 10 |
| Gene-level CNV | 10 / 10 |
| Phosphosite abundance | 10 / 10 |

All ten phosphosite matrices are validated. Earlier transfers of UCEC, LUAD, LSCC, and GBM reached Content-Length while unread pieces were still holes. Those copies contained NUL bytes and were not accepted. The copies now on disk contain no NUL bytes. Sampled byte ranges from those four files match the original S3 objects. The same check matched sampled ranges from BRCA phosphosite, COAD protein, GBM CNV, UCEC isoform RNA, OV phosphosite, and LSCC gene RNA.

SEPEP files are not in this S3 freeze. The six that could be downloaded are stored separately under `phase0_5/acquisition/sepep/`. Their Last-Modified date is 15 April 2026.

## 3. Integrity

The original Phase 0 deliverables were copied to `phase0_5/frozen_phase0/` and locked. They were not edited. The workspace is not a git repository. The freeze record is `phase0_5/frozen_phase0/freeze_record.json`.

AppleDouble files whose names begin with `._` were not read. `.partial` files were not treated as finished matrices. Exact zeros were counted as zeros.

## 4. Reconstructed sample counts

Headers were read in full. Pairing uses the exact case submitter identifier. Row order was not used. The headers contain no aliquot barcode, so the match is case-level.

| Cancer | Recomputed triple | Phase 0 target | Result |
| --- | ---: | ---: | --- |
| BRCA | 121 | 121 | same identifiers |
| CCRCC | 103 | 103 | same identifiers |
| COAD | 96 | 96 | same identifiers |
| GBM | 99 | 99 | same identifiers |
| HNSCC | 108 | 108 | same identifiers |
| LSCC | 108 | 108 | same identifiers |
| LUAD | 110 | 110 | same identifiers |
| OV | 82 | 82 | same identifiers |
| PDAC | 105 | 105 | same identifiers |
| UCEC | 95 | 95 | same identifiers |
| Total | 1,027 | 1,027 | no identifier added or removed |

Sixteen other identifiers appear in at least one matrix and are not in the triple. All sixteen were already outside the Phase 0 triple.

- BRCA `03BR011`: protein, phosphosite, CNV, phenotype, and clinical meta. No RNA.
- OV `01OV029`: protein, phosphosite, phenotype, and clinical meta. No RNA.
- COAD `20CO006`: protein, phosphosite, phenotype, and clinical meta. No RNA.
- COAD RNA without protein: `01CO001`, `05CO005`, `05CO007`, `05CO014`, `05CO055`, `11CO059`, `11CO060`, `14CO002`, `16CO012`, `24CO005`.
- COAD phenotype and clinical meta only: `05CO004`, `13CO001`, `14CO001`.

Protein column counts are higher than the triple where those extra protein cases exist: BRCA 122, COAD 97, OV 83. They were not moved into the analysis set.

## 5. RNA–protein gene intersection

Counts come from the row identifiers, not from file size. In every cancer the gene-RNA matrix has 60,669 rows and the isoform matrix has 108,777 rows.

| Cancer | Protein genes | Shared versioned ENSG | RNA genes with no protein row | Protein genes with no RNA row |
| --- | ---: | ---: | ---: | ---: |
| BRCA | 12,022 | 12,022 | 48,647 | 0 |
| CCRCC | 11,260 | 11,260 | 49,409 | 0 |
| COAD | 9,151 | 9,151 | 51,518 | 0 |
| GBM | 12,100 | 12,100 | 48,569 | 0 |
| HNSCC | 11,561 | 11,561 | 49,108 | 0 |
| LSCC | 12,760 | 12,760 | 47,909 | 0 |
| LUAD | 12,433 | 12,433 | 48,236 | 0 |
| OV | 10,589 | 10,589 | 50,080 | 0 |
| PDAC | 11,404 | 11,404 | 49,265 | 0 |
| UCEC | 11,894 | 11,894 | 48,775 | 0 |

Stripping the version suffix does not recover any additional protein gene. Each gene-RNA matrix also has 45 rows whose identifiers end in `_PAR_`. The protein matrix uses the unsuffixed identifier, and that identifier is already in the versioned intersection.

Every isoform row maps to a GENCODE v34 gene. There are no duplicate transcript identifiers. Of the genes represented, 18,245 have at least two isoform rows. The per-gene complete-observation table has 115,174 cancer-by-gene rows. It is ordered by cancer and gene identifier, not by an outcome.

## 6. Missingness

Gene RNA and isoform RNA have no empty cells in the triple-matched columns. Protein NA fractions are 0.147 (GBM) to 0.234 (COAD). Protein matrices have no exact zeros and no negative values. Isoform matrices have no negative values. Exact zeros are 35–56% of isoform cells, depending on the cancer. Those zeros were left as zeros.

In all ten phosphosite matrices, the NA fraction is 0.514 (GBM) to 0.665 (CCRCC). None of those cells are exact zeros. Nine cancers contain the token `-Inf`. LUAD does not. Those tokens were not recoded. The primary outcome is the global protein matrix. CNV values are log2 ratios, so negative values are part of that measurement. CNV has no NA inside the columns it does contain.

The same four OV cases lack both CNV and WES purity: `01OV007`, `01OV018`, `11OV002`, `17OV013`. Other triple cases missing WES purity, CNV, or both are listed in `cohort_integrity_report.md`. ESTIMATE is numeric for all 1,027 triple-matched cases. No phenotype header contains a batch or plex column.

## 7. SEPEP and peptide evidence

Six tumor SEPEP matrices were downloaded from the public LinkedOmics pages. CCRCC, COAD, GBM, and HNSCC returned HTTP 403 and were not reconstructed. A row is a quantification group named `GENE_SEPEP.<index>_<CODE>`. The file has no peptide sequence and no protein accession. The SEPepQuant code writes classes C1–C5. The files use SI, FD, PD, ND, and MG. Those suffixes were not treated as proof of the C1–C5 mapping.

No row meets Level 3. Level 3 would require a peptide whose sequence is visible and maps uniquely to one isoform. The six matrices can be described only as labeled quantification groups. The global proteome remains one number per gene. Details are in `sepep_feasibility_report.md`.

## 8. Transcript-usage calculation

The isoform values were not summed on the log2 scale, and log2 transcript values were not divided by log2 gene values.

For a stored value y the two inverses that survive the census are:

- Candidate B, inverse of log2(FPKM+1): y = 0 stays 0, and y > 0 becomes 2^y − 1.
- Candidate C, inverse of a clipped or rounded log2(FPKM): y = 0 stays 0, and y > 0 becomes 2^y.

G is the sum of the linear values of the isoform rows for that gene. A gene-case whose observed transcripts are all stored as 0 has G = 0 under both inverses. There are no gene-cases in which every isoform row is missing, because the isoform matrices have no empty cells. G = 0 cases are excluded and are not imputed. The counts are in `gene_protein_intersection.csv`.

The principal transcript is the single GENCODE v34 MANE Select transcript when the gene has one. Otherwise it is the single APPRIS principal-1 transcript. MANE wins when the two tags disagree. 1,001 genes have more than one principal-1 transcript and no MANE transcript; they are not assigned a transcript by file order. P is the principal transcript’s linear value divided by G, and only when G > 0.

Among gene-cases with that fraction defined, the share on which the two inverses differ by more than 0.05 is 0.336 in OV, 0.355 in BRCA, 0.372 in COAD, 0.410 in UCEC, 0.422 in CCRCC, 0.423 in GBM, 0.426 in LSCC, 0.431 in HNSCC and LUAD, and 0.444 in PDAC.

## 9. Statistical analysis plan

The plan is `statistical_analysis_plan.md`. The primary hypothesis is unchanged: a transcript-usage feature from the isoform matrix explains additional variation in gene-level protein abundance beyond a gene-level RNA covariate from that same matrix, tumor purity, and gene-level copy number. It is not a claim about protein function, and it is not a measurement of isoform-specific protein abundance.

M0 is protein ~ log2(G) + WES purity + CNV. M1 adds the principal-transcript fraction. The same cases are used in both. Cancers are fit separately. The reported result is the within-cancer distribution of the partial R², then a random-effects summary. The 1,027 cases are not one regression. Both inverses are fit. The RSEM gene value is a sensitivity covariate, not the primary one.

The one-coefficient large-sample anchor is λ = 7.848880. At the median primary-eligible N, the partial R² detectable at two-sided α = 0.05 and power 0.80 is 0.067–0.091. The historical planning anchor of 0.10 is inside that range for every cancer’s median gene. It is a planning number, not an observed effect. Genes with fewer than 71 complete cases are left out of the primary distribution. That rule is applied to OV as well.

## 10. Cancer-specific eligibility

Primary-eligible means the case has numeric gene RNA, numeric protein, G > 0, a defined principal fraction, numeric WES purity, and numeric CNV, and the gene has at least two isoform rows.

| Cancer | Genes with ≥71 RNA and protein values | Genes with ≥71 primary-eligible cases | Median eligible N | Mean eligible N |
| --- | ---: | ---: | ---: | ---: |
| BRCA | 9,582 | 5,933 | 109 | 90.3 |
| CCRCC | 8,401 | 5,385 | 103 | 83.3 |
| COAD | 6,242 | 3,971 | 95 | 72.8 |
| GBM | 9,801 | 6,189 | 98 | 84.2 |
| HNSCC | 8,959 | 5,678 | 107 | 88.3 |
| LSCC | 10,380 | 6,548 | 107 | 91.6 |
| LUAD | 9,902 | 6,239 | 108 | 90.2 |
| OV | 6,924 | 4,250 | 78 | 61.4 |
| PDAC | 8,236 | 5,260 | 105 | 82.6 |
| UCEC | 8,998 | 5,512 | 86 | 71.5 |

The Phase 0 sample estimate of 61.4 for OV was an approximation of protein completeness. On the full matrix, the mean number of cases with both RNA and protein numeric is 64.9. The mean of the stricter primary-eligible count is 61.4. Those are different quantities. The median primary-eligible N in OV is 78, and 4,250 genes meet the 71-case rule.

## 11. Remaining limitations

- The isoform inverse is not unique, and the choice changes the prespecified fraction in roughly one third to one half of the gene-cases where it is defined.
- No Level 3 isoform-specific peptide evidence is in hand. Four SEPEP files returned HTTP 403.
- The match is case-level. It is not an aliquot-level or same-tissue-block match.
- TMT batch is not a column in the phenotype tables.
- Gene RNA is RSEM. Isoform RNA is a different FPKM quantification. The primary gene covariate is the sum of the inverted isoform rows, not the RSEM value.
- 1,001 genes have no unique principal transcript under the prespecified rule.

## 12. What this phase did not do

No differential expression, differential transcript usage, isoform-switch list, gene ranking, pathway enrichment, survival model, or machine-learning fit was run. TCGA and GTEx were not paired. Confirmation and DIA cohorts were not added. The 1,027-case set was not pooled into one regression. Phase 1 starts only if the project owner authorizes it.

## Decision against the stated requirements

| Requirement | Status |
| --- | --- |
| Validated matrices for the primary RNA and protein comparison | Met for gene RNA, isoform RNA, and protein in all ten cancers |
| All required matrices, including phosphosite | Met: 70 / 70, each with a NUL-free whole-file scan |
| Cohort reconstructed from complete headers | Met: 1,027, with the same identifiers |
| Exact RNA–protein gene intersection | Met: every protein gene is in the RNA matrix |
| One determined transcript-usage transform | Not met: two inverses remain possible and they disagree materially |
| Prespecified model | Met |
| Effective sample size under that model | Met for thousands of genes in every cancer, including OV, at the 71-case gate |
| No fatal barrier to the gene-level comparison | Met |

PARTIAL is the applicable decision: the gene-level comparison is feasible, and the transcript-usage transform remains a material qualification.
