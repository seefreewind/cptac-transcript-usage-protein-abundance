# SEPEP and isoform-specific peptide evidence

Audit date: 2026-11-13. This is an evidence audit. No isoform abundance and no protein-functional-loss claim is made.

## What was retrieved

The LinkedOmicsKB v1.2 S3 URL list does not contain a SEPEP object. The public LinkedOmics pancan pages do. Each page links a tumor file `CPTAC_pancan_<CANCER>_proteomics_SEPEP_log2_ratio_normalized_Tumor.cct`, described there as “Proteome Peptide; Tumor only; log2 ratio.”

Those objects are not the 10 May 2023 S3 freeze. The six files that responded have `Last-Modified: Wed, 15 Apr 2026`. They are stored under `phase0_5/acquisition/sepep/`, separate from the v1.2 matrices.

| Cancer | HTTP | Bytes | Rows | Genes in the row label | Samples | Phase 0 paired cases present | Extra samples |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| BRCA | 200 | 38,342,442 | 32,616 | 12,572 | 122 | 121 / 121 | `03BR011` |
| LSCC | 200 | 35,103,944 | 35,107 | 13,453 | 108 | 108 / 108 | none |
| LUAD | 200 | 33,543,053 | 33,943 | 13,221 | 110 | 110 / 110 | none |
| OV | 200 | 21,284,156 | 27,694 | 11,530 | 83 | 82 / 82 | `01OV029` |
| PDAC | 200 | 26,312,593 | 28,894 | 12,161 | 105 | 105 / 105 | none |
| UCEC | 200 | 28,506,355 | 32,592 | 12,869 | 95 | 95 / 95 | none |
| CCRCC | 403 | — | — | — | — | not retrieved | — |
| COAD | 403 | — | — | — | — | not retrieved | — |
| GBM | 403 | — | — | — | — | not retrieved | — |
| HNSCC | 403 | — | — | — | — | not retrieved | — |

A second request, including the public page as referrer, still returned 403 for CCRCC, COAD, GBM, and HNSCC. Those four files are unavailable. They were not reconstructed.

The 621 Phase 0 case-matched tumors in the six open cohorts are all present as column IDs. The two extra columns are the protein-only cases already excluded in Phase 0. Row order was not used.

Across the six files, about 0.30–0.37 of cells are NA. Exact zeros do not occur. About 0.31–0.35 of cells are negative, which is expected for a log2 ratio. Invalid tokens: 0.

## What a row is

A row identifier has the form `GENE_SEPEP.<index>_<CODE>`, for example `DPM1_SEPEP.1_ND` and `SCYL3_SEPEP.1_FD`. The identifier does not contain a peptide sequence, a protein accession, or a spectrum. Counts of the five suffixes:

| Cancer | SI | FD | PD | ND | MG |
| --- | ---: | ---: | ---: | ---: | ---: |
| BRCA | 3,068 | 3,676 | 10,171 | 7,958 | 7,743 |
| LSCC | 3,350 | 3,980 | 10,947 | 8,531 | 8,299 |
| LUAD | 3,243 | 3,873 | 10,676 | 8,376 | 7,775 |
| OV | 2,751 | 2,846 | 8,507 | 7,205 | 6,385 |
| PDAC | 2,976 | 3,162 | 8,911 | 7,566 | 6,279 |
| UCEC | 3,146 | 3,667 | 10,172 | 8,143 | 7,464 |

SEPepQuant (Dou et al., Nature Communications, 2023; code retrieved 2026-11-13 from `bzhanglab/SEPepQuant`, `scritps/creat-mapping-table-TMT.r`) defines a SEPEP as a set of peptides that connect to the same proteins, and it writes five classes:

- C1: one gene, and the peptides cover that gene’s only protein entry in the database used.
- C2: one gene with more than one protein, and the peptides map to exactly one of those proteins.
- C3: one gene, and the peptides map to more than one but not all of its proteins.
- C4: one gene, and the peptides map to every protein of that gene.
- C5: peptides map to more than one gene.

The code emits `C1`–`C5`. It does not emit `SI`, `FD`, `PD`, `ND`, or `MG`. The initials match the English names in that paper (single isoform, fully discriminative, partially discriminative, non-discriminative, multi-gene), but the LinkedOmics file does not include the mapping table that would show the peptide sequences and protein accessions. This audit therefore does not treat the suffix as a verified crosswalk to C1–C5.

## Evidence level

Level 3 requires a peptide that we can see maps uniquely to one isoform, with identification quality still open to the usual protein-inference limits. None of these rows meets that bar: the sequence and the target protein are not in the file.

| Suffix | How it is recorded here |
| --- | --- |
| SI, FD | Suffix present. Not Level 3. `isoform_peptide_evidence.csv` records these rows as `not_level_3`. |
| PD, MG | Suffix present. If the initials mean what the SEPepQuant prose means, the mapping is ambiguous. The CSV records `LEVEL_2_candidate_label_only`. That is not a verified Level 2 call. |
| ND | Suffix present. If the initials mean what the SEPepQuant prose means, the peptide set is shared across isoforms of one gene. The CSV records `LEVEL_1_candidate_label_only`. That is not a verified Level 1 call. |
| No file (CCRCC, COAD, GBM, HNSCC) | Level 0 for this audit. |

The global proteome matrices remain one number per gene. Phosphosite identifiers carry one representative ENSP and a peptide string, which places the site and does not quantify an isoform.

## PDC

The Phase 0 catalog of 277 PDC UI studies has analytical fractions Proteome, Phosphoproteome, Acetylome, Glycoproteome, and smaller metabolite fractions. It has no peptide fraction. A GraphQL `filesPerStudy` query for PDC000120 returned 127 records whose file name, file type, and data category were null, so no peptide file was retrieved from PDC. Controlled-access PSMs and BAMs were not requested.

## Consequence for the study

Isoform-specific protein abundance cannot be measured from the open tables in hand. The six SEPEP matrices can support a later, separately authorized check only after the mapping tables (peptide sequence, protein accessions, and FDR) are obtained from an open source. Their absence does not block the narrower question of whether transcript usage adds information about gene-level protein abundance.
