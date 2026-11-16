# SEPEP final assessment

Audit date: 2026-11-14. This review uses the six files already stored under `phase0_5/acquisition/sepep/` and the Phase 0.5 retrieval record for the four HTTP 403 responses. No new request was sent to the restricted URLs. No isoform-specific protein abundance is claimed.

## What the open files are

Each downloaded object is a sample-by-feature matrix. The header token is `idx`, and the remaining columns are case identifiers. A fresh read of the first data row in every file returned the identifier `DPM1_SEPEP.1_ND`. That identifier is a gene symbol, a SEPEP index, and a suffix. It is not a peptide sequence, and it is not a protein accession.

| Cancer | Bytes on disk | Header columns, including `idx` | HTTP at retrieval |
| --- | ---: | ---: | --- |
| BRCA | 38,342,442 | 123 | 200 |
| LSCC | 35,103,944 | 109 | 200 |
| LUAD | 33,543,053 | 111 | 200 |
| OV | 21,284,156 | 84 | 200 |
| PDAC | 26,312,593 | 106 | 200 |
| UCEC | 28,506,355 | 96 | 200 |
| CCRCC | — | — | 403 |
| COAD | — | — | 403 |
| GBM | — | — | 403 |
| HNSCC | — | — | 403 |

The column counts match the Phase 0.5 sample counts plus the `idx` column: BRCA 122, LSCC 108, LUAD 110, OV 83, PDAC 105, UCEC 95. These are patient-level measurements of a named quantification group. They are not a peptide-sequence table.

The LinkedOmics page calls the unit “TMT proteomics peptide (log2 ratio)” and the description “Proteome Peptide; Tumor only; log2 ratio.” The filename is `SEPEP_log2_ratio_normalized`. The file body does not carry the peptide sequence or the protein-accession list that would let a row be checked against one isoform.

## Mapping and isoform specificity

SEPepQuant defines a SEPEP as a group of peptides that connect to the same protein set, and its code writes classes C1–C5. The matrices on disk use the suffixes SI, FD, PD, ND, and MG. Those suffixes are not emitted by the retrieved code. Their initials resemble the paper’s English class names. The matrix does not include the mapping table, so the resemblance is not a demonstrated crosswalk.

A Level 3 call would require a peptide sequence that maps uniquely to one isoform, with identification quality still open. None of the six files contains that sequence or that mapping. CCRCC, COAD, GBM, and HNSCC were not retrieved. No public accession in the Phase 0.5 PDC query supplied a peptide file with usable file metadata. Controlled-access searches were not attempted.

## Consequence

The six open matrices can be described as patient-level quantifications of labeled SEPEP groups. They do not identify a unique isoform-specific peptide. The gene-level protein matrix remains the outcome for the primary model.
