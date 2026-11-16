# Transcript-usage sensitivity

Audit date: 2026-11-14. The calculation is in `phase0_6/code/quantification_audit.py`. Cell-level results are in `phase0_6/manifests/transformation_comparison.csv`. No protein model, differential test, isoform-switch ranking, or gene ranking was fit.

The isoform matrix has the same 108,777 transcript rows in every cancer after the BRCA tail repair described below. All 108,777 map to GENCODE v34. There are 60,669 genes, 42,424 with one matrix row and 18,245 with two or more. Duplicate transcript identifiers: 0. Values off the 0.01 grid: 0.

## Annotation rule

The principal transcript is fixed from GENCODE v34 before any abundance model. A gene with exactly one MANE Select transcript uses that transcript. A gene with no MANE Select transcript uses APPRIS principal-1 only when exactly one transcript carries that tag. MANE Select is kept when the two tags disagree. The annotation contains 14,620 genes with one MANE Select transcript, 2,819 with no MANE and exactly one APPRIS principal-1 transcript, and 1,001 genes with no unique principal. Those 1,001 genes are left unassigned.

Of the 17,439 genes with a unique principal, 17,402 have that transcript as a row in the isoform matrix and 37 do not. Of those present, 10,795 genes also have at least two isoform rows. The principal-proportion quantiles below use that gene set. `principal_quantile_observations` in the comparison table is 10,795 times the sample count for literal A, and the count of those gene-samples with a positive literal-B sum for literal B.

## Samples in this audit

The audit uses the intersection of case identifiers in the isoform matrix and the gene-RNA matrix. That intersection equals the Phase 0 protein triple except in COAD, where the RNA intersection has 106 cases and the protein triple has 96. The ten extra COAD identifiers are the RNA-without-protein cases already listed in Phase 0.5. They belong in a reconstruction audit. They do not enter the protein model.

BRCA isoform RNA and UCEC gene RNA were read from `phase0_6/repaired/`. The copies under `phase0_5/acquisition/` still end in a NUL span, and those files were not edited. The repaired BRCA isoform file has 108,777 data rows. The Phase 0.5 missingness row for BRCA isoform reports 108,733 features; that count was taken from the NUL-tailed copy.

## Literal inverses

Literal A sets X = 2^Y, so a stored 0 contributes 1. The gene sum is at least the number of isoform rows, and the usage proportion is defined for every gene-sample.

Literal B sets X = 2^Y − 1, so a stored 0 contributes 0. The gene sum is 0 when every isoform row is 0. Those gene-samples have no usage proportion. Their share of all gene-samples is 0.352 (GBM) to 0.561 (OV).

| Cancer | Samples | Gene-samples | All stored values 0 | One positive row under B | Two or more positive rows under B |
| --- | ---: | ---: | ---: | ---: | ---: |
| BRCA | 121 | 7,340,949 | 3,963,802 | 2,015,765 | 1,361,382 |
| CCRCC | 103 | 6,248,907 | 2,415,990 | 2,473,856 | 1,359,061 |
| COAD | 106 | 6,430,914 | 3,525,914 | 1,702,298 | 1,202,702 |
| GBM | 99 | 6,006,231 | 2,111,413 | 2,529,698 | 1,365,120 |
| HNSCC | 108 | 6,552,252 | 2,466,434 | 2,639,856 | 1,445,962 |
| LSCC | 108 | 6,552,252 | 2,433,374 | 2,652,755 | 1,466,123 |
| LUAD | 110 | 6,673,590 | 2,501,456 | 2,691,391 | 1,480,743 |
| OV | 82 | 4,974,858 | 2,789,700 | 1,288,951 | 896,207 |
| PDAC | 105 | 6,370,245 | 2,390,823 | 2,561,641 | 1,417,781 |
| UCEC | 95 | 5,763,555 | 2,208,679 | 2,294,826 | 1,260,050 |

Under A, every transcript row is positive after the inverse, so “one detected transcript” coincides with the 42,424 single-row genes, in every sample. Under B, the genes with exactly one positive row in every sample number 3,989 (OV) to 11,478 (GBM).

## Principal proportion

For a gene-sample with a unique principal and at least two isoform rows, P is the linear value of that transcript divided by the gene sum. Under A this includes samples whose stored values are all 0: each transcript contributes 1, and P equals 1 over the number of rows. Under B those samples are omitted.

| Cancer | A observations | A p05 | A p50 | A p95 | B observations | B p05 | B p50 | B p95 | Mean \|ΔP\| | Share \|ΔP\| > 0.05 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BRCA | 1,306,195 | 0.084 | 0.500 | 0.939 | 1,208,877 | 0 | 0.670 | 1 | 0.145 | 0.631 |
| CCRCC | 1,111,885 | 0.105 | 0.492 | 0.911 | 1,055,903 | 0 | 0.615 | 1 | 0.138 | 0.629 |
| COAD | 1,144,270 | 0.085 | 0.496 | 0.938 | 1,058,431 | 0 | 0.656 | 1 | 0.146 | 0.630 |
| GBM | 1,068,705 | 0.104 | 0.486 | 0.902 | 1,028,315 | 0 | 0.594 | 1 | 0.132 | 0.614 |
| HNSCC | 1,165,860 | 0.099 | 0.482 | 0.912 | 1,117,966 | 0 | 0.600 | 1 | 0.139 | 0.629 |
| LSCC | 1,165,860 | 0.098 | 0.484 | 0.911 | 1,122,236 | 0 | 0.597 | 1 | 0.134 | 0.616 |
| LUAD | 1,187,450 | 0.102 | 0.492 | 0.911 | 1,139,701 | 0 | 0.614 | 1 | 0.136 | 0.627 |
| OV | 885,190 | 0.081 | 0.498 | 0.942 | 815,409 | 0 | 0.659 | 1 | 0.146 | 0.629 |
| PDAC | 1,133,475 | 0.103 | 0.493 | 0.912 | 1,090,715 | 0 | 0.622 | 1 | 0.138 | 0.639 |
| UCEC | 1,025,525 | 0.100 | 0.491 | 0.911 | 979,310 | 0 | 0.612 | 1 | 0.135 | 0.614 |

The absolute difference is computed only where B is defined. Across cancers, the mean absolute difference is 0.132–0.146, and the share exceeding 0.05 is 0.614–0.639. That is larger than the Phase 0.5 B-versus-C range of 0.336–0.444 because literal A assigns mass 1 to every stored zero, while Phase 0.5 Candidate C assigns those zeros mass 0.

The median of the maximum transcript proportion is 1 under both inverses. Most genes have a single isoform row, and a single row has proportion 1 whenever the sum is positive. This median is a property of the row counts.

## Dominant transcript and entropy

Both inverses are strictly increasing, so they share the unique maximum of Y. The audit found 0 disagreements in every cancer. Gene-samples with two transcripts tied at the maximum stay without a dominant identity. Mean entropy under A, over every gene-sample, is 0.417–0.432 bits. That average includes single-row genes, whose entropy is 0, and all-zero multi-row genes, whose A proportions are equal and whose entropy is the log2 of the row count. Mean entropy under B, over gene-samples with a positive sum, is 0.376–0.402 bits. The two means use different denominators.

## Representation chosen for the model

Five representations were considered. The choice uses the form of the estimand. It does not use an association with protein, because no protein association was calculated.

| Representation | Role |
| --- | --- |
| Raw principal proportion | Primary. One prespecified number, P = X_principal / G. |
| Principal-versus-rest contrast | The same primary contrast. The rest is 1 − P. |
| MANE / APPRIS rule | The rule that names the principal transcript. It is part of the primary contrast, not a second model. |
| Centered log-ratio of the gene’s transcripts | Sensitivity only. The pseudocount is half the smallest positive linear value in that gene, and that value depends on the inverse. |
| Isometric log-ratio | Not the primary model. It replaces the one-degree-of-freedom contrast with a full composition. |

The primary representation is therefore defined. Its numerical value is a pair, one under each literal inverse, until the transform is resolved.
