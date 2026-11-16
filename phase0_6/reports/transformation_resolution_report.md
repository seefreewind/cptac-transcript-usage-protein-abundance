# Transformation resolution

Audit date: 2026-11-14. The stored isoform matrix is `*_RNAseq_isoform_FPKM_log2_Tumor.txt` from LinkedOmicsKB Release 1, `data_freeze_v1.2_reorganized`. This report separates what the source documents say from what the stored numbers can support. The two literal inverses were both calculated. Neither was selected from a correlation with RSEM or with protein.

## What the source documents say

The LinkedOmics portal record in `phase0/raw/api/linkedomics_portal_counts.json` gives the isoform tumor file two fields together:

- description: “Linear mRNA isoforms, upper quartile normalized”
- unit: “Expression (log2 transformed FPKM)”

The filename contains `FPKM_log2`. The same portal record gives the gene-RNA file the unit “Expression (log2 transformed RSEM)”. The project manifest describes that gene matrix as RSEM v1.3.1, with the upper quartile of coding genes scaled to 1,500, then log2.

A different object in the same bucket is named as a log2(x+1) transform. The saved URL list `phase0/raw/api/linkedomicskb_s3_urls.txt` and the saved download page `phase0/raw/html/linkedomicskb_download.html` point to `CPTAC_pan_cancer_RNASeq_TPM_for_Xinpei/*_log2(x+1)_BCM.txt`. Those are gene-level TPM files. They are a different measurement from the isoform FPKM matrix. When that release writes `log2(x+1)` into an object name, the name carries the pseudocount. The isoform filename and the isoform unit string do not.

No processing script for these isoform files was found in the project record, and the saved LinkedOmics manual extract does not state a pseudocount or a floor for this FPKM matrix.

## What the stored values add

The Phase 0.5 missingness table, on the triple-matched columns, records the same pattern in every cancer: no negative values, minimum 0, minimum positive value 0.01, every numeric value on a 0.01 grid, and a point mass at 0. The share of numeric cells in (0, 1) is 0.220 (OV) to 0.345 (HNSCC). Exact zeros are 0.353 (GBM) to 0.555 (OV) of isoform cells.

An unclipped log2 of an FPKM below 1 is negative. The matrices contain values in (0, 1) and contain no negatives, so the stored column is not an unclipped log2(FPKM).

Two written conventions are compatible with that pattern:

- rounded log2(FPKM+1), which is nonnegative for every nonnegative FPKM and equals 0 only at FPKM 0, up to rounding
- log2(FPKM) after a floor at 0 and rounding to two decimals, which puts every FPKM at or below 1 into the stored zero

The documents name the first convention on a different file. They name the second convention’s family, “log2 transformed FPKM”, on this file, and they do not write the floor.

## The two literal inverses

Phase 0.6 uses these names. They are the inverses of the stored number, applied to every cell, including stored zeros.

| Name | Formula | Stored 0 becomes |
| --- | --- | --- |
| Literal A | X = 2^Y | 1 |
| Literal B | X = 2^Y − 1 | 0 |

Literal A is the ordinary inverse of a stored log2 value. On the zero mass it returns 1, which is the top of the interval of FPKM values that a floor at log2 = 0 would have collapsed. Literal B is the ordinary inverse of log2(x+1). On the zero mass it returns 0. Rounding means either formula recovers the center of a small bin rather than a uniquely known preimage, and the zero bin is wider than the 0.01 step in the positive range.

Both functions increase with Y. The transcript with the unique largest stored value is therefore the same transcript under both inverses. Proportions and entropies change, because a stored zero contributes 1 under A and contributes nothing under B.

## How these names sit next to Phase 0.5

Phase 0.5 used a different pair, and both members of that pair keep a stored zero at zero:

- Phase 0.5 Candidate B: Y = 0 stays 0, and Y > 0 becomes 2^Y − 1. This is the same function as literal B.
- Phase 0.5 Candidate C: Y = 0 stays 0, and Y > 0 becomes 2^Y.

Candidate C is a third arithmetic reading. It is the clipped reading of log2(FPKM). It is not literal A. Phase 0.5 reported that B and C differ by more than 0.05 in 0.336–0.444 of principal gene-cases where both fractions are defined. That figure describes B versus C. The literal A versus literal B comparison in this phase is a different pair, and its fractions are in `transcript_usage_sensitivity_report.md`.

## Resolution

The portal unit and the filename support “log2 transformed FPKM”. The separate TPM objects show that this release writes `log2(x+1)` when that is the transform it is naming. The stored pattern rules out unclipped log2(FPKM) and leaves both a floored log2(FPKM) and a rounded log2(FPKM+1) possible. No document in the project states the pseudocount or the floor for this matrix.

The primary transcript-usage calculation therefore keeps both literal inverses. A correlation with RSEM is a scale audit of the reconstructed sum. It is not a rule for choosing the inverse. Phase 1 has not been started.
