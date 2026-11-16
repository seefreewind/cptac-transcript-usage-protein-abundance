# Cohort integrity

Recomputed from complete matrix headers and full matrices in `phase0_5/acquisition`. Phase 0 files were not edited.

GENCODE v34 transcript rows indexed: 228048. Genes with exactly one MANE Select transcript: 14620. Genes with no MANE transcript and exactly one APPRIS principal-1 transcript: 2819. Genes with no unique principal transcript: 1001. A gene in that last set is not assigned a transcript by file order.

The two linear inverses are applied only to describe denominators and the stability of the principal-transcript fraction. They are not used to rank genes or to fit a protein model.

Candidate B, the inverse of log2(FPKM+1): stored 0 stays 0, and any other finite y becomes 2^y − 1.

Candidate C, the inverse of a clipped or rounded log2(FPKM): stored 0 stays 0, and y > 0 becomes 2^y. A negative y blocks Candidate B for that gene-case and is inverted as 2^y under Candidate C.

A gene-case has a zero linear denominator when every observed transcript is stored as 0. Those cells were not recoded as NA.

- BRCA: recomputed triple 121 matches the Phase 0 target. No ID was added or removed.
- BRCA: isoform numeric min 0.0, min positive 0.01, max 16.87, negative cells 0, exact zeros 7032222, values in (0, 1) 3091002 of 13156681 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- BRCA: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- BRCA: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- BRCA: triple cases absent from the CNV header: 01BR009, 01BR042, 01BR043, 01BR044, 03BR006, 05BR004, 09BR004, 11BR027, 604. They remain in the case-level match and are dropped only from models that require CNV.
- BRCA: triple cases with WES_purity missing: 01BR009, 01BR025, 01BR042, 01BR043, 01BR044, 03BR006, 05BR004, 09BR004, 11BR012, 11BR027, 11BR054, 604. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- BRCA: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- CCRCC: recomputed triple 103 matches the Phase 0 target. No ID was added or removed.
- CCRCC: isoform numeric min 0.0, min positive 0.01, max 19.94, negative cells 0, exact zeros 4380621, values in (0, 1) 3628540 of 11204031 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- CCRCC: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- CCRCC: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- CCRCC: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- COAD: recomputed triple 96 matches the Phase 0 target. No ID was added or removed.
- COAD: isoform numeric min 0.0, min positive 0.01, max 16.94, negative cells 0, exact zeros 5555675, values in (0, 1) 2563926 of 10442592 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- COAD: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- COAD: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- COAD: triple cases absent from the CNV header: 05CO041. They remain in the case-level match and are dropped only from models that require CNV.
- COAD: triple cases with WES_purity missing: 05CO041. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- COAD: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- GBM: recomputed triple 99 matches the Phase 0 target. No ID was added or removed.
- GBM: isoform numeric min 0.0, min positive 0.01, max 19.63, negative cells 0, exact zeros 3803243, values in (0, 1) 3674535 of 10768923 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- GBM: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- GBM: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- GBM: triple cases with WES_purity missing: C3L-02955. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- GBM: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- HNSCC: recomputed triple 108 matches the Phase 0 target. No ID was added or removed.
- HNSCC: isoform numeric min 0.0, min positive 0.01, max 19.84, negative cells 0, exact zeros 4490506, values in (0, 1) 4052589 of 11747916 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- HNSCC: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- HNSCC: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- HNSCC: triple cases with WES_purity missing: C3N-02730. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- HNSCC: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- LSCC: recomputed triple 108 matches the Phase 0 target. No ID was added or removed.
- LSCC: isoform numeric min 0.0, min positive 0.01, max 19.85, negative cells 0, exact zeros 4391764, values in (0, 1) 3994178 of 11747916 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- LSCC: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- LSCC: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- LSCC: triple cases with WES_purity missing: C3N-03848. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- LSCC: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- LUAD: recomputed triple 110 matches the Phase 0 target. No ID was added or removed.
- LUAD: isoform numeric min 0.0, min positive 0.01, max 19.89, negative cells 0, exact zeros 4537602, values in (0, 1) 3998944 of 11965470 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- LUAD: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- LUAD: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- LUAD: triple cases absent from the CNV header: C3L-01890. They remain in the case-level match and are dropped only from models that require CNV.
- LUAD: triple cases with WES_purity missing: C3L-01683, C3L-01890. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- OV: recomputed triple 82 matches the Phase 0 target. No ID was added or removed.
- OV: isoform numeric min 0.0, min positive 0.01, max 16.29, negative cells 0, exact zeros 4951816, values in (0, 1) 1958149 of 8919714 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- OV: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- OV: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- OV: triple cases absent from the CNV header: 01OV007, 01OV018, 11OV002, 17OV013. They remain in the case-level match and are dropped only from models that require CNV.
- OV: triple cases with WES_purity missing: 01OV007, 01OV018, 11OV002, 17OV013. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- OV: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- PDAC: recomputed triple 105 matches the Phase 0 target. No ID was added or removed.
- PDAC: isoform numeric min 0.0, min positive 0.01, max 19.71, negative cells 0, exact zeros 4330743, values in (0, 1) 3915869 of 11421585 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- PDAC: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- PDAC: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- PDAC: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.
- UCEC: recomputed triple 95 matches the Phase 0 target. No ID was added or removed.
- UCEC: isoform numeric min 0.0, min positive 0.01, max 19.27, negative cells 0, exact zeros 4022386, values in (0, 1) 3262090 of 10333815 numeric cells, fraction on a 0.01 grid 1.0, duplicate transcript ids 0.
- UCEC: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, so the zero-denominator count does not choose the inverse.
- UCEC: the gene-RNA matrix has 45 additional rows whose identifiers end in _PAR_. The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. Examples: ENSG00000002586.20_PAR_Y, ENSG00000124333.16_PAR_Y, ENSG00000124334.17_PAR_Y, ENSG00000167393.18_PAR_Y.
- UCEC: triple cases with WES_purity missing: C3L-00008, C3L-00769, C3N-00324, C3N-00334, C3N-00383, C3N-00847, C3N-00866, C3N-01003, C3N-01217. They are dropped from M0 and M1. ESTIMATE is recorded separately.
- UCEC: phosphosite matrix contains non-finite tokens -Inf. Those cells are not recoded. The primary outcome is the global protein matrix.

Pairing remains case-level. Matrix headers contain no aliquot barcode.
Per-gene complete-observation rows written: 115174. Rows are ordered by cancer and gene identifier. They are not ordered by an outcome.

The detectable partial R² column uses the one-coefficient large-sample anchor λ = 7.848880 at the median primary-eligible N. It is a planning number, not an observed effect. OV is not given a lower eligibility gate.
