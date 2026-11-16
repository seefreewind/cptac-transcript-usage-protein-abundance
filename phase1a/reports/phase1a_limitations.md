# Phase 1A limitations

- RSS clipping events: 0
- Case matching is case-level.
- S4 CLR was not run: NOT_RUN_UNDERDEFINED.
- CLR sensitivity S4 was not run because the frozen documents define the pseudocount and do not define the one-degree-of-freedom contrast.
- Phase 1B has not been authorized.
- GENCODE v34 is stored as gzip. The 163,390 NUL bytes in that container are compressed-stream bytes. The file was accepted because its SHA256 matches the frozen GENCODE v34 hash. Text matrices with NUL bytes remain a hold.
- Some triple cases are absent from the CNV header. Those cases stay in the cohort and drop out of a gene fit only when CNV is not numeric. BRCA isoform, UCEC gene RNA and PDAC CNV were read from the Phase 0.6 repaired copies. On those repaired copies the Literal B pre-fit counts still matched the published table, so no historical recount was required.
- Isoform row order differs by cancer. Each matrix was required to contain the same 108,777 transcript identifiers as the repaired BRCA matrix.
- The cancer-level estimand is the unweighted mean. Medians are lower than means because the gene-level partial-R² distribution is right-skewed.
- Candidate C is a sensitivity. It is not a third co-primary transformation. S4 CLR was left undefined by the frozen contrast and was not run.
