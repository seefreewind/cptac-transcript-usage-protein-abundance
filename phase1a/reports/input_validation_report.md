# Input validation

Validated at 2026-11-14T10:11:55Z.

Status: PASS. Every used file matches its required hash and contains no NUL bytes.

Repaired copies are used for BRCA isoform RNA, UCEC gene RNA, and PDAC CNV. The Phase 0.5 copies of those three files were not read.
GENCODE v34 is a gzip container. NUL bytes inside that container are compressed-stream bytes. The file is accepted when its SHA256 equals the frozen GENCODE v34 hash. Text matrices are still rejected if they contain NUL bytes.
