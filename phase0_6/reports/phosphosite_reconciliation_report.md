# Phosphosite reconciliation

Audit date: 2026-11-14. Phosphosite matrices are supplementary. The primary model uses gene-level protein abundance.

## Where the two counts came from

The count of six validated phosphosite matrices was an intermediate Phase 0.5 state. At that point BRCA, CCRCC, COAD, HNSCC, OV, and PDAC had passed the size and format checks. A UCEC file whose byte count matched Content-Length contained 6,546,742 NUL bytes and was deleted. LUAD, LSCC, and GBM had not yet been accepted.

The Phase 0.5 final report, written after those four files were retrieved again, states 10 of 10. This reconciliation uses the files now on disk, the acquisition manifest, the Phase 0.5 SHA256 ledger, and a new whole-file NUL scan. The deleted UCEC copy is not counted.

## Authoritative table

Scanned on 2026-11-14. Every current file starts with the header token `idx`, matches the acquisition-manifest byte count, matches the manifest SHA256, and contains zero NUL bytes. `validated` is the value in `phase0_5/manifests/acquisition_manifest.csv`.

| Cancer | Bytes | NUL bytes | Manifest validated | SHA256 |
| --- | ---: | ---: | --- | --- |
| BRCA | 50,633,530 | 0 | yes | `64b2f52e9148747a0380b088892cd5d31cce6ac044c97e8bef88d599e084e190` |
| CCRCC | 46,018,835 | 0 | yes | `e373b43f03df3ddb5c648f5af4b57e4caa61e767c911dfd5aa1b3ded1b772d98` |
| COAD | 30,662,724 | 0 | yes | `dfa526d0ec65755195269148edd5c4f87079a32645ef8cdb78fac7cd07167c4b` |
| GBM | 65,019,693 | 0 | yes | `6dfbe3a23e70fabb428b81b30ec3493a96cbb09ca7cf9c48e0d093f777f6c467` |
| HNSCC | 54,793,177 | 0 | yes | `c793e056d97bb39dd11e0490f3ff44d8af9b2b5e368c609ca99912b38603e861` |
| LSCC | 62,042,879 | 0 | yes | `e14fe8506376ed8f755b6066873756971013d6371092d1b7e9e7644f1e370077` |
| LUAD | 57,949,691 | 0 | yes | `5dbdb71093de80b5f8c4e8893c754ea7b3db065c461d5f10e424aa70559d1868` |
| OV | 29,097,272 | 0 | yes | `509b6daee400feff78896eda5b9ce3ddba13ffb45e94f4165711834b0a0a943e` |
| PDAC | 48,706,052 | 0 | yes | `be6834314f24daaebab9f56ad5736e314af87d7868115963e25b9e82aaf1a89f` |
| UCEC | 57,212,205 | 0 | yes | `0cb6fe7b534b00018cffbc77a64d014f72571baa4b49ad4fc4009cb060c0c70f` |

The authoritative count is 10 of 10. The current UCEC file is the later copy. It is not the rejected copy.

## What the matrices contain

On the Phase 0 triple-matched columns, the Phase 0.5 missingness table records NA fractions from 0.514 (GBM) to 0.665 (CCRCC). Exact zeros are absent. Nine cancers contain the token `-Inf`. LUAD does not. Those tokens were left as stored. The phosphosite identifier carries one representative protein accession and a peptide string for the site. That string places the site. It does not quantify an isoform.

## Separate tail damage in three non-phosphosite files

A whole-file scan of all 70 acquired matrices on 2026-11-14 found NUL bytes in three files that are not phosphosite matrices:

| File | NUL bytes, all in one trailing span | First download SHA256 | SHA256 after the 2026-11-13 11:48 rehash |
| --- | ---: | --- | --- |
| BRCA isoform RNA | 16,671 | `c9e071a8b393e75832673b3fd4d3a1754c4a7d25ac2ac0104c2e03f832f43ca9` | `a69e669488a79c8e2eabf6fcd8b37f3d2771f8387a447bbc84f89ce3bb43d7e2` |
| UCEC gene RNA | 317 | `3d708a346477a51f2aa1fc80479cd60b1888644145b69b1082c31c6405770882` | `473b9f48273abcbf3c693ecd8ba0985d0868355cef284a4fb6d1f8eeb29dd9b0` |
| PDAC CNV | 40,364 | `12f2a560765c9a984251afc60201637eefbe063e09e77b8ac37427ef0da67510` | `024ca58a6bc45c7982bcfa812e1a3835c390aeef86f7e294d24360890a737f64` |

The same byte ranges on the original S3 objects are text and contain no NUL bytes. Replacing only those spans produces files whose SHA256 equals the first completed-download hash. Those repaired copies are in `phase0_6/repaired/`. The Phase 0.5 files were not edited. The isoform audit reads the repaired BRCA isoform matrix and the repaired UCEC gene-RNA matrix. PDAC CNV is not an input to the isoform-inverse audit.
