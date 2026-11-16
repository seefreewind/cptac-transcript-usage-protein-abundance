# Transcript universe validation

GENCODE v34 GTF hash matches `941926bde63e86f4eebbd29a8f747937288750da6b7dda8b21a920980c38b6e1`.

The ten isoform matrices share one transcript-identifier set. Row order differs by cancer. The sorted-set SHA256 is `358c2666fe9833a2f3082733ad079c1b310058be4ded983f4c13bf7d0323b2ea`.

- matrix_rows: observed 108777, expected 108777, match 1
- matrix_genes: observed 60669, expected 60669, match 1
- genes_one_row: observed 42424, expected 42424, match 1
- genes_multi_row: observed 18245, expected 18245, match 1
- unmapped_ids: observed 0, expected 0, match 1
- duplicate_ids: observed 0, expected 0, match 1
- mane: observed 14620, expected 14620, match 1
- appris: observed 2819, expected 2819, match 1
- unassigned_tagged: observed 1001, expected 1001, match 1
- principal: observed 17439, expected 17439, match 1
- principal_absent: observed 37, expected 37, match 1
- candidate: observed 10795, expected 10795, match 1
- frozen_principal_table: observed 0, expected 0, match 1

The principal transcript was assigned by the frozen MANE Select then unique APPRIS principal-1 rule. Partial R², cancer type, and file row order were not used.
