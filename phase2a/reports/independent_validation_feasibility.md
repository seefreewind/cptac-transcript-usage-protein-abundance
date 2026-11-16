# Independent validation feasibility

The search covered GEO, EGA, and the GDC publication page for cohorts that could test transcript usage outside the 1,027 CPTAC cases. No cohort met the locked PASS rule: cancer and control groups, patient identifiers, transcript quantification, GENCODE v34-compatible annotation, at least 50 cases, and no CPTAC sample overlap.

GSE218035 is a public osteosarcoma long-read RNA series with 23 patients and 13 adjacent samples. It is independent of the CPTAC case set and is too small for the locked threshold.

EGAS00001004819 is the Veiga breast long-read resource. Its peptide support uses CPTAC spectra, so it is not an independent protein validation.

GSE270257 contains long-read DNA methylation for 189 POG tumors. The GEO object does not contain transcript quantification. Raw reads are in EGA EGAS00001001159.

The TCGA splicing supplement covers 8,705 patients and has partial sample overlap with TCGA-linked CPTAC cohorts. Its splicing graph is not the frozen GENCODE v34 principal rule.

GSE293661 is a ten-cell-line pancreatic long-read series and has no patient control arm.

An independent RNA cohort can test transcript usage. It cannot validate protein abundance. None of these cohorts is combined with the 1,027 cases, and no DTU test was run.

Validation grade: PARTIAL.
