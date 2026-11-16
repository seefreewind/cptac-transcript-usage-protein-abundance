# Cancer-Specific Heterogeneity in the Incremental Explanatory Value of Transcript Usage for Protein Abundance

Analysis code and frozen derived results supporting the manuscript. The study includes 1,027 case-level matched tumor observations across BRCA, CCRCC, COAD, GBM, HNSCC, LSCC, LUAD, OV, PDAC and UCEC. The outcome is gene-level protein abundance. Both co-primary transcript-expression inverse transformations are retained.

## Contents

- `phase1a` contains the nested-model implementation, protocol, principal-transcript map, gene-fit archives, cancer summaries, sensitivity outputs and random-effects synthesis.
- `phase1b` contains distribution, measurement, influence and heterogeneity audits.
- `phase2a` contains reference transcript/protein annotation linkage and coverage records.
- `phase2b` contains the secondary predicted structural-annotation comparison and its sensitivities.
- `phase0_5` and `phase0_6` contain acquisition/reconstruction scripts, source URLs, checksum metadata and preprocessing audits.
- `phase4c` contains a standalone frozen-result figure exporter, figure PDFs and previews, editable source tables, legends and source manifests.
- `dataset_manifest.csv` records upstream dataset locations and cohort mapping.
- `FILE_MANIFEST.csv` records the size and SHA-256 of every archived payload file.

## Reproduce figure exports

From the repository root, using Python 3.11 or later, run

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python phase4c/code/export_figures.py
```

This entry point reads the archived result files and exports the five main figures and three supplementary figures. It preserves stored estimates and sorts every successful gene-fit value for ECDF visualization. Main-table CSV and XLSX files contain the displayed values. Full-precision fits and uncertainty fields are retained in the phase results.

## Historical analysis entry points

```bash
python phase0_5/code/download_matrices.py
python phase0_5/code/rebuild_cohort.py
python phase0_6/code/quantification_audit.py
python phase1a/code/run_phase1a.py
python phase1b/code/audit_phase1b.py
python phase2a/code/audit_phase2a.py
python phase2b/code/run_phase2b.py
```

These are archived historical entry points. Full model reruns require restoration of the upstream matrices and versioned annotations at their original relative locations, and satisfaction of the historical protocol gates. Inspect the acquisition manifests, input manifests and phase reports first. Raw bulk matrices are downloadable from the recorded LinkedOmicsKB S3 URLs. GENCODE v34 and Ensembl release 100 reference URLs appear in the annotation audit code and input manifests. Archived reports can mention historical local paths and prefit gates. The result-only figure exporter is independently runnable with this repository payload.

## Interpretation

The pooled mean partial R² is approximately 0.0248 under Literal A and 0.0240 under Literal B. Prediction intervals include values below the prespecified 0.02 threshold. Structural comparisons use the separate 0.005 threshold. Case matching is at case level, with aliquot correspondence unresolved. The analysis estimates incremental in-sample explanation of gene-level protein abundance. Product-specific functional consequences require additional evidence.

## Provenance and reuse

Code authored for this project is distributed under MIT. Original derived numerical results and figure/table source data are distributed under CC BY 4.0, https://creativecommons.org/licenses/by/4.0/. Reused annotations and public source metadata retain their upstream terms. Original large expression matrices and controlled-access sequencing data are obtained through their upstream providers.

The historical Phase 2A run log had a preexisting checksum discrepancy. That log is outside this archive. Scientific payload files were copied byte-for-byte from the frozen project outputs. `FILE_MANIFEST.csv` provides self-contained archive checksums.

Zenodo creator metadata use the complete project title as requested. The repository name contains no journal name.
