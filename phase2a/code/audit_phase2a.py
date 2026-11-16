#!/usr/bin/env python3
"""Phase 2A transcript-to-protein annotation qualification.

Reads frozen Phase 1A gene identifiers and GENCODE v34 / Ensembl 100
annotation. It does not refit protein models and does not rank genes.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import os
import platform
import sys
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "phase1a" / "code"))
from common import CANCERS, GTF, GTF_SHA, acquisition_path, load_annotation  # noqa: E402

PHASE = ROOT / "phase2a"
RAW = PHASE / "raw"
GTF_PATH = GTF
SWISSPROT = ROOT / "phase0/raw/annotation/gencode.v34.metadata.SwissProt.gz"
TREMBL = RAW / "gencode.v34.metadata.TrEMBL.gz"
FASTA = RAW / "gencode.v34.pc_translations.fa.gz"
FEATURE = RAW / "protein_feature.txt.gz"
TRANSLATION = RAW / "translation.txt.gz"
ANALYSIS = RAW / "analysis.txt.gz"
FROZEN = ROOT / "phase1a/manifests/frozen_principal_transcripts.csv"
UNIPROT_COUNTS = ROOT / "phase0/raw/api/uniprot_counts.json"

EXPECTED = {
    "matrix_rows": 108777,
    "matrix_genes": 60669,
    "genes_one_row": 42424,
    "genes_multi_row": 18245,
    "unmapped_ids": 0,
    "duplicate_ids": 0,
    "mane": 14620,
    "appris": 2819,
    "unassigned_tagged": 1001,
    "principal": 17439,
    "principal_absent": 37,
    "candidate": 10795,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def log(message: str) -> None:
    line = f"{utc_now()} {message}"
    print(line, flush=True)
    with (PHASE / "logs" / "run.log").open("a") as handle:
        handle.write(line + "\n")


class Transcript:
    __slots__ = (
        "gene", "biotype", "protein_ids", "chrom", "strand", "exons",
        "cds_nt", "cds_min", "cds_max", "has_start", "has_stop",
    )

    def __init__(self):
        self.gene = ""
        self.biotype = ""
        self.protein_ids: set[str] = set()
        self.chrom = ""
        self.strand = ""
        self.exons: list[tuple[int, int]] = []
        self.cds_nt = 0
        self.cds_min = None
        self.cds_max = None
        self.has_start = False
        self.has_stop = False


def attr_values(attr: str, key: str) -> list[str]:
    token = key + ' "'
    values = []
    start = 0
    while True:
        index = attr.find(token, start)
        if index < 0:
            return values
        end = attr.find('"', index + len(token))
        if end < 0:
            return values
        values.append(attr[index + len(token):end])
        start = end + 1


def parse_gtf() -> dict[str, Transcript]:
    records: dict[str, Transcript] = {}
    with gzip.open(GTF_PATH, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9:
                continue
            kind = parts[2]
            if kind not in {"transcript", "exon", "CDS", "start_codon", "stop_codon"}:
                continue
            attr = parts[8]
            txs = attr_values(attr, "transcript_id")
            if not txs:
                continue
            tx_id = txs[0]
            record = records.get(tx_id)
            if record is None:
                record = Transcript()
                records[tx_id] = record
            if kind == "transcript":
                genes = attr_values(attr, "gene_id")
                biotypes = attr_values(attr, "transcript_type")
                record.gene = genes[0] if genes else ""
                record.biotype = biotypes[0] if biotypes else ""
                record.chrom = parts[0]
                record.strand = parts[6]
                record.protein_ids.update(attr_values(attr, "protein_id"))
            else:
                start = int(parts[3])
                end = int(parts[4])
                if kind == "exon":
                    record.exons.append((start, end))
                elif kind == "CDS":
                    record.cds_nt += end - start + 1
                    record.cds_min = start if record.cds_min is None else min(record.cds_min, start)
                    record.cds_max = end if record.cds_max is None else max(record.cds_max, end)
                    record.protein_ids.update(attr_values(attr, "protein_id"))
                elif kind == "start_codon":
                    record.has_start = True
                elif kind == "stop_codon":
                    record.has_stop = True
    return records


def read_matrix_ids(path: Path) -> list[str]:
    identifiers = []
    with path.open("rt") as handle:
        handle.readline()
        for line in handle:
            identifiers.append(line.split("\t", 1)[0].rstrip("\r"))
    return identifiers


def load_crosswalk(path: Path, source: str) -> dict[str, list[tuple[str, str, str]]]:
    found: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            left = parts[1] if len(parts) > 1 else ""
            right = parts[2] if len(parts) > 2 else ""
            found[parts[0]].append((source, left, right))
    return found


def load_fasta() -> dict[str, dict]:
    sequences: dict[str, dict] = {}
    duplicates = 0
    with gzip.open(FASTA, "rt") as handle:
        header = None
        chunks: list[str] = []

        def store() -> None:
            nonlocal duplicates
            if header is None:
                return
            body = "".join(chunks)
            fields = header[1:].split("|")
            ensp = fields[0]
            record = {
                "ensp": ensp,
                "enst": fields[1] if len(fields) > 1 else "",
                "ensg": fields[2] if len(fields) > 2 else "",
                "length": len(body),
                "sha256": sha_text(body),
                "sequence": body,
            }
            if ensp in sequences:
                duplicates += 1
                if sequences[ensp]["sequence"] != body:
                    sequences["__sequence_conflicts__"] = {"count": sequences.get("__sequence_conflicts__", {"count": 0})["count"] + 1}
            else:
                sequences[ensp] = record

        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                header = line
                chunks = []
            else:
                chunks.append(line)
        store()
    sequences["__duplicate_headers__"] = {"count": duplicates}  # type: ignore[assignment]
    return sequences


def split_version(identifier: str) -> tuple[str, int | None]:
    if "." not in identifier:
        return identifier, None
    stable, version = identifier.rsplit(".", 1)
    if version.isdigit():
        return stable, int(version)
    return identifier, None


def is_int(value: str) -> bool:
    return value.lstrip("-").isdigit()


def load_ensembl(needed_stable: set[str]) -> dict:
    analysis_lines = []
    with gzip.open(ANALYSIS, "rt") as handle:
        for line in handle:
            analysis_lines.append(line.rstrip("\n").split("\t"))
    known = {
        "pfam", "smart", "tmhmm", "signalp", "ncoils", "seg", "superfamily",
        "prints", "pirsf", "tigrfam", "gene3d", "hmmpanther", "panther",
    }
    width = max(len(row) for row in analysis_lines[:50])
    logic_col = max(
        range(width),
        key=lambda index: sum(
            1 for row in analysis_lines[:200]
            if index < len(row) and row[index].lower() in known
        ),
    )
    if not any(row[logic_col].lower() in known for row in analysis_lines[:200] if logic_col < len(row)):
        raise SystemExit("Ensembl analysis table has no recognizable logic_name column")
    analysis = {}
    for row in analysis_lines:
        if not row or not is_int(row[0]):
            continue
        logic = row[logic_col].lower() if logic_col < len(row) else ""
        database = row[logic_col + 1] if logic_col + 1 < len(row) else ""
        database_version = row[logic_col + 2] if logic_col + 2 < len(row) else ""
        program = row[logic_col + 4] if logic_col + 4 < len(row) else ""
        program_version = row[logic_col + 5] if logic_col + 5 < len(row) else ""
        analysis[int(row[0])] = {
            "logic_name": logic,
            "db": database,
            "db_version": database_version,
            "program": program,
            "program_version": program_version,
        }
    translations = {}
    stable_rows: dict[str, list[tuple[int, int | None]]] = defaultdict(list)
    with gzip.open(TRANSLATION, "rt") as handle:
        sample = []
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(sample) < 30:
                sample.append(row)
            ensp_cols = [
                index for index, value in enumerate(row)
                if value.startswith("ENSP")
            ]
            if not ensp_cols or not is_int(row[0]):
                continue
            stable = row[ensp_cols[0]]
            version = None
            version_index = ensp_cols[0] + 1
            if version_index < len(row) and is_int(row[version_index]):
                version = int(row[version_index])
            translation_id = int(row[0])
            if stable in needed_stable:
                stable_rows[stable].append((translation_id, version))
                translations[translation_id] = (stable, version)
    if not translations:
        raise SystemExit("Ensembl translation table did not join to any needed ENSP identifier")
    wanted_ids = set(translations)
    features: dict[int, list[tuple[str, str, int | None, int | None]]] = defaultdict(list)
    logic_counts: dict[str, int] = defaultdict(int)
    hit_col = None
    with gzip.open(FEATURE, "rt") as handle:
        for line_number, line in enumerate(handle):
            if line_number > 20000:
                break
            row = line.rstrip("\n").split("\t")
            for index, value in enumerate(row):
                if value.startswith("PF") or value in {"SIGNAL", "TMhelix", "SignalP"}:
                    hit_col = index
                    break
            if hit_col is not None:
                break
    if hit_col is None:
        raise SystemExit("protein_feature hit column was not detected")
    with gzip.open(FEATURE, "rt") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if hit_col >= len(row) or len(row) < 2 or not is_int(row[1]):
                continue
            translation_id = int(row[1])
            if translation_id not in wanted_ids:
                continue
            analysis_index = hit_col + 1
            if analysis_index >= len(row) or not is_int(row[analysis_index]):
                continue
            logic = analysis.get(int(row[analysis_index]), {}).get("logic_name", "")
            logic_counts[logic] += 1
            if logic not in {"pfam", "signalp", "tmhmm", "catalytic_residue", "active_site"}:
                continue
            start = int(row[hit_col - 4]) if hit_col >= 4 and is_int(row[hit_col - 4]) else None
            end = int(row[hit_col - 3]) if hit_col >= 3 and is_int(row[hit_col - 3]) else None
            features[translation_id].append((logic, row[hit_col], start, end))
    return {
        "analysis": analysis,
        "logic_col": logic_col,
        "stable_rows": stable_rows,
        "features": features,
        "logic_counts": dict(logic_counts),
        "hit_col": hit_col,
    }


def overlaps(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return left[0] <= right[1] and right[0] <= left[1]


def has_private_exon(source: list[tuple[int, int]], other: list[tuple[int, int]]) -> bool:
    return any(not any(overlaps(exon, candidate) for candidate in other) for exon in source)


def boundary(record: Transcript, side: str) -> int | None:
    if record.cds_min is None or record.cds_max is None:
        return None
    if record.strand == "+":
        return record.cds_min if side == "start" else record.cds_max
    if record.strand == "-":
        return record.cds_max if side == "start" else record.cds_min
    return None


def pfam_token(hit: str) -> str:
    token = hit.split(".", 1)[0]
    return token if token.startswith("PF") else hit


def supported(start: int | None, end: int | None, length: int) -> bool:
    return start is not None and end is not None and 1 <= start <= end <= length


def feature_status(feature_rows: list[tuple[str, str, int | None, int | None]], length: int) -> dict:
    pfam = set()
    pfam_bad = 0
    pfam_n = 0
    signal = 0
    signal_bad = 0
    tm = 0
    tm_bad = 0
    catalytic = 0
    for logic, hit, start, end in feature_rows:
        ok = supported(start, end, length)
        if logic == "pfam":
            pfam_n += 1
            if ok:
                pfam.add(pfam_token(hit))
            else:
                pfam_bad += 1
        elif logic == "signalp":
            if ok:
                signal += 1
            else:
                signal_bad += 1
        elif logic == "tmhmm":
            if ok:
                tm += 1
            else:
                tm_bad += 1
        elif logic in {"catalytic_residue", "active_site"} and ok:
            catalytic += 1
    return {
        "pfam": pfam,
        "pfam_bad": pfam_bad,
        "pfam_n": pfam_n,
        "signal": signal,
        "signal_bad": signal_bad,
        "tm": tm,
        "tm_bad": tm_bad,
        "catalytic": catalytic,
    }


def domain_class(principal_set: set[str], alternative_set: set[str]) -> str:
    only_principal = principal_set - alternative_set
    only_alternative = alternative_set - principal_set
    if only_principal and only_alternative:
        return "mixed"
    if only_principal:
        return "principal_only"
    if only_alternative:
        return "alternative_only"
    if principal_set:
        return "both_only"
    return "neither"


def map_protein(ensembl: dict, protein_id: str, sequence: dict | None) -> dict:
    stable, version = split_version(protein_id)
    rows = ensembl["stable_rows"].get(stable, [])
    result = {
        "stable": stable,
        "version": version,
        "n_rows": len(rows),
        "version_match": 0,
        "translation_id": "",
        "status": "no_annotation_available",
        "features": [],
    }
    if sequence is None or not rows:
        return result
    matched = [row for row in rows if row[1] == version]
    if len(rows) != 1 or version is None or len(matched) != 1:
        result["status"] = "mapping_uncertain"
        return result
    result["version_match"] = 1
    result["translation_id"] = matched[0][0]
    result["features"] = ensembl["features"].get(matched[0][0], [])
    result["status"] = "version_matched"
    return result


def synthetic_tests() -> None:
    assert split_version("ENSP00000000001.12") == ("ENSP00000000001", 12)
    assert domain_class({"PF1"}, {"PF1", "PF2"}) == "alternative_only"
    assert domain_class({"PF1", "PF3"}, {"PF1", "PF2"}) == "mixed"
    assert domain_class({"PF1"}, {"PF1"}) == "both_only"
    assert domain_class(set(), set()) == "neither"
    assert supported(1, 4, 4)
    assert not supported(1, 5, 4)
    assert sha_text("MKA") == hashlib.sha256(b"MKA").hexdigest()
    (PHASE / "qc").mkdir(parents=True, exist_ok=True)
    (PHASE / "qc" / "synthetic_test_results.txt").write_text("SYNTHETIC_TESTS_PASS\n")


def verify_historical_hashes() -> list[dict]:
    rows = []
    targets = [
        ROOT / "phase0_6/reports/statistical_model_freeze.md",
        ROOT / "phase0_6/reports/transformation_resolution_report.md",
        ROOT / "phase0_6/manifests/final_data_readiness.csv",
        ROOT / "phase1a/reports/phase1a_final_decision.md",
        ROOT / "phase1a/reports/phase1a_primary_results.md",
        ROOT / "phase1b/reports/phase1b_final_decision.md",
        ROOT / "phase1b/reports/phase1b_gene_distribution_report.md",
        ROOT / "phase1b/reports/phase1b_transcript_usage_report.md",
        ROOT / "phase1b/reports/phase1b_transformation_robustness.md",
        GTF_PATH,
    ]
    ledgers = {}
    for ledger_path in (
        ROOT / "phase0_6/SHA256SUMS",
        ROOT / "phase1a/SHA256SUMS",
        ROOT / "phase1b/SHA256SUMS",
    ):
        if not ledger_path.exists():
            continue
        for line in ledger_path.read_text().splitlines():
            if not line.strip():
                continue
            digest, rel = line.split(None, 1)
            ledgers[rel.strip()] = digest
    gtf_digest = sha256_file(GTF_PATH)
    if gtf_digest != GTF_SHA:
        raise SystemExit("GENCODE v34 GTF hash does not match the frozen value")
    for path in targets:
        digest = gtf_digest if path == GTF_PATH else sha256_file(path)
        rel = path.relative_to(ROOT).as_posix()
        recorded = ledgers.get(rel, "")
        rows.append({
            "path": rel,
            "sha256": digest,
            "ledger_sha256": recorded,
            "ledger_match": int(bool(recorded) and recorded == digest),
            "in_ledger": int(bool(recorded)),
        })
    return rows


def read_ok_genes() -> dict[tuple[str, str], set[str]]:
    grouped = {}
    for transform, name in (("A", "primary_literal_A_gene_fits.tsv.gz"), ("B", "primary_literal_B_gene_fits.tsv.gz")):
        with gzip.open(ROOT / "phase1a/results" / name, "rt") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            cancer_i = header.index("cancer")
            gene_i = header.index("gene_id")
            status_i = header.index("fit_status")
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if parts[status_i] != "OK":
                    continue
                grouped.setdefault((parts[cancer_i], transform), set()).add(parts[gene_i])
    return grouped


def expected_success() -> dict[tuple[str, str], tuple[int, int]]:
    found = {}
    with (ROOT / "phase1a/results/eligibility_by_cancer_transform.csv").open() as handle:
        for row in csv.DictReader(handle):
            found[(row["cancer"], row["transformation"])] = (int(row["n_success"]), int(row["n_candidate"]))
    paired = {}
    with (ROOT / "phase1a/results/common_support_summary.csv").open() as handle:
        for row in csv.DictReader(handle):
            paired[row["cancer"]] = int(row["n_paired"])
    return {"success": found, "paired": paired}


def file_row(path: Path, source: str, release: str, retrieved: str, url: str, last_modified: str, kind: str) -> dict:
    return {
        "source": source,
        "release": release,
        "retrieval_date": retrieved,
        "url_or_accession": url,
        "file_size_bytes": path.stat().st_size if path.exists() else 0,
        "sha256": sha256_file(path) if path.exists() else "",
        "annotation_type": kind,
        "species": "Homo sapiens",
        "identifier_versioning": "version suffixes retained; stable id stored separately when a join requires it",
        "source_last_modified": last_modified,
        "used_for_domain_coordinates": int(source.startswith("Ensembl 100")),
    }


def literature_rows() -> list[dict]:
    return [
        {
            "title": "Comprehensive Analysis of Alternative Splicing Across Tumors from 8,705 Patients",
            "authors": "Kahles A, Lehmann KV, Toussaint NC, et al.",
            "journal": "Cancer Cell",
            "year": 2018,
            "doi": "10.1016/j.ccell.2018.07.001",
            "data_resources": "TCGA RNA-seq and WES, 8,705 patients, 32 cancer types; CPTAC spectra used for selected neojunction peptides in breast and ovarian tumors",
            "main_estimand": "Alternative splicing events and tumor-specific exon junctions, plus somatic associations",
            "methods": "RNA and whole-exome reanalysis; splicing-event counts; mass-spectrometry confirmation of selected peptides",
            "main_conclusions": "Tumors carry more splicing events than normal samples, and a subset of neojunction peptides is detectable",
            "direct_overlap": "Pan-cancer splicing and a limited CPTAC peptide confirmation. The estimand is event detection, not incremental partial R2 of a prespecified principal proportion",
            "remaining_gap": "Not established by audited literature",
            "topic_status": "Partially addressed",
        },
        {
            "title": "The Functional Impact of Alternative Splicing in Cancer",
            "authors": "Climente-Gonzalez H, Porta-Pardo E, Godzik A, Eyras E",
            "journal": "Cell Reports",
            "year": 2017,
            "doi": "10.1016/j.celrep.2017.08.012",
            "data_resources": "TCGA tumor RNA-seq",
            "main_estimand": "Predicted functional impact of cancer-associated isoform switches on domains and interactions",
            "methods": "Isoform-switch detection followed by Pfam, Prosite, and disorder annotation",
            "main_conclusions": "A subset of switches changes annotated protein features and protein-interaction potential",
            "direct_overlap": "Predicted domain change from transcript switches. Switches were selected from tumor-versus-normal usage, which this audit does not do",
            "remaining_gap": "Partially addressed",
            "topic_status": "Partially addressed",
        },
        {
            "title": "The Landscape of Isoform Switches in Human Cancers",
            "authors": "Vitting-Seerup K, Sandelin A",
            "journal": "Molecular Cancer Research",
            "year": 2017,
            "doi": "10.1158/1541-7786.MCR-16-0459",
            "data_resources": "TCGA RNA-seq, more than 5,500 patients, 12 solid cancer types",
            "main_estimand": "Isoform switches with predicted functional consequences, including domain loss and survival association",
            "methods": "Isoform-switch testing and functional-consequence annotation",
            "main_conclusions": "Predicted functional switches are common across those 12 cancer types",
            "direct_overlap": "Pan-cancer transcript-usage switching. The outcome is RNA isoform usage, not protein abundance conditional on gene expression, purity, and CNV",
            "remaining_gap": "Partially addressed",
            "topic_status": "Partially addressed",
        },
        {
            "title": "IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences",
            "authors": "Vitting-Seerup K, Sandelin A",
            "journal": "Bioinformatics",
            "year": 2019,
            "doi": "10.1093/bioinformatics/btz247",
            "data_resources": "Software paper; TCGA used as a case study",
            "main_estimand": "A method for genome-wide isoform-switch and consequence annotation",
            "methods": "R package for differential isoform usage and predicted consequences",
            "main_conclusions": "The package systematizes switch detection and consequence prediction",
            "direct_overlap": "Methods for DTU. Phase 2A does not run that test",
            "remaining_gap": "Partially addressed",
            "topic_status": "Partially addressed",
        },
        {
            "title": "Proteogenomics connects somatic mutations to signalling in breast cancer",
            "authors": "Mertins P, Mani DR, Ruggles KV, Gillette MA, et al.",
            "journal": "Nature",
            "year": 2016,
            "doi": "10.1038/nature18003",
            "data_resources": "CPTAC breast tumors, 105 genomically annotated cases, 77 high-quality proteomes",
            "main_estimand": "Gene-level and phosphosite proteomic associations with somatic alterations",
            "methods": "Mass spectrometry integrated with genomic annotation",
            "main_conclusions": "Gene-level proteogenomics connects mutations, copy-number changes, and signalling",
            "direct_overlap": "CPTAC breast protein abundance. The protein measurement is gene-level. The paper also reports splice-associated peptide observations in the broader CPTAC literature, but this article's estimand is not principal-transcript partial R2",
            "remaining_gap": "Partially addressed",
            "topic_status": "Partially addressed",
        },
        {
            "title": "A comprehensive long-read isoform analysis platform and sequencing resource for breast cancer",
            "authors": "Veiga DFT, Nesta A, Zhao Y, Deslattes Mays A, et al.",
            "journal": "Science Advances",
            "year": 2022,
            "doi": "10.1126/sciadv.abg6711",
            "data_resources": "Long-read breast tumor RNA-seq; CPTAC breast mass spectra used for peptide support; EGA EGAS00001004819",
            "main_estimand": "Full-length isoform discovery and tumor-specific splicing events",
            "methods": "Long-read RNA-seq, ORF prediction, spectral matching to CPTAC peptides",
            "main_conclusions": "Long-read sequencing finds unannotated isoforms, and some predicted ORFs have peptide support in CPTAC spectra",
            "direct_overlap": "Breast cancer isoform structure and CPTAC peptide matching. The cohort is not independent of CPTAC proteomics, and the estimand is isoform discovery",
            "remaining_gap": "Partially addressed",
            "topic_status": "Partially addressed",
        },
        {
            "title": "Pan-cancer characterization of protein isoforms uncovers driving alternative splicing events at protein level",
            "authors": "Conference abstract; author list is on the AACR abstract page",
            "journal": "Cancer Research, AACR Annual Meeting abstract 7687",
            "year": 2026,
            "doi": "10.1158/1538-7445.AM2026-7687",
            "data_resources": "Ten CPTAC tumor cohorts and the CPTAC pan-cancer harmonized dataset",
            "main_estimand": "Protein-isoform or isoform-group abundance from peptide-to-isoform graphs",
            "methods": "Tripartite graph grouping of shotgun peptides, then integration with mutation, CNV, methylation, expression, and gene-level protein",
            "main_conclusions": "The abstract reports isoform-group matrices and associations in the same ten CPTAC cohorts",
            "direct_overlap": "Same CPTAC cohorts and a protein-isoform estimand. This is a meeting abstract, not a finished method that replaces the frozen gene-level partial R2",
            "remaining_gap": "Partially addressed",
            "topic_status": "Partially addressed",
        },
        {
            "title": "Frozen project estimand: incremental partial R2 of principal-transcript proportion",
            "authors": "This project",
            "journal": "Not a published paper; recorded so the gap has an explicit row",
            "year": 2026,
            "doi": "",
            "data_resources": "Ten CPTAC cohorts, 1,027 case-level triples, GENCODE v34 principal transcripts",
            "main_estimand": "Partial R2 of principal-transcript proportion beyond reconstructed gene expression, WES purity, and CNV",
            "methods": "Cancer-specific nested OLS, literal inverses A and B, random effects",
            "main_conclusions": "Phase 1A found a small average increment whose magnitude differs by cohort. Phase 1B left the heterogeneity mixed",
            "direct_overlap": "This is the project itself",
            "remaining_gap": "Not established by audited literature",
            "topic_status": "Not established by audited literature",
        },
    ]


def validation_rows() -> list[dict]:
    return [
        {
            "accession": "GSE218035",
            "cancer_type": "osteosarcoma",
            "n_cases": "23 patients, as stated in the GEO summary retrieved 2026-11-14",
            "n_controls": "13 adjacent normal samples, as stated in that GEO summary",
            "platform": "Oxford Nanopore long-read RNA-seq",
            "read_type": "long-read",
            "processed_data": "GEO series is public",
            "transcript_annotation_version": "Not stated as GENCODE v34 in the retrieved summary",
            "dtu_can_be_independently_tested": "possible_from_long_read_rna_but_not_started",
            "patient_identifiers_usable": "patient-level tumor and adjacent tissue are described",
            "independent_of_cptac": "yes_for_the_1027_case_cohort",
            "access_restrictions": "GEO public",
            "protein_validation": "An RNA cohort does not validate protein abundance",
            "notes": "Sample size is below the locked threshold of 50 cases, and the cancer type is outside the ten CPTAC cohorts",
        },
        {
            "accession": "EGAS00001004819",
            "cancer_type": "breast",
            "n_cases": "breast tumor long-read resource described by Veiga et al. 2022",
            "n_controls": "not stated as a matched control arm in the retrieved record",
            "platform": "long-read and short-read RNA-seq",
            "read_type": "long-read",
            "processed_data": "controlled-access study",
            "transcript_annotation_version": "Authors compare isoforms with GENCODE; the freeze version is not established as v34",
            "dtu_can_be_independently_tested": "not_from_open_files",
            "patient_identifiers_usable": "controlled access",
            "independent_of_cptac": "no_protein_arm_uses_CPTAC_spectra",
            "access_restrictions": "EGA controlled access",
            "protein_validation": "Peptide support in that paper uses CPTAC spectra, so it is not an independent protein validation",
            "notes": "Recorded because it is a major breast long-read isoform resource and because it is not independent",
        },
        {
            "accession": "GSE270257",
            "cancer_type": "advanced cancers, multiple types, POG program",
            "n_cases": "189 tumors from 181 patients",
            "n_controls": "41 matched blood samples",
            "platform": "Oxford Nanopore PromethION",
            "read_type": "long-read DNA methylation in the GEO object; accompanying RNA is short-read",
            "processed_data": "methylation frequencies in GEO; sequence reads in EGA EGAS00001001159",
            "transcript_annotation_version": "No transcript quantification matrix in the GEO object",
            "dtu_can_be_independently_tested": "no_from_the_GEO_object",
            "patient_identifiers_usable": "clinical data exist under controlled access",
            "independent_of_cptac": "yes_for_the_1027_case_cohort",
            "access_restrictions": "GEO methylation is public; raw sequence is EGA controlled",
            "protein_validation": "This object does not measure protein abundance",
            "notes": "Useful as a long-read cancer genome resource. It is not a transcript-usage validation set",
        },
        {
            "accession": "GDC PanCanAtlas-Splicing-2018",
            "cancer_type": "32 TCGA cancer types",
            "n_cases": "8705 patients",
            "n_controls": "670 matched normal samples in the Cancer Cell paper",
            "platform": "TCGA short-read RNA-seq",
            "read_type": "short-read",
            "processed_data": "GDC publication page for the splicing supplement",
            "transcript_annotation_version": "Study-specific splicing graph, not the frozen GENCODE v34 principal rule",
            "dtu_can_be_independently_tested": "splicing_calls_exist_but_are_a_different_estimand",
            "patient_identifiers_usable": "TCGA barcodes are public for open-access data",
            "independent_of_cptac": "partial_overlap_with_TCGA_CPTAC_breast_ovarian_colon_endometrial_cases",
            "access_restrictions": "open processed publication data; sequence access follows GDC rules",
            "protein_validation": "The splicing publication's peptide confirmation uses CPTAC spectra for a subset",
            "notes": "Not combined with the 1,027 CPTAC cases",
        },
        {
            "accession": "GSE293661",
            "cancer_type": "pancreatic cancer cell lines",
            "n_cases": "10 cell lines, two replicates",
            "n_controls": "no patient control group",
            "platform": "Oxford Nanopore long-read RNA-seq",
            "read_type": "long-read",
            "processed_data": "GEO public",
            "transcript_annotation_version": "Not stated as GENCODE v34",
            "dtu_can_be_independently_tested": "cell_line_resource_only",
            "patient_identifiers_usable": "no_patient_identifiers",
            "independent_of_cptac": "yes_for_the_1027_case_cohort",
            "access_restrictions": "GEO public",
            "protein_validation": "Cell-line RNA does not validate the CPTAC tumor protein model",
            "notes": "Recorded as a public long-read transcriptome and excluded as a patient validation cohort",
        },
    ]


def overall_label(grades: dict[str, str]) -> str:
    if grades["A"] == "FAIL" or grades["D"] == "FAIL":
        return "STOP_CURRENT_DIRECTION"
    if grades["B"] == "FAIL":
        return "PHASE2B_REQUIRES_REDESIGN"
    if all(grades[key] == "PASS" for key in ("A", "B", "C", "D")):
        return "QUALIFIED_FOR_PHASE2B"
    if grades["A"] == "PASS" and grades["D"] == "PASS" and grades["B"] in {"PASS", "PARTIAL"} and grades["C"] == "PARTIAL":
        return "QUALIFIED_FOR_NARROW_PHASE2B"
    return "PHASE2B_REQUIRES_REDESIGN"


def write_decision(grades: dict[str, str], label: str, audit: str, reason: str) -> None:
    text = f"""# Phase 2A final decision

## A. Transcript annotation completeness

{grades['A']}

## B. Transcript-to-protein mapping

{grades['B']}

## C. Protein domain annotation

{grades['C']}

## D. Phase 1A result linkage

{grades['D']}

## E. Independent validation feasibility

{grades['E']}

## F. Novelty qualification

{grades['F']}

## Overall

{label}

{reason}

No Level 3 isoform-specific peptide evidence is established. The stored SEPEP matrices contain sample-level quantities and do not contain peptide sequences or protein accessions. Gene-level CPTAC protein abundance remains a gene-level measurement.

A structural annotation result is a statement about predicted sequences and annotated coordinates. It is not evidence of protein dysfunction, and it is not isoform-specific protein abundance.

PHASE 2A AUDIT = {audit}

PHASE 2B AUTHORIZATION = FALSE
"""
    (PHASE / "reports" / "phase2a_final_decision.md").write_text(text)


def hash_outputs() -> None:
    lines = []
    for path in sorted(PHASE.rglob("*")):
        if not path.is_file() or path.name.startswith("._") or path.name == "SHA256SUMS":
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        lines.append(f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}")
    (PHASE / "SHA256SUMS").write_text("\n".join(lines) + "\n")


def main() -> None:
    (PHASE / "logs").mkdir(parents=True, exist_ok=True)
    (PHASE / "logs" / "run.log").write_text("")
    log("start")
    synthetic_tests()
    hash_rows = verify_historical_hashes()
    write_csv(
        PHASE / "qc" / "historical_hash_verification.csv",
        hash_rows,
        ["path", "sha256", "ledger_sha256", "ledger_match", "in_ledger"],
    )
    for required in (FASTA, TREMBL, FEATURE, TRANSLATION, ANALYSIS, SWISSPROT):
        if not required.exists() or required.stat().st_size < 1000:
            raise SystemExit(f"required annotation file is missing: {required}")
    log("parse gtf")
    tx_to_gene, principal, rule = load_annotation()
    records = parse_gtf()
    built_mismatch = sum(1 for gene, tx in principal.items() if records.get(tx) is None or records[tx].gene != gene)
    if built_mismatch:
        raise SystemExit(f"principal transcripts missing from the GTF parse: {built_mismatch}")
    log("read matrices")
    matrix_ids = None
    matrix_set = None
    matrix_hash = ""
    for cancer in CANCERS:
        identifiers = read_matrix_ids(acquisition_path(cancer, "rna_isoform"))
        current = set(identifiers)
        if len(current) != len(identifiers):
            raise SystemExit(f"{cancer} isoform matrix has duplicate transcript identifiers")
        if matrix_ids is None:
            matrix_ids = identifiers
            matrix_set = current
            matrix_hash = sha_text("\n".join(sorted(current)))
        elif current != matrix_set:
            raise SystemExit(f"{cancer} isoform identifier set differs from the BRCA reference set")
    assert matrix_ids is not None
    duplicate_ids = len(matrix_ids) - len(set(matrix_ids))
    unmapped = sum(1 for tx in matrix_ids if tx not in tx_to_gene)
    counts: dict[str, int] = defaultdict(int)
    for tx in matrix_ids:
        counts[tx_to_gene.get(tx, "")] += 1
    one_row = sum(1 for value in counts.values() if value == 1)
    multi_row = sum(1 for value in counts.values() if value >= 2)
    mane_n = sum(item == "MANE_Select" for item in rule.values())
    appris_n = sum(item == "APPRIS_principal_1" for item in rule.values())
    tagged = set()
    with gzip.open(GTF_PATH, "rt") as handle:
        for line in handle:
            if not line.startswith("#") and "\ttranscript\t" in line:
                if 'tag "MANE_Select"' in line or 'tag "appris_principal_1"' in line:
                    tagged.add(attr_values(line.split("\t", 8)[8], "gene_id")[0])
    unassigned = len(tagged - set(principal))
    absent = [gene for gene, tx in principal.items() if tx not in set(matrix_ids)]
    candidate_genes = sorted(
        gene for gene, tx in principal.items()
        if tx in set(matrix_ids) and counts.get(gene, 0) >= 2
    )
    observed = {
        "matrix_rows": len(matrix_ids),
        "matrix_genes": len(counts),
        "genes_one_row": one_row,
        "genes_multi_row": multi_row,
        "unmapped_ids": unmapped,
        "duplicate_ids": duplicate_ids,
        "mane": mane_n,
        "appris": appris_n,
        "unassigned_tagged": unassigned,
        "principal": len(principal),
        "principal_absent": len(absent),
        "candidate": len(candidate_genes),
    }
    universe_rows = []
    universe_ok = True
    for key, expected in EXPECTED.items():
        match = int(observed[key] == expected)
        universe_ok = universe_ok and bool(match)
        universe_rows.append({"metric": key, "value": observed[key], "expected": expected, "match": match})
    frozen_mismatch = 0
    frozen_rows = 0
    matrix_set = set(matrix_ids)
    with FROZEN.open() as handle:
        for row in csv.DictReader(handle):
            frozen_rows += 1
            gene = row["gene_id"]
            if principal.get(gene) != row["principal_transcript_id"] or rule.get(gene) != row["principal_rule"]:
                frozen_mismatch += 1
            candidate_flag = int(gene in set(candidate_genes))
            if int(row["in_candidate_universe"]) != candidate_flag:
                frozen_mismatch += 1
            if int(row["in_isoform_matrix"]) != int(row["principal_transcript_id"] in matrix_set):
                frozen_mismatch += 1
    universe_rows.append({
        "metric": "frozen_principal_table",
        "value": frozen_mismatch,
        "expected": 0,
        "match": int(frozen_mismatch == 0 and frozen_rows == 17439),
    })
    universe_ok = universe_ok and frozen_mismatch == 0 and frozen_rows == 17439
    write_csv(PHASE / "manifests" / "transcript_universe.csv", universe_rows, ["metric", "value", "expected", "match"])
    principal_rows = []
    for gene in sorted(principal):
        tx = principal[gene]
        principal_rows.append({
            "gene_id": gene,
            "principal_transcript_id": tx,
            "principal_rule": rule[gene],
            "in_isoform_matrix": int(tx in matrix_set),
            "n_isoform_rows": counts.get(gene, 0),
            "in_candidate_universe": int(gene in set(candidate_genes)),
            "assignment_status": "assigned",
        })
    write_csv(
        PHASE / "manifests" / "principal_transcript_mapping.csv",
        principal_rows,
        ["gene_id", "principal_transcript_id", "principal_rule", "in_isoform_matrix", "n_isoform_rows", "in_candidate_universe", "assignment_status"],
    )
    (PHASE / "reports" / "transcript_universe_validation.md").write_text(
        "# Transcript universe validation\n\n"
        f"GENCODE v34 GTF hash matches `{GTF_SHA}`.\n\n"
        f"The ten isoform matrices share one transcript-identifier set. Row order differs by cancer. The sorted-set SHA256 is `{matrix_hash}`.\n\n"
        + "\n".join(
            f"- {row['metric']}: observed {row['value']}, expected {row['expected']}, match {row['match']}"
            for row in universe_rows
        )
        + "\n\nThe principal transcript was assigned by the frozen MANE Select then unique APPRIS principal-1 rule. "
        "Partial R², cancer type, and file row order were not used.\n"
    )
    if not universe_ok:
        write_decision(
            {"A": "FAIL", "B": "FAIL", "C": "FAIL", "D": "FAIL", "E": "PARTIAL", "F": "PARTIAL"},
            "STOP_CURRENT_DIRECTION",
            "FAILED",
            "The frozen transcript universe did not reproduce. Annotation comparison stopped.",
        )
        hash_outputs()
        log("TECHNICAL_HOLD universe")
        return
    log("sequences and accessions")
    swiss = load_crosswalk(SWISSPROT, "GENCODE_v34_SwissProt")
    trembl = load_crosswalk(TREMBL, "GENCODE_v34_TrEMBL")
    fasta = load_fasta()
    fasta_duplicates = fasta.pop("__duplicate_headers__")["count"]
    fasta_conflicts = fasta.pop("__sequence_conflicts__", {"count": 0})["count"]
    by_gene: dict[str, list[str]] = defaultdict(list)
    for tx in matrix_ids:
        by_gene[tx_to_gene[tx]].append(tx)
    needed = set()
    for gene in candidate_genes:
        for tx in by_gene[gene]:
            record = records.get(tx)
            if record is None:
                continue
            for protein_id in record.protein_ids:
                needed.add(split_version(protein_id)[0])
    log(f"ensembl join for {len(needed)} stable protein ids")
    ensembl = load_ensembl(needed)
    log(f"ensembl logic counts {ensembl['logic_counts']}")

    mapping_fields = [
        "transcript_id", "gene_id", "transcript_biotype", "protein_id", "protein_id_stable",
        "protein_id_version", "uniprot_swissprot", "uniprot_trembl", "annotation_source",
        "annotation_release", "mapping_confidence", "sequence_sha256", "sequence_length",
        "sequence_header_transcript", "mapping_ambiguity", "in_candidate_gene",
    ]
    sequence_rows = {}
    mapping_handle_path = PHASE / "manifests" / "transcript_protein_mapping.csv"
    n_coding = n_coding_one_protein = n_coding_with_sequence = n_multi_sequence = 0
    n_header_conflict = 0
    with mapping_handle_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=mapping_fields, lineterminator="\n")
        writer.writeheader()
        for tx in matrix_ids:
            record = records.get(tx)
            gene = tx_to_gene[tx]
            protein_ids = sorted(record.protein_ids) if record else []
            protein_id = protein_ids[0] if len(protein_ids) == 1 else ""
            stable, version = split_version(protein_id) if protein_id else ("", None)
            sequence = fasta.get(protein_id)
            swiss_acc = sorted({item[1] for item in swiss.get(tx, [])})
            trembl_acc = sorted({item[1] for item in trembl.get(tx, [])})
            ambiguity = []
            if len(protein_ids) > 1:
                ambiguity.append("multiple_protein_ids")
            if len(swiss.get(tx, [])) > 1 or len(trembl.get(tx, [])) > 1:
                ambiguity.append("multiple_uniprot_rows")
            if swiss_acc and trembl_acc and set(swiss_acc) != set(trembl_acc):
                ambiguity.append("swissprot_and_trembl_differ")
            if sequence and sequence["enst"] and sequence["enst"] != tx:
                ambiguity.append("fasta_header_transcript_differs")
                n_header_conflict += 1
            if record and record.biotype == "protein_coding" and gene in set(candidate_genes):
                n_coding += 1
                if len(protein_ids) == 1:
                    n_coding_one_protein += 1
                    if sequence:
                        n_coding_with_sequence += 1
                if len(protein_ids) > 1:
                    n_multi_sequence += 1
            if sequence and protein_id not in sequence_rows:
                sequence_rows[protein_id] = {
                    "protein_id": protein_id,
                    "protein_id_stable": stable,
                    "protein_id_version": version if version is not None else "",
                    "sequence_source": "gencode.v34.pc_translations.fa.gz",
                    "sequence_release": "GENCODE v34",
                    "sequence_length": sequence["length"],
                    "sequence_sha256": sequence["sha256"],
                    "header_transcript_id": sequence["enst"],
                    "header_gene_id": sequence["ensg"],
                    "alignment_method": "exact_string_equality",
                    "alignment_version": "phase2a_exact_v1",
                }
            confidence = "no_protein_id"
            if len(protein_ids) > 1:
                confidence = "ambiguous_protein_id"
            elif protein_id and sequence and (swiss_acc or trembl_acc):
                confidence = "unique_sequence_and_uniprot"
            elif protein_id and sequence:
                confidence = "unique_sequence_no_uniprot"
            elif protein_id:
                confidence = "protein_id_without_sequence"
            writer.writerow({
                "transcript_id": tx,
                "gene_id": gene,
                "transcript_biotype": record.biotype if record else "",
                "protein_id": protein_id,
                "protein_id_stable": stable,
                "protein_id_version": version if version is not None else "",
                "uniprot_swissprot": "|".join(swiss_acc),
                "uniprot_trembl": "|".join(trembl_acc),
                "annotation_source": "GENCODE v34 GTF plus GENCODE v34 SwissProt/TrEMBL metadata plus GENCODE v34 translations",
                "annotation_release": "GENCODE v34; Ensembl 100 coordinate tables are separate",
                "mapping_confidence": confidence,
                "sequence_sha256": sequence["sha256"] if sequence else "",
                "sequence_length": sequence["length"] if sequence else "",
                "sequence_header_transcript": sequence["enst"] if sequence else "",
                "mapping_ambiguity": "|".join(ambiguity),
                "in_candidate_gene": int(gene in set(candidate_genes)),
            })
    sharing = defaultdict(int)
    for row in sequence_rows.values():
        sharing[row["sequence_sha256"]] += 1
    for row in sequence_rows.values():
        row["n_protein_ids_with_this_sequence"] = sharing[row["sequence_sha256"]]
    write_csv(
        PHASE / "manifests" / "protein_sequence_mapping.csv",
        [sequence_rows[key] for key in sorted(sequence_rows)],
        ["protein_id", "protein_id_stable", "protein_id_version", "sequence_source", "sequence_release", "sequence_length", "sequence_sha256", "header_transcript_id", "header_gene_id", "alignment_method", "alignment_version", "n_protein_ids_with_this_sequence"],
    )
    sequence_coverage = n_coding_with_sequence / n_coding_one_protein if n_coding_one_protein else 0
    one_protein_fraction = n_coding_one_protein / n_coding if n_coding else 0

    log("compare principal and alternative transcripts")
    structure_rows = []
    sequence_diff_rows = []
    functional_rows = []
    pair_class_counts = defaultdict(int)
    evidence_counts = defaultdict(int)
    gene_flags = {}
    domain_fields = [
        "gene_id", "principal_transcript_id", "alternative_transcript_id", "principal_protein_id",
        "alternative_protein_id", "principal_seq_sha256", "alternative_seq_sha256",
        "principal_seq_length", "alternative_seq_length", "sequences_identical",
        "alignment_method", "alignment_identity", "unmapped_residue_lower_bound",
        "domain_class", "pfam_shared", "pfam_only_principal", "pfam_only_alternative",
        "signal_principal", "signal_alternative", "tm_principal", "tm_alternative",
        "evidence_level", "coordinate_note",
    ]
    with (PHASE / "manifests" / "domain_mapping_audit.csv").open("w", newline="") as pair_handle:
        pair_writer = csv.DictWriter(pair_handle, fieldnames=domain_fields, lineterminator="\n")
        pair_writer.writeheader()
        candidate_set = set(candidate_genes)
        for gene in candidate_genes:
            ptx = principal[gene]
            principal_tx = records[ptx]
            protein_ids = sorted(principal_tx.protein_ids)
            principal_protein = protein_ids[0] if len(protein_ids) == 1 else ""
            principal_sequence = fasta.get(principal_protein)
            principal_map = map_protein(ensembl, principal_protein, principal_sequence)
            principal_features = feature_status(principal_map["features"], principal_sequence["length"]) if principal_sequence and principal_map["status"] == "version_matched" else None
            n_alt = n_coding_alt = n_nmd = n_noncoding = 0
            n_exon = n_cds = n_start = n_end = 0
            n_with_sequence = n_identical = n_different = n_missing = n_level1 = n_level2 = 0
            swiss_acc = swiss.get(ptx, [])
            trembl_acc = trembl.get(ptx, [])
            for tx in by_gene[gene]:
                if tx == ptx:
                    continue
                alt = records.get(tx)
                if alt is None:
                    continue
                n_alt += 1
                if alt.biotype == "nonsense_mediated_decay":
                    n_nmd += 1
                elif alt.biotype != "protein_coding":
                    n_noncoding += 1
                same_locus = alt.chrom == principal_tx.chrom and alt.strand == principal_tx.strand
                if same_locus and principal_tx.exons and alt.exons and (
                    has_private_exon(principal_tx.exons, alt.exons) or has_private_exon(alt.exons, principal_tx.exons)
                ):
                    n_exon += 1
                if principal_tx.cds_nt and alt.cds_nt and principal_tx.cds_nt != alt.cds_nt:
                    n_cds += 1
                if boundary(principal_tx, "start") and boundary(alt, "start") and boundary(principal_tx, "start") != boundary(alt, "start"):
                    n_start += 1
                if boundary(principal_tx, "end") and boundary(alt, "end") and boundary(principal_tx, "end") != boundary(alt, "end"):
                    n_end += 1
                alt_ids = sorted(alt.protein_ids)
                alt_protein = alt_ids[0] if len(alt_ids) == 1 else ""
                alt_sequence = fasta.get(alt_protein)
                if alt.biotype == "protein_coding" and alt.cds_nt and len(alt_ids) == 1:
                    n_coding_alt += 1
                if (
                    alt.biotype != "protein_coding"
                    or principal_tx.biotype != "protein_coding"
                    or not principal_sequence
                    or not alt_sequence
                ):
                    if alt.biotype == "protein_coding" and principal_protein and alt_protein:
                        n_missing += 1
                    evidence_counts["LEVEL_0"] += 1
                    continue
                identical = principal_sequence["sequence"] == alt_sequence["sequence"]
                n_with_sequence += 1
                if identical:
                    n_identical += 1
                else:
                    n_different += 1
                alt_map = map_protein(ensembl, alt_protein, alt_sequence)
                note = "version_matched"
                klass = "mapping_uncertain"
                pfam_shared = pfam_only_p = pfam_only_a = ""
                signal_p = signal_a = tm_p = tm_a = ""
                if principal_map["status"] != "version_matched" or alt_map["status"] != "version_matched":
                    klass = "no_annotation_available" if principal_map["status"] == "no_annotation_available" or alt_map["status"] == "no_annotation_available" else "mapping_uncertain"
                    note = principal_map["status"] + "|" + alt_map["status"]
                else:
                    alt_features = feature_status(alt_map["features"], alt_sequence["length"])
                    assert principal_features is not None
                    if principal_features["pfam_bad"] or alt_features["pfam_bad"] or principal_features["signal_bad"] or alt_features["signal_bad"] or principal_features["tm_bad"] or alt_features["tm_bad"]:
                        klass = "mapping_uncertain"
                        note = "unsupported_coordinates"
                    elif identical and principal_features["pfam"] != alt_features["pfam"]:
                        klass = "mapping_uncertain"
                        note = "identical_sequence_discordant_pfam"
                    else:
                        klass = domain_class(principal_features["pfam"], alt_features["pfam"])
                        pfam_shared = len(principal_features["pfam"] & alt_features["pfam"])
                        pfam_only_p = len(principal_features["pfam"] - alt_features["pfam"])
                        pfam_only_a = len(alt_features["pfam"] - principal_features["pfam"])
                        signal_p = int(principal_features["signal"] > 0)
                        signal_a = int(alt_features["signal"] > 0)
                        tm_p = principal_features["tm"]
                        tm_a = alt_features["tm"]
                region_diff = False
                if klass in {"principal_only", "alternative_only", "mixed"}:
                    region_diff = True
                if signal_p != "" and (signal_p != signal_a or tm_p != tm_a):
                    region_diff = True
                if identical:
                    level = "LEVEL_0_IDENTICAL_SEQUENCE"
                elif region_diff and klass != "mapping_uncertain":
                    level = "LEVEL_2"
                    n_level2 += 1
                else:
                    level = "LEVEL_1"
                    n_level1 += 1
                evidence_counts[level] += 1
                pair_class_counts[klass] += 1
                pair_writer.writerow({
                    "gene_id": gene,
                    "principal_transcript_id": ptx,
                    "alternative_transcript_id": tx,
                    "principal_protein_id": principal_protein,
                    "alternative_protein_id": alt_protein,
                    "principal_seq_sha256": principal_sequence["sha256"],
                    "alternative_seq_sha256": alt_sequence["sha256"],
                    "principal_seq_length": principal_sequence["length"],
                    "alternative_seq_length": alt_sequence["length"],
                    "sequences_identical": int(identical),
                    "alignment_method": "exact_string_equality",
                    "alignment_identity": 1 if identical else "",
                    "unmapped_residue_lower_bound": abs(principal_sequence["length"] - alt_sequence["length"]),
                    "domain_class": klass,
                    "pfam_shared": pfam_shared,
                    "pfam_only_principal": pfam_only_p,
                    "pfam_only_alternative": pfam_only_a,
                    "signal_principal": signal_p,
                    "signal_alternative": signal_a,
                    "tm_principal": tm_p,
                    "tm_alternative": tm_a,
                    "evidence_level": level,
                    "coordinate_note": note,
                })
            reliable = int(
                principal_sequence is not None
                and principal_map["status"] == "version_matched"
                and principal_features is not None
                and principal_features["pfam_bad"] == 0
                and principal_features["signal_bad"] == 0
                and principal_features["tm_bad"] == 0
            )
            pfam_n = len(principal_features["pfam"]) if principal_features else 0
            flags = {
                "unique_principal": 1,
                "valid_alt_protein_coding": int(n_coding_alt > 0),
                "distinct_alt_sequence": int(n_different > 0),
                "uniprot": int(bool(swiss_acc or trembl_acc)),
                "pfam": int(pfam_n > 0),
                "reliable_coordinates": reliable,
                "interpretable_sequence_difference": int(n_level1 + n_level2 > 0),
                "structure_difference": int(bool(n_exon or n_cds or n_start or n_end or n_nmd)),
                "domain_difference": int(n_level2 > 0),
            }
            gene_flags[gene] = flags
            structure_rows.append({
                "gene_id": gene, "n_matrix_transcripts": counts[gene], "n_alternative_transcripts": n_alt,
                "n_alternative_protein_coding": n_coding_alt, "n_alternative_nmd": n_nmd,
                "n_alternative_noncoding": n_noncoding, "n_exon_composition_difference": n_exon,
                "n_cds_length_difference": n_cds, "n_alternative_start": n_start,
                "n_alternative_termination": n_end,
            })
            sequence_diff_rows.append({
                "gene_id": gene,
                "principal_sequence_length": principal_sequence["length"] if principal_sequence else "",
                "principal_sequence_sha256": principal_sequence["sha256"] if principal_sequence else "",
                "n_alternative_with_sequence": n_with_sequence,
                "n_identical_sequence": n_identical,
                "n_different_sequence": n_different,
                "n_protein_coding_sequence_missing": n_missing,
                "n_level1_pairs": n_level1,
                "n_level2_pairs": n_level2,
            })
            functional_rows.append({
                "gene_id": gene,
                "principal_transcript_id": ptx,
                "principal_protein_id": principal_protein,
                "n_pfam_accessions": pfam_n,
                "signal_peptide": int(bool(principal_features and principal_features["signal"])) if principal_features else "",
                "n_transmembrane_spans": principal_features["tm"] if principal_features else "",
                "n_catalytic_features": principal_features["catalytic"] if principal_features else 0,
                "strongest_pair_evidence": "LEVEL_2" if n_level2 else ("LEVEL_1" if n_level1 else "LEVEL_0"),
                "reliable_domain_coordinates": reliable,
                "uniprot_reviewed_or_unreviewed": int(bool(swiss_acc or trembl_acc)),
            })
    write_csv(PHASE / "results" / "transcript_structure_coverage.csv", structure_rows, list(structure_rows[0]))
    write_csv(PHASE / "results" / "protein_sequence_difference_coverage.csv", sequence_diff_rows, list(sequence_diff_rows[0]))
    write_csv(PHASE / "results" / "functional_annotation_coverage.csv", functional_rows, list(functional_rows[0]))

    version_matched = 0
    principals_with_sequence = 0
    pfam_features_used = 0
    pfam_features_bad = 0
    for gene in candidate_genes:
        record = records[principal[gene]]
        ids = sorted(record.protein_ids)
        if len(ids) != 1 or ids[0] not in fasta:
            continue
        principals_with_sequence += 1
        mapped = map_protein(ensembl, ids[0], fasta[ids[0]])
        if mapped["status"] == "version_matched":
            version_matched += 1
            status = feature_status(mapped["features"], fasta[ids[0]]["length"])
            pfam_features_used += status["pfam_n"]
            pfam_features_bad += status["pfam_bad"]
    version_fraction = version_matched / principals_with_sequence if principals_with_sequence else 0
    coordinate_fraction = 1 - (pfam_features_bad / pfam_features_used) if pfam_features_used else 0

    pfam_analysis = next((item for item in ensembl["analysis"].values() if item["logic_name"] == "pfam"), {})
    signal_analysis = next((item for item in ensembl["analysis"].values() if item["logic_name"] == "signalp"), {})
    tm_analysis = next((item for item in ensembl["analysis"].values() if item["logic_name"] == "tmhmm"), {})
    interpro_present = any(item["logic_name"] == "interpro" for item in ensembl["analysis"].values())
    source_rows = [
        file_row(GTF_PATH, "GENCODE GTF", "v34", "already local", "phase0/raw/annotation/gencode.v34.annotation.gtf.gz", "frozen project copy", "transcript and CDS coordinates"),
        file_row(SWISSPROT, "GENCODE SwissProt metadata", "v34", "already local", "phase0/raw/annotation/gencode.v34.metadata.SwissProt.gz", "bundled with the project GENCODE v34 annotation", "reviewed UniProt accession crosswalk"),
        file_row(TREMBL, "GENCODE TrEMBL metadata", "v34", "2026-11-14", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.metadata.TrEMBL.gz", "Wed, 29 Apr 2020 15:45:23 GMT", "unreviewed UniProt accession crosswalk"),
        file_row(FASTA, "GENCODE protein translations", "v34", "2026-11-14", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.pc_translations.fa.gz", "Wed, 29 Apr 2020 15:45:23 GMT", "predicted protein sequence"),
        file_row(FEATURE, "Ensembl protein_feature", "release 100", "2026-11-14", "https://ftp.ensembl.org/pub/release-100/mysql/homo_sapiens_core_100_38/protein_feature.txt.gz", "Wed, 20 May 2020 10:53:36 GMT", "domain and feature coordinates"),
        file_row(TRANSLATION, "Ensembl translation", "release 100", "2026-11-14", "https://ftp.ensembl.org/pub/release-100/mysql/homo_sapiens_core_100_38/translation.txt.gz", "Wed, 20 May 2020 10:56:16 GMT", "ENSP stable id and version"),
        file_row(ANALYSIS, "Ensembl analysis", "release 100", "2026-11-14", "https://ftp.ensembl.org/pub/release-100/mysql/homo_sapiens_core_100_38/analysis.txt.gz", "Wed, 20 May 2020 10:43:51 GMT", "feature database name and version"),
        file_row(UNIPROT_COUNTS, "UniProt REST count", "2026_04", "2026-11-13", "https://rest.uniprot.org/uniprotkb/search", "Fri, 13 Nov 2026 00:13:56 GMT", "current reviewed-human count only; not joined"),
    ]
    source_rows.append({
        "source": "Pfam through Ensembl 100",
        "release": pfam_analysis.get("db_version", ""),
        "retrieval_date": "2026-11-14",
        "url_or_accession": "logic_name=pfam in the Ensembl 100 analysis table",
        "file_size_bytes": "",
        "sha256": "",
        "annotation_type": f"db={pfam_analysis.get('db', '')}; program={pfam_analysis.get('program', '')} {pfam_analysis.get('program_version', '')}",
        "species": "Homo sapiens",
        "identifier_versioning": "Pfam accession before the dot",
        "source_last_modified": "inherited from the Ensembl 100 analysis table",
        "used_for_domain_coordinates": 1,
    })
    source_rows.append({
        "source": "InterPro",
        "release": "not a separate file",
        "retrieval_date": "2026-11-14",
        "url_or_accession": "",
        "file_size_bytes": "",
        "sha256": "",
        "annotation_type": "present_in_ensembl_logic_names" if interpro_present else "no InterPro logic_name in the Ensembl 100 analysis table",
        "species": "Homo sapiens",
        "identifier_versioning": "",
        "source_last_modified": "",
        "used_for_domain_coordinates": 0,
    })
    write_csv(PHASE / "manifests" / "protein_domain_coverage.csv", source_rows, list(source_rows[0]))

    log("link phase1a")
    ok_genes = read_ok_genes()
    frozen_counts = expected_success()
    linkage_rows = []
    linkage_ok = True

    def count_scope(scope: str, cancer: str, transform: str, genes: list[str]) -> dict:
        rows = [gene_flags[gene] for gene in genes if gene in gene_flags]
        missing = len(genes) - len(rows)
        return {
            "scope": scope,
            "cancer": cancer,
            "transformation": transform,
            "n_genes": len(genes),
            "n_missing_from_candidate_annotation": missing,
            "n_unique_principal": sum(item["unique_principal"] for item in rows),
            "n_valid_alternative_protein_coding": sum(item["valid_alt_protein_coding"] for item in rows),
            "n_distinct_alternative_sequence": sum(item["distinct_alt_sequence"] for item in rows),
            "n_uniprot_mapping": sum(item["uniprot"] for item in rows),
            "n_pfam_annotation": sum(item["pfam"] for item in rows),
            "n_reliable_domain_coordinates": sum(item["reliable_coordinates"] for item in rows),
            "n_interpretable_sequence_difference": sum(item["interpretable_sequence_difference"] for item in rows),
        }

    linkage_rows.append(count_scope("candidate_universe", "ALL", "AB", candidate_genes))
    wide_rows = []
    for cancer in CANCERS:
        per_cancer = {"cancer": cancer}
        sets = {}
        for transform in ("A", "B"):
            genes = ok_genes.get((cancer, transform), set())
            expected_n, expected_candidate = frozen_counts["success"][(cancer, transform)]
            if len(genes) != expected_n or expected_candidate != 10795:
                linkage_ok = False
            row = count_scope("successful_fit", cancer, transform, sorted(genes))
            linkage_rows.append(row)
            sets[transform] = genes
            for key, value in row.items():
                if key not in {"scope", "cancer", "transformation"}:
                    per_cancer[f"{transform}_{key}"] = value
        both = sets["A"] & sets["B"]
        if len(both) != frozen_counts["paired"][cancer]:
            linkage_ok = False
        both_row = count_scope("intersection_A_and_B", cancer, "AB", sorted(both))
        linkage_rows.append(both_row)
        per_cancer["intersection_n_genes"] = both_row["n_genes"]
        per_cancer["intersection_n_distinct_alternative_sequence"] = both_row["n_distinct_alternative_sequence"]
        per_cancer["intersection_n_interpretable_sequence_difference"] = both_row["n_interpretable_sequence_difference"]
        per_cancer["intersection_n_reliable_domain_coordinates"] = both_row["n_reliable_domain_coordinates"]
        wide_rows.append(per_cancer)
    if any(row["n_missing_from_candidate_annotation"] for row in linkage_rows):
        linkage_ok = False
    write_csv(PHASE / "results" / "phase1a_annotation_coverage.csv", linkage_rows, list(linkage_rows[0]))
    write_csv(PHASE / "results" / "annotation_coverage_by_cancer.csv", wide_rows, list(wide_rows[0]))

    grades = {
        "A": "PASS" if universe_ok else "FAIL",
        "B": "PASS" if one_protein_fraction >= 0.99 and sequence_coverage >= 0.99 and n_multi_sequence == 0 else (
            "PARTIAL" if sequence_coverage >= 0.80 else "FAIL"
        ),
        "C": "PASS" if version_fraction >= 0.80 and (pfam_features_used == 0 or coordinate_fraction >= 0.95) else (
            "PARTIAL" if version_fraction >= 0.40 else "FAIL"
        ),
        "D": "PASS" if linkage_ok else "FAIL",
        "E": "PARTIAL",
        "F": "PARTIAL",
    }
    label = overall_label(grades)
    n_structure = sum(gene_flags[gene]["structure_difference"] for gene in candidate_genes)
    n_sequence = sum(gene_flags[gene]["distinct_alt_sequence"] for gene in candidate_genes)
    n_domain = sum(gene_flags[gene]["domain_difference"] for gene in candidate_genes)
    n_uniprot = sum(gene_flags[gene]["uniprot"] for gene in candidate_genes)
    n_pfam = sum(gene_flags[gene]["pfam"] for gene in candidate_genes)
    n_reliable = sum(gene_flags[gene]["reliable_coordinates"] for gene in candidate_genes)
    n_coding_alt_genes = sum(gene_flags[gene]["valid_alt_protein_coding"] for gene in candidate_genes)

    (PHASE / "reports" / "transcript_protein_mapping_report.md").write_text(
        "# Transcript-to-protein mapping\n\n"
        "The join is ENST to ENSG from the GENCODE v34 GTF, ENST to ENSP from the GTF protein identifier, "
        "ENSP to the GENCODE v34 translation, and ENST to SwissProt or TrEMBL accessions from the GENCODE v34 metadata.\n\n"
        f"Candidate measured protein-coding transcripts: {n_coding}. "
        f"Transcripts with one protein identifier: {n_coding_one_protein} ({one_protein_fraction:.4f}). "
        f"Those identifiers found in the translation FASTA: {n_coding_with_sequence} ({sequence_coverage:.4f}). "
        f"Transcripts with more than one protein identifier: {n_multi_sequence}. "
        f"FASTA duplicate protein identifiers: {fasta_duplicates}, of which {fasta_conflicts} had a different sequence. The first identical record was kept. "
        f"Header transcript conflicts: {n_header_conflict}.\n\n"
        f"Candidate genes whose principal transcript has a SwissProt or TrEMBL accession: {n_uniprot} of {len(candidate_genes)}.\n\n"
        "Identical translation strings are counted as the same predicted protein even when the ENSP identifiers differ. "
        "Different accessions are retained and are not treated as proof of a functional difference.\n\n"
        f"Mapping grade: {grades['B']}.\n"
    )
    catalytic_logic = any(
        item["logic_name"] in {"catalytic_residue", "active_site"}
        for item in ensembl["analysis"].values()
    )
    (PHASE / "reports" / "protein_domain_annotation_report.md").write_text(
        "# Protein domain annotation\n\n"
        "Domain coordinates come from Ensembl release 100, the gene build paired with GENCODE v34. "
        f"Pfam in that analysis table is database `{pfam_analysis.get('db', '')}` version `{pfam_analysis.get('db_version', '')}`. "
        f"SignalP version `{signal_analysis.get('db_version', '')}`. TMHMM version `{tm_analysis.get('db_version', '')}`.\n\n"
        f"Candidate principal proteins with a GENCODE translation: {principals_with_sequence}. "
        f"Version-matched Ensembl translations: {version_matched} ({version_fraction:.4f}). "
        f"Pfam features on those rows: {pfam_features_used}. Unsupported Pfam coordinates: {pfam_features_bad}.\n\n"
        f"Candidate genes with at least one supported Pfam accession on the principal protein: {n_pfam}. "
        f"Candidate genes with a version-matched coordinate map and supported feature coordinates: {n_reliable}.\n\n"
        "UniProt 2026_04 is recorded as a current count of reviewed human entries. It is not joined to these coordinates. "
        "InterPro has no domain-coordinate logic name in this analysis table. "
        "Pfam hits were produced by the Ensembl analysis program recorded in the coverage manifest. "
        f"Catalytic-residue logic names are {'present' if catalytic_logic else 'absent'}, so no catalytic-residue difference is called.\n\n"
        "A missing feature is an absent annotation. A predicted composition difference is counted only for version-matched proteins with supported coordinates.\n\n"
        f"Domain grade: {grades['C']}.\n"
        "\nPair classes among compared sequence pairs:\n\n"
        + "\n".join(f"- {key}: {value}" for key, value in sorted(pair_class_counts.items()))
        + "\n"
    )
    candidate_lookup = set(candidate_genes)
    outside_biotype = defaultdict(int)
    outside_genes = defaultdict(set)
    for tx_id, record in records.items():
        if record.gene not in candidate_lookup or tx_id in matrix_set:
            continue
        outside_biotype[record.biotype] += 1
        outside_genes[record.biotype].add(record.gene)
    (PHASE / "reports" / "structural_difference_feasibility.md").write_text(
        "# Structural difference feasibility\n\n"
        "Four evidence layers stay separate.\n\n"
        f"A. Transcript structure. {n_structure} of {len(candidate_genes)} candidate genes have a matrix alternative "
        "with an exon-composition difference, a CDS-length difference, an alternative CDS boundary, or an NMD biotype.\n\n"
        f"B. Predicted protein sequence. {n_sequence} candidate genes have at least one protein-coding matrix alternative "
        "whose GENCODE v34 translation string differs from the principal translation. "
        f"Genes with a protein-coding alternative in the matrix: {n_coding_alt_genes}.\n\n"
        f"C. Predicted functional-region annotation. {n_domain} candidate genes have a Level 2 pair: "
        "the sequences differ and a supported Pfam, signal-peptide, or transmembrane difference is present.\n\n"
        "D. Experimentally supported protein difference. The count is zero. No isoform-specific peptide sequence is in the frozen inputs.\n\n"
        "The sequence method is exact string equality, `phase2a_exact_v1`. "
        "Unequal strings have a blank alignment identity and a length-difference lower bound. "
        "That lower bound is not a residue alignment.\n\n"
        "NMD translations are kept at Level 0. A predicted ORF change is not counted as a stable protein product.\n\n"
        "The comparison uses transcripts present in the isoform matrix. "
        f"For the 10,795 candidate genes, GENCODE v34 also annotates {outside_biotype['nonsense_mediated_decay']} "
        f"NMD transcripts in {len(outside_genes['nonsense_mediated_decay'])} genes and "
        f"{outside_biotype['protein_coding']} protein-coding transcripts that are absent from the matrix. "
        "Those unmeasured transcripts were not compared.\n"
    )
    (PHASE / "reports" / "phase1a_annotation_linkage_report.md").write_text(
        "# Phase 1A annotation linkage\n\n"
        "The linkage uses every successful Phase 1A fit. Partial R² was not read. "
        "The OK counts were checked against `eligibility_by_cancer_transform.csv`, and the A-and-B intersection was checked against `n_paired`.\n\n"
        "Annotation does not change by cancer. Cancer-to-cancer differences in these tables are differences in which genes have a successful fit.\n\n"
        f"Linkage grade: {grades['D']}.\n\n"
        "The eight counts are in `phase2a/results/phase1a_annotation_coverage.csv`.\n"
    )
    criteria = [
        {"level": "LEVEL_0", "definition": "Transcript annotation, a missing translation, an NMD biotype, or an identical pair of translations.", "assigned_events": evidence_counts.get("LEVEL_0", 0) + evidence_counts.get("LEVEL_0_IDENTICAL_SEQUENCE", 0)},
        {"level": "LEVEL_1", "definition": "Two protein-coding GENCODE v34 translations differ, and no supported functional-region difference is called.", "assigned_events": evidence_counts.get("LEVEL_1", 0)},
        {"level": "LEVEL_2", "definition": "A Level 1 sequence difference plus a supported Pfam, signal-peptide, or transmembrane difference on version-matched coordinates.", "assigned_events": evidence_counts.get("LEVEL_2", 0)},
        {"level": "LEVEL_3", "definition": "An isoform-specific peptide or other protein identification distinguishes the isoform. Gene-level abundance does not qualify.", "assigned_events": 0},
    ]
    write_csv(PHASE / "manifests" / "evidence_level_criteria.csv", criteria, ["level", "definition", "assigned_events"])
    (PHASE / "reports" / "evidence_hierarchy_report.md").write_text(
        "# Evidence hierarchy\n\n"
        "Level 0 is transcript annotation. Level 1 is a predicted protein-sequence difference. "
        "Level 2 is a predicted functional-region difference on supported coordinates. "
        "Level 3 is independent isoform-specific protein evidence.\n\n"
        + "\n".join(f"- {row['level']}: {row['assigned_events']} pair events. {row['definition']}" for row in criteria)
        + f"\n\nOf the Level 0 events, {evidence_counts.get('LEVEL_0_IDENTICAL_SEQUENCE', 0)} are identical protein-coding translations and {evidence_counts.get('LEVEL_0', 0)} are matrix alternatives without a paired protein-coding translation.\n\n"
        "The current SEPEP files are sample-level quantitative groupings. Six accessible files were assessed in Phase 0.5, four returned HTTP 403, and none of the accessible files contain peptide sequences or protein accessions. "
        "This audit assigns no Level 3 events and does not promote gene-level CPTAC abundance to isoform evidence.\n"
    )
    write_csv(
        PHASE / "manifests" / "independent_transcript_validation_candidates.csv",
        validation_rows(),
        list(validation_rows()[0]),
    )
    (PHASE / "reports" / "independent_validation_feasibility.md").write_text(
        "# Independent validation feasibility\n\n"
        "The search covered GEO, EGA, and the GDC publication page for cohorts that could test transcript usage outside the 1,027 CPTAC cases. "
        "No cohort met the locked PASS rule: cancer and control groups, patient identifiers, transcript quantification, GENCODE v34-compatible annotation, at least 50 cases, and no CPTAC sample overlap.\n\n"
        "GSE218035 is a public osteosarcoma long-read RNA series with 23 patients and 13 adjacent samples. "
        "It is independent of the CPTAC case set and is too small for the locked threshold.\n\n"
        "EGAS00001004819 is the Veiga breast long-read resource. Its peptide support uses CPTAC spectra, so it is not an independent protein validation.\n\n"
        "GSE270257 contains long-read DNA methylation for 189 POG tumors. The GEO object does not contain transcript quantification. Raw reads are in EGA EGAS00001001159.\n\n"
        "The TCGA splicing supplement covers 8,705 patients and has partial sample overlap with TCGA-linked CPTAC cohorts. Its splicing graph is not the frozen GENCODE v34 principal rule.\n\n"
        "GSE293661 is a ten-cell-line pancreatic long-read series and has no patient control arm.\n\n"
        "An independent RNA cohort can test transcript usage. It cannot validate protein abundance. None of these cohorts is combined with the 1,027 cases, and no DTU test was run.\n\n"
        "Validation grade: PARTIAL.\n"
    )
    write_csv(
        PHASE / "manifests" / "literature_competition_audit.csv",
        literature_rows(),
        list(literature_rows()[0]),
    )
    (PHASE / "reports" / "novelty_and_competition_report.md").write_text(
        "# Novelty and competition\n\n"
        "The literature audit is a targeted reading of primary papers and one 2026 meeting abstract. It is not a claim that every related paper was found.\n\n"
        "Pan-cancer isoform switching and predicted domain change are already published for TCGA RNA-seq, including Climente-González et al. 2017 and Vitting-Seerup and Sandelin 2017. "
        "Kahles et al. 2018 quantify splicing events across 8,705 tumors and confirm a subset of neojunction peptides with CPTAC spectra. "
        "Mertins et al. 2016 connect somatic alterations to gene-level breast-cancer proteomes. "
        "Veiga et al. 2022 add long-read breast isoforms and CPTAC peptide matches. "
        "A 2026 AACR abstract describes peptide-to-isoform graphs in the ten CPTAC cohorts.\n\n"
        "Those studies leave the frozen estimand unestablished: the partial R² of a prespecified GENCODE v34 principal-transcript proportion beyond reconstructed gene expression, tumor purity, and CNV, estimated inside each cancer and then summarized across cancers. "
        "The audit does not treat that gap as a reason to change the Phase 1A result.\n\n"
        "Novelty grade: PARTIAL.\n"
    )
    reason = (
        f"Sequence coverage for one-protein candidate transcripts is {sequence_coverage:.4f}. "
        f"Version-matched principal proteins are {version_fraction:.4f} of principals with a translation. "
        f"The structural-study label is {label}. A functional-loss study remains unqualified because Level 3 evidence is absent "
        "and because a predicted coordinate difference is not a demonstrated protein product."
    )
    write_decision(grades, label, "COMPLETED", reason)
    env = (
        f"python {platform.python_version()}\n"
        f"implementation {platform.python_implementation()}\n"
        f"platform {platform.platform()}\n"
        "numpy imported through phase1a/code/common.py for the frozen principal parser\n"
        "sequence comparison phase2a_exact_v1\n"
        f"run_finished {utc_now()}\n"
    )
    (PHASE / "manifests" / "software_environment.txt").write_text(env)
    input_rows = source_rows + hash_rows
    write_csv(
        PHASE / "manifests" / "input_manifest.csv",
        [{"path": row.get("path", row.get("source", "")), "sha256": row.get("sha256", ""), "role": row.get("annotation_type", row.get("ledger_match", ""))} for row in input_rows],
        ["path", "sha256", "role"],
    )
    (PHASE / "qc" / "decision_grades.txt").write_text(
        "\n".join(f"{key} {value}" for key, value in grades.items()) + f"\nOVERALL {label}\n"
    )
    hash_outputs()
    log(f"PHASE2A_DONE {grades['A']} {grades['B']} {grades['C']} {grades['D']} {label}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        (PHASE / "logs").mkdir(parents=True, exist_ok=True)
        (PHASE / "logs" / "run_error.txt").write_text(traceback.format_exc())
        if not (PHASE / "reports" / "phase2a_final_decision.md").exists():
            write_decision(
                {"A": "FAIL", "B": "FAIL", "C": "FAIL", "D": "FAIL", "E": "FAIL", "F": "FAIL"},
                "STOP_CURRENT_DIRECTION",
                "FAILED",
                "The audit stopped on an execution error. See logs/run_error.txt.",
            )
            hash_outputs()
        raise
