"""Shared paths and readers for Phase 1A. Does not fit protein models."""
from __future__ import annotations

import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PHASE1A = ROOT / "phase1a"
CANCERS = ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"]
EXPECTED = {
    "BRCA": 121, "CCRCC": 103, "COAD": 96, "GBM": 99, "HNSCC": 108,
    "LSCC": 108, "LUAD": 110, "OV": 82, "PDAC": 105, "UCEC": 95,
}
GTF = ROOT / "phase0/raw/annotation/gencode.v34.annotation.gtf.gz"
GTF_SHA = "941926bde63e86f4eebbd29a8f747937288750da6b7dda8b21a920980c38b6e1"
REPAIRED_SHA = {
    ("BRCA", "rna_isoform"): "c9e071a8b393e75832673b3fd4d3a1754c4a7d25ac2ac0104c2e03f832f43ca9",
    ("UCEC", "rna_gene"): "3d708a346477a51f2aa1fc80479cd60b1888644145b69b1082c31c6405770882",
    ("PDAC", "cnv"): "12f2a560765c9a984251afc60201637eefbe063e09e77b8ac37427ef0da67510",
}
LAYER_SUFFIX = {
    "rna_gene": "RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "rna_isoform": "RNAseq_isoform_FPKM_log2_Tumor.txt",
    "protein": "proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "cnv": "WES_CNV_gene_ratio_log2.txt",
    "phenotype": "phenotype.txt",
}
NA_TOKENS = {b"", b"NA", b"NaN", b"nan", b"Na"}
LOG2E = float(np.log(2.0))


def acquisition_path(cancer: str, layer: str) -> Path:
    if layer == "phenotype":
        name = f"{cancer}_phenotype.txt"
    else:
        name = f"{cancer}_{LAYER_SUFFIX[layer]}"
    repaired = ROOT / "phase0_6" / "repaired" / cancer / name
    if (cancer, layer) in REPAIRED_SHA:
        return repaired
    return ROOT / "phase0_5" / "acquisition" / cancer / name


def historical_path(cancer: str, layer: str) -> Path:
    if layer == "phenotype":
        name = f"{cancer}_phenotype.txt"
    else:
        name = f"{cancer}_{LAYER_SUFFIX[layer]}"
    return ROOT / "phase0_5" / "acquisition" / cancer / name


def sha256_and_nul(path: Path) -> tuple[str, int, int, int]:
    digest = hashlib.sha256()
    size = nul = newlines = 0
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
            size += len(block)
            nul += block.count(b"\x00")
            newlines += block.count(b"\n")
    return digest.hexdigest(), size, nul, newlines


def load_ledger(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(maxsplit=1)
        out[name.strip()] = digest
    return out


def header_names(path: Path) -> list[str]:
    with path.open("rb") as handle:
        parts = handle.readline().rstrip(b"\r\n").split(b"\t")
    if not parts or parts[0] != b"idx":
        raise SystemExit(f"header does not start with idx: {path}")
    names = [part.decode() for part in parts[1:]]
    if any(name.startswith("._") for name in names):
        raise SystemExit(f"AppleDouble column in {path}")
    return names


def load_annotation() -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    tx_to_gene: dict[str, str] = {}
    bare_genes: dict[str, set[str]] = defaultdict(set)
    mane: dict[str, set[str]] = defaultdict(set)
    appris: dict[str, set[str]] = defaultdict(set)
    with gzip.open(GTF, "rb") as handle:
        for line in handle:
            if line.startswith(b"#"):
                continue
            parts = line.split(b"\t", 8)
            if len(parts) < 9 or parts[2] != b"transcript":
                continue
            attr = parts[8]
            tx = gene = None
            for field in attr.split(b";"):
                field = field.strip()
                if field.startswith(b"transcript_id "):
                    tx = field.split(b'"')[1].decode()
                elif field.startswith(b"gene_id "):
                    gene = field.split(b'"')[1].decode()
            if not tx or not gene:
                continue
            tx_to_gene[tx] = gene
            bare_genes[tx.split(".", 1)[0]].add(gene)
            if b'tag "MANE_Select"' in attr:
                mane[gene].add(tx)
            if b'tag "appris_principal_1"' in attr:
                appris[gene].add(tx)
    principal: dict[str, str] = {}
    rule: dict[str, str] = {}
    for gene, txs in mane.items():
        if len(txs) == 1:
            principal[gene] = next(iter(txs))
            rule[gene] = "MANE_Select"
    for gene, txs in appris.items():
        if gene in mane:
            continue
        if len(txs) == 1:
            principal[gene] = next(iter(txs))
            rule[gene] = "APPRIS_principal_1"
    return tx_to_gene, principal, rule


def parse_numeric_row(parts: list[bytes], columns: list[int | None]) -> np.ndarray:
    row = np.empty(len(columns), dtype=np.float64)
    width = len(parts)
    for j, col in enumerate(columns):
        if col is None or col + 1 >= width:
            row[j] = np.nan
            continue
        token = parts[col + 1].strip()
        if token in NA_TOKENS:
            row[j] = np.nan
        else:
            row[j] = float(token)
    return row


def load_gene_matrix(path: Path, case_ids: list[str], allow_missing: bool = False) -> dict[str, np.ndarray]:
    names = header_names(path)
    locate = {name: i for i, name in enumerate(names)}
    if len(locate) != len(names):
        raise SystemExit(f"duplicate sample ids in {path}")
    missing = [case for case in case_ids if case not in locate]
    if missing and not allow_missing:
        raise SystemExit(f"{path} is missing {len(missing)} triple cases, first {missing[:3]}")
    columns = [locate.get(case) for case in case_ids]
    out: dict[str, np.ndarray] = {}
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip() or b"\x00" in line:
                if b"\x00" in line:
                    raise SystemExit(f"NUL byte inside a data line of {path}")
                continue
            parts = line.rstrip(b"\r\n").split(b"\t")
            gene = parts[0].decode()
            if gene in out:
                raise SystemExit(f"duplicate feature row {gene} in {path}")
            out[gene] = parse_numeric_row(parts, columns)
    return out


def load_phenotype(path: Path, case_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
    with path.open("rb") as handle:
        header = [part.decode() for part in handle.readline().rstrip(b"\r\n").split(b"\t")]
        purity_i = header.index("WES_purity")
        estimate_i = header.index("ESTIMATE_ESTIMATEScore")
        table: dict[str, tuple[float, float]] = {}
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip(b"\r\n").split(b"\t")
            sample = parts[0].decode()
            if sample in table:
                raise SystemExit(f"duplicate phenotype id {sample} in {path}")

            def value(index: int) -> float:
                token = parts[index].strip() if index < len(parts) else b""
                if token in NA_TOKENS:
                    return np.nan
                return float(token)

            table[sample] = (value(purity_i), value(estimate_i))
    purity = np.array([table.get(case, (np.nan, np.nan))[0] for case in case_ids])
    estimate = np.array([table.get(case, (np.nan, np.nan))[1] for case in case_ids])
    return purity, estimate


def triple_from_headers(cancer: str) -> list[str]:
    gene = set(header_names(acquisition_path(cancer, "rna_gene")))
    isoform = set(header_names(acquisition_path(cancer, "rna_isoform")))
    protein = set(header_names(acquisition_path(cancer, "protein")))
    return sorted(gene & isoform & protein)


def phase0_triple_ids() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {cancer: set() for cancer in CANCERS}
    path = ROOT / "sample_pairing_audit.csv"
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["pairing_call"] == "case_level_gene_isoform_protein":
                out[row["cancer"]].add(row["sample_id"])
    return out


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
