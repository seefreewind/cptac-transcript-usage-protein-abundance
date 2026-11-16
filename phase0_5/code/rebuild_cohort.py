#!/usr/bin/env python3
"""Rebuild the case intersection and RNA–protein coverage from complete matrices.

Reads only files already validated into phase0_5/acquisition. Does not modify
Phase 0 files. Does not test transcript usage against protein abundance.
"""
from __future__ import annotations

import csv
import gzip
import math
import sys
from array import array
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ACQ = ROOT / "phase0_5" / "acquisition"
MAN = ROOT / "phase0_5" / "manifests"
REP = ROOT / "phase0_5" / "reports"
GTF = ROOT / "phase0/raw/annotation/gencode.v34.annotation.gtf.gz"
PHASE0_PAIR = ROOT / "sample_pairing_audit.csv"
CANCERS = ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"]
EXPECTED = {
    "BRCA": 121, "CCRCC": 103, "COAD": 96, "GBM": 99, "HNSCC": 108,
    "LSCC": 108, "LUAD": 110, "OV": 82, "PDAC": 105, "UCEC": 95,
}
LAYERS = {
    "rna_gene": "RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
    "rna_isoform": "RNAseq_isoform_FPKM_log2_Tumor.txt",
    "protein": "proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "phospho": "phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "cnv": "WES_CNV_gene_ratio_log2.txt",
    "meta": "_meta.txt",
    "phenotype": "_phenotype.txt",
}
LAMBDA = (1.959964 + 0.841621) ** 2
NUMERIC = {"negative", "zero", "positive"}


def path_for(cancer: str, layer: str) -> Path | None:
    name = f"{cancer}{LAYERS[layer]}" if layer in {"meta", "phenotype"} else f"{cancer}_{LAYERS[layer]}"
    p = ACQ / cancer / name
    if not p.exists() or p.name.startswith("._"):
        return None
    return p


def header_ids(path: Path) -> list[str]:
    with path.open("rb") as handle:
        line = handle.readline()
    parts = line.rstrip(b"\r\n").split(b"\t")
    if parts[0] != b"idx":
        raise SystemExit(f"unexpected header in {path}")
    return [part.decode() for part in parts[1:]]


def row_ids(path: Path, skip_descriptor: bool) -> list[str]:
    ids = []
    with path.open("rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(b"\t")
        names = [part.decode() for part in header]
        key = 0
        for i, name in enumerate(names):
            if name.lower() in {"case_id", "case", "sample", "sample_id"}:
                key = i
                break
        first = True
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip(b"\r\n").split(b"\t")
            if first and skip_descriptor:
                first = False
                token = cols[key] if key < len(cols) else b""
                if token in (b"data_type", b"") or (len(cols) > 1 and cols[1] == b"data_type"):
                    continue
            elif first:
                first = False
            ids.append(cols[key].decode() if key < len(cols) else "")
    return ids


def phase0_sets() -> dict[str, dict[str, set[str]]]:
    out = {cancer: defaultdict(set) for cancer in CANCERS}
    with PHASE0_PAIR.open(newline="") as handle:
        for row in csv.DictReader(handle):
            cancer = row["cancer"]
            if cancer in out:
                out[cancer][row["pairing_call"]].add(row["sample_id"])
    return out


def classify_token(raw: bytes) -> str:
    token = raw.strip()
    if token in (b"", b"NA", b"NaN", b"nan", b"Na"):
        return "na"
    try:
        value = float(token)
    except ValueError:
        return "invalid"
    if math.isinf(value):
        return "infinite"
    if value < 0:
        return "negative"
    if value == 0.0:
        return "zero"
    return "positive"


def bare(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def unique_bare(identifiers: set[str]) -> dict[str, str]:
    grouped: dict[str, set[str]] = defaultdict(set)
    for identifier in identifiers:
        grouped[bare(identifier)].add(identifier)
    return {key: next(iter(values)) for key, values in grouped.items() if len(values) == 1}


def duplicates(values: list[str]) -> list[str]:
    counts = Counter(values)
    return sorted(key for key, count in counts.items() if count > 1)


def partial_r2(n: float | str) -> str:
    if isinstance(n, str) or n <= 0:
        return ""
    f2 = LAMBDA / float(n)
    return f"{f2 / (1.0 + f2):.4f}"


def mean(values: list[int]) -> float | str:
    if not values:
        return ""
    return round(sum(values) / len(values), 4)


def median(values: list[int]) -> float | str:
    if not values:
        return ""
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def load_annotation() -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, int]]:
    """Map GENCODE v34 transcripts to genes and choose one principal transcript per gene.

    MANE Select wins over APPRIS principal-1. A bare transcript id is stored only
    when it points at one gene.
    """
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
            bare_genes[bare(tx)].add(gene)
            if b'tag "MANE_Select"' in attr:
                mane[gene].add(tx)
            if b'tag "appris_principal_1"' in attr:
                appris[gene].add(tx)
    principal: dict[str, str] = {}
    ambiguous = 0
    for gene, txs in mane.items():
        if len(txs) == 1:
            principal[gene] = next(iter(txs))
        else:
            ambiguous += 1
    appris_only = 0
    for gene, txs in appris.items():
        if gene in mane:
            continue
        if len(txs) == 1:
            principal[gene] = next(iter(txs))
            appris_only += 1
        else:
            ambiguous += 1
    tx_bare = {key: next(iter(genes)) for key, genes in bare_genes.items() if len(genes) == 1}
    counts = {
        "transcripts": len(tx_to_gene),
        "genes_with_principal": len(principal),
        "genes_with_mane": sum(len(txs) == 1 for txs in mane.values()),
        "genes_with_appris_only": appris_only,
        "genes_with_ambiguous_principal": ambiguous,
        "bare_transcripts_with_one_gene": len(tx_bare),
    }
    return tx_to_gene, tx_bare, principal, counts


def phenotype_table(path: Path) -> tuple[dict[str, bool], dict[str, bool], list[str]]:
    with path.open("rb") as handle:
        header = [part.decode() for part in handle.readline().rstrip(b"\r\n").split(b"\t")]
        if "WES_purity" not in header or "ESTIMATE_ESTIMATEScore" not in header:
            raise SystemExit(f"purity columns missing in {path}")
        purity_i = header.index("WES_purity")
        estimate_i = header.index("ESTIMATE_ESTIMATEScore")
        batch_columns = [name for name in header if "batch" in name.lower() or "plex" in name.lower()]
        purity: dict[str, bool] = {}
        estimate: dict[str, bool] = {}
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip(b"\r\n").split(b"\t")
            sample = cols[0].decode()
            purity[sample] = classify_token(cols[purity_i] if purity_i < len(cols) else b"") in NUMERIC
            estimate[sample] = classify_token(cols[estimate_i] if estimate_i < len(cols) else b"") in NUMERIC
    return purity, estimate, batch_columns


def scan(path: Path, case_order: list[str], store_bits: bool) -> dict:
    with path.open("rb") as handle:
        names = [part.decode() for part in handle.readline().rstrip(b"\r\n").split(b"\t")[1:]]
        dup_samples = duplicates(names)
        position = {name: i for i, name in enumerate(names)}
        columns = [position.get(sample) for sample in case_order]
        counts: Counter[str] = Counter()
        feature_ids: list[str] = []
        seen: set[str] = set()
        dup_features: list[str] = []
        bits: dict[str, int] = {}
        value_min = None
        value_max = None
        min_positive = None
        hundredths = 0
        numeric = 0
        open_unit = 0
        nonfinite: list[str] = []
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip(b"\r\n").split(b"\t")
            feature = cols[0].decode()
            if feature in seen:
                dup_features.append(feature)
            seen.add(feature)
            feature_ids.append(feature)
            mask = 0
            for bit, col in enumerate(columns):
                if col is None:
                    continue
                raw = cols[col + 1] if col + 1 < len(cols) else b""
                kind = classify_token(raw)
                counts[kind] += 1
                if kind == "infinite" and len(nonfinite) < 3 and raw.strip().decode() not in nonfinite:
                    nonfinite.append(raw.strip().decode())
                if kind not in NUMERIC:
                    continue
                value = float(raw.strip())
                numeric += 1
                if value_min is None or value < value_min:
                    value_min = value
                if value_max is None or value > value_max:
                    value_max = value
                if value > 0 and (min_positive is None or value < min_positive):
                    min_positive = value
                if 0 < value < 1:
                    open_unit += 1
                if abs(value * 100 - round(value * 100)) < 1e-6:
                    hundredths += 1
                mask |= 1 << bit
            if store_bits and feature not in bits:
                bits[feature] = mask
    cells = sum(counts.values())
    return {
        "ids": feature_ids,
        "n_features": len(feature_ids),
        "n_columns_present": sum(col is not None for col in columns),
        "counts": counts,
        "cells": cells,
        "bits": bits,
        "duplicate_samples": dup_samples,
        "duplicate_features": dup_features[:20],
        "n_duplicate_features": len(dup_features),
        "min": value_min,
        "max": value_max,
        "min_positive": min_positive,
        "open_unit": open_unit,
        "numeric": numeric,
        "fraction_on_hundredths": (hundredths / numeric) if numeric else None,
        "nonfinite_tokens": nonfinite,
    }


def miss_row(cancer: str, layer: str, scanned: dict) -> dict:
    cells = scanned["cells"] or 1
    counts = scanned["counts"]
    return {
        "cancer": cancer,
        "layer": layer,
        "features": scanned["n_features"],
        "cases_in_triple_present_in_this_matrix": scanned["n_columns_present"],
        "cells": scanned["cells"],
        "na_fraction": round(counts["na"] / cells, 6),
        "zero_fraction": round(counts["zero"] / cells, 6),
        "negative_fraction": round(counts["negative"] / cells, 6),
        "positive_fraction": round(counts["positive"] / cells, 6),
        "open_unit_interval_fraction_of_numeric": round(scanned["open_unit"] / scanned["numeric"], 6) if scanned["numeric"] else "",
        "infinite_fraction": round(counts["infinite"] / cells, 6),
        "invalid_fraction": round(counts["invalid"] / cells, 6),
        "min_numeric": scanned["min"],
        "min_positive": scanned["min_positive"],
        "max_numeric": scanned["max"],
        "fraction_of_numeric_values_on_0.01_grid": "" if scanned["fraction_on_hundredths"] is None else round(scanned["fraction_on_hundredths"], 6),
        "duplicate_feature_rows": scanned["n_duplicate_features"],
        "nonfinite_tokens": ";".join(scanned["nonfinite_tokens"]),
    }


def lookup_gene(tx: str, tx_to_gene: dict[str, str], tx_bare: dict[str, str]) -> str | None:
    gene = tx_to_gene.get(tx)
    if gene:
        return gene
    return tx_bare.get(bare(tx))


def usage_census(path: Path, case_order: list[str], tx_to_gene: dict[str, str], tx_bare: dict[str, str], principal: dict[str, str]) -> dict:
    """Recover both candidate linear scales and count zero denominators.

    Candidate B is X = 0 when y = 0 and X = 2^y - 1 otherwise.
    Candidate C is X = 0 when y = 0 and X = 2^y otherwise.
    Neither proportion is tested against protein.
    """
    principal_tx = {tx: gene for gene, tx in principal.items()}

    with path.open("rb") as handle:
        names = [part.decode() for part in handle.readline().rstrip(b"\r\n").split(b"\t")[1:]]
        position = {name: i for i, name in enumerate(names)}
        columns = [position.get(sample) for sample in case_order]
        n = len(case_order)
        gene_index: dict[str, int] = {}
        ntx: list[int] = []
        nobs: list[array] = []
        sum_b: list[array] = []
        sum_c: list[array] = []
        pobs: list[array] = []
        psum_b: list[array] = []
        psum_c: list[array] = []
        b_bad: list[array] = []
        mapped = 0
        unmapped = 0
        negatives = 0
        rows = 0
        seen_tx: set[str] = set()
        duplicate_tx = 0

        def slot(gene: str) -> int:
            gi = gene_index.get(gene)
            if gi is not None:
                return gi
            gi = len(ntx)
            gene_index[gene] = gi
            ntx.append(0)
            nobs.append(array("H", [0]) * n)
            sum_b.append(array("d", [0.0]) * n)
            sum_c.append(array("d", [0.0]) * n)
            pobs.append(array("H", [0]) * n)
            psum_b.append(array("d", [0.0]) * n)
            psum_c.append(array("d", [0.0]) * n)
            b_bad.append(array("B", [0]) * n)
            return gi

        for line in handle:
            if not line.strip():
                continue
            rows += 1
            cols = line.rstrip(b"\r\n").split(b"\t")
            tx = cols[0].decode()
            if tx in seen_tx:
                duplicate_tx += 1
            seen_tx.add(tx)
            gene = lookup_gene(tx, tx_to_gene, tx_bare)
            if not gene:
                unmapped += 1
                continue
            mapped += 1
            gi = slot(gene)
            ntx[gi] += 1
            chosen = principal.get(gene)
            is_principal = False
            if chosen is not None:
                if tx == chosen:
                    is_principal = True
                elif tx not in tx_to_gene and bare(tx) == bare(chosen):
                    is_principal = True
            for i, col in enumerate(columns):
                if col is None:
                    continue
                raw = cols[col + 1] if col + 1 < len(cols) else b""
                kind = classify_token(raw)
                if kind not in NUMERIC:
                    continue
                y = float(raw.strip())
                if not math.isfinite(y):
                    continue
                if y < 0:
                    negatives += 1
                    xb = math.pow(2.0, y)
                    xc = xb
                    b_bad[gi][i] = 1
                elif y == 0.0:
                    xb = 0.0
                    xc = 0.0
                else:
                    pow2 = math.pow(2.0, y)
                    if not math.isfinite(pow2):
                        continue
                    xb = pow2 - 1.0
                    xc = pow2
                nobs[gi][i] += 1
                sum_b[gi][i] += xb
                sum_c[gi][i] += xc
                if is_principal:
                    pobs[gi][i] += 1
                    psum_b[gi][i] += xb
                    psum_c[gi][i] += xc
            if rows % 40000 == 0:
                print(f"  isoform rows {rows}", flush=True)

    by_gene = {}
    g0 = gpos = all_na = b_blocked = p_defined = p_diff = 0
    genes_ge2 = 0
    for gene, gi in gene_index.items():
        n_rows = ntx[gi]
        if n_rows >= 2:
            genes_ge2 += 1
        a_na = a_g0 = a_gpos = a_block = a_def = a_diff = 0
        usable = 0
        for i in range(n):
            if nobs[gi][i] == 0:
                a_na += 1
                continue
            if b_bad[gi][i]:
                a_block += 1
                continue
            if sum_b[gi][i] <= 0.0 or sum_c[gi][i] <= 0.0:
                a_g0 += 1
                continue
            a_gpos += 1
            if n_rows < 2 or pobs[gi][i] == 0:
                continue
            a_def += 1
            pb = psum_b[gi][i] / sum_b[gi][i]
            pc = psum_c[gi][i] / sum_c[gi][i]
            if abs(pb - pc) > 0.05:
                a_diff += 1
            usable |= 1 << i
        by_gene[gene] = {
            "n_tx": n_rows,
            "all_na": a_na,
            "g0": a_g0,
            "gpos": a_gpos,
            "b_blocked": a_block,
            "p_defined": a_def,
            "p_diff": a_diff,
            "usable": usable,
        }
        all_na += a_na
        g0 += a_g0
        gpos += a_gpos
        b_blocked += a_block
        p_defined += a_def
        p_diff += a_diff
    return {
        "by_gene": by_gene,
        "mapped_rows": mapped,
        "unmapped_rows": unmapped,
        "matrix_rows": rows,
        "negative_cells": negatives,
        "genes": len(gene_index),
        "genes_ge2": genes_ge2,
        "gene_cases_all_na": all_na,
        "gene_cases_linear_zero": g0,
        "gene_cases_linear_positive": gpos,
        "gene_cases_inverse_b_blocked": b_blocked,
        "gene_cases_principal_defined": p_defined,
        "gene_cases_principal_fraction_differs_by_more_than_0.05": p_diff,
        "duplicate_transcript_ids": duplicate_tx,
    }


def mask_from(bits: dict[str, int], feature: str, feature_bare: dict[str, str]) -> tuple[int, str]:
    if feature in bits:
        return bits[feature], "versioned"
    stripped = feature_bare.get(bare(feature))
    if stripped and stripped in bits:
        return bits[stripped], "unversioned"
    return 0, "absent"


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    MAN.mkdir(parents=True, exist_ok=True)
    REP.mkdir(parents=True, exist_ok=True)
    print("loading GENCODE v34", flush=True)
    tx_to_gene, tx_bare, principal, ann_counts = load_annotation()
    print(f"annotation {ann_counts}", flush=True)
    old = phase0_sets()
    pair_rows: list[dict] = []
    inter_rows: list[dict] = []
    miss_rows: list[dict] = []
    cov_rows: list[dict] = []
    notes: list[str] = []
    per_gene_path = MAN / "per_gene_complete_observations.csv"
    per_gene_handle = per_gene_path.open("w", newline="")
    per_gene_writer = None
    per_gene_n = 0

    for cancer in CANCERS:
        files = {layer: path_for(cancer, layer) for layer in LAYERS}
        missing_files = [layer for layer, path in files.items() if path is None]
        if not files["rna_gene"] or not files["rna_isoform"] or not files["protein"]:
            notes.append(
                f"{cancer}: gene RNA, isoform RNA, and protein are not all validated "
                f"({', '.join(missing_files) or 'unread'}). Intersection not recomputed."
            )
            continue
        print(f"{cancer} headers", flush=True)
        headers = {
            "rna_gene": header_ids(files["rna_gene"]),
            "rna_isoform": header_ids(files["rna_isoform"]),
            "protein": header_ids(files["protein"]),
        }
        if files["phospho"]:
            headers["phospho"] = header_ids(files["phospho"])
        if files["cnv"]:
            headers["cnv"] = header_ids(files["cnv"])
        if files["phenotype"]:
            headers["phenotype"] = row_ids(files["phenotype"], skip_descriptor=False)
        if files["meta"]:
            headers["meta"] = row_ids(files["meta"], skip_descriptor=True)
        for layer, ids in headers.items():
            dup = duplicates(ids)
            if dup:
                notes.append(f"{cancer}: duplicate {layer} identifiers {dup[:8]}. Those columns were not used to force a match.")
        gene_ids = headers["rna_gene"]
        iso_ids = headers["rna_isoform"]
        prot_ids = headers["protein"]
        gene_set, iso_set, prot_set = set(gene_ids), set(iso_ids), set(prot_ids)
        triple = gene_set & iso_set & prot_set
        case_order = sorted(triple)
        old_triple = old[cancer].get("case_level_gene_isoform_protein", set())
        gained = sorted(triple - old_triple)
        lost = sorted(old_triple - triple)
        if gained or lost or len(triple) != EXPECTED[cancer]:
            notes.append(
                f"{cancer}: recomputed triple {len(triple)}; Phase 0 target {EXPECTED[cancer]}; "
                f"only in new headers {gained or 'none'}; missing from new headers {lost or 'none'}."
            )
        else:
            notes.append(f"{cancer}: recomputed triple {len(triple)} matches the Phase 0 target. No ID was added or removed.")
        universe = set()
        for key in ("rna_gene", "rna_isoform", "protein", "phospho", "cnv", "phenotype", "meta"):
            universe |= set(headers.get(key, []))
        phospho_set = set(headers.get("phospho", []))
        cnv_set = set(headers.get("cnv", []))
        phenotype_set = set(headers.get("phenotype", []))
        meta_set = set(headers.get("meta", []))
        for sample in sorted(universe):
            pair_rows.append({
                "cancer": cancer,
                "sample_id": sample,
                "in_rna_gene_tumor": int(sample in gene_set),
                "in_rna_isoform_tumor": int(sample in iso_set),
                "in_protein_tumor": int(sample in prot_set),
                "in_phospho_tumor": int(sample in phospho_set),
                "in_cnv": int(sample in cnv_set),
                "in_phenotype": int(sample in phenotype_set),
                "in_clinical_meta": int(sample in meta_set),
                "in_phase0_triple": int(sample in old_triple),
                "pairing_call": "case_level_gene_isoform_protein" if sample in triple else "not_in_triple",
                "identifier_level": "case_submitter_id",
                "same_aliquot_barcode_in_rna_and_protein": "no_aliquot_barcode_in_matrix_header",
            })

        print(f"{cancer} scanning gene RNA and protein", flush=True)
        gene_scan = scan(files["rna_gene"], case_order, store_bits=True)
        prot_scan = scan(files["protein"], case_order, store_bits=True)
        print(f"{cancer} scanning isoform values and linear denominators", flush=True)
        iso_miss = scan(files["rna_isoform"], case_order, store_bits=False)
        usage = usage_census(files["rna_isoform"], case_order, tx_to_gene, tx_bare, principal)
        purity_mask = 0
        estimate_mask = 0
        batch_columns: list[str] = []
        if files["phenotype"]:
            purity, estimate, batch_columns = phenotype_table(files["phenotype"])
            for i, sample in enumerate(case_order):
                if purity.get(sample):
                    purity_mask |= 1 << i
                if estimate.get(sample):
                    estimate_mask |= 1 << i
        cnv_scan = None
        if files["cnv"]:
            print(f"{cancer} scanning CNV", flush=True)
            cnv_scan = scan(files["cnv"], case_order, store_bits=True)
        if files["phospho"]:
            print(f"{cancer} scanning phosphoproteome", flush=True)
            phospho_scan = scan(files["phospho"], case_order, store_bits=False)
            miss_rows.append(miss_row(cancer, "phospho", phospho_scan))

        rna_ids = set(gene_scan["ids"])
        protein_ids = set(prot_scan["ids"])
        shared = rna_ids & protein_ids
        rna_bare_all = {bare(x) for x in rna_ids}
        protein_bare_all = {bare(x) for x in protein_ids}
        shared_bare = rna_bare_all & protein_bare_all
        par_rows = sorted(x for x in rna_ids if "_PAR_" in x)
        usage_bare = unique_bare(set(usage["by_gene"]))
        cnv_bits = cnv_scan["bits"] if cnv_scan else {}
        cnv_bare = unique_bare(set(cnv_bits))
        both_counts = []
        eligible_counts = []
        eligible_ge71 = 0
        eligible_without_cnv_ge71 = 0
        versioned_isoform = 0
        bare_isoform_only = 0
        no_isoform_gene = 0
        for gid in sorted(shared):
            rna_mask, _ = mask_from(gene_scan["bits"], gid, {})
            prot_mask, _ = mask_from(prot_scan["bits"], gid, {})
            both = rna_mask & prot_mask
            both_n = both.bit_count()
            both_counts.append(both_n)
            usage_gene = gid if gid in usage["by_gene"] else None
            isoform_match = "versioned" if usage_gene else "absent"
            if usage_gene is None:
                stripped = usage_bare.get(bare(gid))
                if stripped and stripped in usage["by_gene"]:
                    usage_gene = stripped
                    isoform_match = "unversioned"
                    bare_isoform_only += 1
                else:
                    no_isoform_gene += 1
            else:
                versioned_isoform += 1
            info = usage["by_gene"].get(usage_gene or "", {})
            usable = info.get("usable", 0)
            cnv_mask, cnv_match = mask_from(cnv_bits, gid, cnv_bare) if cnv_bits else (0, "cnv_file_absent" if cnv_scan is None else "absent")
            with_cnv = both & purity_mask & usable & (cnv_mask if cnv_bits else 0)
            without_cnv = both & purity_mask & usable
            elig_n = with_cnv.bit_count() if cnv_bits else ""
            if cnv_bits and elig_n >= 71 and info.get("n_tx", 0) >= 2:
                eligible_ge71 += 1
            if without_cnv.bit_count() >= 71 and info.get("n_tx", 0) >= 2:
                eligible_without_cnv_ge71 += 1
            if isinstance(elig_n, int) and info.get("n_tx", 0) >= 2 and usage_gene in principal:
                eligible_counts.append(elig_n)
            row = {
                "cancer": cancer,
                "gene_id": gid,
                "isoform_gene_match": isoform_match,
                "n_cases_rna_and_protein_numeric": both_n,
                "n_isoform_rows": info.get("n_tx", 0),
                "n_cases_all_isoform_rows_missing": info.get("all_na", ""),
                "n_cases_linear_sum_zero": info.get("g0", ""),
                "n_cases_linear_sum_positive": info.get("gpos", ""),
                "n_cases_principal_fraction_defined": info.get("p_defined", ""),
                "n_cases_primary_eligible_with_purity_and_cnv": elig_n,
                "cnv_id_match": cnv_match,
            }
            if per_gene_writer is None:
                per_gene_writer = csv.DictWriter(per_gene_handle, fieldnames=list(row.keys()))
                per_gene_writer.writeheader()
            per_gene_writer.writerow(row)
            per_gene_n += 1

        complete_ge71 = sum(n >= 71 for n in both_counts)
        inter_rows.append({
            "cancer": cancer,
            "triple_cases": len(triple),
            "phase0_target": EXPECTED[cancer],
            "id_delta_gained": ";".join(gained),
            "id_delta_lost": ";".join(lost),
            "rna_gene_columns": len(gene_ids),
            "isoform_columns": len(iso_ids),
            "protein_columns": len(prot_ids),
            "rna_genes": len(rna_ids),
            "protein_genes": len(protein_ids),
            "shared_genes_versioned": len(shared),
            "rna_not_in_protein": len(rna_ids - protein_ids),
            "protein_not_in_rna": len(protein_ids - rna_ids),
            "shared_genes_unversioned": len(shared_bare),
            "protein_ids_unmatched_even_without_version": len(protein_bare_all - rna_bare_all),
            "rna_rows_with_par_suffix": len(par_rows),
            "isoform_rows": usage["matrix_rows"],
            "isoform_rows_in_gencode_v34": usage["mapped_rows"],
            "isoform_rows_not_in_gencode_v34": usage["unmapped_rows"],
            "genes_with_at_least_one_isoform_row": usage["genes"],
            "genes_with_at_least_two_isoform_rows": usage["genes_ge2"],
            "shared_genes_with_versioned_isoform_gene": versioned_isoform,
            "shared_genes_with_only_unversioned_isoform_gene": bare_isoform_only,
            "shared_genes_with_no_isoform_gene": no_isoform_gene,
            "mean_cases_rna_and_protein_numeric": mean(both_counts),
            "median_cases_rna_and_protein_numeric": median(both_counts),
            "genes_with_at_least_71_rna_and_protein_observations": complete_ge71,
            "genes_with_at_least_71_primary_eligible_cases": eligible_ge71 if cnv_bits else "",
            "genes_with_at_least_71_eligible_before_cnv": eligible_without_cnv_ge71,
            "mean_primary_eligible_cases": mean(eligible_counts),
            "median_primary_eligible_cases": median(eligible_counts) if eligible_counts else "",
            "detectable_partial_r2_at_median_eligible_n": partial_r2(median(eligible_counts) if eligible_counts else ""),
            "gene_cases_all_transcripts_missing": usage["gene_cases_all_na"],
            "gene_cases_linear_denominator_zero": usage["gene_cases_linear_zero"],
            "gene_cases_linear_denominator_positive": usage["gene_cases_linear_positive"],
            "gene_cases_negative_value_blocks_log2_plus_1": usage["gene_cases_inverse_b_blocked"],
            "gene_cases_principal_fraction_defined": usage["gene_cases_principal_defined"],
            "gene_cases_where_two_inverses_change_principal_fraction_by_more_than_0.05": usage["gene_cases_principal_fraction_differs_by_more_than_0.05"],
            "isoform_negative_cells_in_triple": usage["negative_cells"],
            "triple_cases_with_numeric_wes_purity": purity_mask.bit_count(),
            "triple_cases_with_numeric_estimate_score": estimate_mask.bit_count(),
            "phenotype_batch_or_plex_columns": ";".join(batch_columns),
        })
        miss_rows.append(miss_row(cancer, "rna_gene", gene_scan))
        miss_rows.append(miss_row(cancer, "rna_isoform", iso_miss))
        miss_rows.append(miss_row(cancer, "protein", prot_scan))
        if cnv_scan:
            miss_rows.append(miss_row(cancer, "cnv", cnv_scan))
        cov_rows.append({
            "cancer": cancer,
            "triple_cases": len(triple),
            "wes_purity_numeric": purity_mask.bit_count(),
            "estimate_score_numeric": estimate_mask.bit_count(),
            "cnv_header_overlap": len(triple & cnv_set),
            "batch_or_plex_columns": ";".join(batch_columns),
        })
        iso_counts = iso_miss["counts"]
        notes.append(
            f"{cancer}: isoform numeric min {iso_miss['min']}, min positive {iso_miss['min_positive']}, "
            f"max {iso_miss['max']}, negative cells {iso_counts['negative']}, "
            f"exact zeros {iso_counts['zero']}, values in (0, 1) {iso_miss['open_unit']} of {iso_miss['numeric']} numeric cells, "
            f"fraction on a 0.01 grid {iso_miss['fraction_on_hundredths']}, "
            f"duplicate transcript ids {usage['duplicate_transcript_ids']}."
        )
        if iso_counts["negative"] == 0 and iso_miss["open_unit"] > 0:
            notes.append(
                f"{cancer}: an unclipped log2(FPKM) is incompatible with stored values in (0, 1) and no negative cells. "
                "log2(FPKM+1) and a clipped or rounded log2(FPKM) both remain possible. Exact zeros are linear zeros under both, "
                "so the zero-denominator count does not choose the inverse."
            )
        elif iso_counts["negative"] > 0:
            notes.append(f"{cancer}: negative isoform values are present, so log2(FPKM+1) is not the transform.")
        if par_rows:
            notes.append(
                f"{cancer}: the gene-RNA matrix has {len(par_rows)} additional rows whose identifiers end in _PAR_. "
                "The protein matrix uses the unsuffixed identifier. The versioned intersection already contains that unsuffixed identifier. "
                f"Examples: {', '.join(par_rows[:4])}."
            )
        if files["cnv"]:
            missing_cnv = sorted(triple - cnv_set)
            if missing_cnv:
                notes.append(f"{cancer}: triple cases absent from the CNV header: {', '.join(missing_cnv)}. They remain in the case-level match and are dropped only from models that require CNV.")
        if files["phenotype"]:
            missing_purity = [sample for i, sample in enumerate(case_order) if ((purity_mask >> i) & 1) == 0]
            if missing_purity:
                notes.append(f"{cancer}: triple cases with WES_purity missing: {', '.join(missing_purity)}. They are dropped from M0 and M1. ESTIMATE is recorded separately.")
        if files["phospho"] and phospho_scan["counts"]["infinite"]:
            notes.append(
                f"{cancer}: phosphosite matrix contains non-finite tokens {', '.join(phospho_scan['nonfinite_tokens'])}. "
                "Those cells are not recoded. The primary outcome is the global protein matrix."
            )
        del usage, gene_scan, prot_scan, iso_miss

    per_gene_handle.close()
    write_csv(MAN / "reconstructed_pairing.csv", pair_rows)
    write_csv(MAN / "gene_protein_intersection.csv", inter_rows)
    write_csv(MAN / "missingness_report.csv", miss_rows)
    write_csv(MAN / "covariate_availability.csv", cov_rows)
    lines = [
        "# Cohort integrity",
        "",
        "Recomputed from complete matrix headers and full matrices in `phase0_5/acquisition`. Phase 0 files were not edited.",
        "",
        f"GENCODE v34 transcript rows indexed: {ann_counts['transcripts']}. Genes with exactly one MANE Select transcript: {ann_counts['genes_with_mane']}. Genes with no MANE transcript and exactly one APPRIS principal-1 transcript: {ann_counts['genes_with_appris_only']}. Genes with no unique principal transcript: {ann_counts['genes_with_ambiguous_principal']}. A gene in that last set is not assigned a transcript by file order.",
        "",
        "The two linear inverses are applied only to describe denominators and the stability of the principal-transcript fraction. They are not used to rank genes or to fit a protein model.",
        "",
        "Candidate B, the inverse of log2(FPKM+1): stored 0 stays 0, and any other finite y becomes 2^y − 1.",
        "",
        "Candidate C, the inverse of a clipped or rounded log2(FPKM): stored 0 stays 0, and y > 0 becomes 2^y. A negative y blocks Candidate B for that gene-case and is inverted as 2^y under Candidate C.",
        "",
        "A gene-case has a zero linear denominator when every observed transcript is stored as 0. Those cells were not recoded as NA.",
        "",
    ]
    if not inter_rows:
        lines.append("No cancer yet has gene RNA, isoform RNA, and protein on disk together. The intersection was not invented.")
    lines.extend(f"- {note}" for note in notes)
    lines.append("")
    lines.append("Pairing remains case-level. Matrix headers contain no aliquot barcode.")
    lines.append(f"Per-gene complete-observation rows written: {per_gene_n}. Rows are ordered by cancer and gene identifier. They are not ordered by an outcome.")
    lines.append("")
    lines.append("The detectable partial R² column uses the one-coefficient large-sample anchor λ = 7.848880 at the median primary-eligible N. It is a planning number, not an observed effect. OV is not given a lower eligibility gate.")
    (REP / "cohort_integrity_report.md").write_text("\n".join(lines) + "\n")
    print(f"cancers_with_triple_files {len(inter_rows)} pairing_rows {len(pair_rows)} per_gene_rows {per_gene_n}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--annotation-only":
        _, _, _, counts = load_annotation()
        print(counts)
    else:
        main()
