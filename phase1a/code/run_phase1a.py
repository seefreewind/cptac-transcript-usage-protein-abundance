#!/usr/bin/env python3
"""Phase 1A frozen-model execution. Protein regressions start only after the gates pass."""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import csv
import gzip
import hashlib
import platform
import sys
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    CANCERS,
    EXPECTED,
    GTF,
    GTF_SHA,
    LOG2E,
    PHASE1A,
    REPAIRED_SHA,
    ROOT,
    acquisition_path,
    header_names,
    historical_path,
    load_annotation,
    load_gene_matrix,
    load_ledger,
    NA_TOKENS,
    load_phenotype,
    parse_numeric_row,
    phase0_triple_ids,
    sha256_and_nul,
    triple_from_headers,
    write_csv,
)
from models import (  # noqa: E402
    NumericalFailure,
    benjamini_hochberg,
    fit_nested,
    jackknife_mean,
    partial_r2_fast,
    reml_meta,
    run_synthetic_tests,
)

PRIMARY_FIELDS = [
    "cancer", "gene_id", "principal_transcript_id", "principal_rule", "n_isoform_rows",
    "transformation", "n_candidate_cases", "n_complete", "n_stored_all_zero",
    "fraction_stored_all_zero", "rank_m0", "rank_m1", "df_residual_m1", "condition_number",
    "cnv_vif", "r2_m0", "r2_m1", "adjusted_r2_m0", "adjusted_r2_m1", "partial_r2",
    "expected_null_partial_r2", "null_excess_partial_r2", "signed_partial_r",
    "beta_p_standardized", "se_beta_p", "f_statistic", "p_value", "bh_q_value", "fit_status",
]
HIST_GE71 = {
    "BRCA": 5933, "CCRCC": 5385, "COAD": 3971, "GBM": 6189, "HNSCC": 5678,
    "LSCC": 6548, "LUAD": 6239, "OV": 4250, "PDAC": 5260, "UCEC": 5512,
}
HIST_MEAN = {
    "BRCA": 90.3395, "CCRCC": 83.2773, "COAD": 72.8404, "GBM": 84.1828, "HNSCC": 88.33,
    "LSCC": 91.6051, "LUAD": 90.2311, "OV": 61.4289, "PDAC": 82.6167, "UCEC": 71.511,
}
HIST_MEDIAN = {
    "BRCA": 109, "CCRCC": 103, "COAD": 95, "GBM": 98, "HNSCC": 107,
    "LSCC": 107, "LUAD": 108, "OV": 78, "PDAC": 105, "UCEC": 86,
}
STARTED = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (float, np.floating)):
        if not np.isfinite(value):
            return ""
        return f"{float(value):.10g}"
    return str(value)


def transcript_ids(path: Path) -> list[str]:
    ids = []
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip():
                continue
            if b"\x00" in line:
                raise SystemExit(f"NUL byte in transcript table {path}")
            ids.append(line.split(b"\t", 1)[0].decode())
    return ids


def accumulate(path: Path, case_ids: list[str], expected_ids: list[str], gene_index: dict[str, int],
               tx_to_gene: dict[str, str], principal: dict[str, str]) -> dict:
    n_cases = len(case_ids)
    n_genes = len(gene_index)
    sum_a = np.zeros((n_genes, n_cases))
    sum_b = np.zeros((n_genes, n_cases))
    sum_c = np.zeros((n_genes, n_cases))
    part_a = np.zeros((n_genes, n_cases))
    part_b = np.zeros((n_genes, n_cases))
    part_c = np.zeros((n_genes, n_cases))
    n_pos = np.zeros((n_genes, n_cases), dtype=np.int16)
    names = header_names(path)
    locate = {name: i for i, name in enumerate(names)}
    columns = [locate[case] for case in case_ids]
    expected = set(expected_ids)
    seen: set[str] = set()
    rows = 0
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip():
                continue
            if b"\x00" in line:
                raise SystemExit(f"NUL byte in {path}")
            parts = line.rstrip(b"\r\n").split(b"\t")
            tx = parts[0].decode()
            if tx in seen or tx not in expected:
                raise SystemExit(f"unexpected or duplicate transcript {tx} in {path}")
            seen.add(tx)
            rows += 1
            gene = tx_to_gene[tx]
            gi = gene_index.get(gene)
            if gi is None:
                if rows % 40000 == 0:
                    print(f"    isoform rows {rows}", flush=True)
                continue
            y = parse_numeric_row(parts, columns)
            if np.any(~np.isfinite(y)):
                raise SystemExit(f"non-finite isoform value in {path} at {tx}")
            power = np.exp2(y)
            sum_a[gi] += power
            sum_b[gi] += np.expm1(y * LOG2E)
            positive = y > 0
            sum_c[gi] += np.where(positive, power, 0.0)
            n_pos[gi] += positive.astype(np.int16)
            if principal.get(gene) == tx:
                part_a[gi] += power
                part_b[gi] += np.expm1(y * LOG2E)
                part_c[gi] += np.where(positive, power, 0.0)
            if rows % 40000 == 0:
                print(f"    isoform rows {rows}", flush=True)
    if seen != expected:
        raise SystemExit(f"{path} transcript set differs from the BRCA isoform matrix by {len(seen ^ expected)} ids")
    if np.any(part_a <= 0):
        raise SystemExit(f"principal row missing from an accumulated gene in {path}")
    log_b = np.full(sum_b.shape, np.nan)
    log_c = np.full(sum_c.shape, np.nan)
    p_b = np.full(sum_b.shape, np.nan)
    p_c = np.full(sum_c.shape, np.nan)
    np.log2(sum_b, out=log_b, where=sum_b > 0)
    np.log2(sum_c, out=log_c, where=sum_c > 0)
    np.divide(part_b, sum_b, out=p_b, where=sum_b > 0)
    np.divide(part_c, sum_c, out=p_c, where=sum_c > 0)
    return {
        "logG_A": np.log2(sum_a),
        "P_A": part_a / sum_a,
        "all_zero": n_pos == 0,
        "logG_B": log_b,
        "P_B": p_b,
        "logG_C": log_c,
        "P_C": p_c,
    }


def validate_inputs() -> list[dict]:
    ledger = load_ledger(ROOT / "phase0_5" / "SHA256SUMS")
    rows = []
    hold = []
    targets = [("ALL", "gencode_v34", GTF, GTF_SHA)]
    for cancer in CANCERS:
        for layer in ("rna_gene", "rna_isoform", "protein", "cnv", "phenotype"):
            path = acquisition_path(cancer, layer)
            expected = REPAIRED_SHA.get((cancer, layer))
            targets.append((cancer, layer, path, expected))
    for cancer, layer, path, expected in targets:
        if path.name.startswith("._") or path.suffix == ".partial" or ".partial" in path.name:
            hold.append(f"{path} is an AppleDouble or partial file")
            continue
        if not path.exists():
            hold.append(f"missing {path}")
            continue
        digest, size, nul, newlines = sha256_and_nul(path)
        rel = str(path.relative_to(ROOT))
        if expected is not None and digest != expected:
            hold.append(f"unexpected hash for {rel}")
        if (cancer, layer) not in REPAIRED_SHA and layer != "gencode_v34":
            ledger_key = f"phase0_5/acquisition/{cancer}/{path.name}"
            if ledger.get(ledger_key) != digest:
                hold.append(f"hash differs from the Phase 0.5 ledger for {ledger_key}")
        gzip_container = path.suffix == ".gz"
        if nul and not gzip_container:
            hold.append(f"NUL bytes in used file {rel}: {nul}")
        names = [] if layer == "gencode_v34" else header_names(path)
        rows.append({
            "cancer": cancer,
            "layer": layer,
            "path": rel,
            "sha256": digest,
            "bytes": size,
            "nul_bytes": nul,
            "newlines": newlines,
            "n_header_columns": len(names),
            "header_status": "gzip_container" if gzip_container else "idx",
            "expected_hash_checked": "yes" if expected else "phase0_5_ledger",
        })
        print(f"validated {rel} nul={nul}", flush=True)
    write_csv(PHASE1A / "manifests" / "input_manifest.csv", rows,
              ["cancer", "layer", "path", "sha256", "bytes", "nul_bytes", "newlines",
               "n_header_columns", "header_status", "expected_hash_checked"])
    (PHASE1A / "manifests" / "input_SHA256SUMS").write_text(
        "".join(f"{row['sha256']}  {row['path']}\n" for row in rows)
    )
    lines = ["# Input validation", "", f"Validated at {now()}.", ""]
    if hold:
        lines.append("Status: TECHNICAL_HOLD")
        lines.extend(f"- {item}" for item in hold)
    else:
        lines.append("Status: PASS. Every used file matches its required hash and contains no NUL bytes.")
    lines.append("")
    lines.append("Repaired copies are used for BRCA isoform RNA, UCEC gene RNA, and PDAC CNV. The Phase 0.5 copies of those three files were not read.")
    lines.append("GENCODE v34 is a gzip container. NUL bytes inside that container are compressed-stream bytes. The file is accepted when its SHA256 equals the frozen GENCODE v34 hash. Text matrices are still rejected if they contain NUL bytes.")
    (PHASE1A / "reports" / "input_validation_report.md").write_text("\n".join(lines) + "\n")
    if hold:
        raise SystemExit("TECHNICAL_HOLD\n" + "\n".join(hold))
    return rows


def build_principal(tx_to_gene, principal, rule, transcript_list):
    counts: dict[str, int] = {}
    id_set = set(transcript_list)
    for tx in transcript_list:
        gene = tx_to_gene.get(tx)
        if gene is None:
            raise SystemExit(f"unmapped transcript {tx}")
        counts[gene] = counts.get(gene, 0) + 1
    candidates = []
    missing = []
    multi_missing = []
    rows = []
    for gene in sorted(principal):
        tx = principal[gene]
        present = tx in id_set
        n_tx = counts.get(gene, 0)
        in_universe = present and n_tx >= 2
        if in_universe:
            candidates.append(gene)
        if not present:
            missing.append(gene)
            if n_tx >= 2:
                multi_missing.append(gene)
        rows.append({
            "gene_id": gene,
            "principal_transcript_id": tx,
            "principal_rule": rule[gene],
            "in_isoform_matrix": int(present),
            "n_isoform_rows": n_tx,
            "in_candidate_universe": int(in_universe),
        })
    write_csv(
        PHASE1A / "manifests" / "frozen_principal_transcripts.csv",
        rows,
        ["gene_id", "principal_transcript_id", "principal_rule", "in_isoform_matrix", "n_isoform_rows", "in_candidate_universe"],
    )
    reproduction = [
        {"metric": "genes_with_principal", "value": len(principal), "expected": 17439, "match": int(len(principal) == 17439)},
        {"metric": "principal_transcripts_absent", "value": len(missing), "expected": 37, "match": int(len(missing) == 37)},
        {"metric": "candidate_genes", "value": len(candidates), "expected": 10795, "match": int(len(candidates) == 10795)},
        {"metric": "isoform_rows", "value": len(transcript_list), "expected": 108777, "match": int(len(transcript_list) == 108777)},
        {"metric": "genes_with_mane", "value": sum(item == "MANE_Select" for item in rule.values()), "expected": 14620, "match": int(sum(item == "MANE_Select" for item in rule.values()) == 14620)},
    ]
    write_csv(PHASE1A / "qc" / "principal_transcript_reproduction.csv", reproduction, ["metric", "value", "expected", "match"])
    if any(not row["match"] for row in reproduction):
        raise SystemExit("principal-transcript reproduction does not match the frozen counts")
    return candidates, {gene: counts[gene] for gene in candidates}, multi_missing


def reconstruct_cohort(phase0_ids) -> dict[str, list[str]]:
    rows = []
    triples = {}
    problems = []
    for cancer in CANCERS:
        triple = triple_from_headers(cancer)
        triples[cancer] = triple
        got = set(triple)
        old = phase0_ids[cancer]
        if len(triple) != EXPECTED[cancer] or got != old:
            problems.append(f"{cancer}: header triple {len(triple)} expected {EXPECTED[cancer]}; symmetric difference {sorted(got ^ old)[:8]}")
        cnv_names = set(header_names(acquisition_path(cancer, "cnv")))
        purity, estimate = load_phenotype(acquisition_path(cancer, "phenotype"), triple)
        for i, sample in enumerate(triple):
            rows.append({
                "cancer": cancer,
                "sample_id": sample,
                "in_triple": 1,
                "in_phase0_record": int(sample in old),
                "in_cnv_header": int(sample in cnv_names),
                "wes_purity_numeric": int(np.isfinite(purity[i])),
                "estimate_numeric": int(np.isfinite(estimate[i])),
                "same_aliquot": "no_aliquot_barcode_in_matrix_header",
            })
    write_csv(
        PHASE1A / "qc" / "cohort_reconstruction.csv",
        rows,
        ["cancer", "sample_id", "in_triple", "in_phase0_record", "in_cnv_header", "wes_purity_numeric", "estimate_numeric", "same_aliquot"],
    )
    text = ["# Cohort reconstruction", "", f"Recomputed from complete headers at {now()}.", ""]
    if problems:
        text.append("Status: MISMATCH")
        text.extend(f"- {item}" for item in problems)
    else:
        text.append("Status: PASS. Each cancer reproduces the Phase 0 case set. The total is 1,027 case-level matches. COAD contributes its 96 protein-overlapping cases. The match is not an aliquot match.")
    text.append("")
    for cancer in CANCERS:
        text.append(f"- {cancer}: {len(triples[cancer])}")
    (PHASE1A / "reports" / "cohort_reconstruction_report.md").write_text("\n".join(text) + "\n")
    if problems:
        raise SystemExit("cohort reconstruction mismatch")
    return triples


def prepare_cancer(cancer, case_ids, expected_ids, gene_ids, gene_index, n_tx, tx_to_gene, principal, rule):
    print(f"{cancer} loading covariates", flush=True)
    protein = load_gene_matrix(acquisition_path(cancer, "protein"), case_ids)
    cnv = load_gene_matrix(acquisition_path(cancer, "cnv"), case_ids, allow_missing=True)
    rsem = load_gene_matrix(acquisition_path(cancer, "rna_gene"), case_ids)
    purity, estimate = load_phenotype(acquisition_path(cancer, "phenotype"), case_ids)
    print(f"{cancer} accumulating isoforms", flush=True)
    iso = accumulate(acquisition_path(cancer, "rna_isoform"), case_ids, expected_ids, gene_index, tx_to_gene, principal)
    n_cases = len(case_ids)
    n_genes = len(gene_ids)
    protein_m = np.full((n_genes, n_cases), np.nan)
    cnv_m = np.full((n_genes, n_cases), np.nan)
    rsem_m = np.full((n_genes, n_cases), np.nan)
    in_protein = np.zeros(n_genes, dtype=bool)
    in_cnv = np.zeros(n_genes, dtype=bool)
    in_rsem = np.zeros(n_genes, dtype=bool)
    for i, gene in enumerate(gene_ids):
        if gene in protein:
            protein_m[i] = protein[gene]
            in_protein[i] = True
        if gene in cnv:
            cnv_m[i] = cnv[gene]
            in_cnv[i] = True
        if gene in rsem:
            rsem_m[i] = rsem[gene]
            in_rsem[i] = True
    pack = {
        "cancer": cancer,
        "case_ids": case_ids,
        "gene_ids": gene_ids,
        "n_tx": np.array([n_tx[gene] for gene in gene_ids], dtype=np.int32),
        "principal_tx": [principal[gene] for gene in gene_ids],
        "principal_rule": [rule[gene] for gene in gene_ids],
        "protein": protein_m,
        "cnv": cnv_m,
        "rsem": rsem_m,
        "purity": purity,
        "estimate": estimate,
        "in_protein": in_protein,
        "in_cnv": in_cnv,
        "in_rsem": in_rsem,
        "protein_ids": set(protein),
        "rsem_ids": set(rsem),
        **iso,
    }
    return pack


def _first_numeric_masks(path: Path, case_ids: list[str]) -> dict[str, np.ndarray]:
    """Numeric-observation mask from a historical matrix. Lines containing NUL are skipped."""
    names = header_names(path)
    locate = {name: i for i, name in enumerate(names)}
    columns = [locate.get(case) for case in case_ids]
    masks: dict[str, np.ndarray] = {}
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip() or b"\x00" in line:
                continue
            parts = line.rstrip(b"\r\n").split(b"\t")
            gene = parts[0].decode()
            if gene in masks:
                continue
            mask = np.zeros(len(case_ids), dtype=bool)
            width = len(parts)
            for i, col in enumerate(columns):
                if col is None or col + 1 >= width:
                    continue
                token = parts[col + 1].strip()
                if token in NA_TOKENS:
                    continue
                try:
                    value = float(token)
                except ValueError:
                    continue
                if np.isfinite(value):
                    mask[i] = True
            masks[gene] = mask
    return masks


def _historical_usable(path: Path, case_ids: list[str], tx_to_gene: dict[str, str], principal: dict[str, str]):
    """Phase 0.5 Literal B usable mask. NUL lines are not used as observations."""
    names = header_names(path)
    locate = {name: i for i, name in enumerate(names)}
    columns = [locate[case] for case in case_ids]
    n_cases = len(case_ids)
    n_tx: dict[str, int] = {}
    n_pos: dict[str, np.ndarray] = {}
    principal_seen: dict[str, np.ndarray] = {}
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip() or b"\x00" in line:
                continue
            parts = line.rstrip(b"\r\n").split(b"\t")
            gene = tx_to_gene.get(parts[0].decode())
            if gene is None:
                continue
            n_tx[gene] = n_tx.get(gene, 0) + 1
            if gene not in n_pos:
                n_pos[gene] = np.zeros(n_cases, dtype=np.int16)
                principal_seen[gene] = np.zeros(n_cases, dtype=bool)
            y = parse_numeric_row(parts, columns)
            finite = np.isfinite(y)
            n_pos[gene] += ((y > 0) & finite).astype(np.int16)
            if principal.get(gene) == parts[0].decode():
                principal_seen[gene] |= finite
    usable = {gene: (n_pos[gene] > 0) & principal_seen[gene] for gene in n_pos}
    return n_tx, usable


def recount_historical(cancer: str, case_ids: list[str], tx_to_gene, principal) -> tuple[int, float, float]:
    """Read-only recount on the Phase 0.5 files. Damaged lines are skipped and are not model inputs."""
    print(f"{cancer} historical recount", flush=True)
    rna = _first_numeric_masks(historical_path(cancer, "rna_gene"), case_ids)
    protein = _first_numeric_masks(historical_path(cancer, "protein"), case_ids)
    cnv = _first_numeric_masks(historical_path(cancer, "cnv"), case_ids)
    purity, _estimate = load_phenotype(historical_path(cancer, "phenotype"), case_ids)
    purity_ok = np.isfinite(purity)
    n_tx, usable = _historical_usable(historical_path(cancer, "rna_isoform"), case_ids, tx_to_gene, principal)
    values = []
    ge71 = 0
    for gene in rna.keys() & protein.keys():
        rows = n_tx.get(gene, 0)
        use = usable.get(gene)
        if use is None or gene not in cnv:
            n_complete = 0
        else:
            n_complete = int(np.sum(rna[gene] & protein[gene] & cnv[gene] & purity_ok & use))
        if rows >= 2 and n_complete >= 71:
            ge71 += 1
        if rows >= 2 and gene in principal:
            values.append(n_complete)
    ordered = sorted(values)
    mid = len(ordered) // 2
    median = ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2
    mean = sum(values) / len(values)
    return ge71, mean, float(median)


def historical_literal_b(pack, multi_missing) -> tuple[int, float, float]:
    values = []
    ge71 = 0
    finite_rsem = np.isfinite(pack["rsem"])
    finite_protein = np.isfinite(pack["protein"])
    finite_cnv = np.isfinite(pack["cnv"])
    purity = np.isfinite(pack["purity"])
    usable = np.isfinite(pack["logG_B"]) & np.isfinite(pack["P_B"])
    for i, gene in enumerate(pack["gene_ids"]):
        if not pack["in_protein"][i] or not pack["in_rsem"][i]:
            continue
        n_complete = int(np.sum(finite_protein[i] & finite_rsem[i] & finite_cnv[i] & purity & usable[i]))
        values.append(n_complete)
        if n_complete >= 71:
            ge71 += 1
    for gene in multi_missing:
        if gene in pack["protein_ids"] and gene in pack["rsem_ids"]:
            values.append(0)
    ordered = sorted(values)
    mid = len(ordered) // 2
    median = ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2
    mean = sum(values) / len(values)
    return ge71, mean, float(median)


def fit_rows(pack, transform: str, clips: list) -> tuple[list[dict], np.ndarray, list[np.ndarray]]:
    if transform == "A":
        log_g_all, prop_all, all_zero = pack["logG_A"], pack["P_A"], pack["all_zero"]
    elif transform == "B":
        log_g_all, prop_all, all_zero = pack["logG_B"], pack["P_B"], np.zeros_like(pack["all_zero"])
    else:
        log_g_all, prop_all, all_zero = pack["logG_C"], pack["P_C"], np.zeros_like(pack["all_zero"])
    purity_ok = np.isfinite(pack["purity"])
    n_genes = len(pack["gene_ids"])
    partial = np.full(n_genes, np.nan)
    membership = [np.empty(0, dtype=np.int32) for _ in range(n_genes)]
    stored = [None] * n_genes
    p_values = np.full(n_genes, np.nan)
    for i in range(n_genes):
        mask = (
            np.isfinite(pack["protein"][i])
            & np.isfinite(pack["cnv"][i])
            & purity_ok
            & np.isfinite(log_g_all[i])
            & np.isfinite(prop_all[i])
        )
        idx = np.flatnonzero(mask).astype(np.int32)
        n_complete = int(idx.size)
        n_zero = int(all_zero[i, idx].sum()) if n_complete else 0
        record = {
            "fit_status": "BELOW_N",
            "n_complete": n_complete,
            "n_stored_all_zero": n_zero,
        }
        if not pack["in_protein"][i]:
            record["fit_status"] = "NOT_IN_PROTEIN_MATRIX"
        elif not pack["in_cnv"][i]:
            record["fit_status"] = "NOT_IN_CNV_MATRIX"
        elif n_complete >= 71:
            result = fit_nested(
                log_g_all[i, idx], prop_all[i, idx], pack["purity"][idx], pack["cnv"][i, idx], pack["protein"][i, idx], clips
            )
            record.update(result)
            record["n_complete"] = n_complete
            record["n_stored_all_zero"] = n_zero
            if result["fit_status"] == "OK":
                partial[i] = result["partial_r2"]
                membership[i] = idx
                p_values[i] = result["p_value"]
        stored[i] = record
    q_values = benjamini_hochberg(p_values)
    rows = []
    for i, gene in enumerate(pack["gene_ids"]):
        result = stored[i]
        n_complete = result.get("n_complete", 0)
        n_zero = result.get("n_stored_all_zero", 0)
        row = {field: "" for field in PRIMARY_FIELDS}
        row.update({
            "cancer": pack["cancer"],
            "gene_id": gene,
            "principal_transcript_id": pack["principal_tx"][i],
            "principal_rule": pack["principal_rule"][i],
            "n_isoform_rows": int(pack["n_tx"][i]),
            "transformation": transform,
            "n_candidate_cases": len(pack["case_ids"]),
            "n_complete": n_complete,
            "n_stored_all_zero": n_zero,
            "fraction_stored_all_zero": fmt(n_zero / n_complete) if n_complete else "",
            "fit_status": result["fit_status"],
        })
        for key in (
            "rank_m0", "rank_m1", "df_residual_m1", "condition_number", "cnv_vif", "r2_m0", "r2_m1",
            "adjusted_r2_m0", "adjusted_r2_m1", "partial_r2", "expected_null_partial_r2",
            "null_excess_partial_r2", "signed_partial_r", "beta_p_standardized", "se_beta_p",
            "f_statistic", "p_value",
        ):
            if key in result:
                row[key] = fmt(result[key])
        if result["fit_status"] == "OK":
            row["bh_q_value"] = fmt(q_values[i])
        rows.append(row)
    return rows, partial, membership


def cancer_summary(rows, jack=None) -> dict:
    ok = [row for row in rows if row["fit_status"] == "OK"]
    partial = np.array([float(row["partial_r2"]) for row in ok], dtype=float)
    expected = np.array([float(row["expected_null_partial_r2"]) for row in ok], dtype=float)
    excess = np.array([float(row["null_excess_partial_r2"]) for row in ok], dtype=float)
    qvals = np.array([float(row["bh_q_value"]) for row in ok], dtype=float)
    counts = Counter(row["fit_status"] for row in rows)
    prefit_status = {"OK", "RANK_DEFICIENT", "NO_RESIDUAL_P_VARIATION", "ZERO_PROTEIN_VARIANCE", "M0_RSS_ZERO"}
    out = {
        "cancer": rows[0]["cancer"],
        "transformation": rows[0]["transformation"],
        "n_candidate": len(rows),
        "n_prefit": sum(row["fit_status"] in prefit_status for row in rows),
        "n_success": len(ok),
        "n_rank_deficient": counts["RANK_DEFICIENT"],
        "n_no_residual_p": counts["NO_RESIDUAL_P_VARIATION"],
        "mean_partial_r2": float(partial.mean()) if partial.size else np.nan,
        "median_partial_r2": float(np.median(partial)) if partial.size else np.nan,
        "sd_partial_r2": float(partial.std(ddof=1)) if partial.size > 1 else np.nan,
        "q25": float(np.quantile(partial, 0.25)) if partial.size else np.nan,
        "q75": float(np.quantile(partial, 0.75)) if partial.size else np.nan,
        "p90": float(np.quantile(partial, 0.90)) if partial.size else np.nan,
        "p95": float(np.quantile(partial, 0.95)) if partial.size else np.nan,
        "mean_expected_null": float(expected.mean()) if expected.size else np.nan,
        "mean_null_excess": float(excess.mean()) if excess.size else np.nan,
        "prop_ge_0.01": float(np.mean(partial >= 0.01)) if partial.size else np.nan,
        "prop_ge_0.02": float(np.mean(partial >= 0.02)) if partial.size else np.nan,
        "prop_ge_0.05": float(np.mean(partial >= 0.05)) if partial.size else np.nan,
        "prop_ge_0.10": float(np.mean(partial >= 0.10)) if partial.size else np.nan,
        "prop_q_lt_0.05": float(np.mean(qvals < 0.05)) if qvals.size else np.nan,
        "prop_q_and_r2": float(np.mean((qvals < 0.05) & (partial >= 0.02))) if partial.size else np.nan,
    }
    if jack:
        out.update({
            "jackknife_n": jack["n_jackknife"],
            "jackknife_se": jack["se"],
            "ci_low": jack["ci_low"],
            "ci_high": jack["ci_high"],
            "ci_low_clipped": np.clip(jack["ci_low"], 0, 1),
            "ci_high_clipped": np.clip(jack["ci_high"], 0, 1),
            "jackknife_min_success": jack["min_success"],
            "jackknife_max_success": jack["max_success"],
            "jackknife_genes_always_ok": jack["genes_always_ok"],
            "jackknife_unstable": int(jack["unstable"]),
        })
    return out


def write_primary(path: Path, rows: list[dict]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PRIMARY_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sensitivity_fit(pack, idx_list, transform, mode: str, clips: list) -> list[float]:
    values = []
    if transform == "A":
        log_g, prop, all_zero = pack["logG_A"], pack["P_A"], pack["all_zero"]
    elif transform == "B":
        log_g, prop, all_zero = pack["logG_B"], pack["P_B"], None
    else:
        log_g, prop, all_zero = pack["logG_C"], pack["P_C"], None
    for i, idx in enumerate(idx_list):
        if idx.size == 0:
            continue
        use = idx
        this_log = log_g[i]
        this_purity = pack["purity"]
        this_cnv = pack["cnv"][i]
        if mode == "drop_all_zero":
            use = idx[~all_zero[i, idx]]
        elif mode == "rsem":
            keep = np.isfinite(pack["rsem"][i, idx])
            use = idx[keep]
            this_log = pack["rsem"][i]
        elif mode == "estimate":
            keep = np.isfinite(pack["estimate"][idx])
            use = idx[keep]
            this_purity = pack["estimate"]
        elif mode == "no_cnv":
            y = pack["protein"][i, use]
            result = fit_nested(this_log[use], prop[i, use], this_purity[use], np.zeros(use.size), y, clips)
            # A constant CNV makes the design rank-deficient. Refit by passing a
            # non-constant dummy only when the dedicated no-CNV design is used below.
            result = None
        if mode == "no_cnv":
            if use.size < 71:
                continue
            y = pack["protein"][i, use]
            x_log = this_log[use]
            x_p = prop[i, use]
            x_pur = this_purity[use]
            # Three-predictor sensitivity: drop CNV by giving fit_nested a column that
            # is a copy of purity, then accept only if we instead call a reduced design.
            reduced = _fit_without_cnv(x_log, x_p, x_pur, y, clips)
            if reduced is not None:
                values.append(reduced)
            continue
        if use.size < 71:
            continue
        result = fit_nested(this_log[use], prop[i, use], this_purity[use], this_cnv[use], pack["protein"][i, use], clips)
        if result["fit_status"] == "OK":
            values.append(result["partial_r2"])
    return values


def _fit_without_cnv(log_g, proportion, purity, protein, clips) -> float | None:
    from models import _design, _svd_fit
    x0 = _design([log_g, purity])
    x1 = _design([log_g, proportion, purity])
    if x0 is None or x1 is None:
        return None
    fit0 = _svd_fit(x0, protein)
    fit1 = _svd_fit(x1, protein)
    if not fit0["ok"] or not fit1["ok"] or fit0["rss"] <= 0:
        return None
    gap = fit0["rss"] - fit1["rss"]
    if gap < -CLIP_TOL:
        raise NumericalFailure("sensitivity RSS inversion")
    return max(gap, 0.0) / fit0["rss"]


CLIP_TOL = 1e-8


def paired_stats(partial_a, partial_b) -> dict:
    both = np.isfinite(partial_a) & np.isfinite(partial_b)
    a = partial_a[both]
    b = partial_b[both]
    if a.size < 3:
        return {"n_paired": int(a.size)}
    diff = a - b
    return {
        "n_paired": int(a.size),
        "spearman": float(stats.spearmanr(a, b).statistic),
        "pearson": float(stats.pearsonr(a, b).statistic),
        "mean_diff_a_minus_b": float(diff.mean()),
        "median_diff_a_minus_b": float(np.median(diff)),
        "mean_abs_diff": float(np.mean(np.abs(diff))),
        "prop_abs_diff_gt_0.01": float(np.mean(np.abs(diff) > 0.01)),
        "prop_abs_diff_gt_0.02": float(np.mean(np.abs(diff) > 0.02)),
        "both_ge_0.02": int(np.sum((a >= 0.02) & (b >= 0.02))),
        "only_a_ge_0.02": int(np.sum((a >= 0.02) & (b < 0.02))),
        "only_b_ge_0.02": int(np.sum((a < 0.02) & (b >= 0.02))),
        "neither_ge_0.02": int(np.sum((a < 0.02) & (b < 0.02))),
    }


def classify_science(meta_a, meta_b, paired_a, paired_b, loo_a, loo_b) -> str:
    def above(meta):
        return bool(np.isfinite(meta["ci_low"]) and meta["ci_low"] > 0.02)

    if above(meta_a) != above(meta_b):
        return "TRANSFORMATION_DEPENDENT"
    if above(meta_a) and above(meta_b):
        loo_ok = all(row["ci_low"] > 0.02 for row in loo_a + loo_b)
        consistent = above(paired_a) and above(paired_b) and meta_a["pi_low"] > 0.02 and meta_b["pi_low"] > 0.02 and loo_ok
        return "ROBUST_PRACTICAL_INCREMENT" if consistent else "CANCER_SPECIFIC_OR_HETEROGENEOUS"
    return "SMALL_OR_UNCERTAIN_INCREMENT"


def write_reports(science, execution, notes: list[str]) -> None:
    final = f"""# Phase 1A final decision

## A. Execution status

{execution}

## B. Scientific result

{science}

## C. Evidence boundary

The 1,027 tumors are matched at case level. The protein outcome is the stored gene-level abundance. Literal A and Literal B are co-primary transformations, and neither was selected from the protein fits. This analysis does not identify an isoform-specific protein abundance and does not demonstrate protein functional loss. Partial R² is incremental statistical explanation. M0 and M1 use RNA-derived variables from the same isoform matrix. Under Literal A, a gene–case whose stored isoform values are all zero remains in the model with equal transcript weights; that behavior is tested by sensitivity S1.

## D. Next-phase authorization

PHASE 1B AUTHORIZATION = FALSE
"""
    (PHASE1A / "reports" / "phase1a_final_decision.md").write_text(final)
    limit = ["# Phase 1A limitations", ""]
    limit.extend(f"- {note}" for note in notes)
    limit.append("- CLR sensitivity S4 was not run because the frozen documents define the pseudocount and do not define the one-degree-of-freedom contrast.")
    limit.append("- Phase 1B has not been authorized.")
    (PHASE1A / "reports" / "phase1a_limitations.md").write_text("\n".join(limit) + "\n")


def main() -> None:
    for folder in ("protocol", "code", "logs", "qc", "results", "reports", "manifests"):
        (PHASE1A / folder).mkdir(parents=True, exist_ok=True)
    env = [
        f"start_utc {STARTED}",
        f"system {platform.platform()}",
        f"python {sys.version}",
        f"numpy {np.__version__}",
        f"scipy {stats.__file__}",
        "threads OMP=1 OPENBLAS=1 MKL=1 VECLIB=1",
        "blas accelerate",
        "standardize ddof=0 predictors only",
        "seed permutation_audit_only 20261114",
    ]
    import scipy
    env[4] = f"scipy {scipy.__version__}"
    (PHASE1A / "manifests" / "software_environment.txt").write_text("\n".join(env) + "\n")
    prefit = []
    for path in sorted((PHASE1A / "protocol").glob("*")) + sorted((PHASE1A / "code").glob("*.py")):
        if path.name.startswith("._"):
            continue
        digest, _, _, _ = sha256_and_nul(path)
        prefit.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
    (PHASE1A / "manifests" / "protocol_prefit_SHA256SUMS").write_text("\n".join(prefit) + "\n")
    print("protocol hashed", now(), flush=True)
    validate_inputs()
    print("inputs validated", now(), flush=True)
    tx_to_gene, principal, rule = load_annotation()
    expected_ids = transcript_ids(acquisition_path("BRCA", "rna_isoform"))
    gene_ids, n_tx, multi_missing = build_principal(tx_to_gene, principal, rule, expected_ids)
    gene_index = {gene: i for i, gene in enumerate(gene_ids)}
    print(f"candidate genes {len(gene_ids)} multi_missing_principal {len(multi_missing)}", flush=True)
    triples = reconstruct_cohort(phase0_triple_ids())
    packs = {}
    reproduction_rows = []
    unresolved = []
    for cancer in CANCERS:
        pack = prepare_cancer(cancer, triples[cancer], expected_ids, gene_ids, gene_index, n_tx, tx_to_gene, principal, rule)
        packs[cancer] = pack
        ge71, mean, median = historical_literal_b(pack, multi_missing)
        ge_ok = ge71 == HIST_GE71[cancer]
        mean_ok = abs(round(mean, 4) - HIST_MEAN[cancer]) <= 0.00015
        median_ok = abs(median - HIST_MEDIAN[cancer]) <= 1e-6
        note = "matched_repaired_inputs"
        hist_ge = hist_mean = hist_median = ""
        if not (ge_ok and mean_ok and median_ok):
            if cancer in {"BRCA", "UCEC", "PDAC"}:
                h_ge, h_mean, h_med = recount_historical(cancer, pack["case_ids"], tx_to_gene, principal)
                hist_ge, hist_mean, hist_median = h_ge, round(h_mean, 4), h_med
                hist_ok = (
                    h_ge == HIST_GE71[cancer]
                    and abs(round(h_mean, 4) - HIST_MEAN[cancer]) <= 0.00015
                    and abs(h_med - HIST_MEDIAN[cancer]) <= 1e-6
                )
                if hist_ok:
                    note = "repaired_input_differs_historical_recount_matches"
                else:
                    note = "MISMATCH_after_historical_recount"
                    unresolved.append(cancer)
                print(f"{cancer} historical recount >=71 {h_ge} mean {h_mean:.4f} median {h_med}", flush=True)
            else:
                note = "MISMATCH"
                unresolved.append(cancer)
        reproduction_rows.append({
            "cancer": cancer,
            "literal_b_ge71": ge71,
            "expected_ge71": HIST_GE71[cancer],
            "ge71_match": int(ge_ok),
            "mean_n": round(mean, 4),
            "expected_mean_n": HIST_MEAN[cancer],
            "mean_match": int(mean_ok),
            "median_n": median,
            "expected_median_n": HIST_MEDIAN[cancer],
            "median_match": int(median_ok),
            "historical_recount_ge71": hist_ge,
            "historical_recount_mean_n": hist_mean,
            "historical_recount_median_n": hist_median,
            "note": note,
        })
        print(f"{cancer} literal B >=71 {ge71} expected {HIST_GE71[cancer]} mean {mean:.4f} median {median}", flush=True)
    write_csv(
        PHASE1A / "qc" / "literal_b_prefit_reproduction.csv",
        reproduction_rows,
        ["cancer", "literal_b_ge71", "expected_ge71", "ge71_match", "mean_n", "expected_mean_n", "mean_match", "median_n", "expected_median_n", "median_match", "historical_recount_ge71", "historical_recount_mean_n", "historical_recount_median_n", "note"],
    )
    if unresolved:
        raise SystemExit("PREMODEL_REPRODUCTION_FAILURE " + ",".join(unresolved))
    print("synthetic tests", flush=True)
    synthetic = run_synthetic_tests()
    (PHASE1A / "qc" / "synthetic_test_results.txt").write_text("\n".join(synthetic) + "\n")
    print("protein models start", now(), flush=True)
    clips: list[float] = []
    primary_rows = {"A": [], "B": []}
    partials = {}
    memberships = {}
    summaries = []
    exclusions = []
    for cancer in CANCERS:
        partials[cancer] = {}
        memberships[cancer] = {}
        for transform in ("A", "B"):
            print(f"fit {cancer} {transform}", now(), flush=True)
            rows, partial, membership = fit_rows(packs[cancer], transform, clips)
            primary_rows[transform].extend(rows)
            partials[cancer][transform] = partial
            memberships[cancer][transform] = membership
            for row in rows:
                if row["fit_status"] != "OK":
                    exclusions.append({
                        "cancer": cancer,
                        "gene_id": row["gene_id"],
                        "transformation": transform,
                        "reason": row["fit_status"],
                        "n_complete": row["n_complete"],
                    })
    write_primary(PHASE1A / "results" / "primary_literal_A_gene_fits.tsv.gz", primary_rows["A"])
    write_primary(PHASE1A / "results" / "primary_literal_B_gene_fits.tsv.gz", primary_rows["B"])
    write_csv(PHASE1A / "results" / "exclusion_reasons.csv", exclusions, ["cancer", "gene_id", "transformation", "reason", "n_complete"])
    print("jackknife", now(), flush=True)
    jack_rows = []
    jack_by = {}
    for cancer in CANCERS:
        pack = packs[cancer]
        jack_by[cancer] = {}
        for transform in ("A", "B"):
            partial = partials[cancer][transform]
            membership = memberships[cancer][transform]
            log_g = pack["logG_A"] if transform == "A" else pack["logG_B"]
            prop = pack["P_A"] if transform == "A" else pack["P_B"]

            def refit(gene_index_i, dropped, log_g=log_g, prop=prop, membership=membership, pack=pack):
                idx = membership[gene_index_i]
                keep = idx[idx != dropped]
                if keep.size < 5:
                    return None
                return partial_r2_fast(log_g[gene_index_i, keep], prop[gene_index_i, keep], pack["purity"][keep], pack["cnv"][gene_index_i, keep], pack["protein"][gene_index_i, keep])

            print(f"  jackknife {cancer} {transform}", flush=True)
            jack = jackknife_mean(partial, membership, refit)
            jack_by[cancer][transform] = jack
            jack_rows.append({"cancer": cancer, "transformation": transform, **{key: fmt(value) if not isinstance(value, (int, np.integer, bool)) else value for key, value in jack.items()}})
    write_csv(PHASE1A / "qc" / "jackknife_diagnostics.csv", jack_rows, list(jack_rows[0].keys()) if jack_rows else ["cancer"])
    summary_rows = []
    by_cancer_rows = {"A": {}, "B": {}}
    for transform in ("A", "B"):
        rows_by_cancer = {}
        for row in primary_rows[transform]:
            rows_by_cancer.setdefault(row["cancer"], []).append(row)
        for cancer in CANCERS:
            summary = cancer_summary(rows_by_cancer[cancer], jack_by[cancer][transform])
            summary_rows.append(summary)
            by_cancer_rows[transform][cancer] = summary
    write_csv(PHASE1A / "results" / "cancer_level_summary.csv", summary_rows, list(summary_rows[0].keys()))
    write_csv(
        PHASE1A / "results" / "eligibility_by_cancer_transform.csv",
        [{key: row[key] for key in ("cancer", "transformation", "n_candidate", "n_prefit", "n_success", "n_rank_deficient", "n_no_residual_p")} for row in summary_rows],
        ["cancer", "transformation", "n_candidate", "n_prefit", "n_success", "n_rank_deficient", "n_no_residual_p"],
    )
    paired_rows = []
    common_partial = {"A": {}, "B": {}}
    for cancer in CANCERS:
        stats_row = paired_stats(partials[cancer]["A"], partials[cancer]["B"])
        stats_row["cancer"] = cancer
        paired_rows.append(stats_row)
        both = np.isfinite(partials[cancer]["A"]) & np.isfinite(partials[cancer]["B"])
        common_partial["A"][cancer] = np.where(both, partials[cancer]["A"], np.nan)
        common_partial["B"][cancer] = np.where(both, partials[cancer]["B"], np.nan)
    write_csv(PHASE1A / "results" / "transformation_paired_comparison.csv", paired_rows, list(paired_rows[0].keys()))
    print("common-support jackknife", now(), flush=True)
    common_jack = {"A": {}, "B": {}}
    for cancer in CANCERS:
        pack = packs[cancer]
        for transform in ("A", "B"):
            log_g = pack["logG_A"] if transform == "A" else pack["logG_B"]
            prop = pack["P_A"] if transform == "A" else pack["P_B"]
            membership = memberships[cancer][transform]

            def refit(gene_index_i, dropped, log_g=log_g, prop=prop, membership=membership, pack=pack):
                idx = membership[gene_index_i]
                keep = idx[idx != dropped]
                if keep.size < 5:
                    return None
                return partial_r2_fast(log_g[gene_index_i, keep], prop[gene_index_i, keep], pack["purity"][keep], pack["cnv"][gene_index_i, keep], pack["protein"][gene_index_i, keep])

            common_jack[transform][cancer] = jackknife_mean(common_partial[transform][cancer], membership, refit)
    common_rows = []
    for cancer in CANCERS:
        common_rows.append({
            "cancer": cancer,
            "n_paired": paired_rows[CANCERS.index(cancer)]["n_paired"],
            "mean_A": common_jack["A"][cancer]["full_mean"],
            "se_A": common_jack["A"][cancer]["se"],
            "ci_low_A": common_jack["A"][cancer]["ci_low"],
            "ci_high_A": common_jack["A"][cancer]["ci_high"],
            "mean_B": common_jack["B"][cancer]["full_mean"],
            "se_B": common_jack["B"][cancer]["se"],
            "ci_low_B": common_jack["B"][cancer]["ci_low"],
            "ci_high_B": common_jack["B"][cancer]["ci_high"],
        })
    write_csv(PHASE1A / "results" / "common_support_summary.csv", common_rows, list(common_rows[0].keys()))

    def meta_from(summary_map, transform, cancers, se_key="jackknife_se", mean_key="mean_partial_r2"):
        theta = np.array([summary_map[cancer][mean_key] for cancer in cancers], dtype=float)
        se = np.array([summary_map[cancer][se_key] for cancer in cancers], dtype=float)
        return reml_meta(theta, se)

    # Primary metas use cancer summaries. Common-support metas use common_jack.
    primary_meta = {}
    paired_meta = {}
    no_ov_meta = {}
    loo_rows = []
    meta_rows = []
    for transform in ("A", "B"):
        smap = {cancer: by_cancer_rows[transform][cancer] for cancer in CANCERS}
        primary_meta[transform] = meta_from(smap, transform, CANCERS)
        no_ov_meta[transform] = meta_from(smap, transform, [cancer for cancer in CANCERS if cancer != "OV"])
        theta = np.array([common_jack[transform][cancer]["full_mean"] for cancer in CANCERS])
        se = np.array([common_jack[transform][cancer]["se"] for cancer in CANCERS])
        paired_meta[transform] = reml_meta(theta, se)
        for label, meta in (("all_ten", primary_meta[transform]), ("exclude_OV", no_ov_meta[transform]), ("paired_support", paired_meta[transform])):
            meta_rows.append({"analysis": label, "transformation": transform, "method": "REML_HartungKnapp", **meta})
        for left_out in CANCERS:
            kept = [cancer for cancer in CANCERS if cancer != left_out]
            meta = meta_from(smap, transform, kept)
            loo_rows.append({"analysis": "all_ten", "transformation": transform, "left_out": left_out, "mu": meta["mu"], "ci_low": meta["ci_low"], "ci_high": meta["ci_high"]})
    write_csv(PHASE1A / "results" / "random_effects_summary.csv", meta_rows, list(meta_rows[0].keys()))
    write_csv(PHASE1A / "results" / "leave_one_cancer_out.csv", loo_rows, list(loo_rows[0].keys()))
    print("sensitivities", now(), flush=True)
    sens_rows = []
    for cancer in CANCERS:
        pack = packs[cancer]
        for transform in ("A", "B"):
            membership = memberships[cancer][transform]
            if transform == "A":
                dropped = []
                purity_ok = np.isfinite(pack["purity"])
                for i in range(len(pack["gene_ids"])):
                    if not pack["in_protein"][i] or not pack["in_cnv"][i]:
                        continue
                    mask = (
                        np.isfinite(pack["protein"][i]) & np.isfinite(pack["cnv"][i]) & purity_ok
                        & np.isfinite(pack["logG_A"][i]) & np.isfinite(pack["P_A"][i]) & ~pack["all_zero"][i]
                    )
                    idx = np.flatnonzero(mask)
                    if idx.size < 71:
                        continue
                    result = fit_nested(pack["logG_A"][i, idx], pack["P_A"][i, idx], pack["purity"][idx], pack["cnv"][i, idx], pack["protein"][i, idx], clips)
                    if result["fit_status"] == "OK":
                        dropped.append(result["partial_r2"])
                sens_rows.append({"sensitivity": "S1_drop_stored_all_zero", "cancer": cancer, "transformation": "A", "n_genes": len(dropped), "mean_partial_r2": float(np.mean(dropped)) if dropped else np.nan, "median_partial_r2": float(np.median(dropped)) if dropped else np.nan})
            rsem_values = sensitivity_fit(pack, membership, transform, "rsem", clips)
            sens_rows.append({"sensitivity": "S2_rsem", "cancer": cancer, "transformation": transform, "n_genes": len(rsem_values), "mean_partial_r2": float(np.mean(rsem_values)) if rsem_values else np.nan, "median_partial_r2": float(np.median(rsem_values)) if rsem_values else np.nan})
            est_values = sensitivity_fit(pack, membership, transform, "estimate", clips)
            sens_rows.append({"sensitivity": "S5_estimate", "cancer": cancer, "transformation": transform, "n_genes": len(est_values), "mean_partial_r2": float(np.mean(est_values)) if est_values else np.nan, "median_partial_r2": float(np.median(est_values)) if est_values else np.nan})
            # S6 and S7 need the stored VIF and residuals. Refit from primary rows.
            high_vif = []
            outlier_values = []
            rows = [row for row in primary_rows[transform] if row["cancer"] == cancer and row["fit_status"] == "OK"]
            gene_pos = {gene: i for i, gene in enumerate(pack["gene_ids"])}
            log_g = pack["logG_A"] if transform == "A" else pack["logG_B"]
            prop = pack["P_A"] if transform == "A" else pack["P_B"]
            for row in rows:
                gi = gene_pos[row["gene_id"]]
                idx = membership[gi]
                vif = float(row["cnv_vif"]) if row["cnv_vif"] else np.nan
                if np.isfinite(vif) and vif > 10:
                    value = _fit_without_cnv(log_g[gi, idx], prop[gi, idx], pack["purity"][idx], pack["protein"][gi, idx], clips)
                    if value is not None:
                        high_vif.append(value)
                result = fit_nested(log_g[gi, idx], prop[gi, idx], pack["purity"][idx], pack["cnv"][gi, idx], pack["protein"][gi, idx], clips)
                if result["fit_status"] != "OK":
                    continue
                keep = np.abs(result["pearson_m0"]) <= 4
                if keep.all():
                    outlier_values.append(result["partial_r2"])
                elif int(keep.sum()) >= 71:
                    value = partial_r2_fast(log_g[gi, idx][keep], prop[gi, idx][keep], pack["purity"][idx][keep], pack["cnv"][gi, idx][keep], pack["protein"][gi, idx][keep])
                    if value is not None:
                        outlier_values.append(value)
            sens_rows.append({"sensitivity": "S6_cnv_vif_gt_10", "cancer": cancer, "transformation": transform, "n_genes": len(high_vif), "mean_partial_r2": float(np.mean(high_vif)) if high_vif else np.nan, "median_partial_r2": float(np.median(high_vif)) if high_vif else np.nan})
            sens_rows.append({"sensitivity": "S7_m0_pearson_gt_4", "cancer": cancer, "transformation": transform, "n_genes": len(outlier_values), "mean_partial_r2": float(np.mean(outlier_values)) if outlier_values else np.nan, "median_partial_r2": float(np.median(outlier_values)) if outlier_values else np.nan})
        c_rows, c_partial, _c_mem = fit_rows(pack, "C", clips)
        c_ok = [float(row["partial_r2"]) for row in c_rows if row["fit_status"] == "OK"]
        sens_rows.append({"sensitivity": "S3_candidate_C", "cancer": cancer, "transformation": "C", "n_genes": len(c_ok), "mean_partial_r2": float(np.mean(c_ok)) if c_ok else np.nan, "median_partial_r2": float(np.median(c_ok)) if c_ok else np.nan})
    for transform in ("A", "B"):
        sens_rows.append({"sensitivity": "S4_clr", "cancer": "ALL", "transformation": transform, "n_genes": 0, "mean_partial_r2": "", "median_partial_r2": "", "note": "NOT_RUN_UNDERDEFINED"})
        sens_rows.append({"sensitivity": "S8_meta_without_OV", "cancer": "ALL_EXCEPT_OV", "transformation": transform, "n_genes": 9, "mean_partial_r2": no_ov_meta[transform]["mu"], "median_partial_r2": "", "ci_low": no_ov_meta[transform]["ci_low"], "ci_high": no_ov_meta[transform]["ci_high"]})
    write_csv(PHASE1A / "results" / "sensitivity_summary.csv", sens_rows, sorted({key for row in sens_rows for key in row}))
    print("spot checks", now(), flush=True)
    spot_rows = spot_check(packs, partials, memberships, expected_ids, tx_to_gene, principal)
    write_csv(PHASE1A / "qc" / "numerical_equivalence_checks.csv", spot_rows, list(spot_rows[0].keys()))
    loo_a = [row for row in loo_rows if row["transformation"] == "A"]
    loo_b = [row for row in loo_rows if row["transformation"] == "B"]
    science = classify_science(primary_meta["A"], primary_meta["B"], paired_meta["A"], paired_meta["B"], loo_a, loo_b)
    notes = [
        f"RSS clipping events: {len(clips)}",
        "Case matching is case-level.",
        "S4 CLR was not run: NOT_RUN_UNDERDEFINED.",
    ]
    if any(jack_by[cancer][transform]["unstable"] for cancer in CANCERS for transform in ("A", "B")):
        notes.append("At least one jackknife interval was flagged unstable.")
    write_reports(science, "EXECUTION_PASS", notes)
    write_narrative(summary_rows, paired_rows, meta_rows, sens_rows, science)
    (PHASE1A / "logs" / "run.log").write_text(f"start {STARTED}\nend {now()}\nscience {science}\nexecution EXECUTION_PASS\n")
    write_run_manifest("EXECUTION_PASS", science)
    hash_outputs()
    print("PHASE1A_DONE", science, now(), flush=True)


def spot_check(packs, partials, memberships, expected_ids, tx_to_gene, principal) -> list[dict]:
    rows = []
    for cancer in CANCERS:
        pack = packs[cancer]
        chosen = {}
        for transform in ("A", "B"):
            ok = np.flatnonzero(np.isfinite(partials[cancer][transform]))
            order = sorted(ok, key=lambda i: hashlib.sha256(pack["gene_ids"][i].encode()).hexdigest())
            chosen[transform] = order[:10]
        needed = {pack["gene_ids"][i] for transform in ("A", "B") for i in chosen[transform]}
        recomputed = recompute_selected(acquisition_path(cancer, "rna_isoform"), pack["case_ids"], expected_ids, needed, tx_to_gene, principal)
        for transform in ("A", "B"):
            log_key = "logG_A" if transform == "A" else "logG_B"
            p_key = "P_A" if transform == "A" else "P_B"
            for gi in chosen[transform]:
                gene = pack["gene_ids"][gi]
                idx = memberships[cancer][transform][gi]
                fresh_g, fresh_p = recomputed[gene][transform]
                g_ok = np.allclose(pack[log_key][gi, idx], fresh_g[idx], rtol=1e-8, atol=1e-8)
                p_ok = np.allclose(pack[p_key][gi, idx], fresh_p[idx], rtol=1e-8, atol=1e-8)
                again = partial_r2_fast(pack[log_key][gi, idx], pack[p_key][gi, idx], pack["purity"][idx], pack["cnv"][gi, idx], pack["protein"][gi, idx])
                r_ok = again is not None and abs(again - partials[cancer][transform][gi]) <= 1e-8
                rows.append({
                    "cancer": cancer,
                    "gene_id": gene,
                    "transformation": transform,
                    "n_complete": int(idx.size),
                    "logG_match": int(g_ok),
                    "proportion_match": int(p_ok),
                    "partial_r2_match": int(r_ok),
                })
                if not (g_ok and p_ok and r_ok):
                    raise NumericalFailure(f"spot check failed {cancer} {gene} {transform}")
    return rows


def recompute_selected(path, case_ids, expected_ids, genes, tx_to_gene, principal):
    n_cases = len(case_ids)
    acc = {gene: {"a": np.zeros(n_cases), "b": np.zeros(n_cases), "pa": np.zeros(n_cases), "pb": np.zeros(n_cases)} for gene in genes}
    names = header_names(path)
    columns = [{name: i for i, name in enumerate(names)}[case] for case in case_ids]
    rows = 0
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip(b"\r\n").split(b"\t")
            tx = parts[0].decode()
            rows += 1
            gene = tx_to_gene[tx]
            if gene not in acc:
                continue
            y = parse_numeric_row(parts, columns)
            power = np.exp2(y)
            acc[gene]["a"] += power
            acc[gene]["b"] += np.expm1(y * LOG2E)
            if principal.get(gene) == tx:
                acc[gene]["pa"] += power
                acc[gene]["pb"] += np.expm1(y * LOG2E)
    out = {}
    for gene, item in acc.items():
        out[gene] = {
            "A": (np.log2(item["a"]), item["pa"] / item["a"]),
            "B": (np.where(item["b"] > 0, np.log2(item["b"]), np.nan), np.where(item["b"] > 0, item["pb"] / item["b"], np.nan)),
        }
    return out


def write_narrative(summary_rows, paired_rows, meta_rows, sens_rows, science) -> None:
    def cell(value):
        return fmt(value) if not isinstance(value, str) else value

    lines = ["# Phase 1A primary results", "", f"Scientific classification: {science}", ""]
    lines.append("Cancer-level means are unweighted across successful genes. Intervals are delete-one-case jackknife intervals on the raw partial-R² scale.")
    lines.append("")
    lines.append("| Cancer | Transform | Successful genes | Mean partial R² | 95% interval | Proportion ≥ 0.02 |")
    lines.append("| --- | --- | ---: | ---: | --- | ---: |")
    for row in summary_rows:
        lines.append(f"| {row['cancer']} | {row['transformation']} | {row['n_success']} | {cell(row['mean_partial_r2'])} | {cell(row['ci_low'])}, {cell(row['ci_high'])} | {cell(row['prop_ge_0.02'])} |")
    lines.append("")
    lines.append("Cross-cancer REML / Hartung–Knapp summaries are in `results/random_effects_summary.csv`. OV remains in its own rows above.")
    (PHASE1A / "reports" / "phase1a_primary_results.md").write_text("\n".join(lines) + "\n")
    sens = ["# Phase 1A sensitivity results", "", "Every prespecified sensitivity was attempted. S4 is recorded as NOT_RUN_UNDERDEFINED.", ""]
    sens.append("| Sensitivity | Cancer | Transform | Genes | Mean partial R² |")
    sens.append("| --- | --- | --- | ---: | ---: |")
    for row in sens_rows:
        sens.append(f"| {row.get('sensitivity','')} | {row.get('cancer','')} | {row.get('transformation','')} | {row.get('n_genes','')} | {cell(row.get('mean_partial_r2',''))} |")
    (PHASE1A / "reports" / "phase1a_sensitivity_results.md").write_text("\n".join(sens) + "\n")
    exe = ["# Phase 1A execution report", "", "Execution status: EXECUTION_PASS", "", f"Started {STARTED}. Completed {now()}.", "", f"Scientific result: {science}", "", "Gates passed: input hashes, NUL scans, 1,027-case reconstruction, GENCODE v34 principal-transcript counts, and the Literal B pre-fit eligibility counts.", ""]
    exe.append(f"RSS values clipped for floating-point negativity are listed by count in the limitations report.")
    exe.append("")
    exe.append("PHASE 1B AUTHORIZATION = FALSE")
    (PHASE1A / "reports" / "phase1a_execution_report.md").write_text("\n".join(exe) + "\n")


def write_run_manifest(status: str, science: str = "") -> None:
    write_csv(
        PHASE1A / "manifests" / "run_manifest.csv",
        [
            {"item": "start_utc", "value": STARTED},
            {"item": "end_utc", "value": now()},
            {"item": "execution_status", "value": status},
            {"item": "scientific_result", "value": science},
            {"item": "python", "value": sys.version.split()[0]},
            {"item": "numpy", "value": np.__version__},
            {"item": "threads", "value": "1"},
            {"item": "seed", "value": "20261114"},
            {"item": "phase_1b_authorization", "value": "FALSE"},
        ],
        ["item", "value"],
    )


def write_closed_decision(reason: str) -> None:
    path = PHASE1A / "reports" / "phase1a_final_decision.md"
    if path.exists():
        return
    path.write_text(
        "# Phase 1A final decision\n\n"
        "## A. Execution status\n\nEXECUTION_FAIL\n\n"
        "## B. Scientific result\n\nTECHNICAL_NO_RESULT\n\n"
        "## C. Evidence boundary\n\n"
        f"The frozen primary models did not finish. {reason}\n\n"
        "The 1,027 tumors are matched at case level. No isoform-specific protein or functional-loss claim is made.\n\n"
        "## D. Next-phase authorization\n\nPHASE 1B AUTHORIZATION = FALSE\n"
    )
    write_run_manifest("EXECUTION_FAIL", "TECHNICAL_NO_RESULT")


def hash_outputs() -> None:
    lines = []
    for path in sorted(PHASE1A.rglob("*")):
        if not path.is_file() or path.name.startswith("._") or path.name == "SHA256SUMS":
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        digest, _, _, _ = sha256_and_nul(path)
        lines.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
    (PHASE1A / "SHA256SUMS").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as exc:
        write_closed_decision(str(exc))
        raise
    except Exception:
        traceback.print_exc()
        (PHASE1A / "logs" / "run_error.txt").write_text(traceback.format_exc())
        write_closed_decision("An uncaught exception stopped the run. See logs/run_error.txt.")
        raise
