#!/usr/bin/env python3
"""Technical audit of the two literal isoform inverses.

Candidate A is X = 2^Y, so a stored zero becomes 1.
Candidate B is X = 2^Y - 1, so a stored zero stays 0.
Neither inverse is chosen from a protein association. Phase 0 files are not modified.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from array import array
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1].parent / "phase0_5" / "code"))
from rebuild_cohort import CANCERS, load_annotation, path_for  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUTM = ROOT / "phase0_6" / "manifests"
OUTL = ROOT / "phase0_6" / "logs"
# GENCODE v34 principal transcript present in the isoform matrix, gene has >=2 rows.
PRINCIPAL_MULTI_IN_MATRIX = 10795
for folder in (OUTM, OUTL):
    folder.mkdir(parents=True, exist_ok=True)

LUT = tuple(2.0 ** (centi / 100.0) for centi in range(0, 4001))
LOG2 = math.log2


def centi(token: bytes) -> int | None:
    token = token.strip()
    if token in (b"", b"NA", b"NaN", b"nan"):
        return None
    if token in (b"0", b"0.0", b"0.00"):
        return 0
    if b"." in token:
        whole, frac = token.split(b".", 1)
        if len(frac) > 2:
            return int(round(float(token) * 100))
        frac = (frac + b"00")[:2]
        sign = -1 if whole.startswith(b"-") else 1
        if whole.startswith(b"-"):
            whole = whole[1:]
        return sign * (int(whole or b"0") * 100 + int(frac))
    return int(token) * 100


def quantile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    pos = q * (len(ordered) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return ordered[lo]
    weight = pos - lo
    return ordered[lo] * (1 - weight) + ordered[hi] * weight


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    n = len(values)
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = average
        i = j + 1
    return out


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    num = dx = dy = 0.0
    for a, b in zip(x, y):
        da = a - mx
        db = b - my
        num += da * db
        dx += da * da
        dy += db * db
    if dx <= 0 or dy <= 0:
        return float("nan")
    return num / math.sqrt(dx * dy)


def spearman(x: list[float], y: list[float]) -> float:
    if len(x) < 3:
        return float("nan")
    return pearson(ranks(x), ranks(y))


def fmt(value: float) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return f"{value:.6g}"


def load_rsem(path: Path, samples: list[str]) -> dict[str, array]:
    with path.open("rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(b"\t")
        names = [part.decode() for part in header[1:]]
        locate = {name: i for i, name in enumerate(names)}
        columns = [locate.get(sample) for sample in samples]
        out: dict[str, array] = {}
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip(b"\r\n").split(b"\t")
            row = array("d")
            for col in columns:
                if col is None or col + 1 >= len(cols):
                    row.append(float("nan"))
                    continue
                token = cols[col + 1].strip()
                if token in (b"", b"NA", b"NaN", b"nan"):
                    row.append(float("nan"))
                else:
                    row.append(float(token))
            out[cols[0].decode()] = row
    return out


def ntx_bin(k: int) -> str:
    if k <= 1:
        return "1"
    if k == 2:
        return "2"
    if k <= 5:
        return "3-5"
    if k <= 10:
        return "6-10"
    return "11+"


def subset_corr(r_values: list[float], g_values: list[float], mask: list[bool]) -> tuple[float, float, int]:
    rr = [a for a, keep in zip(r_values, mask) if keep]
    gg = [b for b, keep in zip(g_values, mask) if keep]
    return spearman(gg, rr), pearson(gg, rr), len(rr)


def source_matrix(cancer: str, layer: str) -> Path | None:
    original = path_for(cancer, layer)
    if original is None:
        return None
    repaired = ROOT / "phase0_6" / "repaired" / cancer / original.name
    return repaired if repaired.exists() else original


def audit_cancer(cancer: str, tx_to_gene: dict[str, str], tx_bare: dict[str, str], principal: dict[str, str]) -> dict:
    isoform = source_matrix(cancer, "rna_isoform")
    gene_path = source_matrix(cancer, "rna_gene")
    if isoform is None or gene_path is None:
        raise SystemExit(f"missing matrix for {cancer}")
    with isoform.open("rb") as handle:
        iso_names = [part.decode() for part in handle.readline().rstrip(b"\r\n").split(b"\t")[1:]]
    with gene_path.open("rb") as handle:
        gene_names = {part.decode() for part in handle.readline().rstrip(b"\r\n").split(b"\t")[1:]}
    samples = [name for name in iso_names if name in gene_names]
    n = len(samples)
    rsem = load_rsem(gene_path, samples)
    principal_tx = {tx: gene for gene, tx in principal.items()}

    genes: list[str] = []
    txs: list[str] = []
    values = array("H")
    unmapped = offgrid = duplicate = rows = 0
    seen: set[str] = set()
    sample_pos = {name: i for i, name in enumerate(iso_names)}
    columns = [sample_pos[sample] + 1 for sample in samples]
    with isoform.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip(b"\r\n").split(b"\t")
            tx = cols[0].decode()
            rows += 1
            if tx in seen:
                duplicate += 1
            seen.add(tx)
            gene = tx_to_gene.get(tx) or tx_bare.get(tx.split(".", 1)[0])
            if not gene:
                unmapped += 1
                continue
            genes.append(gene)
            txs.append(tx)
            width = len(cols)
            for col in columns:
                token = cols[col] if col < width else b"NA"
                coded = centi(token)
                if coded is None:
                    raise SystemExit(f"{cancer} non-numeric isoform token {token!r}")
                if coded < 0 or coded >= len(LUT):
                    raise SystemExit(f"{cancer} isoform centi out of range {coded}")
                if b"." in token and len(token.split(b".", 1)[1]) > 2:
                    offgrid += 1
                values.append(coded)
            if rows % 40000 == 0:
                print(f"  {cancer} isoform rows {rows}", flush=True)
    n_tx = len(txs)
    if len(values) != n_tx * n:
        raise SystemExit(f"{cancer} value grid {len(values)} != {n_tx}*{n}")

    groups: dict[str, list[int]] = defaultdict(list)
    for i, gene in enumerate(genes):
        groups[gene].append(i)
    tx_index = {tx: i for i, tx in enumerate(txs)}

    structural_one = structural_multi = 0
    one_b = multi_b = 0
    always_one = always_multi = 0
    all_zero = dominant_defined = dominant_mismatch = 0
    principal_defined = delta_gt = 0
    delta_sum = entropy_a_sum = entropy_b_sum = 0.0
    entropy_a_n = entropy_b_n = 0
    p_principal_a: list[float] = []
    p_principal_b: list[float] = []
    p_max_a: list[float] = []
    p_max_b: list[float] = []

    r_a: list[float] = []
    g_a: list[float] = []
    k_a: list[int] = []
    max_a: list[int] = []
    r_b: list[float] = []
    g_b: list[float] = []
    k_b: list[int] = []
    max_b: list[int] = []
    rsem_zero_when_iso_zero = iso_zero_with_rsem = 0

    for gene, idxs in groups.items():
        k = len(idxs)
        if k == 1:
            structural_one += 1
        else:
            structural_multi += 1
        chosen = principal.get(gene)
        principal_i = None
        if chosen is not None:
            principal_i = tx_index.get(chosen)
            if principal_i is None:
                bare_chosen = chosen.split(".", 1)[0]
                for i in idxs:
                    if txs[i].split(".", 1)[0] == bare_chosen and txs[i] not in tx_to_gene:
                        principal_i = i
                        break
        rsem_row = rsem.get(gene)
        sample_one = sample_multi = 0
        for s in range(n):
            gA = gB = 0.0
            npos = 0
            max_c = -1
            for i in idxs:
                coded = values[i * n + s]
                xA = LUT[coded]
                gA += xA
                if coded:
                    gB += xA - 1.0
                    npos += 1
                if coded > max_c:
                    max_c = coded
            if npos == 0:
                all_zero += 1
            elif npos == 1:
                one_b += 1
                sample_one += 1
            else:
                multi_b += 1
                sample_multi += 1
            # Both inverses are strictly increasing, so the first maximum of Y
            # is the first maximum of X. A mismatch would mean the lookup failed.
            max_xA = max_xB = -1.0
            for i in idxs:
                coded = values[i * n + s]
                xA = LUT[coded]
                xB = 0.0 if coded == 0 else xA - 1.0
                if xA > max_xA:
                    max_xA = xA
                if xB > max_xB:
                    max_xB = xB
            y_at_first = LUT[max_c] if max_c >= 0 else 0.0
            if npos > 0 and (abs(max_xA - y_at_first) > 1e-6 or abs(max_xB - (y_at_first - 1.0)) > 1e-6):
                dominant_mismatch += 1
            if npos > 0 and sum(1 for i in idxs if values[i * n + s] == max_c) == 1:
                dominant_defined += 1
            hA = hB = 0.0
            for i in idxs:
                coded = values[i * n + s]
                pA = LUT[coded] / gA
                hA -= pA * LOG2(pA)
                if coded and gB > 0:
                    pB = (LUT[coded] - 1.0) / gB
                    hB -= pB * LOG2(pB)
            entropy_a_sum += hA
            entropy_a_n += 1
            p_max_a.append(max_xA / gA)
            if gB > 0:
                entropy_b_sum += hB
                entropy_b_n += 1
                p_max_b.append(max_xB / gB)
            if principal_i is not None and k >= 2:
                coded = values[principal_i * n + s]
                pA = LUT[coded] / gA
                p_principal_a.append(pA)
                if gB > 0:
                    pB = 0.0 if coded == 0 else (LUT[coded] - 1.0) / gB
                    p_principal_b.append(pB)
                    principal_defined += 1
                    diff = abs(pA - pB)
                    delta_sum += diff
                    if diff > 0.05:
                        delta_gt += 1
            if rsem_row is not None and math.isfinite(rsem_row[s]):
                r_value = rsem_row[s]
                r_a.append(r_value)
                g_a.append(LOG2(gA))
                k_a.append(k)
                max_a.append(max_c if max_c > 0 else 0)
                if npos == 0:
                    iso_zero_with_rsem += 1
                    if r_value == 0.0:
                        rsem_zero_when_iso_zero += 1
                if gB > 0:
                    r_b.append(r_value)
                    g_b.append(LOG2(gB))
                    k_b.append(k)
                    max_b.append(max_c)
        if sample_one == n and n:
            always_one += 1
        if sample_multi == n and n:
            always_multi += 1

    def corr_table(r_values: list[float], g_values: list[float], k_values: list[int], max_values: list[int]) -> dict:
        sp, pe, n_pairs = subset_corr(r_values, g_values, [True] * len(r_values))
        residuals = [g - r for g, r in zip(g_values, r_values)]
        ordered_r = sorted(r_values)
        if ordered_r:
            q1 = quantile(ordered_r, 0.25)
            q2 = quantile(ordered_r, 0.50)
            q3 = quantile(ordered_r, 0.75)
        else:
            q1 = q2 = q3 = float("nan")
        bins = {
            "q1": [r <= q1 for r in r_values],
            "q2": [q1 < r <= q2 for r in r_values],
            "q3": [q2 < r <= q3 for r in r_values],
            "q4": [r > q3 for r in r_values],
            "ntx_1": [k == 1 for k in k_values],
            "ntx_2": [k == 2 for k in k_values],
            "ntx_3_5": [3 <= k <= 5 for k in k_values],
            "ntx_6_10": [6 <= k <= 10 for k in k_values],
            "ntx_11": [k >= 11 for k in k_values],
            "low_max_y": [m <= 100 for m in max_values],
            "high_max_y": [m > 100 for m in max_values],
        }
        out = {
            "spearman": sp,
            "pearson": pe,
            "n": n_pairs,
            "median_residual": quantile(residuals, 0.5) if residuals else float("nan"),
            "residual_p05": quantile(residuals, 0.05) if residuals else float("nan"),
            "residual_p95": quantile(residuals, 0.95) if residuals else float("nan"),
        }
        for name, mask in bins.items():
            sp_b, pe_b, n_b = subset_corr(r_values, g_values, mask)
            out[f"spearman_{name}"] = sp_b
            out[f"pearson_{name}"] = pe_b
            out[f"n_{name}"] = n_b
        return out

    print(f"  {cancer} correlations A={len(r_a)} B={len(r_b)}", flush=True)
    corr_a = corr_table(r_a, g_a, k_a, max_a)
    corr_b = corr_table(r_b, g_b, k_b, max_b)
    return {
        "cancer": cancer,
        "samples_shared": n,
        "isoform_rows_read": rows,
        "isoform_rows_mapped": n_tx,
        "unmapped": unmapped,
        "duplicate_transcript_ids": duplicate,
        "offgrid_cells": offgrid,
        "genes": len(groups),
        "structural_one": structural_one,
        "structural_multi": structural_multi,
        "gene_samples": len(groups) * n,
        "all_zero_gene_samples": all_zero,
        "zero_denominator_A": 0,
        "zero_denominator_B": all_zero,
        "detected_one_B": one_b,
        "detected_multi_B": multi_b,
        "detected_one_A": structural_one * n,
        "detected_multi_A": structural_multi * n,
        "genes_always_one_detected_B": always_one,
        "genes_always_multi_detected_B": always_multi,
        "dominant_defined": dominant_defined,
        "dominant_mismatch": dominant_mismatch,
        "principal_pairs": principal_defined,
        "delta_gt_0.05": delta_gt,
        "delta_mean": delta_sum / principal_defined if principal_defined else float("nan"),
        "delta_fraction_gt_0.05": delta_gt / principal_defined if principal_defined else float("nan"),
        "entropy_mean_A_all": entropy_a_sum / entropy_a_n if entropy_a_n else float("nan"),
        "entropy_mean_B_defined": entropy_b_sum / entropy_b_n if entropy_b_n else float("nan"),
        "p_principal_A_p05": quantile(p_principal_a, 0.05),
        "p_principal_A_p50": quantile(p_principal_a, 0.50),
        "p_principal_A_p95": quantile(p_principal_a, 0.95),
        "p_principal_B_p05": quantile(p_principal_b, 0.05),
        "p_principal_B_p50": quantile(p_principal_b, 0.50),
        "p_principal_B_p95": quantile(p_principal_b, 0.95),
        "p_max_A_p50": quantile(p_max_a, 0.50),
        "p_max_B_p50": quantile(p_max_b, 0.50),
        "iso_zero_with_rsem": iso_zero_with_rsem,
        "rsem_zero_when_iso_zero": rsem_zero_when_iso_zero,
        "rsem_zero_fraction_when_iso_zero": (
            rsem_zero_when_iso_zero / iso_zero_with_rsem if iso_zero_with_rsem else float("nan")
        ),
        "corr_A": corr_a,
        "corr_B": corr_b,
    }


def write_tables(results: list[dict]) -> None:
    comp_fields = [
        "cancer", "inverse", "formula", "zero_stored_becomes", "gene_samples",
        "zero_denominator_gene_samples", "zero_denominator_fraction",
        "genes_with_one_matrix_row", "genes_with_multiple_matrix_rows",
        "gene_samples_one_detected", "gene_samples_multiple_detected",
        "genes_with_one_detected_in_every_sample",
        "principal_contrast_gene_samples", "principal_quantile_observations",
        "principal_proportion_p05", "principal_proportion_p50", "principal_proportion_p95",
        "max_proportion_p50", "entropy_mean_bits",
        "dominant_unique_and_positive", "dominant_identity_mismatches",
        "mean_abs_delta_principal_proportion", "fraction_abs_delta_principal_gt_0.05",
    ]
    with (OUTM / "transformation_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=comp_fields)
        writer.writeheader()
        for row in results:
            base = row["gene_samples"]
            writer.writerow({
                "cancer": row["cancer"], "inverse": "A", "formula": "2^Y", "zero_stored_becomes": 1,
                "gene_samples": base, "zero_denominator_gene_samples": row["zero_denominator_A"],
                "zero_denominator_fraction": 0,
                "genes_with_one_matrix_row": row["structural_one"],
                "genes_with_multiple_matrix_rows": row["structural_multi"],
                "gene_samples_one_detected": row["detected_one_A"],
                "gene_samples_multiple_detected": row["detected_multi_A"],
                "genes_with_one_detected_in_every_sample": row["structural_one"],
                "principal_contrast_gene_samples": row["principal_pairs"],
                "principal_quantile_observations": PRINCIPAL_MULTI_IN_MATRIX * row["samples_shared"],
                "principal_proportion_p05": fmt(row["p_principal_A_p05"]),
                "principal_proportion_p50": fmt(row["p_principal_A_p50"]),
                "principal_proportion_p95": fmt(row["p_principal_A_p95"]),
                "max_proportion_p50": fmt(row["p_max_A_p50"]),
                "entropy_mean_bits": fmt(row["entropy_mean_A_all"]),
                "dominant_unique_and_positive": row["dominant_defined"],
                "dominant_identity_mismatches": row["dominant_mismatch"],
                "mean_abs_delta_principal_proportion": fmt(row["delta_mean"]),
                "fraction_abs_delta_principal_gt_0.05": fmt(row["delta_fraction_gt_0.05"]),
            })
            writer.writerow({
                "cancer": row["cancer"], "inverse": "B", "formula": "2^Y-1", "zero_stored_becomes": 0,
                "gene_samples": base, "zero_denominator_gene_samples": row["zero_denominator_B"],
                "zero_denominator_fraction": fmt(row["zero_denominator_B"] / base if base else float("nan")),
                "genes_with_one_matrix_row": row["structural_one"],
                "genes_with_multiple_matrix_rows": row["structural_multi"],
                "gene_samples_one_detected": row["detected_one_B"],
                "gene_samples_multiple_detected": row["detected_multi_B"],
                "genes_with_one_detected_in_every_sample": row["genes_always_one_detected_B"],
                "principal_contrast_gene_samples": row["principal_pairs"],
                "principal_quantile_observations": row["principal_pairs"],
                "principal_proportion_p05": fmt(row["p_principal_B_p05"]),
                "principal_proportion_p50": fmt(row["p_principal_B_p50"]),
                "principal_proportion_p95": fmt(row["p_principal_B_p95"]),
                "max_proportion_p50": fmt(row["p_max_B_p50"]),
                "entropy_mean_bits": fmt(row["entropy_mean_B_defined"]),
                "dominant_unique_and_positive": row["dominant_defined"],
                "dominant_identity_mismatches": row["dominant_mismatch"],
                "mean_abs_delta_principal_proportion": fmt(row["delta_mean"]),
                "fraction_abs_delta_principal_gt_0.05": fmt(row["delta_fraction_gt_0.05"]),
            })
    audit_fields = ["cancer", "samples_shared", "isoform_rows_mapped", "unmapped", "genes",
                    "structural_one_transcript_genes", "structural_multi_transcript_genes",
                    "all_zero_gene_samples", "dominant_identity_mismatches",
                    "fraction_abs_delta_principal_gt_0.05", "mean_abs_delta_principal_proportion",
                    "rsem_zero_fraction_when_all_isoform_rows_are_zero"]
    for label in ("A", "B"):
        audit_fields.extend([
            f"spearman_{label}", f"pearson_{label}", f"n_{label}", f"median_residual_{label}",
            f"residual_p05_{label}", f"residual_p95_{label}",
        ])
        for name in ("q1", "q2", "q3", "q4", "ntx_1", "ntx_2", "ntx_3_5", "ntx_6_10", "ntx_11", "low_max_y", "high_max_y"):
            audit_fields.append(f"spearman_{label}_{name}")
            audit_fields.append(f"n_{label}_{name}")
    with (OUTM / "cancer_specific_quantification_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields)
        writer.writeheader()
        for row in results:
            out = {
                "cancer": row["cancer"],
                "samples_shared": row["samples_shared"],
                "isoform_rows_mapped": row["isoform_rows_mapped"],
                "unmapped": row["unmapped"],
                "genes": row["genes"],
                "structural_one_transcript_genes": row["structural_one"],
                "structural_multi_transcript_genes": row["structural_multi"],
                "all_zero_gene_samples": row["all_zero_gene_samples"],
                "dominant_identity_mismatches": row["dominant_mismatch"],
                "fraction_abs_delta_principal_gt_0.05": fmt(row["delta_fraction_gt_0.05"]),
                "mean_abs_delta_principal_proportion": fmt(row["delta_mean"]),
                "rsem_zero_fraction_when_all_isoform_rows_are_zero": fmt(row["rsem_zero_fraction_when_iso_zero"]),
            }
            for label, corr in (("A", row["corr_A"]), ("B", row["corr_B"])):
                out[f"spearman_{label}"] = fmt(corr["spearman"])
                out[f"pearson_{label}"] = fmt(corr["pearson"])
                out[f"n_{label}"] = corr["n"]
                out[f"median_residual_{label}"] = fmt(corr["median_residual"])
                out[f"residual_p05_{label}"] = fmt(corr["residual_p05"])
                out[f"residual_p95_{label}"] = fmt(corr["residual_p95"])
                for name in ("q1", "q2", "q3", "q4", "ntx_1", "ntx_2", "ntx_3_5", "ntx_6_10", "ntx_11", "low_max_y", "high_max_y"):
                    out[f"spearman_{label}_{name}"] = fmt(corr[f"spearman_{name}"])
                    out[f"n_{label}_{name}"] = corr[f"n_{name}"]
            writer.writerow(out)
    (OUTL / "quantification_audit.json").write_text(json.dumps(results, indent=2) + "\n")


def main() -> None:
    selected = sys.argv[1:] or CANCERS
    print("loading GENCODE v34", flush=True)
    tx_to_gene, tx_bare, principal, counts = load_annotation()
    print(counts, flush=True)
    results = []
    for cancer in selected:
        print(cancer, flush=True)
        results.append(audit_cancer(cancer, tx_to_gene, tx_bare, principal))
    write_tables(results)
    print("wrote quantification tables", flush=True)


if __name__ == "__main__":
    main()
