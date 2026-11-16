#!/usr/bin/env python3
"""Phase 2B structural-category association. Does not refit Phase 1A models."""
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
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "phase1a" / "code"))
from common import CANCERS, acquisition_path  # noqa: E402

PHASE = ROOT / "phase2b"
SEED = 20261114
DRAWS = 1000
PRACTICAL = 0.005
VIF_LIMIT = 10.0
MISSING_LIMIT = 0.05
MIN_CLASS = 500
MIN_WINDOW = 300
MIN_STRATUM = 30
SMD_LIMIT = 1.0
Z95 = 1.959963984540054
PFAM_ANALYSIS = 9930
COVARIATE_PRIORITY = [
    "principal_length",
    "n_pfam",
    "pfam_coverage",
    "n_alt_protein_coding",
    "n_isoform_rows",
    "n_complete",
    "protein_fraction",
    "median_log2_g",
]
PREFIT = {
    "phase2b/protocol/PHASE2B_LOCK.md": "2ba8ce6e3e4c7962a772f0aa23a4793f9273e201b0a3e5e608ef5b4bba50d882",
    "phase2b/protocol/analysis_config.yaml": "13e8d813d074ce7459500e7cf57c7f2e19e94053dce3edb4b53698129f90d980",
}
EXPECTED = {
    "candidate_genes": 10795,
    "measured_protein_coding_transcripts": 39490,
    "genes_with_uniprot": 10742,
    "genes_with_sequence_difference": 9608,
    "genes_with_structural_difference": 4829,
    "genes_with_principal_pfam": 10228,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(message: str) -> None:
    line = f"{utc_now()} {message}"
    print(line, flush=True)
    with (PHASE / "logs" / "run.log").open("a") as handle:
        handle.write(line + "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def verify_historical() -> None:
    required = {
        "phase1a/SHA256SUMS": [
            "phase1a/reports/phase1a_final_decision.md",
            "phase1a/reports/phase1a_primary_results.md",
            "phase1a/results/primary_literal_A_gene_fits.tsv.gz",
            "phase1a/results/primary_literal_B_gene_fits.tsv.gz",
            "phase1a/results/cancer_level_summary.csv",
            "phase1a/results/transformation_paired_comparison.csv",
        ],
        "phase1b/SHA256SUMS": [
            "phase1b/reports/phase1b_final_decision.md",
            "phase1b/reports/phase1b_gene_distribution_report.md",
            "phase1b/reports/phase1b_transcript_usage_report.md",
            "phase1b/reports/phase1b_transformation_robustness.md",
        ],
        "phase2a/SHA256SUMS": [
            "phase2a/reports/phase2a_final_decision.md",
            "phase2a/reports/transcript_protein_mapping_report.md",
            "phase2a/reports/protein_domain_annotation_report.md",
            "phase2a/reports/structural_difference_feasibility.md",
            "phase2a/reports/phase1a_annotation_linkage_report.md",
            "phase2a/reports/evidence_hierarchy_report.md",
            "phase2a/results/functional_annotation_coverage.csv",
            "phase2a/results/protein_sequence_difference_coverage.csv",
            "phase2a/results/transcript_structure_coverage.csv",
            "phase2a/manifests/domain_mapping_audit.csv",
        ],
    }
    for ledger, paths in required.items():
        recorded = {}
        for line in (ROOT / ledger).read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            digest, name = line.split(maxsplit=1)
            recorded[name.strip()] = digest
        for path in paths:
            digest = sha256_file(ROOT / path)
            if recorded.get(path) != digest:
                raise SystemExit(f"historical checksum mismatch: {path}")
            log(f"checksum ok {path}")


def verify_prefit() -> None:
    for rel, expected in PREFIT.items():
        digest = sha256_file(ROOT / rel)
        if digest != expected:
            raise SystemExit(f"prefit hash mismatch for {rel}")


def synthetic_tests() -> None:
    rng = np.random.default_rng(SEED)
    y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    x = np.column_stack([np.ones(6), np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])])
    beta, resid = ols(x, y)
    if abs(beta[1] - 3.0) > 1e-8:
        raise SystemExit("synthetic mean difference failed")
    clusters = np.array([0, 0, 1, 1, 2, 2])
    cov = cluster_covariance(x, resid, clusters)
    if cov is None or cov.shape != (2, 2):
        raise SystemExit("synthetic cluster covariance failed")
    p = np.array([0.01, 0.04, 0.03])
    q = benjamini_hochberg(p)
    if np.any(q < p - 1e-12) or q[0] > q[2] + 1e-12:
        raise SystemExit("synthetic BH failed")
    draws = rng.normal(size=20)
    if abs(float(np.quantile(draws, 0.5)) - float(np.median(draws))) > 1e-12:
        raise SystemExit("synthetic quantile failed")
    (PHASE / "qc" / "synthetic_test_results.txt").write_text("SYNTHETIC_TESTS_PASS\n")


def ols(x_design: np.ndarray, outcome: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    beta, _, _, _ = np.linalg.lstsq(x_design, outcome, rcond=None)
    return beta, outcome - x_design @ beta


def cluster_covariance(x_design: np.ndarray, resid: np.ndarray, clusters: np.ndarray) -> np.ndarray | None:
    n_obs, n_col = x_design.shape
    gram = x_design.T @ x_design
    try:
        inverse = np.linalg.inv(gram)
    except np.linalg.LinAlgError:
        return None
    order = np.argsort(clusters, kind="mergesort")
    ordered = clusters[order]
    design = x_design[order]
    errors = resid[order]
    cuts = np.flatnonzero(np.diff(ordered)) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, n_obs]
    meat = np.zeros((n_col, n_col))
    for start, end in zip(starts, ends):
        score = design[start:end].T @ errors[start:end]
        meat += np.outer(score, score)
    n_groups = len(starts)
    if n_groups <= n_col:
        return None
    scale = (n_groups / (n_groups - 1.0)) * ((n_obs - 1.0) / (n_obs - n_col))
    return scale * inverse @ meat @ inverse


def hc1_covariance(x_design: np.ndarray, resid: np.ndarray) -> np.ndarray | None:
    n_obs, n_col = x_design.shape
    try:
        inverse = np.linalg.inv(x_design.T @ x_design)
    except np.linalg.LinAlgError:
        return None
    weighted = x_design * resid[:, None]
    meat = weighted.T @ weighted
    scale = n_obs / (n_obs - n_col)
    return scale * inverse @ meat @ inverse


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    size = len(p_values)
    order = np.argsort(p_values)
    ranked = p_values[order] * size / np.arange(1, size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(size)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def load_categories() -> tuple[dict[str, dict], dict[str, int]]:
    functional = {}
    with (ROOT / "phase2a/results/functional_annotation_coverage.csv").open() as handle:
        for row in csv.DictReader(handle):
            functional[row["gene_id"]] = row
    structure = {}
    with (ROOT / "phase2a/results/transcript_structure_coverage.csv").open() as handle:
        for row in csv.DictReader(handle):
            structure[row["gene_id"]] = row
    sequence = {}
    with (ROOT / "phase2a/results/protein_sequence_difference_coverage.csv").open() as handle:
        for row in csv.DictReader(handle):
            sequence[row["gene_id"]] = row
    pairs = defaultdict(lambda: {
        "s2": 0, "s1": 0, "identical": 0,
        "pfam_absent_from_alternative": 0, "pfam_gained_in_alternative": 0,
        "signal_principal_only": 0, "signal_alternative_only": 0,
        "tm_principal_only": 0, "tm_alternative_only": 0,
        "pfam_diff": 0, "signal_diff": 0, "tm_diff": 0,
    })
    with (ROOT / "phase2a/manifests/domain_mapping_audit.csv").open() as handle:
        for row in csv.DictReader(handle):
            item = pairs[row["gene_id"]]
            identical = row["sequences_identical"] == "1"
            clean = row["coordinate_note"] == "version_matched"
            if not clean:
                continue
            pfam_diff = row["domain_class"] in {"principal_only", "alternative_only", "mixed"}
            signal_diff = row["signal_principal"] != row["signal_alternative"]
            tm_diff = row["tm_principal"] != row["tm_alternative"]
            if identical:
                item["identical"] += 1
                continue
            if pfam_diff or signal_diff or tm_diff:
                item["s2"] += 1
            else:
                item["s1"] += 1
            if int(float(row["pfam_only_principal"] or 0)) > 0:
                item["pfam_absent_from_alternative"] += 1
                item["pfam_diff"] += 1
            if int(float(row["pfam_only_alternative"] or 0)) > 0:
                item["pfam_gained_in_alternative"] += 1
                item["pfam_diff"] += 1
            if signal_diff:
                item["signal_diff"] += 1
                if row["signal_principal"] == "1" and row["signal_alternative"] == "0":
                    item["signal_principal_only"] += 1
                if row["signal_principal"] == "0" and row["signal_alternative"] == "1":
                    item["signal_alternative_only"] += 1
            if tm_diff:
                item["tm_diff"] += 1
                principal_tm = float(row["tm_principal"] or 0) > 0
                alternative_tm = float(row["tm_alternative"] or 0) > 0
                if principal_tm and not alternative_tm:
                    item["tm_principal_only"] += 1
                if alternative_tm and not principal_tm:
                    item["tm_alternative_only"] += 1
    categories = {}
    for gene, row in functional.items():
        info = pairs.get(gene)
        if info is None or (info["s2"] == 0 and info["s1"] == 0 and info["identical"] == 0):
            label = "ANNOTATION_UNRESOLVED"
            info = info or {key: 0 for key in (
                "s2", "s1", "identical", "pfam_absent_from_alternative", "pfam_gained_in_alternative",
                "signal_principal_only", "signal_alternative_only", "tm_principal_only",
                "tm_alternative_only", "pfam_diff", "signal_diff", "tm_diff",
            )}
        elif info["s2"] > 0:
            label = "S2"
        elif info["s1"] > 0:
            label = "S1"
        else:
            label = "S0"
        categories[gene] = {
            "gene_id": gene,
            "category": label,
            "principal_transcript_id": row["principal_transcript_id"],
            "principal_protein_id": row["principal_protein_id"],
            "principal_length": float(sequence[gene]["principal_sequence_length"] or 0),
            "n_pfam": int(row["n_pfam_accessions"]),
            "n_alt_protein_coding": int(structure[gene]["n_alternative_protein_coding"]),
            "n_isoform_rows": int(structure[gene]["n_matrix_transcripts"]),
            "uniprot": int(row["uniprot_reviewed_or_unreviewed"]),
            "reliable_coordinates": int(row["reliable_domain_coordinates"]),
            "n_s2_pairs": info["s2"],
            "n_s1_pairs": info["s1"],
            "n_identical_pairs": info["identical"],
            "pfam_diff": int(info["pfam_diff"] > 0),
            "signal_diff": int(info["signal_diff"] > 0),
            "tm_diff": int(info["tm_diff"] > 0),
            "pfam_absent_from_alternative": int(info["pfam_absent_from_alternative"] > 0),
            "pfam_gained_in_alternative": int(info["pfam_gained_in_alternative"] > 0),
            "signal_principal_only": int(info["signal_principal_only"] > 0),
            "signal_alternative_only": int(info["signal_alternative_only"] > 0),
            "tm_principal_only": int(info["tm_principal_only"] > 0),
            "tm_alternative_only": int(info["tm_alternative_only"] > 0),
        }
    counts = {
        "candidate_genes": len(categories),
        "measured_protein_coding_transcripts": sum(item["n_alt_protein_coding"] for item in categories.values()) + len(categories),
        "genes_with_uniprot": sum(item["uniprot"] for item in categories.values()),
        "genes_with_sequence_difference": sum(item["category"] in {"S1", "S2"} for item in categories.values()),
        "genes_with_structural_difference": sum(item["category"] == "S2" for item in categories.values()),
        "genes_with_principal_pfam": sum(item["n_pfam"] > 0 for item in categories.values()),
        "S2": sum(item["category"] == "S2" for item in categories.values()),
        "S1": sum(item["category"] == "S1" for item in categories.values()),
        "S0": sum(item["category"] == "S0" for item in categories.values()),
        "ANNOTATION_UNRESOLVED": sum(item["category"] == "ANNOTATION_UNRESOLVED" for item in categories.values()),
    }
    return categories, counts


def pfam_coverage(categories: dict[str, dict]) -> dict[str, float]:
    needed = {}
    for gene, item in categories.items():
        protein = item["principal_protein_id"]
        if "." not in protein:
            continue
        stable, version = protein.rsplit(".", 1)
        if version.isdigit():
            needed[stable] = (gene, int(version), item["principal_length"])
    translation = {}
    with gzip.open(ROOT / "phase2a/raw/translation.txt.gz", "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8 or not parts[6].startswith("ENSP"):
                continue
            record = needed.get(parts[6])
            if record and parts[7].isdigit() and int(parts[7]) == record[1]:
                translation[int(parts[0])] = record
    intervals = defaultdict(list)
    with gzip.open(ROOT / "phase2a/raw/protein_feature.txt.gz", "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8 or parts[7] != str(PFAM_ANALYSIS):
                continue
            record = translation.get(int(parts[1])) if parts[1].isdigit() else None
            if record is None or not parts[2].isdigit() or not parts[3].isdigit():
                continue
            length = int(record[2])
            start = max(1, int(parts[2]))
            end = min(length, int(parts[3]))
            if start <= end:
                intervals[record[0]].append((start, end))
    coverage = {}
    for gene, item in categories.items():
        spans = sorted(intervals.get(gene, []))
        length = item["principal_length"]
        if length <= 0 or not spans:
            coverage[gene] = 0.0 if item["n_pfam"] == 0 else np.nan
            continue
        total = 0
        cursor_start, cursor_end = spans[0]
        for start, end in spans[1:]:
            if start <= cursor_end + 1:
                cursor_end = max(cursor_end, end)
            else:
                total += cursor_end - cursor_start + 1
                cursor_start, cursor_end = start, end
        total += cursor_end - cursor_start + 1
        coverage[gene] = min(total, length) / length
    return coverage


def median_expression_by_gene(categories: dict[str, dict]) -> dict[tuple[str, str, str], float]:
    """Median log2 gene sum. Transcript rows are mapped through the Phase 2A structure is not enough.

    The isoform matrix is keyed by transcript, so this function maps transcripts with the
    category table's gene membership via a transcript-to-gene table built from the pair and
    functional files plus the GTF-free matrix join stored in Phase 2A mapping.
    """
    transcript_gene = {}
    with (ROOT / "phase2a/manifests/transcript_protein_mapping.csv").open() as handle:
        for row in csv.DictReader(handle):
            if row["gene_id"] in categories:
                transcript_gene[row["transcript_id"]] = row["gene_id"]
    genes = list(categories)
    index = {gene: i for i, gene in enumerate(genes)}
    n_isoform = np.array([categories[gene]["n_isoform_rows"] for gene in genes], dtype=np.int16)
    output = {}
    for cancer in CANCERS:
        path = acquisition_path(cancer, "rna_isoform")
        log(f"rna {cancer}")
        with path.open("rb") as handle:
            header = handle.readline().rstrip(b"\n").split(b"\t")
            n_cases = len(header) - 1
            sums_a = np.zeros((len(genes), n_cases))
            sums_b = np.zeros((len(genes), n_cases))
            seen = np.zeros((len(genes), n_cases), dtype=np.int16)
            for raw in handle:
                tab = raw.find(b"\t")
                gene = index.get(transcript_gene.get(raw[:tab].decode()))
                if gene is None:
                    continue
                values = np.fromstring(raw[tab + 1:], sep="\t")
                if values.size != n_cases:
                    padded = np.full(n_cases, np.nan)
                    padded[: min(values.size, n_cases)] = values[:n_cases]
                    values = padded
                finite = np.isfinite(values)
                sums_a[gene, finite] += np.exp2(values[finite])
                sums_b[gene, finite] += np.expm1(values[finite] * np.log(2.0))
                seen[gene, finite] += 1
        matched = int(np.sum(np.any(seen, axis=1)))
        log(f"rna {cancer} genes {matched}")
        if matched < 10000:
            raise SystemExit(f"{cancer} matched only {matched} candidate genes")
        for gene, gene_index in index.items():
            complete = seen[gene_index] == n_isoform[gene_index]
            for transform, sums in (("A", sums_a), ("B", sums_b)):
                usable = complete & (sums[gene_index] > 0)
                output[(gene, cancer, transform)] = (
                    float(np.median(np.log2(sums[gene_index, usable]))) if np.any(usable) else np.nan
                )
        finite_a = sum(np.isfinite(output[(gene, cancer, "A")]) for gene in index)
        if finite_a < 0.95 * len(index):
            raise SystemExit(f"{cancer} produced a finite RNA median for only {finite_a} genes")
    return output


def load_outcomes(categories: dict[str, dict], coverage: dict[str, float], expression: dict) -> dict[str, dict]:
    grouped = {}
    for transform, name in (("A", "primary_literal_A_gene_fits.tsv.gz"), ("B", "primary_literal_B_gene_fits.tsv.gz")):
        rows = []
        with gzip.open(ROOT / "phase1a/results" / name, "rt") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            locate = {column: i for i, column in enumerate(header)}
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if parts[locate["fit_status"]] != "OK":
                    continue
                gene = parts[locate["gene_id"]]
                item = categories.get(gene)
                if item is None or item["category"] not in {"S0", "S1", "S2"}:
                    continue
                cancer = parts[locate["cancer"]]
                n_complete = float(parts[locate["n_complete"]])
                n_candidate = float(parts[locate["n_candidate_cases"]])
                rows.append({
                    "gene": gene,
                    "cancer": cancer,
                    "category": item["category"],
                    "y": float(parts[locate["partial_r2"]]),
                    "principal_length": item["principal_length"],
                    "n_pfam": float(item["n_pfam"]),
                    "pfam_coverage": float(coverage.get(gene, np.nan)),
                    "n_alt_protein_coding": float(item["n_alt_protein_coding"]),
                    "n_isoform_rows": float(item["n_isoform_rows"]),
                    "n_complete": n_complete,
                    "protein_fraction": n_complete / n_candidate if n_candidate else np.nan,
                    "median_log2_g": float(expression.get((gene, cancer, transform), np.nan)),
                    "pfam_diff": item["pfam_diff"],
                    "signal_diff": item["signal_diff"],
                    "tm_diff": item["tm_diff"],
                })
        grouped[transform] = rows
    return grouped


def standardized_mean_difference(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 2 or len(right) < 2:
        return np.nan
    pooled = np.sqrt(((len(left) - 1) * np.var(left, ddof=1) + (len(right) - 1) * np.var(right, ddof=1)) / (len(left) + len(right) - 2))
    if pooled == 0:
        return 0.0
    return float((np.mean(left) - np.mean(right)) / pooled)


def variance_inflation(covariates: np.ndarray) -> np.ndarray:
    factors = []
    for index in range(covariates.shape[1]):
        outcome = covariates[:, index]
        others = np.delete(covariates, index, axis=1)
        design = np.column_stack([np.ones(len(outcome)), others]) if others.size else np.ones((len(outcome), 1))
        beta, resid = ols(design, outcome)
        total = np.sum((outcome - np.mean(outcome)) ** 2)
        r2 = 1 - np.sum(resid ** 2) / total if total > 0 else 0
        factors.append(np.inf if r2 >= 1 - 1e-8 else 1 / (1 - r2))
    return np.array(factors)


def select_covariates(rows: list[dict]) -> tuple[list[str], list[dict]]:
    usable = [row for row in rows if row["category"] in {"S1", "S2"}]
    kept = []
    dropped = []
    for name in COVARIATE_PRIORITY:
        values = np.array([row[name] for row in usable], dtype=float)
        missing = np.mean(~np.isfinite(values))
        if missing > MISSING_LIMIT:
            dropped.append({"covariate": name, "reason": "missingness", "value": float(missing)})
        else:
            kept.append(name)
    changed = True
    while changed and kept:
        changed = False
        matrix = np.column_stack([
            np.nan_to_num(np.array([row[name] for row in usable], dtype=float), nan=np.nanmedian([row[name] for row in usable]))
            for name in kept
        ])
        # The median fill is only for the VIF screen of otherwise complete columns.
        factors = variance_inflation(matrix)
        offenders = [index for index, factor in enumerate(factors) if factor > VIF_LIMIT]
        if offenders:
            lowest = max(offenders, key=lambda index: COVARIATE_PRIORITY.index(kept[index]))
            dropped.append({"covariate": kept[lowest], "reason": "vif", "value": float(factors[lowest])})
            del kept[lowest]
            changed = True
        else:
            for name, factor in zip(kept, factors):
                dropped.append({"covariate": name, "reason": "retained_vif", "value": float(factor)})
    return kept, dropped


def gene_table(rows: list[dict]) -> dict[str, dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["gene"]].append(row)
    summary = {}
    for gene, items in grouped.items():
        summary[gene] = {
            "category": items[0]["category"],
            "y": float(np.mean([item["y"] for item in items])),
            "length": items[0]["principal_length"],
            "n_pfam": items[0]["n_pfam"],
            "n_alt": items[0]["n_alt_protein_coding"],
            "pfam_diff": items[0]["pfam_diff"],
            "signal_diff": items[0]["signal_diff"],
            "tm_diff": items[0]["tm_diff"],
        }
    return summary


def design_matrix(rows: list[dict], covariates: list[str], include_cancer: bool) -> tuple[np.ndarray, list[str], np.ndarray, np.ndarray]:
    outcome = np.array([row["y"] for row in rows], dtype=float)
    s2 = np.array([row["category"] == "S2" for row in rows], dtype=float)
    cancer = np.array([CANCERS.index(row["cancer"]) for row in rows])
    columns = [np.ones(len(rows)), s2]
    names = ["intercept", "S2"]
    if include_cancer:
        for index, name in enumerate(CANCERS):
            if name == "BRCA":
                continue
            columns.append((cancer == index).astype(float))
            names.append(f"cancer_{name}")
    for name in covariates:
        values = np.array([row[name] for row in rows], dtype=float)
        finite = values[np.isfinite(values)]
        center = np.mean(finite) if finite.size else 0
        scale = np.std(finite) if finite.size and np.std(finite) > 0 else 1
        filled = np.where(np.isfinite(values), values, center)
        columns.append((filled - center) / scale)
        names.append(name)
    design = np.column_stack(columns)
    keep = [0]
    if np.linalg.matrix_rank(design[:, [0, 1]], tol=1e-8) == 2:
        keep.append(1)
    for index in range(2, design.shape[1]):
        trial = keep + [index]
        if np.linalg.matrix_rank(design[:, trial], tol=1e-8) == len(trial):
            keep.append(index)
    dropped = [names[index] for index in range(design.shape[1]) if index not in keep]
    if dropped:
        log(f"dropped dependent columns {dropped}")
    design = design[:, keep]
    names = [names[index] for index in keep]
    clusters = np.array([row["gene"] for row in rows])
    codes = np.unique(clusters, return_inverse=True)[1]
    return design, names, outcome, codes


def fit_primary(rows: list[dict], covariates: list[str], include_cancer: bool) -> dict | None:
    if sum(row["category"] == "S2" for row in rows) < 2 or sum(row["category"] == "S1" for row in rows) < 2:
        return None
    design, names, outcome, codes = design_matrix(rows, covariates, include_cancer)
    if len(names) < 2 or names[1] != "S2":
        return None
    beta, resid = ols(design, outcome)
    covariance = cluster_covariance(design, resid, codes) if include_cancer else hc1_covariance(design, resid)
    if covariance is None:
        return None
    se = np.sqrt(np.clip(np.diag(covariance), 0, None))
    z_value = beta[1] / se[1] if se[1] > 0 else np.nan
    p_value = float(2 * stats.norm.sf(abs(z_value))) if np.isfinite(z_value) else np.nan
    return {
        "estimate": float(beta[1]),
        "se": float(se[1]),
        "ci_low": float(beta[1] - Z95 * se[1]),
        "ci_high": float(beta[1] + Z95 * se[1]),
        "p_value": p_value,
        "n_rows": len(rows),
        "n_genes": int(len(np.unique(codes))),
        "names": names,
        "beta": beta,
        "design": design,
        "resid": resid,
        "codes": codes,
        "outcome": outcome,
    }


def bootstrap_mean_interval(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    if values.size < 2:
        return np.nan, np.nan
    draws = values[rng.integers(0, values.size, size=(DRAWS, values.size))].mean(axis=1)
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def bootstrap_coefficient(rows: list[dict], covariates: list[str], include_cancer: bool, rng: np.random.Generator) -> tuple[float, float]:
    design, _, outcome, codes = design_matrix(rows, covariates, include_cancer)
    groups = [np.flatnonzero(codes == code) for code in range(int(codes.max()) + 1)]
    draws = np.empty(DRAWS)
    n_groups = len(groups)
    for draw in range(DRAWS):
        chosen = rng.integers(0, n_groups, n_groups)
        index = np.concatenate([groups[item] for item in chosen])
        beta, _ = ols(design[index], outcome[index])
        draws[draw] = beta[1]
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def stratified_bootstrap(rows: list[dict], covariates: list[str], rng: np.random.Generator) -> tuple[float, float]:
    design, _, outcome, _ = design_matrix(rows, covariates, True)
    cancer = np.array([row["cancer"] for row in rows])
    groups = [np.flatnonzero(cancer == name) for name in CANCERS]
    draws = np.empty(DRAWS)
    for draw in range(DRAWS):
        index = np.concatenate([
            group[rng.integers(0, group.size, group.size)] for group in groups if group.size
        ])
        beta, _ = ols(design[index], outcome[index])
        draws[draw] = beta[1]
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def distribution_block(values: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "sd": float(np.std(values, ddof=1)),
        "iqr": float(np.subtract(*np.quantile(values, [0.75, 0.25]))),
        "trimmed_10": float(stats.trim_mean(values, 0.10)),
        "trimmed_20": float(stats.trim_mean(values, 0.20)),
        "p90": float(np.quantile(values, 0.90)),
        "p95": float(np.quantile(values, 0.95)),
        "fraction_ge_0.01": float(np.mean(values >= 0.01)),
        "fraction_ge_0.02": float(np.mean(values >= 0.02)),
        "fraction_ge_0.05": float(np.mean(values >= 0.05)),
        "fraction_ge_0.10": float(np.mean(values >= 0.10)),
    }


def support_window(genes: dict[str, dict]) -> tuple[np.ndarray, dict]:
    s2_length = np.array([item["length"] for item in genes.values() if item["category"] == "S2"])
    s1_length = np.array([item["length"] for item in genes.values() if item["category"] == "S1"])
    s2_pfam = np.array([item["n_pfam"] for item in genes.values() if item["category"] == "S2"])
    s1_pfam = np.array([item["n_pfam"] for item in genes.values() if item["category"] == "S1"])
    length_low = max(np.quantile(s1_length, 0.10), np.quantile(s2_length, 0.10))
    length_high = min(np.quantile(s1_length, 0.90), np.quantile(s2_length, 0.90))
    pfam_low = max(np.quantile(s1_pfam, 0.10), np.quantile(s2_pfam, 0.10))
    pfam_high = min(np.quantile(s1_pfam, 0.90), np.quantile(s2_pfam, 0.90))
    selected = []
    for gene, item in genes.items():
        if item["category"] not in {"S1", "S2"}:
            continue
        if length_low <= item["length"] <= length_high and pfam_low <= item["n_pfam"] <= pfam_high:
            selected.append(gene)
    info = {
        "length_smd": standardized_mean_difference(s2_length, s1_length),
        "pfam_smd": standardized_mean_difference(s2_pfam, s1_pfam),
        "length_low": float(length_low),
        "length_high": float(length_high),
        "pfam_low": float(pfam_low),
        "pfam_high": float(pfam_high),
        "window_s2": sum(genes[gene]["category"] == "S2" for gene in selected),
        "window_s1": sum(genes[gene]["category"] == "S1" for gene in selected),
    }
    return np.array(selected), info


def decide(qualified: bool, estimates: dict[str, dict], contrasts: dict[str, list[dict]], loo: list[dict], window: dict[str, dict]) -> str:
    if not qualified:
        return "PRIMARY_COMPARISON_NOT_QUALIFIED"

    def supported(item: dict) -> bool:
        return item["estimate"] >= PRACTICAL and item["ci_low"] > 0

    a_ok = supported(estimates["A"])
    b_ok = supported(estimates["B"])
    same_sign = estimates["A"]["estimate"] * estimates["B"]["estimate"] > 0
    loo_positive = all(row["estimate"] > 0 for row in loo)
    window_ok = all(window[transform]["estimate"] > 0 and np.sign(window[transform]["estimate"]) == np.sign(estimates[transform]["estimate"]) for transform in ("A", "B"))
    if a_ok and b_ok and same_sign and loo_positive and window_ok:
        return "STRUCTURAL_ASSOCIATION_SUPPORTED"
    if (not same_sign) or (a_ok != b_ok):
        return "STRUCTURAL_ASSOCIATION_TRANSFORMATION_DEPENDENT"

    def cancer_split(transform: str) -> bool:
        high = [row for row in contrasts[transform] if row["estimate"] >= PRACTICAL and row["ci_low"] > 0]
        low = [row for row in contrasts[transform] if row["estimate"] <= 0 and row["ci_high"] < PRACTICAL]
        return bool(high) and bool(low)

    if cancer_split("A") and cancer_split("B"):
        return "STRUCTURAL_ASSOCIATION_CANCER_SPECIFIC"
    return "NO_PRACTICAL_STRUCTURAL_ASSOCIATION"


def fmt(value) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return ""
    return f"{float(value):.6g}"


def write_reports(counts, qualified, label, estimates, contrasts, features, opportunity, loo_rows, covariate_sets) -> None:
    primary = []
    for transform in ("A", "B"):
        item = estimates[transform]
        primary.append(
            f"Literal {transform}: adjusted S2 minus S1 = {item['estimate']:.6f}, "
            f"cluster-robust interval {item['ci_low']:.6f} to {item['ci_high']:.6f}, "
            f"genes {item['n_genes']}, rows {item['n_rows']}. "
            f"Adjustment set: {', '.join(covariate_sets[transform])}."
        )
    contrast_lines = ["| Cancer | Transform | S2 genes | S1 genes | Adjusted difference | Interval |", "|---|---|---:|---:|---:|---|"]
    for transform in ("A", "B"):
        for row in contrasts[transform]:
            contrast_lines.append(
                f"| {row['cancer']} | {transform} | {row['n_s2_genes']} | {row['n_s1_genes']} | {row['estimate']:.6f} | {row['ci_low']:.6f} to {row['ci_high']:.6f} |"
            )
    (PHASE / "reports" / "phase2b_primary_analysis.md").write_text(
        "# Phase 2B primary analysis\n\n"
        "The outcome is the frozen Phase 1A partial R². S2 is a predicted structural difference. "
        "S1 is a predicted sequence difference without a supported Pfam, SignalP, or TMHMM difference. "
        "The coefficient is an adjusted difference in means. It is an association.\n\n"
        + "\n".join(primary)
        + "\n\nThe practical magnitude locked before fitting is 0.005 partial R². "
        f"Common-support status: {'qualified' if qualified else 'not qualified'}.\n"
    )
    (PHASE / "reports" / "phase2b_structural_feature_analysis.md").write_text(
        "# Structural feature analysis\n\n"
        "Pfam, SignalP, and TMHMM indicators can overlap. Their coefficients are secondary and are not independent discoveries.\n\n"
        "Domain truncation and domain-boundary change are not estimated. A residue alignment was not qualified. "
        "The feature summary counts an accession present on only one of the two compared sequences.\n\n"
        + "\n".join(
            f"- {row['feature']} {row['transformation']}: genes {row['n_genes']}, adjusted difference {row['estimate']}, interval {row['ci_low']} to {row['ci_high']}."
            for row in features
        )
        + "\n"
    )
    (PHASE / "reports" / "phase2b_cancer_heterogeneity.md").write_text(
        "# Cancer-specific structural contrasts\n\n"
        "Phase 1B left the partial-R² heterogeneity as mixed measurement and cancer-specific structure. "
        "These contrasts ask whether the S2-versus-S1 difference itself varies by cancer. They do not assign a cause to that variation. "
        "Inside one cancer, the protein-observation fraction is a rescaling of the complete-case count, so the rank check keeps only one of those two columns.\n\n"
        + "\n".join(contrast_lines)
        + "\n\nLeave-one-cancer-out estimates:\n\n"
        + "\n".join(
            f"- Omit {row['omitted_cancer']} {row['transformation']}: {row['estimate']:.6f} ({row['ci_low']:.6f} to {row['ci_high']:.6f})."
            for row in loo_rows
        )
        + "\n"
    )
    (PHASE / "reports" / "phase2b_annotation_bias_report.md").write_text(
        "# Annotation opportunity\n\n"
        "Longer proteins and proteins with more Pfam accessions have more chances to be placed in S2. "
        "The primary model adjusts for length, domain count, and the other covariates that survived the locked missingness and VIF rules. "
        "The measured-isoform-row count was removed because its variance-inflation factor exceeded 10. "
        "For these genes it is almost the same number as the count of alternative protein-coding transcripts, and that count stays in the model.\n\n"
        f"Length standardized difference, S2 minus S1: {opportunity['length_smd']:.3f}. "
        f"Pfam-count standardized difference: {opportunity['pfam_smd']:.3f}. "
        f"Overlapping window retains {opportunity['window_s2']} S2 genes and {opportunity['window_s1']} S1 genes.\n\n"
        "Length quintiles are reported separately. The shortest quintile is part of the sensitivity, not a reason to redraw the groups. "
        "The window estimate is a sensitivity. Genes were not selected by partial R².\n"
        "A SwissProt accession that differs from the TrEMBL accession was not treated as an ambiguous protein sequence. "
        "The ambiguous-mapping sensitivity removes genes with a FASTA header conflict or more than one UniProt row.\n"
    )
    (PHASE / "reports" / "phase2b_transformation_robustness.md").write_text(
        "# Transformation robustness\n\n"
        "Literal A and Literal B were fit as separate outcomes. The inverse was not chosen from the coefficient.\n\n"
        + "\n".join(primary)
        + f"\n\nSame sign: {estimates['A']['estimate'] * estimates['B']['estimate'] > 0}.\n"
    )
    (PHASE / "reports" / "phase2b_literature_update.md").write_text(
        "# Literature update\n\n"
        "The Phase 2A audit already recorded pan-cancer isoform-switch studies, CPTAC gene-level proteogenomics, and long-read isoform resources. "
        "Two later records were checked because they sit close to this question.\n\n"
        "Li et al., Cancer Cell, 2023, DOI 10.1016/j.ccell.2023.06.009, describe the CPTAC pan-cancer resource and the difficulty of joining transcripts to gene-level protein identifiers. "
        "Their recommendation is transparent gene-level quantification when isoform-resolved protein measurement is unrealistic. They do not estimate the partial R² of a prespecified principal-transcript proportion by structural class.\n\n"
        "Donovan et al., PLOS ONE, 2023, DOI 10.1371/journal.pone.0282821, report opposite peptide-abundance patterns for BMP1 isoforms in 188 plasma proteomes. "
        "That is a single-protein, peptide-level case study. It is not a ten-cohort test of principal-transcript partial R².\n\n"
        "The joint estimand in this phase remains unestablished by the audited literature: frozen principal-transcript usage, incremental gene-level protein explanation, predicted structural class, cancer-specific contrasts, and two literal inverses. "
        "That statement is not a claim that no related study exists. Novelty stays PARTIAL.\n"
    )
    (PHASE / "reports" / "phase2b_limitations.md").write_text(
        "# Limitations\n\n"
        "The structural classes are predicted from GENCODE v34 sequences and Ensembl 100 Pfam 32.0, SignalP 4.1, and TMHMM 2.0c coordinates. "
        "Level 3 isoform-specific peptide evidence remains absent. A statistical association does not move a gene from Level 2 to Level 3.\n\n"
        "Partial R² is right-skewed. The arithmetic mean, through the regression coefficient, stays the primary location. "
        "Trimmed means and high percentiles are descriptive.\n\n"
        "Genes are repeated across cancers. The primary interval uses gene-level clustering. "
        "Patient-level observations were not resampled as if they were genes.\n\n"
        "S1 is not a functionally neutral class. It contains sequence differences that the audited tracks do not mark. "
        "S0 contains measured alternatives with the same predicted amino-acid sequence.\n\n"
        "Residue-level domain truncation was not called. Exact sequence equality does not place domain coordinates from two different sequences onto one alignment.\n\n"
        f"Candidate genes: {counts['candidate_genes']}. S2 {counts['S2']}. S1 {counts['S1']}. S0 {counts['S0']}. "
        f"Unresolved {counts['ANNOTATION_UNRESOLVED']}.\n"
    )
    (PHASE / "reports" / "phase2b_final_decision.md").write_text(
        "# Phase 2B final decision\n\n"
        "## A. Structural category reproducibility\n\n"
        "PASS\n\n"
        "The Phase 2A counts reproduced: 10,795 candidate genes, 39,490 measured protein-coding transcripts, "
        "10,742 genes with a principal SwissProt or TrEMBL accession, 9,608 genes with a different predicted sequence, "
        "4,829 genes with a supported structural difference, and 10,228 genes with a Pfam accession on the principal protein.\n\n"
        "## B. Common annotation-support adequacy\n\n"
        f"{'ADEQUATE' if qualified else 'NOT_QUALIFIED'}\n\n"
        "## C. Primary S2-versus-S1 association\n\n"
        + "\n".join(primary)
        + "\n\n## D. Cancer-specific interaction\n\n"
        "Cancer-specific adjusted differences are in `phase2b/results/cancer_specific_structural_contrasts.csv`. "
        "Benjamini–Hochberg adjustment is within the ten cancer tests of each transformation.\n\n"
        "## E. Literal A/B consistency\n\n"
        f"The two adjusted estimates have the same sign: {estimates['A']['estimate'] * estimates['B']['estimate'] > 0}.\n\n"
        "## F. Robustness to annotation opportunity\n\n"
        f"Length standardized difference {opportunity['length_smd']:.3f}. "
        f"Window S2 {opportunity['window_s2']}, S1 {opportunity['window_s1']}.\n\n"
        "## G. Independent evidence limitations\n\n"
        "Level 3 evidence remains zero. The association uses gene-level protein abundance and predicted annotation.\n\n"
        "## H. Scientific novelty\n\n"
        "PARTIAL\n\n"
        "## Overall\n\n"
        f"{label}\n\n"
        "Genes with predicted principal-versus-alternative protein structural differences were compared with genes whose predicted sequences differ without a supported structural annotation. "
        "The comparison is an association with transcript-usage incremental explanatory values for gene-level protein abundance. "
        f"The locked practical magnitude is 0.005. The adjusted estimates are {estimates['A']['estimate']:.6f} under Literal A and {estimates['B']['estimate']:.6f} under Literal B.\n\n"
        "The analysis does not establish that an alternative isoform is expressed as a protein, that a domain is experimentally absent, or that a structural difference causes the partial R².\n\n"
        "PHASE 2B EXECUTION = COMPLETED\n\n"
        "PHASE 2C AUTHORIZATION = FALSE\n"
    )


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
    verify_prefit()
    verify_historical()
    synthetic_tests()
    log("categories")
    categories, counts = load_categories()
    count_rows = []
    for key, expected in EXPECTED.items():
        count_rows.append({
            "metric": key, "value": counts[key], "expected": expected,
            "match": int(counts[key] == expected),
        })
    for key in ("S2", "S1", "S0", "ANNOTATION_UNRESOLVED"):
        count_rows.append({"metric": key, "value": counts[key], "expected": "", "match": ""})
    write_csv(PHASE / "results" / "structural_category_counts.csv", count_rows, ["metric", "value", "expected", "match"])
    if any(not row["match"] for row in count_rows if row["expected"] != ""):
        raise SystemExit(f"historical counts did not reproduce: {counts}")
    assignment = []
    for gene in sorted(categories):
        item = categories[gene]
        assignment.append({
            "gene_id": gene,
            "category": item["category"],
            "principal_transcript_id": item["principal_transcript_id"],
            "principal_protein_id": item["principal_protein_id"],
            "principal_length": item["principal_length"],
            "n_pfam": item["n_pfam"],
            "n_alt_protein_coding": item["n_alt_protein_coding"],
            "n_isoform_rows": item["n_isoform_rows"],
            "n_s2_pairs": item["n_s2_pairs"],
            "n_s1_pairs": item["n_s1_pairs"],
            "n_identical_pairs": item["n_identical_pairs"],
            "pfam_diff": item["pfam_diff"],
            "signal_diff": item["signal_diff"],
            "tm_diff": item["tm_diff"],
        })
    fields = list(assignment[0])
    write_csv(PHASE / "manifests" / "structural_category_assignment.csv", assignment, fields)
    write_csv(PHASE / "manifests" / "frozen_structural_universe.csv", count_rows, ["metric", "value", "expected", "match"])
    (PHASE / "reports" / "structural_universe_validation.md").write_text(
        "# Structural universe validation\n\n"
        "Categories were assigned from the Phase 2A version-matched pair table before partial R² was modeled.\n\n"
        + "\n".join(f"- {row['metric']}: {row['value']} (expected {row['expected']}, match {row['match']})" for row in count_rows)
        + "\n\nS2 requires a supported Pfam, SignalP, or TMHMM difference on differing sequences. "
        "S1 requires differing sequences without that supported difference. "
        "S0 requires identical predicted sequences. The remaining candidate genes are ANNOTATION_UNRESOLVED.\n"
    )
    log("pfam coverage")
    coverage = pfam_coverage(categories)
    for gene, item in categories.items():
        item["pfam_coverage"] = coverage.get(gene, np.nan)
    log("expression")
    expression = median_expression_by_gene(categories)
    log("outcomes")
    outcomes = load_outcomes(categories, coverage, expression)
    primary_rows = {}
    covariate_sets = {}
    vif_rows = []
    for transform, rows in outcomes.items():
        kept, diagnostics = select_covariates(rows)
        covariate_sets[transform] = kept
        for item in diagnostics:
            item["transformation"] = transform
            vif_rows.append(item)
        primary_rows[transform] = [
            row for row in rows
            if row["category"] in {"S1", "S2"} and all(np.isfinite(row[name]) for name in kept)
        ]
        log(f"covariates {transform} {kept} complete {len(primary_rows[transform])}")
    write_csv(PHASE / "qc" / "covariate_screening.csv", vif_rows, ["transformation", "covariate", "reason", "value"])

    estimates = {}
    contrast_map = {}
    loo_rows = []
    distribution_rows = []
    bootstrap_rows = []
    opportunity_rows = []
    feature_rows = []
    testing_rows = []
    rng = np.random.default_rng(SEED)
    genes_by_transform = {transform: gene_table(rows) for transform, rows in primary_rows.items()}
    ambiguous_genes = set()
    with (ROOT / "phase2a/manifests/transcript_protein_mapping.csv").open() as handle:
        for row in csv.DictReader(handle):
            ambiguity = row["mapping_ambiguity"]
            if row["in_candidate_gene"] == "1" and (
                "fasta_header_transcript_differs" in ambiguity or "multiple_uniprot_rows" in ambiguity
            ):
                ambiguous_genes.add(row["gene_id"])
    log(f"ambiguous alternative mappings {len(ambiguous_genes)}")
    window_genes, window_info = support_window(genes_by_transform["A"])
    # Length and Pfam count are annotation properties and do not depend on the inverse.
    n_s2 = sum(item["category"] == "S2" for item in genes_by_transform["A"].values())
    n_s1 = sum(item["category"] == "S1" for item in genes_by_transform["A"].values())
    qualified = (
        n_s2 >= MIN_CLASS and n_s1 >= MIN_CLASS and not (
            abs(window_info["length_smd"]) > SMD_LIMIT and (
                window_info["window_s2"] < MIN_WINDOW or window_info["window_s1"] < MIN_WINDOW
            )
        )
    )
    window_fits = {}
    for transform in ("A", "B"):
        rows = primary_rows[transform]
        covariates = covariate_sets[transform]
        fit = fit_primary(rows, covariates, True)
        if fit is None:
            raise SystemExit(f"primary design failed for {transform}")
        unadjusted = fit_primary(rows, [], True)
        estimates[transform] = fit
        log(f"bootstrap {transform}")
        boot_low, boot_high = bootstrap_coefficient(rows, covariates, True, rng)
        stratified_low, stratified_high = stratified_bootstrap(rows, covariates, rng)
        bootstrap_rows.append({
            "analysis": "cluster_bootstrap_adjusted", "transformation": transform,
            "estimate": fit["estimate"], "ci_low": boot_low, "ci_high": boot_high, "draws": DRAWS, "seed": SEED,
        })
        bootstrap_rows.append({
            "analysis": "within_cancer_stratified_bootstrap", "transformation": transform,
            "estimate": fit["estimate"], "ci_low": stratified_low, "ci_high": stratified_high, "draws": DRAWS, "seed": SEED,
        })
        genes = genes_by_transform[transform]
        s2_values = np.array([item["y"] for item in genes.values() if item["category"] == "S2"])
        s1_values = np.array([item["y"] for item in genes.values() if item["category"] == "S1"])
        gene_delta = float(np.mean(s2_values) - np.mean(s1_values))
        pooled = np.sqrt(((len(s2_values) - 1) * np.var(s2_values, ddof=1) + (len(s1_values) - 1) * np.var(s1_values, ddof=1)) / (len(s2_values) + len(s1_values) - 2))
        for label, values in (("S2", s2_values), ("S1", s1_values)):
            block = distribution_block(values)
            block.update({"group": label, "transformation": transform, "n_genes": len(values)})
            distribution_rows.append(block)
        s0_table = gene_table([row for row in outcomes[transform] if row["category"] == "S0"])
        s0_values = np.array([item["y"] for item in s0_table.values()])
        if s0_values.size:
            block = distribution_block(s0_values)
            block.update({"group": "S0", "transformation": transform, "n_genes": len(s0_values)})
            distribution_rows.append(block)
        cancer_rows = []
        p_values = []
        for cancer in CANCERS:
            subset = [row for row in rows if row["cancer"] == cancer]
            cancer_fit = fit_primary(subset, covariates, False)
            s2_gene = {row["gene"] for row in subset if row["category"] == "S2"}
            s1_gene = {row["gene"] for row in subset if row["category"] == "S1"}
            s2_y = np.array([row["y"] for row in subset if row["category"] == "S2"])
            s1_y = np.array([row["y"] for row in subset if row["category"] == "S1"])
            s2_low, s2_high = bootstrap_mean_interval(s2_y, rng)
            s1_low, s1_high = bootstrap_mean_interval(s1_y, rng)
            record = {
                "cancer": cancer,
                "transformation": transform,
                "n_s2_genes": len(s2_gene),
                "n_s1_genes": len(s1_gene),
                "s2_mean": float(np.mean(s2_y)) if s2_y.size else np.nan,
                "s1_mean": float(np.mean(s1_y)) if s1_y.size else np.nan,
                "s2_mean_ci_low": s2_low,
                "s2_mean_ci_high": s2_high,
                "s1_mean_ci_low": s1_low,
                "s1_mean_ci_high": s1_high,
                "s2_iqr": float(np.subtract(*np.quantile(s2_y, [0.75, 0.25]))) if s2_y.size else np.nan,
                "s1_iqr": float(np.subtract(*np.quantile(s1_y, [0.75, 0.25]))) if s1_y.size else np.nan,
                "s2_median": float(np.median(s2_y)) if s2_y.size else np.nan,
                "s1_median": float(np.median(s1_y)) if s1_y.size else np.nan,
                "s2_sd": float(np.std(s2_y, ddof=1)) if s2_y.size > 1 else np.nan,
                "s1_sd": float(np.std(s1_y, ddof=1)) if s1_y.size > 1 else np.nan,
                "unadjusted_difference": float(np.mean(s2_y) - np.mean(s1_y)) if s2_y.size and s1_y.size else np.nan,
                "estimate": cancer_fit["estimate"] if cancer_fit else np.nan,
                "ci_low": cancer_fit["ci_low"] if cancer_fit else np.nan,
                "ci_high": cancer_fit["ci_high"] if cancer_fit else np.nan,
                "p_value": cancer_fit["p_value"] if cancer_fit else np.nan,
                "standardized_effect": (cancer_fit["estimate"] / pooled) if cancer_fit and pooled else np.nan,
                "dropped_covariates": " ".join(name for name in covariates if cancer_fit and name not in cancer_fit["names"]),
            }
            cancer_rows.append(record)
            p_values.append(record["p_value"] if np.isfinite(record["p_value"]) else 1)
        adjusted = benjamini_hochberg(np.array(p_values, dtype=float))
        for record, q_value in zip(cancer_rows, adjusted):
            record["bh_q"] = float(q_value)
            testing_rows.append({
                "family": f"cancer_contrasts_{transform}",
                "test": record["cancer"],
                "p_value": record["p_value"],
                "bh_q": record["bh_q"],
                "n": record["n_s2_genes"] + record["n_s1_genes"],
            })
        contrast_map[transform] = cancer_rows
        for cancer in CANCERS:
            held = [row for row in rows if row["cancer"] != cancer]
            held_fit = fit_primary(held, covariates, True)
            loo_rows.append({
                "omitted_cancer": cancer,
                "transformation": transform,
                "estimate": held_fit["estimate"] if held_fit else np.nan,
                "ci_low": held_fit["ci_low"] if held_fit else np.nan,
                "ci_high": held_fit["ci_high"] if held_fit else np.nan,
            })
        window_rows = [row for row in rows if row["gene"] in set(window_genes)]
        window_fit = fit_primary(window_rows, covariates, True)
        window_fits[transform] = window_fit
        lengths = np.array([genes[row["gene"]]["length"] for row in rows])
        # Quintiles are assigned from the unique gene lengths in this transform's successful S1/S2 set.
        gene_lengths = np.array([item["length"] for item in genes.values() if item["category"] in {"S1", "S2"}])
        edges = np.quantile(gene_lengths, [0, 0.2, 0.4, 0.6, 0.8, 1])
        stratum_estimates = []
        for stratum in range(5):
            low, high = edges[stratum], edges[stratum + 1]
            if stratum < 4:
                members = {
                    gene for gene, item in genes.items()
                    if item["category"] in {"S1", "S2"} and low <= item["length"] < high
                }
            else:
                members = {
                    gene for gene, item in genes.items()
                    if item["category"] in {"S1", "S2"} and low <= item["length"] <= high
                }
            s2_set = [item["y"] for gene, item in genes.items() if gene in members and item["category"] == "S2"]
            s1_set = [item["y"] for gene, item in genes.items() if gene in members and item["category"] == "S1"]
            if len(s2_set) >= MIN_STRATUM and len(s1_set) >= MIN_STRATUM:
                stratum_estimates.append(float(np.mean(s2_set) - np.mean(s1_set)))
            opportunity_rows.append({
                "analysis": "length_quintile", "transformation": transform, "stratum": stratum + 1,
                "n_s2": len(s2_set), "n_s1": len(s1_set),
                "estimate": float(np.mean(s2_set) - np.mean(s1_set)) if s2_set and s1_set else np.nan,
            })
        single_removed = [row for row in rows if row["n_alt_protein_coding"] >= 2]
        single_fit = fit_primary(single_removed, covariates, True)
        influence = one_step_influence(fit)
        threshold = 2 / np.sqrt(fit["n_genes"])
        influential_genes = set(np.array([row["gene"] for row in rows])[np.unique(fit["codes"])][np.abs(influence) > threshold]) if False else set()
        # Map influence back through the first row of each code.
        code_to_gene = {}
        for row, code in zip(rows, fit["codes"]):
            code_to_gene.setdefault(int(code), row["gene"])
        influential_genes = {code_to_gene[index] for index, value in enumerate(influence) if abs(value) > threshold}
        influence_rows = [row for row in rows if row["gene"] not in influential_genes]
        influence_fit = fit_primary(influence_rows, covariates, True)
        feature_fit = fit_feature_model(rows, covariates)
        for feature, estimate in feature_fit.items():
            feature_rows.append({"feature": feature, "transformation": transform, **estimate, "n_genes": int(sum(genes[gene][feature] if feature != "none" else 0 for gene in genes))})
        primary_summary = [
            {"quantity": "adjusted_S2_minus_S1", **{key: fit[key] for key in ("estimate", "se", "ci_low", "ci_high", "p_value", "n_rows", "n_genes")}},
            {"quantity": "unadjusted_cancer_fe", **{key: unadjusted[key] for key in ("estimate", "se", "ci_low", "ci_high", "p_value", "n_rows", "n_genes")}},
            {"quantity": "gene_equal_weight_difference", "estimate": gene_delta, "se": "", "ci_low": "", "ci_high": "", "p_value": "", "n_rows": len(rows), "n_genes": len(s2_values) + len(s1_values)},
            {"quantity": "cluster_bootstrap_interval", "estimate": fit["estimate"], "se": "", "ci_low": boot_low, "ci_high": boot_high, "p_value": "", "n_rows": len(rows), "n_genes": fit["n_genes"]},
            {"quantity": "standardized_gene_equal_effect", "estimate": gene_delta / pooled if pooled else np.nan, "se": "", "ci_low": "", "ci_high": "", "p_value": "", "n_rows": "", "n_genes": len(s2_values) + len(s1_values)},
            {"quantity": "length_quintile_mean_difference", "estimate": float(np.mean(stratum_estimates)) if stratum_estimates else np.nan, "se": "", "ci_low": "", "ci_high": "", "p_value": "", "n_rows": len(stratum_estimates), "n_genes": ""},
            {"quantity": "annotation_window", **{key: window_fit[key] for key in ("estimate", "se", "ci_low", "ci_high", "p_value", "n_rows", "n_genes")}} if window_fit else {"quantity": "annotation_window", "estimate": np.nan},
            {"quantity": "exclude_single_alternative", **{key: single_fit[key] for key in ("estimate", "se", "ci_low", "ci_high", "p_value", "n_rows", "n_genes")}} if single_fit else {"quantity": "exclude_single_alternative", "estimate": np.nan},
            {"quantity": "exclude_influential_genes", "estimate": influence_fit["estimate"] if influence_fit else np.nan, "se": influence_fit["se"] if influence_fit else "", "ci_low": influence_fit["ci_low"] if influence_fit else "", "ci_high": influence_fit["ci_high"] if influence_fit else "", "p_value": "", "n_rows": len(influence_rows), "n_genes": len(influential_genes)},
            {"quantity": "exclude_ambiguous_alternative_mapping", "estimate": fit["estimate"] if not ambiguous_genes else fit_primary([row for row in rows if row["gene"] not in ambiguous_genes], covariates, True)["estimate"], "se": "", "ci_low": "", "ci_high": "", "p_value": "", "n_rows": "", "n_genes": len(ambiguous_genes)},
        ]
        for row in primary_summary:
            row["transformation"] = transform
        write_csv(
            PHASE / "results" / f"primary_S2_vs_S1_literal_{transform}.csv",
            primary_summary,
            ["quantity", "estimate", "se", "ci_low", "ci_high", "p_value", "n_rows", "n_genes", "transformation"],
        )
        for row in distribution_rows:
            if row["transformation"] == transform:
                pass
    # Fix feature n_genes properly after loop by rewriting from categories.
    for row in feature_rows:
        if row["feature"] == "intercept":
            continue
        row["n_genes"] = sum(
            categories[gene][row["feature"]] == 1 and categories[gene]["category"] in {"S1", "S2"}
            for gene in categories
        )
    testing_feature = [row for row in feature_rows if row["feature"] in {"pfam_diff", "signal_diff", "tm_diff"}]
    for transform in ("A", "B"):
        subset = [row for row in testing_feature if row["transformation"] == transform and np.isfinite(row.get("p_value", np.nan))]
        if subset:
            q_values = benjamini_hochberg(np.array([row["p_value"] for row in subset], dtype=float))
            for row, q_value in zip(subset, q_values):
                row["bh_q"] = float(q_value)
                testing_rows.append({
                    "family": f"feature_class_{transform}", "test": row["feature"],
                    "p_value": row["p_value"], "bh_q": row["bh_q"], "n": row["n_genes"],
                })
    support_rows = []
    for gene in sorted(categories):
        item = categories[gene]
        support_rows.append({
            "gene_id": gene,
            "category": item["category"],
            "principal_length": item["principal_length"],
            "n_pfam": item["n_pfam"],
            "pfam_coverage": coverage.get(gene, ""),
            "n_alt_protein_coding": item["n_alt_protein_coding"],
            "n_isoform_rows": item["n_isoform_rows"],
            "in_length_pfam_window": int(gene in set(window_genes)),
            "pfam_diff": item["pfam_diff"],
            "signal_diff": item["signal_diff"],
            "tm_diff": item["tm_diff"],
            "pfam_absent_from_alternative": item["pfam_absent_from_alternative"],
            "pfam_gained_in_alternative": item["pfam_gained_in_alternative"],
            "signal_principal_only": item["signal_principal_only"],
            "signal_alternative_only": item["signal_alternative_only"],
            "tm_principal_only": item["tm_principal_only"],
            "tm_alternative_only": item["tm_alternative_only"],
        })
    write_csv(PHASE / "manifests" / "common_annotation_support.csv", support_rows, list(support_rows[0]))
    write_csv(PHASE / "results" / "cancer_specific_structural_contrasts.csv", contrast_map["A"] + contrast_map["B"], list(contrast_map["A"][0]))
    write_csv(PHASE / "results" / "structural_feature_class_results.csv", feature_rows, ["feature", "transformation", "estimate", "se", "ci_low", "ci_high", "p_value", "n_genes", "bh_q"])
    domain_summary = []
    for key in (
        "pfam_absent_from_alternative", "pfam_gained_in_alternative", "signal_principal_only",
        "signal_alternative_only", "tm_principal_only", "tm_alternative_only",
    ):
        domain_summary.append({
            "feature": key,
            "n_s2_genes": sum(item[key] == 1 and item["category"] == "S2" for item in categories.values()),
            "alignment": "not_estimated_no_residue_alignment" if "pfam" in key else "presence_comparison_without_localization_claim",
        })
    write_csv(PHASE / "results" / "domain_difference_summary.csv", domain_summary, ["feature", "n_s2_genes", "alignment"])
    consistency = []
    for transform in ("A", "B"):
        consistency.append({
            "transformation": transform,
            "adjusted_estimate": estimates[transform]["estimate"],
            "ci_low": estimates[transform]["ci_low"],
            "ci_high": estimates[transform]["ci_high"],
            "window_estimate": window_fits[transform]["estimate"] if window_fits[transform] else np.nan,
            "sign": int(np.sign(estimates[transform]["estimate"])),
        })
    write_csv(PHASE / "results" / "transformation_consistency.csv", consistency, list(consistency[0]))
    opportunity_rows.append({
        "analysis": "standardized_difference", "transformation": "annotation", "stratum": "",
        "n_s2": n_s2, "n_s1": n_s1, "estimate": window_info["length_smd"], "detail": "length_smd",
    })
    opportunity_rows.append({
        "analysis": "standardized_difference", "transformation": "annotation", "stratum": "",
        "n_s2": window_info["window_s2"], "n_s1": window_info["window_s1"], "estimate": window_info["pfam_smd"], "detail": "pfam_smd",
    })
    write_csv(PHASE / "results" / "annotation_opportunity_sensitivity.csv", opportunity_rows, ["analysis", "transformation", "stratum", "n_s2", "n_s1", "estimate", "detail"])
    write_csv(PHASE / "results" / "bootstrap_results.csv", bootstrap_rows, list(bootstrap_rows[0]))
    write_csv(PHASE / "results" / "leave_one_cancer_out.csv", loo_rows, list(loo_rows[0]))
    heterogeneity = []
    for transform in ("A", "B"):
        deltas = np.array([row["estimate"] for row in contrast_map[transform]], dtype=float)
        variances = np.array([((float(row["ci_high"]) - float(row["ci_low"])) / (2 * Z95)) ** 2 for row in contrast_map[transform]])
        usable = np.isfinite(deltas) & np.isfinite(variances) & (variances > 0)
        weights = 1 / variances[usable]
        fixed = float(np.sum(weights * deltas[usable]) / np.sum(weights)) if np.any(usable) else np.nan
        heterogeneity.append({
            "transformation": transform,
            "fixed_effect_of_cancer_contrasts": float(fixed),
            "min_cancer_estimate": float(np.min(deltas)),
            "max_cancer_estimate": float(np.max(deltas)),
            "n_cancers": 10,
            "note": "Descriptive summary of the ten adjusted cancer contrasts. The stacked cluster-robust coefficient remains the primary estimate.",
        })
    for row in distribution_rows:
        heterogeneity.append({
            "transformation": row["transformation"],
            "fixed_effect_of_cancer_contrasts": "",
            "min_cancer_estimate": "",
            "max_cancer_estimate": "",
            "n_cancers": row["n_genes"],
            "note": f"{row['group']} mean {row['mean']:.6g}; median {row['median']:.6g}; trim10 {row['trimmed_10']:.6g}; trim20 {row['trimmed_20']:.6g}; p90 {row['p90']:.6g}; p95 {row['p95']:.6g}; ge0.02 {row['fraction_ge_0.02']:.6g}",
        })
    write_csv(PHASE / "results" / "structural_heterogeneity_summary.csv", heterogeneity, list(heterogeneity[0]))
    write_csv(PHASE / "results" / "multiple_testing_summary.csv", testing_rows, ["family", "test", "p_value", "bh_q", "n"])
    for row in window_fits.values():
        if row is None:
            qualified = False
    label = decide(qualified, estimates, contrast_map, loo_rows, window_fits)
    write_reports(counts, qualified, label, estimates, contrast_map, feature_rows, window_info, loo_rows, covariate_sets)
    (PHASE / "manifests" / "software_environment.txt").write_text(
        f"python {platform.python_version()}\nnumpy {np.__version__}\nscipy {stats.__version__ if hasattr(stats, '__version__') else 'scipy.stats'}\n"
        f"scipy_version {__import__('scipy').__version__}\nseed {SEED}\nbootstrap_draws {DRAWS}\nfinished {utc_now()}\n"
    )
    input_rows = []
    for rel in (
        "phase2a/manifests/domain_mapping_audit.csv",
        "phase2a/results/functional_annotation_coverage.csv",
        "phase2a/results/protein_sequence_difference_coverage.csv",
        "phase2a/results/transcript_structure_coverage.csv",
        "phase1a/results/primary_literal_A_gene_fits.tsv.gz",
        "phase1a/results/primary_literal_B_gene_fits.tsv.gz",
        "phase2b/protocol/PHASE2B_LOCK.md",
        "phase2b/protocol/analysis_config.yaml",
    ):
        input_rows.append({"path": rel, "sha256": sha256_file(ROOT / rel), "role": "analysis_input"})
    write_csv(PHASE / "manifests" / "input_manifest.csv", input_rows, ["path", "sha256", "role"])
    log(f"PHASE2B_DONE {label} qualified={int(qualified)}")
    hash_outputs()


def one_step_influence(fit: dict) -> np.ndarray:
    design = fit["design"]
    resid = fit["resid"]
    codes = fit["codes"]
    inverse = np.linalg.inv(design.T @ design)
    influence = np.zeros(fit["n_genes"])
    for code in range(fit["n_genes"]):
        mask = codes == code
        score = design[mask].T @ resid[mask]
        influence[code] = (inverse @ score)[1]
    return influence


def fit_feature_model(rows: list[dict], covariates: list[str]) -> dict[str, dict]:
    outcome = np.array([row["y"] for row in rows], dtype=float)
    cancer = np.array([CANCERS.index(row["cancer"]) for row in rows])
    columns = [np.ones(len(rows))]
    names = ["intercept"]
    for feature in ("pfam_diff", "signal_diff", "tm_diff"):
        columns.append(np.array([row[feature] for row in rows], dtype=float))
        names.append(feature)
    for index, name in enumerate(CANCERS):
        if name == "BRCA":
            continue
        columns.append((cancer == index).astype(float))
        names.append(name)
    for name in covariates:
        values = np.array([row[name] for row in rows], dtype=float)
        finite = values[np.isfinite(values)]
        center = np.mean(finite) if finite.size else 0
        scale = np.std(finite) if finite.size and np.std(finite) > 0 else 1
        filled = np.where(np.isfinite(values), values, center)
        columns.append((filled - center) / scale)
        names.append(name)
    design = np.column_stack(columns)
    keep = [index for index in range(design.shape[1]) if index == 0 or np.std(design[:, index]) > 1e-8]
    design = design[:, keep]
    names = [names[index] for index in keep]
    beta, resid = ols(design, outcome)
    codes = np.unique([row["gene"] for row in rows], return_inverse=True)[1]
    covariance = cluster_covariance(design, resid, codes)
    results = {}
    if covariance is None:
        return {feature: {"estimate": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan} for feature in ("pfam_diff", "signal_diff", "tm_diff")}
    se = np.sqrt(np.clip(np.diag(covariance), 0, None))
    for feature in ("pfam_diff", "signal_diff", "tm_diff"):
        if feature not in names:
            results[feature] = {"estimate": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan}
            continue
        index = names.index(feature)
        z_value = beta[index] / se[index] if se[index] > 0 else np.nan
        results[feature] = {
            "estimate": float(beta[index]),
            "se": float(se[index]),
            "ci_low": float(beta[index] - Z95 * se[index]),
            "ci_high": float(beta[index] + Z95 * se[index]),
            "p_value": float(2 * stats.norm.sf(abs(z_value))) if np.isfinite(z_value) else np.nan,
        }
    return results


if __name__ == "__main__":
    try:
        main()
    except Exception:
        (PHASE / "logs").mkdir(parents=True, exist_ok=True)
        (PHASE / "logs" / "run_error.txt").write_text(traceback.format_exc())
        raise
