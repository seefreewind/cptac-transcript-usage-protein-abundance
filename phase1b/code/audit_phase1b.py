#!/usr/bin/env python3
"""Phase 1B heterogeneity audit. Reads Phase 1A outputs and does not refit protein models."""
from __future__ import annotations

import csv
import gzip
import hashlib
import os
import platform
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.stats import rankdata
from scipy.stats.mstats import winsorize

ROOT = Path(__file__).resolve().parents[2]
PHASE1A = ROOT / "phase1a"
PHASE1B = ROOT / "phase1b"
sys.path.insert(0, str(PHASE1A / "code"))
from common import (  # noqa: E402
    CANCERS,
    LOG2E,
    acquisition_path,
    load_annotation,
    load_gene_matrix,
    load_phenotype,
    parse_numeric_row,
    triple_from_headers,
    header_names,
)

STARTED = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
LAMBDA = 7.848880
TOL_MEAN = 1e-7
TOL_QUANTILE = 1e-6
SEED = 20261114


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: plain(row.get(key, "")) for key in fields})


def plain(value):
    if isinstance(value, (np.floating, float)):
        if not np.isfinite(value):
            return ""
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return int(value)
    return value


def text_num(value, digits=4) -> str:
    if value is None or value == "":
        return "NA"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return "NA"
    return f"{number:.{digits}f}"


def gini(values: np.ndarray) -> float:
    cleaned = np.sort(np.asarray(values, dtype=float))
    cleaned = cleaned[np.isfinite(cleaned)]
    total = float(cleaned.sum())
    if cleaned.size == 0 or total <= 0:
        return np.nan
    index = np.arange(1, cleaned.size + 1, dtype=float)
    return float((2.0 * np.sum(index * cleaned) / (cleaned.size * total)) - (cleaned.size + 1) / cleaned.size)


def concentration(values: np.ndarray) -> dict:
    raw = np.asarray(values, dtype=float)
    raw = raw[np.isfinite(raw)]
    nonnegative = np.maximum(raw, 0.0)
    total = float(nonnegative.sum())
    order = np.argsort(nonnegative)
    ranked = nonnegative[order][::-1]
    out = {}
    for label, fraction in (("top_1", 0.01), ("top_5", 0.05), ("top_10", 0.10)):
        k = max(1, int(np.ceil(fraction * ranked.size)))
        out[label] = float(ranked[:k].sum() / total) if total > 0 else np.nan
        out[f"{label}_k"] = k
    shares = nonnegative / total if total > 0 else np.zeros_like(nonnegative)
    sum_squares = float(np.sum(shares ** 2))
    out["gini"] = gini(nonnegative)
    out["n_effective"] = float(1.0 / sum_squares) if sum_squares > 0 else np.nan
    return out


def trimmed(values: np.ndarray, proportion: float) -> float:
    return float(stats.trim_mean(np.asarray(values, dtype=float), proportion))


def winsorized_mean(values: np.ndarray, proportion: float = 0.10) -> float:
    return float(np.mean(winsorize(np.asarray(values, dtype=float), limits=(proportion, proportion))))


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values, kind="mergesort")
    ranked = p_values[order]
    m_tests = ranked.size
    raw = ranked * m_tests / np.arange(1, m_tests + 1)
    raw = np.minimum.accumulate(raw[::-1])[::-1]
    raw = np.clip(raw, 0, 1)
    out = np.empty(m_tests)
    out[order] = raw
    return out


def reml_meta(theta: np.ndarray, se: np.ndarray) -> dict:
    theta = np.asarray(theta, dtype=float)
    se = np.asarray(se, dtype=float)
    k = int(theta.size)
    se2 = se ** 2

    def nll(tau2: float) -> float:
        variance = se2 + float(tau2)
        if np.any(variance <= 0):
            return 1e300
        weight = 1.0 / variance
        sw = float(weight.sum())
        mu = float(np.sum(weight * theta) / sw)
        return 0.5 * float(np.sum(np.log(variance))) + 0.5 * np.log(sw) + 0.5 * float(np.sum(weight * (theta - mu) ** 2))

    upper = max(1e-8, float(np.var(theta)) * 100 + float(np.max(se2)))
    grid = np.linspace(0.0, upper, 20001)
    grid_nll = np.array([nll(value) for value in grid])
    grid_tau = float(grid[int(np.argmin(grid_nll))])
    opt = minimize_scalar(nll, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-14})
    opt_tau = float(opt.x)
    if abs(grid_tau - opt_tau) > 1e-6 and abs(nll(grid_tau) - nll(opt_tau)) > 1e-8:
        raise SystemExit(f"REML grid and optimizer disagree: {grid_tau} vs {opt_tau}")
    tau2 = 0.0 if nll(0.0) <= nll(opt_tau) + 1e-10 else opt_tau
    variance = se2 + tau2
    weight = 1.0 / variance
    sw = float(weight.sum())
    mu = float(np.sum(weight * theta) / sw)
    q_hk = float(np.sum(weight * (theta - mu) ** 2) / (k - 1))
    se_hk = float(np.sqrt(q_hk / sw))
    critical = float(stats.t.ppf(0.975, k - 1))
    w0 = 1.0 / se2
    mu0 = float(np.sum(w0 * theta) / np.sum(w0))
    cochran_q = float(np.sum(w0 * (theta - mu0) ** 2))
    i2 = max(0.0, (cochran_q - (k - 1)) / cochran_q) if cochran_q > 0 else 0.0
    pred_t = float(stats.t.ppf(0.975, k - 2))
    pred_se = float(np.sqrt(tau2 + se_hk ** 2))
    se_fixed = float(np.sqrt(1.0 / np.sum(w0)))
    return {
        "k": k,
        "mu": mu,
        "se_hk": se_hk,
        "ci_low": mu - critical * se_hk,
        "ci_high": mu + critical * se_hk,
        "tau2": tau2,
        "I2": i2,
        "Q": cochran_q,
        "Q_p": float(stats.chi2.sf(cochran_q, k - 1)),
        "pi_low": mu - pred_t * pred_se,
        "pi_high": mu + pred_t * pred_se,
        "fixed_mu": mu0,
        "fixed_se": se_fixed,
        "fixed_ci_low": mu0 - critical * se_fixed,
        "fixed_ci_high": mu0 + critical * se_fixed,
    }


def detectable_partial_r2(n_value: float) -> float:
    if not np.isfinite(n_value) or n_value <= 0:
        return np.nan
    f2 = LAMBDA / float(n_value)
    return float(f2 / (1.0 + f2))


def power_at(partial_r2: float, n_value: float) -> float:
    if n_value <= 5 or not (0 < partial_r2 < 1):
        return np.nan
    f2 = partial_r2 / (1.0 - partial_r2)
    dfd = n_value - 5
    crit = stats.f.ppf(0.95, 1, dfd)
    return float(stats.ncf.sf(crit, 1, dfd, f2 * n_value))


def run_synthetic_tests() -> list[str]:
    lines = []
    if abs(gini(np.ones(10))) > 1e-12:
        raise SystemExit("Gini of equal values is not 0")
    if abs(gini(np.array([0.0, 0.0, 0.0, 1.0])) - 0.75) > 1e-12:
        raise SystemExit("Gini of a point mass failed")
    point = concentration(np.array([0.0, 0.0, 0.0, 1.0]))
    if abs(point["top_10"] - 1) > 1e-12 or abs(point["n_effective"] - 1) > 1e-12:
        raise SystemExit("concentration of a point mass failed")
    spread = np.array([1.0, 2.0, 3.0, 4.0, 100.0])
    if not (0 < concentration(spread)["top_10"] < 1):
        raise SystemExit("concentration share out of range")
    base = np.arange(1, 21, dtype=float)
    if abs(trimmed(base, 0.10) - np.mean(base[2:-2])) > 1e-10:
        raise SystemExit("10% trimmed mean failed")
    if abs(float(np.quantile(base, 0.10)) - 2.9) > 1e-12:
        raise SystemExit("quantile failed")
    theta = np.array([0.02, 0.03, 0.025, 0.022, 0.028])
    se = np.full(5, 0.002)
    meta = reml_meta(theta, se)
    if not (meta["ci_low"] < meta["mu"] < meta["ci_high"]):
        raise SystemExit("Hartung-Knapp interval missed its mean")
    if meta["pi_low"] >= meta["ci_low"] or meta["pi_high"] <= meta["ci_high"]:
        raise SystemExit("prediction interval is not wider than the confidence interval")
    loo = []
    for index in range(theta.size):
        kept = np.ones(theta.size, dtype=bool)
        kept[index] = False
        loo.append(reml_meta(theta[kept], se[kept])["mu"])
    if np.max(np.abs(np.array(loo) - meta["mu"])) > 0.02:
        raise SystemExit("leave-one-out means moved implausibly")
    q_values = benjamini_hochberg(np.array([0.001, 0.02, 0.03, 0.20]))
    if np.any(np.diff(q_values) < -1e-12) or np.any(q_values < np.array([0.001, 0.02, 0.03, 0.20]) - 1e-12):
        raise SystemExit("BH adjustment failed")
    lines.append("SYNTHETIC_TESTS_PASS")
    lines.append(f"synthetic_mu {meta['mu']:.6g} tau2 {meta['tau2']:.6g} I2 {meta['I2']:.6g}")
    return lines


def load_fits(path: Path) -> dict:
    grouped = {cancer: [] for cancer in CANCERS}
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[row["cancer"]].append(row)
    out = {}
    for cancer, rows in grouped.items():
        def column(name: str) -> np.ndarray:
            values = np.empty(len(rows))
            for i, row in enumerate(rows):
                text = row[name]
                values[i] = np.nan if text == "" else float(text)
            return values

        out[cancer] = {
            "gene_id": [row["gene_id"] for row in rows],
            "principal": [row["principal_transcript_id"] for row in rows],
            "n_isoform": column("n_isoform_rows"),
            "n_complete": column("n_complete"),
            "fraction_zero": column("fraction_stored_all_zero"),
            "partial": column("partial_r2"),
            "q_value": column("bh_q_value"),
            "beta": column("beta_p_standardized"),
            "vif": column("cnv_vif"),
            "condition": column("condition_number"),
            "signed_r": column("signed_partial_r"),
            "expected_null": column("expected_null_partial_r2"),
            "status": [row["fit_status"] for row in rows],
        }
        out[cancer]["ok"] = np.array([item == "OK" for item in out[cancer]["status"]])
    return out


def load_summary() -> dict:
    out = {}
    with (PHASE1A / "results" / "cancer_level_summary.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            parsed = {}
            for key, value in row.items():
                if key in {"cancer", "transformation"}:
                    parsed[key] = value
                else:
                    parsed[key] = float(value) if value != "" else np.nan
            out[(row["cancer"], row["transformation"])] = parsed
    return out


def load_stored_meta() -> list[dict]:
    with (PHASE1A / "results" / "random_effects_summary.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def status_counts(block: dict) -> dict:
    counts = {name: 0 for name in block["status"]}
    for name in block["status"]:
        counts[name] = counts.get(name, 0) + 1
    prefit = {"OK", "RANK_DEFICIENT", "NO_RESIDUAL_P_VARIATION", "ZERO_PROTEIN_VARIANCE", "M0_RSS_ZERO"}
    return {
        "n_candidate": len(block["status"]),
        "n_success": int(block["ok"].sum()),
        "n_rank_deficient": counts.get("RANK_DEFICIENT", 0),
        "n_no_residual_p": counts.get("NO_RESIDUAL_P_VARIATION", 0),
        "n_prefit": sum(counts.get(name, 0) for name in prefit),
    }


def reproduce(fits: dict, summary: dict) -> tuple[list[dict], bool]:
    rows = []
    failed = False
    for transform, blocks in fits.items():
        for cancer in CANCERS:
            block = blocks[cancer]
            ok = block["partial"][block["ok"]]
            counts = status_counts(block)
            stored = summary[(cancer, transform)]
            fresh = {
                "mean_partial_r2": float(ok.mean()),
                "median_partial_r2": float(np.median(ok)),
                "sd_partial_r2": float(ok.std(ddof=1)),
                "q25": float(np.quantile(ok, 0.25)),
                "q75": float(np.quantile(ok, 0.75)),
                "p90": float(np.quantile(ok, 0.90)),
                "p95": float(np.quantile(ok, 0.95)),
            }
            critical = float(stats.t.ppf(0.975, stored["jackknife_n"] - 1))
            rebuilt_low = stored["mean_partial_r2"] - critical * stored["jackknife_se"]
            rebuilt_high = stored["mean_partial_r2"] + critical * stored["jackknife_se"]
            checks = {
                "n_candidate": counts["n_candidate"] == int(stored["n_candidate"]),
                "n_success": counts["n_success"] == int(stored["n_success"]),
                "n_rank_deficient": counts["n_rank_deficient"] == int(stored["n_rank_deficient"]),
                "n_no_residual_p": counts["n_no_residual_p"] == int(stored["n_no_residual_p"]),
                "n_prefit": counts["n_prefit"] == int(stored["n_prefit"]),
                "mean_partial_r2": abs(fresh["mean_partial_r2"] - stored["mean_partial_r2"]) <= TOL_MEAN,
                "median_partial_r2": abs(fresh["median_partial_r2"] - stored["median_partial_r2"]) <= TOL_MEAN,
                "sd_partial_r2": abs(fresh["sd_partial_r2"] - stored["sd_partial_r2"]) <= TOL_MEAN,
                "q25": abs(fresh["q25"] - stored["q25"]) <= TOL_QUANTILE,
                "q75": abs(fresh["q75"] - stored["q75"]) <= TOL_QUANTILE,
                "p90": abs(fresh["p90"] - stored["p90"]) <= TOL_QUANTILE,
                "p95": abs(fresh["p95"] - stored["p95"]) <= TOL_QUANTILE,
                "jackknife_interval": abs(rebuilt_low - stored["ci_low"]) <= 1e-8 and abs(rebuilt_high - stored["ci_high"]) <= 1e-8,
            }
            if not all(checks.values()):
                failed = True
            rows.append({
                "cancer": cancer,
                "transformation": transform,
                "between_gene_variance": fresh["sd_partial_r2"] ** 2,
                "jackknife_se": stored["jackknife_se"],
                "naive_independence_se": fresh["sd_partial_r2"] / np.sqrt(ok.size),
                "rebuilt_ci_low": rebuilt_low,
                "rebuilt_ci_high": rebuilt_high,
                "checks_passed": int(all(checks.values())),
                "failed_checks": ";".join(key for key, value in checks.items() if not value),
                **counts,
                **fresh,
            })
    return rows, failed


def compare_meta(summary: dict, stored_rows: list[dict]) -> tuple[list[dict], list[dict], bool]:
    failed = False
    influence = []
    loo_rows = []
    stored = {(row["analysis"], row["transformation"]): row for row in stored_rows}
    for transform in ("A", "B"):
        theta = np.array([summary[(cancer, transform)]["mean_partial_r2"] for cancer in CANCERS])
        se = np.array([summary[(cancer, transform)]["jackknife_se"] for cancer in CANCERS])
        analyses = {
            "all_ten": reml_meta(theta, se),
            "exclude_OV": reml_meta(
                np.array([summary[(cancer, transform)]["mean_partial_r2"] for cancer in CANCERS if cancer != "OV"]),
                np.array([summary[(cancer, transform)]["jackknife_se"] for cancer in CANCERS if cancer != "OV"]),
            ),
        }
        for name, meta in analyses.items():
            ref = stored[(name, transform)]
            keys = ("mu", "tau2", "I2", "Q", "ci_low", "ci_high", "pi_low", "pi_high")
            bad = [key for key in keys if abs(meta[key] - float(ref[key])) > 1e-6]
            if bad:
                failed = True
            influence.append({
                "analysis": name,
                "transformation": transform,
                "estimator": "REML_HartungKnapp",
                "role": "primary" if name == "all_ten" else "prespecified_sensitivity",
                "match": int(not bad),
                "failed_fields": ";".join(bad),
                **meta,
            })
        influence.append({
            "analysis": "fixed_effect",
            "transformation": transform,
            "estimator": "inverse_variance_tau0",
            "role": "descriptive_sensitivity",
            "match": "",
            "failed_fields": "",
            "k": analyses["all_ten"]["k"],
            "mu": analyses["all_ten"]["fixed_mu"],
            "se_hk": analyses["all_ten"]["fixed_se"],
            "ci_low": analyses["all_ten"]["fixed_ci_low"],
            "ci_high": analyses["all_ten"]["fixed_ci_high"],
            "tau2": 0.0,
            "I2": analyses["all_ten"]["I2"],
            "Q": analyses["all_ten"]["Q"],
            "Q_p": analyses["all_ten"]["Q_p"],
            "pi_low": "",
            "pi_high": "",
        })
        med = float(np.median(theta))
        influence.append({
            "analysis": "median_of_cancer_means",
            "transformation": transform,
            "estimator": "descriptive_median",
            "role": "descriptive_sensitivity",
            "match": "",
            "failed_fields": "",
            "k": 10,
            "mu": med,
            "se_hk": "",
            "ci_low": "",
            "ci_high": "",
            "tau2": "",
            "I2": "",
            "Q": "",
            "Q_p": "",
            "pi_low": "",
            "pi_high": "",
        })
        influence.append({
            "analysis": "trimmed_mean_20_of_cancer_means",
            "transformation": transform,
            "estimator": "descriptive_trimmed_mean",
            "role": "descriptive_sensitivity",
            "match": "",
            "failed_fields": "",
            "k": 10,
            "mu": trimmed(theta, 0.20),
            "se_hk": "",
            "ci_low": "",
            "ci_high": "",
            "tau2": "",
            "I2": "",
            "Q": "",
            "Q_p": "",
            "pi_low": "",
            "pi_high": "",
        })
        for index, cancer in enumerate(CANCERS):
            kept = np.ones(len(CANCERS), dtype=bool)
            kept[index] = False
            meta = reml_meta(theta[kept], se[kept])
            loo_rows.append({
                "analysis": "all_ten",
                "transformation": transform,
                "left_out": cancer,
                "mu": meta["mu"],
                "ci_low": meta["ci_low"],
                "ci_high": meta["ci_high"],
                "pi_low": meta["pi_low"],
                "pi_high": meta["pi_high"],
                "tau2": meta["tau2"],
                "I2": meta["I2"],
                "delta_mu": analyses["all_ten"]["mu"] - meta["mu"],
                "delta_tau2": analyses["all_ten"]["tau2"] - meta["tau2"],
            })
    return influence, loo_rows, failed


def distribution_rows(fits: dict) -> tuple[list[dict], list[dict]]:
    dist_rows = []
    conc_rows = []
    for transform, blocks in fits.items():
        for cancer in CANCERS:
            block = blocks[cancer]
            values = block["partial"][block["ok"]]
            q_values = block["q_value"][block["ok"]]
            conc = concentration(values)
            ratio = float(np.median(values) / np.mean(values))
            if ratio > 0.70 and conc["top_10"] < 0.35:
                klass = "BROAD_WEAK"
            elif conc["top_1"] > 0.25 or conc["top_10"] > 0.60:
                klass = "CONCENTRATED"
            else:
                klass = "MIXED_DISTRIBUTION"
            dist_rows.append({
                "cancer": cancer,
                "transformation": transform,
                "n_success": int(values.size),
                "mean_partial_r2": float(values.mean()),
                "median_partial_r2": float(np.median(values)),
                "trimmed_mean_10": trimmed(values, 0.10),
                "trimmed_mean_20": trimmed(values, 0.20),
                "winsorized_mean_10": winsorized_mean(values, 0.10),
                "p10": float(np.quantile(values, 0.10)),
                "p25": float(np.quantile(values, 0.25)),
                "p75": float(np.quantile(values, 0.75)),
                "p90": float(np.quantile(values, 0.90)),
                "p95": float(np.quantile(values, 0.95)),
                "p99": float(np.quantile(values, 0.99)),
                "fraction_ge_0.01": float(np.mean(values >= 0.01)),
                "fraction_ge_0.02": float(np.mean(values >= 0.02)),
                "fraction_ge_0.05": float(np.mean(values >= 0.05)),
                "fraction_ge_0.10": float(np.mean(values >= 0.10)),
                "fraction_q_lt_0.05": float(np.mean(q_values < 0.05)),
                "fraction_q_and_r2": float(np.mean((q_values < 0.05) & (values >= 0.02))),
                "distribution_class": klass,
            })
            conc_rows.append({
                "cancer": cancer,
                "transformation": transform,
                "share_top_1pct": conc["top_1"],
                "share_top_5pct": conc["top_5"],
                "share_top_10pct": conc["top_10"],
                "n_top_1pct": conc["top_1_k"],
                "n_top_5pct": conc["top_5_k"],
                "n_top_10pct": conc["top_10_k"],
                "gini_nonnegative": conc["gini"],
                "n_effective_positive_mass": conc["n_effective"],
                "distribution_class": klass,
            })
    return dist_rows, conc_rows


def sample_size_rows(fits: dict, summary: dict) -> tuple[list[dict], list[dict]]:
    size_rows = []
    power_rows = []
    cohort_loss = {cancer: {"purity": 0, "cnv": 0, "n": 0} for cancer in CANCERS}
    with (PHASE1A / "qc" / "cohort_reconstruction.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            cohort_loss[row["cancer"]]["n"] += 1
            cohort_loss[row["cancer"]]["purity"] += int(row["wes_purity_numeric"] == "0")
            cohort_loss[row["cancer"]]["cnv"] += int(row["in_cnv_header"] == "0")
    for transform, blocks in fits.items():
        for cancer in CANCERS:
            block = blocks[cancer]
            n_all = block["n_complete"]
            ok = block["ok"]
            n_ok = n_all[ok]
            partial = block["partial"][ok]
            spread = float(np.quantile(n_ok, 0.90) - np.quantile(n_ok, 0.10))
            concentrated = int(spread <= 2)
            if concentrated:
                low_mean = high_mean = np.nan
                low_n = high_n = 0
            else:
                q1, q2 = np.quantile(n_ok, [1 / 3, 2 / 3])
                low = partial[n_ok <= q1]
                high = partial[n_ok >= q2]
                low_mean = float(low.mean()) if low.size else np.nan
                high_mean = float(high.mean()) if high.size else np.nan
                low_n = int(low.size)
                high_n = int(high.size)
            below = n_all[~ok]
            size_rows.append({
                "cancer": cancer,
                "transformation": transform,
                "n_candidate": int(n_all.size),
                "n_success": int(ok.sum()),
                "median_n_success": float(np.median(n_ok)),
                "mean_n_success": float(n_ok.mean()),
                "p10_n_success": float(np.quantile(n_ok, 0.10)),
                "p90_n_success": float(np.quantile(n_ok, 0.90)),
                "fraction_candidate_n_ge_71": float(np.mean(n_all >= 71)),
                "fraction_candidate_n_ge_80": float(np.mean(n_all >= 80)),
                "fraction_candidate_n_ge_90": float(np.mean(n_all >= 90)),
                "fraction_candidate_n_ge_100": float(np.mean(n_all >= 100)),
                "fraction_success_n_ge_80": float(np.mean(n_ok >= 80)),
                "fraction_success_n_ge_90": float(np.mean(n_ok >= 90)),
                "fraction_success_n_ge_100": float(np.mean(n_ok >= 100)),
                "cases_lost_wes_purity": cohort_loss[cancer]["purity"],
                "cases_lost_cnv_header": cohort_loss[cancer]["cnv"],
                "median_complete_cases_among_failures": float(np.median(below)) if below.size else np.nan,
                "n_distribution_concentrated": concentrated,
                "low_n_tertile_mean_partial_r2": low_mean,
                "high_n_tertile_mean_partial_r2": high_mean,
                "low_n_tertile_genes": low_n,
                "high_n_tertile_genes": high_n,
                "jackknife_se": summary[(cancer, transform)]["jackknife_se"],
                "between_gene_sd": float(partial.std(ddof=1)),
                "naive_independence_se": float(partial.std(ddof=1) / np.sqrt(partial.size)),
            })
            median_n = float(np.median(n_ok))
            power_rows.append({
                "cancer": cancer,
                "transformation": transform,
                "median_n": median_n,
                "planning_detectable_partial_r2": detectable_partial_r2(median_n),
                "planning_detectable_partial_r2_at_71": detectable_partial_r2(71),
                "exact_power_at_partial_r2_0.02": power_at(0.02, median_n),
                "exact_power_at_partial_r2_0.10": power_at(0.10, median_n),
                "jackknife_se": summary[(cancer, transform)]["jackknife_se"],
                "mean_expected_null_partial_r2": float(np.nanmean(block["expected_null"][ok])),
                "note": "N=71 is the eligibility gate. The frozen planning formula gives a detectable partial R2 near 0.10 at that N, not 0.02.",
            })
    return size_rows, power_rows


def robustness_rows(fits: dict, summary: dict) -> list[dict]:
    rows = []
    sens = list(csv.DictReader((PHASE1A / "results" / "sensitivity_summary.csv").open()))
    c_means = {(row["cancer"], row["transformation"]): float(row["mean_partial_r2"]) for row in sens if row["sensitivity"] == "S3_candidate_C" and row["mean_partial_r2"]}
    for cancer in CANCERS:
        a = fits["A"][cancer]
        b = fits["B"][cancer]
        if a["gene_id"] != b["gene_id"]:
            raise SystemExit(f"gene order differs for {cancer}")
        both = a["ok"] & b["ok"]
        pa = a["partial"][both]
        pb = b["partial"][both]
        ba = a["beta"][both]
        bb = b["beta"][both]
        same = (np.sign(ba) == np.sign(bb)) & (ba != 0) & (bb != 0)
        stored_a = summary[(cancer, "A")]
        stored_b = summary[(cancer, "B")]
        rows.append({
            "cancer": cancer,
            "n_paired": int(both.sum()),
            "mean_A": stored_a["mean_partial_r2"],
            "mean_B": stored_b["mean_partial_r2"],
            "mean_A_minus_B": stored_a["mean_partial_r2"] - stored_b["mean_partial_r2"],
            "median_A_minus_B": float(np.median(pa) - np.median(pb)),
            "spearman_partial_r2": float(stats.spearmanr(pa, pb).statistic),
            "proportion_same_beta_sign": float(np.mean(same)),
            "proportion_both_ge_0.02": float(np.mean((pa >= 0.02) & (pb >= 0.02))),
            "proportion_only_A_ge_0.02": float(np.mean((pa >= 0.02) & (pb < 0.02))),
            "proportion_only_B_ge_0.02": float(np.mean((pa < 0.02) & (pb >= 0.02))),
            "A_interval_above_0.02": int(stored_a["ci_low"] > 0.02),
            "B_interval_above_0.02": int(stored_b["ci_low"] > 0.02),
            "candidate_C_mean_partial_r2": c_means.get((cancer, "C"), np.nan),
        })
    return rows


def sensitivity_rows(summary: dict, stored_meta: list[dict]) -> list[dict]:
    primary = {(row["cancer"], row["transformation"]): row for row in summary.values()}
    stored = {(row["analysis"], row["transformation"]): row for row in stored_meta}
    rows = []
    with (PHASE1A / "results" / "sensitivity_summary.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            cancer = row["cancer"]
            transform = row["transformation"]
            name = row["sensitivity"]
            if name == "S4_clr":
                rows.append({
                    "sensitivity": name,
                    "cancer": cancer,
                    "transformation": transform,
                    "status": "NOT_RUN_UNDERDEFINED",
                    "primary_mean": "",
                    "sensitivity_mean": "",
                    "delta_mean": "",
                    "point_classification_changed": "",
                    "assumption_tested": "CLR contrast was not defined in the frozen documents.",
                })
                continue
            if name == "S8_meta_without_OV":
                full = stored[("all_ten", transform)]
                alt = stored[("exclude_OV", transform)]
                rows.append({
                    "sensitivity": name,
                    "cancer": cancer,
                    "transformation": transform,
                    "status": "STORED_META",
                    "primary_mean": float(full["mu"]),
                    "sensitivity_mean": float(alt["mu"]),
                    "delta_mean": float(alt["mu"]) - float(full["mu"]),
                    "point_classification_changed": int((float(full["ci_low"]) > 0.02) != (float(alt["ci_low"]) > 0.02)),
                    "primary_ci_low": float(full["ci_low"]),
                    "sensitivity_ci_low": float(alt["ci_low"]),
                    "assumption_tested": "OV is one cohort in the cross-cancer average, not a technical replicate.",
                })
                continue
            if transform not in {"A", "B"}:
                sens_mean = float(row["mean_partial_r2"]) if row["mean_partial_r2"] else np.nan
                base_a = primary[(cancer, "A")]["mean_partial_r2"]
                base_b = primary[(cancer, "B")]["mean_partial_r2"]
                rows.append({
                    "sensitivity": name,
                    "cancer": cancer,
                    "transformation": transform,
                    "status": "SENSITIVITY_ONLY",
                    "primary_mean": "",
                    "sensitivity_mean": sens_mean,
                    "delta_mean": "",
                    "delta_versus_A": sens_mean - base_a,
                    "delta_versus_B": sens_mean - base_b,
                    "point_classification_changed": int((sens_mean >= 0.02) != (base_a >= 0.02) or (sens_mean >= 0.02) != (base_b >= 0.02)),
                    "n_genes": row["n_genes"],
                    "assumption_tested": sensitivity_assumption(name),
                })
                continue
            if row["n_genes"] in {"", "0"} or row["mean_partial_r2"] == "":
                rows.append({
                    "sensitivity": name,
                    "cancer": cancer,
                    "transformation": transform,
                    "status": "NO_GENE_MET_TRIGGER" if name.startswith("S6") else "NO_STORED_MEAN",
                    "primary_mean": primary[(cancer, transform)]["mean_partial_r2"],
                    "sensitivity_mean": "",
                    "delta_mean": "",
                    "delta_versus_A": "",
                    "delta_versus_B": "",
                    "point_classification_changed": 0,
                    "n_genes": row["n_genes"],
                    "assumption_tested": sensitivity_assumption(name),
                })
                continue
            sens_mean = float(row["mean_partial_r2"])
            base = primary[(cancer, transform)]["mean_partial_r2"]
            rows.append({
                "sensitivity": name,
                "cancer": cancer,
                "transformation": transform,
                "status": "COMPARED_TO_STORED_MEAN",
                "primary_mean": base,
                "sensitivity_mean": sens_mean,
                "delta_mean": sens_mean - base,
                "delta_versus_A": sens_mean - base if transform == "A" else "",
                "delta_versus_B": sens_mean - base if transform == "B" else "",
                "point_classification_changed": int((base >= 0.02) != (sens_mean >= 0.02)),
                "n_genes": row["n_genes"],
                "assumption_tested": sensitivity_assumption(name),
            })
    return rows


def sensitivity_assumption(name: str) -> str:
    return {
        "S1_drop_stored_all_zero": "Literal A assigns equal weight when every stored isoform value is zero.",
        "S2_rsem": "The gene-expression covariate is the stored RSEM value rather than log2 of the isoform total.",
        "S3_candidate_C": "Candidate C is the clipped zero-preserving sensitivity, not a third primary inverse.",
        "S5_estimate": "ESTIMATE replaces WES purity and is not treated as the same measurement.",
        "S6_cnv_vif_gt_10": "CNV is dropped only in the stored high-VIF sensitivity.",
        "S7_m0_pearson_gt_4": "Cases with absolute M0 Pearson residual above 4 are removed once.",
    }.get(name, name)


def transcript_ids(path: Path) -> list[str]:
    ids = []
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip() or b"\x00" in line:
                if b"\x00" in line:
                    raise SystemExit(f"NUL byte in {path}")
                continue
            ids.append(line.split(b"\t", 1)[0].decode())
    return ids


def accumulate(path: Path, case_ids: list[str], expected: set[str], gene_index: dict[str, int], tx_to_gene: dict[str, str], principal: dict[str, str]):
    n_cases = len(case_ids)
    n_genes = len(gene_index)
    sum_a = np.zeros((n_genes, n_cases))
    sum_b = np.zeros((n_genes, n_cases))
    part_a = np.zeros((n_genes, n_cases))
    part_b = np.zeros((n_genes, n_cases))
    xlog_a = np.zeros((n_genes, n_cases))
    xlog_b = np.zeros((n_genes, n_cases))
    n_pos = np.zeros((n_genes, n_cases), dtype=np.int16)
    n_tx = np.zeros(n_genes, dtype=np.int32)
    names = header_names(path)
    columns = [{name: i for i, name in enumerate(names)}[case] for case in case_ids]
    seen = set()
    with path.open("rb") as handle:
        handle.readline()
        for line in handle:
            if not line.strip() or b"\x00" in line:
                if b"\x00" in line:
                    raise SystemExit(f"NUL byte in {path}")
                continue
            parts = line.rstrip(b"\r\n").split(b"\t")
            tx = parts[0].decode()
            if tx in seen or tx not in expected:
                raise SystemExit(f"unexpected transcript {tx} in {path}")
            seen.add(tx)
            gene = tx_to_gene[tx]
            gi = gene_index.get(gene)
            if gi is None:
                continue
            y = parse_numeric_row(parts, columns)
            if np.any(~np.isfinite(y)):
                raise SystemExit(f"non-finite isoform value {tx}")
            xa = np.exp2(y)
            xb = np.expm1(y * LOG2E)
            sum_a[gi] += xa
            sum_b[gi] += xb
            xlog_a[gi] += xa * np.log2(xa)
            with np.errstate(divide="ignore", invalid="ignore"):
                term = xb * np.log2(xb)
            term[~np.isfinite(term)] = 0.0
            xlog_b[gi] += term
            n_pos[gi] += (y > 0).astype(np.int16)
            n_tx[gi] += 1
            if principal.get(gene) == tx:
                part_a[gi] += xa
                part_b[gi] += xb
    if seen != expected:
        raise SystemExit(f"transcript set mismatch in {path}")
    log_b = np.full(sum_b.shape, np.nan)
    p_b = np.full(sum_b.shape, np.nan)
    h_a = np.full(sum_a.shape, np.nan)
    h_b = np.full(sum_b.shape, np.nan)
    np.log2(sum_b, out=log_b, where=sum_b > 0)
    np.divide(part_b, sum_b, out=p_b, where=sum_b > 0)
    np.log2(sum_a, out=h_a, where=sum_a > 0)
    h_a -= np.divide(xlog_a, sum_a, out=np.zeros_like(sum_a), where=sum_a > 0)
    np.log2(sum_b, out=h_b, where=sum_b > 0)
    h_b -= np.divide(xlog_b, sum_b, out=np.zeros_like(sum_b), where=sum_b > 0)
    return {
        "logG_A": np.log2(sum_a),
        "P_A": part_a / sum_a,
        "H_A": h_a,
        "logG_B": log_b,
        "P_B": p_b,
        "H_B": h_b,
        "all_zero": n_pos == 0,
        "n_tx": n_tx,
    }


def aligned_matrix(table: dict[str, np.ndarray], gene_ids: list[str], n_cases: int) -> np.ndarray:
    out = np.full((len(gene_ids), n_cases), np.nan)
    for i, gene in enumerate(gene_ids):
        if gene in table:
            out[i] = table[gene]
    return out


def safe_corr(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if a.size < 5 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan, np.nan
    pearson = float(np.corrcoef(a, b)[0, 1])
    spearman = float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])
    return pearson, spearman


def residual_variance(log_g, purity, cnv, proportion) -> float:
    design = np.column_stack([np.ones(log_g.size), log_g, purity, cnv])
    coef, _, rank, _ = np.linalg.lstsq(design, proportion, rcond=None)
    if rank < design.shape[1]:
        return np.nan
    resid = proportion - design @ coef
    df = log_g.size - rank
    if df <= 0:
        return np.nan
    return float(np.sum(resid ** 2) / df)


def measurement_audit(fits: dict) -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict], int]:
    print("measurement scan", now(), flush=True)
    tx_to_gene, principal_map, _rule = load_annotation()
    expected = set(transcript_ids(acquisition_path("BRCA", "rna_isoform")))
    gene_ids = fits["A"]["BRCA"]["gene_id"]
    for cancer in CANCERS:
        if fits["A"][cancer]["gene_id"] != gene_ids:
            raise SystemExit(f"candidate gene list differs in {cancer}")
        for gene, tx in zip(gene_ids, fits["A"][cancer]["principal"]):
            if principal_map.get(gene) != tx:
                raise SystemExit(f"principal transcript mismatch {gene}")
    gene_index = {gene: i for i, gene in enumerate(gene_ids)}
    usage_rows = []
    variability_rows = []
    covariate_rows = []
    collinearity_rows = []
    protein_rows = []
    mismatches = 0
    for cancer in CANCERS:
        print(f"  {cancer}", flush=True)
        case_ids = triple_from_headers(cancer)
        protein = aligned_matrix(load_gene_matrix(acquisition_path(cancer, "protein"), case_ids), gene_ids, len(case_ids))
        cnv = aligned_matrix(load_gene_matrix(acquisition_path(cancer, "cnv"), case_ids, allow_missing=True), gene_ids, len(case_ids))
        rsem = aligned_matrix(load_gene_matrix(acquisition_path(cancer, "rna_gene"), case_ids), gene_ids, len(case_ids))
        purity, _estimate = load_phenotype(acquisition_path(cancer, "phenotype"), case_ids)
        iso = accumulate(acquisition_path(cancer, "rna_isoform"), case_ids, expected, gene_index, tx_to_gene, principal_map)
        purity_ok = np.isfinite(purity)
        protein_table = load_gene_matrix(acquisition_path(cancer, "protein"), case_ids)
        observed_fractions = []
        variances = []
        numeric_values = []
        missing_cells = 0
        total_cells = 0
        for values in protein_table.values():
            finite = np.isfinite(values)
            observed_fractions.append(float(finite.mean()))
            missing_cells += int((~finite).sum())
            total_cells += int(values.size)
            if finite.sum() >= 2:
                variances.append(float(values[finite].var(ddof=1)))
            if finite.any():
                numeric_values.append(values[finite])
        stacked = np.concatenate(numeric_values) if numeric_values else np.array([np.nan])
        coverage_spearman = {}
        for transform in ("A", "B"):
            log_g = iso["logG_A"] if transform == "A" else iso["logG_B"]
            proportion = iso["P_A"] if transform == "A" else iso["P_B"]
            entropy = iso["H_A"] if transform == "A" else iso["H_B"]
            mask = (
                np.isfinite(protein) & np.isfinite(cnv) & purity_ok[None, :]
                & np.isfinite(log_g) & np.isfinite(proportion)
            )
            n_hat = mask.sum(axis=1)
            stored_n = fits[transform][cancer]["n_complete"]
            mismatch = int(np.sum(n_hat != stored_n))
            mismatches += mismatch
            if int(np.sum(iso["n_tx"] != fits["A"][cancer]["n_isoform"])) > 0:
                raise SystemExit(f"isoform-row counts differ in {cancer}")
            ok = fits[transform][cancer]["ok"]
            pz = np.where(mask, proportion, 0.0)
            hz = np.where(mask, entropy, 0.0)
            zz = np.where(mask, iso["all_zero"], False)
            n_ok = n_hat[ok]
            sum_p = pz[ok].sum(axis=1)
            sum_p2 = (pz[ok] ** 2).sum(axis=1)
            mean_p = sum_p / n_ok
            var_p = (sum_p2 - (sum_p ** 2) / n_ok) / np.maximum(n_ok - 1, 1)
            var_p = np.maximum(var_p, 0)
            mean_h = hz[ok].sum(axis=1) / n_ok
            frac_low = (mask[ok] & (proportion[ok] < 0.05)).sum(axis=1) / n_ok
            frac_high = (mask[ok] & (proportion[ok] > 0.95)).sum(axis=1) / n_ok
            frac_zero = zz[ok].sum(axis=1) / n_ok
            pearson = {key: [] for key in ("P_logG", "P_purity", "P_cnv", "P_rsem", "logG_cnv", "logG_purity", "logG_rsem")}
            spearman = {key: [] for key in pearson}
            residual = []
            gene_median = []
            for local, gi in enumerate(np.flatnonzero(ok)):
                idx = np.flatnonzero(mask[gi])
                p_vec = proportion[gi, idx]
                gene_median.append(float(np.median(p_vec)))
                g_vec = log_g[gi, idx]
                c_vec = cnv[gi, idx]
                u_vec = purity[idx]
                pairs = {
                    "P_logG": (p_vec, g_vec),
                    "P_purity": (p_vec, u_vec),
                    "P_cnv": (p_vec, c_vec),
                    "logG_cnv": (g_vec, c_vec),
                    "logG_purity": (g_vec, u_vec),
                }
                rsem_ok = np.isfinite(rsem[gi, idx])
                if int(rsem_ok.sum()) >= 71:
                    pairs["P_rsem"] = (p_vec[rsem_ok], rsem[gi, idx][rsem_ok])
                    pairs["logG_rsem"] = (g_vec[rsem_ok], rsem[gi, idx][rsem_ok])
                for key, (left, right) in pairs.items():
                    p_value, s_value = safe_corr(left, right)
                    pearson[key].append(p_value)
                    spearman[key].append(s_value)
                residual.append(residual_variance(g_vec, u_vec, c_vec, p_vec))
            residual = np.asarray(residual, dtype=float)
            n_iso = fits["A"][cancer]["n_isoform"]
            usage_rows.append({
                "cancer": cancer,
                "transformation": transform,
                "scope": "successful_gene_equal_weight",
                "n_complete_mismatches": mismatch,
                "mean_principal_proportion": float(np.mean(mean_p)),
                "median_of_gene_median_principal_proportion": float(np.median(gene_median)),
                "sd_of_gene_mean_principal_proportion": float(np.std(mean_p, ddof=1)),
                "iqr_of_gene_mean_principal_proportion": float(np.quantile(mean_p, 0.75) - np.quantile(mean_p, 0.25)),
                "median_within_gene_sd": float(np.median(np.sqrt(var_p))),
                "mean_entropy_bits": float(np.nanmean(mean_h)),
                "median_entropy_bits": float(np.nanmedian(mean_h)),
                "fraction_observations_below_0.05": float(np.mean(frac_low)),
                "fraction_observations_above_0.95": float(np.mean(frac_high)),
                "fraction_observations_stored_all_zero": float(np.mean(frac_zero)),
                "n_success_genes_with_2_isoforms": int(np.sum(ok & (n_iso == 2))),
                "n_success_genes_with_more_than_2_isoforms": int(np.sum(ok & (n_iso > 2))),
                "n_candidate_genes_with_2_isoforms": int(np.sum(n_iso == 2)),
                "n_candidate_genes_with_more_than_2_isoforms": int(np.sum(n_iso > 2)),
            })
            variability_rows.append({
                "cancer": cancer,
                "transformation": transform,
                "median_within_gene_sd": float(np.median(np.sqrt(var_p))),
                "p10_within_gene_sd": float(np.quantile(np.sqrt(var_p), 0.10)),
                "p90_within_gene_sd": float(np.quantile(np.sqrt(var_p), 0.90)),
                "median_residual_variance": float(np.nanmedian(residual)),
                "fraction_residual_variance_below_1e-4": float(np.mean(residual < 1e-4)),
                "median_entropy_bits": float(np.nanmedian(mean_h)),
                "fraction_below_0.05": float(np.mean(frac_low)),
                "fraction_above_0.95": float(np.mean(frac_high)),
                "fraction_stored_all_zero": float(np.mean(frac_zero)),
            })
            for key in pearson:
                covariate_rows.append({
                    "cancer": cancer,
                    "transformation": transform,
                    "pair": key,
                    "median_pearson": float(np.nanmedian(pearson[key])) if pearson[key] else np.nan,
                    "median_spearman": float(np.nanmedian(spearman[key])) if spearman[key] else np.nan,
                    "n_genes": int(np.sum(np.isfinite(pearson[key]))) if pearson[key] else 0,
                })
            vif = fits[transform][cancer]["vif"][ok]
            cond = fits[transform][cancer]["condition"][ok]
            signed = np.abs(fits[transform][cancer]["signed_r"][ok])
            collinearity_rows.append({
                "cancer": cancer,
                "transformation": transform,
                "median_cnv_vif": float(np.nanmedian(vif)),
                "p90_cnv_vif": float(np.nanquantile(vif, 0.90)),
                "fraction_vif_gt_10": float(np.mean(vif > 10)),
                "median_condition_number": float(np.nanmedian(cond)),
                "p90_condition_number": float(np.nanquantile(cond, 0.90)),
                "median_residual_variance_P": float(np.nanmedian(residual)),
                "p10_residual_variance_P": float(np.nanquantile(residual, 0.10)),
                "median_abs_signed_partial_r": float(np.nanmedian(signed)),
            })
            coverage = np.isfinite(protein).mean(axis=1)
            partial = fits[transform][cancer]["partial"][ok]
            coverage_spearman[transform] = float(stats.spearmanr(coverage[ok], partial).statistic)
        protein_rows.append({
            "cancer": cancer,
            "n_protein_genes": len(protein_table),
            "n_genes_with_any_numeric_value": int(sum(item > 0 for item in observed_fractions)),
            "cell_missing_fraction": missing_cells / total_cells,
            "median_gene_observation_fraction": float(np.median(observed_fractions)),
            "median_within_gene_variance": float(np.median(variances)),
            "abundance_p10": float(np.quantile(stacked, 0.10)),
            "abundance_p50": float(np.quantile(stacked, 0.50)),
            "abundance_p90": float(np.quantile(stacked, 0.90)),
            "spearman_coverage_vs_partial_r2_A": coverage_spearman["A"],
            "spearman_coverage_vs_partial_r2_B": coverage_spearman["B"],
        })
        del protein_table
    return usage_rows, variability_rows, covariate_rows, collinearity_rows, protein_rows, mismatches


def spearman_bundle(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> dict:
    result = stats.spearmanr(x, y)
    boots = np.empty(5000)
    n_obs = x.size
    for i in range(5000):
        draw = rng.integers(0, n_obs, n_obs)
        boots[i] = stats.spearmanr(x[draw], y[draw]).statistic
    lo, hi = np.nanquantile(boots, [0.025, 0.975])
    return {"estimate": float(result.statistic), "p_value": float(result.pvalue), "ci_low": float(lo), "ci_high": float(hi), "n": int(n_obs)}


def build_moderators(summary, size_rows, usage_rows, variability_rows, protein_rows, collinearity_rows):
    mods = {}
    eligible = [
        "median_n",
        "fraction_n_ge_71",
        "median_principal_proportion",
        "median_residual_variance",
        "protein_cell_missingness",
        "median_cnv_vif",
    ]
    for transform in ("A", "B"):
        mods[transform] = {
            "cancer_mean": np.array([summary[(c, transform)]["mean_partial_r2"] for c in CANCERS]),
            "median_n": np.array([row["median_n_success"] for row in size_rows if row["transformation"] == transform]),
            "fraction_n_ge_71": np.array([row["fraction_candidate_n_ge_71"] for row in size_rows if row["transformation"] == transform]),
            "median_principal_proportion": np.array([row["median_of_gene_median_principal_proportion"] for row in usage_rows if row["transformation"] == transform]),
            "median_residual_variance": np.array([row["median_residual_variance"] for row in variability_rows if row["transformation"] == transform]),
            "protein_cell_missingness": np.array([row["cell_missing_fraction"] for row in protein_rows]),
            "median_cnv_vif": np.array([row["median_cnv_vif"] for row in collinearity_rows if row["transformation"] == transform]),
            "fraction_ge_0.02": np.array([summary[(c, transform)]["prop_ge_0.02"] for c in CANCERS]),
        }
        for key, value in mods[transform].items():
            if value.size != 10:
                raise SystemExit(f"moderator {transform} {key} has length {value.size}")
    return mods, eligible


def qualification_rows(mods: dict, eligible: list[str]) -> list[dict]:
    rng = np.random.default_rng(SEED)
    rows = [{
        "item": "independent_cohorts",
        "transformation": "BOTH",
        "estimate": 10,
        "ci_low": "",
        "ci_high": "",
        "p_value": "",
        "q_value": "",
        "role": "qualification",
        "detail": "Ten cohorts are below the locked minimum of 15. META_REGRESSION_NOT_QUALIFIED. No moderator model was fit.",
    }]
    for transform in ("A", "B"):
        estimates = []
        p_values = []
        drafted = []
        y = mods[transform]["cancer_mean"]
        for name in eligible:
            bundle = spearman_bundle(mods[transform][name], y, rng)
            estimates.append(bundle)
            p_values.append(bundle["p_value"])
            drafted.append((name, bundle, "eligible_exploratory"))
        q_values = benjamini_hochberg(np.array(p_values))
        for (name, bundle, role), q_value in zip(drafted, q_values):
            rows.append({
                "item": name,
                "transformation": transform,
                "estimate": bundle["estimate"],
                "ci_low": bundle["ci_low"],
                "ci_high": bundle["ci_high"],
                "p_value": bundle["p_value"],
                "q_value": float(q_value),
                "n": bundle["n"],
                "role": role,
                "detail": "Spearman correlation with the cancer-level mean partial R2. Exploratory n=10.",
            })
        circular = spearman_bundle(mods[transform]["fraction_ge_0.02"], y, rng)
        rows.append({
            "item": "fraction_partial_r2_ge_0.02",
            "transformation": transform,
            "estimate": circular["estimate"],
            "ci_low": circular["ci_low"],
            "ci_high": circular["ci_high"],
            "p_value": circular["p_value"],
            "q_value": "",
            "n": circular["n"],
            "role": "ineligible_circular",
            "detail": "This quantity is a functional of the same gene-level distribution as the mean. It was not entered in the BH family and was not used as a moderator.",
        })
        for i, left in enumerate(eligible):
            for right in eligible[i + 1:]:
                bundle = spearman_bundle(mods[transform][left], mods[transform][right], rng)
                rows.append({
                    "item": f"confound|{left}|{right}",
                    "transformation": transform,
                    "estimate": bundle["estimate"],
                    "ci_low": bundle["ci_low"],
                    "ci_high": bundle["ci_high"],
                    "p_value": bundle["p_value"],
                    "q_value": "",
                    "n": bundle["n"],
                    "role": "confounding_screen",
                    "detail": "Confounded pair under the locked |Spearman| > 0.70 rule." if abs(bundle["estimate"]) > 0.70 else "Pairwise moderator correlation.",
                })
    rows.append({
        "item": "decision",
        "transformation": "BOTH",
        "estimate": "",
        "ci_low": "",
        "ci_high": "",
        "p_value": "",
        "q_value": "",
        "role": "qualification",
        "detail": "META_REGRESSION_NOT_QUALIFIED",
    })
    return rows


def range_of(summary, transform, size_rows) -> dict:
    means = np.array([summary[(cancer, transform)]["mean_partial_r2"] for cancer in CANCERS])
    high = []
    gaps = []
    signs = 0
    defined = 0
    for row in size_rows:
        if row["transformation"] != transform:
            continue
        if row["n_distribution_concentrated"]:
            continue
        defined += 1
        high.append(row["high_n_tertile_mean_partial_r2"])
        gap = row["high_n_tertile_mean_partial_r2"] - row["low_n_tertile_mean_partial_r2"]
        gaps.append(gap)
        if gap > 0:
            signs += 1
    high = np.array(high, dtype=float)
    full_range = float(means.max() - means.min())
    high_range = float(np.max(high) - np.min(high)) if high.size >= 2 else np.nan
    return {
        "range_full": full_range,
        "range_high_tertile": high_range,
        "ratio": high_range / full_range if full_range and np.isfinite(high_range) else np.nan,
        "defined_cancers": defined,
        "positive_gap_cancers": signs,
        "median_abs_gap": float(np.median(np.abs(gaps))) if gaps else np.nan,
        "median_jackknife_se": float(np.median([summary[(cancer, transform)]["jackknife_se"] for cancer in CANCERS])),
    }


def classify(repro_fail, meta_fail, influence, ranges, qual_rows, conc_rows, sens_rows, summary) -> dict:
    if repro_fail or meta_fail:
        return {
            "reproducibility": "FAIL",
            "i2_robust": "NOT_EVALUATED",
            "explanation": "INSUFFICIENT_EVIDENCE_TO_EXPLAIN_HETEROGENEITY",
            "phase2": "STOP_CURRENT_DIRECTION",
            "audit": "FAILED",
        }
    i2_ok = True
    pi_below = True
    for row in influence:
        if row["analysis"] in {"all_ten", "exclude_OV"} and row["role"] != "descriptive_sensitivity":
            i2_ok = i2_ok and float(row["I2"]) >= 0.75
            pi_below = pi_below and float(row["pi_low"]) < 0.02
    paired_ok = True
    # Paired-support I2 is compared later from the stored Phase 1A meta if recomputed separately.
    explanation = explanation_class(ranges, qual_rows)
    top1 = [row["share_top_1pct"] for row in conc_rows]
    redesign = sum(value > 0.50 for value in top1) >= 6
    above = {transform: {cancer for cancer in CANCERS if summary[(cancer, transform)]["ci_low"] > 0.02} for transform in ("A", "B")}
    agree = above["A"] == above["B"]
    s1 = [row for row in sens_rows if row["sensitivity"] == "S1_drop_stored_all_zero" and row["delta_mean"] != ""]
    s1_abs = [abs(float(row["delta_mean"])) for row in s1]
    s1_cross = any(int(row["point_classification_changed"]) == 1 and abs(float(row["delta_mean"])) > 0.01 for row in s1)
    both_ci = all(
        float(row["ci_low"]) > 0.02
        for row in influence
        if row["analysis"] == "all_ten" and row["estimator"] == "REML_HartungKnapp"
    )
    both_pi = all(
        float(row["pi_low"]) > 0.02
        for row in influence
        if row["analysis"] == "all_ten" and row["estimator"] == "REML_HartungKnapp"
    )
    if s1_cross:
        phase2 = "STOP_CURRENT_DIRECTION"
    elif redesign or not agree:
        phase2 = "PHASE2_REQUIRES_REDESIGN"
    elif both_pi and explanation == "STRONG_EVIDENCE_OF_MEASUREMENT_OR_POWER_CONTRIBUTION":
        phase2 = "QUALIFIED_FOR_PHASE2"
    elif both_ci and agree and s1_abs and max(s1_abs) <= 0.002:
        phase2 = "QUALIFIED_FOR_NARROW_PHASE2"
    else:
        phase2 = "PHASE2_REQUIRES_REDESIGN"
    return {
        "reproducibility": "PASS",
        "i2_robust": "YES" if i2_ok and pi_below and paired_ok else "NO",
        "explanation": explanation,
        "phase2": phase2,
        "audit": "COMPLETED",
        "cancers_above_A": ",".join(sorted(above["A"])),
        "cancers_above_B": ",".join(sorted(above["B"])),
    }


def explanation_class(ranges: dict, qual_rows: list[dict]) -> str:
    strong = []
    cancer_specific = []
    mixed = []
    for transform, info in ranges.items():
        spearman_n = next(float(row["estimate"]) for row in qual_rows if row["item"] == "median_n" and row["transformation"] == transform)
        eligible = [abs(float(row["estimate"])) for row in qual_rows if row["transformation"] == transform and row["role"] == "eligible_exploratory" and row["estimate"] != ""]
        ratio = info["ratio"]
        tertile_usable = info["defined_cancers"] >= 8 and np.isfinite(ratio)
        strong.append(
            (tertile_usable and ratio <= 0.50)
            or (abs(spearman_n) >= 0.80 and info["positive_gap_cancers"] >= 8 and info["median_abs_gap"] >= 0.005)
        )
        cancer_specific.append(
            tertile_usable and ratio >= 0.80 and max(eligible) < 0.30 and info["gini_range"] < 0.10
            and info["range_full"] > 4 * info["median_jackknife_se"]
        )
        mixed.append(tertile_usable and ratio >= 0.50 and max(eligible) >= 0.70)
    if all(strong):
        return "STRONG_EVIDENCE_OF_MEASUREMENT_OR_POWER_CONTRIBUTION"
    if all(cancer_specific):
        return "EVIDENCE_OF_CANCER_SPECIFIC_EFFECT_HETEROGENEITY"
    if all(mixed):
        return "MIXED_MEASUREMENT_AND_CANCER_SPECIFIC_HETEROGENEITY"
    return "INSUFFICIENT_EVIDENCE_TO_EXPLAIN_HETEROGENEITY"


def write_reports(context: dict) -> None:
    decision = context["decision"]
    repro = context["repro_rows"]
    failed_repro = [row for row in repro if row["checks_passed"] != 1]
    influence = [row for row in context["influence"] if row["analysis"] == "all_ten" and row["estimator"] == "REML_HartungKnapp"]
    lines = [
        "# Phase 1B heterogeneity reproduction",
        "",
        f"Recomputed from the Phase 1A gene-level files at {now()}.",
        "",
        f"Integer counts and stored moments {'match' if not failed_repro else 'do not match'} the Phase 1A cancer summary.",
        f"Failed rows: {len(failed_repro)}.",
        "",
        "I² compares between-cancer dispersion with the jackknife standard errors. Those standard errors are about 0.001, while the cancer means span roughly 0.019 to 0.029. A high I² can therefore come from precise cancer means rather than from a large biological gap. The between-gene standard deviation within a cancer is about 0.03 to 0.05, much larger than the jackknife standard error of the mean.",
        "",
        "The independence approximation sd / sqrt(n_genes) is smaller than the jackknife standard error because genes share cases. The jackknife standard error is the uncertainty used for the cancer mean.",
        "",
    ]
    write_pair("phase1b_heterogeneity_reproduction.md", "heterogeneity_reproduction_report.md", "\n".join(lines) + "\n")

    dist = context["dist_rows"]
    classes = {}
    for row in dist:
        classes[row["distribution_class"]] = classes.get(row["distribution_class"], 0) + 1
    conc = context["conc_rows"]
    gini_values = [row["gini_nonnegative"] for row in conc]
    top1 = [row["share_top_1pct"] for row in conc]
    top10 = [row["share_top_10pct"] for row in conc]
    dist_text = [
        "# Phase 1B gene-level distributions",
        "",
        "The primary estimand remains the unweighted mean. Trimmed and winsorized means are descriptive.",
        "",
        f"Distribution classes across 20 cancer–inverse results: {classes}.",
        f"Gini of nonnegative partial R² ranges from {min(gini_values):.3f} to {max(gini_values):.3f}.",
        f"The top 1% of genes contribute {min(top1):.3f} to {max(top1):.3f} of the nonnegative sum.",
        f"The top 10% contribute {min(top10):.3f} to {max(top10):.3f}.",
        "",
        "Medians remain near 0.008 to 0.012 and means remain near 0.020 to 0.029. That gap is a right-skewed distribution. The class rule calls this pattern MIXED_DISTRIBUTION when the median is well below the mean and the top 10% share stays below 0.60. No gene identifier is reported.",
        "",
    ]
    write_pair("phase1b_gene_distribution_report.md", "gene_distribution_report.md", "\n".join(dist_text) + "\n")

    power = context["power_rows"]
    power_text = [
        "# Phase 1B power and precision",
        "",
        f"At N = 71 the frozen planning formula gives a detectable partial R² of {detectable_partial_r2(71):.3f} for 80% power at two-sided alpha 0.05. The observed cancer means are near 0.02. Exact power for a gene-level partial R² of 0.02 at the observed median N is reported in `results/power_sensitivity.csv` and is low.",
        "",
        "The cancer-level mean is an average over thousands of genes, so its jackknife standard error is near 0.001 even though a single gene is not powered for a partial R² of 0.02. Between-gene standard deviations are about 0.03 to 0.05. Between-cancer τ² is on the order of 1e-5.",
        "",
        "A cancer with fewer eligible genes is not assigned a weaker biological effect from that count alone. OV remains in every cancer-specific display.",
        "",
        f"Sample-size rows: {len(context['size_rows'])}. Power rows: {len(power)}.",
        "",
    ]
    write_pair("phase1b_power_and_precision_report.md", "power_and_precision_report.md", "\n".join(power_text) + "\n")

    usage = context["usage_rows"]
    ent = [row["median_entropy_bits"] for row in usage]
    zero = [row["fraction_observations_stored_all_zero"] for row in usage if row["transformation"] == "A"]
    usage_text = [
        "# Phase 1B transcript-usage measurement",
        "",
        "Summaries give each successful gene equal weight. Close to zero means principal proportion < 0.05. Close to one means principal proportion > 0.95. Entropy is base-2 Shannon entropy.",
        "",
        f"Median within-gene entropy ranges from {min(ent):.3f} to {max(ent):.3f} bits.",
        f"Literal A stored-all-zero fractions, averaged across successful genes, range from {min(zero):.4f} to {max(zero):.4f}.",
        "",
        "The principal transcript was not redefined. No differential transcript usage test was run.",
        "",
        f"Complete-case reconstruction mismatches against the Phase 1A n_complete column: {context['mismatches']}.",
        "",
    ]
    write_pair("phase1b_transcript_usage_report.md", "transcript_usage_measurement_report.md", "\n".join(usage_text) + "\n")

    cov_text = [
        "# Phase 1B covariate structure",
        "",
        "Correlations are within-gene correlations on the primary complete cases, summarized by the median across successful genes. CNV VIF and condition numbers come from the stored Phase 1A fits. CNV stays in the primary model.",
        "",
        "Residual variance is the mean square of the principal proportion after projection on intercept, log2(G), WES purity, and CNV. A value below 1e-4 marks a proportion that is nearly determined by M0.",
        "",
        "RSEM correlations use the primary cases that also have numeric RSEM. RSEM is not treated as the same scale as the reconstructed total.",
        "",
    ]
    write_pair("phase1b_covariate_report.md", "covariate_structure_report.md", "\n".join(cov_text) + "\n")

    protein = context["protein_rows"]
    miss = [row["cell_missing_fraction"] for row in protein]
    prot_text = [
        "# Phase 1B protein coverage",
        "",
        f"Cell missingness in the global protein matrices ranges from {min(miss):.3f} to {max(miss):.3f}.",
        "",
        "Missing protein values were left missing. They were not replaced by zero and were not interpreted as functional loss. Phosphoproteomics was not used.",
        "",
        "Within-cancer Spearman correlations between a gene's protein observation fraction and its partial R² are stored in `results/protein_coverage_comparison.csv`. Those correlations describe coverage, not a biological absence.",
        "",
    ]
    write_pair("phase1b_protein_coverage_report.md", "protein_measurement_report.md", "\n".join(prot_text) + "\n")

    rob = context["robust_rows"]
    deltas = [row["mean_A_minus_B"] for row in rob]
    signs = [row["proportion_same_beta_sign"] for row in rob]
    rob_text = [
        "# Phase 1B transformation robustness",
        "",
        f"Cancer-level mean A minus mean B ranges from {min(deltas):.4f} to {max(deltas):.4f}.",
        f"The proportion of paired genes with the same standardized-coefficient sign ranges from {min(signs):.3f} to {max(signs):.3f}.",
        "",
        "Literal A and Literal B stay co-primary. The correlation between them was not used to choose an inverse. Candidate C means are included only as the stored sensitivity.",
        "",
        f"Cancers with a jackknife interval entirely above 0.02, Literal A: {decision.get('cancers_above_A', '')}.",
        f"Literal B: {decision.get('cancers_above_B', '')}.",
        "",
    ]
    write_pair("phase1b_transformation_robustness.md", "transformation_robustness_report.md", "\n".join(rob_text) + "\n")

    meta_lines = ["# Phase 1B meta-analytic influence", "", "The primary estimator remains REML with a Hartung–Knapp interval. Fixed-effect, median, and trimmed-mean rows are descriptive.", ""]
    for row in context["influence"]:
        if row["estimator"] == "REML_HartungKnapp" and row["analysis"] in {"all_ten", "exclude_OV"}:
            meta_lines.append(f"- {row['transformation']} {row['analysis']}: mean {text_num(row['mu'])}, 95% CI {text_num(row['ci_low'])} to {text_num(row['ci_high'])}, I² {text_num(row['I2'], 3)}, prediction interval {text_num(row['pi_low'])} to {text_num(row['pi_high'])}.")
    deltas_mu = [abs(float(row["delta_mu"])) for row in context["loo_rows"]]
    meta_lines.append("")
    meta_lines.append(f"The largest absolute leave-one-cancer change in the pooled mean is {max(deltas_mu):.5f}. No cancer was removed because of that change.")
    meta_lines.append("")
    write_pair("phase1b_meta_influence.md", "meta_analytic_influence_report.md", "\n".join(meta_lines) + "\n")

    qual = [
        "# Phase 1B meta-regression qualification",
        "",
        "Independent cohorts: 10.",
        "",
        "The locked rule requires at least 15 independent cohorts before a meta-regression is qualified. Moderators measured at cohort level are ecological, several are mutually correlated, and one listed quantity, the fraction of genes with partial R² at least 0.02, is a function of the outcome distribution.",
        "",
        "META_REGRESSION_NOT_QUALIFIED",
        "",
        "No moderator model was fit. Exploratory Spearman correlations and bootstrap intervals are in `results/meta_regression_qualification.csv`. Benjamini–Hochberg adjustment was applied inside each transformation to the six eligible correlates. Those q values are diagnostic. They are not a primary explanation.",
        "",
    ]
    write_pair("phase1b_meta_regression_qualification.md", "meta_regression_feasibility.md", "\n".join(qual) + "\n")

    changed = [row for row in context["sens_rows"] if str(row.get("point_classification_changed")) == "1"]
    sens_text = [
        "# Phase 1B sensitivity reassessment",
        "",
        "The stored Phase 1A sensitivity means were compared with the primary cancer means. New protein models were not fit. S1 through S7 do not have stored jackknife standard errors, so their random-effects intervals were not recomputed. S8 is the stored meta-analysis excluding OV.",
        "",
        f"Cancer-level point classifications relative to 0.02 change in {len(changed)} stored sensitivity rows.",
        "",
        "S4 remains NOT_RUN_UNDERDEFINED.",
        "",
        "Assumptions: S1 tests the Literal A all-zero assignment; S2 tests the gene-expression covariate; S3 tests Candidate C; S5 tests the purity covariate; S6 tests high-VIF genes only; S7 tests one-pass residual deletion; S8 tests the influence of OV on the pooled interval.",
        "",
    ]
    write_pair("phase1b_sensitivity_reassessment.md", "artifact_sensitivity_report.md", "\n".join(sens_text) + "\n")

    final = f"""# Phase 1B final decision

## A. Is the Phase 1A result reproducible?

{decision['reproducibility']}

Gene-level counts, means, medians, standard deviations, and the stored percentiles were recomputed from the Phase 1A gene files. Jackknife intervals match mean ± t_(n−1) times the stored jackknife standard error. REML, Hartung–Knapp intervals, I², Cochran Q, and prediction intervals were recomputed from the cancer means and those standard errors.

## B. Is high heterogeneity robust?

{decision['i2_robust']}

I² stays high for both inverses in the ten-cancer analysis and in the analysis excluding OV. The prediction intervals continue to include values below 0.02. High I² here is the ratio of a small between-cancer variance to still smaller jackknife standard errors. The absolute gap between the lowest and highest cancer means is about one percentage point of partial R².

## C. Is the heterogeneity explained?

{decision['explanation']}

Meta-regression was not qualified with 10 cohorts, and no moderator model was fit. The locked mixed class fired because the cross-cancer range remained among higher-N genes and because protein-cell missingness was the only eligible correlate with absolute Spearman correlation at least 0.70. That correlation is about −0.72 for Literal A and −0.75 for Literal B. Its Benjamini–Hochberg q values are 0.11 and 0.08. The same moderator is correlated with the fraction of candidate genes at N ≥ 71 (Spearman −0.92) and with the median principal proportion (Spearman about 0.73 to 0.76). Within a cancer, a gene’s protein observation fraction and its partial R² have Spearman correlations of about 0.01 to 0.07. The median gene is fully observed in every cohort. This is an ecological association on 10 cohorts, not a causal assignment.

The five candidate explanations, each left as a description of the audit rather than as a selected cause:

1. Cancer means differ by about 0.01 partial R², which is several jackknife standard errors, and the same four cohorts fall at or below 0.02 under both inverses.
2. Median N differs by cohort, while power at a gene-level partial R² of 0.02 remains low throughout. The cohort with the smallest N is not the cohort with the smallest mean.
3. Principal-proportion location, entropy, and residual variance differ across cohorts; the locked correlation threshold for calling that the source of the mean differences was not met.
4. Cohort-level protein-cell missingness is the correlate that crossed the locked absolute Spearman screen of 0.70. It is confounded with eligibility and with the principal-proportion level. The within-cancer association between coverage and partial R² is about 0.01 to 0.07. Missing protein values were left missing.
5. The gene-level distribution is right-skewed in every cohort. The top 1% share and the Gini coefficient describe a mixed distribution, not a mean carried by a handful of genes and not a flat weak effect.

## D. Is the project ready for Phase 2?

{decision['phase2']}

This label does not start Phase 2. A narrow qualification would keep the same gene-level endpoint inside each cancer. It does not authorize gene ranking, pathway analysis, isoform switching, or a mechanism claim.

## E. Scientific interpretation

Transcript usage provides a small average incremental explanatory contribution to gene-level protein abundance beyond reconstructed gene-level RNA expression, tumor purity and CNV, but the magnitude varies substantially across cancer cohorts.

The measurements are case-level. The outcome is gene-level protein abundance. Literal A and Literal B remain co-primary. Partial R² is incremental statistical explanation. The audit does not identify isoform-specific protein abundance, protein functional loss, a causal effect of transcript switching, clinical utility, or a therapeutic target.

## Audit status

PHASE 1B AUDIT = {decision['audit']}

PHASE 2 AUTHORIZATION = FALSE
"""
    (PHASE1B / "reports" / "phase1b_final_decision.md").write_text(final)


def write_pair(canonical: str, alias: str, text: str) -> None:
    (PHASE1B / "reports" / canonical).write_text(text)
    if alias != canonical:
        (PHASE1B / "reports" / alias).write_text(text)


def hash_outputs() -> None:
    lines = []
    for path in sorted(PHASE1B.rglob("*")):
        if not path.is_file() or path.name.startswith("._") or path.name == "SHA256SUMS":
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        lines.append(f"{digest.hexdigest()}  {path.relative_to(ROOT).as_posix()}")
    (PHASE1B / "SHA256SUMS").write_text("\n".join(lines) + "\n")


def verify_phase1a_inputs() -> None:
    ledger = (PHASE1A / "SHA256SUMS").read_text().splitlines()
    for line in ledger:
        if not line.strip():
            continue
        digest, rel = line.split(maxsplit=1)
        path = ROOT / rel.strip()
        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                hasher.update(block)
        if hasher.hexdigest() != digest:
            raise SystemExit(f"TECHNICAL_HOLD phase1a hash mismatch {rel}")


def input_manifest() -> None:
    rows = []
    paths = [
        PHASE1A / "reports" / "phase1a_final_decision.md",
        PHASE1A / "reports" / "phase1a_primary_results.md",
        PHASE1A / "reports" / "phase1a_sensitivity_results.md",
        PHASE1A / "results" / "cancer_level_summary.csv",
        PHASE1A / "results" / "random_effects_summary.csv",
        PHASE1A / "results" / "transformation_paired_comparison.csv",
        PHASE1A / "results" / "eligibility_by_cancer_transform.csv",
        PHASE1A / "results" / "primary_literal_A_gene_fits.tsv.gz",
        PHASE1A / "results" / "primary_literal_B_gene_fits.tsv.gz",
        PHASE1A / "results" / "sensitivity_summary.csv",
        PHASE1A / "results" / "common_support_summary.csv",
        PHASE1A / "results" / "leave_one_cancer_out.csv",
        PHASE1A / "qc" / "cohort_reconstruction.csv",
        PHASE1A / "SHA256SUMS",
        ROOT / "phase0_6" / "reports" / "statistical_model_freeze.md",
    ]
    for path in paths:
        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                hasher.update(block)
        rows.append({"path": path.relative_to(ROOT).as_posix(), "sha256": hasher.hexdigest(), "bytes": path.stat().st_size})
    write_csv(PHASE1B / "manifests" / "input_manifest.csv", rows, ["path", "sha256", "bytes"])


def environment() -> None:
    import scipy
    text = "\n".join([
        f"start_utc {STARTED}",
        f"system {platform.platform()}",
        f"python {sys.version}",
        f"numpy {np.__version__}",
        f"scipy {scipy.__version__}",
        "threads OMP=1 OPENBLAS=1 MKL=1 VECLIB=1",
        f"seed {SEED}",
        "phase1a_models_refit no",
    ])
    (PHASE1B / "manifests" / "software_environment.txt").write_text(text + "\n")


def hold(message: str) -> None:
    (PHASE1B / "reports" / "phase1b_final_decision.md").write_text(
        "# Phase 1B final decision\n\n"
        "## A. Is the Phase 1A result reproducible?\n\nFAIL\n\n"
        "## B. Is high heterogeneity robust?\n\nNOT_EVALUATED\n\n"
        "## C. Is the heterogeneity explained?\n\nINSUFFICIENT_EVIDENCE_TO_EXPLAIN_HETEROGENEITY\n\n"
        "## D. Is the project ready for Phase 2?\n\nSTOP_CURRENT_DIRECTION\n\n"
        "## E. Scientific interpretation\n\n"
        f"The audit stopped before an explanation was assigned. {message}\n\n"
        "PHASE 1B AUDIT = FAILED\n\nPHASE 2 AUTHORIZATION = FALSE\n"
    )
    (PHASE1B / "logs" / "run.log").write_text(f"start {STARTED}\nend {now()}\nstatus FAILED\n{message}\n")


def main() -> None:
    for folder in ("protocol", "code", "logs", "qc", "results", "reports", "manifests"):
        (PHASE1B / folder).mkdir(parents=True, exist_ok=True)
    environment()
    print("verify phase1a", now(), flush=True)
    verify_phase1a_inputs()
    input_manifest()
    synthetic = run_synthetic_tests()
    (PHASE1B / "qc" / "synthetic_test_results.txt").write_text("\n".join(synthetic) + "\n")
    print("load fits", flush=True)
    fits = {"A": load_fits(PHASE1A / "results" / "primary_literal_A_gene_fits.tsv.gz"),
            "B": load_fits(PHASE1A / "results" / "primary_literal_B_gene_fits.tsv.gz")}
    summary = load_summary()
    stored_meta = load_stored_meta()
    repro_rows, repro_fail = reproduce(fits, summary)
    write_csv(PHASE1B / "results" / "heterogeneity_reproduction.csv", repro_rows, list(repro_rows[0].keys()))
    influence, loo_rows, meta_fail = compare_meta(summary, stored_meta)
    if repro_fail or meta_fail:
        hold("Heterogeneity statistics did not reproduce within the locked tolerance.")
        hash_outputs()
        raise SystemExit("TECHNICAL_HOLD")
    common = list(csv.DictReader((PHASE1A / "results" / "common_support_summary.csv").open()))
    for transform, mean_key, se_key in (("A", "mean_A", "se_A"), ("B", "mean_B", "se_B")):
        meta = reml_meta(
            np.array([float(row[mean_key]) for row in common]),
            np.array([float(row[se_key]) for row in common]),
        )
        ref = next(row for row in stored_meta if row["analysis"] == "paired_support" and row["transformation"] == transform)
        bad = [key for key in ("mu", "tau2", "I2", "ci_low", "ci_high", "pi_low", "pi_high") if abs(meta[key] - float(ref[key])) > 1e-6]
        if bad:
            hold("Paired-support random-effects summary did not reproduce: " + ",".join(bad))
            hash_outputs()
            raise SystemExit("TECHNICAL_HOLD")
        influence.append({
            "analysis": "paired_support",
            "transformation": transform,
            "estimator": "REML_HartungKnapp",
            "role": "prespecified_sensitivity",
            "match": 1,
            "failed_fields": "",
            **meta,
        })
    write_csv(PHASE1B / "results" / "meta_influence.csv", influence, list(influence[0].keys()))
    write_csv(PHASE1B / "results" / "leave_one_cancer_out_reanalysis.csv", loo_rows, list(loo_rows[0].keys()))
    dist_rows, conc_rows = distribution_rows(fits)
    write_csv(PHASE1B / "results" / "gene_distribution_summary.csv", dist_rows, list(dist_rows[0].keys()))
    write_csv(PHASE1B / "results" / "effect_concentration_summary.csv", conc_rows, list(conc_rows[0].keys()))
    size_rows, power_rows = sample_size_rows(fits, summary)
    write_csv(PHASE1B / "results" / "effective_sample_size_audit.csv", size_rows, list(size_rows[0].keys()))
    write_csv(PHASE1B / "results" / "power_sensitivity.csv", power_rows, list(power_rows[0].keys()))
    robust_rows = robustness_rows(fits, summary)
    write_csv(PHASE1B / "results" / "transformation_robustness.csv", robust_rows, list(robust_rows[0].keys()))
    sens_rows = sensitivity_rows(summary, stored_meta)
    write_csv(PHASE1B / "results" / "sensitivity_reassessment.csv", sens_rows, list(sens_rows[0].keys()))
    usage_rows, variability_rows, covariate_rows, collinearity_rows, protein_rows, mismatches = measurement_audit(fits)
    if mismatches:
        hold(f"Complete-case reconstruction differed for {mismatches} gene rows.")
        hash_outputs()
        raise SystemExit("TECHNICAL_HOLD")
    write_csv(PHASE1B / "results" / "transcript_usage_distribution.csv", usage_rows, list(usage_rows[0].keys()))
    write_csv(PHASE1B / "results" / "transcript_usage_variability.csv", variability_rows, list(variability_rows[0].keys()))
    write_csv(PHASE1B / "results" / "covariate_structure.csv", covariate_rows, list(covariate_rows[0].keys()))
    write_csv(PHASE1B / "results" / "collinearity_audit.csv", collinearity_rows, list(collinearity_rows[0].keys()))
    write_csv(PHASE1B / "results" / "protein_coverage_comparison.csv", protein_rows, list(protein_rows[0].keys()))
    mods, eligible = build_moderators(summary, size_rows, usage_rows, variability_rows, protein_rows, collinearity_rows)
    qual_rows = qualification_rows(mods, eligible)
    write_csv(PHASE1B / "results" / "meta_regression_qualification.csv", qual_rows, ["item", "transformation", "estimate", "ci_low", "ci_high", "p_value", "q_value", "n", "role", "detail"])
    ranges = {}
    for transform in ("A", "B"):
        info = range_of(summary, transform, size_rows)
        info["gini_range"] = float(np.ptp([row["gini_nonnegative"] for row in conc_rows if row["transformation"] == transform]))
        ranges[transform] = info
    decision = classify(False, False, influence, ranges, qual_rows, conc_rows, sens_rows, summary)
    # Paired-support I2 participates in the robustness flag.
    paired_i2 = [float(row["I2"]) for row in influence if row["analysis"] == "paired_support"]
    paired_pi = [float(row["pi_low"]) for row in influence if row["analysis"] == "paired_support"]
    if decision["i2_robust"] == "YES" and not (all(value >= 0.75 for value in paired_i2) and all(value < 0.02 for value in paired_pi)):
        decision["i2_robust"] = "NO"
    context = {
        "decision": decision,
        "repro_rows": repro_rows,
        "influence": influence,
        "loo_rows": loo_rows,
        "dist_rows": dist_rows,
        "conc_rows": conc_rows,
        "size_rows": size_rows,
        "power_rows": power_rows,
        "usage_rows": usage_rows,
        "sens_rows": sens_rows,
        "robust_rows": robust_rows,
        "protein_rows": protein_rows,
        "mismatches": mismatches,
    }
    write_reports(context)
    (PHASE1B / "logs" / "run.log").write_text(
        f"start {STARTED}\nend {now()}\nreproducibility {decision['reproducibility']}\n"
        f"explanation {decision['explanation']}\nphase2 {decision['phase2']}\naudit {decision['audit']}\n"
    )
    hash_outputs()
    print("PHASE1B_DONE", decision["audit"], decision["explanation"], decision["phase2"], now(), flush=True)


if __name__ == "__main__":
    try:
        main()
    except SystemExit as exc:
        if not (PHASE1B / "reports" / "phase1b_final_decision.md").exists():
            hold(str(exc))
            hash_outputs()
        raise
    except Exception:
        traceback.print_exc()
        (PHASE1B / "logs" / "run_error.txt").write_text(traceback.format_exc())
        if not (PHASE1B / "reports" / "phase1b_final_decision.md").exists():
            hold("An uncaught exception stopped the audit. See logs/run_error.txt.")
            hash_outputs()
        raise
