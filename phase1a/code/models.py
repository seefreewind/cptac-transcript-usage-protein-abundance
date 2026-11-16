"""OLS, jackknife, and random-effects helpers. Protein files are not read here."""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar

EPS = np.finfo(float).eps
CROSSCHECK_TOL = 1e-8
CLIP_TOL = 1e-8


class NumericalFailure(RuntimeError):
    pass


def _svd_fit(x_matrix: np.ndarray, y: np.ndarray) -> dict | None:
    u_mat, singular, vt_mat = np.linalg.svd(x_matrix, full_matrices=False)
    if singular.size == 0 or singular[0] == 0:
        return {"ok": False, "rank": 0, "cond": np.inf}
    tolerance = max(x_matrix.shape) * EPS * singular[0]
    rank = int(np.sum(singular > tolerance))
    if rank < x_matrix.shape[1]:
        cond = np.inf if rank == 0 else float(singular[0] / singular[rank - 1])
        return {"ok": False, "rank": rank, "cond": cond}
    inv_s = 1.0 / singular
    beta = vt_mat.T @ (inv_s * (u_mat.T @ y))
    resid = y - x_matrix @ beta
    rss = float(resid @ resid)
    cond = float(singular[0] / singular[-1])
    xtx_inv_diag = np.sum((vt_mat.T ** 2) * (inv_s ** 2), axis=1)
    return {
        "ok": True,
        "beta": beta,
        "resid": resid,
        "rss": rss,
        "rank": rank,
        "cond": cond,
        "xtx_inv_diag": xtx_inv_diag,
    }


def _design(columns: list[np.ndarray]) -> np.ndarray | None:
    raw = np.column_stack(columns)
    sd = raw.std(axis=0, ddof=0)
    scale = np.maximum(1.0, np.max(np.abs(raw), axis=0))
    if np.any(sd <= 1e-8 * scale) or np.any(~np.isfinite(sd)):
        return None
    scaled = (raw - raw.mean(axis=0)) / sd
    return np.column_stack([np.ones(raw.shape[0]), scaled])


def fit_nested(log_g: np.ndarray, proportion: np.ndarray, purity: np.ndarray,
               cnv: np.ndarray, protein: np.ndarray, clips: list | None = None) -> dict:
    n_obs = protein.shape[0]
    if np.nanstd(protein) == 0 or np.unique(protein).size < 2:
        return {"fit_status": "ZERO_PROTEIN_VARIANCE", "n_complete": n_obs}
    x0 = _design([log_g, purity, cnv])
    x1 = _design([log_g, proportion, purity, cnv])
    if x0 is None or x1 is None:
        return {"fit_status": "NO_RESIDUAL_P_VARIATION", "n_complete": n_obs}
    fit0 = _svd_fit(x0, protein)
    fit1 = _svd_fit(x1, protein)
    if not fit0["ok"] or not fit1["ok"]:
        return {
            "fit_status": "RANK_DEFICIENT",
            "n_complete": n_obs,
            "rank_m0": fit0["rank"],
            "rank_m1": fit1["rank"],
            "condition_number": fit1["cond"],
        }
    p_fit = _svd_fit(x0, proportion)
    if not p_fit["ok"] or p_fit["rss"] <= 1e-18 * max(1.0, float(proportion @ proportion)):
        return {"fit_status": "NO_RESIDUAL_P_VARIATION", "n_complete": n_obs, "rank_m0": fit0["rank"], "rank_m1": fit1["rank"]}
    rss0 = fit0["rss"]
    rss1 = fit1["rss"]
    gap = rss0 - rss1
    if gap < 0:
        scale = max(rss0, 1.0)
        if gap < -CLIP_TOL * scale:
            raise NumericalFailure(f"RSS_M1 exceeds RSS_M0 by {-gap}")
        if clips is not None:
            clips.append(gap)
        gap = 0.0
    if rss0 <= 0:
        return {"fit_status": "M0_RSS_ZERO", "n_complete": n_obs, "rank_m0": fit0["rank"], "rank_m1": fit1["rank"]}
    partial = gap / rss0
    resid_r = float(np.corrcoef(fit0["resid"], p_fit["resid"])[0, 1])
    if not np.isfinite(resid_r):
        raise NumericalFailure("partial correlation is not finite")
    if abs(resid_r ** 2 - partial) > CROSSCHECK_TOL:
        raise NumericalFailure(f"partial R2 cross-check failed: {resid_r ** 2} vs {partial}")
    df1 = n_obs - fit1["rank"]
    df0 = n_obs - fit0["rank"]
    if df1 <= 0 or df0 <= 0:
        return {"fit_status": "RANK_DEFICIENT", "n_complete": n_obs, "rank_m0": fit0["rank"], "rank_m1": fit1["rank"]}
    if rss1 == 0:
        f_stat = np.inf
        p_value = 0.0
    else:
        f_stat = (gap / 1.0) / (rss1 / df1)
        p_value = float(stats.f.sf(f_stat, 1, df1))
    sigma2 = rss1 / df1
    se_beta = float(np.sqrt(sigma2 * fit1["xtx_inv_diag"][2]))
    tss = float(np.sum((protein - protein.mean()) ** 2))
    r2_0 = 1.0 - rss0 / tss
    r2_1 = 1.0 - rss1 / tss
    adj0 = 1.0 - (1.0 - r2_0) * (n_obs - 1) / df0
    adj1 = 1.0 - (1.0 - r2_1) * (n_obs - 1) / df1
    cnv_col = x1[:, 4]
    others = np.column_stack([x1[:, 0], x1[:, 1], x1[:, 2], x1[:, 3]])
    cnv_fit = _svd_fit(others, cnv_col)
    if not cnv_fit["ok"] or cnv_fit["rss"] <= 0:
        vif = np.inf
    else:
        tss_cnv = float(np.sum((cnv_col - cnv_col.mean()) ** 2))
        r2_cnv = 1.0 - cnv_fit["rss"] / tss_cnv if tss_cnv > 0 else 1.0
        vif = np.inf if r2_cnv >= 1 else 1.0 / (1.0 - r2_cnv)
    expected = 1.0 / (df1 + 1.0)
    pearson = fit0["resid"] / np.sqrt(rss0 / df0)
    return {
        "fit_status": "OK",
        "n_complete": n_obs,
        "rank_m0": fit0["rank"],
        "rank_m1": fit1["rank"],
        "df_residual_m1": df1,
        "condition_number": fit1["cond"],
        "cnv_vif": float(vif),
        "r2_m0": r2_0,
        "r2_m1": r2_1,
        "adjusted_r2_m0": adj0,
        "adjusted_r2_m1": adj1,
        "partial_r2": partial,
        "expected_null_partial_r2": expected,
        "null_excess_partial_r2": partial - expected,
        "signed_partial_r": resid_r,
        "beta_p_standardized": float(fit1["beta"][2]),
        "se_beta_p": se_beta,
        "f_statistic": float(f_stat),
        "p_value": p_value,
        "pearson_m0": pearson,
    }


def partial_r2_fast(log_g, proportion, purity, cnv, protein) -> float | None:
    x0 = _design([log_g, purity, cnv])
    x1 = _design([log_g, proportion, purity, cnv])
    if x0 is None or x1 is None:
        return None
    fit0 = _svd_fit(x0, protein)
    fit1 = _svd_fit(x1, protein)
    if not fit0["ok"] or not fit1["ok"] or fit0["rss"] <= 0:
        return None
    gap = fit0["rss"] - fit1["rss"]
    if gap < -CLIP_TOL * max(fit0["rss"], 1.0):
        raise NumericalFailure("jackknife RSS inversion")
    if gap < 0:
        gap = 0.0
    return gap / fit0["rss"]


def partial_r2_only(log_g, proportion, purity, cnv, protein, clips: list | None = None) -> float | None:
    result = fit_nested(log_g, proportion, purity, cnv, protein, clips)
    if result["fit_status"] != "OK":
        return None
    return float(result["partial_r2"])


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    q_values = np.full(p_values.shape, np.nan)
    finite = np.isfinite(p_values)
    if finite.sum() == 0:
        return q_values
    idx = np.flatnonzero(finite)
    ordered = np.argsort(p_values[idx], kind="mergesort")
    ranked = p_values[idx][ordered]
    m_tests = ranked.size
    raw = ranked * m_tests / np.arange(1, m_tests + 1)
    raw = np.minimum.accumulate(raw[::-1])[::-1]
    raw = np.clip(raw, 0, 1)
    placed = np.empty(m_tests)
    placed[ordered] = raw
    q_values[idx] = placed
    return q_values


def jackknife_mean(partial: np.ndarray, membership: list[np.ndarray], refit) -> dict:
    """refit(gene_index, dropped_case) -> float or None."""
    ok = np.isfinite(partial)
    genes = np.flatnonzero(ok)
    base = partial[genes]
    full_mean = float(base.mean()) if base.size else np.nan
    users: list[list[int]] = []
    n_cases = 0
    if membership:
        n_cases = max((int(idx.max()) if idx.size else -1) for idx in membership) + 1
        users = [[] for _ in range(n_cases)]
        for local, gene in enumerate(genes):
            for case in membership[gene]:
                users[int(case)].append(local)
    active = [i for i, group in enumerate(users) if group]
    n_jack = len(active)
    thetas = np.empty(n_jack)
    success = np.empty(n_jack, dtype=int)
    failures = np.empty(n_jack, dtype=int)
    fail_by_gene = np.zeros(genes.size, dtype=int)
    touch_by_gene = np.zeros(genes.size, dtype=int)
    total = float(base.sum())
    for slot, case in enumerate(active):
        if slot % 25 == 0:
            print(f"    jackknife replicate {slot + 1}/{n_jack}", flush=True)
        delta = 0.0
        n_ok = genes.size
        fails = 0
        for local in users[case]:
            touch_by_gene[local] += 1
            updated = refit(int(genes[local]), case)
            if updated is None or not np.isfinite(updated):
                fails += 1
                fail_by_gene[local] += 1
                n_ok -= 1
                delta -= base[local]
            else:
                delta += float(updated) - base[local]
        success[slot] = n_ok
        failures[slot] = fails
        thetas[slot] = (total + delta) / n_ok if n_ok else np.nan
    if n_jack < 2 or np.any(~np.isfinite(thetas)):
        se = np.nan
        ci = (np.nan, np.nan)
    else:
        center = float(thetas.mean())
        se = float(np.sqrt((n_jack - 1) / n_jack * np.sum((thetas - center) ** 2)))
        critical = float(stats.t.ppf(0.975, n_jack - 1))
        ci = (full_mean - critical * se, full_mean + critical * se)
    touched = touch_by_gene > 0
    always_ok = int(np.sum(fail_by_gene == 0))
    unstable_genes = 0
    if n_jack:
        unstable_genes = int(np.sum((fail_by_gene / np.maximum(touch_by_gene, 1) >= 0.05) & touched))
    unstable = bool(genes.size and unstable_genes / genes.size > 0.01)
    return {
        "full_mean": full_mean,
        "n_jackknife": n_jack,
        "se": se,
        "ci_low": ci[0],
        "ci_high": ci[1],
        "min_success": int(success.min()) if n_jack else 0,
        "max_success": int(success.max()) if n_jack else 0,
        "genes_always_ok": always_ok,
        "unstable": unstable,
        "max_failures": int(failures.max()) if n_jack else 0,
        "sum_failures": int(failures.sum()) if n_jack else 0,
    }


def reml_meta(theta: np.ndarray, se: np.ndarray) -> dict:
    theta = np.asarray(theta, dtype=float)
    se = np.asarray(se, dtype=float)
    k = theta.size
    se2 = se ** 2

    def nll(tau2: float) -> float:
        variance = se2 + tau2
        if np.any(variance <= 0):
            return 1e300
        weight = 1.0 / variance
        sw = float(weight.sum())
        mu = float(np.sum(weight * theta) / sw)
        return 0.5 * float(np.sum(np.log(variance))) + 0.5 * np.log(sw) + 0.5 * float(np.sum(weight * (theta - mu) ** 2))

    upper = max(1e-8, float(np.var(theta)) * 100 + float(np.max(se2)))
    opt = minimize_scalar(nll, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-14})
    tau2 = 0.0 if nll(0.0) <= nll(float(opt.x)) + 1e-10 else float(opt.x)
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
    q_p = float(stats.chi2.sf(cochran_q, k - 1))
    pred_se = float(np.sqrt(tau2 + se_hk ** 2))
    pred_t = float(stats.t.ppf(0.975, k - 2))
    return {
        "k": k,
        "mu": mu,
        "se_hk": se_hk,
        "ci_low": mu - critical * se_hk,
        "ci_high": mu + critical * se_hk,
        "tau2": tau2,
        "I2": i2,
        "Q": cochran_q,
        "Q_p": q_p,
        "pi_low": mu - pred_t * pred_se,
        "pi_high": mu + pred_t * pred_se,
    }


def run_synthetic_tests() -> list[str]:
    lines = []
    rng = np.random.default_rng(20261114)
    n_obs = 80
    log_g = rng.normal(size=n_obs)
    purity = rng.normal(size=n_obs)
    cnv = rng.normal(size=n_obs)
    proportion = rng.normal(size=n_obs)
    protein = 1.2 + 0.4 * log_g - 0.2 * purity + 0.1 * cnv + 0.5 * proportion
    positive = fit_nested(log_g, proportion, purity, cnv, protein)
    negative = fit_nested(log_g, proportion, purity, cnv, -protein)
    lines.append(f"positive_beta {positive['beta_p_standardized']:.6g} partial {positive['partial_r2']:.6g}")
    lines.append(f"negative_beta {negative['beta_p_standardized']:.6g} sign {negative['signed_partial_r']:.6g}")
    if positive["fit_status"] != "OK" or positive["beta_p_standardized"] <= 0 or positive["signed_partial_r"] <= 0:
        raise NumericalFailure("synthetic positive coefficient failed")
    if negative["beta_p_standardized"] >= 0 or negative["signed_partial_r"] >= 0:
        raise NumericalFailure("synthetic negative coefficient failed")
    if abs(positive["partial_r2"] - negative["partial_r2"]) > 1e-10:
        raise NumericalFailure("sign flip changed partial R2")
    noise = 0.05 * rng.normal(size=n_obs)
    protein0 = 0.3 + 0.8 * log_g + 0.2 * purity - 0.1 * cnv + noise
    zero_add = fit_nested(log_g, proportion, purity, cnv, protein0)
    lines.append(f"near_zero_increment partial {zero_add['partial_r2']:.6g} expected {zero_add['expected_null_partial_r2']:.6g}")
    if zero_add["partial_r2"] > 0.2:
        raise NumericalFailure("orthogonal increment is too large for this draw")
    collinear = fit_nested(log_g, log_g.copy(), purity, cnv, protein0)
    lines.append(f"collinear_status {collinear['fit_status']}")
    if collinear["fit_status"] not in {"RANK_DEFICIENT", "NO_RESIDUAL_P_VARIATION"}:
        raise NumericalFailure("collinear P was accepted")
    flat = fit_nested(log_g, np.full(n_obs, 0.4), purity, cnv, protein0)
    lines.append(f"flat_P_status {flat['fit_status']}")
    if flat["fit_status"] != "NO_RESIDUAL_P_VARIATION":
        raise NumericalFailure("constant P was accepted")
    missing = protein.copy()
    missing[:5] = np.nan
    keep = np.isfinite(missing)
    same_n = keep.sum()
    held = fit_nested(log_g[keep], proportion[keep], purity[keep], cnv[keep], missing[keep])
    lines.append(f"complete_case_n {held['n_complete']} expected {int(same_n)}")
    if held["n_complete"] != int(same_n):
        raise NumericalFailure("complete-case count changed between M0 and M1")
    p = np.array([0.001, 0.04, 0.2, 0.8])
    q = benjamini_hochberg(p)
    lines.append("bh " + " ".join(f"{value:.6g}" for value in q))
    if np.any(np.diff(q) < -1e-12):
        raise NumericalFailure("BH q values decreased with increasing p")
    # Jackknife on a fixed eight-observation problem with a known leave-one-out mean.
    small_p = np.array([0.2, 0.4, 0.1, 0.5, 0.3, 0.6, 0.15, 0.45])
    membership = [np.array([i]) for i in range(8)]

    def refit(gene, dropped):
        return float(small_p[gene] + (0.01 if dropped == gene else 0.0))

    jack = jackknife_mean(small_p, membership, refit)
    manual = np.mean([np.mean(np.r_[small_p[:i], small_p[i] + 0.01, small_p[i + 1:]]) for i in range(8)])
    # Each deletion changes only gene i, so the mean shifts by 0.01/8.
    expected = float(small_p.mean() + 0.01 / 8)
    lines.append(f"jackknife_center_check reported_mean {jack['full_mean']:.6g} loo_shift_target {expected:.6g} manual {manual:.6g}")
    if abs(jack["full_mean"] - float(small_p.mean())) > 1e-12 or abs(expected - manual) > 1e-12:
        raise NumericalFailure("jackknife diagnostic failed")
    meta = reml_meta(np.array([0.03, 0.04, 0.02, 0.05]), np.array([0.01, 0.01, 0.01, 0.01]))
    lines.append(f"reml_mu {meta['mu']:.6g} tau2 {meta['tau2']:.6g} ci {meta['ci_low']:.6g} {meta['ci_high']:.6g}")
    if not (meta["ci_low"] < meta["mu"] < meta["ci_high"]):
        raise NumericalFailure("Hartung-Knapp interval does not cover its mean")
    lines.append("SYNTHETIC_TESTS_PASS")
    return lines
