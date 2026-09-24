"""Donor-level prediction-powered inference (PPI++) for cell-type composition.

Estimand: for a cell type k, the mean over donors of the donor's proportion of
cells of type k, and its difference between two donor groups. Donors are the
sampling units; a random subset of donors in each group has gold-standard
labels, the rest only have model predictions.

    theta_hat = mean_L(Y) + lam * (mean_U(f) - mean_L(f))
    Var       = Var(Y - lam f)/n + lam^2 Var(f)/N

Y = gold donor proportion, f = predicted donor proportion, L = labeled donors,
U = unlabeled donors. The labeled donors MUST be a random subset of donors in
the group (not chosen by prediction confidence), otherwise the rectifier is
biased. With few labeled donors the t reference (df = n_labeled - 1) is used.

Baselines here are deliberately simple and labelled as proxies:
  * naive_plugin_diff       : predicted proportions treated as truth
  * cell_pooled_diff        : pooled cells, donors ignored (pseudoreplication)
  * confusion_corrected_diff: confusion-matrix inversion with donor bootstrap.
    This is a stand-in for the idea behind DCATS, NOT a run of DCATS itself.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def donor_proportions(labels, donors, n_classes: int, donor_ids=None):
    """(donor_ids, matrix n_donors x n_classes) of per-donor class proportions."""
    labels = np.asarray(labels).astype(int)
    donors = np.asarray(donors)
    ids = np.unique(donors) if donor_ids is None else np.asarray(donor_ids)
    P = np.zeros((ids.size, n_classes))
    for i, d in enumerate(ids):
        m = donors == d
        if m.any():
            P[i] = np.bincount(labels[m], minlength=n_classes) / m.sum()
    return ids, P


def _crit(alpha: float, df: int, dist: str) -> float:
    if dist == "t":
        return float(stats.t.ppf(1 - alpha / 2, max(df, 1)))
    return float(stats.norm.ppf(1 - alpha / 2))


def ppi_mean(y_lab, f_lab, f_unlab, alpha: float = 0.05, lam="auto", dist: str = "t") -> dict:
    """PPI++ mean estimate with a (1-alpha) interval over donors."""
    y = np.asarray(y_lab, float)
    fl = np.asarray(f_lab, float)
    fu = np.asarray(f_unlab, float)
    n, N = y.size, fu.size
    if n < 3 or N < 2:
        raise ValueError("need >=3 labeled and >=2 unlabeled donors")
    if lam == "auto":
        vf = np.var(np.concatenate([fl, fu]), ddof=1)
        cov = np.cov(y, fl, ddof=1)[0, 1]
        lam_v = 0.0 if vf <= 1e-15 else float(np.clip(cov / ((1 + n / N) * vf), 0.0, 1.0))
    else:
        lam_v = float(lam)
    est = y.mean() + lam_v * (fu.mean() - fl.mean())
    var = np.var(y - lam_v * fl, ddof=1) / n + lam_v ** 2 * np.var(fu, ddof=1) / N
    se = float(np.sqrt(var))
    c = _crit(alpha, n - 1, dist)
    return {"estimate": float(est), "se": se, "ci": (float(est - c * se), float(est + c * se)),
            "lambda": lam_v, "n_labeled": int(n), "n_unlabeled": int(N)}


def ppi_group_diff(g1, g0, alpha: float = 0.05, dist: str = "t", lam="auto") -> dict:
    """Difference of two groups; each g is a dict with y_lab, f_lab, f_unlab."""
    a = ppi_mean(g1["y_lab"], g1["f_lab"], g1["f_unlab"], alpha, lam, dist)
    b = ppi_mean(g0["y_lab"], g0["f_lab"], g0["f_unlab"], alpha, lam, dist)
    est = a["estimate"] - b["estimate"]
    se = float(np.hypot(a["se"], b["se"]))
    df = min(a["n_labeled"], b["n_labeled"]) - 1
    c = _crit(alpha, df, dist)
    return {"estimate": est, "se": se, "ci": (est - c * se, est + c * se)}


def naive_plugin_diff(f1_all, f0_all, alpha: float = 0.05) -> dict:
    """Welch interval on predicted donor proportions, predictions taken as truth."""
    a, b = np.asarray(f1_all, float), np.asarray(f0_all, float)
    va, vb = a.var(ddof=1) / a.size, b.var(ddof=1) / b.size
    est = a.mean() - b.mean()
    se = float(np.sqrt(va + vb))
    df = (va + vb) ** 2 / (va ** 2 / (a.size - 1) + vb ** 2 / (b.size - 1)) if se > 0 else 1
    c = float(stats.t.ppf(1 - alpha / 2, df))
    return {"estimate": float(est), "se": se, "ci": (float(est - c * se), float(est + c * se))}


def cell_pooled_diff(pred1, pred0, k: int, alpha: float = 0.05) -> dict:
    """Wald interval for pooled-cell proportions of predicted class k (donors ignored)."""
    p1 = float(np.mean(np.asarray(pred1) == k))
    p0 = float(np.mean(np.asarray(pred0) == k))
    se = float(np.sqrt(p1 * (1 - p1) / len(pred1) + p0 * (1 - p0) / len(pred0)))
    c = float(stats.norm.ppf(1 - alpha / 2))
    est = p1 - p0
    return {"estimate": est, "se": se, "ci": (est - c * se, est + c * se)}


def _confusion(true_lab, pred_lab, K):
    M = np.zeros((K, K))
    np.add.at(M, (np.asarray(true_lab).astype(int), np.asarray(pred_lab).astype(int)), 1.0)
    rs = M.sum(axis=1, keepdims=True)
    M = np.where(rs > 0, M / np.where(rs == 0, 1, rs), np.eye(K))
    return M


def confusion_corrected_diff(lab_true, lab_pred, lab_donor, unl1_props, unl0_props,
                             k: int, K: int, alpha: float = 0.05, n_boot: int = 300,
                             seed: int = 0) -> dict:
    """Confusion-matrix correction of mean predicted proportions, donor bootstrap CI.

    lab_* : cell-level arrays for the labeled donors (pooled over both groups).
    unl*_props : (n_donors x K) predicted donor proportions of unlabeled donors.
    """
    rng = np.random.default_rng(seed)
    lab_true, lab_pred, lab_donor = map(np.asarray, (lab_true, lab_pred, lab_donor))
    ids = np.unique(lab_donor)
    cells = {d: np.nonzero(lab_donor == d)[0] for d in ids}

    def correct(pmean, M):
        sol, *_ = np.linalg.lstsq(M.T, pmean, rcond=None)
        sol = np.clip(sol, 0, None)
        s = sol.sum()
        return sol / s if s > 0 else pmean

    def once(d_ids, u1, u0):
        idx = np.concatenate([cells[d] for d in d_ids])
        M = _confusion(lab_true[idx], lab_pred[idx], K)
        return correct(u1.mean(axis=0), M)[k] - correct(u0.mean(axis=0), M)[k]

    est = once(ids, unl1_props, unl0_props)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        d_b = rng.choice(ids, size=ids.size, replace=True)
        u1 = unl1_props[rng.integers(0, len(unl1_props), len(unl1_props))]
        u0 = unl0_props[rng.integers(0, len(unl0_props), len(unl0_props))]
        boots[b] = once(d_b, u1, u0)
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    return {"estimate": float(est), "se": float(boots.std(ddof=1)), "ci": (float(lo), float(hi))}
