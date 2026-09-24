"""Label-free, donor-level shift screen.

The test asks whether the target donors' summary statistics look exchangeable
with the calibration donors' summaries. Donors, not cells, are permuted, so the
p-value is valid under donor exchangeability of the summaries.

What a pass means
-----------------
A large p-value is the ABSENCE of detected shift, not a certificate. The screen
is label-free, so it cannot see pure concept shift (P(label | features) changes
while features and confidences do not). Power is limited by donor numbers; the
smallest attainable p-value is reported so a user can see when the test cannot
possibly reject. Treat "pass" as "no evidence against", never as "guaranteed".
"""
from __future__ import annotations

import itertools
import math

import numpy as np


def donor_summaries(probs, donors, embeddings=None, centroid=None, cov_inv=None):
    """One row per donor: mean max-probability, mean entropy, [embedding distance].

    embeddings + centroid (+ optional inverse covariance) add the mean distance
    of the donor's cells to the calibration centroid.
    """
    probs = np.asarray(probs, dtype=float)
    donors = np.asarray(donors)
    uniq, inv = np.unique(donors, return_inverse=True)
    cnt = np.bincount(inv)
    maxp = np.bincount(inv, weights=probs.max(axis=1)) / cnt
    ent = np.bincount(inv, weights=-(probs * np.log(np.clip(probs, 1e-12, 1))).sum(axis=1)) / cnt
    cols = [maxp, ent]
    names = ["mean_max_prob", "mean_entropy"]
    if embeddings is not None and centroid is not None:
        diff = np.asarray(embeddings, dtype=float) - np.asarray(centroid, dtype=float)
        if cov_inv is None:
            dist = np.sqrt((diff ** 2).sum(axis=1))
        else:
            dist = np.sqrt(np.einsum("ij,jk,ik->i", diff, cov_inv, diff))
        cols.append(np.bincount(inv, weights=dist) / cnt)
        names.append("mean_centroid_distance")
    return uniq, np.column_stack(cols), names


def _stat(a: np.ndarray, b: np.ndarray, scale: np.ndarray) -> float:
    d = (a.mean(axis=0) - b.mean(axis=0)) / scale
    return float(d @ d)


def donor_permutation_test(cal_summ, tgt_summ, n_perm: int = 9999,
                           max_exact: int = 20000, seed: int = 0) -> dict:
    """Permutation test on donor summaries (squared standardised mean difference)."""
    A = np.atleast_2d(np.asarray(cal_summ, dtype=float))
    B = np.atleast_2d(np.asarray(tgt_summ, dtype=float))
    if A.shape[1] != B.shape[1]:
        raise ValueError("summary dimensions differ")
    nA, nB = A.shape[0], B.shape[0]
    if nA < 2 or nB < 1:
        raise ValueError("need >=2 calibration donors and >=1 target donor")
    pool = np.vstack([A, B])
    scale = pool.std(axis=0, ddof=1)
    scale[scale == 0] = 1.0
    obs = _stat(A, B, scale)
    n = nA + nB
    total = math.comb(n, nB)
    if total <= max_exact:
        cnt = 0
        for idx in itertools.combinations(range(n), nB):
            mask = np.zeros(n, bool)
            mask[list(idx)] = True
            if _stat(pool[~mask], pool[mask], scale) >= obs - 1e-12:
                cnt += 1
        p = cnt / total
        exact, min_p = True, 1.0 / total
    else:
        rng = np.random.default_rng(seed)
        cnt = 0
        for _ in range(n_perm):
            perm = rng.permutation(n)
            if _stat(pool[perm[:nA]], pool[perm[nA:]], scale) >= obs - 1e-12:
                cnt += 1
        p = (cnt + 1) / (n_perm + 1)
        exact, min_p = False, 1.0 / (n_perm + 1)
    return {"statistic": obs, "p_value": float(p), "exact": exact,
            "min_attainable_p": float(min_p), "n_cal_donors": int(nA), "n_target_donors": int(nB)}


def shift_screen(cal_summ, tgt_summ, level: float = 0.05, **kw) -> dict:
    """Flag the target dataset when the donor-level test rejects at ``level``."""
    res = donor_permutation_test(cal_summ, tgt_summ, **kw)
    res["level"] = level
    res["flagged"] = bool(res["p_value"] <= level)
    res["can_reject"] = bool(res["min_attainable_p"] <= level)
    res["note"] = (
        "Label-free: cannot detect concept shift. Not flagged means no detected "
        "evidence against exchangeability, not a coverage certificate."
    )
    return res
