"""Target-adaptive donor-level risk control (framing B).

Setting: calibrated on SOURCE donors, deployed on TARGET donors from a different dataset, with a few
(k) labeled target donors available. Target donors are given weight omega per donor relative to weight 1
for source donors.

Weighted donor-level conformal risk control
-------------------------------------------
    lam_hat = inf{ lam : ( sum_d w_d L_d(lam) + w_test ) / ( sum_d w_d + w_test ) <= alpha }

With all weights 1 and w_test = 1 this is exactly the unweighted donor-level CRC. Weighted-exchangeability
theory (Tibshirani et al. 2019; Angelopoulos et al. CRC) gives E[L_test(lam_hat)] <= alpha when the weights
are the true likelihood ratios of the test donor's distribution to each calibration donor's, with w_test the
ratio at the test donor. Those ratios are UNKNOWN here, so omega is a tuning parameter and the guarantee is
NOT claimed: the method is an empirically validated heuristic that is exact only in the limiting cases
omega = 0 (source-only CRC) and omega = 1 with target = source.

Choosing omega
--------------
Leave-one-target-donor-out (LOTO): for each candidate omega and each labeled target donor j, calibrate
with the remaining donors and record donor j's miss rate. Among omegas whose mean LOTO miss <= alpha pick the
smallest threshold (smallest sets); if none qualifies pick the omega with the smallest LOTO miss. With k of 3
to 6 donors this selection is noisy and optimistic; that is a stated limitation, not a solved problem.
"""
from __future__ import annotations

import math

import numpy as np

from .scores import score_matrix, true_label_scores

DEFAULT_OMEGAS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)


def donor_risk_matrix(scores: np.ndarray, donors: np.ndarray, cand: np.ndarray):
    """(donor_ids, matrix n_donors x n_cand) of the fraction of each donor's cells with score > cand."""
    uniq, inv = np.unique(donors, return_inverse=True)
    R = np.empty((uniq.size, cand.size))
    for d in range(uniq.size):
        sd = np.sort(scores[inv == d])
        R[d] = 1.0 - np.searchsorted(sd, cand, side="right") / sd.size
    return uniq, R


def weighted_crc_threshold(R: np.ndarray, w: np.ndarray, w_test: float, alpha: float, cand: np.ndarray) -> float:
    num = (w[:, None] * R).sum(axis=0) + w_test
    ok = np.nonzero(num / (w.sum() + w_test) <= alpha)[0]
    return float(cand[ok[0]]) if ok.size else math.inf


class WeightedDonorCRC:
    """Donor-level CRC with per-donor weights. weights: dict donor_id -> weight (missing = 1)."""

    def __init__(self, alpha: float = 0.1, score: str = "lac"):
        self.alpha, self.score = alpha, score
        self.threshold_ = None

    def fit(self, probs, labels, donors, weights: dict | None = None, test_weight: float = 1.0):
        s = true_label_scores(probs, labels, self.score)
        cand = np.unique(s)
        ids, R = donor_risk_matrix(s, np.asarray(donors), cand)
        w = np.array([(weights or {}).get(i, 1.0) for i in ids], dtype=float)
        self.threshold_ = weighted_crc_threshold(R, w, test_weight, self.alpha, cand)
        return self

    def predict_sets(self, probs) -> np.ndarray:
        return score_matrix(probs, self.score) <= self.threshold_


class TargetAdaptiveCRC:
    """Source donors + k labeled target donors; omega chosen by leave-one-target-donor-out."""

    def __init__(self, alpha: float = 0.1, score: str = "lac", omegas=DEFAULT_OMEGAS,
                 finite_sample_correction: bool = True):
        """finite_sample_correction=False drops the w_test term: a plug-in estimate with smaller sets
        and NO finite-sample statement, even in the limiting cases."""
        self.alpha, self.score, self.omegas = alpha, score, tuple(omegas)
        self.correct = finite_sample_correction
        self.threshold_ = None

    def fit(self, src_probs, src_labels, src_donors, tgt_probs=None, tgt_labels=None, tgt_donors=None):
        s_src = true_label_scores(src_probs, src_labels, self.score)
        k_donors = 0 if tgt_donors is None else np.unique(tgt_donors).size
        if k_donors:
            s_tgt = true_label_scores(tgt_probs, tgt_labels, self.score)
            cand = np.unique(np.concatenate([s_src, s_tgt]))
            _, R_t = donor_risk_matrix(s_tgt, np.asarray(tgt_donors), cand)
        else:
            cand = np.unique(s_src)
            R_t = np.zeros((0, cand.size))
        _, R_s = donor_risk_matrix(s_src, np.asarray(src_donors), cand)
        ns, k = R_s.shape[0], R_t.shape[0]
        self.n_source_, self.n_target_ = ns, k
        table = []
        if k == 0:
            self.omega_ = 0.0
            self.threshold_ = weighted_crc_threshold(R_s, np.ones(ns), 1.0 if self.correct else 0.0, self.alpha, cand)
            self.loto_table_ = table
            return self
        for om in self.omegas:
            miss, thr = [], []
            for j in range(k):
                keep = np.arange(k) != j
                if om == 0.0:
                    R, w, wt = R_s, np.ones(ns), (1.0 if self.correct else 0.0)
                else:
                    R = np.vstack([R_s, R_t[keep]])
                    w = np.concatenate([np.ones(ns), np.full(k - 1, om)])
                    wt = om if self.correct else 0.0
                lam = weighted_crc_threshold(R, w, wt, self.alpha, cand)
                thr.append(lam)
                miss.append(0.0 if math.isinf(lam) else float(R_t[j][np.searchsorted(cand, lam)]))
            finite = [1.0 if math.isinf(t) else t for t in thr]
            table.append({"omega": om, "loto_miss": float(np.mean(miss)), "mean_threshold": float(np.mean(finite))})
        feas = [t for t in table if t["loto_miss"] <= self.alpha]
        best = (min(feas, key=lambda t: (t["mean_threshold"], t["omega"])) if feas
                else min(table, key=lambda t: (t["loto_miss"], -t["mean_threshold"])))
        self.omega_ = best["omega"]
        self.loto_table_ = table
        if self.omega_ == 0.0:
            self.threshold_ = weighted_crc_threshold(R_s, np.ones(ns), 1.0 if self.correct else 0.0, self.alpha, cand)
        else:
            self.threshold_ = weighted_crc_threshold(
                np.vstack([R_s, R_t]), np.concatenate([np.ones(ns), np.full(k, self.omega_)]),
                self.omega_ if self.correct else 0.0, self.alpha, cand)
        return self

    def predict_sets(self, probs) -> np.ndarray:
        if self.threshold_ is None:
            raise RuntimeError("call fit first")
        return score_matrix(probs, self.score) <= self.threshold_
