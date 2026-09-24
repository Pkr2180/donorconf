"""Cell-level baseline and donor-level conformal risk control.

Donor-level guarantee (Angelopoulos et al., conformal risk control)
-------------------------------------------------------------------
Let each calibration donor d have loss L_d(lam) = fraction of its cells whose
true label is NOT in the prediction set at threshold lam. L_d is non-increasing
in lam and bounded by 1. If a new donor is exchangeable with the calibration
donors, choosing

    lam_hat = inf{ lam : n/(n+1) * mean_d L_d(lam) + 1/(n+1) <= alpha }

gives E[L_new(lam_hat)] <= alpha. The expectation is over donors with EQUAL
weight, so the guarantee is about a typical donor, not a typical cell. It is a
statement in expectation, not a per-donor or per-cell-conditional guarantee.
It needs alpha >= 1/(n+1); otherwise the only valid threshold is infinity
(all classes) and the fit records that.
"""
from __future__ import annotations

import math

import numpy as np

from .scores import score_matrix, true_label_scores


def conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """Split-conformal threshold: the ceil((n+1)(1-alpha))-th smallest score."""
    s = np.sort(np.asarray(scores, dtype=float))
    n = s.size
    if n == 0:
        return math.inf
    k = math.ceil((n + 1) * (1.0 - alpha))
    return math.inf if k > n else float(s[k - 1])


class CellLevelConformal:
    """Standard split conformal that pools all calibration cells.

    Baseline only. It treats cells as exchangeable, which cells from a shared
    donor are not, so its guarantee does not apply to new donors as stated.
    """

    def __init__(self, alpha: float = 0.1, score: str = "lac"):
        self.alpha = alpha
        self.score = score
        self.threshold_: float | None = None

    def fit(self, probs, labels, donors=None):
        s = true_label_scores(probs, labels, self.score)
        self.threshold_ = conformal_quantile(s, self.alpha)
        self.n_cells_ = int(s.size)
        return self

    def predict_sets(self, probs) -> np.ndarray:
        if self.threshold_ is None:
            raise RuntimeError("call fit first")
        return score_matrix(probs, self.score) <= self.threshold_


class ClasswiseConformal:
    """Class-conditional (Mondrian) split conformal over pooled cells: one threshold per class.

    A stand-in for the classwise option of published scRNA conformal annotation, NOT that package. A class
    with too few calibration cells for the requested level gets an infinite threshold (always included).
    """

    def __init__(self, alpha: float = 0.1, score: str = "lac"):
        self.alpha, self.score = alpha, score
        self.thresholds_: np.ndarray | None = None

    def fit(self, probs, labels, donors=None):
        s = true_label_scores(probs, labels, self.score)
        labels = np.asarray(labels).astype(int)
        K = np.asarray(probs).shape[1]
        self.thresholds_ = np.array([conformal_quantile(s[labels == c], self.alpha) for c in range(K)])
        return self

    def predict_sets(self, probs) -> np.ndarray:
        if self.thresholds_ is None:
            raise RuntimeError("call fit first")
        return score_matrix(probs, self.score) <= self.thresholds_[None, :]


class DonorWeightedQuantile:
    """Pooled-cell quantile with each cell weighted 1/(cells in its donor).

    Removes the bias that arises when donor size is related to difficulty, and
    does not pay the 1/(n+1) finite-sample penalty of DonorLevelCRC. PLUG-IN
    ESTIMATOR: it has NO finite-sample guarantee, only the asymptotic target
    "donor-equal-weight mixture quantile". Use as an efficient companion to
    DonorLevelCRC, never as a replacement where a guarantee is claimed.
    """

    def __init__(self, alpha: float = 0.1, score: str = "lac"):
        self.alpha = alpha
        self.score = score
        self.threshold_: float | None = None

    def fit(self, probs, labels, donors):
        s = true_label_scores(probs, labels, self.score)
        donors = np.asarray(donors)
        _, inv, cnt = np.unique(donors, return_inverse=True, return_counts=True)
        w = 1.0 / cnt[inv]
        order = np.argsort(s, kind="stable")
        cw = np.cumsum(w[order]) / w.sum()
        k = int(np.searchsorted(cw, 1.0 - self.alpha - 1e-12, side="left"))
        self.threshold_ = float(s[order][min(k, s.size - 1)])
        self.n_donors_ = int(cnt.size)
        return self

    def predict_sets(self, probs) -> np.ndarray:
        if self.threshold_ is None:
            raise RuntimeError("call fit first")
        return score_matrix(probs, self.score) <= self.threshold_


class DonorLevelCRC:
    """Donor-level conformal risk control with the donor as exchangeable unit."""

    def __init__(self, alpha: float = 0.1, score: str = "lac"):
        self.alpha = alpha
        self.score = score
        self.threshold_: float | None = None

    def fit(self, probs, labels, donors):
        donors = np.asarray(donors)
        s = true_label_scores(probs, labels, self.score)
        uniq, inv = np.unique(donors, return_inverse=True)
        n = uniq.size
        if n < 2:
            raise ValueError("donor-level calibration needs at least 2 donors")
        cand = np.unique(s)
        risk = np.zeros(cand.size)
        for d in range(n):
            sd = np.sort(s[inv == d])
            risk += 1.0 - np.searchsorted(sd, cand, side="right") / sd.size
        risk /= n
        adjusted = n / (n + 1.0) * risk + 1.0 / (n + 1.0)
        ok = np.nonzero(adjusted <= self.alpha)[0]
        self.threshold_ = float(cand[ok[0]]) if ok.size else math.inf
        self.n_donors_ = int(n)
        self.min_alpha_ = 1.0 / (n + 1.0)
        self.vacuous_ = bool(math.isinf(self.threshold_))
        self.calibration_risk_ = float(risk[ok[0]]) if ok.size else float("nan")
        return self

    def predict_sets(self, probs) -> np.ndarray:
        if self.threshold_ is None:
            raise RuntimeError("call fit first")
        return score_matrix(probs, self.score) <= self.threshold_

    def to_dict(self) -> dict:
        return {
            "method": "donor_level_crc",
            "alpha": self.alpha,
            "score": self.score,
            "threshold": None if math.isinf(self.threshold_) else self.threshold_,
            "vacuous_all_classes": self.vacuous_,
            "n_calibration_donors": self.n_donors_,
            "min_valid_alpha": self.min_alpha_,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "DonorLevelCRC":
        m = cls(alpha=d["alpha"], score=d["score"])
        m.threshold_ = math.inf if d["threshold"] is None else float(d["threshold"])
        m.n_donors_ = d["n_calibration_donors"]
        m.min_alpha_ = d["min_valid_alpha"]
        m.vacuous_ = d["vacuous_all_classes"]
        return m
