"""Nonconformity scores computed from class-probability matrices.

Every function returns an (n_cells, n_classes) matrix holding the score the
cell would receive if each candidate class were its true label. A prediction
set at threshold ``lam`` is ``score_matrix(...) <= lam``, so sets grow
monotonically with ``lam`` (required by conformal risk control).
"""
from __future__ import annotations

import numpy as np

SCORE_KINDS = ("lac", "aps")


def _check_probs(probs: np.ndarray) -> np.ndarray:
    p = np.asarray(probs, dtype=float)
    if p.ndim != 2:
        raise ValueError("probs must be a 2-D (cells x classes) array")
    if not np.all(np.isfinite(p)):
        raise ValueError("probs contains non-finite values")
    if p.min() < -1e-9 or not np.allclose(p.sum(axis=1), 1.0, atol=1e-4):
        raise ValueError("each row of probs must be a probability vector")
    return np.clip(p, 0.0, 1.0)


def score_matrix(probs: np.ndarray, kind: str = "lac") -> np.ndarray:
    """Score of every candidate label for every cell.

    lac: 1 - p_y (least ambiguous classifier; smallest average sets).
    aps: cumulative probability mass of all classes ranked at or above y
         (adaptive prediction sets, non-randomised, hence conservative).
    """
    p = _check_probs(probs)
    if kind == "lac":
        return 1.0 - p
    if kind == "aps":
        order = np.argsort(-p, axis=1, kind="stable")
        cum = np.cumsum(np.take_along_axis(p, order, axis=1), axis=1)
        out = np.empty_like(cum)
        np.put_along_axis(out, order, cum, axis=1)
        return np.clip(out, 0.0, 1.0)
    raise ValueError(f"unknown score kind {kind!r}; choose from {SCORE_KINDS}")


def true_label_scores(probs: np.ndarray, labels: np.ndarray, kind: str = "lac") -> np.ndarray:
    labels = np.asarray(labels)
    S = score_matrix(probs, kind)
    if labels.shape[0] != S.shape[0]:
        raise ValueError("labels and probs disagree on the number of cells")
    if labels.min() < 0 or labels.max() >= S.shape[1]:
        raise ValueError("labels must be integer class indices within probs columns")
    return S[np.arange(S.shape[0]), labels.astype(int)]
