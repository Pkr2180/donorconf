"""Coverage metrics that keep the donor, not the cell, as the reporting unit."""
from __future__ import annotations

import numpy as np


def per_donor_coverage(sets: np.ndarray, labels: np.ndarray, donors: np.ndarray):
    """Return (donor_ids, coverage, mean_set_size) with one entry per donor."""
    donors = np.asarray(donors)
    uniq, inv = np.unique(donors, return_inverse=True)
    hit = sets[np.arange(sets.shape[0]), np.asarray(labels).astype(int)].astype(float)
    size = sets.sum(axis=1).astype(float)
    cov = np.bincount(inv, weights=hit) / np.bincount(inv)
    sz = np.bincount(inv, weights=size) / np.bincount(inv)
    return uniq, cov, sz


def coverage_summary(sets, labels, donors, alpha: float, tol: float = 0.10) -> dict:
    """Donor-equal-weight and cell-pooled coverage side by side.

    tol: a donor counts as 'poorly covered' when its coverage < 1 - alpha - tol.
    """
    _, cov, sz = per_donor_coverage(sets, labels, donors)
    hit = sets[np.arange(sets.shape[0]), np.asarray(labels).astype(int)]
    return {
        "donor_mean_coverage": float(cov.mean()),
        "cell_pooled_coverage": float(hit.mean()),
        "donor_min_coverage": float(cov.min()),
        "donor_q10_coverage": float(np.quantile(cov, 0.10)),
        "frac_donors_poorly_covered": float((cov < 1 - alpha - tol).mean()),
        "mean_set_size": float(sz.mean()),
        "frac_empty_sets": float((sets.sum(axis=1) == 0).mean()),
        "n_donors": int(cov.size),
        "n_cells": int(sets.shape[0]),
    }


def donor_bootstrap_coverage_ci(sets, labels, donors, level: float = 0.95,
                                n_boot: int = 2000, seed: int = 0):
    """Percentile CI for donor-mean coverage by resampling donors.

    Uses labeled target donors only. The interval describes realised coverage
    on those donors; it does not repair a failed calibration.
    """
    _, cov, _ = per_donor_coverage(sets, labels, donors)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, cov.size, size=(n_boot, cov.size))
    means = cov[idx].mean(axis=1)
    lo, hi = np.quantile(means, [(1 - level) / 2, 1 - (1 - level) / 2])
    return float(cov.mean()), float(lo), float(hi)
