"""Hierarchical single-cell simulator with donor effects.

SIMULATION ONLY. Output of this module is never biological evidence. It exists
to (a) test that the statistical code meets its stated guarantees under known
ground truth and (b) probe how the methods behave when assumptions are broken.

Generative model
----------------
class means  mu_k (fixed by world_seed, length class_sep)
donor d      random effect  b_d ~ N(0, sd_d^2 I) added to every cell
             composition    pi_d ~ Dirichlet(comp_conc * prior)
             cell count     n_d  ~ lognormal(median, sigma)
             sd_d = donor_sd * exp(size_difficulty * z_d), z_d the standardised
             log cell count; size_difficulty != 0 correlates donor size with
             classification difficulty, the case in which pooling cells across
             donors biases the calibration quantile.
cell         x = mu_y + b_d + offset + N(0, cell_noise^2 I)

Dataset shift is an ``offset`` added to target donors; label shift a different
``prior``; concept shift moves class-0's mean toward class-1's in the target.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from sklearn.linear_model import LogisticRegression


@dataclass
class SimConfig:
    n_classes: int = 6
    dim: int = 10
    class_sep: float = 2.5
    cell_noise: float = 1.0
    donor_sd: float = 0.8
    cells_median: int = 200
    cells_sigma: float = 0.6
    size_difficulty: float = 0.0
    comp_conc: float = 20.0
    world_seed: int = 0
    prior: tuple = field(default_factory=tuple)

    def base_prior(self) -> np.ndarray:
        if len(self.prior):
            p = np.asarray(self.prior, float)
            return p / p.sum()
        return np.full(self.n_classes, 1.0 / self.n_classes)


@dataclass
class Cohort:
    X: np.ndarray
    y: np.ndarray
    donor: np.ndarray

    @property
    def n_donors(self) -> int:
        return int(np.unique(self.donor).size)

    def subset_donors(self, donor_ids) -> "Cohort":
        m = np.isin(self.donor, donor_ids)
        return Cohort(self.X[m], self.y[m], self.donor[m])


def class_means(cfg: SimConfig) -> np.ndarray:
    rng = np.random.default_rng(cfg.world_seed)
    mu = rng.normal(size=(cfg.n_classes, cfg.dim))
    return cfg.class_sep * mu / np.linalg.norm(mu, axis=1, keepdims=True)


def shifted_means(cfg: SimConfig, concept: float = 0.0) -> np.ndarray:
    mu = class_means(cfg).copy()
    if concept:
        mu[0] = (1 - concept) * mu[0] + concept * mu[1]
    return mu


def make_cohort(cfg: SimConfig, n_donors: int, rng: np.random.Generator,
                prior=None, offset=None, means=None, first_id: int = 0) -> Cohort:
    mu = class_means(cfg) if means is None else means
    prior = cfg.base_prior() if prior is None else np.asarray(prior, float) / np.sum(prior)
    off = np.zeros(cfg.dim) if offset is None else np.asarray(offset, float)
    Xs, ys, ds = [], [], []
    for i in range(n_donors):
        z = rng.normal()
        n_cells = max(20, int(round(cfg.cells_median * np.exp(cfg.cells_sigma * z))))
        sd = cfg.donor_sd * np.exp(cfg.size_difficulty * z)
        b = rng.normal(0.0, sd, size=cfg.dim)
        pi = rng.dirichlet(cfg.comp_conc * prior)
        y = rng.choice(cfg.n_classes, size=n_cells, p=pi)
        X = mu[y] + b + off + rng.normal(0.0, cfg.cell_noise, size=(n_cells, cfg.dim))
        Xs.append(X)
        ys.append(y)
        ds.append(np.full(n_cells, first_id + i))
    return Cohort(np.vstack(Xs), np.concatenate(ys), np.concatenate(ds))


def fit_head(cohort: Cohort, C: float = 1.0) -> LogisticRegression:
    """Multinomial logistic head on embeddings (used for simulated and real data)."""
    return LogisticRegression(C=C, max_iter=1000).fit(cohort.X, cohort.y)


def head_probs(head: LogisticRegression, X: np.ndarray, n_classes: int) -> np.ndarray:
    """Probabilities with a column for every class, even if training missed one."""
    p = head.predict_proba(X)
    if p.shape[1] == n_classes:
        return p
    out = np.zeros((p.shape[0], n_classes))
    out[:, head.classes_] = p
    return out
