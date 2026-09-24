"""A fixed, prespecified classification head for embedding inputs.

standardise (reference statistics) -> PCA -> multinomial logistic regression.
The recipe is fixed in advance and identical for every embedding; it is not tuned
on calibration or target donors. Probabilities always have one column per class.
"""
from __future__ import annotations

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression


class EmbeddingHead:
    def __init__(self, n_components: int = 64, C: float = 1.0, seed: int = 0):
        self.n_components, self.C, self.seed = n_components, C, seed

    def fit(self, X: np.ndarray, y: np.ndarray, n_classes: int):
        X = np.asarray(X, dtype=np.float32)
        self.mu_ = X.mean(axis=0)
        self.sd_ = X.std(axis=0) + 1e-6
        Z = (X - self.mu_) / self.sd_
        k = min(self.n_components, Z.shape[1], Z.shape[0] - 1)
        self.pca_ = PCA(n_components=k, svd_solver="randomized", random_state=self.seed).fit(Z)
        F = self.pca_.transform(Z)
        self.clf_ = LogisticRegression(C=self.C, max_iter=500).fit(F, y)
        self.n_classes_ = n_classes
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        return self.pca_.transform((np.asarray(X, dtype=np.float32) - self.mu_) / self.sd_)

    def predict_proba(self, X: np.ndarray, feats: np.ndarray | None = None) -> np.ndarray:
        F = self.transform(X) if feats is None else feats
        p = self.clf_.predict_proba(F)
        out = np.zeros((p.shape[0], self.n_classes_))
        out[:, self.clf_.classes_] = p
        return out
