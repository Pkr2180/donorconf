import numpy as np
import pytest

from donorconf import CellLevelConformal
from donorconf.lineage import lung_coarse


def _need_torchcp():
    pytest.importorskip("torch")
    from donorconf.baseline_torchcp import _import_torchcp
    try:
        _import_torchcp()          # inserts a torchsort stub when torchsort is not installed
    except ImportError:
        pytest.skip("torchcp not installed")


def _exchangeable(n=3000, k=5, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(n, k)) * 1.5
    y = rng.integers(0, k, n)
    z[np.arange(n), y] += 2
    P = np.exp(z)
    P /= P.sum(1, keepdims=True)
    return P, y


def test_torchcp_standard_thr_matches_pooled_cell():
    _need_torchcp()
    from donorconf.baseline_torchcp import torchcp_sets
    P, y = _exchangeable()
    ours = CellLevelConformal(0.1).fit(P[:2000], y[:2000]).predict_sets(P[2000:])
    theirs = torchcp_sets(P[:2000], y[:2000], P[2000:], 0.1, "tcp_standard_thr")
    assert theirs.dtype == bool and theirs.shape == ours.shape
    assert (theirs != ours).mean() < 0.002       # same rule; only float32 ties may differ


def test_torchcp_variants_cover_on_exchangeable_data():
    _need_torchcp()
    from donorconf.baseline_torchcp import VARIANTS, torchcp_sets
    P, y = _exchangeable(seed=1)
    for v in VARIANTS:
        S = torchcp_sets(P[:2000], y[:2000], P[2000:], 0.1, v)
        cov = S[np.arange(1000), y[2000:]].mean()
        assert 0.85 < cov < 0.95, (v, cov)      # regression: a 0/1 tensor once read as index lists gave 0.38


def test_lung_coarse_rules():
    assert lung_coarse("type II pneumocyte") == "AT2"
    assert lung_coarse("alveolar macrophage") == "macrophage"
    assert lung_coarse("CD4-positive, alpha-beta T cell") == "T"
    assert lung_coarse("mast cell") is None          # must not match the T-cell rule
    assert lung_coarse("pancreatic A cell") is None
    assert lung_coarse("unknown") is None
