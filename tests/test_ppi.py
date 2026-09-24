import numpy as np
import pytest

from donorconf.ppi import (cell_pooled_diff, confusion_corrected_diff, donor_proportions,
                           naive_plugin_diff, ppi_group_diff, ppi_mean)


def test_donor_proportions_hand_computed():
    ids, P = donor_proportions([0, 0, 1, 1, 1, 0], [7, 7, 7, 9, 9, 9], 2)
    assert list(ids) == [7, 9]
    assert np.allclose(P, [[2 / 3, 1 / 3], [1 / 3, 2 / 3]])


def test_ppi_needs_enough_donors():
    with pytest.raises(ValueError):
        ppi_mean([0.1, 0.2], [0.1, 0.2], [0.1, 0.2, 0.3])


def _draw_group(rng, n_donors, truth, bias, noise=0.05):
    y = np.clip(rng.normal(truth, 0.08, n_donors), 0, 1)
    f = np.clip(y + bias + rng.normal(0, noise, n_donors), 0, 1)
    return y, f


def test_ppi_covers_truth_when_predictor_is_biased_but_plugin_does_not():
    rng = np.random.default_rng(0)
    truth, bias, R = 0.30, 0.10, 400
    cover_ppi = cover_plug = 0
    for _ in range(R):
        y, f = _draw_group(rng, 40, truth, bias)
        idx = rng.permutation(40)
        L, U = idx[:10], idx[10:]
        r = ppi_mean(y[L], f[L], f[U], alpha=0.05)
        cover_ppi += r["ci"][0] <= truth <= r["ci"][1]
        se = f.std(ddof=1) / np.sqrt(40)
        cover_plug += abs(f.mean() - truth) <= 1.96 * se
    assert cover_ppi / R >= 0.92          # near nominal 0.95
    assert cover_plug / R < 0.5           # plug-in is badly biased


def test_ppi_group_difference_covers_true_difference():
    rng = np.random.default_rng(1)
    R, cover = 300, 0
    for _ in range(R):
        g = {}
        for name, truth in (("g1", 0.35), ("g0", 0.25)):
            y, f = _draw_group(rng, 30, truth, 0.08)
            idx = rng.permutation(30)
            L, U = idx[:8], idx[8:]
            g[name] = {"y_lab": y[L], "f_lab": f[L], "f_unlab": f[U]}
        r = ppi_group_diff(g["g1"], g["g0"], alpha=0.05)
        cover += r["ci"][0] <= 0.10 <= r["ci"][1]
    assert cover / R >= 0.90


def test_baselines_return_intervals():
    rng = np.random.default_rng(2)
    a, b = rng.uniform(0.2, 0.4, 20), rng.uniform(0.1, 0.3, 20)
    r = naive_plugin_diff(a, b)
    assert r["ci"][0] < r["estimate"] < r["ci"][1]
    r = cell_pooled_diff(rng.integers(0, 3, 500), rng.integers(0, 3, 500), 1)
    assert r["ci"][0] < r["ci"][1]


def test_confusion_correction_recovers_known_mixing():
    rng = np.random.default_rng(5)
    K = 3
    true = rng.integers(0, K, 3000)
    pred = np.where(rng.random(3000) < 0.2, (true + 1) % K, true)  # 20% one-step errors
    donor = np.repeat(np.arange(30), 100)
    u1 = np.tile([0.5, 0.3, 0.2], (10, 1))
    u0 = np.tile([0.3, 0.3, 0.4], (10, 1))
    r = confusion_corrected_diff(true, pred, donor, u1, u0, k=0, K=K, n_boot=50)
    assert np.isfinite(r["estimate"]) and r["ci"][0] <= r["ci"][1]
