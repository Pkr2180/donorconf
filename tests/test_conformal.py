import math

import numpy as np
import pytest

from donorconf import CellLevelConformal, DonorLevelCRC, conformal_quantile, score_matrix
from donorconf.metrics import coverage_summary
from donorconf.simulate import SimConfig, fit_head, head_probs, make_cohort


def test_conformal_quantile_known_values():
    s = np.arange(1, 11, dtype=float)  # 1..10
    assert conformal_quantile(s, 0.2) == 9.0  # ceil(11*0.8)=9
    assert math.isinf(conformal_quantile(s, 0.05))  # ceil(11*0.95)=11 > 10
    assert math.isinf(conformal_quantile(np.array([]), 0.1))


def test_score_matrices_hand_computed():
    p = np.array([[0.6, 0.3, 0.1]])
    assert np.allclose(score_matrix(p, "lac"), [[0.4, 0.7, 0.9]])
    # aps: class 0 -> 0.6, class 1 -> 0.9, class 2 -> 1.0
    assert np.allclose(score_matrix(p, "aps"), [[0.6, 0.9, 1.0]])


def test_bad_probs_rejected():
    with pytest.raises(ValueError):
        score_matrix(np.array([[0.5, 0.2]]))
    with pytest.raises(ValueError):
        score_matrix(np.array([[np.nan, 1.0]]))


def _tiny_calibration(n_donors=10, seed=0):
    cfg = SimConfig(cells_median=60)
    rng = np.random.default_rng(seed)
    ref = make_cohort(cfg, 10, rng)
    head = fit_head(ref)
    cal = make_cohort(cfg, n_donors, rng, first_id=100)
    return cfg, head, cal


def test_sets_grow_with_threshold():
    cfg, head, cal = _tiny_calibration()
    P = head_probs(head, cal.X, cfg.n_classes)
    S = score_matrix(P)
    small = S <= 0.2
    big = S <= 0.6
    assert np.all(big[small])


def test_donor_crc_vacuous_when_alpha_below_one_over_n_plus_one():
    cfg, head, cal = _tiny_calibration(n_donors=5)
    P = head_probs(head, cal.X, cfg.n_classes)
    m = DonorLevelCRC(alpha=0.10).fit(P, cal.y, cal.donor)  # 1/(5+1)=0.167 > 0.10
    assert m.vacuous_ and math.isinf(m.threshold_)
    assert m.predict_sets(P).all()


def test_donor_crc_needs_two_donors():
    cfg, head, cal = _tiny_calibration(n_donors=3)
    P = head_probs(head, cal.X, cfg.n_classes)
    one = cal.donor == cal.donor[0]
    with pytest.raises(ValueError):
        DonorLevelCRC().fit(P[one], cal.y[one], cal.donor[one])


def test_donor_crc_finite_sample_guarantee_monte_carlo():
    """E[new-donor miscoverage] <= alpha under donor exchangeability (equal weights).

    Fixed classifier, fresh calibration and test donors each replication. Allow
    3 Monte-Carlo standard errors of slack; a broken implementation misses by far more.
    """
    alpha, n_cal, R = 0.10, 15, 300
    cfg = SimConfig(cells_median=80, donor_sd=1.0)
    rng = np.random.default_rng(1)
    head = fit_head(make_cohort(cfg, 25, rng))
    miss = []
    for r in range(R):
        cal = make_cohort(cfg, n_cal, rng, first_id=0)
        te = make_cohort(cfg, 10, rng, first_id=1000)
        m = DonorLevelCRC(alpha).fit(head_probs(head, cal.X, cfg.n_classes), cal.y, cal.donor)
        sets = m.predict_sets(head_probs(head, te.X, cfg.n_classes))
        miss.append(1 - coverage_summary(sets, te.y, te.donor, alpha)["donor_mean_coverage"])
    miss = np.array(miss)
    se = miss.std(ddof=1) / np.sqrt(R)
    assert miss.mean() <= alpha + 3 * se, (miss.mean(), alpha, se)


def test_cell_level_baseline_runs_and_is_not_worse_than_trivial():
    cfg, head, cal = _tiny_calibration(n_donors=10)
    P = head_probs(head, cal.X, cfg.n_classes)
    m = CellLevelConformal(0.1).fit(P, cal.y)
    sets = m.predict_sets(P)
    assert sets[np.arange(len(cal.y)), cal.y].mean() >= 0.85


def test_donor_crc_is_tight_with_many_donors_and_conservative_with_few():
    """Validity alone is satisfied by all-class sets, so also pin how tight it is.

    Measured at development time (150 reps): miss = 0.045 at 15 donors, 0.086 at
    60 donors for alpha=0.10; 0.182 at 60 and 0.194 at 120 donors for alpha=0.20.
    """
    alpha, R = 0.20, 100
    cfg = SimConfig(cells_median=80, donor_sd=1.0)
    rng = np.random.default_rng(7)
    head = fit_head(make_cohort(cfg, 25, rng))

    def mean_miss(n_cal):
        out = []
        for _ in range(R):
            cal = make_cohort(cfg, n_cal, rng)
            te = make_cohort(cfg, 10, rng, first_id=1000)
            m = DonorLevelCRC(alpha).fit(head_probs(head, cal.X, cfg.n_classes), cal.y, cal.donor)
            s = m.predict_sets(head_probs(head, te.X, cfg.n_classes))
            out.append(1 - coverage_summary(s, te.y, te.donor, alpha)["donor_mean_coverage"])
        return float(np.mean(out))

    few, many = mean_miss(8), mean_miss(60)
    assert many <= alpha + 0.015           # valid
    assert many >= alpha - 0.05            # and not vacuous
    assert few < many                      # conservatism shrinks as donors grow


def test_donor_weighted_quantile_hand_example():
    from donorconf import DonorWeightedQuantile
    # donor A: 1 cell, LAC score 0.9 (total weight 1); donor B: 9 cells, score 0.1 (total weight 1)
    probs = np.array([[0.1, 0.9]] + [[0.9, 0.1]] * 9)
    labels = np.zeros(10, dtype=int)
    donors = np.array([0] + [1] * 9)
    # weighted mass at 0.1 is exactly 0.5, so alpha=0.5 -> 0.1 but alpha=0.4 needs mass 0.6 -> 0.9
    assert abs(DonorWeightedQuantile(alpha=0.5).fit(probs, labels, donors).threshold_ - 0.1) < 1e-9
    assert abs(DonorWeightedQuantile(alpha=0.4).fit(probs, labels, donors).threshold_ - 0.9) < 1e-9
    # pooling cells lets the big donor dominate: 9/10 of the cells score 0.1
    assert CellLevelConformal(alpha=0.4).fit(probs, labels).threshold_ <= 0.1 + 1e-9
