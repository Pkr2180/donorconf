import math

import numpy as np

from donorconf import (ClasswiseConformal, DonorLevelCRC, TargetAdaptiveCRC, WeightedDonorCRC, score_matrix)
from donorconf.metrics import coverage_summary
from donorconf.simulate import SimConfig, fit_head, head_probs, make_cohort, shifted_means


def _world(seed=0, cells=80):
    cfg = SimConfig(cells_median=cells, donor_sd=0.8)
    rng = np.random.default_rng(seed)
    head = fit_head(make_cohort(cfg, 25, rng))
    return cfg, rng, head


def test_unit_weights_reproduce_unweighted_donor_crc():
    cfg, rng, head = _world()
    cal = make_cohort(cfg, 15, rng)
    P = head_probs(head, cal.X, cfg.n_classes)
    a = DonorLevelCRC(0.1).fit(P, cal.y, cal.donor)
    b = WeightedDonorCRC(0.1).fit(P, cal.y, cal.donor)
    assert a.threshold_ == b.threshold_


def test_no_target_donors_falls_back_to_source_only():
    cfg, rng, head = _world(1)
    cal = make_cohort(cfg, 15, rng)
    P = head_probs(head, cal.X, cfg.n_classes)
    m = TargetAdaptiveCRC(0.1).fit(P, cal.y, cal.donor)
    assert m.omega_ == 0.0 and m.threshold_ == DonorLevelCRC(0.1).fit(P, cal.y, cal.donor).threshold_


def test_weighted_threshold_grows_when_target_donors_are_harder():
    """Up-weighting hard target donors must not shrink the threshold."""
    cfg, rng, head = _world(2)
    src = make_cohort(cfg, 15, rng)
    tgt = make_cohort(cfg, 4, rng, means=shifted_means(cfg, 0.7), first_id=500)
    Ps, Pt = head_probs(head, src.X, cfg.n_classes), head_probs(head, tgt.X, cfg.n_classes)
    P = np.vstack([Ps, Pt]); y = np.concatenate([src.y, tgt.y]); d = np.concatenate([src.donor, tgt.donor])
    t0 = WeightedDonorCRC(0.1).fit(P, y, d, weights={i: 0.0 for i in np.unique(tgt.donor)}).threshold_
    t8 = WeightedDonorCRC(0.1).fit(P, y, d, weights={i: 8.0 for i in np.unique(tgt.donor)}, test_weight=8.0).threshold_
    assert t8 >= t0


def test_adaptive_improves_coverage_under_concept_shift():
    """Average over replications: with 4 labeled target donors coverage moves toward target vs source-only."""
    alpha, R = 0.10, 40
    cfg, rng, head = _world(3)
    cov_src, cov_ad = [], []
    for _ in range(R):
        src = make_cohort(cfg, 20, rng)
        tgt = make_cohort(cfg, 14, rng, means=shifted_means(cfg, 0.7), first_id=500)
        ids = np.unique(tgt.donor)
        lab = np.isin(tgt.donor, ids[:4]); ev = ~lab
        Ps, Pt = head_probs(head, src.X, cfg.n_classes), head_probs(head, tgt.X, cfg.n_classes)
        base = DonorLevelCRC(alpha).fit(Ps, src.y, src.donor)
        ad = TargetAdaptiveCRC(alpha).fit(Ps, src.y, src.donor, Pt[lab], tgt.y[lab], tgt.donor[lab])
        cov_src.append(coverage_summary(base.predict_sets(Pt[ev]), tgt.y[ev], tgt.donor[ev], alpha)["donor_mean_coverage"])
        cov_ad.append(coverage_summary(ad.predict_sets(Pt[ev]), tgt.y[ev], tgt.donor[ev], alpha)["donor_mean_coverage"])
    assert np.mean(cov_ad) >= np.mean(cov_src) - 0.01        # never materially worse
    assert abs(np.mean(cov_ad) - 0.9) <= abs(np.mean(cov_src) - 0.9) + 0.02


def test_classwise_gives_per_class_coverage_on_exchangeable_data():
    cfg, rng, head = _world(4, cells=150)
    cal = make_cohort(cfg, 30, rng); te = make_cohort(cfg, 30, rng, first_id=900)
    Pc, Pt = head_probs(head, cal.X, cfg.n_classes), head_probs(head, te.X, cfg.n_classes)
    m = ClasswiseConformal(0.1).fit(Pc, cal.y)
    sets = m.predict_sets(Pt)
    hit = sets[np.arange(len(te.y)), te.y]
    for c in range(cfg.n_classes):
        if (te.y == c).sum() > 100:
            assert hit[te.y == c].mean() >= 0.85
