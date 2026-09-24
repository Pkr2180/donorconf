import json

import numpy as np

from donorconf.cli import main
from donorconf.simulate import SimConfig, fit_head, head_probs, make_cohort


def test_simulator_is_reproducible_and_shaped():
    cfg = SimConfig()
    a = make_cohort(cfg, 4, np.random.default_rng(0))
    b = make_cohort(cfg, 4, np.random.default_rng(0))
    assert np.array_equal(a.X, b.X) and a.n_donors == 4
    assert a.X.shape[1] == cfg.dim and a.X.shape[0] == a.y.size == a.donor.size


def test_size_difficulty_changes_donor_effects():
    rng = np.random.default_rng(0)
    c0 = make_cohort(SimConfig(size_difficulty=0.0), 30, rng)
    c1 = make_cohort(SimConfig(size_difficulty=1.0), 30, np.random.default_rng(0))
    assert not np.allclose(c0.X.std(), c1.X.std())


def test_cli_round_trip(tmp_path):
    cfg = SimConfig(cells_median=60)
    rng = np.random.default_rng(0)
    head = fit_head(make_cohort(cfg, 10, rng))
    cal = make_cohort(cfg, 12, rng, first_id=100)
    tgt = make_cohort(cfg, 6, rng, first_id=500)
    Pc, Pt = head_probs(head, cal.X, cfg.n_classes), head_probs(head, tgt.X, cfg.n_classes)
    np.savez(tmp_path / "cal.npz", probs=Pc, labels=cal.y, donors=cal.donor)
    np.savez(tmp_path / "tgt.npz", probs=Pt, labels=tgt.y, donors=tgt.donor)
    assert main(["calibrate", "--input", str(tmp_path / "cal.npz"), "--alpha", "0.2",
                 "--out", str(tmp_path / "m.json")]) == 0
    assert main(["predict", "--model", str(tmp_path / "m.json"), "--input", str(tmp_path / "tgt.npz"),
                 "--out", str(tmp_path / "sets.npz")]) == 0
    assert np.load(tmp_path / "sets.npz")["sets"].shape == Pt.shape
    assert main(["screen", "--cal", str(tmp_path / "cal.npz"), "--target", str(tmp_path / "tgt.npz"),
                 "--out", str(tmp_path / "s.json")]) == 0
    assert "p_value" in json.load(open(tmp_path / "s.json"))
    assert main(["report", "--model", str(tmp_path / "m.json"), "--input", str(tmp_path / "tgt.npz"),
                 "--out", str(tmp_path / "r.json")]) == 0
    assert "donor_mean_coverage_ci95" in json.load(open(tmp_path / "r.json"))


def test_lineage_mapping_regressions():
    from donorconf.lineage import blood_fine, oral_coarse
    assert oral_coarse("mast cell") == "mast"            # substring "t cell" once mapped this to T
    assert oral_coarse("T cell") == "T" and oral_coarse("regulatory T cell") == "T"
    assert oral_coarse("mature NK T cell") is None and oral_coarse("plasmacytoid dendritic cell") is None
    assert oral_coarse("plasmablast") == "plasma" and oral_coarse("lymphocyte") is None
    assert blood_fine("non-classical monocyte") == "mono_nonclassical"
    assert blood_fine("classical monocyte") == "mono_classical"
    assert blood_fine("mast cell") is None
    assert blood_fine("CD16-negative, CD56-bright natural killer cell, human") == "NK"
