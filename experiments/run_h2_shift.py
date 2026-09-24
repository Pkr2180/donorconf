"""H2 kill test: does the label-free donor-level screen flag the shifts that break coverage?

Scenarios: no shift (type-I error), dataset shift (feature offset), label shift
(class-prior tilt), concept shift (class 0 drifts toward class 1). For each, the
donor-level CRC is calibrated on calibration donors and applied to target donors.
Reported per scenario: screen rejection rate; coverage without repair; coverage
after recalibrating with k labeled target donors; and how often coverage fails
given the screen did or did not flag.

Concept shift is expected to be INVISIBLE to the label-free screen. It is
included on purpose so the limitation is measured, not asserted.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from _common import write_results
from donorconf import DonorLevelCRC
from donorconf.metrics import coverage_summary
from donorconf.shift import donor_summaries, shift_screen
from donorconf.simulate import SimConfig, fit_head, head_probs, make_cohort, shifted_means


def scenario_specs(cfg):
    K = cfg.n_classes
    specs = [("none", 0.0)]
    specs += [("dataset_offset", m) for m in (0.5, 1.0, 2.0)]
    specs += [("label_shift", t) for t in (2.0, 5.0)]
    specs += [("concept_shift", c) for c in (0.3, 0.6)]
    return specs


def run(alpha, reps, n_cal, n_tgt, k_recal, seed):
    cfg = SimConfig()
    K = cfg.n_classes
    rng = np.random.default_rng(seed)
    head = fit_head(make_cohort(cfg, 30, rng))
    direction = rng.normal(size=cfg.dim)
    direction /= np.linalg.norm(direction)
    rows = []
    for kind, mag in scenario_specs(cfg):
        res = []
        for r in range(reps):
            cal = make_cohort(cfg, n_cal, rng, first_id=0)
            kw = {}
            if kind == "dataset_offset":
                kw["offset"] = mag * direction
            elif kind == "label_shift":
                pr = np.ones(K)
                pr[0] = mag
                kw["prior"] = pr
            elif kind == "concept_shift":
                kw["means"] = shifted_means(cfg, mag)
            tgt = make_cohort(cfg, n_tgt, rng, first_id=10_000, **kw)
            Pc, Pt = head_probs(head, cal.X, K), head_probs(head, tgt.X, K)
            centroid = cal.X.mean(axis=0)
            _, sc, _ = donor_summaries(Pc, cal.donor, cal.X, centroid)
            _, st, _ = donor_summaries(Pt, tgt.donor, tgt.X, centroid)
            flag = shift_screen(sc, st, level=0.05, n_perm=499, max_exact=0,
                                seed=int(rng.integers(1e9)))["flagged"]

            donors = np.unique(tgt.donor)
            lab = rng.choice(donors, size=k_recal, replace=False)
            m_lab = np.isin(tgt.donor, lab)
            ev = ~m_lab
            base = DonorLevelCRC(alpha).fit(Pc, cal.y, cal.donor)
            cov0 = coverage_summary(base.predict_sets(Pt[ev]), tgt.y[ev], tgt.donor[ev], alpha)["donor_mean_coverage"]
            rec = DonorLevelCRC(alpha).fit(
                np.vstack([Pc, Pt[m_lab]]), np.concatenate([cal.y, tgt.y[m_lab]]),
                np.concatenate([cal.donor, tgt.donor[m_lab]]))
            cov1 = coverage_summary(rec.predict_sets(Pt[ev]), tgt.y[ev], tgt.donor[ev], alpha)["donor_mean_coverage"]
            res.append((flag, cov0, cov1))
        a = np.array(res, dtype=float)
        fail = a[:, 1] < (1 - alpha - 0.05)
        flag = a[:, 0].astype(bool)
        rows.append({
            "shift": kind, "magnitude": mag, "target": 1 - alpha,
            "screen_rejection_rate": flag.mean(),
            "coverage_no_repair": a[:, 1].mean(),
            "coverage_after_recal": a[:, 2].mean(),
            "p_fail_no_repair": fail.mean(),
            "p_fail_given_flagged": fail[flag].mean() if flag.any() else np.nan,
            "p_fail_given_not_flagged": fail[~flag].mean() if (~flag).any() else np.nan,
        })
        print(rows[-1], flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=150)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--n-cal", type=int, default=20)
    ap.add_argument("--n-tgt", type=int, default=12)
    ap.add_argument("--k-recal", type=int, default=3)
    ap.add_argument("--seed", type=int, default=2027)
    a = ap.parse_args()
    df = run(a.alpha, a.reps, a.n_cal, a.n_tgt, a.k_recal, a.seed)
    out = write_results("h2_shift", df, vars(a))
    print(df.round(3).to_string(index=False))
    print("written to", out)


if __name__ == "__main__":
    main()
