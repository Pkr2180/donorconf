"""Simulation check of target-adaptive donor-level CRC (framing B). SIMULATION ONLY.

Shifts: none, dataset offset, concept shift. k labeled target donors are used for adaptation; all methods are
scored on the OTHER target donors. Methods:
  source_crc        donor-level CRC on source donors only
  source_weighted   donor-weighted plug-in quantile on source only
  target_only_k     donor-weighted plug-in quantile on the k labeled target donors only
  pooled_plugin     donor-weighted plug-in quantile on source + labeled target donors (equal donor weight)
  adaptive_crc      B: weighted donor CRC, omega chosen by leave-one-target-donor-out
  adaptive_plugin   B without the finite-sample term: smaller sets, no finite-sample statement
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd

from _common import write_results
from donorconf import DonorLevelCRC, DonorWeightedQuantile, TargetAdaptiveCRC
from donorconf.metrics import coverage_summary
from donorconf.simulate import SimConfig, fit_head, head_probs, make_cohort, shifted_means


def run(alpha, reps, n_src, n_tgt, ks, seed):
    cfg = SimConfig()
    K = cfg.n_classes
    rng = np.random.default_rng(seed)
    head = fit_head(make_cohort(cfg, 30, rng))
    direction = rng.normal(size=cfg.dim); direction /= np.linalg.norm(direction)
    scen = [("none", 0.0), ("dataset_offset", 1.0), ("dataset_offset", 2.0), ("concept_shift", 0.4), ("concept_shift", 0.7)]
    rows = []
    for (kind, mag), k in itertools.product(scen, ks):
        res = {m: [] for m in ("source_crc", "source_weighted", "target_only_k", "pooled_plugin", "adaptive_crc", "adaptive_plugin")}
        omegas = []
        for _ in range(reps):
            src = make_cohort(cfg, n_src, rng)
            kw = {}
            if kind == "dataset_offset": kw["offset"] = mag * direction
            if kind == "concept_shift": kw["means"] = shifted_means(cfg, mag)
            tgt = make_cohort(cfg, n_tgt, rng, first_id=10_000, **kw)
            ids = rng.permutation(np.unique(tgt.donor))
            lab = np.isin(tgt.donor, ids[:k]); ev = ~lab
            Ps, Pt = head_probs(head, src.X, K), head_probs(head, tgt.X, K)
            fits = {
                "source_crc": DonorLevelCRC(alpha).fit(Ps, src.y, src.donor),
                "source_weighted": DonorWeightedQuantile(alpha).fit(Ps, src.y, src.donor),
                "target_only_k": DonorWeightedQuantile(alpha).fit(Pt[lab], tgt.y[lab], tgt.donor[lab]),
                "pooled_plugin": DonorWeightedQuantile(alpha).fit(
                    np.vstack([Ps, Pt[lab]]), np.concatenate([src.y, tgt.y[lab]]), np.concatenate([src.donor, tgt.donor[lab]])),
            }
            ad = TargetAdaptiveCRC(alpha).fit(Ps, src.y, src.donor, Pt[lab], tgt.y[lab], tgt.donor[lab])
            fits["adaptive_crc"] = ad
            fits["adaptive_plugin"] = TargetAdaptiveCRC(alpha, finite_sample_correction=False).fit(
                Ps, src.y, src.donor, Pt[lab], tgt.y[lab], tgt.donor[lab])
            omegas.append(ad.omega_)
            for m, f in fits.items():
                cs = coverage_summary(f.predict_sets(Pt[ev]), tgt.y[ev], tgt.donor[ev], alpha)
                res[m].append((cs["donor_mean_coverage"], cs["mean_set_size"]))
        for m, v in res.items():
            a = np.array(v)
            rows.append({"shift": kind, "magnitude": mag, "k_labeled_target": k, "method": m, "target": 1 - alpha,
                         "coverage": a[:, 0].mean(), "cov_sd": a[:, 0].std(ddof=1), "set_size": a[:, 1].mean(),
                         "p_under5": float((a[:, 0] < 1 - alpha - 0.05).mean()),
                         "mean_omega": float(np.mean(omegas)) if m == "adaptive_crc" else np.nan})
        print(kind, mag, "k=", k, "done", flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--n-src", type=int, default=20)
    ap.add_argument("--n-tgt", type=int, default=20)
    ap.add_argument("--ks", type=int, nargs="*", default=[3, 6])
    ap.add_argument("--seed", type=int, default=2029)
    a = ap.parse_args()
    df = run(a.alpha, a.reps, a.n_src, a.n_tgt, a.ks, a.seed)
    write_results("sim_b_adaptive", df, vars(a))
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
