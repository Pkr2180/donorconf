"""Real-data evaluation of target-adaptive donor-level CRC (framing B).

Per replication: fit a head on reference donors, take calibration (source) donors and target donors from
different datasets, label k random target donors, and score every method on the REMAINING target donors.
Exploratory; see results/real/ANALYSIS_LOCK.md (Addendum B written before this run).
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd

from _common import ROOT, stamp
from donorconf import (CellLevelConformal, DonorLevelCRC, DonorWeightedQuantile, TargetAdaptiveCRC)
from donorconf.head import EmbeddingHead
from donorconf.metrics import coverage_summary
from real_common import EMBEDDINGS, prepare_pair, save_meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--sets", nargs=3, required=True)
    ap.add_argument("--scheme", required=True)
    ap.add_argument("--permute-roles", action="store_true")
    ap.add_argument("--n-ref", type=int, default=30)
    ap.add_argument("--n-cal", type=int, default=30)
    ap.add_argument("--n-tgt", type=int, default=30)
    ap.add_argument("--min-cells", type=int, default=100)
    ap.add_argument("--min-donors", type=int, default=5)
    ap.add_argument("--ks", type=int, nargs="*", default=[3, 6])
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--embeddings", nargs="*", default=list(EMBEDDINGS))
    a = ap.parse_args()

    roles = list(itertools.permutations(a.sets)) if a.permute_roles else [tuple(a.sets)]
    rows, info = [], None
    for emb in a.embeddings:
        X, y, donors, dset, classes, info = prepare_pair(a.sets, a.scheme, a.min_cells, a.min_donors, emb)
        K = len(classes)
        for r in range(a.reps):
            rng = np.random.default_rng(r)
            ref_t, cal_t, tgt_t = roles[r % len(roles)]
            pick = lambda t, n: rng.choice(np.unique(donors[dset == t]), size=min(n, np.unique(donors[dset == t]).size), replace=False)
            ref_d, cal_d, tgt_d = pick(ref_t, a.n_ref), pick(cal_t, a.n_cal), pick(tgt_t, a.n_tgt)
            m_ref, m_cal, m_tgt = (np.isin(donors, s) for s in (ref_d, cal_d, tgt_d))
            head = EmbeddingHead(64, seed=r).fit(X[m_ref], y[m_ref], K)
            Pc, Pt = head.predict_proba(X[m_cal]), head.predict_proba(X[m_tgt])
            yc, dc, yt, dt = y[m_cal], donors[m_cal], y[m_tgt], donors[m_tgt]
            perm = rng.permutation(tgt_d)
            for k in a.ks:
                if k >= len(tgt_d) - 1:
                    continue
                lab = np.isin(dt, perm[:k]); ev = ~lab
                fits = {"cell_pooled": CellLevelConformal(a.alpha).fit(Pc, yc),
                        "donor_weighted": DonorWeightedQuantile(a.alpha).fit(Pc, yc, dc),
                        "donor_crc": DonorLevelCRC(a.alpha).fit(Pc, yc, dc),
                        "target_only_k": DonorWeightedQuantile(a.alpha).fit(Pt[lab], yt[lab], dt[lab])}
                ad = TargetAdaptiveCRC(a.alpha).fit(Pc, yc, dc, Pt[lab], yt[lab], dt[lab])
                fits["adaptive_crc"] = ad
                fits["adaptive_plugin"] = TargetAdaptiveCRC(a.alpha, finite_sample_correction=False).fit(
                    Pc, yc, dc, Pt[lab], yt[lab], dt[lab])
                for m, f in fits.items():
                    cs = coverage_summary(f.predict_sets(Pt[ev]), yt[ev], dt[ev], a.alpha)
                    rows.append({"embedding": emb, "rep": r, "k": k, "reference": ref_t, "calibration": cal_t,
                                 "target": tgt_t, "method": m, **cs,
                                 "omega": ad.omega_ if m == "adaptive_crc" else np.nan})
        print(f"{emb} done", flush=True)
    out = ROOT / "results" / "real" / a.name
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "per_rep.csv", index=False)
    save_meta(out / "provenance.json", stamp(vars(a)) | {"simulation_only": False, "exploratory": True, "data": info})
    print("written to", out)


if __name__ == "__main__":
    main()
