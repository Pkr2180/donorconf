"""C: benchmark of reported-vs-realised coverage under dataset shift across blood datasets and embeddings.

Ordered role triples (reference, calibration, target) over blood datasets A-E; a fixed-seed subset is run, with
several replications of donor sampling each. Exploratory; see results/real/ANALYSIS_LOCK.md Addendum B.
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd

from _common import ROOT, stamp
from real_common import EMBEDDINGS, evaluate_split, prepare_pair, save_meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", default=["blood_A", "blood_B", "blood_C", "blood_D", "blood_E"])
    ap.add_argument("--scheme", default="blood_fine")
    ap.add_argument("--n-triples", type=int, default=24)
    ap.add_argument("--reps-per-triple", type=int, default=2)
    ap.add_argument("--n-donors", type=int, default=30)
    ap.add_argument("--min-cells", type=int, default=100)
    ap.add_argument("--min-donors", type=int, default=5)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--embeddings", nargs="*", default=list(EMBEDDINGS))
    ap.add_argument("--name", default="c_benchmark_blood")
    ap.add_argument("--torchcp", action="store_true", help="also run the published-package (torchCP) predictors")
    a = ap.parse_args()

    triples = list(itertools.permutations(a.datasets, 3))
    order = np.random.default_rng(0).permutation(len(triples))[: a.n_triples]
    triples = [triples[i] for i in order]
    rows, screens, infos = [], [], {}
    for emb in a.embeddings:
        for ti, (ref_t, cal_t, tgt_t) in enumerate(triples):
            try:
                X, y, donors, dset, classes, info = prepare_pair([ref_t, cal_t, tgt_t], a.scheme, a.min_cells, a.min_donors, emb)
            except ValueError as e:
                print("skip", ref_t, cal_t, tgt_t, e, flush=True)
                continue
            infos[f"{ref_t}|{cal_t}|{tgt_t}"] = info
            for r in range(a.reps_per_triple):
                seed = ti * 100 + r
                rng = np.random.default_rng(seed)
                pick = lambda t: rng.choice(np.unique(donors[dset == t]), size=min(a.n_donors, np.unique(donors[dset == t]).size), replace=False)
                res, sc = evaluate_split(X, y, donors, pick(ref_t), pick(cal_t), pick(tgt_t), a.alpha, len(classes), seed=seed, torchcp=a.torchcp)
                key = {"embedding": emb, "triple": ti, "rep": r, "reference": ref_t, "calibration": cal_t, "target": tgt_t,
                       "n_classes": len(classes)}
                rows += [{**key, **row} for row in res]
                screens.append({**key, **sc})
            print(f"{emb} triple {ti + 1}/{len(triples)} done", flush=True)
    out = ROOT / "results" / "real" / a.name
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "per_rep.csv", index=False)
    pd.DataFrame(screens).to_csv(out / "screens.csv", index=False)
    save_meta(out / "provenance.json", stamp(vars(a)) | {"simulation_only": False, "exploratory": True, "data": infos})
    print("written to", out)


if __name__ == "__main__":
    main()
