"""Leave-dataset-out experiments on real Census data (R2 blood, R3 oral). See results/real/ANALYSIS_LOCK.md.

Reference, calibration and target donors come from DIFFERENT datasets, so target donors differ from
calibration donors by study, protocol and (for R3) disease: a genuine dataset shift.

    python real_crossdataset.py --name r2_blood --sets blood_A blood_B blood_C --permute-roles \
        --n-ref 30 --n-cal 30 --n-tgt 30 --min-cells 100
    python real_crossdataset.py --name r3_oral --sets gingiva oral_normal_multisite oral_scc \
        --n-ref 20 --n-cal 14 --n-tgt 10 --min-cells 40 --min-donors 4
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd

from _common import ROOT, stamp
from donorconf.lineage import SCHEMES, map_labels
from real_common import EMBEDDINGS, encode, evaluate_split, load, save_meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--sets", nargs=3, required=True, help="tags for reference / calibration / target")
    ap.add_argument("--permute-roles", action="store_true", help="cycle all 6 role assignments")
    ap.add_argument("--n-ref", type=int, default=30)
    ap.add_argument("--n-cal", type=int, default=30)
    ap.add_argument("--n-tgt", type=int, default=30)
    ap.add_argument("--min-cells", type=int, default=100, help="per dataset, for a class to enter the vocabulary")
    ap.add_argument("--min-donors", type=int, default=5, help="per dataset")
    ap.add_argument("--scheme", default="none", choices=sorted(SCHEMES),
                    help="coarse-lineage mapping applied to every dataset (see config/lineage_audit.csv)")
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--embeddings", nargs="*", default=list(EMBEDDINGS))
    a = ap.parse_args()

    data = {t: load(t) for t in a.sets}
    for z in data.values():
        z["labels"] = map_labels(z["labels"], a.scheme)     # '' = excluded (unmapped); never a class
    vocab = None
    for t, z in data.items():
        ok = set()
        for c in np.unique(z["labels"]):
            if c == "":
                continue
            m = z["labels"] == c
            if m.sum() >= a.min_cells and np.unique(z["donors"][m]).size >= a.min_donors:
                ok.add(c)
        vocab = ok if vocab is None else vocab & ok
    classes = np.array(sorted(vocab))
    if classes.size < 3:
        raise SystemExit(f"only {classes.size} shared cell types: {classes.tolist()}; cannot run")
    info = {"classes": classes.tolist(), "n_classes": int(classes.size), "sets": a.sets, "per_set": {}}
    X_all = {e: [] for e in a.embeddings}
    y_all, d_all, ds_all = [], [], []
    for t, z in data.items():
        y = encode(z["labels"], classes)
        k = y >= 0
        info["per_set"][t] = {"cells": int(len(y)), "kept": int(k.sum()),
                              "donors": int(np.unique(z["donors"][k]).size)}
        for e in a.embeddings:
            X_all[e].append(z[e][k])
        y_all.append(y[k]); d_all.append(z["donors"][k]); ds_all.append(np.full(k.sum(), t))
    y, donors, dset = np.concatenate(y_all), np.concatenate(d_all), np.concatenate(ds_all)
    Xs = {e: np.vstack(v) for e, v in X_all.items()}
    print(info, flush=True)

    role_sets = list(itertools.permutations(a.sets)) if a.permute_roles else [tuple(a.sets)]
    rows, screens = [], []
    for e in a.embeddings:
        for r in range(a.reps):
            rng = np.random.default_rng(r)
            ref_t, cal_t, tgt_t = role_sets[r % len(role_sets)]
            pick = lambda t, n: rng.choice(np.unique(donors[dset == t]), size=min(n, np.unique(donors[dset == t]).size), replace=False)
            res, sc = evaluate_split(Xs[e], y, donors, pick(ref_t, a.n_ref), pick(cal_t, a.n_cal),
                                     pick(tgt_t, a.n_tgt), a.alpha, len(classes), seed=r)
            for row in res:
                rows.append({"embedding": e, "rep": r, "reference": ref_t, "calibration": cal_t, "target": tgt_t, **row})
            screens.append({"embedding": e, "rep": r, "reference": ref_t, "calibration": cal_t, "target": tgt_t, **sc})
        print(f"{e} done", flush=True)
    out = ROOT / "results" / "real" / a.name
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "per_rep.csv", index=False)
    pd.DataFrame(screens).to_csv(out / "screens.csv", index=False)
    save_meta(out / "provenance.json", stamp(vars(a)) | {"simulation_only": False, "exploratory": True, "data": info})
    print("written to", out)


if __name__ == "__main__":
    main()
