"""R1: real gingiva atlas, repeated donor-level splits (see results/real/ANALYSIS_LOCK.md)."""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from _common import ROOT, stamp
from real_common import EMBEDDINGS, choose_classes, encode, evaluate_split, load, save_meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--embeddings", nargs="*", default=list(EMBEDDINGS))
    a = ap.parse_args()

    z = load("gingiva")
    labels, donors = z["labels"], z["donors"]
    classes = choose_classes(labels, donors, min_cells=150, min_donors=12)
    y = encode(labels, classes)
    keep = y >= 0
    info = {"classes": classes.tolist(), "n_classes": int(classes.size),
            "cells_total": int(len(y)), "cells_dropped_outside_vocabulary": int((~keep).sum()),
            "n_donors": int(np.unique(donors).size)}
    print(info, flush=True)
    ud = np.unique(donors[keep])
    rows, screens = [], []
    for emb in a.embeddings:
        X = z[emb][keep]
        yy, dd = y[keep], donors[keep]
        for r in range(a.reps):
            rng = np.random.default_rng(r)
            perm = rng.permutation(ud)
            ref, cal, tgt = perm[:12], perm[12:24], perm[24:34]
            res, sc = evaluate_split(X, yy, dd, ref, cal, tgt, a.alpha, len(classes), seed=r)
            for row in res:
                rows.append({"embedding": emb, "rep": r, **row})
            screens.append({"embedding": emb, "rep": r, **sc})
        print(f"{emb} done", flush=True)
    out = ROOT / "results" / "real" / "r1_gingiva"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "per_rep.csv", index=False)
    pd.DataFrame(screens).to_csv(out / "screens.csv", index=False)
    save_meta(out / "provenance.json", stamp(vars(a)) | {"simulation_only": False, "exploratory": True,
                                                        "data": info, "lock": "results/real/ANALYSIS_LOCK.md"})
    print("written to", out)


if __name__ == "__main__":
    main()
