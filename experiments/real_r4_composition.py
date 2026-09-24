"""R4: real gingiva atlas, normal (20 donors) vs periodontitis (10 donors) composition differences.

Ground reference = the difference in mean per-donor proportions computed from ALL donors' author-derived
Census labels (a finite-population quantity). Predictions for every donor come from a head trained only on
other donors (5-fold donor cross-fitting), so no donor is scored by a model that saw it. Each replication
draws a random labeled subset of donors per group; the rest keep predicted labels only.

Design-based check, not an inference about periodontitis biology: the finite-population reference makes
intervals mildly conservative (no finite-population correction), and 30 donors is small. Exploratory.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from _common import ROOT, stamp
from donorconf.head import EmbeddingHead
from donorconf.ppi import (cell_pooled_diff, confusion_corrected_diff, donor_proportions,
                           naive_plugin_diff, ppi_group_diff)
from real_common import choose_classes, encode, load, save_meta


def crossfit_predictions(X, y, donors, n_classes, folds=5, seed=0):
    ud = np.unique(donors)
    rng = np.random.default_rng(seed)
    order = rng.permutation(ud)
    fold_of = {d: i % folds for i, d in enumerate(order)}
    f = np.array([fold_of[d] for d in donors])
    pred = np.zeros(len(y), dtype=int)
    for k in range(folds):
        tr, te = f != k, f == k
        head = EmbeddingHead(n_components=64, seed=seed).fit(X[tr], y[tr], n_classes)
        pred[te] = head.predict_proba(X[te]).argmax(axis=1)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=150)
    ap.add_argument("--n-lab-normal", type=int, default=6)
    ap.add_argument("--n-lab-perio", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=60)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--embeddings", nargs="*", default=["emb_scvi", "emb_tf-sapiens"])
    a = ap.parse_args()

    z = load("gingiva")
    labels, donors, disease = z["labels"], z["donors"], z["disease"]
    classes = choose_classes(labels, donors, 150, 12)
    y = encode(labels, classes)
    keep = (y >= 0) & np.isin(disease, ["normal", "periodontitis"])
    K = len(classes)
    dis_of = {d: disease[donors == d][0] for d in np.unique(donors[keep])}
    rows = []
    for emb in a.embeddings:
        X = z[emb][keep]
        yy, dd = y[keep], donors[keep]
        pred = crossfit_predictions(X, yy, dd, K)
        ids, P_true = donor_proportions(yy, dd, K)
        _, P_pred = donor_proportions(pred, dd, K, donor_ids=ids)
        grp = np.array([dis_of[d] for d in ids])
        g1, g0 = np.where(grp == "periodontitis")[0], np.where(grp == "normal")[0]
        truth = P_true[g1].mean(axis=0) - P_true[g0].mean(axis=0)
        acc = float((pred == yy).mean())
        print(f"{emb}: cross-fit accuracy {acc:.3f}; donors normal={len(g0)} perio={len(g1)}", flush=True)
        cell_pred_by_donor = {d: pred[dd == d] for d in ids}
        for r in range(a.reps):
            rng = np.random.default_rng(1000 + r)
            lab1 = rng.choice(g1, a.n_lab_perio, replace=False)
            lab0 = rng.choice(g0, a.n_lab_normal, replace=False)
            unl1, unl0 = np.setdiff1d(g1, lab1), np.setdiff1d(g0, lab0)
            lab_all = np.concatenate([lab1, lab0])
            m_lab = np.isin(dd, ids[lab_all])
            pooled1 = np.concatenate([cell_pred_by_donor[d] for d in ids[g1]])
            pooled0 = np.concatenate([cell_pred_by_donor[d] for d in ids[g0]])
            for k in range(K):
                res = {
                    "naive_plugin": naive_plugin_diff(P_pred[g1, k], P_pred[g0, k], a.alpha),
                    "cell_pooled": cell_pooled_diff(pooled1, pooled0, k, a.alpha),
                    "confusion_proxy": confusion_corrected_diff(
                        yy[m_lab], pred[m_lab], dd[m_lab], P_pred[unl1], P_pred[unl0], k, K, a.alpha,
                        a.n_boot, seed=int(rng.integers(1e9))),
                    "ppi_donor": ppi_group_diff(
                        {"y_lab": P_true[lab1, k], "f_lab": P_pred[lab1, k], "f_unlab": P_pred[unl1, k]},
                        {"y_lab": P_true[lab0, k], "f_lab": P_pred[lab0, k], "f_unlab": P_pred[unl0, k]}, a.alpha),
                }
                for m, rr in res.items():
                    rows.append({"embedding": emb, "rep": r, "cell_type": classes[k], "method": m,
                                 "truth": truth[k], "estimate": rr["estimate"], "lo": rr["ci"][0], "hi": rr["ci"][1],
                                 "covered": rr["ci"][0] <= truth[k] <= rr["ci"][1], "width": rr["ci"][1] - rr["ci"][0]})
        print(f"{emb} done", flush=True)
    out = ROOT / "results" / "real" / "r4_composition"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "per_rep.csv", index=False)
    save_meta(out / "provenance.json", stamp(vars(a)) | {"simulation_only": False, "exploratory": True,
              "classes": classes.tolist(), "note": "finite-population design-based check; DCATS not run"})
    print("written to", out)


if __name__ == "__main__":
    main()
