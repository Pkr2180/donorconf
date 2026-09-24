"""H3 kill test: do donor-level PPI intervals for a composition difference cover the truth?

Two donor groups differ in true mean composition by construction (class 0 up,
class 1 down, other classes equal). Cell types are labelled by a classification
head trained on reference donors; donor effects make its predictions biased.
Target: mean donor-level proportion difference, whose truth is the prior
difference. A random subset of donors per group has gold labels.

Methods: naive plug-in, cell-pooled Wald (donors ignored), confusion-matrix
correction (a DCATS-style PROXY, not DCATS), donor-level PPI++.
Classes 0 and 1 have a true difference; class 2 is a true null.
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd

from _common import write_results
from donorconf.ppi import (cell_pooled_diff, confusion_corrected_diff, donor_proportions,
                           naive_plugin_diff, ppi_group_diff)
from donorconf.simulate import SimConfig, fit_head, head_probs, make_cohort

EVAL_CLASSES = (0, 1, 2)
DELTA = 0.08


def one_scenario(donor_sd, n_donors, n_lab, alpha, reps, n_boot, seed):
    cfg = SimConfig(donor_sd=donor_sd, comp_conc=30.0)
    K = cfg.n_classes
    rng = np.random.default_rng(seed)
    head = fit_head(make_cohort(cfg, 30, rng))
    p0 = cfg.base_prior()
    p1 = p0.copy()
    p1[0] += DELTA
    p1[1] -= DELTA
    truth = {k: p1[k] - p0[k] for k in EVAL_CLASSES}
    hits = {(m, k): [] for m in ("naive_plugin", "cell_pooled", "confusion_proxy", "ppi_donor")
            for k in EVAL_CLASSES}
    width = {key: [] for key in hits}
    for _ in range(reps):
        groups = {}
        for g, pr, fid in ((1, p1, 0), (0, p0, 5000)):
            co = make_cohort(cfg, n_donors, rng, prior=pr, first_id=fid)
            pred = head_probs(head, co.X, K).argmax(axis=1)
            ids, Pf = donor_proportions(pred, co.donor, K)
            _, Py = donor_proportions(co.y, co.donor, K, donor_ids=ids)
            lab = rng.permutation(ids.size)[:n_lab]
            unl = np.setdiff1d(np.arange(ids.size), lab)
            groups[g] = dict(co=co, pred=pred, ids=ids, Pf=Pf, Py=Py, lab=lab, unl=unl)
        g1, g0 = groups[1], groups[0]
        lab_cells = np.concatenate([
            np.nonzero(np.isin(g["co"].donor, g["ids"][g["lab"]]))[0] + off
            for g, off in ((g1, 0), (g0, len(g1["pred"])))])
        all_true = np.concatenate([g1["co"].y, g0["co"].y])
        all_pred = np.concatenate([g1["pred"], g0["pred"]])
        all_donor = np.concatenate([g1["co"].donor, g0["co"].donor])
        for k in EVAL_CLASSES:
            res = {
                "naive_plugin": naive_plugin_diff(g1["Pf"][:, k], g0["Pf"][:, k], alpha),
                "cell_pooled": cell_pooled_diff(g1["pred"], g0["pred"], k, alpha),
                "confusion_proxy": confusion_corrected_diff(
                    all_true[lab_cells], all_pred[lab_cells], all_donor[lab_cells],
                    g1["Pf"][g1["unl"]], g0["Pf"][g0["unl"]], k, K, alpha, n_boot,
                    seed=int(rng.integers(1e9))),
                "ppi_donor": ppi_group_diff(
                    {"y_lab": g1["Py"][g1["lab"], k], "f_lab": g1["Pf"][g1["lab"], k], "f_unlab": g1["Pf"][g1["unl"], k]},
                    {"y_lab": g0["Py"][g0["lab"], k], "f_lab": g0["Pf"][g0["lab"], k], "f_unlab": g0["Pf"][g0["unl"], k]},
                    alpha),
            }
            for m, r in res.items():
                hits[(m, k)].append(r["ci"][0] <= truth[k] <= r["ci"][1])
                width[(m, k)].append(r["ci"][1] - r["ci"][0])
    rows = []
    for (m, k), h in hits.items():
        rows.append({"donor_sd": donor_sd, "donors_per_group": n_donors, "labeled_per_group": n_lab,
                     "cell_type": k, "true_diff": round(truth[k], 3), "method": m,
                     "ci_coverage": float(np.mean(h)), "mean_ci_width": float(np.mean(width[(m, k)])),
                     "nominal": 1 - alpha})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=150)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--n-boot", type=int, default=100)
    ap.add_argument("--seed", type=int, default=2028)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    grid = dict(donor_sd=[0.5, 1.0], n_donors=[12, 24, 48], n_lab=[4, 8])
    if a.quick:
        grid = dict(donor_sd=[1.0], n_donors=[24], n_lab=[8])
    frames = []
    for i, (sd, n, nl) in enumerate(itertools.product(*grid.values())):
        frames.append(one_scenario(sd, n, nl, a.alpha, a.reps, a.n_boot, a.seed + i))
        print(f"scenario {i + 1} done: donor_sd={sd} donors/group={n} labeled={nl}", flush=True)
    df = pd.concat(frames, ignore_index=True)
    out = write_results("h3_composition", df, vars(a) | {"grid": grid})
    print(df.round(3).to_string(index=False))
    print("written to", out)


if __name__ == "__main__":
    main()
