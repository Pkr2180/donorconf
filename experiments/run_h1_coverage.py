"""H1 kill test: does cell-pooled conformal misbehave on NEW donors, and does donor-level CRC fix it?

For each scenario a classification head is trained once on reference donors; then
per replication fresh calibration donors and fresh, independent test donors are
drawn. Compared, all at the same nominal alpha:

  cell_pooled     split conformal over pooled calibration cells
  donor_weighted  pooled quantile, cells weighted 1/(donor size); plug-in, no finite-sample guarantee
  donor_crc       donor-level conformal risk control (donor = unit)

and for the cell-pooled sets a "reported" coverage as a within-dataset random
cell split would measure it (held-out cells of the SAME donors) versus the
"realised" coverage on new donors.

Reading the result honestly: pooled-cell calibration is approximately correct
on average when donor size is unrelated to difficulty (size_difficulty = 0);
the harm is expected where donor size and difficulty are related, and in the
gap between reported and realised coverage. If neither appears, the hypothesis
is not supported.
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd

from _common import write_results
from donorconf import CellLevelConformal, DonorLevelCRC, DonorWeightedQuantile
from donorconf.metrics import coverage_summary
from donorconf.simulate import SimConfig, fit_head, head_probs, make_cohort


def one_scenario(donor_sd, size_difficulty, n_cal, alpha, reps, n_test, seed):
    cfg = SimConfig(donor_sd=donor_sd, size_difficulty=size_difficulty)
    rng = np.random.default_rng(seed)
    head = fit_head(make_cohort(cfg, 30, rng))
    K = cfg.n_classes
    rows = []
    for r in range(reps):
        cal = make_cohort(cfg, n_cal, rng, first_id=0)
        te = make_cohort(cfg, n_test, rng, first_id=10_000)
        Pc, Pt = head_probs(head, cal.X, K), head_probs(head, te.X, K)

        # within-dataset random cell split on the calibration donors (what is usually reported)
        perm = rng.permutation(len(cal.y))
        half = len(perm) // 2
        a, b = perm[:half], perm[half:]
        leak = CellLevelConformal(alpha).fit(Pc[a], cal.y[a])
        s_b = leak.predict_sets(Pc[b])
        reported = float(s_b[np.arange(len(b)), cal.y[b]].mean())

        for name, model in (("cell_pooled", CellLevelConformal(alpha)), ("donor_weighted", DonorWeightedQuantile(alpha)),
                            ("donor_crc", DonorLevelCRC(alpha))):
            model.fit(Pc, cal.y, cal.donor)
            cs = coverage_summary(model.predict_sets(Pt), te.y, te.donor, alpha, tol=0.10)
            rows.append({"method": name, "rep": r, **cs,
                         "reported_random_split_coverage": reported if name == "cell_pooled" else np.nan})
    df = pd.DataFrame(rows)
    df["realised_below_target_by_5pts"] = df["donor_mean_coverage"] < (1 - alpha - 0.05)
    g = df.groupby("method")
    out = g.agg(donor_mean_coverage=("donor_mean_coverage", "mean"),
                donor_mean_coverage_sd=("donor_mean_coverage", "std"),
                cell_pooled_coverage=("cell_pooled_coverage", "mean"),
                frac_donors_poorly_covered=("frac_donors_poorly_covered", "mean"),
                mean_set_size=("mean_set_size", "mean"),
                p_rep_undercovers_5pts=("realised_below_target_by_5pts", "mean"),
                reported_random_split=("reported_random_split_coverage", "mean")).reset_index()
    out.insert(0, "n_cal_donors", n_cal)
    out.insert(0, "size_difficulty", size_difficulty)
    out.insert(0, "donor_sd", donor_sd)
    out.insert(0, "target", 1 - alpha)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--n-test", type=int, default=30)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--quick", action="store_true", help="tiny grid for smoke testing")
    a = ap.parse_args()
    grid = dict(donor_sd=[0.4, 0.8, 1.2], size_difficulty=[-0.5, 0.0, 0.5], n_cal=[6, 12, 24])
    if a.quick:
        grid = dict(donor_sd=[0.8], size_difficulty=[0.0], n_cal=[12])
    frames = []
    for i, (sd, sz, n) in enumerate(itertools.product(*grid.values())):
        frames.append(one_scenario(sd, sz, n, a.alpha, a.reps, a.n_test, a.seed + i))
        print(f"scenario {i + 1} done: donor_sd={sd} size_difficulty={sz} n_cal={n}", flush=True)
    df = pd.concat(frames, ignore_index=True)
    out = write_results("h1_coverage", df, vars(a) | {"grid": grid})
    print(df.round(3).to_string(index=False))
    print("written to", out)


if __name__ == "__main__":
    main()
