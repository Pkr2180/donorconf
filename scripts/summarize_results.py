"""Generate results/SUMMARY.md from the experiment CSVs. No number is typed by hand.

    python scripts/summarize_results.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"


def md(df: pd.DataFrame, nd: int = 3) -> str:
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = [" / ".join(str(x) for x in c if str(x)).replace("donor_mean_coverage", "cov")
                      .replace("mean_set_size", "size").replace("p_rep_undercovers_5pts", "p_rep_under")
                      .replace("frac_donors_poorly_covered", "frac_poor").replace("ci_coverage", "cover")
                      .replace("mean_ci_width", "width") for c in df.columns]
    df = df.reset_index()
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(f"{v:.{nd}f}" if isinstance(v, float) else str(v) for v in row) + " |")
    return "\n".join(lines)


def prov(name: str) -> str:
    p = R / name / "provenance.json"
    if not p.exists():
        return "_no provenance file_"
    j = json.loads(p.read_text())
    c = j["config"]
    return f"generated {j['utc_time']} | reps={c.get('reps')} | alpha={c.get('alpha')} | SIMULATION ONLY"


def main():
    out = ["# Development-time simulation results (auto-generated)", "",
           "**All numbers below come from simulated data with known ground truth.** They test whether the "
           "code meets its stated properties and show where methods break. They are not biological evidence "
           "and were seen before the real-data protocol was locked (see PROTOCOL.md).", ""]

    h1 = R / "h1_coverage" / "summary.csv"
    if h1.exists():
        d = pd.read_csv(h1)
        keys = ["donor_mean_coverage", "mean_set_size", "p_rep_undercovers_5pts", "frac_donors_poorly_covered"]
        out += ["## H1 coverage on new donors (target = 1 - alpha)", f"_{prov('h1_coverage')}_", "",
                "Columns: donor_mean_coverage = coverage averaged with equal donor weight; "
                "p_rep_undercovers_5pts = share of replications whose realised coverage fell >5 points below target; "
                "frac_donors_poorly_covered = share of donors covered >10 points below target.", "",
                "### Overall (mean over all scenarios)", md(d.groupby("method")[keys].mean()), ""]
        for f in ("size_difficulty", "n_cal_donors", "donor_sd"):
            out += [f"### By {f}", md(d.groupby([f, "method"])[keys].mean().unstack("method").round(3)), ""]
        rs = d.dropna(subset=["reported_random_split"])
        out += ["### What a within-dataset random cell split would report",
                f"reported coverage range across scenarios: {rs.reported_random_split.min():.3f} to "
                f"{rs.reported_random_split.max():.3f}; realised (new-donor) coverage of the same sets ranged "
                f"{d[d.method == 'cell_pooled'].donor_mean_coverage.min():.3f} to "
                f"{d[d.method == 'cell_pooled'].donor_mean_coverage.max():.3f}.", ""]
        for m in d.method.unique():
            s = d[d.method == m]
            out.append(f"- {m}: worst scenario coverage {s.donor_mean_coverage.min():.3f}; "
                       f"scenarios below target-0.05: {(s.donor_mean_coverage < s.target - 0.05).sum()} of {len(s)}")
        out.append("")

    h2 = R / "h2_shift" / "summary.csv"
    if h2.exists():
        d = pd.read_csv(h2)
        out += ["## H2 label-free donor-level shift screen", f"_{prov('h2_shift')}_", "",
                "p_fail = probability that realised donor-mean coverage fell >5 points below target. "
                "A useful screen has p_fail_given_flagged clearly above p_fail_given_not_flagged.", "",
                md(d.set_index(["shift", "magnitude"]).round(3)), ""]

    h3 = R / "h3_composition" / "summary.csv"
    if h3.exists():
        d = pd.read_csv(h3)
        out += ["## H3 composition-difference intervals (nominal 0.95)", f"_{prov('h3_composition')}_", "",
                "confusion_proxy is a confusion-matrix correction standing in for the idea behind DCATS; "
                "DCATS itself has NOT been run.", "",
                "### Overall", md(d.groupby("method")[["ci_coverage", "mean_ci_width"]].mean()), ""]
        for f in ("donors_per_group", "labeled_per_group", "donor_sd"):
            out += [f"### By {f}",
                    md(d.groupby([f, "method"])[["ci_coverage", "mean_ci_width"]].mean().unstack("method").round(3)), ""]
        out += ["### By cell type (0 and 1 have a true difference; 2 is a true null)",
                md(d.groupby(["cell_type", "method"]).ci_coverage.mean().unstack("method").round(3)), ""]

    geo = R / "geo_check.json"
    if geo.exists():
        out += ["## GEO record check (scripts/verify_geo_accessions.py)", ""]
        for r in json.loads(geo.read_text()):
            out.append(f"- {r['accession']}: found={r.get('found')} samples={r.get('n_samples')} "
                       f"taxon={r.get('taxon')} | {(r.get('title') or '')[:80]}")
        out.append("")
    (R / "SUMMARY.md").write_text("\n".join(out), encoding="utf-8")
    print("wrote", R / "SUMMARY.md")


if __name__ == "__main__":
    main()
