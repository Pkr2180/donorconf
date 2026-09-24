"""Generate results/real/SUMMARY_REAL.md from the real-data CSVs. No number is typed by hand."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1] / "results" / "real"


def md(df: pd.DataFrame, nd: int = 3) -> str:
    df = df.reset_index()
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(f"{v:.{nd}f}" if isinstance(v, float) else str(v) for v in row) + " |")
    return "\n".join(lines)


def coverage_block(name: str, title: str, alpha: float = 0.10) -> list[str]:
    d = R / name / "per_rep.csv"
    if not d.exists():
        return [f"## {title}", "_not run yet_", ""]
    df = pd.read_csv(d)
    prov = json.loads((R / name / "provenance.json").read_text())
    tgt = 1 - alpha
    df["under5"] = df.donor_mean_coverage < tgt - 0.05
    g = df.groupby(["embedding", "method"]).agg(
        reps=("rep", "nunique"), coverage=("donor_mean_coverage", "mean"), cov_sd=("donor_mean_coverage", "std"),
        set_size=("mean_set_size", "mean"), p_rep_under5=("under5", "mean"),
        frac_poor_donors=("frac_donors_poorly_covered", "mean"),
        reported_random_split=("reported_random_split_coverage", "mean"))
    out = [f"## {title}", f"_data: {json.dumps({k: v for k, v in prov.get('data', {}).items() if k != 'classes'}, default=str)[:400]}_",
           f"_classes: {prov.get('data', {}).get('classes')}_", "",
           f"Target coverage {tgt:.2f}. coverage = donor-mean coverage on TARGET donors. p_rep_under5 = share of "
           f"replications more than 5 points below target.", "", md(g), ""]
    cp = df[df.method == "cell_pooled"]
    out.append("Reported (within-calibration random split) minus realised (target donors) coverage for cell_pooled, by embedding:")
    gap = (cp.groupby("embedding").reported_random_split_coverage.mean() - cp.groupby("embedding").donor_mean_coverage.mean())
    out += [md(gap.rename("reported_minus_realised").to_frame()), ""]
    s = R / name / "screens.csv"
    if s.exists():
        sc = pd.read_csv(s)
        rows = []
        for meth in ("cell_pooled", "donor_crc"):
            m = df[df.method == meth][["embedding", "rep", "donor_mean_coverage"]].merge(sc, on=["embedding", "rep"])
            m["fail"] = m.donor_mean_coverage < tgt - 0.05
            for emb, g in m.groupby("embedding"):
                fl = g.screen_flagged.astype(bool)
                rows.append({"embedding": emb, "method_failing": meth, "screen_flag_rate": fl.mean(),
                             "p_fail": g.fail.mean(),
                             "p_fail_given_flagged": g.fail[fl].mean() if fl.any() else np.nan,
                             "p_fail_given_not_flagged": g.fail[~fl].mean() if (~fl).any() else np.nan,
                             "head_accuracy_target": g.head_accuracy_target.mean()})
        out += ["Shift screen vs coverage failure (fail = donor-mean coverage >5 points below target). "
                "A useful screen has p_fail_given_flagged well above p_fail_given_not_flagged:",
                md(pd.DataFrame(rows).set_index(["embedding", "method_failing"])), ""]
    return out


def main():
    out = ["# Real-data results, CELLxGENE Census 2025-11-08 (auto-generated)", "",
           "**Exploratory.** Analysis parameters were fixed in `ANALYSIS_LOCK.md` before results, but the protocol is not "
           "externally registered. Embeddings are hosted by the Census; TranscriptFormer and scVI were trained on Census cells "
           "that probably include these donors, so evaluation is not independent of pretraining. Labels are the Census "
           "`cell_type` (curated mapping of the authors' labels), not model output.", ""]
    out += coverage_block("r1_gingiva", "R1 gingiva: random donor splits (12 ref / 12 cal / 10 target donors)")
    out += coverage_block("r2_blood", "R2 blood: leave-dataset-out (30 donors per role, 7 coarse types)")
    out += coverage_block("r3_oral", "R3 oral: gingiva -> oropharynx -> oral SCC")
    out += coverage_block("r3b_oral", "R3b oral: gingiva -> multi-site normal oral -> oral SCC")
    p = R / "r4_composition" / "per_rep.csv"
    out.append("## R4 gingiva composition: normal vs periodontitis (nominal 0.95)")
    if p.exists():
        d = pd.read_csv(p)
        out += ["Reference = difference in mean per-donor proportions from all donors' Census labels. Cross-fitted predictions. "
                "confusion_proxy is NOT DCATS.", "",
                "### Overall", md(d.groupby(["embedding", "method"]).agg(coverage=("covered", "mean"), width=("width", "mean"))), ""]
        top = d[d.rep == 0].assign(a=lambda x: x.truth.abs()).drop_duplicates(["embedding", "cell_type"]).sort_values("a", ascending=False)
        big = set(top.head(5).cell_type)
        out += [f"### Five cell types with the largest reference difference: {sorted(big)}",
                md(d[d.cell_type.isin(big)].groupby(["embedding", "method"]).agg(coverage=("covered", "mean"), width=("width", "mean"))), ""]
    else:
        out.append("_not run yet_")
    (R / "SUMMARY_REAL.md").write_text("\n".join(out), encoding="utf-8")
    print("wrote", R / "SUMMARY_REAL.md")


if __name__ == "__main__":
    main()
