"""results/real/SUMMARY_FINAL.md: published-package baselines (blood) and the second tissue (lung). No hand-typed numbers."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from summarize_bc import EMB, md

R = Path(__file__).resolve().parents[1] / "results" / "real"
ALPHA = 0.10


def load(name):
    p = R / name / "per_rep.csv"
    if not p.exists():
        return None, None
    d = pd.read_csv(p)
    d["emb"] = d.embedding.map(EMB)
    d["fail5"] = d.donor_mean_coverage < 1 - ALPHA - 0.05
    return d, pd.read_csv(R / name / "screens.csv")


def method_table(d):
    return md(d.groupby("method").agg(reps=("rep", "size"), coverage=("donor_mean_coverage", "mean"),
                                      cov_sd=("donor_mean_coverage", "std"), set_size=("mean_set_size", "mean"),
                                      p_fail5=("fail5", "mean"), frac_poor_donors=("frac_donors_poorly_covered", "mean")))


def signal_block(d, s):
    cp = d[d.method == "cell_pooled"].merge(s, on=["embedding", "triple", "rep"], suffixes=("", "_s"))
    cp["gap"] = cp.reported_random_split_coverage - cp.donor_mean_coverage
    cp["shortfall"] = (1 - ALPHA) - cp.donor_mean_coverage
    g = cp.groupby("triple")[["conf_shift", "entropy_shift", "screen_p", "head_accuracy_calibration", "shortfall"]].mean()
    rows = []
    for c in ("conf_shift", "entropy_shift", "screen_p", "head_accuracy_calibration"):
        r, p = stats.spearmanr(g[c], g.shortfall)
        rows.append({"signal": c, "spearman_per_triple": r, "p_value": p, "n": len(g)})
    head = (f"reported {cp.reported_random_split_coverage.mean():.3f} vs realised {cp.donor_mean_coverage.mean():.3f}; "
            f"|gap| > 0.03 in {(cp.gap.abs() > 0.03).mean():.1%} of replications "
            f"(over-reported {(cp.gap > 0.03).mean():.1%}, under-reported {(cp.gap < -0.03).mean():.1%})")
    return head, md(pd.DataFrame(rows).set_index("signal"))


def main():
    out = ["# Published-package baselines and second tissue (auto-generated)", "",
           "**Exploratory** (`ANALYSIS_LOCK.md` Addendum C). tcp_* rows are the torchCP 1.0.2 predictors used by the "
           "conformalized single-cell annotator (Bioinformatics 2025, btaf521) on this study's head and splits; the annotator's "
           "OOD detector and neural classifier are not included. Embeddings are Census-hosted and probably saw these cells.", ""]
    for name, title in (("c_benchmark_blood_tcp", "Blood (5 datasets, 24 triples) with published-package baselines"),
                        ("c_benchmark_lung", "Lung (6 datasets: C, E, F, G, H, I), second tissue")):
        d, s = load(name)
        out.append(f"## {title}")
        if d is None:
            out += ["_not run yet_", ""]
            continue
        head, sig = signal_block(d, s)
        out += [f"Rows {len(d)}; triples {d.triple.nunique()}; target coverage {1 - ALPHA:.2f}.", "",
                "### Methods (mean over triples, replications, embeddings)", method_table(d), "",
                "### By embedding (coverage / set size)",
                md(d.groupby(["emb", "method"]).agg(coverage=("donor_mean_coverage", "mean"), set_size=("mean_set_size", "mean"),
                                                    p_fail5=("fail5", "mean"))), "",
                "### Reported vs realised (cell_pooled)", head, "",
                "### Label-free signals vs shortfall (per triple; kill rule |Spearman| >= 0.4)", sig, ""]
        # does the published package differ from this study's pooled-cell baseline?
        tcp = d[d.method == "tcp_standard_thr"][["embedding", "triple", "rep", "donor_mean_coverage"]]
        pooled = d[d.method == "cell_pooled"][["embedding", "triple", "rep", "donor_mean_coverage"]]
        if len(tcp):
            m = tcp.merge(pooled, on=["embedding", "triple", "rep"], suffixes=("_tcp", "_pooled"))
            out += [f"Cross-check: torchCP standard-THR vs this package's pooled-cell coverage, mean absolute difference "
                    f"{(m.donor_mean_coverage_tcp - m.donor_mean_coverage_pooled).abs().mean():.4f} over {len(m)} runs.", ""]
    (R / "SUMMARY_FINAL.md").write_text("\n".join(out), encoding="utf-8")
    print("wrote", R / "SUMMARY_FINAL.md")


if __name__ == "__main__":
    main()
