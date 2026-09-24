"""Generate results/real/SUMMARY_BC.md and figures for framings B and C. No number is typed by hand."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parents[1] / "results" / "real"
FIG = R / "figures"
EMB = {"emb_tf-sapiens": "TF-Sapiens", "emb_tf-exemplar-human": "TF-Exemplar", "emb_scvi": "scVI"}
INK, MUTED = "#0b0b0b", "#52514e"
PAL = {"TF-Sapiens": "#2a78d6", "TF-Exemplar": "#eb6834", "scVI": "#1baf7a"}   # validated slots 1-3, all-pairs
MARK = {"TF-Sapiens": "o", "TF-Exemplar": "s", "scVI": "^"}                     # secondary encoding (contrast relief)


def md(df: pd.DataFrame, nd: int = 3) -> str:
    df = df.reset_index()
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(f"{v:.{nd}f}" if isinstance(v, float) else str(v) for v in row) + " |")
    return "\n".join(lines)


def b_block(name, title, alpha=0.10):
    p = R / name / "per_rep.csv"
    if not p.exists():
        return [f"### {title}", "_not run yet_", ""]
    d = pd.read_csv(p)
    d["under5"] = d.donor_mean_coverage < 1 - alpha - 0.05
    g = d.groupby(["k", "method"]).agg(reps=("rep", "nunique"), coverage=("donor_mean_coverage", "mean"),
                                       cov_sd=("donor_mean_coverage", "std"), set_size=("mean_set_size", "mean"),
                                       p_rep_under5=("under5", "mean"), mean_omega=("omega", "mean"))
    return [f"### {title}",
            f"Target {1 - alpha:.2f}; scored on target donors NOT used for adaptation; pooled over embeddings.", "", md(g), ""]


def c_block(alpha=0.10):
    p, s = R / "c_benchmark_blood" / "per_rep.csv", R / "c_benchmark_blood" / "screens.csv"
    if not p.exists():
        return ["## C benchmark", "_not run yet_", ""], None
    d, sc = pd.read_csv(p), pd.read_csv(s)
    d["emb"] = d.embedding.map(EMB)
    d["under5"] = d.donor_mean_coverage < 1 - alpha - 0.05
    out = ["## C: blood benchmark (reported vs realised coverage under dataset shift)",
           f"Ordered dataset triples: {d.groupby(['reference', 'calibration', 'target']).ngroups}; "
           f"rows: {len(d)}. Target coverage {1 - alpha:.2f}. Exploratory.", "",
           "### Coverage and set size by method (mean over triples, replications, embeddings)",
           md(d.groupby("method").agg(coverage=("donor_mean_coverage", "mean"), cov_sd=("donor_mean_coverage", "std"),
                                      set_size=("mean_set_size", "mean"), p_rep_under5=("under5", "mean"),
                                      frac_poor_donors=("frac_donors_poorly_covered", "mean"))), "",
           "### By embedding",
           md(d.groupby(["emb", "method"]).agg(coverage=("donor_mean_coverage", "mean"),
                                               set_size=("mean_set_size", "mean"), p_rep_under5=("under5", "mean"))), ""]
    cp = d[d.method == "cell_pooled"].merge(sc, on=["embedding", "triple", "rep"], suffixes=("", "_s"))
    cp["gap"] = cp.reported_random_split_coverage - cp.donor_mean_coverage
    cp["shortfall"] = (1 - alpha) - cp.donor_mean_coverage
    out += ["### Reported minus realised coverage (cell_pooled)",
            f"mean {cp.gap.mean():.3f}, sd {cp.gap.std():.3f}; share with |gap| > 0.03: {(cp.gap.abs() > 0.03).mean():.3f}; "
            f"share over-reported (reported > realised + 0.03): {(cp.gap > 0.03).mean():.3f}; "
            f"share under-reported (reported < realised - 0.03): {(cp.gap < -0.03).mean():.3f}", ""]
    rows = []
    for col in ("conf_shift", "entropy_shift", "screen_p", "head_accuracy_calibration"):
        r, pv = stats.spearmanr(cp[col], cp.shortfall)
        rows.append({"label_free_signal": col, "spearman_vs_shortfall": r, "p_value": pv})
    agg_rows = []
    for label, keys in (("per embedding x triple (n=72)", ["embedding", "triple"]), ("per triple, embeddings averaged (n=24)", ["triple"])):
        g = cp.groupby(keys)[["conf_shift", "entropy_shift", "screen_p", "head_accuracy_calibration", "shortfall"]].mean()
        for col in ("conf_shift", "entropy_shift", "screen_p", "head_accuracy_calibration"):
            r, pv = stats.spearmanr(g[col], g.shortfall)
            agg_rows.append({"unit": label, "label_free_signal": col, "spearman_vs_shortfall": r, "p_value": pv, "n": len(g)})
    fail = (cp.donor_mean_coverage < 1 - alpha - 0.05).astype(int)
    if fail.nunique() == 2:
        from sklearn.metrics import roc_auc_score
        rows.append({"label_free_signal": "screen flagged: AUROC for >5pt failure",
                     "spearman_vs_shortfall": roc_auc_score(fail, cp.screen_flagged.astype(int)), "p_value": np.nan})
    out += ["### Can a label-free signal predict shortfall? (kill rule: |Spearman| >= 0.4)",
            "Correlations are across replications that reuse the same five datasets, so p-values are optimistic.",
            md(pd.DataFrame(rows).set_index("label_free_signal"), 3), "",
            "Same correlations with replications averaged first (the fairer unit; the same five datasets still recur):",
            md(pd.DataFrame(agg_rows).set_index(["unit", "label_free_signal"]), 3), ""]
    return out, (d, cp)


def figures(d, cp, alpha=0.10):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                         "axes.spines.right": False})
    cp = cp.assign(emb=cp.embedding.map(EMB))
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for e, g in cp.groupby("emb"):
        ax.scatter(g.reported_random_split_coverage, g.donor_mean_coverage, s=22, c=PAL[e], marker=MARK[e],
                   edgecolors="#fcfcfb", linewidths=0.6, alpha=0.85, label=e)
    lo = min(cp.donor_mean_coverage.min(), 0.7)
    ax.plot([lo, 1], [lo, 1], color=MUTED, lw=1, ls="--")
    ax.set_xlim(0.85, 0.95)
    ax.axhline(1 - alpha, color="#c3c2b7", lw=0.8)
    ax.axvline(1 - alpha, color="#c3c2b7", lw=0.8)
    ax.set_xlabel("Coverage reported by a within-dataset random split")
    ax.set_ylabel("Realised coverage on target-dataset donors")
    ax.set_title("Reported coverage is ~0.90 whatever the realised coverage", loc="left", fontsize=11)
    ax.legend(frameon=False, title="Embedding", loc="upper left", bbox_to_anchor=(0.0, 0.62))
    fig.tight_layout()
    fig.savefig(FIG / "fig1_reported_vs_realised.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for e, g in cp.groupby("emb"):
        ax.scatter(g.conf_shift, g.shortfall, s=22, c=PAL[e], marker=MARK[e], edgecolors="#fcfcfb",
                   linewidths=0.6, alpha=0.85, label=e)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("Target minus calibration mean confidence (label-free)")
    ax.set_ylabel("Coverage shortfall (target minus realised)")
    ax.set_title("Does a label-free signal predict shortfall?", loc="left", fontsize=11)
    ax.legend(frameon=False, title="Embedding")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_shortfall_vs_confidence_shift.png", dpi=200)
    plt.close(fig)

    d = d.assign(emb=d.embedding.map(EMB))
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for ax, (e, g) in zip(axes, d.groupby("emb")):
        m = g.groupby("method").agg(cov=("donor_mean_coverage", "mean"), size=("mean_set_size", "mean"))
        for k in ("cell_pooled", "classwise", "donor_weighted", "donor_crc"):
            if k in m.index:
                ax.scatter(m.loc[k, "size"], m.loc[k, "cov"], s=46, c=PAL[e], edgecolors="#fcfcfb", linewidths=0.8, zorder=3)
                off = {"cell_pooled": (7, -12), "donor_weighted": (6, 5), "classwise": (6, -3), "donor_crc": (-6, 6)}[k]
                ax.annotate(k, (m.loc[k, "size"], m.loc[k, "cov"]), textcoords="offset points", xytext=off,
                            fontsize=8, color=MUTED, ha="right" if off[0] < 0 else "left")
        ax.axhline(1 - alpha, color="#c3c2b7", lw=0.8)
        ax.set_title(e, loc="left", fontsize=10)
        ax.set_xlabel("Mean prediction-set size")
    axes[0].set_ylabel("Realised coverage (grey line: target)")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_coverage_vs_size.png", dpi=200)
    plt.close(fig)


def main():
    out = ["# Framings B and C, real data (auto-generated)", "",
           "**Exploratory.** Written after R1-R4; see `ANALYSIS_LOCK.md` Addendum B. Embedding and pretraining caveats "
           "from `SUMMARY_REAL.md` apply. Nothing here is a confirmatory result.", "",
           "## B: target-adaptive calibration with k labeled target donors", ""]
    out += b_block("b_blood", "Blood, leave-dataset-out (A, B, C, permuted roles)")
    out += b_block("b_oral", "Oral: gingiva -> oropharynx -> oral SCC")
    out += b_block("b_oral_multisite", "Oral: gingiva -> multi-site oral -> oral SCC")
    c, obj = c_block()
    out += c
    if obj is not None:
        figures(*obj)
        out += ["Figures: `figures/fig1_reported_vs_realised.png`, `fig2_shortfall_vs_confidence_shift.png`, "
                "`fig3_coverage_vs_size.png`.", ""]
    (R / "SUMMARY_BC.md").write_text("\n".join(out), encoding="utf-8")
    print("wrote", R / "SUMMARY_BC.md")


if __name__ == "__main__":
    main()
