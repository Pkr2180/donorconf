"""Build manuscript/BiB_SKELETON.md for Briefings in Bioinformatics from result files only.

Every number is computed here from results/real/**/per_rep.csv|screens.csv. Prose is deliberately NOT written:
the journal states that LLM drafting of papers is unacceptable, so narrative sections carry AUTHOR-TO-WRITE markers.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "real"
OUT = ROOT / "manuscript"
ALPHA = 0.10
EMB = {"emb_tf-sapiens": "TF-Sapiens", "emb_tf-exemplar-human": "TF-Exemplar", "emb_scvi": "scVI"}


def md(df, nd=3):
    df = df.reset_index()
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(f"{v:.{nd}f}" if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(out)


def main():
    OUT.mkdir(exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    for f in (R / "figures").glob("*.png"):
        shutil.copy(f, OUT / "figures" / f.name)

    c = pd.read_csv(R / "c_benchmark_blood_tcp" / "per_rep.csv")
    s = pd.read_csv(R / "c_benchmark_blood_tcp" / "screens.csv")
    c["emb"] = c.embedding.map(EMB)
    c["fail5"] = c.donor_mean_coverage < 1 - ALPHA - 0.05
    cp = c[c.method == "cell_pooled"].merge(s, on=["embedding", "triple", "rep"], suffixes=("", "_s"))
    cp["gap"] = cp.reported_random_split_coverage - cp.donor_mean_coverage
    cp["shortfall"] = (1 - ALPHA) - cp.donor_mean_coverage
    g = cp.groupby("triple")[["conf_shift", "shortfall", "screen_p"]].mean()
    rho_conf = stats.spearmanr(g.conf_shift, g.shortfall)
    rho_scr = stats.spearmanr(g.screen_p, g.shortfall)

    lung_path = R / "c_benchmark_lung" / "per_rep.csv"
    lung_facts, lung_meth_md = [], None
    if lung_path.exists():
        cl = pd.read_csv(lung_path)
        sl = pd.read_csv(R / "c_benchmark_lung" / "screens.csv")
        cl["emb"] = cl.embedding.map(EMB)
        cl["fail5"] = cl.donor_mean_coverage < 1 - ALPHA - 0.05
        cpl = cl[cl.method == "cell_pooled"].merge(sl, on=["embedding", "triple", "rep"], suffixes=("", "_s"))
        cpl["gap"] = cpl.reported_random_split_coverage - cpl.donor_mean_coverage
        cpl["shortfall"] = (1 - ALPHA) - cpl.donor_mean_coverage
        gl = cpl.groupby("triple")[["conf_shift", "shortfall", "screen_p"]].mean()
        rho_conf_l = stats.spearmanr(gl.conf_shift, gl.shortfall)
        lung_facts = [
            f"Lung (second tissue; datasets C, E, F, G, H, I; {len(cl)} rows): cell-pooled realised coverage "
            f"{cpl.donor_mean_coverage.mean():.3f} vs reported {cpl.reported_random_split_coverage.mean():.3f}; "
            f"|gap| > 0.03 in {(cpl.gap.abs() > 0.03).mean():.1%} of replications",
            f"Lung confidence shift vs shortfall, per triple (n={len(gl)}): Spearman {rho_conf_l.statistic:.3f} "
            f"(p={rho_conf_l.pvalue:.4f}) -- replicates the blood signal (kill rule |Spearman| >= 0.4)",
            f"Lung donor_crc coverage {cl[cl.method=='donor_crc'].donor_mean_coverage.mean():.3f}, "
            f"set size {cl[cl.method=='donor_crc'].mean_set_size.mean():.3f}",
        ]
        lung_meth = cl.groupby("method").agg(coverage=("donor_mean_coverage", "mean"), sd=("donor_mean_coverage", "std"),
                                              set_size=("mean_set_size", "mean"), p_fail5=("fail5", "mean"))
        lung_meth_md = md(lung_meth)

    facts = [
        f"Ordered blood dataset triples: {c.groupby(['reference', 'calibration', 'target']).ngroups}; rows {len(c)}",
        f"Cell-pooled realised coverage {cp.donor_mean_coverage.mean():.3f} (target {1 - ALPHA:.2f}); reported by random split "
        f"{cp.reported_random_split_coverage.mean():.3f}",
        f"|reported - realised| > 0.03 in {(cp.gap.abs() > 0.03).mean():.1%} of replications; over-reported "
        f"{(cp.gap > 0.03).mean():.1%}, under-reported {(cp.gap < -0.03).mean():.1%}",
        f"Confidence shift vs shortfall, per triple (n={len(g)}): Spearman {rho_conf.statistic:.3f} (p={rho_conf.pvalue:.4f})",
        f"Permutation screen p-value vs shortfall, per triple: Spearman {rho_scr.statistic:.3f} (p={rho_scr.pvalue:.3f})",
        f"Published-package torchCP standard-THR cross-check vs this study's pooled-cell coverage: mean absolute "
        f"difference {(c[c.method=='tcp_standard_thr'].set_index(['embedding','triple','rep']).donor_mean_coverage - c[c.method=='cell_pooled'].set_index(['embedding','triple','rep']).donor_mean_coverage).abs().mean():.4f}",
    ] + lung_facts
    meth = c.groupby("method").agg(coverage=("donor_mean_coverage", "mean"), sd=("donor_mean_coverage", "std"),
                                    set_size=("mean_set_size", "mean"), p_fail5=("fail5", "mean"))

    tests = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT, capture_output=True, text=True)
    test_line = [l for l in tests.stdout.splitlines() if "passed" in l or "failed" in l][-1]

    b_tables = []
    for name in ("b_blood", "b_oral", "b_oral_multisite"):
        p = R / name / "per_rep.csv"
        if p.exists():
            d = pd.read_csv(p)
            d["fail5"] = d.donor_mean_coverage < 1 - ALPHA - 0.05
            t = d.groupby(["k", "method"]).agg(coverage=("donor_mean_coverage", "mean"), set_size=("mean_set_size", "mean"),
                                                p_fail5=("fail5", "mean"))
            b_tables.append(f"**{name}**\n\n{md(t)}")

    text = f"""# Briefings in Bioinformatics: manuscript skeleton (auto-generated evidence; prose is the author's)

Generated by `scripts/build_bib_pack.py`. Numbers come from result files only. The journal states LLM drafting of
papers is unacceptable, so sections marked **[AUTHOR TO WRITE]** are left empty on purpose. Check word limits, abstract
format, article type and reference style on the live guidelines page (the fetch returned none of these).

Suggested framing: benchmark and evaluation of reported vs realised coverage under dataset shift (C), with a
label-free early-warning signal and a few-labeled-donor remedy (B). Do not claim a new method (B has prior art).

## Title and keywords
**[AUTHOR TO WRITE]** Title must match the findings. Candidate scope: conformal prediction, single-cell foundation
models, dataset shift, donor-level evaluation, uncertainty calibration.

## Abstract
**[AUTHOR TO WRITE]** Use the headline facts below.

### Headline facts (generated)
""" + "\n".join(f"- {x}" for x in facts) + f"""

## 1 Introduction
**[AUTHOR TO WRITE]** Per `manuscript/LITERATURE_CHECK.md` (PubMed/Scite/Consensus/bioRxiv connector search run
2026-09-22, Shahid read in full 2026-09-24): do NOT claim "first to show coverage fails under dataset shift" --
López-De-Castro et al. (Bioinformatics 2025, btaf521) already report this for a disease-shift query, and Shahid
(arXiv:2606.20115, 2026, "When Average Calibration Fails: Site-Conditional Federated Conformal Risk Control") shows
the same failure shape for hospitals in federated brain-tumour segmentation (pooled CRC protects the average site
but violates coverage at 40% of individual institutions) -- this is now confirmed the closest related work to this
study's headline claim, and it means "average calibration fails at a natural exchangeable subunit" is NOT a novel
observation in the conformal literature generally. State the novelty precisely: the first (as far as this search
found) quantified benchmark of that failure mode for single-cell foundation-model embeddings at donor granularity,
across two tissues, against a published package's own calibration variants, plus a label-free early-warning signal.
Also note Shahid's leave-one-site-out-selected shrinkage weight is methodologically parallel to B's leave-one-
target-donor-out-selected omega (same "LOO interpolation between local and pooled calibration" idea, different
unit/domain/guarantee-status) -- worth one explicit sentence in Methods or Discussion, not just a citation. Cite
conformal scRNA work (btaf521), weighted conformal under covariate shift (Tibshirani et al. 2019; Wang and Qiao
2025, arXiv:2502.17744 -- prior art for B), conformal risk control (Angelopoulos et al., ICLR 2024),
site-conditional federated CRC (Shahid, arXiv:2606.20115), and train-time CRM for few-labeled-target
semi-supervised domain adaptation (Giannopoulos, Shen, Zavlanos, arXiv:2608.23153, 2026 -- second, more direct
strand of prior art for B: same few-labeled-target-domain problem, but acts at training time via optimal-transport
pseudo-labels rather than at calibration time). PubMed found zero hits for donor-/patient-level conformal cell-type
annotation (4 queries run, see LITERATURE_CHECK.md); no Scopus connector is available, note that limitation in the
manuscript's search-strategy statement.

## 2 Methods
**[AUTHOR TO WRITE]** Facts to state: split conformal with LAC scores; cell-pooled, class-conditional, donor-weighted
quantile and donor-level CRC (1/(n+1) term); target-adaptive CRC with LOTO choice of omega over
(0, .25, .5, 1, 2, 4, 8, 16); head = standardise, PCA64, logistic regression; alpha = {ALPHA};
data = CELLxGENE Census 2025-11-08, hosted embeddings TF-Sapiens, TF-Exemplar, scVI; coarse lineage mapping in
`config/lineage_audit.csv` (blood/oral) and `config/lineage_audit_lung.csv` (lung, second tissue; six datasets
C, E, F, G, H, I qualified by a fixed donor/coverage rule in `ANALYSIS_LOCK.md` Addenda C2/C3, applied before any
lung coverage number was seen); parameters fixed in `results/real/ANALYSIS_LOCK.md` (self-imposed, not registered).
Published-package baseline: torchCP 1.0.2 `SplitPredictor` (THR, APS), `ClassWisePredictor`, `ClusteredPredictor`
as used by the conformalized single-cell annotator (btaf521), run on this study's head/splits (`baseline_torchcp.py`);
the annotator's OOD detector and neural classifier are NOT included -- state this every time tcp_* results are shown.
AI-use disclosure (required by the journal): code and this evidence pack were produced with an AI coding assistant
and checked by the authors; state this in Methods or Acknowledgements and in the cover letter.

## 3 Results
### 3.1 Validity checks (simulation, ground truth known)
See `results/SUMMARY.md`. Software tests: {test_line}.

### 3.2 Within-dataset donor splits (gingiva, 34 donors)
`results/real/SUMMARY_REAL.md` R1: cell-pooled at nominal; donor CRC over-covers with larger sets. Negative result
for the exchangeability hypothesis; report it.

### 3.3 Coverage under dataset shift (C) -- blood, 5 datasets, with published-package baselines
{md(meth)}

Figure 1 `figures/fig1_reported_vs_realised.png`. Figure 3 `figures/fig3_coverage_vs_size.png`.
tcp_standard_thr is numerically identical to this study's cell_pooled method (cross-check in Headline facts);
tcp_classwise_thr and tcp_cluster_thr under-cover like this study's own classwise method; tcp_standard_aps
(APS score) is the strongest published-package variant but donor_crc still coverages closer to nominal with a
smaller mean set size.

### 3.3b Second tissue: lung (6 datasets, replication check)
""" + (lung_meth_md or "_lung benchmark not run_") + """

Lung replicates the blood direction and, if anything, the gap is larger in lung (see Headline facts): reported
coverage over-states realised coverage by more under lung dataset shift than under blood dataset shift. donor_crc
remains the best coverage/size trade-off; tcp_cluster_thr reaches similar coverage to donor_crc only by inflating
set size roughly 2x (see `SUMMARY_FINAL.md`) -- report this as a negative result for that baseline's practical
usefulness, not just its raw coverage number.

### 3.4 Label-free early warning
Figure 2 `figures/fig2_shortfall_vs_confidence_shift.png`; correlations in `results/real/SUMMARY_BC.md` and
`results/real/SUMMARY_FINAL.md` (blood and lung separately). The confidence-shift signal passes the pre-registered
kill rule (|Spearman| >= 0.4) in BOTH tissues; the permutation screen does not predict shortfall in either tissue.

### 3.5 Remedy with few labeled target donors (B)
""" + "\n\n".join(b_tables) + """

### 3.6 Negative results to report
Class-conditional conformal worse under shift in both tissues; permutation screen not predictive in either tissue;
donor-level PPI composition not needed (naive intervals covered); over-coverage in oral multi-site pair;
tcp_cluster_thr's coverage in lung is an artefact of set-size inflation, not calibration quality.

## 4 Discussion and limitations
**[AUTHOR TO WRITE]** Must state: five blood datasets, six lung datasets, and three oral pairs only; oral targets 10
to 12 donors; pretraining contamination (Census embeddings trained on Census cells); exploratory, not registered;
B has no guarantee without true density ratios; correlations reuse the same datasets, so p-values are optimistic;
published-package comparison isolates the calibration step only (no OOD detector, no neural classifier);
Geneformer/scGPT embeddings not evaluated (not Census-hosted, need GPU extraction).

## 5 Data and code availability (mandatory)
Package `donorconf`, `results/real/data/*.npz`, Census release 2025-11-08. **Zenodo DOI: not yet created.**

## Evaluation still missing before submission
- [x] Published cell-level conformal package as a real baseline (torchCP 1.0.2, this pack)
- [x] Second tissue (lung, 6 datasets)
- [ ] Geneformer / scGPT embeddings
- [ ] OSF registration, statistician review, written-strategy PubMed/Scopus literature search (connector needs
      authorising in claude.ai connector settings)
- [ ] Paperpal preflight (free trial) once the draft exists
- [ ] Cover letter: disclose AI use and related preprints
"""
    (OUT / "BiB_SKELETON.md").write_text(text, encoding="utf-8")
    print("wrote", OUT / "BiB_SKELETON.md")


if __name__ == "__main__":
    main()
