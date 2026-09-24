# donorconf

Donor-level conformal calibration and donor-level composition inference for transferred
single-cell annotations. **Pre-release research code (v0.0.1).**

Working paper title: *Donor-Level Conformal Calibration for Single-Cell Foundation Models Under Dataset Shift.*

## Status — read this first
* **Real data are retrieved and analysed (exploratory).** Real datasets from CELLxGENE Census 2025-11-08
  (gingiva, other oral sites, oral and oropharyngeal cancer, blood, lung) with hosted TranscriptFormer and scVI
  embeddings were pulled on Modal (`modal/modal_census.py`, workspace `pradeepaiperio`). Results:
  `results/real/SUMMARY_REAL.md`, `results/real/SUMMARY_FINAL.md`. Parameters were fixed in
  `results/real/ANALYSIS_LOCK.md` first, but the protocol is **not** externally registered.
* **The real data narrowed the claim again.** Within one dataset with random donor splits, pooled-cell
  conformal was already at nominal coverage, so "cells are not exchangeable" is not shown to matter there.
  Miscoverage appeared under **dataset shift**, in either direction; see `PROTOCOL.md` §1.
* Embedding extraction is not implemented (the Census hosts the embeddings used). Foundation-model pretraining
  probably included these cells, so evaluation is not independent of pretraining.
* **Framings B and C were then run** (`results/real/SUMMARY_BC.md`, figures in `results/real/figures/`): a
  target-adaptive calibration using a few labeled target donors (B) and a 24-triple blood benchmark of reported vs
  realised coverage with a label-free early-warning signal (C). Both exploratory; prior art for B is substantial
  (see WORKPLAN.md).
* **A second tissue (lung, 6 datasets) replicates C.** Reported-vs-realised gap is larger in lung than in blood
  (71.5% vs 47.9% of replications with |gap| > 0.03); the confidence-shift early warning replicates (Spearman
  −0.776 in lung vs −0.732 in blood, both pass the |Spearman| ≥ 0.4 kill rule); the permutation screen again does
  not predict shortfall. See `results/real/SUMMARY_FINAL.md`.
* **Published-package baseline run.** The conformalized single-cell annotator's own torchCP 1.0.2 calibration
  (standard/THR, standard/APS, classwise, clustered) was run on this study's head/splits
  (`src/donorconf/baseline_torchcp.py`) — not a re-implementation. Its standard-THR predictor is numerically
  identical to this study's cell-pooled method (cross-checked, <0.002 disagreement); its classwise/clustered
  variants under-cover the same way this study's own classwise method does. Does NOT include the annotator's OOD
  detector or neural classifier.
* The package is tested (32 tests). Simulation results (`results/SUMMARY.md`) are exploratory too.
* **Manuscript evidence pack:** `manuscript/BiB_SKELETON.md` (auto-generated from result files by
  `scripts/build_bib_pack.py`; prose sections are marked **[AUTHOR TO WRITE]** because the journal states LLM
  drafting is unacceptable). Literature check (web search only, PubMed connector not authorised):
  `manuscript/LITERATURE_CHECK.md`.

## The idea
Cells from one donor are not independent, so a guarantee that treats cells as exchangeable does not
describe what happens on a new donor. `donorconf` makes the **donor** the unit:

| Component | What it gives | Caveat |
|---|---|---|
| `DonorLevelCRC` | Finite-sample bound on expected per-donor miscoverage (equal donor weights) | Conservative with few donors; vacuous when α < 1/(n+1) |
| `DonorWeightedQuantile` | Efficient, removes size–difficulty bias | Plug-in estimate, **no** finite-sample guarantee |
| `CellLevelConformal` | Baseline for comparison | Cell exchangeability is not satisfied across donors |
| `shift_screen` | Label-free donor-level permutation test of donor summaries | Cannot see concept shift; "not flagged" is not a certificate |
| `TargetAdaptiveCRC` | Weighted donor CRC using k labeled target donors, weight chosen by leave-one-target-donor-out | Guarantee only under the true density-ratio weights, which are unknown; needs k labeled donors |
| `ClasswiseConformal` | Class-conditional baseline (stand-in for the published classwise option) | Not the published package |
| `ppi_group_diff` | Donor-level prediction-powered interval for composition differences | Needs gold labels on a random subset of donors |

## Install and test
```bash
pip install -e .            # or: pip install -r requirements.txt
python -m pytest -q         # 32 tests
```

## Use (real data you have obtained and reviewed)
```bash
# 1. probabilities from a head trained on reference donors
python scripts/prepare_inputs.py --input atlas.h5ad --embedding X_scGPT --reference DS_A \
    --calibration DS_B --target DS_C --out data/run1
# 2. calibrate, screen, report
donorconf calibrate --input data/run1/calibration.npz --alpha 0.1 --out data/run1/model.json
donorconf screen    --cal data/run1/calibration.npz --target data/run1/target.npz --out data/run1/screen.json
donorconf report    --model data/run1/model.json --input data/run1/target.npz --out data/run1/report.json
```
Python API:
```python
from donorconf import DonorLevelCRC, DonorWeightedQuantile
m = DonorLevelCRC(alpha=0.1).fit(probs_cal, labels_cal, donors_cal)
sets = m.predict_sets(probs_new)      # boolean (cells x classes)
```

## Reproduce the simulations
```bash
cd experiments
python run_h1_coverage.py --reps 100      # new-donor coverage, 27 scenarios
python run_h2_shift.py    --reps 150      # shift screen
python run_h3_composition.py --reps 150   # composition intervals
cd .. && python scripts/summarize_results.py   # writes results/SUMMARY.md from the CSVs
```

## Layout
```
src/donorconf/     package (scores, conformal, adaptive, shift, ppi, metrics, simulate, head, lineage, cli,
                    baseline_torchcp: published-package torchCP calibration baseline)
modal/             modal_census.py: survey / pull real Census data on Modal (pradeepaiperio); --tissue-general filter
experiments/       simulation kill tests (run_h*.py) and real-data runs (real_r1..r4, real_crossdataset,
                    real_benchmark.py --torchcp for C, real_common.py)
scripts/           verify_geo_accessions, prepare_inputs, summarize_results, summarize_real, summarize_bc,
                    summarize_final (blood tcp + lung), build_bib_pack (manuscript skeleton), audit_lung_fghi,
                    pull_lung*.sh, retry_pull.sh
config/            datasets_manifest.csv, lineage_audit.csv (blood/oral), lineage_audit_lung.csv (lung)
tests/             32 tests (incl. tests/test_baseline_lung.py: torchCP cross-check + lung_coarse rules)
results/           SUMMARY.md (simulation), real/SUMMARY_REAL.md, real/SUMMARY_BC.md, real/SUMMARY_FINAL.md
                    (published-package baseline + lung), real/ANALYSIS_LOCK.md (addenda A-C3)
manuscript/        BiB_SKELETON.md (auto-generated evidence pack), LITERATURE_CHECK.md
PROTOCOL.md        hypotheses, kill criteria, data rules (draft, not registered)
WORKPLAN.md        phases, decisions for the author, current bottleneck
MANUSCRIPT_OUTLINE.md
HANDOFF_2026-09-21.md   session handoff notes (superseded once lung/tcp work below is read)
```

## Reproduce the real-data runs
```bash
set MODAL_PROFILE=pradeepaiperio
python -m modal run modal/modal_census.py::survey_main          # donor structure of oral / blood datasets
python -m modal run modal/modal_census.py::pull_main --dataset-id <census dataset id> --tag gingiva
python -m modal volume get donorconf-data gingiva.npz results/real/data/
cd experiments && python real_r1_gingiva.py && python real_crossdataset.py --help
cd .. && python scripts/summarize_real.py
```
Run the heavy real-data jobs one at a time: three in parallel exhausted 16 GB of RAM.

## Limits
Guarantees are expectations over donors, assume target donors are exchangeable with calibration donors,
and cover only the calibration label space. Not for clinical use.
