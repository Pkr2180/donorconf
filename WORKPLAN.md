# Work plan

Legend: DONE = built and tested here (simulation only). TODO = needs your decision or real data.

## Phase 0 — Package and simulation kill tests  (DONE)
* `src/donorconf/`: scores, cell-level baseline, donor-level CRC, donor-weighted quantile, donor-level shift
  screen, donor-level PPI++ and two baselines, hierarchical simulator, CLI. 23 tests pass, including a
  Monte-Carlo check of the finite-sample guarantee and a tightness check.
* `experiments/`: H1, H2, H3 simulations; `scripts/summarize_results.py` generates `results/SUMMARY.md`.

## Phase 1 — Decide the framing  (TODO — needs the author)
The simulations narrowed the claim. Pick one before real data:

| Option | Contribution | Risk |
|---|---|---|
| A. Validity/efficiency trade-off | Donor-level CRC (guaranteed) vs donor-weighted quantile (efficient, no guarantee), reported-vs-realised gap, donor-level composition inference | Moderate; honest and defensible |
| B. Add a middle method | A donor-bootstrap upper-confidence threshold giving a probability statement about realised coverage with smaller sets than CRC | Higher novelty, needs a proof/simulation of its own |
| C. Benchmark-only | Leave-donor-out and leave-dataset-out coverage benchmark across foundation models | Lowest risk, less methodological novelty |

Recommendation at the time: A. **Superseded by the real-data run** (see "Status after real-data run 1" below): A alone is
no longer sufficient.

## Phase 2 — Real data  (DONE for run 1; see results/real/SUMMARY_REAL.md)
1. Run `scripts/verify_geo_accessions.py` (done once; see `results/geo_check.json`). The four datasets
   used by the prior conformal scRNA paper have only 3 to 16 GEO samples each, except one lung-cancer set
   (40, human and mouse mixed). **They cannot support donor-level calibration.**
2. Find multi-donor sources: CELLxGENE Census filtered to donors ≥ 30 per tissue; periodontitis and healthy
   gingiva scRNA (not yet identified); oral cancer / OPMD scRNA (not yet identified).
3. Record donor counts from sample metadata (not GEO sample counts) in the manifest.
4. Extract foundation-model embeddings with each model's own tooling (scGPT, Geneformer; not implemented
   here), then `scripts/prepare_inputs.py` → `donorconf calibrate/screen/report`.

Go/no-go: if no oral dataset has ≥ 10 donors with independent author labels, the oral case study is dropped
and the paper uses non-oral multi-donor atlases; the oral link then moves to future work.

## Phase 3 — Baselines to run for real  (torchCP baseline DONE 2026-09-21/22; DCATS and graph-structured method not run)
* ~~cell-level conformal from the 2025 *Bioinformatics* paper~~ DONE: its torchCP 1.0.2 calibration step
  (standard/THR, standard/APS, classwise, clustered) was run on this study's head/splits
  (`src/donorconf/baseline_torchcp.py`, `results/real/c_benchmark_blood_tcp/`). Isolates calibration only —
  the annotator's OOD detector and neural classifier are NOT included. standard-THR is numerically identical
  to this study's own cell_pooled method; classwise/clustered under-cover like this study's classwise method;
  standard-APS is the strongest published-package variant but donor_crc still wins on coverage/size trade-off.
* DCATS (R package) for composition — not run (composition likely out of scope, see below).
* The graph-structured conformal method — not run.
* Foundation models: at least two plus a PCA baseline — still open (Geneformer/scGPT not Census-hosted).

## Phase 4 — Package and release
* `pip install`-able; add a tutorial notebook on a public multi-donor atlas; scverse/AnnData integration
  (`prepare_inputs.py` already reads `.h5ad`); tag a release, deposit on Zenodo, add repository URL to
  `CITATION.cff`.
* Choose the package name after a collision search (the name `cellconformal` is taken on GitHub).

## Phase 5 — Manuscript
See `MANUSCRIPT_OUTLINE.md`. Target *Bioinformatics* or *Briefings in Bioinformatics*. Lock the protocol
(OSF) before touching real target data. A statistician should review the guarantee statements.

## Open items and known limits of this scaffold
* Foundation-model embedding extraction is **not implemented**.
* Real-data results: one exploratory run exists (`results/real/SUMMARY_REAL.md`); nothing is confirmatory.
* The donor-level shift-screen's usefulness is unproven (see PROTOCOL.md §1).
* No literature check beyond web searches has been done; PubMed/Scite connectors need authorisation.
* Modal (CPU) is used for Census retrieval; GPU would be needed only to extract embeddings the Census does not host (Geneformer, scGPT).

## Status after real-data run 1 (2026-09-21)
* Retrieved and analysed 7 real Census datasets on Modal (`pradeepaiperio`). Gingiva atlas has 34 donors
  (20 normal, 10 periodontitis, 4 gingivitis). See PROTOCOL.md §1b for outcomes against the kill criteria.
* Recommended next decision: the cell-exchangeability headline is not supported, so choose between
  (B) a less conservative guarantee under dataset shift, for example calibrating with a few LABELED target
  donors plus a donor-bootstrap confidence bound, or (C) a benchmark paper on reported-vs-realised coverage
  under dataset shift across foundation-model embeddings. Framing A alone is no longer enough.
* Still to do before any claim: run DCATS and the published cell-level conformal method as real baselines;
  add more dataset pairs (only three were run); register the protocol; a statistician review.

## Status after B and C (2026-09-21)
* B and C were both run on real Census data (Modal `pradeepaiperio` for retrieval); see PROTOCOL.md section 1c.
* **Prior art for B is substantial**: weighted conformal with few labeled target samples, and federated
  conformal risk control with risk-curve shrinkage and equal-site weighting (arXiv 2606.20115). Do not present
  B as a new method. Cite these and position B as the donor-level, scRNA-seq application with a real-data study.
* **C is the stronger paper**: a dataset-shift benchmark of reported vs realised coverage for foundation-model
  embeddings, with a label-free early-warning signal (confidence shift) and a practical remedy (B with k=3 to 6
  labeled donors). Recommended framing: C as the headline, B as the remedy, the donor-level guarantee as a footnote.
* Open before any submission: the published cell-level conformal package as a real baseline; more dataset pairs
  and tissues (only blood and three oral pairs); Geneformer and scGPT embeddings (not Census-hosted for this
  version); DCATS only if composition stays in scope (probably drop it); OSF registration; statistician review;
  a proper literature check on PubMed (connector unauthorised).

## Status after torchCP baseline + second tissue (2026-09-21/22)
* Published-package torchCP baseline run on blood (see Phase 3 above); result strengthens the paper — the
  finding is not an artefact of this study's own conformal code, since the paper's own predictor shows the same
  pattern (standard THR under-covers under shift; classwise/clustered under-cover more; APS is the best
  published-package variant but still loses to donor_crc).
* **Second tissue (lung) run and replicates C.** Six lung datasets qualified a fixed, pre-declared inclusion rule
  (`ANALYSIS_LOCK.md` Addenda C2/C3: ≥15 donors with ≥20 mapped cells each, ≥60% cells mapped, lung tissue only;
  applied to the audit BEFORE any coverage number was computed). Reported-vs-realised gap is *larger* in lung
  (71.5% of replications with |gap|>0.03, vs 47.9% in blood); the label-free confidence-shift signal replicates
  (Spearman −0.776 in lung vs −0.732 in blood); the permutation screen again fails to predict shortfall. This is
  the strongest evidence yet that C is a real, generalisable phenomenon and not a blood-specific artefact.
* Negative finding worth reporting: tcp_cluster_thr reaches donor_crc-like coverage in lung only by inflating
  mean prediction-set size ~2x — coverage number alone would be misleading without reporting set size alongside it.
* Remaining before submission: Geneformer/scGPT embeddings, OSF registration, statistician review, a
  written-strategy PubMed/Scopus search (PubMed connector needs authorising in claude.ai connector settings —
  cannot be done from this non-interactive session), Zenodo DOI, Paperpal preflight, AI-use disclosure in the
  cover letter. See `manuscript/BiB_SKELETON.md` for the auto-generated evidence pack and
  `manuscript/LITERATURE_CHECK.md` for the (web-search-only) literature check.
