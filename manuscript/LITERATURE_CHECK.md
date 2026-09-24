# Literature check

**Method (updated 2026-09-22).** The 2026-09-21 pass was web search only. On 2026-09-22 the user authorised the
PubMed, bioRxiv, Scite, and Consensus connectors and this pass re-ran with them directly: PubMed `search_articles`
(controlled-vocabulary translation, shown per query below), Scite `search_literature` (220M+ papers, full-text/
citation-graph aware), Consensus `search` (title/abstract across Semantic Scholar/PubMed/Scopus/arXiv), and
bioRxiv `get_preprint`/`search_preprints` (date/category only -- **no bioRxiv keyword search exists**, so bioRxiv
preprints were reached only when a DOI was already known from Scite/Consensus). This is closer to a written-strategy
search than the 2026-09-21 pass but is still not a full systematic review (no deduplication protocol, no second
reviewer, no PRISMA-style log). Queries actually run are listed below so they can be reproduced or extended.

## PubMed queries run (2026-09-22)
1. `conformal prediction single-cell RNA sequencing cell type annotation coverage` -> 1 hit: PMID 40973204
   (López-De-Castro et al., already known -- see Read #1 below).
2. `donor-level conformal prediction OR "conformal risk control" single cell` -> 0 hits.
3. `single-cell foundation model dataset shift benchmark cross-dataset annotation reliability` -> 0 hits.
4. `conformal prediction cell type annotation exchangeability patient donor` -> 0 hits.

**Consequence:** PubMed indexes no article combining donor-/patient-level exchangeability with conformal cell-type
annotation. The one relevant PubMed hit is the paper already treated as prior art for "coverage fails under shift."

## Read in full (via Scite metadata + PubMed record; abstracts and citation context checked)
1. **López-De-Castro M, García-Galindo A, González-Gomariz J, Armañanzas R. Conformal inference for reliable single
   cell RNA-seq annotation.** *Bioinformatics* 41(10):btaf521 (2025). PMID 40973204, PMC12506889, DOI
   [10.1093/bioinformatics/btaf521](https://doi.org/10.1093/bioinformatics/btaf521). Gold OA, CC-BY, no retraction/
   correction/erratum on record (checked via Scite `editorialNotices`).
   - Cell-level conformal annotation with an OOD detector; standard, classwise and clustered taxonomies via torchCP.
     Calibration set = 40% of the training samples (cells), not donors.
   - **Already reports** that exchangeability is violated under disease shift (lung tumour query) and coverage is
     closer to nominal for a healthy query; uses leave-one-subject-out on pancreas. So "coverage fails under dataset
     shift in scRNA-seq" is NOT new. What this study adds: a quantified reported-vs-realised gap across many
     dataset pairs and foundation-model embeddings (blood + lung, 11 datasets total), a label-free early-warning
     signal, and a head-to-head against this paper's own torchCP calibration variants in both tissues.
   - Code: https://github.com/digital-medicine-research-group-UNAV/conformalized_single_cell_annotator (not pip
     installable; not licensed in the README). Its conformal step is torchCP 1.0.2 -- matches this study's baseline.
2. **Wang B, Qiao X. Conformal Prediction Under Generalized Covariate Shift with Posterior Drift.**
   arXiv:[2502.17744](https://arxiv.org/abs/2502.17744) (2025). Confirmed via Scite still arXiv-only (no journal
   version found), OA green, CC-BY. A weighted conformal classifier combining abundant source data with scarce
   labeled target samples, with a coverage guarantee. **Prior art for framing B**; B is an application to
   donor-level scRNA-seq with a leave-one-target-donor-out choice of the weight, not a new method.

## New prior art found this pass (2026-09-22) -- important for framing
3. **Giannopoulos M, Shen Y, Zavlanos MM. Conformal Risk Minimization for Semi-Supervised Domain Adaptation via
   Optimal Transport.** arXiv:[2608.23153](https://arxiv.org/abs/2608.23153) (2026, preprint only). Found via
   Consensus and confirmed via Scite (OA green, CC-BY). Integrates the conformal-risk-control objective directly
   into semi-supervised domain-adaptation *training* (not post-hoc calibration) using optimal-transport pseudo-
   labels, specifically for the "limited-labeled-target-data regime" -- evaluated on Office-Home and a skin-lesion
   (HAM10000 -> ISIC/PH2/Derm7pt) transfer setting, 1-10 labeled target samples per class.
   - **This is closer prior art for B than previously found.** It targets exactly the "few labeled target examples
     + distribution shift + conformal guarantee" problem B addresses, though it acts at training time (retraining
     the classifier) rather than at calibration time (reweighting an already-trained model, as B does), and is
     imaging/domain-adaptation, not single-cell. Framing B's contribution must now be stated relative to *both*
     Wang & Qiao (2025, post-hoc weighted CP under covariate shift) and Giannopoulos et al. (2026, train-time CRM
     for SSDA) -- B is a post-hoc, donor-level instantiation of an idea with prior art on both the "weighted
     calibration" side and the "few-labeled-target-domain" side.
4. **Shahid NF. When Average Calibration Fails: Site-Conditional Federated Conformal Risk Control.**
   arXiv:[2606.20115](https://arxiv.org/abs/2606.20115) (2026, preprint only). Confirmed via Scite (title corrected
   from the 2026-09-21 pass, which had only a paraphrase -- "federated CRC risk-curve shrinkage" -- not the real
   title). OA green. **Read in full 2026-09-24** (Scite `read_fulltext` returned no indexed body text for this
   arXiv record, so the full text was read via the arXiv abstract/PDF page directly).
   - **Setting:** 20 hospitals, federated brain-tumour MRI segmentation (FeTS-2022, 1,251 subjects), predicting
     sets that control false-negative rate. Genuinely federated: "no patient-level images, masks, or per-volume
     scores leave any site" -- only summary risk-curve statistics are shared between sites.
   - **Their core empirical finding is structurally the closest match to this study's C headline finding found in
     the whole search:** "naive pooled CRC protects the average hospital but violates coverage at 40% of
     individual institutions," with the worst site 7.8 percentage points above the target false-negative rate.
     This is the same shape of result as this study's reported-vs-realised gap (47.9% of blood replications,
     71.5% of lung replications with |gap| > 0.03) -- pooled/marginal calibration looking fine on average while
     failing badly for a specific exchangeable unit (their hospital; this study's donor).
   - **Their remedy ("risk-curve shrinkage") is methodologically closer to B than first appreciated:** each site
     transmits its empirical risk curve plus a hyperparameter n0 that interpolates between site-specific local
     calibration and sample-size-weighted pooled calibration, with n0 chosen by **leave-one-site-out** sensitivity
     analysis (n0=19 in their experiments). B chooses its interpolation weight omega between source-only and
     target-adaptive weighting by **leave-one-target-donor-out** -- the same "leave-one-unit-out selection of an
     interpolation weight between local and pooled calibration" idea, applied to a different unit (site vs.
     donor) and a different task (segmentation FNR control vs. cell-type set coverage). State this parallel
     explicitly in Methods/Discussion rather than only listing Shahid as a citation.
   - **Key differences to state in the manuscript:** (i) genuinely federated/privacy-constrained (data never
     pooled) vs. this study's centrally-pooled Census data (no federation constraint, full access to all cells);
     (ii) site = institution (large n, heterogeneity from scanners/protocols) vs. donor = individual (much smaller
     per-unit n, heterogeneity from disease/biological state); (iii) their guarantee is finite-sample and
     per-site by construction (the shrinkage estimator is designed for it); this study's C is descriptive/
     benchmarking (measuring how badly a *marginal* guarantee is violated per donor), while B is the
     guarantee-bearing remedy, evaluated only empirically here, not proved; (iv) domain: medical imaging
     segmentation vs. single-cell transcriptomic classification.
   - **Consequence for novelty claim:** "average calibration can fail conditional on a natural exchangeable
     subunit of the data" is now confirmed NOT a novel observation in the conformal-prediction literature at
     large -- Shahid (2026) shows it for hospitals in federated segmentation. This study's contribution must be
     stated precisely as: the first (as far as this search found) quantified benchmark of this failure mode for
     single-cell foundation-model embeddings at donor granularity, across two tissues and against a published
     package's own calibration variants, plus a label-free early-warning signal -- not as "discovering" the
     general phenomenon.
5. **Mehrtens H, Bucher T, Brinker TJ. Pitfalls of Conformal Predictions for Medical Image Classification.**
   arXiv:[2506.18162](https://arxiv.org/abs/2506.18162) (2025, preprint only). Confirmed via Scite, OA green, CC-BY.
   General cautionary piece on conformal prediction failure modes in a clinical imaging setting; relevant as a
   "conformal predictions can mislead in practice" citation for the Introduction/Discussion, not a direct
   competitor.
6. **Bhattacharyya A, Barber RF. Group-weighted conformal prediction.** *Electronic Journal of Statistics* (2024).
   Found via Consensus. Handles covariate shift driven by a *known, finite number of groups* (e.g. stratified
   sampling) -- conceptually close to treating donors as groups, but assumes the group-shift structure is known in
   advance and does not address measuring a reported-vs-realised gap empirically. Worth one related-work sentence;
   not a direct competitor to C or B.

## Broader covariate-shift/conformal-risk-control theory confirmed present (not direct competitors, standard cites)
Confirmed via Consensus/Scite as existing, correctly attributed literature that the Methods/Introduction should
cite where relevant: Tibshirani, Barber, Candès, Ramdas (NeurIPS 2019, weighted conformal under covariate shift);
Barber, Candès, Ramdas, Tibshirani (*Annals of Statistics* 2022, "Conformal prediction beyond exchangeability");
Angelopoulos et al. (Conformal Risk Control, 2022/ICLR 2024); Gibbs & Candès (conditional guarantees, *JRSS-B*
2023); Cauchois, Gupta, Duchi (robust validation under f-divergence shift, *JASA* 2020); Yang, Kuchibhotla, Tchetgen
Tchetgen (doubly robust calibration under covariate shift, *JRSS-B* 2022); Xu et al. (Wasserstein-regularized CP
under general distribution shift, arXiv 2025). None of these apply to single-cell data or benchmark a
reported-vs-realised gap; they are general conformal-under-shift theory that motivates this study's setting.

## Single-cell foundation-model benchmark landscape (adjacent, not conformal)
Several 2025-2026 preprints/papers benchmark scGPT/Geneformer/UCE/etc. cross-dataset generalisation and zero-shot
performance (Kedzierska et al., *Genome Biology* 2025, "Zero-shot evaluation reveals limitations of single-cell
foundation models"; VCBench, bioRxiv 2026; "Harmonised benchmarking of foundation models for single-cell and
spatial transcriptomics," 2026; Gaballa et al., bioRxiv 2026, "Benchmarking single-cell foundation models in a
zero-shot setting"). None of their abstracts mention conformal prediction or a coverage guarantee -- they benchmark
point-prediction accuracy/embedding quality, not uncertainty calibration. Cite 1-2 of these in the Introduction to
motivate "foundation-model embeddings are known to generalise unevenly across datasets," which is the premise this
study's coverage benchmark depends on.

## Correction to the 2026-09-21 pass
- **CellBench-LS (previously cited as bioRxiv 2026.04.01.714123)**: this DOI does **not resolve** via the bioRxiv
  connector on either the biorxiv or medrxiv server (`get_preprint` returns "No preprint found"). The 2026-09-21
  web-search pass had this from a search snippet, not a verified page. **Do not cite this DOI in the manuscript
  until it is independently re-found and verified** (it may be mistyped, unpublished, or from a different
  server/aggregator). Treat as unverified, not as confirmed related work.

## Consequence for the manuscript
- Claim: a quantified benchmark of coverage under dataset shift for conformal annotation on foundation-model
  embeddings (blood + lung, 11 datasets), with a label-free early warning. Do not claim "first to show coverage
  fails" -- Shahid (2026) shows the same failure shape for federated hospital sites, so state novelty precisely as
  the first quantified benchmark of this failure mode *for single-cell foundation-model embeddings at donor
  granularity*, not as discovering the general phenomenon.
- Claim about B: an empirical evaluation of an existing idea, now with **three** strands of prior art to cite --
  Wang & Qiao (2025, post-hoc weighted CP), Giannopoulos et al. (2026, train-time CRM for few-labeled-target SSDA),
  and Shahid (2026, LOSO-selected shrinkage weight for federated per-site CRC -- methodologically parallel to B's
  LOTO-selected omega) -- plus Tibshirani et al. (2019) for the foundational weighted-conformal-under-covariate-
  shift result.
- **Done (2026-09-24):** Shahid (arXiv:2606.20115) read in full via WebFetch (Scite had no indexed full text for
  this record); the explicit donor-vs-site comparison paragraph is written above and in `build_bib_pack.py`'s
  Introduction guidance.
- Still to do before submission: verify every reference against the publisher/arXiv page one more time at
  submission time (three of the newly found sources are 2026 preprints and may be updated or published between now
  and submission); Scopus was not searched (no Scopus connector available) -- note this limitation in the
  manuscript's search-strategy statement.
