# Analysis lock — real-data run 1 (self-imposed, NOT an external registration)

Written before any real-data result was computed. Status: EXPLORATORY, because PROTOCOL.md is not
externally registered and the simulation results were already seen.

## Data
Source: CELLxGENE Census 2025-11-08 (primary cells), pulled with modal/modal_census.py on Modal
workspace pradeepaiperio. Embeddings hosted by the Census: TranscriptFormer TF-Sapiens (2048-d),
TranscriptFormer TF-Exemplar-Human (2048-d), scVI (50-d, trained with batch/donor/dataset correction).
Labels: Census `cell_type` (curated ontology mapping of the data authors' labels), not model output.
Caveat stated in advance: these models were trained on Census cells that probably include these
datasets, so evaluation is not independent of pretraining.

## Experiment R1 — gingiva, donor-level splits (single dataset, 34 donors: 20 normal / 10 periodontitis / 4 gingivitis)
* Label vocabulary: cell types with >= 150 cells and present in >= 12 donors; other cells dropped and counted.
* Each replication: random donor partition into 12 reference / 12 calibration / 10 target donors.
* Head (fixed for all embeddings): standardise on reference -> PCA(64) -> logistic regression C=1.
* alpha = 0.10; methods: cell_pooled, donor_weighted, donor_crc; reported = within-calibration random half split.
* 100 replications per embedding; seeds 0..99. Metric: donor-mean coverage on target donors, set size,
  fraction of poorly covered donors (< 0.80). Label-free shift screen p-value recorded (all splits are
  exchangeable by construction, so it should flag at about the nominal 5% rate).

## Experiment R2 — blood, leave-dataset-out (three datasets, <= 80 donors each, <= 150 cells/donor)
* Label vocabulary: cell types present in all three datasets with >= 100 cells in each.
* Each replication: reference = 30 random donors of dataset A; calibration = 30 donors of dataset B;
  target = 30 donors of dataset C (all six role permutations across replications). alpha = 0.10.

## Kill criteria (from PROTOCOL.md)
H1 not supported if cell_pooled realised donor-mean coverage stays within +-0.02 of 0.90 in R1 and R2
and its reported-vs-realised gap is under 0.03. H2 screen claim only if P(fail|flagged) - P(fail|not) >= 0.15.

Results are reported whatever they are, including failures.

## Addendum A (2026-09-21) — label harmonisation for cross-dataset runs
Written after inspecting the raw Census label NAMES of the datasets (needed because a name-level
intersection left one shared type) and BEFORE any cross-dataset coverage result was computed.
* `src/donorconf/lineage.py` maps source labels to coarse types by fixed substring rules; every mapping
  and every excluded label is in `config/lineage_audit.csv`. A regression bug (substring "t cell" matching
  "mast cell") was found from that audit and fixed before running.
* R2 blood uses scheme `blood_fine` (CD4_T, CD8_T, NK, B, mono_classical, mono_nonclassical, cDC).
* R3 oral uses scheme `oral_coarse`, reference gingiva -> calibration oropharynx -> target oral SCC
  (n_ref 20, n_cal 10, n_tgt 10, min 40 cells and 4 donors per class per dataset), 60 replications.
  R3b: gingiva -> oral_normal_multisite -> oral_scc, same settings.
* The oropharynx calibration set contains both normal and tumour donors; the oral_scc target is entirely
  tumour. That disease shift is deliberate.
* Not specified in advance and therefore exploratory: which cross-dataset pairs to run beyond these.

## Addendum B (2026-09-21) — framing B (target-adaptive calibration) and framing C (benchmark)
Written after real runs R1-R4 were seen and after the simulation study of B, BEFORE any real-data result for B or
for the extended benchmark. B and C are therefore a response to R1-R4 and are exploratory.
* B methods: donor_crc (source only), donor_weighted (source), target_only_k, adaptive_crc, adaptive_plugin
  (omega grid 0..16 chosen by leave-one-target-donor-out). k in {3, 6} labeled target donors; all methods are scored
  on the remaining target donors. alpha = 0.10. Real pairs: R2 blood (permuted roles, n=30 donors per role),
  R3 (gingiva -> oropharynx -> oral SCC), R3b (gingiva -> multi-site oral -> oral SCC); 30 replications.
* C benchmark: blood datasets A-E, 24 role triples chosen by a fixed-seed draw, 2 replications each, methods
  cell_pooled, classwise (Mondrian), donor_weighted, donor_crc. Label-free covariates recorded: screen p-value,
  mean-confidence shift, calibration accuracy. Kill rule for a "label-free predictor of shortfall" claim:
  Spearman correlation between confidence shift and coverage shortfall must exceed 0.4 in absolute value.
* Not specified in advance: any further pair or method beyond the above.

## Addendum C (2026-09-21) - second tissue (lung) and published-package baseline
Written after the blood B/C results were seen and BEFORE any lung result was computed.
* Second tissue: human lung, five Census 2025-11-08 datasets chosen from the survey by donor count and a
  mostly-normal disease composition (350237e0, 6725ee8e, eb499fd8, 1b350d0a, d8da613f), tags lung_A..E, at most 150
  cells and 40 donors per dataset, same three hosted embeddings. Chosen before looking at any lung label or coverage.
* Labels: scheme `lung_coarse` in `src/donorconf/lineage.py`; the rules were written before the lung data were
  inspected and may be amended only by an addendum after the audit CSV is produced.
* Design: identical to C blood (24 fixed-seed role triples, 2 replications, 30 donors per role, alpha 0.10, kill rules
  unchanged: confidence-shift predictor needs |Spearman| >= 0.4).
* Published-package baseline: the torchCP 1.0.2 SplitPredictor (THR and APS), ClassWisePredictor and
  ClusteredPredictor used by the conformalized single-cell annotator (Bioinformatics 2025, btaf521), run on the same
  head's probabilities and the same splits. The annotator's OOD detector and neural classifier are NOT included.
  Cross-check: torchCP standard THR reproduced this package's pooled-cell sets exactly on exchangeable simulated data.
* Blood benchmark re-run with the same seeds plus the torchCP variants: `results/real/c_benchmark_blood_tcp/`.
* Not done and not claimed: Geneformer and scGPT embeddings (not hosted by the Census; would need GPU extraction).

## Addendum C2 (2026-09-21) - lung dataset inclusion, written after the lung label audit and BEFORE any lung coverage result
* The audit (`config/lineage_audit_lung.csv`) showed the first pull kept whole datasets, and lung_B and lung_D are
  multi-organ atlases (lung_B is mostly intestine and pancreas; lung_D mostly spleen and lymph nodes); lung_A is
  100% lung tissue but only 24% of its cells map to the shared vocabulary (labels such as "unknown", "stem cell").
* Inclusion rule, fixed now: a lung dataset is used only if, restricted to lung tissue, it has >= 15 donors with
  >= 20 mapped cells each and >= 60% of its cells map to the shared vocabulary. By this rule lung_A, lung_B and
  lung_D are EXCLUDED; lung_C and lung_E qualify.
* Two datasets cannot form reference/calibration/target triples, so four more large lung datasets were chosen from the
  same survey by donor count only (9f222629, f14bc322, d68a8b48, 1e6a6ef9) and pulled with tissue_general == lung
  (tags lung_F..I). The same inclusion rule applies to them. Some contain disease donors (COVID, COPD, ILD, BPD, tumour):
  disease shift between roles is part of the dataset shift being measured, as in the oral pairs.

## Addendum C3 (2026-09-22) - lung_F..I audit result, written BEFORE any lung coverage result
* `scripts/audit_lung_fghi.py` applied the Addendum C2 rule (>= 15 donors with >= 20 mapped cells each, >= 60%
  cells mapped) to lung_F..I. Result: all four qualify (frac_mapped 0.642/0.927/0.888/0.683; donors with
  >= 20 mapped cells 38/38/21/38 out of 40/40/21/40 kept). Rows appended to `config/lineage_audit_lung.csv`.
* The lung benchmark therefore uses six datasets: lung_C, lung_E, lung_F, lung_G, lung_H, lung_I. lung_A, lung_B,
  lung_D remain excluded per Addendum C2.
