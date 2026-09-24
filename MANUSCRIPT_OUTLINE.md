# Manuscript outline (skeleton — no results are claimed)

Working title: *Donor-Level Conformal Calibration for Single-Cell Foundation Models Under Dataset Shift*
Alternates: *Cells Are Not Exchangeable: Donor-Aware Uncertainty for Transferred Single-Cell Annotation*.
Rule from PROTOCOL.md: the final title must match what the real-data results support.

Target: *Bioinformatics* (original paper + software) or *Briefings in Bioinformatics*.
Every number in the final manuscript must be generated from result files, never typed.

## Abstract (structure)
Motivation (donors, not cells, are the independent unit) · Results (filled from real-data outputs only) ·
Availability (package, benchmark, DOI) · Contact.

## 1 Introduction
* Foundation-model annotations are increasingly transferred across datasets and used for downstream
  statistics; uncertainty is usually reported per cell.
* Prior conformal work for scRNA-seq calibrates cell by cell within a dataset; its own limitations
  section names disease and batch shift. State the gap precisely and cite it fairly.

## 2 Methods
1. Setting and notation; donor as exchangeable unit.
2. Cell-pooled split conformal (baseline). Donor-weighted quantile. Donor-level conformal risk control
   with the finite-sample statement and its 1/(n+1) requirement.
3. Label-free donor-level shift screen and its stated limits.
4. Donor-level PPI++ for composition; baselines including DCATS.
5. Simulator (for validity checks only).

## 3 Results (each panel fed by a script)
| Figure | Content | Source |
|---|---|---|
| 1 | Schema: donors, cells, the three calibrations | diagram |
| 2 | Reported (random split) vs realised (new-donor) coverage, by donor heterogeneity and size–difficulty correlation | `run_h1_coverage.py` then real data |
| 3 | Coverage vs set size trade-off across donor counts (valid-but-large vs efficient) | H1 |
| 4 | Shift screen: power by shift type; failure-prediction value | H2 |
| 5 | Composition-interval coverage and width vs donors and labeled donors | H3 |
| 6 | Real-data leave-dataset-out coverage across foundation models and PCA | real data |
| Supp | Failure cases; minimum-donor table; all scenarios | all |

## 4 Discussion
Validity versus efficiency; when few donors make guarantees expensive; what a pass on the screen does not
mean; exchangeability assumptions; label-space limits; why composition inference needs some gold labels.

## 5 Availability
Package, tagged release with DOI, benchmark splits, one-command rerun, exact accessions and hashes.

## Reporting checklist before submission
Datasets verified and donor counts recorded · protocol locked before target data · DCATS run · two
foundation models plus baseline · failure cases shown · statistician review · literature check on PubMed.
