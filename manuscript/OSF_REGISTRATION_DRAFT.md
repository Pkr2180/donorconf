# OSF pre-registration draft

**Status of this draft.** This is a *post-hoc* write-up of a self-imposed analysis lock
(`results/real/ANALYSIS_LOCK.md`) that was time-stamped only in this repository's git history, not on an
external registry. Registering it on OSF now will NOT make the blood/lung/B/C analyses reported in the
manuscript prospectively pre-registered — those results were already seen before this text is submitted.
Say so explicitly in the registration itself (OSF has a field for "as-planned vs. as-run" and for
retrospective registrations) and in the manuscript's Methods. What OSF registration still buys you here:
(a) a fixed, timestamped, third-party record of the analysis rules going forward — for a companion or
revision analysis, or a Registered Report second study — and (b) a place to point reviewers who ask about
pre-registration. Decide before submitting whether to register this as a *retrospective* protocol (labelled
as such) or reserve OSF for a genuinely prospective follow-up.

## How to submit (you do this; ~10 minutes)
1. Create an OSF account at osf.io (or log in) with the same identity as your ORCID (0000-0002-6653-4123).
2. New project → "Register" → choose a template (a generic "OSF Preregistration" template is fine for a
   retrospective protocol; do not use "Prospective" wording).
3. Paste the sections below into the matching fields.
4. Set visibility as you prefer (private with embargo, or public immediately) and submit.
5. Send me the resulting OSF DOI/URL and I will add it to the manuscript's Methods and Data/code availability
   section.

---

## Title
Donor-level conformal calibration for single-cell foundation-model annotation under dataset shift: a
benchmark of reported vs. realised coverage (retrospective protocol record)

## Research questions / estimands
- **E1 (donor coverage).** For a new donor drawn from the target population, the expected fraction of its
  cells whose true label lies in the conformal prediction set, averaged with equal weight per donor.
  Nominal target 1 − α (α = 0.10).
- **E2 (reported-vs-realised gap).** The difference between coverage reported by a within-dataset random
  cell-level split and coverage realised on held-out donors/datasets.
- **E3 (label-free early warning).** Whether a confidence-shift statistic, computed without target labels,
  predicts the coverage shortfall.

## Hypotheses and kill criteria
Reproduced verbatim from `PROTOCOL.md` §2 and `results/real/ANALYSIS_LOCK.md` — do not paraphrase; copy those
sections directly into the OSF form so the record matches the repository exactly. Include H1 through H5
(H5 = the lung second-tissue replication, added 2026-09-21/22).

## Data
CELLxGENE Census 2025-11-08, pulled via `modal/modal_census.py` (Modal workspace `pradeepaiperio`). Datasets,
inclusion rules, and exclusions are listed in `results/real/ANALYSIS_LOCK.md` Addenda A–C3 and
`config/lineage_audit.csv` / `config/lineage_audit_lung.csv`. State plainly in the OSF record that dataset
inclusion/exclusion for lung (Addenda C2/C3) was decided from an audit of label coverage, not from any
coverage result, and that this order is verifiable from file timestamps / git commit order.

## Analysis plan
Split conformal (LAC/THR and APS scores); cell-pooled, class-conditional (Mondrian), donor-weighted quantile,
and donor-level conformal risk control (CRC) baselines; target-adaptive weighted CRC (weight chosen by
leave-one-target-donor-out) as the remedy (framing B); published-package torchCP 1.0.2 baselines
(`src/donorconf/baseline_torchcp.py`) isolating the calibration step of López-De-Castro et al. (2025,
*Bioinformatics*, btaf521). Head: standardise → PCA(64) → logistic regression, fixed across embeddings.
α = 0.10 throughout.

## Deviations already known (disclose these up front)
- The protocol was fixed in stages via dated addenda as new tissues/baselines were added, not in one
  prospective document — this is why the registration is retrospective, not prospective.
- Simulation results (Phase 0) were seen before any real-data hypothesis was finalised.
- Composition inference (donor-level PPI, H3) and DCATS were part of the original plan but were not carried
  through to the final manuscript scope (state whichever is true by the time you submit).
