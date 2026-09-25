# Zenodo DOI — prep and steps

## Status (2026-09-25): repo re-synced with the 5-embedding results — Zenodo connection is what's left

Re-synced 2026-09-25: the release repo was two rounds of results behind (Geneformer/scGPT embeddings for
all 11 datasets, then the C benchmark re-run across all 5 embeddings). Copied the current working copy
over again, dropped ephemeral debug logs (raw stdout from the Modal GPU runs — not scientific artifacts)
via an expanded `.gitignore`, verified 32/32 tests pass independently in the release repo, committed and
pushed (`bf63d56`). `CITATION.cff` now also has an author email. Everything below still applies — only the
GitHub OAuth / release-tagging steps remain, and those need your login, not mine.

## Status (2026-09-24): repo extracted and pushed — Zenodo connection is what's left

The blocker found on 2026-09-22 (`donorconf/` lived inside one large personal monorepo rooted at the
whole `D:\` drive, no GitHub remote) is resolved:

- Extracted `donorconf/` (minus the ~900 MB of raw `.npz` Census pulls, which are regenerable via
  `modal/modal_census.py` and not suited to git) into a fresh local repo at
  `C:\Users\PRADEEP KUMAR\donorconf-release`.
- Verified it is fully self-contained: `python -m pytest -q` passes 32/32 there, independent of the
  original monorepo.
- Created **https://github.com/Pkr2180/donorconf** (currently **private**, per your choice) and pushed both
  commits (initial release + the `CITATION.cff` URL fix).
- `CITATION.cff`'s `repository-code` field is filled in, in both the release repo and the original
  `d:\pending work\bioinformatics - top q1\donorconf\CITATION.cff`.

## What's still yours to do

1. **Make the repo public when you're ready.** Zenodo only archives and mints DOIs from *public* GitHub
   repositories — https://github.com/Pkr2180/donorconf/settings → Danger Zone → Change visibility. There's
   no rush; do this whenever you're comfortable with the code being publicly visible (e.g. once you're
   closer to submission, or whenever you'd like the DOI to exist).
2. **Connect Zenodo to GitHub**: zenodo.org → Settings → GitHub → toggle the `donorconf` repository on.
   Needs your Zenodo account, logged in via GitHub OAuth — this step is yours, not mine.
3. **Tag a release on GitHub** (e.g. `v0.1.0`, matching `CITATION.cff`'s `version: 0.0.1` — bump that field
   to match whatever tag you use, e.g. `0.1.0`). Zenodo automatically archives the release and mints a DOI
   within a few minutes of the tag being pushed.
4. **Send me the DOI.** I'll add the Zenodo badge/DOI to `README.md`, `CITATION.cff`, and the manuscript's
   Data and code availability section (all sourced from files already in the repo, no new numbers invented).

### Suggested v0.1.0 release notes (draft, unchanged from 2026-09-22)
> Pre-release research code for donor-level conformal calibration of single-cell foundation-model
> annotations under dataset shift. Includes: donor-level conformal risk control and donor-weighted quantile
> calibrators; a target-adaptive weighted-CRC remedy using a few labeled target donors; a benchmark of
> reported-vs-realised coverage across blood (5 datasets) and lung (6 datasets) tissue-shift scenarios on
> CELLxGENE Census 2025-11-08 data across five embeddings (three Census-hosted: TF-Sapiens, TF-Exemplar,
> scVI; two foundation-model: Geneformer-V2-104M, scGPT-human); a
> published-package (torchCP 1.0.2) calibration baseline; 32 passing tests. See `results/real/SUMMARY_FINAL.md`
> and `manuscript/BiB_SKELETON.md` for the full evidence pack. Exploratory research code — not yet
> peer-reviewed; see `results/real/ANALYSIS_LOCK.md` for what was pre-declared versus exploratory.

## Notes on what's in the release repo vs. the working copy
- The release repo (`C:\Users\PRADEEP KUMAR\donorconf-release`, pushed to GitHub) is a **snapshot**, not a
  live mirror of `d:\pending work\bioinformatics - top q1\donorconf\`. Any further changes here (new
  results, updated manuscript drafts) need to be copied over and re-committed before the next release tag —
  ask me to do this before tagging, so the tagged release matches what the manuscript actually cites.
- `results/real/data/*.npz` (the raw Census pulls) were deliberately excluded from git (a `.gitignore` rule
  plus a `results/real/data/README.md` explaining why and how to regenerate them). If you want the raw data
  archived too, Zenodo can host it as an additional upload attached to the same DOI record, separate from
  the GitHub-archived code — that's a manual step on zenodo.org after the DOI exists, not something GitHub
  release tagging does automatically.
- A `LICENSE` file (MIT, matching `CITATION.cff`) was already present in the working copy and carried over
  into the release repo — confirmed, not just assumed.
