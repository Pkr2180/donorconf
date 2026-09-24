# AI-use disclosure — draft text

*Briefings in Bioinformatics* requires AI-use disclosure in both the cover letter and the manuscript
(Methods or Acknowledgements). Below are drafts for both; edit the specifics (tool name/version, what it did
and did not do) to match what actually happened by the time you submit, and to your own voice.

## For the cover letter (paste and edit)
> During preparation of this manuscript, the corresponding author used an AI coding assistant (Claude,
> Anthropic) to help build and run the data-retrieval and analysis pipeline (Modal-based CELLxGENE Census
> retrieval, conformal-prediction implementations, and result-summarisation scripts), and to prepare a
> structured evidence pack (tables and figures) from the resulting output files. The assistant did not draft
> the manuscript text; all narrative sections were written by the author(s), who take full responsibility for
> the accuracy of the reported results and interpretation. All code is available at [repository URL] and all
> numerical results in the manuscript are generated directly from the archived result files.

## For Methods or Acknowledgements (paste and edit)
> **AI-use statement.** Data retrieval, the conformal-prediction pipeline, and the generation of summary
> tables and figures from result files were implemented and executed with the assistance of an AI coding
> assistant (Claude, Anthropic, [model + approximate dates]). Manuscript text was written by the author(s).
> Code and a full record of the analysis pipeline, including a pre-analysis lock document with dated
> addenda, are available at [repository URL / DOI].

## Things to double-check before submitting
- [ ] Confirm the exact tool/model name and rough date range you want disclosed (the pipeline work spanned
      several sessions in September 2026).
- [ ] Confirm the repository URL / Zenodo DOI is live before the disclosure references it.
- [ ] Check whether the journal wants this only in the cover letter, only in Methods, or both — the
      guidelines page linked earlier said "cover letter and Methods or Acknowledgements," so both are safest.
- [ ] If a statistician reviewed the guarantee statements, disclose that separately (see
      `STATISTICIAN_REVIEW_PACKET.md`) — it is a human review, not an AI-use item, but reviewers will ask
      about both.
