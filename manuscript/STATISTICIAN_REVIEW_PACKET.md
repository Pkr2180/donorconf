# Statistician review packet

**Purpose.** A short packet for a statistician (colleague, paid reviewer, or co-author) to check the
guarantee statements and the correlational claims before submission. Not written by the reviewer — this is
the material to hand them, pulled directly from the code's own docstrings and the result files, so nothing
here is a summary that could drift from what the code actually does.

## What to ask the reviewer to check

### 1. The finite-sample guarantee claims (highest priority)
Three conformal calibrators are compared. Their own docstrings state the guarantee differences:

- `DonorLevelCRC` (`src/donorconf/conformal.py:118`): "Donor-level conformal risk control with the donor as
  exchangeable unit." — has the finite-sample guarantee, subject to α ≥ 1/(n+1) donors.
- `DonorWeightedQuantile` (`src/donorconf/conformal.py:85`): "PLUG-IN ESTIMATOR: it has NO finite-sample
  guarantee, only the asymptotic target 'donor-equal-weight mixture quantile'. Use as an efficient companion
  to DonorLevelCRC, never as a replacement where a guarantee is claimed."
- `TargetAdaptiveCRC` (`src/donorconf/adaptive.py:71-77`): "Source donors + k labeled target donors; omega
  chosen by leave-one-target-donor-out." with a `finite_sample_correction` flag: "finite_sample_correction=
  False drops the w_test term: a plug-in estimate with smaller sets and NO finite-sample statement, even in
  the limiting cases."

**Ask:** does the manuscript state, everywhere a coverage number from `donor_weighted` or the
`adaptive_plugin` variant is reported, that it is a plug-in estimate without the finite-sample guarantee? Is
the distinction between "guaranteed but conservative" (CRC) and "efficient but unguaranteed" (weighted
quantile / plug-in adaptive) stated clearly enough that a reader could not mistake one for the other?

### 2. The α ≥ 1/(n+1) floor
`README.md`: "Conservative with few donors; vacuous when α < 1/(n+1)." **Ask:** for every reported
`donor_crc` result, was the calibration-donor count checked against this floor? (Blood/lung calibration
donor counts and α = 0.10 imply a floor of 9 calibration donors; check this was respected in every run,
not just asserted.)

### 3. The correlational claims (confidence-shift signal)
Result files: `results/real/SUMMARY_FINAL.md`. Blood: confidence-shift vs. shortfall Spearman ρ = −0.732
(p < 0.001, n = 24 triples); lung: ρ = −0.776 (p < 0.001, n = 24 triples). Kill rule pre-declared in
`results/real/ANALYSIS_LOCK.md` Addendum B: |ρ| ≥ 0.4.

**Ask:**
- Is n = 24 triples (2 replications × role permutations of the same underlying datasets) an appropriate
  denominator, or does dataset reuse across triples make the effective sample size much smaller than 24 and
  the reported p-values overstated? (This is already flagged as a limitation in `PROTOCOL.md`/discussion
  notes — ask the reviewer whether the flag is sufficient or whether a cluster-robust correction is needed.)
- Is Spearman the right correlation given the shortfall metric's distribution (see `per_rep.csv` for the
  raw values), or should this be reported with a regression + prediction interval instead of a single ρ?

### 4. The permutation screen negative result
Screen p-value vs. shortfall: blood ρ = −0.297 (p = 0.159, not significant); lung ρ = 0.106 (p = 0.622, not
significant, and opposite sign to blood). **Ask:** is reporting these as a consistent "not predictive"
negative result defensible given the sign flip between tissues, or does the sign flip itself need
discussion (e.g., as evidence the screen is picking up noise rather than a real but weak signal)?

### 5. Multiple comparisons
Four label-free signals were tested against shortfall (`conf_shift`, `entropy_shift`, `screen_p`,
`head_accuracy_calibration`) in each of two tissues — eight tests total, only one pre-declared with a kill
rule (`conf_shift`). **Ask:** does `entropy_shift`'s strong correlation (ρ = 0.730 blood, 0.788 lung) need a
multiplicity caveat since it was not the pre-declared signal, even though it moves in the expected direction
and closely tracks `conf_shift`?

### 6. Published-package cross-check
`baseline_torchcp.py` reproduces the published annotator's own torchCP calibration. Cross-check: mean
absolute coverage difference between `tcp_standard_thr` and this study's `cell_pooled` is 0.0000 (blood
simulation regression test) and 0.0000–0.0013 (real data). **Ask:** is this tight enough to say "isolates
the same calibration rule," or should the residual 0.0013 in real data be explained (likely float32 ties;
worth a sentence either way)?

## What NOT to ask them to check
- Data provenance / dataset inclusion rules (Addenda A–C3 in `ANALYSIS_LOCK.md`) — that's a domain
  (single-cell biology) question, not a statistics one; a different reviewer should look at that.
- Code correctness at the implementation level — covered by the 32-test suite (`python -m pytest -q`); a
  statistician reviewing guarantee *statements* doesn't need to re-derive the code, only check the claims
  against the stated assumptions.
