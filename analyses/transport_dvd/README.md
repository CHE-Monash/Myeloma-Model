# DVd Line 2 — Calibrated Transport

Case study for the **Calibrated Transport** method: a transportability approach for *pre-market* health technology assessment (HTA) that calibrates a trial treatment effect to real-world practice using a comparator (Vd) observed in **both** the trial (CASTOR) and the registry (MRDR). Applied to second-line daratumumab + bortezomib + dexamethasone (DVd) in multiple myeloma.

This README documents the **code and outputs** for the analysis. Manuscript status and the conceptual write-up live with the paper (vault: `myeloma model/papers/2026 dvd-transport/`), not here.

## Method in brief

Conventional transportability needs real-world data on the novel treatment, which does not exist before funding. Calibrated Transport gets around this: Vd is observed in both CASTOR and the MRDR, so it acts as an efficacy–effectiveness anchor (within-treatment, across-setting) that lets the trial DVd effect be calibrated to the Australian setting and expressed as an *absolute* real-world prediction.

**Comparator: Vd, across all scenarios** — matching the comparator the PBAC actually used in its DVd decision. The decision-relevant contrast is always DVd vs Vd; the scenarios differ only in how the second-line outcomes are derived.

Empirical hook: Vd was only **2.5%** of second-line treatment in the MRDR (71 of 2,791, all registry years; the earlier 1.5% used an inconsistent denominator). This does not change the comparator (the decision was made against Vd) — it is carried as context and a limitation. Vd is also thin at second line (only **71** patients ever receive it at L2 across the whole registry; Vd is mostly a first-line regimen), so the transport anchor pools relapsed-setting Vd across lines (n = 135; see *Underlying models*).

## Scenarios

The analysis compares three ways of predicting real-world second-line outcomes. Each is a **separate dispatcher run**, selected by the `$scenario` global, and saves to its own `simulated/<scenario>/` folder.

All three scenarios estimate the same decision-relevant contrast — **DVd vs Vd** — and differ only in how the second-line outcomes are derived. Each is a separate dispatcher run.

| `$scenario` | Prediction approach | How DVd & Vd outcomes are derived |
|------------------------|------------------------|------------------------|
| `A_trial` | Traditional trial-based | CASTOR trial BCR applied directly |
| `B_transport` | **Calibrated Transport** | DVd calibrated to the Australian setting via the Vd anchor |
| `C_mrdr` | Observed (validation target) | Observed MRDR outcomes |

## Underlying models

The comparison is **DVd vs Vd** (the PBAC comparator) in all three scenarios. Vd is observed in both CASTOR and the MRDR, which is what makes it usable as the efficacy–effectiveness anchor for Calibrated Transport. Observed L2 DVd (n = 533 evaluable; 565 starts − 32 missing BCR) is the validation benchmark.

**A line-adjusted ordered logistic regression on a pooled relapsed-Vd anchor:** `ologit BCR MRDR DVd line3 line4 line5`, fit on observed data (no imputation). The anchor pools **all relapsed-setting registry Vd** (regimen code 5, second line onward; n = 135 — L2 71, L3 28, L4 14, L5 17, L6 5; MRDR = 1, DVd = 0), with the resampled CASTOR arms — CASTOR Vd (n = 231; MRDR = 0, DVd = 0, the reference cell) and CASTOR DVd (n = 238; MRDR = 0, DVd = 1). The fit yields `bL2_BCR_T = (β_MRDR, β_DVd, line3, line4, line5 | 5 cutpoints)`. **β_MRDR** is the registry-vs-trial cross-setting shift that carries the whole trial→real-world calibration; **β_DVd ≈ −1.158** is the CASTOR DVd-vs-Vd effect (the consistency check against the published trial).

**Why a pooled anchor with line indicators.** CASTOR enrolled relapsed/refractory patients with at least one prior line (median 2, range 1–9; Palumbo 2016) — it is not an L2-only population. Restricting the anchor to L2 Vd (n = 71) both under-powers it and mismatches CASTOR's line mix, which bakes a line difference into β_MRDR. Pooling all relapsed-setting Vd with a line-of-therapy indicator corrects the mismatch and roughly doubles the anchor; first-line Vd (n = 240, not CASTOR-eligible) is excluded. CASTOR's line of therapy is generated per pseudo-patient from the Palumbo 2016 prior-line shares (DVd 48.6/27.9/14.7/8.8%; Vd 45.7/30.0/13.0/11.3%), reference line L2 and L5+ collapsed into `line5`; the line effect is identified from the registry, while the CASTOR line draw balances the line mix for the comparison.

**Prediction is at L2**, where `line3 = line4 = line5 = 0`: the real-world DVd cell is (MRDR = 1, DVd = 1) and real-world Vd is (MRDR = 1, DVd = 0), so DVd − Vd in scenario B is exactly β_DVd. The override (`outcomes/sim_bcr_override.do`) builds a two-column design (MRDR, DVd) and reads the first two coefficients plus the five cutpoints, so the line coefficients are skipped at prediction.

**Patient-covariate adjustment (age, sex, ISS) was tested and rejected.** Three covariate-adjusted variants all worsened out-of-sample DVd BCR calibration, and CASTOR reports only age, sex and ISS (no ECOG, no R-ISS), so a faithful match to the base `BCR_L2` covariates is not feasible. The line-of-therapy indicator above is a correction to the *anchor's* line composition, not a patient-covariate adjustment, and is retained. The pooled Vd anchor (n = 135) is still modest, so the cross-setting test has limited power — *no evidence of* a residual effect, not *evidence of none*.

## How to run

The dispatcher `transport_dvd.do` runs **one scenario** (one or two arms) per invocation. Set the configuration globals at the top — `$scenario`, `$coeffs` (`dvd_pre` / `dvd_post`), `$int1`/`$int0`, `$line`, bootstrap flags, and the PDF report flag `$report` — then run it. The pipeline is `load_patients → mata_setup → simulation → process_data → export_results`, followed by validation and (when `$report = 1`) the PDF report. CSV export runs by default — no flag.

**The BCR-transport generator** `calibrated_transport.do` is separate from the dispatcher: it fits the calibrated-transport ordered logit (β_MRDR / β_DVd, observed `m=0` data — no `mi`) that feeds `outcomes/sim_bcr_override.do`. `do calibrated_transport.do 0` runs the deterministic point estimate **and** writes the MRDR Vd anchor baseline for Table 1 → `results/mrdr_vd_baseline.csv`; `do calibrated_transport.do 1 <min> <max>` runs the bootstrap. The Table 1 baseline tabulation runs **only on the deterministic call** (`0`), not during the bootstrap.

To produce the full comparison, run all three scenarios, then run the cross-scenario aggregation (below).

## Outputs

Output follows the model-wide convention in the root `README.md` (Result Exports for Downstream Access). Three tiers:

**Tier 1 — engine-level, per run (`core/`).** `core/export_results.do` writes the CSVs every analysis needs (per-patient summary, BCR distribution, mean cost/QALY/LY) into `simulated/<scenario>/`. It runs by default as part of the simulation pipeline (immediately after `process_data`, once per arm; skipped during bootstrap) — no flag. Not specific to this analysis; it currently emits `bcr_<stub>.csv`, `econ_<stub>.csv`, and `patients_<stub>.csv` (where `<stub> = <int>_<line>_<data>_<min_id>_<max_id>`).

**Tier 2 — analysis-level, per scenario (`analyses/transport_dvd/export_tables.do`).** *(planned)* The DVd-transport-specific tables that no other analysis shares — the treatment-mix / 2.5% table and the transport-model regression coefficients (β_MRDR = the cross-setting shift, β_DVd ≈ −1.158 = the CASTOR contrast). Called from the dispatcher at the end of a scenario run; writes into `simulated/<scenario>/`.

**Tier 3 — cross-scenario aggregation (`analyses/transport_dvd/bootstrap_summary.do`).** Written; run once after the 6×500 bootstrap exists. **All three scenarios are bootstrapped** (A = resampled CASTOR BCR, B = calibrated-transport bootstrap, C = resampled observed cohort), so every CI reflects the same sources and is comparable. It reads `simulated/<scenario>/bootstrap/{dvd,vd}_2_predicted_B<b>.dta`, computes the **MAE paired against a common C_b** each iteration (so the A-vs-B reduction CI isn't inflated by benchmark noise) with the reduction (abs + %) and 95% CI, the side-by-side DVd-vs-Vd ICER / inc-cost / inc-QALY CIs, and the per-scenario BCR distributions; writes `bcr_distributions.csv`, `mae_comparison.csv`, `icer_comparison.csv`, `results.md` (and `bootstrap_iterations.dta`) into `results/`. This tier cannot live in a single run, since each scenario is run separately.

### Read surface: `results/`

`analyses/transport_dvd/results/` is the canonical place a downstream consumer reads from. **It is git-ignored** — these are manuscript/appendix outputs, regenerated locally by running the analysis and `bootstrap_summary.do`, and are not distributed in the repo:

```         
results/
├── treatment_mix.csv       # pre-funding L2 mix; the 1.5% finding
├── bcr_distributions.csv   # predicted vs observed BCR × 6 categories, all scenarios
├── mae_comparison.csv      # MAE per scenario with bootstrap 95% CIs
├── icer_comparison.csv     # DVd vs Vd ICER, per scenario
├── mrdr_vd_baseline.csv    # MRDR Vd anchor baseline → Table 1 (from `calibrated_transport.do 0`)
└── results.md              # narrates and labels the key figures for downstream use
```

Intermediate per-scenario CSVs stay in `simulated/<scenario>/`; only the final, cross-scenario CSVs plus `results.md` belong in `results/`.

## Current structure

```         
analyses/transport_dvd/
├── README.md                # this file
│   # 1. Cohort construction
├── cohort_pool.do           # build the L2-entry pool (case-mix) once; all_l2 = windowed (0) vs all-era (1)
├── ce_cohort.do             # draw the production cohort (size N) from the pool -> canonical patient file
│   # 2. Analysis
├── calibrated_transport.do  # BCR-transport generator: `0`=deterministic (+ Table 1 baseline), `1`=bootstrap
├── transport_dvd.do         # dispatcher: runs ONE scenario (deterministic + bootstrap)
├── bootstrap_summary.do     # Tier 3 cross-scenario aggregation of the bootstrap output
│   # 3. Sample-size justification (methods; runs independently of production)
├── ce_precision.do          # per-patient sigma_pp + TSD 15 convergence figures
├── ce_sample_size.do        # parameter SD (from the PSA) + required N + appendix figure
├── coefficients/            # coefficients_dvd_pre / coefficients_dvd_post (+ bootstrap/)
├── outcomes/                # sim_bcr_override.do (+ per-scenario subfolders)
├── patients/                # pool + drawn cohort .dta
├── results/                 # cross-scenario CSVs, figures, results.md (git-ignored; local)
└── simulated/
    ├── A_trial/             # per-scenario .dta + Tier-1/2 CSVs
    ├── B_transport/
    ├── C_mrdr/
    └── report/              # per-run PDF reports
```

Planned additions: `export_tables.do`. `bootstrap_summary.do` is written; `results/` is created when it runs.

## Implementation notes

- `core/generate_report.do` already computes the BCR-by-regimen, cost, QALY and ICER quantities the CSVs need, but as `putpdf` tables and using older variable names (`cTotald`, `qTotald`) that have drifted from `process_data.do` (`cost_total_d`, `qaly_total_d`). The CSV export must read the actual `process_data` variables; ideally each summary is computed once and fed to both the PDF and the CSV so they cannot disagree.
- The dispatcher's inline ICER calculation only fires in the two-arm case. The cross-scenario ICER comparison belongs in `bootstrap_summary.do`, not the dispatcher.

## Data

- **Vd comparator / anchor:** MRDR relapsed-setting Vd, **pooled across second line onward (n = 135: L2 71, L3 28, L4 14, L5 17, L6 5)** with a line-of-therapy indicator, to match CASTOR's relapsed line mix; first-line Vd (n = 240, not CASTOR-eligible) is excluded. Plus aggregate CASTOR DVd arm characteristics.
- **Non-Vd reference:** MRDR second-line non-Vd, **same all-years window (n = 2,720)** (= 2,791 total L2 − 71 Vd). One window for both treatment groups removes the era asymmetry; a calendar-period covariate is optional.
- **DVd validation (held out):** MRDR second-line DVd, **n = 533 evaluable** (565 starts − 32 missing BCR). The DVd prediction uses no observed DVd outcomes — these are reserved for out-of-sample validation (the paper's holdout claim, replacing the earlier "pre-funding data only" framing). This raw observed distribution (CR 11.4 / VGPR 23.1 / PR 33.4 / MR 7.7 / SD 15.4 / PD 9.0 %) is the **BCR-accuracy benchmark**.

MRDR patient data access is via the MRDR Steering Committee ([mrdr.net.au](https://www.mrdr.net.au/)).

## Status

In development. Manuscript Introduction and Methods complete (Methods currently reflect the earlier bare model and need re-reconciling to the line-adjusted model — see Next steps); Results and Discussion drafted (tracked changes, 9 Jun 2026) — three tables + Figure 1, Discussion across all six subsections. Abstract and Key points remain.

**Done.** Tier-1 CSV export (`core/export_results.do`) live and wired into the dispatcher; all three scenarios (A_trial / B_transport / C_mrdr) run for the DVd and Vd arms with point-estimate CSVs in `simulated/<scenario>/`. All outcomes reported from L2 onwards.

**Method (implemented): the line-adjusted `ologit BCR MRDR DVd line3 line4 line5` on a pooled relapsed-Vd anchor.** The anchor pools all relapsed-setting registry Vd (n = 135) with a line-of-therapy indicator to match CASTOR's line mix; prediction is at L2 (line dummies = 0), so DVd − Vd = β_DVd. Patient-covariate adjustment (age/sex/ISS) was tested in three forms and worsened out-of-sample calibration, so it is not used. `calibrated_transport.do` + `outcomes/sim_bcr_override.do` implement this.

**Results (matched cohort n = 50,180; DVd vs Vd, discounted; deterministic + 500-iteration bootstrap).** Source: `results/` via `bootstrap_summary.do`. *(These figures were recorded against the earlier L2-only anchor; reconfirm them against the current line-adjusted fit before use.)*

- **BCR accuracy** (benchmark = raw observed DVd, n = 533, fixed): DVd MAE falls from **8.4 pp** (Traditional) to **5.4 pp** (Calibrated Transport); bootstrap 8.5 → 5.4 pp, reduction 36% (95% CI −0.3–6.2 pp). The reduction is **not significant on the symmetric CI**, but Calibrated Transport is more accurate in **96.0% of bootstrap replicates** (≥VGPR closer in 96.0%) — report this exceedance probability rather than significance.
- **Economics:** ICER A **$456k** → B **$691k** → C **$2.64M** (bootstrap $459k / $711k / $2.40M; B 95% CI $487k–1.52M). ΔQALY A 0.285, B 0.171, C 0.039 (C crosses zero → C's ICER CI is uninterpretable; use its inc-cost/inc-QALY CIs). Calibrated Transport moves the ICER off the over-optimistic trial value toward observed reality — the policy-relevant result, independent of MAE significance.
- All three arms bootstrapped for comparability (so the reduction CI is conservative). Vd cost cap (Next steps #2) still to apply before final ICERs.

## Next steps (pick up here)

1.  **Re-fit the risk equations using all 71 L2 Vd and all 2,720 non-Vd** (same all-years window for both; drop the pre-funding restriction; keep the DVd holdout). This stabilises the Vd category in `bL2_BCR` and clears the sparse-Vd bootstrap failures (which were in the unused TXR_L2 regression anyway). Then re-run all three scenarios + bootstrap.
2.  **Finish the Vd cost in `core/process_data.do`.** Vd is now costed (`cVd = 724`, regimen code 5) — done. Remaining: per the MSAG 2022 guideline (Table 8), CASTOR Vd is **8 × 21-day cycles** (bortezomib D1,4,8,11), so switch the Vd line from `/28` (uncapped) to a 21-day basis with an 8-cycle cap, mirroring VCd: `cVd * min(8, TXD_L2 * 30.4375 / 21)` — confirming `cVd` is a per-21-day-cycle figure. Also tidy the `cOther` comment (no longer includes Vd).
3.  **Bootstrap + Tier 3.** Run the 6×500 bootstrap (all three scenarios, both arms), then run `bootstrap_summary.do` (written) → `results/`. All three scenarios bootstrapped; MAE paired against a common C_b; ICER CIs side-by-side and comparable.
4.  **Tier 2 (`export_tables.do`).** Treatment-mix / 2.5% table and the transport-model coefficients (β_MRDR = cross-setting shift; β_DVd ≈ −1.158 = CASTOR contrast).
5.  **Manuscript framing follow-through** (vault `todo.md`): **reconcile Methods/Assumptions/Limitations to the line-adjusted transport model** (pooled relapsed-Vd anchor + line indicator) and add the covariate-test robustness result; comparator is Vd across all scenarios; pooled-SoC removed entirely; "pre-funding only" replaced by the DVd-holdout framing; relapsed Vd pooled across lines (n = 135).
6.  Draft Results from `results/`, then Discussion → Abstract → Conclusion.
