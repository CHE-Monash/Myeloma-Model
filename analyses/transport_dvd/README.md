# DVd Line 2: Calibrated Transport

Case study for the **Calibrated Transport** method: a transportability approach for *pre-market* health technology assessment (HTA) that calibrates a trial treatment effect to real-world practice using a comparator (Vd) observed in **both** the trial (CASTOR) and the registry (MRDR). Applied to second-line daratumumab + bortezomib + dexamethasone (DVd) in multiple myeloma.

This README documents the **code and outputs**. The method, its assumptions and the results are reported in the associated manuscript.

> **Status: manuscript submitted (28 October 2026), under review.** The analysis is frozen at the submitted result set (run 10 August 2026 on the v3.0 engine). Code and fitted coefficients are in the repository; the simulation cohorts, simulated output, `results/` and the sample-size scripts are git-ignored and regenerate locally.

## Method in brief

Covariate-based transportability methods typically adjust trial results for measured case mix; comparator-based calibration instead combines a trial treatment contrast with outcomes under the comparator in the target setting, under explicit transportability and representativeness assumptions. Calibrated Transport applies this: Vd is observed in both CASTOR and the MRDR, so the registry Vd patients form a **calibration cohort** (the same regimen in two settings) that lets the trial DVd effect be calibrated to the Australian setting and expressed as an *absolute* real-world prediction.

**Comparator: Vd, across all scenarios**, matching the comparator the PBAC used in its DVd decision. The decision-relevant contrast is always DVd vs Vd; the scenarios differ only in how the second-line outcomes are derived.

Vd was only **2.5%** of second-line treatment in the MRDR (71 of 2,791, all registry years). This does not change the comparator; it is carried as context and a limitation. Because Vd is thin at second line, the calibration cohort pools relapsed-setting Vd across lines (135 episodes; see *Underlying models*).

## Scenarios

Three ways of predicting real-world second-line outcomes, all estimating the same **DVd vs Vd** contrast and differing only in how the outcomes are derived. Each is a separate dispatcher run, selected by `$scenario`, saved to its own `simulated/<scenario>/` folder.

| `$scenario` | Prediction approach | How DVd and Vd outcomes are derived |
|----|----|----|
| `A_trial` | Trial-based | CASTOR trial BCR applied directly to both arms |
| `B_transport` | **Calibrated Transport** | DVd calibrated to the Australian setting via the Vd calibration cohort; Vd from the registry equation |
| `C_mrdr` | Registry-based (economic reference) | DVd and Vd arms from the core registry second-line response equation (patient-specific, MI-fitted on all second-line records with a regimen term). Not the raw observed DVd marginal, which is the accuracy benchmark only |

## Underlying models

The comparison is **DVd vs Vd** (the PBAC comparator) in all three scenarios. The recorded second-line DVd response distribution (n = 399 of 565 with a recorded response) is the accuracy benchmark (see Data).

**A line-adjusted ordered logistic regression:** `ologit BCR MRDR DVd line3 line4 line5`, fitted on observed data (no imputation). The calibration cohort is **all relapsed-setting registry Vd** (regimen code 5, second line onward; 135 episodes, of which n = 93 with a recorded response enter the complete-case fit; MRDR = 1, DVd = 0), stacked with the resampled CASTOR arms: CASTOR Vd (n = 231; MRDR = 0, DVd = 0, the reference cell) and CASTOR DVd (n = 238; MRDR = 0, DVd = 1). The fit yields `bL2_BCR_T = (β_MRDR, β_DVd, line3, line4, line5 | 5 cutpoints)`. **β_MRDR = 0.665** is the **calibration parameter**, the registry-versus-trial difference that carries the whole efficacy–effectiveness adjustment; **β_DVd = −1.132, conditional on the line indicators**, is consistent with the CASTOR DVd-versus-Vd contrast (the consistency check against the published trial).

**Why a pooled cohort with line indicators.** CASTOR enrolled relapsed/refractory patients with at least one prior line (median 2, range 1–9; Palumbo 2016); it is not a second-line-only population. Restricting the calibration cohort to second-line Vd (n = 71) both under-powers it and mismatches CASTOR's line mix, which bakes a line difference into β_MRDR. Pooling all relapsed-setting Vd with line-of-therapy indicators corrects the mismatch and roughly doubles the cohort; first-line Vd (not CASTOR-eligible) is excluded. CASTOR's line of therapy is generated per pseudo-patient from the Palumbo 2016 prior-line shares (DVd 48.6/27.9/14.7/8.8%; Vd 45.7/30.0/13.0/11.3%), with reference line L2 and L5+ collapsed into `line5`; the line effect is identified from the registry, while the CASTOR line draw balances the line mix for the comparison.

**Prediction is at L2**, where `line3 = line4 = line5 = 0`: the real-world DVd cell is (MRDR = 1, DVd = 1). The override is dispatched for the DVd arm only; scenario B's Vd arm is the core registry second-line equation (the same Vd arm as C), so DVd − Vd in scenario B is **not** β_DVd. Scenario A overrides both arms with the CASTOR distributions, assigned by prognostic rank; B's DVd arm uses a common predicted distribution with stochastic CRN draws. The override (`outcomes/sim_bcr_override.do`) builds a two-column design (MRDR, DVd) and reads the first two coefficients plus the five cutpoints, so the line coefficients are skipped at prediction.

**Patient-covariate adjustment (age, sex, ISS) was tested and rejected.** Three covariate-adjusted variants all worsened out-of-sample DVd BCR calibration, and CASTOR reports only age, sex and ISS (no ECOG, no R-ISS), so a faithful match to the base `BCR_L2` covariates is not feasible. The line-of-therapy indicators are a correction to the calibration cohort's line composition, not a patient-covariate adjustment, and are retained. The calibration cohort is still modest, so the cross-setting test has limited power: *no evidence of* a residual effect, not *evidence of none*.

**Proportional odds is pre-specified.** `outcomes/calibrated_transport_prop.do` holds the partial proportional-odds diagnostic, kept as a record; it is not part of the operative pipeline.

## Data

- **Calibration cohort:** MRDR relapsed-setting Vd, pooled across second line onward (135 episodes, 93 with a recorded response) with line-of-therapy indicators, to match CASTOR's relapsed line mix; first-line Vd is excluded. Plus the published aggregate CASTOR response distributions for both arms.
- **DVd validation (held out):** recorded second-line DVd benchmark, **n = 399 of 565** with a recorded response (CR 8.77%, VGPR 22.81%, PR 34.84%, MR 9.27%, SD 13.53%, PD 10.78%). The DVd prediction uses no observed DVd outcomes; these are reserved for out-of-sample validation and are the **BCR-accuracy benchmark**, held fixed.

MRDR patient data is restricted; access is via the MRDR Steering Committee ([mrdr.net.au](https://www.mrdr.net.au/)). Its location is set per machine in `config.do` (`$data_dir`) and never committed. Nothing patient-level is in this folder: the repository carries code, fitted coefficients and the published CASTOR aggregates only.

## How to run

**Prerequisite: `config.do`.** The data-prep scripts read the restricted MRDR extract via `$data_dir`, set in a git-ignored `config.do` at the repo root. Copy `config.example.do` → `config.do` and set your paths once. (Without it, `$data_dir` is empty and the data-load fails with a plain "file not found".) Run everything from the repository root.

1.  **Fit the transport coefficients** with `outcomes/calibrated_transport.do` (separate from the dispatcher). It fits the ordered logit that feeds `outcomes/sim_bcr_override.do`:
    - `do "analyses/transport_dvd/outcomes/calibrated_transport.do" 0` gives the deterministic point estimate; it also writes the calibration-cohort and DVd baseline tables, the transport coefficients and the observed benchmarks into `results/` (these tabulations run **only** on this `0` call).
    - `do "analyses/transport_dvd/outcomes/calibrated_transport.do" 1 <min> <max>` runs the bootstrap.
2.  **Run a scenario.** `simulate.do` runs **one scenario** (one or two arms) per invocation. Config globals at the top set the run (`$scenario`, `$coeffs` = `transport_dvd`, `$int1`/`$int0`, `$line`, bootstrap flags, `$report`); `$scenario` is also overridable as the 4th positional arg (`do "analyses/transport_dvd/simulate.do" <boot> <min_bs> <max_bs> <scenario>`), which is how `run.do` loops the three scenarios. The pipeline is `load_patients → mata_setup → simulation → process_data → export_results`, then validation and (when `$report = 1`) the PDF. CSV export runs by default.
3.  **Aggregate.** Run all three scenarios, then the cross-scenario aggregation `simulated/bootstrap_summary.do` (below).

The full ordered pipeline (deterministic locally, plus the bootstrap HPC plumbing) is recorded in `run.do`.

## Repository structure

```         
analyses/transport_dvd/
├── README.md                        # this file
├── run.do                           # analysis runbook: deterministic (local) + bootstrap (HPC) pipeline
├── simulate.do                      # dispatcher: runs ONE scenario (deterministic + bootstrap)
├── draw_figure1.R                   # response figure (manuscript Figure 2) from results/
├── draw_figure1_poster.py           # poster version of the same figure, same CSVs
├── outcomes/
│   ├── calibrated_transport.do      # transport generator: `0` = deterministic, `1` = bootstrap
│   ├── calibrated_transport_prop.do # partial proportional-odds diagnostic (record only)
│   ├── sim_bcr_override.do          # L2 BCR override, branching on scenario and arm
│   ├── sim_txr_override.do          # forces the L2 regimen per arm
│   ├── txr_transport_dvd.do         # per-line regimen lists for the coefficient set
│   ├── A_trial/                     # CASTOR aggregate response distributions (.mmat)
│   └── B_transport/                 # fitted transport coefficients (.mmat)
├── coefficients/                    # coefficients_transport_dvd.mmat (bootstrap/ is local)
├── patients/                        # cohort .dta (git-ignored)
│   └── sample_size/                 # cohort and Monte Carlo precision scripts (git-ignored, local)
├── results/                         # cross-scenario CSVs, figures, results.md (git-ignored, local)
└── simulated/
    ├── bootstrap_summary.do         # cross-scenario aggregation of the bootstrap output
    ├── A_trial/                     # per-scenario .dta + CSVs (git-ignored)
    ├── B_transport/
    ├── C_mrdr/
    └── report/                      # per-run PDF reports
```

## Outputs

Output follows the model-wide convention in the root `README.md` (Result Exports).

**Engine-level, per run (`core/`).** `core/export_results.do` writes the CSVs every analysis needs (per-patient summary, BCR distribution, mean cost/QALY/LY) into `simulated/<scenario>/`. It runs by default as part of the pipeline (immediately after `process_data`, once per arm; skipped during bootstrap). It emits `bcr_<stub>.csv`, `econ_<stub>.csv`, and `patients_<stub>.csv` (where `<stub> = <int>_<line>_<data>_<min_id>_<max_id>`).

**Cross-scenario aggregation (`simulated/bootstrap_summary.do`).** Run once after the 6×500 bootstrap exists. All three scenarios are bootstrapped (A = resampled CASTOR BCR, B = Calibrated Transport bootstrap, C = the registry equation re-estimated on each resampled registry), so every interval reflects the same sources and is comparable. It reads `simulated/<scenario>/bootstrap/{dvd,vd}_2_predicted_B<b>.dta`, computes MAE for A and B against the **fixed raw observed DVd benchmark (n = 399)** each iteration, paired within replicate, the side-by-side DVd-versus-Vd ICER, incremental cost and incremental QALY intervals, and the per-scenario BCR distributions. This step cannot live in a single run, since each scenario is run separately.

### Read surface: `results/`

`results/` is the canonical place a downstream consumer reads from. **It is git-ignored**: these are manuscript and appendix outputs, regenerated locally.

```         
results/
├── bcr_distributions.csv    # predicted vs observed BCR × 6 categories, all scenarios
├── mae_comparison.csv       # MAE per scenario with bootstrap 95% intervals
├── icer_comparison.csv      # DVd vs Vd ICER, per scenario
├── transport_dvd_coefs.csv  # transport regression coefficients
├── mrdr_vd_baseline.csv     # calibration-cohort baseline (from `calibrated_transport.do 0`)
├── mrdr_dvd_baseline.csv    # DVd validation-cohort baseline
├── observed_dvd_bcr.csv     # observed DVd benchmark
├── observed_vd_bcr.csv      # observed Vd response in the calibration cohort
└── results.md               # narrates and labels the key figures
```

Intermediate per-scenario CSVs stay in `simulated/<scenario>/`; only the final cross-scenario CSVs plus `results.md` belong in `results/`.

## Implementation notes

- `core/generate_report.do` computes the BCR-by-regimen, cost, QALY and ICER quantities the CSVs need, but as `putpdf` tables and with older variable names (`cTotald`, `qTotald`) that have drifted from `process_data.do` (`cost_total_d`, `qaly_total_d`). The CSV export reads the actual `process_data` variables.
- The dispatcher's inline ICER calculation only fires in the two-arm case; the cross-scenario ICER comparison belongs in `bootstrap_summary.do`, not the dispatcher.
- The cohort is fixed at N = 50,000 per arm by a Monte Carlo precision rule. The scripts that build the line-entry pool, draw the cohort and size it (`patients/sample_size/`) are kept local and are not in the repository; the method is described in the manuscript's Supplementary Information and in `docs/technical_review.md`.
- Six of the 500 bootstrap replicates failed at the risk-equation stage and were replaced by spare replicates generated above index 500 (`hpc/fill_gaps.sh`).

## Results

**Matched cohort N = 50,000 per arm; DVd vs Vd, discounted; 2025 cost year; line-adjusted model with common random numbers; deterministic plus 500-replicate bootstrap.** Source: `results/` via `bootstrap_summary.do` (run 10 August 2026).

- **BCR accuracy** (benchmark = recorded second-line DVd, n = 399, fixed): bootstrap MAE 9.35 → 4.12 percentage points; mean paired reduction 56.3% (95% interval 16.2 to 83.1%); Calibrated Transport closer to observed in 99.4% of replicates.
- **Economics:** ICERs for A / B / C of A$708,306 / A$1,016,845 / A$1,386,469 (incremental QALYs 0.272 / 0.177 / 0.126, all intervals excluding zero). C is the fitted registry-based reference, not observed lifetime economics.
- All three scenarios are bootstrapped (500 replicates) for comparability.
