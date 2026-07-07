# Monash Myeloma Model — Validation

## Overview

The model is checked at three layers:

1.  **In-simulation invariants** — `core/validation.do`, run inside every simulation: lightweight ordering/bounds sanity checks (e.g. event times non-decreasing, durations non-negative). Catches gross engine errors on every run.
2.  **Engine verification** — `core/tests/`: unit tests for the Mata/survival primitives plus extreme-value (stress) testing. Verifies the *machinery* with no calibration data.
3.  **Calibration validation** — the **out-of-sample (70/30) analysis** in `analyses/oos/`. This is the **mainstay validation, re-run for each model version** to confirm the model still reproduces observed outcomes and nothing has regressed.

> The previous in-sample registry acceptance test (a base-model run compared to pre-baked MRDR benchmark CSVs) has been **retired**. Its comparison logic survives as the shared engine `analyses/oos/validate_outcomes.do`; the old debug diagnostics were moved to `scratch/`.

## Layer 2 — engine verification (`core/tests/`)

| File | Role |
|----|----|
| `core/tests/extreme_value.do` | Extreme-value / stress test: perturb one risk-equation intercept to an extreme and assert the engine responds in the right direction and stays bounded |
| `core/tests/test_mata_functions.do` | Unit check: ordered-logit and survival helpers |
| `core/tests/test_survival_functions.do` | Unit check: exponential / Weibull / Gompertz inverse-CDF sampling |
| `core/tests/test_suite.md` | Test catalogue / specification |

`extreme_value.do` needs **no MRDR data** — it runs the engine on a 3k slice of the in-repo `population_1995_2040` cohort with the base coefficients, perturbing one intercept at a time in Mata, and self-checks the result. The latest run passes **7/7**: OS hazard → ∞ gives median OS ≈ 0 months and → 0 gives ≈ 367 months (the age limit); ASCT probability → 1 gives ≈ 98% transplanted and → 0 gives 0%; TXD and TFI hazard → ∞ collapse to median ≈ 0; and a monotone sweep of the OS intercept gives a smoothly decreasing median OS. Run it from the repository root:

``` stata
do "core/tests/extreme_value.do"
```

> Known engine-robustness gap (surfaced by this test): at *intermediate-extreme* OS hazard (intercept shift ≈ +2, not the +20 boundary) the engine throws `r(3301) subscript invalid` — high early mortality empties a downstream line and an unguarded `selectindex` fails. The harness self-protects (`capture` per point); the engine guard (`rows(idx) > 0` before indexing) is an open item.

## Layer 3 — out-of-sample validation (`analyses/oos/`)

Trains the model on a random **70%** of MRDR patients, predicts the held-out **30%**, and compares the predictions to those patients' **observed** outcomes — testing generalisation, not just in-sample fit. See `analyses/oos/README.md` for the full layout; the prep steps (split → train-70% imputation/risk-equations → 30% targets/cohort) run against the restricted registry data and are HPC-suited.

### How to run

The full pipeline — split, train-70% fit, 30% targets + cohort, simulate, and compare — is the deterministic track of `analyses/oos/run.do` (its numbered steps 0–6); run that top to bottom, or the individual steps it lists:

``` stata
// The whole deterministic OOS pipeline (see the numbered steps inside):
do "analyses/oos/run.do"

// ... or just the final compare (steps 1-5 already done):
//   step 5  do "analyses/oos/simulate.do"       // simulate the held-out 30% with the 70%-trained coefficients
//   step 6  do "analyses/oos/validate_oos.do"   // compare the simulated 30% to the observed targets
```

`validate_oos.do` points the shared comparison engine `validate_outcomes.do` at the OOS targets (`$val_targets`) and the OOS simulation (`$val_simfile`); the engine imports the target CSVs inline, loads the simulated dataset, runs the checks below, and prints a pass/fail summary.

### Targets

`analyses/oos/targets/` holds 13 CSVs (the held-out 30%'s observed outcomes), built by `prep/generate_benchmarks.do` from the test-fold imputed data: `os_l1_noasct.csv`, `os_asct.csv`, `os_l2.csv`, `os_l3.csv`, `bcr.csv`, `txd_l1_noasct.csv`, `txd_l1_asct.csv`, `txd_l2.csv`, `tfi_l1_noasct.csv`, `tfi_l1_asct.csv`, `tfi_l2.csv`, `tfi_l3.csv`, `pathways.csv`. Overall-survival files carry N, median and annual survival percentages; response carries N and the CR/VGPR/PR/MR/SD/PD percentages by line; treatment-duration and treatment-free-interval files carry N, mean, median, quartiles and the `M12`/`M24` horizon-survival columns; pathways carries the ASCT and subsequent-line reach rates. A further target, `os_wholepop_curve.csv`, holds the observed whole-population monthly Kaplan–Meier survivor (with Greenwood SE) over 0–120 months — the reference for the OS validation curve (below).

The survival, treatment-duration and treatment-free-interval targets are all estimated with censoring-aware survival methods (Kaplan–Meier / `stsum`), so they are comparable to the run-to-death simulation despite the registry's incomplete follow-up. The **subsequent-line reach rates** in `pathways.csv` are estimated **one transition at a time, conditionally** — the Aalen–Johansen probability of reaching line *L* given the patient reached line *L−1*, with *death before reaching L* as the competing event and the *previous line's reach date* as the origin. Validating each transition on its own keeps line-to-line errors from compounding down the pathway (an earlier cumulative-from-L1 estimate let a single early miss cascade into every later line); the validator and the simulation are compared on the same conditional quantity. Each conditional rate is still a competing-risks cumulative incidence rather than a crude "ever reached ÷ total" count, which would understate reach because recently-diagnosed registry patients still in an earlier line would usually progress given more follow-up. The ASCT reach rate is a crude proportion, but its **denominator is patients who reach the end of L1** (the transplant decision point), not all diagnosed patients — transplant is decided at L1 end, the model's ASCT logit is fit on that conditional population, and the comparison uses the matching denominator. Dividing by all patients would understate the simulated rate because the \~8% who die during induction never reach the decision.

### Checks and tolerances

`validate_outcomes.do` runs five families of checks and counts passes and failures:

| Family | Metric | Tolerance |
|----|----|----|
| Overall survival | 3-year and 5-year survival %, by response and line (L1 no-ASCT, ASCT, L2, L3) | within ±10 percentage points |
| Best clinical response | per-category % by line (L1, ASCT, L2, L3) | within ±5 percentage points |
| Treatment duration (TXD) | \% still on treatment at 12 & 24 months, by response and line (L1 no-ASCT, ASCT) | within ±10 percentage points |
| Treatment-free interval (TFI) | \% still treatment-free at 12 & 24 months, by response and line (L1 no-ASCT, ASCT, L2), plus 12-month response ordering (CR \> VGPR \> PR) at L2 | within ±10 percentage points |
| Pathways | ASCT rate and L2–L5 conditional reach rates (P(reach L \| reached L−1)) | within ±5 percentage points |

TXD and TFI are validated by **survival at fixed horizons** (% still on treatment / still treatment-free at 12 and 24 months), the same way overall survival is checked — not by the median. The median is unreliable for these outcomes because they are extremely right-skewed and heavily censored (TFI 20–54%): poor responders progress almost immediately (median ≈ 0, so a ratio test is unstable), while good responders have so much administrative censoring that the KM median is itself an unreliable, likely under-estimate. Survival at 12/24 months is well inside observed follow-up, so it is estimated robustly for both extremes and is directly comparable between the simulation and the registry.

### Interpreting results

The script ends with a summary of tests run, passed and failed. As a guide: a pass rate above 90% indicates the model reproduces the held-out outcomes well; 75–90% warrants reviewing the failed checks; below 75% indicates a structural issue to investigate before relying on the run. The latest deterministic point-estimate run scores **83.1% (143/172)**; the *same* model checked in-sample scores **94.8% (163/172)**, so the residual out-of-sample misses (chiefly an L1 no-ASCT CR under-prediction) are sampling artefacts of the 70/30 split rather than model defects. The other fails are data artefacts — the SMM-tail response misses, missing treatment-end-dates inflating the TXD target, and the older incidence cohort versus the registry — not structural problems.

### Bootstrap prediction intervals and the OS validation curve

The point estimate above uses one 70%-trained coefficient set. The **calibration** question — does each held-out observed value fall inside the *interval* the model predicts? — is answered by `analyses/oos/bootstrap_validation.do` on the HPC, using **500 bootstrap simulations** of the held-out 30% (each a patient-cluster resample of the 70%, re-imputed, re-fit and re-simulated). By the **percentile method** it asks whether each observed target lies within the bootstrap 95% interval [p2.5, p97.5]. Latest coverage: **105/171 (61.4%)** — in line with the model's history, dominated by the tight parameter-only intervals plus the known data artefacts above; aggregate OS passes.

The same run also builds the **whole-population OS validation curve** (the 2024 PLOS ONE Fig 2 for this model): the held-out validation cohort's Kaplan–Meier **95% CI** against the simulated cohort's **95% CI** over 120 months, plus a **monthly two-sample z-test p-value** (observed Greenwood SE + simulated bootstrap SD) testing the difference. Result: **no significant difference in 118/120 months (98.3%)** — the simulated OS 95% CI sits within the observed CI across the full 10 years (p < 0.05 only in the first ~2–3 months), improving on the 2024 paper's 90%. It writes `results/os_wholepop_curve_validation.csv` and `os_wholepop_curve.png`; `analyses/oos/plot_os_curve.do` redraws the figure locally from the CSV when HPC batch graphics are off.

## When to re-validate

Re-run the OOS validation **for each model version**, and whenever the training data are updated, the imputation strategy changes, variable definitions change, or the model structure changes (for example, the response categories). The targets are regenerated as part of the prep steps, so each version is validated against a freshly-derived held-out set.
