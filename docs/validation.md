# Monash Myeloma Model — Validation

## Overview

Validation compares a base-model simulation against benchmarks derived from the MRDR registry to confirm the model reproduces observed outcomes within agreed tolerances. Benchmarks are CSVs held in `validation/benchmarks/`, generated from the imputed MRDR data by `prep/generate_benchmarks.do` (run against the restricted registry data — see `prep/README.md`).

## How to run

Validation reads the base-model output, so run the base-model dispatcher first, then the validation script — both from the repository root:

```stata
// 1. Produce the base-model simulation
do "analyses/base_model/base_model.do"
//    → analyses/base_model/simulated/all_0_population_1_101212.dta

// 2. Validate it against the MRDR benchmarks
do "validation/validate_simulation.do"
```

`validate_simulation.do` imports the benchmark CSVs inline (there is no separate loader), loads the simulated dataset, runs the checks below, and prints a pass/fail summary.

## Files

| File | Role |
|---|---|
| `validation/validate_simulation.do` | Headline acceptance test: simulation vs MRDR benchmarks |
| `validation/benchmarks/*.csv` | Pre-baked MRDR benchmark targets (13 files) |
| `core/validation.do` | Lighter invariant checks run inside every simulation (distinct from the acceptance test) |
| `validation/validate_vectors.do` | Developer check: vectorised representation reproduces the old matrix representation |
| `validation/test_mata_functions.do` | Unit check: ordered-logit and survival helpers |
| `validation/test_survival_functions.do` | Unit check: exponential / Weibull / Gompertz inverse-CDF sampling |
| `validation/test_SCT_DN.do` | Loop-vs-vector equivalence check for ASCT-at-diagnosis |

## Benchmarks

`validation/benchmarks/` contains 13 CSVs: `os_l1_noasct.csv`, `os_asct.csv`, `os_l2.csv`, `os_l3.csv`, `bcr.csv`, `txd_l1_noasct.csv`, `txd_l1_asct.csv`, `txd_l2.csv`, `tfi_l1_noasct.csv`, `tfi_l1_asct.csv`, `tfi_l2.csv`, `tfi_l3.csv`, and `pathways.csv`. Overall-survival files carry N, median and annual survival percentages; response carries N and the CR/VGPR/PR/MR/SD/PD percentages by line; treatment-duration and treatment-free-interval files carry N, mean, median and quartiles; pathways carries the ASCT and subsequent-line reach rates.

The survival, treatment-duration and treatment-free-interval benchmarks are all estimated with censoring-aware survival methods (Kaplan–Meier / `stsum`), so they are comparable to the run-to-death simulation despite the registry's incomplete follow-up. The **subsequent-line reach rates** in `pathways.csv` (L2–L8) are estimated as a **competing-risks cumulative incidence** — the Aalen–Johansen probability of reaching each line by end of follow-up, with *death before reaching it* as the competing event and L1 start as the origin. This estimates the lifetime reach probability, which is what the run-to-death simulation reports; a crude "ever reached ÷ total" count would understate it, because recently-diagnosed registry patients who are still alive in an earlier line would usually progress given more follow-up. The ASCT reach rate is left as a crude proportion (transplant is decided during L1 and is observed for essentially every patient, so it is not subject to that cumulative follow-up censoring).

## Checks and tolerances

`validate_simulation.do` runs five families of checks and counts passes and failures:

| Family | Metric | Tolerance |
|---|---|---|
| Overall survival | 3-year and 5-year survival %, by response and line (L1 no-ASCT, ASCT, L2, L3) | within ±10 percentage points |
| Best clinical response | per-category % by line (L1, ASCT, L2, L3) | within ±5 percentage points |
| Treatment duration (TXD) | median (simulated ÷ benchmark) | within ±20% |
| Treatment-free interval (TFI) | median (simulated ÷ benchmark), plus response ordering (CR > VGPR > PR > MR > SD) at L2 | within ±20% |
| Pathways | ASCT and L2–L5 reach rates | within ±5 percentage points |

The wider tolerances on treatment duration and treatment-free interval reflect parametric underestimation against censored training data and sparser cells in poor-responder and later-line groups.

## Interpreting results

The script ends with a summary of tests run, passed and failed. As a guide: a pass rate above 90% indicates the model reproduces the benchmarks well; 75–90% warrants reviewing the failed checks; below 75% indicates a structural issue to investigate before relying on the run.

## Developer equivalence checks

`validate_vectors.do` and the `test_*.do` scripts are standalone developer tools, not part of the acceptance test. They assert that the vectorised implementation reproduces the previous loop-based results and that the Mata helpers match their analytic formulae, at machine tolerance. They are run individually as needed during engine changes.

> Note: `validation/test_SCT_DN.do` currently loads coefficients from a path that no longer exists (`analyses/base_model/data/coefficients/…`); it needs updating to the current `analyses/base_model/coefficients/` layout before it will run.

## When to refresh benchmarks

Regenerate the benchmark CSVs (re-run `prep/generate_benchmarks.do`) when the training data are updated, the imputation strategy changes, variable definitions change, or the model structure changes (for example, the response categories). Commit the refreshed benchmarks so model performance can be tracked across versions and regressions detected.

> Note: the pathways benchmark estimator was changed from a crude reach count to a competing-risks cumulative incidence (see *Benchmarks* above) so that the L2–L8 targets reflect lifetime reach probability and are comparable to the run-to-death simulation. The CSV schema is unchanged, so `validate_simulation.do` is unaffected. Re-running `generate_benchmarks.do` is required for the new pathways targets to take effect.
