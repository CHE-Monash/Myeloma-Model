# prep/ — MRDR → model inputs

Builds the model's inputs from the Australia & New Zealand Myeloma and Related Diseases Registry (MRDR). These scripts run against the **restricted** registry data on the shared university drive (located per machine via `$data_dir` in `config.do`) and are kept local (git-ignored) — they can't run without MRDR access. Their **outputs** (the coefficient `.mmat` sets, the population cohorts, and the validation benchmarks) are what the rest of the model consumes and what ship.

Run everything from the repository root. The data cut is currently fixed by a hardcoded `local Date 251128` (28 Nov 2025) in each script.

## Pipeline

```         
data_extraction.do      raw MRDR tbl_*.dta (registry drive)      ->  "MRDR Long.dta"
multiple_imputation.do  "MRDR Long.dta"                          ->  "MRDR Long MI.dta"
                                                                     (+ "MRDR Wide MI.dta",
                                                                      + bootstrap/MRDR Long MI B<b>.dta)
    ├─ risk_equations.do       "MRDR Long MI.dta" + txr_<coeffs>.do -> analyses/<a>/coefficients/coefficients_<coeffs>.mmat
    ├─ generate_benchmarks.do  "MRDR Long MI.dta"                   -> scratch/benchmarks/*.csv  (18 files)
    └─ population_1995_2040.do  Forecast.xlsx + "MRDR Wide MI.dta"  -> patients/population_1995_2040_1..10.dta
```

Steps 1→2 are sequential; once `MRDR Long MI.dta` exists, `risk_equations` / `generate_benchmarks` / `population_1995_2040` are independent.

## Scripts

**`data_extraction.do`** — reads the raw MRDR registry tables (`tbl_Diagnosis`, `tbl_Patient`, `tbl_Chemotherapy`, `tbl_MaintenanceTherapy`, `tbl_Review`, `tbl_ASCT`, …), cleans/labels each via the `sub/` helpers, and reshapes into a long event-history with derived lines of therapy, BCR, ISS/R-ISS, treatment durations (TXD), treatment-free intervals (TFI) and model covariates → **`MRDR Long.dta`**. No arguments; run standalone.

**`multiple_imputation.do`** — multiply-imputes missing diagnosis covariates and BCR in `MRDR Long.dta` → **`MRDR Long MI.dta`** (plus `MRDR Wide MI.dta`, the per-imputation covariate seed for the population cohorts). Args: `imp` (number of imputations), `boot` (`0` = single dataset; `1` = HPC bootstrap branch), `min_bs`, `max_bs`. The bootstrap branch writes `bootstrap/MRDR Long MI B<b>.dta` per resample.

**`risk_equations.do`** — fits the 50 risk equations (OS, ASCT, regimen multinomials, BCR ordered logits, TXD/TFI parametric survival, by line) from `MRDR Long MI.dta` → a Mata coefficient set **`analyses/<analysis>/coefficients/coefficients_<coeffs>.mmat`**. Invoked per analysis from that analysis's `run.do`. Args (8): `analysis`, `coeffs`, `min_year`, `max_year`, `boot`, `min_bs`, `max_bs`, and an optional `sample` (`""` = standard MI data; `train`/`test` select the OOS fold). It loads the analysis's **regimen definition** via `do "analyses/$analysis/outcomes/txr_$coeffs.do"` (below).

**`generate_benchmarks.do`** — extracts observed validation targets (OS, BCR, TXD, TFI, pathways) from `MRDR Long MI.dta` → **18 benchmark CSVs** (default `scratch/benchmarks/`; the OOS analysis redirects them to `analyses/oos/targets/`) that the shared comparison engine `analyses/oos/validate_outcomes.do` checks the simulation against. Args (optional): `bench_in`, `bench_out`, `sample`.

**`population_1995_2040.do`** — generates the synthetic incident-MM simulation cohorts (incidence 1995–2020 AIHW + 2021–2040 Daffodil Centre, covariates seeded from `MRDR Wide MI.dta`) → **`patients/population_1995_2040_1..10.dta`**, the cohorts `core/load_patients.do` reads for `$data = population`. No arguments.

## Regimen definitions — `txr_<coeffs>.do`

`risk_equations.do` learns which regimens to model per line by running `analyses/$analysis/outcomes/txr_$coeffs.do`. For example `analyses/base_model/outcomes/txr_base_model.do` maps MRDR `Regimen` codes onto `TXR_L1`–`TXR_L9` (everything else → 0 = "other"): **L1** VCd(4)/VRd(31); **L2** Rd(7)/DVd(80); **L3** Kd(49)/Rd(7); **L4** Kd(49)/Pd(56); **L5–L9** unset (all "other"). Each coefficient set has a matching `txr_<coeffs>.do` defining its regimen exposure.

## `sub/` — registry cleaning helpers

`do`-ne by `data_extraction.do`:

- **`sub/TRU/`** — the registry team's generic label/clean code for the raw `tbl_*` tables (`Clean*` / `Label*` for Diagnosis, Patient, Chemotherapy, Maintenance, Review; `CalcDiagnostics.do` derives ISS / R-ISS / FISH risk).
- **`sub/MRDR/`** — project-specific derivations: `MRDR Chemo Regimen.do` (standardised regimen codes), `MRDR Maint Regimen.do`, `MRDR EQ5D AU.do` (Australian EQ-5D utility), `MRDR Manual Data Cleaning.do` (per-patient fixes). `MRDR/_archive/` holds superseded copies.

## Downstream consumers

- **Coefficients** (`coefficients_<coeffs>.mmat`) → loaded by each analysis dispatcher via `mata matuse "$coefficients_path/coefficients_$coeffs"` (`base_model`, `transport_dvd`, `oos`, …).
- **Population cohorts** → `core/load_patients.do` (`use "patients/population_1995_2040_<n>.dta"`).
- **Benchmarks** → `analyses/oos/validate_outcomes.do` (imports the target CSVs as comparison matrices).

## Paths & config

All MRDR/machine paths go through the git-ignored `config.do` at the repo root (see `config.example.do`). Each prep script does `capture run "config.do"` and reads:

- **`$data_dir`** — the dated EpiMAP working dir (`MRDR Long.dta`, `MRDR Long MI.dta`, `MRDR Wide MI.dta`, `bootstrap/`). Used by all five scripts.
- **`$mrdr_raw_dir`** — the raw MRDR registry source (the `tbl_*.dta` tables; a different drive branch). Used by `data_extraction.do` and `sub/TRU/CalcDiagnostics.do`.
- **`$epimap_dir`** — the EpiMAP project base on the drive (`$data_dir = ${epimap_dir}/Data/${data_cut}`).
- **`$data_cut`** — the data-cut date (e.g. `251128`); **`$scratch_dir`** — bootstrap/scratch output (HPC).

Run from the repository root. `data_extraction.do` sources the `sub/` helpers from the repo (`prep/sub/…`), not the drive; the stale `cd`s in `generate_benchmarks.do` / `population_1995_2040.do` are gone (outputs resolve relative to the repo root); `CalcDiagnostics.do` now uses `$mrdr_raw_dir` (was pinned to an older `2024/241129` cut); and `multiple_imputation.do`'s working `temp/` is created relative to the current directory.

> `population_1995_2040.do` reads its incidence inputs from `patients/population_forecast.csv` (Daffodil projections, 2010–2043) and `patients/population_historical.csv` (AIHW, 1995–2020) — both git-ignored. Imported with `case(preserve)` so the `Sex/AgeGroup/Year/Incidence` names survive (the files carry a UTF-8 BOM, which Stata strips on import).
