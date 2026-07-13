# prep/ — model inputs

Builds the model's inputs. Two independent groups, with **different access requirements**:

1. **Clinical inputs (MRDR-gated)** — risk equations, synthetic cohorts and validation benchmarks, built from the Australia & New Zealand Myeloma and Related Diseases Registry (MRDR), starting from `MRDR Long.dta` (the cleaned long event-history). These run against the **restricted** registry data on the shared university drive (located per machine via `$data_dir` in `config.do`) — they need MRDR access to run, though the code ships. The upstream step that builds `MRDR Long.dta` from the raw `tbl_*` registry tables is a separate, restricted extraction layer kept **local (git-ignored)**. Their **outputs** (the coefficient `.mmat` sets, the synthetic cohorts, and the validation benchmarks) are what the rest of the model consumes and what ship.

2. **Cost inputs (public; no MRDR needed)** — see [Cost inputs](#cost-inputs) below. These build from published sources (the PBS Schedule, ABS CPI, Yap 2025), so **anyone with the repo can run them and reproduce the cost engine end to end**. Both the scripts and their outputs are committed.

Run everything from the repository root. The MRDR data cut is currently fixed by a hardcoded `local Date 251128` (28 Nov 2025) in each of the clinical scripts.

## Pipeline

```         
# "MRDR Long.dta" is built upstream from the raw registry by a local (git-ignored) extraction step.
multiple_imputation.do  "MRDR Long.dta"                          ->  "MRDR Long MI.dta"
                                                                     (+ "MRDR Wide MI.dta",
                                                                      + bootstrap/MRDR Long MI B<b>.dta)
    ├─ risk_equations.do       "MRDR Long MI.dta" + txr_<coeffs>.do -> analyses/<a>/coefficients/coefficients_<coeffs>.mmat
    ├─ generate_benchmarks.do  "MRDR Long MI.dta"                   -> scratch/benchmarks/*.csv  (18 files)
    └─ synthetic_1995_2040.do  Forecast.xlsx + "MRDR Wide MI.dta"  -> patients/synthetic_1995_2040_1..10.dta
```

Steps 1→2 are sequential; once `MRDR Long MI.dta` exists, `risk_equations` / `generate_benchmarks` / `synthetic_1995_2040` are independent.

## Scripts

**`multiple_imputation.do`** — multiply-imputes missing diagnosis covariates and BCR in `MRDR Long.dta` → **`MRDR Long MI.dta`** (plus `MRDR Wide MI.dta`, the per-imputation covariate seed for the population cohorts). Args: `imp` (number of imputations), `boot` (`0` = single dataset; `1` = HPC bootstrap branch), `min_bs`, `max_bs`. The bootstrap branch writes `bootstrap/MRDR Long MI B<b>.dta` per resample.

**`risk_equations.do`** — fits the 50 risk equations (OS, ASCT, regimen multinomials, BCR ordered logits, TXD/TFI parametric survival, by line) from `MRDR Long MI.dta` → a Mata coefficient set **`analyses/<analysis>/coefficients/coefficients_<coeffs>.mmat`**. Invoked per analysis from that analysis's `run.do`. Args (8): `analysis`, `coeffs`, `min_year`, `max_year`, `boot`, `min_bs`, `max_bs`, and an optional `sample` (`""` = standard MI data; `train`/`test` select the OOS fold). It loads the analysis's **regimen definition** via `do "analyses/$analysis/outcomes/txr_$coeffs.do"` (below).

**`generate_benchmarks.do`** — extracts observed validation targets (OS, BCR, TXD, TFI, pathways) from `MRDR Long MI.dta` → **18 benchmark CSVs** (default `scratch/benchmarks/`; the OOS analysis redirects them to `analyses/default/targets/`) that the shared comparison engine `analyses/default/validate_outcomes.do` checks the simulation against. Args (optional): `bench_in`, `bench_out`, `sample`.

**`synthetic_1995_2040.do`** — generates the synthetic incident-MM simulation cohorts (incidence 1995–2020 AIHW + 2021–2040 Daffodil Centre, covariates seeded from `MRDR Wide MI.dta`) → **`patients/synthetic_1995_2040_1..10.dta`**, the cohorts `core/load_patients.do` reads for `$data = synthetic`. No arguments.

## Cost inputs

**No MRDR access required.** These build from public sources and are committed with their outputs, so the cost engine is reproducible by anyone with the repo. Derivations, prices and per-cycle figures: `docs/economic_inputs.md`.

```
build_cost_index.do          ABS CPI 6401.0            -> inputs/cost_index.csv
extract_pbs_costs.do         dated PBS Schedule extract-> inputs/pbs_{prices,fees,markups,copayments}_<year>.csv
extract_pbs_restrictions.do    "      "        "       -> inputs/pbs_restrictions.csv   (reference only)
    └─ treatment_costs.do    [year] [oral_policy] [net_copay]
                             + inputs/treatment_regimens.csv (dosing spec)
                             + inputs/other_costs.csv        (non-drug)
                                                       -> inputs/treatment_costs_<year>.csv
```

**`build_cost_index.do`** — the ABS CPI series (6401.0, All groups), hardcoded so the deflator lives in code → `inputs/cost_index.csv`. It inflates the **non-drug** costs from their source year to the target year. Drug prices are **not** inflated: they use actual dated PBS values.

**`extract_pbs_costs.do`** — filters a dated full PBS Schedule extract (`$pbs_src`, ~40 MB, git-ignored, in the sibling `data/` folder) down to the modelled myeloma drugs, and commits the small subset: AEMPs, fees, mark-up bands and co-payments per year. The raw extract stays out of the repo; the script records exactly how the committed subset was produced, so the chain is reproducible from a re-downloadable dated source.

**`extract_pbs_restrictions.do`** — the PBS restriction map (which regimens are reimbursed at which line) → `inputs/pbs_restrictions.csv`. An **eligibility reference only**; it is not wired into the cost calculation.

**`treatment_costs.do`** — builds per-cycle regimen drug costs from first principles as the PBS **Dispensed Price for Maximum Quantity (DPMQ)**, plus the inflated non-drug costs → **`inputs/treatment_costs_<year>.csv`**, which `core/process_data.do` reads at simulate time. Three optional args:

| Arg | Default | Meaning |
|---|---|---|
| `year` | `2026` | Target price year — selects the `pbs_*_<year>.csv` inputs. |
| `oral_policy` | `wholepack` | `wholepack` costs whole packs per cycle (ceiling to the pack boundary, so **dispensing wastage is costed** — the PBAC convention). `prorata` costs only units consumed, as a sensitivity. Injectables are always fewest-vials. |
| `net_copay` | off | Nets the patient co-payment off the DPMQ, giving a cost-to-government view. The model's default perspective is the **Australian health system**, so the full DPMQ is used with the co-payment **not** netted off. |

> **`$cost_year` is a trap.** If no cost file exists for the requested year, `process_data.do` silently falls back to the **latest** available file. 2026 generic price disclosure moved drug costs materially, so re-running an older analysis now can reprice it without saying so. **Check the cost year before re-running any near-submission analysis.**

Non-treatment (phase-of-care) costs come from Yap 2025 and are applied in `core/process_data.do`; the initial phase is netted of the transplant admission to avoid double-counting it against the ASCT cost.

## Regimen definitions — `txr_<coeffs>.do`

`risk_equations.do` learns which regimens to model per line by running `analyses/$analysis/outcomes/txr_$coeffs.do`. For example `analyses/default/outcomes/txr_full.do` maps MRDR `Regimen` codes onto `TXR_L1`–`TXR_L9` (everything else → 0 = "other"): **L1** VCd(4)/VRd(31); **L2** Rd(7)/DVd(80); **L3** Kd(49)/Rd(7); **L4** Kd(49)/Pd(56); **L5–L9** unset (all "other"). Each coefficient set has a matching `txr_<coeffs>.do` defining its regimen exposure.

## Downstream consumers

- **Coefficients** (`coefficients_<coeffs>.mmat`) → loaded by each analysis dispatcher via `mata matuse "$coefficients_path/coefficients_$coeffs"` (`default` loads `full`/`train`; `transport_dvd`; …).
- **Synthetic cohorts** → `core/load_patients.do` (`use "patients/synthetic_1995_2040_<n>.dta"`).
- **Benchmarks** → `analyses/default/validate_outcomes.do` (imports the target CSVs as comparison matrices).

## Paths & config

All MRDR/machine paths go through the git-ignored `config.do` at the repo root (see `config.example.do`). Each prep script does `capture run "config.do"` and reads:

- **`$data_dir`** — the dated EpiMAP working dir (`MRDR Long.dta`, `MRDR Long MI.dta`, `MRDR Wide MI.dta`, `bootstrap/`). Used by all four scripts above.
- **`$mrdr_raw_dir`** — the raw MRDR registry source (the `tbl_*.dta` tables; a different drive branch). Used by the local (git-ignored) registry-extraction step that produces `MRDR Long.dta`.
- **`$epimap_dir`** — the EpiMAP project base on the drive (`$data_dir = ${epimap_dir}/Data/${data_cut}`).
- **`$data_cut`** — the data-cut date (e.g. `251128`); **`$scratch_dir`** — bootstrap/scratch output (HPC).

Run from the repository root. The stale `cd`s in `generate_benchmarks.do` / `synthetic_1995_2040.do` are gone (outputs resolve relative to the repo root), and `multiple_imputation.do`'s working `temp/` is created relative to the current directory.

> `synthetic_1995_2040.do` reads its incidence inputs from `patients/population_forecast.csv` (Daffodil projections, 2010–2043) and `patients/population_historical.csv` (AIHW, 1995–2020) — both git-ignored. Imported with `case(preserve)` so the `Sex/AgeGroup/Year/Incidence` names survive (the files carry a UTF-8 BOM, which Stata strips on import).
