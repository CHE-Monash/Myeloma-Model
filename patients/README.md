# Synthetic Cohorts

## Overview

This directory contains 10 independent realisations of the synthetic incident-MM cohort — patients diagnosed in Australia between 1995 and 2040. These are base datasets that can be used across different analyses without modification. Each dataset contains 101,212 patients.

## Files

**Cohorts** (the simulation inputs, read by `core/load_patients.do`):

- `synthetic_1995_2040_1.dta` to `synthetic_1995_2040_10.dta` — the 10 cohort realisations. **Git-ignored** (regenerated locally by `prep/synthetic_1995_2040.do`), not shipped in the public repo.

**Incidence inputs** (read *by* `prep/synthetic_1995_2040.do` to build the cohorts — despite the name, these are *not* cohort files and must not be deleted with them):

- `population_historical.csv` / `.xlsx` — observed incidence, AIHW, 1995–2020.
- `population_forecast.csv` — projected incidence, Daffodil Centre, 2010–2043.

The superseded `population_1995_2040_*.dta` cohorts were deleted in July 2026: they carried the retired ordinal comorbidity score (`CMc`) and the unused `CM_LVR` / `CM_PNR` / `CM_MLG` flags, and their covariates came from an imputation model that included them. The `synthetic_*` files replace them; there is no `population` cohort token any more.

## Cohort Characteristics

Each cohort contains 101,212 patients with:

| Variable | Description | Range/Categories |
|----------|-------------|------------------|
| **ID** | Patient identifier | 1 to N |
| **Age** | Age at diagnosis | 18-100 years |
| **Male** | Sex indicator | 0=Female, 1=Male |
| **ECOGcc** | ECOG performance status | 0, 1, 2 |
| **RISS** | Revised International Staging | 1, 2, 3 |
| **CM_CKD** | Renal impairment (chronic kidney disease) flag | 0, 1 |
| **CM_CRD** | Cardiac comorbidity flag | 0, 1 |
| **CM_PLM** | Pulmonary comorbidity flag | 0, 1 |
| **CM_DBT** | Diabetes flag | 0, 1 |
| **State** | Initial disease state | Usually 1 |
| **DateDN** | Date of diagnosis | Date format |

## Usage in Analyses

### Option 1: Use specific population
```stata
use "patients/synthetic_1995_2040_3.dta", clear
```

### Option 2: Programmatic selection
```stata
local pop_number = 1
use "patients/synthetic_1995_2040_`pop_number'.dta", clear
```

### Option 3: Parameter-driven selection
```stata
// In analysis file, use global parameter
global PopulationNumber = 5
use "patients/synthetic_1995_2040_${PopulationNumber}.dta", clear
```

## Why ten realisations

Each file is a different random realisation of the same underlying demographic and clinical distributions. `$data = "synthetic"` loads file 1 and is the default for a single projection run (`analyses/default`, `$scenario ""`).

The other nine exist for the **line-entry cohort pools**: a line-specific decision analysis needs enough patients who *reach* line L within a fixed case-mix window, and one incident cohort does not supply them. `patients/cohort_pool.do` (in `template`, `car_t` and `transport_dvd`) loops `$data = "synthetic_1..10"` and pools the line-L entrants. Widen the sample count, never the case-mix window.

## Maintenance

- Regenerate (`prep/synthetic_1995_2040.do`) whenever the incidence inputs, the covariate imputation model or the cohort schema change. The schema must stay in step with `core/load_patients.do` and `core/mata_setup.do`.
- Regenerating changes the realisation, so any cohort pool built from these files must be rebuilt with them.
