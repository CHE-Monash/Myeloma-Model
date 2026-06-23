**********
* EpiMAP Myeloma - DVd L2 Cohort Pool (build once, reuse)
*
* Purpose: Build the reusable second-line ENTRY POOL for the DVd-vs-Vd analysis:
*          all patients (across the independent population samples) who are alive
*          at end of L1 and reach L2 within the case-mix window. This is the
*          expensive step - it simulates each incident population from diagnosis.
*          Build it once; downstream scripts draw cohorts from it without
*          re-simulating:
*            - ce_cohort.do       : draws the production cohort (size N) from the pool
*            - ce_precision.do    : the convergence / sigma_pp study (reads the pool)
*            - transport_dvd.do   : production runs read the drawn cohort
*
* Case-mix: defined by the window (WHO reaches L2), held fixed. Do not widen the
*           window to change size - that shifts case-mix. Size is a downstream
*           sampling choice, not a property of the pool.
*
* Independent samples: loops data="population_`s'" over the n_samples independent
*           population files (NOT data="population", which always loads file 1).
*
* Output: analyses/transport_dvd/patients/cohort_pool_<line>.dta
*
* Usage:  do cohort_pool.do
*
* Author: Adam Irving
* Date: June 2026
**********

clear all
set more off
macro drop _all

cd "/Users/adami/Documents/Monash/Vault/research/models/myeloma model/repo"

**********
* Configuration
**********

local n_samples  = 10          // independent population files to pool (max 10)
local all_l2     = 0           // 1 = ALL L2-reachers (large DISTINCT pool for the sample-size study)
                               // 0 = restrict to start_year-end_year (production case-mix)
local start_year = 2020        // case-mix window, used only when all_l2 == 0
local end_year   = 2025

global analysis  "transport_dvd"
global line      "2"
global coeffs    "dvd_post"
global cost_year "2020"
global drate     "0.05"
global boot      "0"

global coefficients_path  "analyses/$analysis/coefficients"
global outcomes_path      "analyses/$analysis/outcomes"
global patients_path      "analyses/$analysis/patients"
global simulated_path     "analyses/$analysis/simulated"

// all-L2 pool (sample-size study) saved separately from the windowed production pool
if (`all_l2') local pool_file "$patients_path/cohort_pool_all_${line}.dta"
else          local pool_file "$patients_path/cohort_pool_${line}.dta"

**********
* Load core programs (incl. the shared engine pass)
**********

run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"
run "core/run_pipeline.do"

**********
* Build pool: simulate each independent population -> L2 entry -> filter
**********

di as text _n "=== Building cohort pool from `n_samples' independent populations ==="

forval s = 1/`n_samples' {

	mata: mata clear

	global int      "dvd"            // natural-history setting only; L2 regimen blanked below
	global data     "population_`s'" // INDEPENDENT population file s
	global min_year "1995"
	global max_year "2020"
	global min_id   "1"
	global max_id   "101212"

	qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"

	run_pipeline

	// Filters: alive at end of L1; restrict to the case-mix window only when all_l2 == 0
	keep if MOR_L`= ${line} - 1'E == 0
	if (`all_l2' == 0) keep if YearL${line} >= `start_year' & YearL${line} <= `end_year'

	// Clean for simulation (blank the L2 regimen so the cohort is arm-agnostic)
	replace State = ${line} * 2
	replace Age_L${line}S = .
	replace TXR_L${line} = .
	replace ID = _n
	replace DateDN = td(1jan2020) - (TSD_L${line}S * 12)
	replace YearDN = yofd(DateDN)
	drop DateL* DateMOR YearL* YearMOR c* q* OC_TIME_L TSD_*_ref
	cap drop DateSCT YearSCT

	gen Sample = `s'
	save "$patients_path/pool_part_`s'.dta", replace
}

// Combine the parts into the pool
use "$patients_path/pool_part_1.dta", clear
erase "$patients_path/pool_part_1.dta"
forval s = 2/`n_samples' {
	append using "$patients_path/pool_part_`s'.dta", nolabel
	erase "$patients_path/pool_part_`s'.dta"
}

replace ID = _n
save "`pool_file'", replace

di as text "Pool built: " as result _N as text " L2-entry patients -> `pool_file'"
tab Sample, missing
summarize TSD_L${line}S, detail
