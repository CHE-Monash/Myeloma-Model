**********
* EpiMAP Myeloma - DVd L2 Cost-Effectiveness Cohort
*
* Purpose: Draw the production cost-effectiveness cohort for the DVd-vs-Vd
*          analysis: a fixed-size, fixed-seed sample of the 2020-2025 windowed
*          L2-entry pool, saved as the canonical patient file the dispatcher and
*          PSA read.
*
*          Supersedes the calendar-linked BI-cohort approach for the CE
*          analysis. The 2020-2025 window now defines the CASE-MIX only; the
*          cohort SIZE is a Monte Carlo precision choice - N = 35,000 gives a
*          Monte Carlo SD of the incremental QALY equal to 5% of the parameter
*          SD (see results/ce_sample_size_appendix.md). The budget-impact count
*          is a separate question that still uses the full windowed pool.
*
* Prerequisite: the windowed pool exists, i.e. cohort_pool.do has been run
*               with all_l2 = 0 (start_year 2020, end_year 2025), producing
*               patients/cohort_pool_<line>.dta.
*
* Output: analyses/transport_dvd/patients/patients_transport_dvd_<line>.dta
*         (the canonical file load_patients reads).
*
* Usage:  do ce_cohort.do
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
global analysis "transport_dvd"
global line     "2"
local  target_n = 50000          // CE cohort size (5% MC rule; see ce_sample_size_appendix.md)
local  seed     = 20260616       // fixed -> reproducible cohort

global patients_path "analyses/$analysis/patients"
local pool_file   "$patients_path/cohort_pool_${line}.dta"           // 2020-2025 windowed pool (cohort_pool all_l2 = 0)
local cohort_file "$patients_path/patients_${analysis}_${line}.dta"  // canonical file (load_patients reads this)

**********
* Require the windowed pool (and never read and overwrite the same file)
**********
capture confirm file "`pool_file'"
if (_rc) {
    di as error "Windowed pool not found: `pool_file'"
    di as error "Run cohort_pool.do with all_l2 = 0 (2020-2025) first."
    exit 601
}
if ("`pool_file'" == "`cohort_file'") {
    di as error "pool_file and cohort_file are the same - point pool_file at the pool, not the output."
    exit 198
}

**********
* Draw the CE cohort
**********
use "`pool_file'", clear
local pool_n = _N
di as text "Windowed pool: " as result `pool_n' as text " L2-entry patients (`pool_file')."

set seed `seed'
if (`target_n' <= `pool_n') {
    sample `target_n', count                          // distinct patients, no replacement
    di as text "Drew " as result `target_n' as text " distinct patients (seed `seed')."
}
else {
    // target exceeds the pool: stochastic replication at the fixed case-mix
    di as text "Target `target_n' > pool `pool_n': expanding (stochastic replication)."
    local base  = floor(`target_n' / `pool_n')
    local extra = `target_n' - `base' * `pool_n'
    gen long _exp = `base'
    if (`extra' > 0) {
        gen double _r = runiform()
        sort _r
        replace _exp = _exp + 1 in 1/`extra'
        drop _r
    }
    expand _exp
    drop _exp
}

replace ID = _n                                       // canonical row order (required by mata_setup / CRN)
save "`cohort_file'", replace

di as text _n "Saved CE cohort: " as result _N as text " patients -> `cohort_file'"
summarize TSD_L${line}S, detail
