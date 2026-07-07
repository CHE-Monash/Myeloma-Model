**********
* Monash Myeloma Model - Draw Cohort (template)
*
* Purpose: second step of the line-specific cohort workflow (see cohort_pool.do). Draw a fixed-size,
*          fixed-seed sample of the line-entry pool and save it as the canonical patient file that the
*          dispatcher's load_patients reads.
* Usage:   do draw_cohort.do -- edit the Configuration block first. Needs cohort_pool.do's output
*          patients/cohort_pool_<line>.dta; writes patients/patients_<analysis>_<line>.dta.
* Notes:   Cohort SIZE is a precision choice, independent of case-mix: pick N large enough that Monte
*          Carlo error is small vs parameter uncertainty (transport_dvd used N=50,000 via a 5% rule; see
*          its results/ce_sample_size_appendix.md). Fixed seed -> reproducible draw. Worked example:
*          analyses/transport_dvd/patients/sample_size/ce_cohort.do.
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)
macro drop _all

**********
* Configuration  -- EDIT THESE
**********
global analysis "template"        // <- your analysis
global line     "2"               // <- the decision line L (must match cohort_pool.do)
local  target_n = 50000           // <- cohort size (precision choice; see ce_sample_size for a rule)
local  seed     = 20260616        // <- fixed -> reproducible cohort

global patients_path "analyses/$analysis/patients"
local pool_file   "$patients_path/cohort_pool_${line}.dta"
local cohort_file "$patients_path/patients_${analysis}_${line}.dta"   // canonical file (load_patients reads this)

**********
* Require the pool (and never read and overwrite the same file)
**********
capture confirm file "`pool_file'"
if (_rc) {
    di as error "Pool not found: `pool_file'"
    di as error "Run cohort_pool.do first."
    exit 601
}
if ("`pool_file'" == "`cohort_file'") {
    di as error "pool_file and cohort_file are the same - point pool_file at the pool, not the output."
    exit 198
}

**********
* Draw the cohort
**********
use "`pool_file'", clear
local pool_n = _N
di as text "Pool: " as result `pool_n' as text " line-${line} entrants (`pool_file')."

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

di as text _n "Saved cohort: " as result _N as text " patients -> `cohort_file'"
summarize TSD_L${line}S, detail
