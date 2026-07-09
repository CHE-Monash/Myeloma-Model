**********
* Monash Myeloma Model - Reproducibility verification
*
* Purpose: Guard determinism. With a fixed seed and a fixed cohort the engine MUST produce
*          byte-identical output. Two checks, no golden file (which would need updating on every
*          intended model change):
*            (1) the common-random-number matrix mRN is identical across two mata_setup builds
*                with the same seed;
*            (2) a full run_pipeline pass, repeated in the same session with the same seed, yields
*                an identical Stata dataset signature (checksum of the whole output).
*          Catches non-determinism regressions -- unstable sorts in the load, unordered Mata
*          operations, or RNG drift -- the class of bug that silently breaks CRN variance reduction.
* Usage:   Run from the repository root. No MRDR drive required (uses the in-repo synthetic cohort).
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

**********
* Configuration (small cohort for speed; fixed seed)
**********
global analysis     "default"
global int          "all"
global line         "0"
global coeffs       "full"
global data         "synthetic"
global min_year     "1995"
global max_year     "2040"
global min_id       "1"
global max_id       "3000"          // reduced cohort -- enough to exercise the pipeline
global boot         "0"
global min_bs       ""
global max_bs       ""
global cost_year    "2025"
global drate        "0.05"
global report       "0"
global scenario     ""
global crn_seed_base 20260615       // fixed CRN seed (mata_setup's default; pinned here explicitly)
global coefficients_path "analyses/$analysis/coefficients"
global outcomes_path     "analyses/$analysis/outcomes"
global patients_path     "analyses/$analysis/patients"
global simulated_path    "analyses/$analysis/simulated"

**********
* Load engine
**********
run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"
run "core/mata_functions.do"

qui mata: mata matuse "$coefficients_path/coefficients_$coeffs", replace

* One full deterministic pipeline pass, leaving the output dataset in memory
capture program drop rp_sim
program define rp_sim
	qui load_patients
	qui mata_setup
	qui simulation
	qui process_data
end

* PASS/FAIL tally; second arg is a pre-evaluated 0/1
global rp_run  = 0
global rp_pass = 0
global rp_fail = 0
capture program drop rp_assert
program define rp_assert
	args ok label
	global rp_run = $rp_run + 1
	if `ok' {
		global rp_pass = $rp_pass + 1
		di as text "  " %-54s "`label'" "   " as result "PASS"
	}
	else {
		global rp_fail = $rp_fail + 1
		di as text "  " %-54s "`label'" "   " as error "FAIL"
	}
end

di _n "{hline 90}"
di "REPRODUCIBILITY VERIFICATION  (cohort N <= $max_id, seed $crn_seed_base)"
di "{hline 90}"

**********
* Check 1: the CRN matrix is identical across two seeded builds
**********
qui load_patients
qui mata_setup
mata: mRN_ref = mRN
qui load_patients
qui mata_setup
mata: st_numscalar("crn_diff", max(vec(abs(mRN - mRN_ref))))
mata: st_numscalar("crn_dim",  (rows(mRN)==rows(mRN_ref)) & (cols(mRN)==cols(mRN_ref)))
mata: mata drop mRN_ref
local crn_ok = (crn_dim == 1 & crn_diff == 0)
rp_assert `crn_ok' "CRN matrix mRN identical across two seeded builds"

**********
* Check 2: full run_pipeline output is identical when repeated with the same seed
**********
rp_sim
qui datasignature
local sigA = r(datasignature)
rp_sim
qui datasignature
local sigB = r(datasignature)
di as text "    run A: `sigA'"
di as text "    run B: `sigB'"
local sig_ok = (`"`sigA'"' == `"`sigB'"')
rp_assert `sig_ok' "run_pipeline output identical on re-run (same seed)"

**********
* Summary
**********
di _n "{hline 90}"
di "REPRODUCIBILITY SUMMARY:  " $rp_pass "/" $rp_run " passed"
if $rp_fail == 0  di as result "All reproducibility checks passed -- fixed seed -> byte-identical output."
else              di as error  "$rp_fail check(s) FAILED -- the engine is NON-DETERMINISTIC; investigate."
di "{hline 90}"
