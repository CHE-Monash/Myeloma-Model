**********
* EpiMAP Myeloma - DVd L2 Cost-Effectiveness Cohort-Size Convergence Sweep
*
* Purpose: For a grid of simulated sample sizes, draw replicate cohorts from the
*          existing pool, simulate both arms (DVd, Vd) on each, and report the
*          incremental cost/QALY and their Monte Carlo SD ACROSS REPLICATES - the
*          empirical "MC error vs simulated sample size" curve.
*
* Sizing (sample vs expand):
*   - size <= pool : drawn as DISTINCT patients (sample, no replacement).
*   - size >  pool : built with `expand` - each pool patient is simulated several
*       times. The engine never re-seeds, so a fresh per-rep seed flows into the
*       simulation and each copy gets INDEPENDENT stochastic draws. This is
*       legitimate Monte Carlo replication of the STOCHASTIC component at a FIXED
*       case-mix. The across-rep SD therefore measures stochastic MC error - the
*       dominant component here, given the large per-patient increment SD - NOT
*       case-mix-sampling error, so do not read the >pool points as "N
*       independent patients".
*
* Prerequisite: build_cohort_pool.do has been run (creates cohort_pool_<line>.dta).
*               The scenario's deterministic BCR artefacts must exist
*               (e.g. outcomes/B_transport/transport_dvd.mmat).
*
* Reading it: judge on incremental COST and QALY, not the ICER (the ratio is
*   unstable when incremental QALY is near zero). dC_mcsd / dQ_mcsd are the Monte
*   Carlo SD across replicates at each size; pick the smallest size where they
*   are small relative to your bootstrap CI width.
*
* Output: results/ce_convergence_raw.dta   (per replicate)
*         results/ce_convergence.csv        (table)
*         results/ce_convergence_mcsd.png/.pdf  (empirical MC SD vs simulated size)
*
* Usage:  do ce_convergence.do
*
* NB: heavy (2 lifetime simulations x reps x sizes x scenarios) and UNTESTED
*     against the live engine here - check the first run's log.
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

local sizes         "12500 25000 50000 100000 200000"  // <=pool = distinct; >pool = expand
local reps          = 10                         // replicates per size; SD precision ~ 1/sqrt(2(reps-1)), lower if compute-bound
local scenario_list "B_transport"                // A_trial / B_transport / C_mrdr
local seedbase      = 20260612

global analysis  "transport_dvd"
global line      "2"
global coeffs    "dvd_pre"
global cost_year "2020"
global drate     "0.05"
global boot      "0"

global coefficients_path  "analyses/$analysis/coefficients"
global outcomes_path      "analyses/$analysis/outcomes"
global patients_path      "analyses/$analysis/patients"
global simulated_path     "analyses/$analysis/simulated"

local pool_file   "$patients_path/cohort_pool_all_${line}.dta"   // all-L2 pool (build_cohort_pool all_l2=1); windowed = cohort_pool_${line}.dta
local cohort_file "$patients_path/patients_${analysis}_${line}.dta"   // canonical (read by load_patients)
cap mkdir "analyses/$analysis/results"

**********
* Require the pool (built once by build_cohort_pool.do)
**********

capture confirm file "`pool_file'"
if (_rc) {
	di as error "Pool not found: `pool_file'."
	di as error "Run build_cohort_pool.do first (it is the expensive, build-once step)."
	exit 601
}
qui use "`pool_file'", clear
local pool_n = _N
di as text "Using pool: " as result `pool_n' as text " patients (`pool_file')."

**********
* Load core programs (incl. the shared engine pass)
**********

run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"
run "core/run_pipeline.do"

// Outcome-simulation settings (cohort read from the canonical predicted file)
global data     "predicted"
global min_year "1995"
global max_year "2040"
global min_id   "1"
global max_id   "100000000"      // do not truncate the drawn cohort

**********
* Sweep sizes x replicates
**********

tempname P
postfile `P' str16 scenario size rep double dC double dQ double icer ///
	using "analyses/$analysis/results/ce_convergence_raw.dta", replace

local sd = `seedbase'        // running seed: a distinct value per (size, rep, arm-set)

foreach scen of local scenario_list {
	global scenario "`scen'"
	di as text _n "=== Scenario: `scen' ==="

	foreach n of local sizes {
		forvalues r = 1/`reps' {

			// draw a size-n cohort from the pool; fresh seed feeds the simulation too
			local ++sd
			use "`pool_file'", clear
			set seed `sd'
			if (`n' <= `pool_n') {
				sample `n', count                       // distinct patients (no replacement)
			}
			else {
				// expand: each patient `base' times, plus one extra copy for `extra' random patients,
				// so the cohort is exactly `n' rows of the SAME case-mix (stochastic replication)
				local base  = floor(`n' / `pool_n')
				local extra = `n' - `base' * `pool_n'
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
			replace ID = _n
			save "`cohort_file'", replace

			// intervention arm (DVd)
			global int "dvd"
			qui mata: mata clear
			qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
			run_pipeline
			qui summarize cost_total_d, meanonly
			local cost1 = r(mean)
			qui summarize qaly_total_d, meanonly
			local qaly1 = r(mean)

			// comparator arm (Vd)
			global int "vd"
			qui mata: mata clear
			qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
			run_pipeline
			qui summarize cost_total_d, meanonly
			local cost0 = r(mean)
			qui summarize qaly_total_d, meanonly
			local qaly0 = r(mean)

			local dC = `cost1' - `cost0'
			local dQ = `qaly1' - `qaly0'
			local icer = cond(`dQ' == 0, ., `dC' / `dQ')
			post `P' ("`scen'") (`n') (`r') (`dC') (`dQ') (`icer')

			di as text "  `scen'  N=`n'  rep `r'/`reps'" ///
				as result "   dC=" %9.0fc `dC' "   dQ=" %6.3f `dQ' "   ICER=" %9.0fc `icer'
		}
	}
}
postclose `P'

**********
* Convergence table
**********

use "analyses/$analysis/results/ce_convergence_raw.dta", clear
collapse (mean) dC_mean=dC dQ_mean=dQ ///
         (sd)   dC_mcsd=dC dQ_mcsd=dQ ///
         (count) reps=dC, by(scenario size)

gen icer_mean = dC_mean / dQ_mean
gen cv_dQ = dQ_mcsd / dQ_mean
gen cv_dC = dC_mcsd / dC_mean

format dC_mean dC_mcsd icer_mean %12.0fc
format dQ_mean dQ_mcsd %8.4f
format cv_dQ cv_dC %6.3f
order scenario size reps dC_mean dC_mcsd cv_dC dQ_mean dQ_mcsd cv_dQ icer_mean
sort scenario size

di as text _n "=== DVd vs Vd: cohort-size convergence ==="
list, clean noobs
export delimited using "analyses/$analysis/results/ce_convergence.csv", replace

**********
* Figure: empirical Monte Carlo SD of incremental QALY vs simulated size
*   dashed line marks the pool size (distinct patients to the left, expand to the right)
**********

qui summarize dQ_mcsd
local ytop = r(max)

twoway (connected dQ_mcsd size, sort msymbol(O) lcolor(navy) mcolor(navy)) ///
	, ///
	xline(`pool_n', lpattern(dash) lcolor(gs8)) ///
	xscale(log) ///
	xlabel(12500 25000 50000 100000 200000, angle(45) labsize(small)) ///
	ylabel(, angle(0) format(%6.4f) labsize(small)) ///
	xtitle("Simulated sample size") ///
	ytitle("Monte Carlo SD of incremental QALY (across reps)") ///
	text(`ytop' `pool_n' "pool: distinct -> expand", place(e) size(small) color(gs8)) ///
	legend(off) graphregion(color(white)) plotregion(margin(medium))

graph export "analyses/$analysis/results/ce_convergence_mcsd.png", replace width(2200)
cap graph export "analyses/$analysis/results/ce_convergence_mcsd.pdf", replace
di as text "Saved figure -> results/ce_convergence_mcsd.png (.pdf)"

di as text _n "NOTE: the canonical cohort file now holds the last drawn sample." ///
	" Draw a production cohort (or rebuild from the pool) before a production run."
