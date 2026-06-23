**********
* EpiMAP Myeloma - CE Monte Carlo Precision (multi-scenario)
*
* Purpose: For EACH scenario, from ONE deterministic two-arm run on a shared
*          cohort, produce the DVd-vs-Vd per-patient SD of the increment (sigma_pp)
*          and the TSD 15 convergence figures. The Monte Carlo SD of a size-N mean
*          is MCSD(N) = sigma_pp / sqrt(N); sigma_pp is a population property, so
*          one run per scenario fixes it for every N.
*
*          Writes results/ce_precision_sigma.csv (scenario, sd_pp_dQ, sd_pp_dC, ...)
*          which ce_sample_size.do reads to size N against each scenario's
*          bootstrap parameter SD. Run once, after the cohort pool exists.
*
* Prerequisite: cohort_pool.do has been run (cohort_pool_<line>.dta), and each
*   scenario's deterministic BCR artefacts exist (the .mmat files in outcomes/<scenario>/).
*
* NB: heavy (scenarios x 2 lifetime sims on a large cohort) and UNTESTED against
*     the live engine here - check the first run's log.
*
* Author: Adam Irving
* Date: June 2026
**********

clear all
set more off
macro drop _all

cap cd "/Users/adami/Documents/Monash/Vault/research/models/myeloma model/repo"   // local (Mac)
cap cd "~/em76/adam"                                                              // HPC repo root

**********
* Configuration
**********
local scenarios    "A_trial B_transport C_mrdr"   // scenarios to evaluate
local run_n        = 150000          // cohort size for the convergence run (0 = full pool; >pool = expand)
local seed         = 20260612

local grid         "10000 20000 50000 100000 150000"   // analytic MCSD reference points (printed)
local target_dQ    = 0.005           // absolute MC-SD target for the printed reqN
local batch_sizes  "2000 5000 10000 20000 50000 75000"   // empirical batch sizes (bs<=run_n/2 used)

local chosen_n     = 50000           // cohort size marked on the figures
local curve_min    = 1000            // figure x-range (analytic curve)
local curve_max    = 200000

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

local pool_file   "$patients_path/cohort_pool_${line}.dta"            // 2020-2025 windowed pool (cohort_pool all_l2=0)
local prec_file   "$patients_path/cohort_precision_${line}.dta"       // separate convergence cohort (NOT the production file)
global cohort_file "`prec_file'"                                      // load_patients reads this instead of the canonical patient file
local R           "analyses/$analysis/results"
cap mkdir "`R'"

**********
* Require the pool
**********
capture confirm file "`pool_file'"
if (_rc) {
	di as error "Pool not found: `pool_file'. Run cohort_pool.do first."
	exit 601
}

**********
* Load core programs (once)
**********
run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"
run "core/run_pipeline.do"

**********
* Draw the convergence cohort ONCE (shared by all scenarios) -> cohort_precision_<line>.dta
*   This is a SEPARATE file from the production cohort (load_patients reads it via
*   $cohort_file), so ce_precision never overwrites patients_<analysis>_<line>.dta -
*   the production cohort and this convergence cohort stay independent and reproducible.
**********
use "`pool_file'", clear
local pool_ceiling = _N
if (`run_n' > 0 & `run_n' != `pool_ceiling') {
	set seed `seed'
	if (`run_n' < `pool_ceiling') {
		sample `run_n', count
	}
	else {
		local base  = floor(`run_n' / `pool_ceiling')
		local extra = `run_n' - `base' * `pool_ceiling'
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
}
replace ID = _n
save "`prec_file'", replace

global data     "predicted"
global min_year "1995"
global max_year "2040"
global min_id   "1"
global max_id   "100000000"

* Guard: the scenario list must be populated. If it is empty (e.g. the file was run
*   as a code selection rather than the whole do-file), the loop would silently do
*   nothing - this turns that into a clear error.
if ("`scenarios'" == "") {
	di as error "Config 'scenarios' is empty - run the WHOLE do-file (do ce_precision.do), not a selection."
	exit 198
}

**********
* sigma_pp table (one row per scenario) -> read by ce_sample_size.do
**********
tempname SIG
postfile `SIG' str16 scenario double dQ dC icer sd_pp_dQ sd_pp_dC run_N reqN_dQ ///
	using "`R'/ce_precision_sigma.dta", replace

**********
* Loop scenarios
**********
foreach scen of local scenarios {

	global scenario "`scen'"
	di as text _n "{hline 60}" _n "=== Scenario: `scen' ===" _n "{hline 60}"

	// Intervention arm (DVd)
	global int "dvd"
	qui mata: mata clear
	qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
	run_pipeline
	keep ID cost_total_d qaly_total_d
	rename cost_total_d c_dvd
	rename qaly_total_d q_dvd
	tempfile dvd
	qui save `dvd'

	// Comparator arm (Vd)
	global int "vd"
	qui mata: mata clear
	qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
	run_pipeline
	keep ID cost_total_d qaly_total_d
	rename cost_total_d c_vd
	rename qaly_total_d q_vd

	// Pair by patient ID (common-random-number increment)
	merge 1:1 ID using `dvd', keep(match) nogenerate
	gen double dC_i = c_dvd - c_vd
	gen double dQ_i = q_dvd - q_vd

	// Point estimates and per-patient SD
	qui count
	local N = r(N)
	qui summarize dC_i
	local mC = r(mean)
	local sC = r(sd)
	qui summarize dQ_i
	local mQ = r(mean)
	local sQ = r(sd)
	local icer = `mC' / `mQ'
	qui correlate dC_i dQ_i, covariance
	local cov_CQ = r(cov_12)
	local reqN = ceil((`sQ'/`target_dQ')^2)

	di as text "dC=" as result %12.0fc `mC' as text "  dQ=" %8.4f `mQ' ///
		as text "  ICER=" %12.0fc `icer' as text "  sigma_pp(dQ)=" %7.4f `sQ'

	// record sigma_pp for this scenario (before the figures, so it's safe)
	post `SIG' ("`scen'") (`mQ') (`mC') (`icer') (`sQ') (`sC') (`N') (`reqN')

	// Analytic MCSD over the reference grid (printed)
	di as text "Analytic MCSD(N)=sigma_pp/sqrt(N):  N / MCSD(dQ)"
	foreach n of local grid {
		di as text "   " %9.0fc `n' as result "   " %7.5f (`sQ'/sqrt(`n'))
	}

	// Empirical batch-means (shuffle, split into non-overlapping batches)
	set seed `seed'
	gen double _u = runiform()
	sort _u
	gen long _row = _n

	tempname B M
	postfile `B' double N double mcsd_emp_dQ double mcsd_emp_dC double nbatch ///
		using "`R'/ce_precision_batch_`scen'.dta", replace
	postfile `M' double bs double dQ_batch double dC_batch ///
		using "`R'/ce_precision_estimates_`scen'.dta", replace
	foreach bs of local batch_sizes {
		local k = floor(`N' / `bs')
		if (`k' >= 2) {
			preserve
			keep if _row <= `k' * `bs'
			gen long _batch = ceil(_row / `bs')
			collapse (mean) dQ_i dC_i, by(_batch)
			qui summarize dQ_i
			local seq = r(sd)
			qui summarize dC_i
			local sec = r(sd)
			forvalues j = 1/`=_N' {
				post `M' (`bs') (dQ_i[`j']) (dC_i[`j'])
			}
			restore
			post `B' (`bs') (`seq') (`sec') (`k')
		}
	}
	postclose `B'
	postclose `M'

	**********
	* Figure 1: MCSD(dQ) = sigma_pp/sqrt(N), analytic curve + empirical batch means
	**********
	clear
	local npts = 250
	set obs `npts'
	gen double N = round(exp(ln(`curve_min') + (_n - 1)/(`npts' - 1) * (ln(`curve_max') - ln(`curve_min'))))
	gen double mcsd = `sQ' / sqrt(N)
	gen double mcsd_solid = mcsd if N <= `pool_ceiling'
	gen double mcsd_dash  = mcsd if N >= `pool_ceiling'
	local ychosen = `sQ' / sqrt(`chosen_n')
	append using "`R'/ce_precision_batch_`scen'.dta"
	sort N
	twoway ///
	    (line mcsd_solid N, lcolor(navy) lwidth(medthick)) ///
	    (line mcsd_dash  N, lcolor(navy) lpattern(dash)) ///
	    (scatter mcsd_emp_dQ N, msymbol(Oh) mcolor(maroon) msize(medium)) ///
	    (scatteri `ychosen' `chosen_n', msymbol(O) mcolor(black) msize(medium)) ///
	    , xscale(log) ///
	    xlabel(1000 2000 5000 10000 20000 50000 100000, angle(45) labsize(small)) ///
	    ylabel(, angle(0) format(%5.3f) labsize(small)) ///
	    xtitle("Simulated cohort size, N") ytitle("Monte Carlo SD of incremental QALY") ///
	    text(`ychosen' `chosen_n' "  N = `chosen_n'", place(e) size(small)) ///
	    title("`scen': MC SD of incremental QALY vs N", size(medium)) ///
	    legend(order(1 "Analytic (within pool)" 2 "Analytic extrapolation" 3 "Empirical batch means") ///
	           pos(2) ring(0) cols(1) size(small) region(lstyle(none))) ///
	    graphregion(color(white)) plotregion(margin(medium))
	graph export "`R'/ce_precision_mcsd_`scen'.png", replace width(2200)

	**********
	* Figure 2: TSD 15 convergence funnels (cost / QALY / ICER) vs N
	**********
	clear
	local npts = 200
	set obs `npts'
	gen double bs = round(exp(ln(`curve_min') + (_n - 1)/(`npts' - 1) * (ln(`curve_max') - ln(`curve_min'))))
	gen double q_hi = `mQ' + 1.96 * `sQ' / sqrt(bs)
	gen double q_lo = `mQ' - 1.96 * `sQ' / sqrt(bs)
	gen double c_hi = `mC' + 1.96 * `sC' / sqrt(bs)
	gen double c_lo = `mC' - 1.96 * `sC' / sqrt(bs)
	gen double icer_se = (1/abs(`mQ')) * sqrt( max(0, (`sC'^2)/bs + (`icer'^2)*(`sQ'^2)/bs - 2*`icer'*`cov_CQ'/bs) )
	gen double i_hi = `icer' + 1.96 * icer_se
	gen double i_lo = `icer' - 1.96 * icer_se
	append using "`R'/ce_precision_estimates_`scen'.dta"
	gen double icer_batch = dC_batch / dQ_batch
	sort bs

	twoway (rarea c_hi c_lo bs, color(navy%12) lwidth(none)) ///
	       (scatter dC_batch bs, msymbol(oh) mcolor(maroon) msize(vsmall)) ///
	       , yline(`mC', lcolor(navy)) xline(`chosen_n', lpattern(dot) lcolor(gs9)) ///
	       xscale(log) xlabel(2000 5000 20000 50000 100000, angle(45) labsize(vsmall)) ///
	       ylabel(, angle(0) format(%9.0fc) labsize(vsmall)) ///
	       xtitle("Simulated cohort size, N", size(vsmall)) ytitle("Incremental cost", size(small)) ///
	       title("A. Incremental cost", size(small)) legend(off) name(g_cost, replace) graphregion(color(white))
	twoway (rarea q_hi q_lo bs, color(navy%12) lwidth(none)) ///
	       (scatter dQ_batch bs, msymbol(oh) mcolor(maroon) msize(vsmall)) ///
	       , yline(`mQ', lcolor(navy)) xline(`chosen_n', lpattern(dot) lcolor(gs9)) ///
	       xscale(log) xlabel(2000 5000 20000 50000 100000, angle(45) labsize(vsmall)) ///
	       ylabel(, angle(0) format(%5.3f) labsize(vsmall)) ///
	       xtitle("Simulated cohort size, N", size(vsmall)) ytitle("Incremental QALY", size(small)) ///
	       title("B. Incremental QALY", size(small)) legend(off) name(g_qaly, replace) graphregion(color(white))
	twoway (rarea i_hi i_lo bs, color(navy%12) lwidth(none)) ///
	       (scatter icer_batch bs, msymbol(oh) mcolor(maroon) msize(vsmall)) ///
	       , yline(`icer', lcolor(navy)) xline(`chosen_n', lpattern(dot) lcolor(gs9)) ///
	       xscale(log) xlabel(2000 5000 20000 50000 100000, angle(45) labsize(vsmall)) ///
	       ylabel(, angle(0) format(%9.0fc) labsize(vsmall)) ///
	       xtitle("Simulated cohort size, N", size(vsmall)) ytitle("Cost per QALY (ICER)", size(small)) ///
	       title("C. Cost per QALY (ICER)", size(small)) legend(off) name(g_icer, replace) graphregion(color(white))
	graph combine g_cost g_qaly g_icer, rows(1) xsize(11) ysize(4) ///
	    title("Convergence of incremental outcomes with cohort size (`scen')", size(small)) ///
	    subtitle("Full-sample estimate (line), +/- 1.96 MC SE band, batch estimates (points); dotted N = `chosen_n'", size(vsmall)) ///
	    note("After NICE DSU TSD 15 (Davis et al. 2014), Figures 5-7. ICER band is delta-method.", size(vsmall)) ///
	    graphregion(color(white))
	graph export "`R'/ce_precision_convergence_`scen'.png", replace width(3200)
	di as text "  saved figures -> ce_precision_mcsd_`scen'.png, ce_precision_convergence_`scen'.png"
}
postclose `SIG'

**********
* Export the sigma_pp table
**********
use "`R'/ce_precision_sigma.dta", clear
format dQ %8.4f
format sd_pp_dQ %8.4f
format dC icer sd_pp_dC reqN_dQ %12.0fc
export delimited using "`R'/ce_precision_sigma.csv", replace
di as text _n "Saved sigma_pp table -> `R'/ce_precision_sigma.csv"
list scenario dQ dC icer sd_pp_dQ sd_pp_dC reqN_dQ, clean noobs
