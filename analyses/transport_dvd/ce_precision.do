**********
* EpiMAP Myeloma - DVd L2 Cost-Effectiveness Monte Carlo Precision (single run)
*
* Purpose: From ONE run on the cohort pool, produce the DVd-vs-Vd incremental
*          point estimates AND their Monte Carlo precision, without a size sweep
*          or replicates. The Monte Carlo SD of a size-N mean is
*
*              MCSD(N) = sigma_pp / sqrt(N)
*
*          where sigma_pp is the per-patient SD of the incremental outcome (a
*          population property). We estimate sigma_pp from the per-patient
*          increments in one run, then report MCSD for any N analytically and the
*          N implied by a target precision.
*
* Replaces ce_convergence.do: no reps, no grid of simulations. Reps were only a
*   brute-force way to measure the same run-to-run SD; sigma_pp/sqrt(N) gives it
*   from a single run, and avoids the finite-pool overlap bias that deflates
*   rep-based SDs when N approaches the pool size.
*
* Prerequisite: build_cohort_pool.do has been run (cohort_pool_<line>.dta), and
*   the scenario's deterministic BCR artefacts exist
*   (e.g. outcomes/B_transport/transport_dvd.mmat).
*
* Pairing: sigma_pp uses the PAIRED per-patient increment (dQ_i = QALY_DVd,i -
*   QALY_Vd,i), so both arms are simulated on the same cohort and merged by ID.
*   If the arms are common-random-number aligned this is the (small) paired SD;
*   if not, it is a slightly conservative SD - either way it is the correct
*   per-patient SD of the increment and needs only this one run.
*
* Output: analyses/transport_dvd/results/ce_precision.csv  (point estimates, sigma_pp, MCSD curve)
*
* Usage:  do ce_precision.do
*
* NB: UNTESTED against the live engine here - check the first run's log.
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

local scenario_run "B_transport"     // scenario to evaluate (A_trial / B_transport / C_mrdr)
local run_n        = 0                // 0 = full pool (~296k distinct); batch points honest to pool/2 (~148k)
                                      //   set e.g. 280000 if the full-pool run is too heavy (still 2 batches at 140k)
local seed         = 20260612

// Reference sizes for the analytic MCSD curve, and a target precision for dQ
local grid         "10000 20000 40000 80000 140000"
local target_dQ    = 0.005           // desired Monte Carlo SD of the incremental QALY

// Empirical batch-means convergence (SAME single run): split the run into
//   non-overlapping random batches of each size; SD across batch means = MCSD(size).
local batch_sizes  "2000 5000 10000 20000 50000 100000 140000"   // sizes with <2 batches are skipped automatically

// Publication figure (appendix): MCSD(dQ) = sigma_pp/sqrt(N)
local chosen_n     = 50000           // cohort size chosen for the bootstrap (marked on the figure)
local param_sd_dQ  = 0.06            // bootstrap parameter SD of dQ for the reference line
                                     //   UPDATE once the corrected-cohort bootstrap is run
local curve_min    = 1000            // figure x-range
local curve_max    = 140000          // analytic line/band x-range

global analysis  "transport_dvd"
global line      "2"
global coeffs    "dvd_pre"
global cost_year "2020"
global drate     "0.05"
global boot      "0"
global scenario  "`scenario_run'"

global coefficients_path  "analyses/$analysis/coefficients"
global outcomes_path      "analyses/$analysis/outcomes"
global patients_path      "analyses/$analysis/patients"
global simulated_path     "analyses/$analysis/simulated"

local pool_file   "$patients_path/cohort_pool_all_${line}.dta"   // all-L2 pool (build_cohort_pool all_l2=1); windowed = cohort_pool_${line}.dta
local cohort_file "$patients_path/patients_${analysis}_${line}.dta"
cap mkdir "analyses/$analysis/results"

**********
* Require the pool
**********

capture confirm file "`pool_file'"
if (_rc) {
	di as error "Pool not found: `pool_file'. Run build_cohort_pool.do first."
	exit 601
}

**********
* Load core programs
**********

run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"
run "core/run_pipeline.do"

**********
* Draw the run cohort (full pool by default) -> canonical file
**********

use "`pool_file'", clear
local pool_ceiling = _N              // independent-patient ceiling (solid/dashed split on the figure)
if (`run_n' > 0 & `run_n' != `pool_ceiling') {
	set seed `seed'
	if (`run_n' < `pool_ceiling') {
		sample `run_n', count                       // subsample distinct patients
	}
	else {
		// expand beyond the pool: stochastic replication at fixed case-mix
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
save "`cohort_file'", replace

global data     "predicted"
global min_year "1995"
global max_year "2040"
global min_id   "1"
global max_id   "100000000"

**********
* Simulate both arms on the same cohort; keep per-patient outcomes
**********

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

**********
* Point estimates and per-patient SD (sigma_pp)
**********

qui count
local N = r(N)
qui summarize dC_i
local mC = r(mean)
local sC = r(sd)
qui summarize dQ_i
local mQ = r(mean)
local sQ = r(sd)
local icer = `mC' / `mQ'

di as text _n "=== DVd vs Vd (`scenario_run'), single run on N = `N' patients ==="
di as text "Incremental cost   dC = " as result %12.0fc `mC' as text "   (per-patient SD " %9.0fc `sC' ")"
di as text "Incremental QALY   dQ = " as result %8.4f `mQ' as text "       (per-patient SD " %7.3f `sQ' ")"
di as text "ICER                  = " as result %12.0fc `icer'
di as text _n "Monte Carlo SD at the run size (sigma_pp/sqrt(N)):"
di as text "   MCSD(dC) = " as result %9.1f (`sC'/sqrt(`N')) as text "    MCSD(dQ) = " as result %7.5f (`sQ'/sqrt(`N'))

**********
* Analytic MCSD curve over the reference grid + required N for the target
**********

tempname Q
postfile `Q' double size double mcsd_dC double mcsd_dQ ///
	using "analyses/$analysis/results/ce_precision_raw.dta", replace
di as text _n "Analytic MCSD(N) = sigma_pp/sqrt(N):"
di as text "        N        MCSD(dC)     MCSD(dQ)"
foreach n of local grid {
	local mc = `sC'/sqrt(`n')
	local mq = `sQ'/sqrt(`n')
	di as text %9.0fc `n' as result "   " %9.1f `mc' "     " %7.5f `mq'
	post `Q' (`n') (`mc') (`mq')
}
postclose `Q'

// N needed so that MCSD(dQ) <= target_dQ
local reqN = ceil((`sQ'/`target_dQ')^2)
di as text _n "To get MCSD(dQ) <= " as result %6.4f `target_dQ' as text ", need N >= " as result %9.0fc `reqN'

**********
* Empirical batch-means convergence (same single run, no reps, no re-simulation)
*   Shuffle, then for each batch size bs split the run into floor(N/bs) non-
*   overlapping batches; the SD of the per-batch incremental means estimates
*   MCSD(bs) directly. These points should sit on the analytic sigma_pp/sqrt(N)
*   line - that agreement is the cross-check. (If the run was expanded beyond the
*   pool, batches at large bs reuse patients, so they read as stochastic-only.)
**********

set seed `seed'
gen double _u = runiform()
sort _u
gen long _row = _n

tempname B M
postfile `B' double N double mcsd_emp_dQ double mcsd_emp_dC double nbatch ///
	using "analyses/$analysis/results/ce_precision_batch.dta", replace
postfile `M' double bs double dQ_batch double dC_batch ///
	using "analyses/$analysis/results/ce_precision_estimates.dta", replace
di as text _n "Empirical batch-means MCSD:"
di as text "       bs     nbatch     MCSD(dQ)"
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
		forvalues j = 1/`=_N' {                 // keep each batch estimate for the QALY-vs-size plot
			post `M' (`bs') (dQ_i[`j']) (dC_i[`j'])
		}
		restore
		post `B' (`bs') (`seq') (`sec') (`k')
		di as text %9.0fc `bs' "   " %6.0f `k' as result "    " %7.5f `seq'
	}
}
postclose `B'
postclose `M'

**********
* Publication figure: MCSD(dQ) = sigma_pp/sqrt(N)
*   solid within the independent pool, dashed (analytic extrapolation) beyond it,
*   horizontal reference at the bootstrap parameter SD, marker at the chosen N.
**********

clear
local npts = 250
set obs `npts'
gen double N = round(exp(ln(`curve_min') + (_n - 1)/(`npts' - 1) * (ln(`curve_max') - ln(`curve_min'))))
gen double mcsd = `sQ' / sqrt(N)
gen double mcsd_solid = mcsd if N <= `pool_ceiling'
gen double mcsd_dash  = mcsd if N >= `pool_ceiling'    // overlap by 1 point so the lines join

local ychosen = `sQ' / sqrt(`chosen_n')

// overlay the empirical batch-means points on the analytic curve
append using "analyses/$analysis/results/ce_precision_batch.dta"
sort N

twoway ///
    (line mcsd_solid N, lcolor(navy) lwidth(medthick)) ///
    (line mcsd_dash  N, lcolor(navy) lpattern(dash)) ///
    (scatter mcsd_emp_dQ N, msymbol(Oh) mcolor(maroon) msize(medium)) ///
    (scatteri `ychosen' `chosen_n', msymbol(O) mcolor(black) msize(medium)) ///
    , ///
    yline(`param_sd_dQ', lcolor(red) lpattern(shortdash)) ///
    xscale(log) ///
    xlabel(1000 2000 5000 10000 20000 50000 100000, angle(45) labsize(small)) ///
    ylabel(, angle(0) format(%5.3f) labsize(small)) ///
    xtitle("Simulated cohort size, N") ///
    ytitle("Monte Carlo SD of incremental QALY") ///
    text(`param_sd_dQ' `curve_min' "Bootstrap parameter SD", place(ne) size(small) color(red)) ///
    text(`ychosen' `chosen_n' "  N = `chosen_n'", place(e) size(small)) ///
    legend(order(1 "Analytic (within pool)" 2 "Analytic extrapolation" 3 "Empirical batch means") ///
           pos(2) ring(0) cols(1) size(small) region(lstyle(none))) ///
    graphregion(color(white)) plotregion(margin(medium))

graph export "analyses/$analysis/results/ce_precision_mcsd.png", replace width(2200)
cap graph export "analyses/$analysis/results/ce_precision_mcsd.pdf", replace

// curve data, for re-plotting elsewhere if a different look is wanted
keep N mcsd mcsd_emp_dQ
gen double pool_ceiling = `pool_ceiling'
gen double chosen_n     = `chosen_n'
gen double param_sd_dQ  = `param_sd_dQ'
export delimited using "analyses/$analysis/results/ce_precision_curve.csv", replace
di as text "Saved figure -> results/ce_precision_mcsd.png (.pdf) and curve -> results/ce_precision_curve.csv"

**********
* Figure 2: incremental QALY estimate vs batch size (convergence funnel)
*   batch estimates scatter, the full-sample mean, and a +/- 1.96*sigma_pp/sqrt(N)
*   Monte Carlo band that narrows with N.
**********

clear
local npts = 200
set obs `npts'
gen double bs = round(exp(ln(`curve_min') + (_n - 1)/(`npts' - 1) * (ln(`curve_max') - ln(`curve_min'))))
gen double band_hi = `mQ' + 1.96 * `sQ' / sqrt(bs)
gen double band_lo = `mQ' - 1.96 * `sQ' / sqrt(bs)
append using "analyses/$analysis/results/ce_precision_estimates.dta"
sort bs

twoway ///
    (rarea band_hi band_lo bs, color(navy%12) lwidth(none)) ///
    (scatter dQ_batch bs, msymbol(oh) mcolor(maroon) msize(small)) ///
    , ///
    yline(`mQ', lcolor(navy)) ///
    xline(`chosen_n', lpattern(dot) lcolor(gs9)) ///
    xscale(log) ///
    xlabel(2000 5000 10000 20000 50000 100000, angle(45) labsize(small)) ///
    ylabel(, angle(0) format(%5.3f) labsize(small)) ///
    xtitle("Batch size (simulated patients per estimate)") ///
    ytitle("Incremental QALY, DVd vs Vd") ///
    text(`mQ' `curve_min' "full-sample mean", place(ne) size(small) color(navy)) ///
    legend(order(2 "Batch estimates" 1 "Mean +/- 1.96 MC SE") ///
           pos(1) ring(0) cols(1) size(small) region(lstyle(none))) ///
    graphregion(color(white)) plotregion(margin(medium))

graph export "analyses/$analysis/results/ce_precision_estimate.png", replace width(2200)
cap graph export "analyses/$analysis/results/ce_precision_estimate.pdf", replace
di as text "Saved figure -> results/ce_precision_estimate.png (.pdf)"

**********
* Save a tidy summary
**********

clear
set obs 1
gen scenario   = "`scenario_run'"
gen run_N      = `N'
gen dC         = `mC'
gen dQ         = `mQ'
gen icer       = `icer'
gen sd_pp_dC   = `sC'
gen sd_pp_dQ   = `sQ'
gen mcsd_dC_runN = `sC'/sqrt(`N')
gen mcsd_dQ_runN = `sQ'/sqrt(`N')
gen target_dQ  = `target_dQ'
gen reqN_dQ    = `reqN'
export delimited using "analyses/$analysis/results/ce_precision.csv", replace

di as text _n "Saved analyses/$analysis/results/ce_precision.csv"
