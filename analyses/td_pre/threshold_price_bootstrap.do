**********
* EpiMAP Myeloma - Threshold Price Curve
*
* Calculates the maximum cost-effective total cost increment
* at various WTP thresholds, with 95% confidence intervals
**********

**********
* Configuration
**********

global line "2"

use "analyses/td_pre/simulated/bootstrap_td_pre_${line}.dta", clear

// WTP thresholds to evaluate (A$/QALY)
local wtp_min = 0
local wtp_max = 150000
local wtp_step = 5000
local n_bs = _N

**********
* Incremental Outcomes (for each iteration)
**********

gen inc_cost = cTotald1 - cTotald0
gen inc_qaly = qTotald1 - qTotald0
gen inc_os = (os1 - os0)/12
gen inc_cTXd = cTXd1 - cTXd0
gen inc_cNTd = cNTd1 - cNTd0
gen cycles = txd_l${line}_1 * 30.4375 / 28

**********
* Threshold Price (for each iteration)
**********

local wtp_values ""
forval w = `wtp_min'(`wtp_step')`wtp_max' {
    local wtp_values "`wtp_values' `w'"
    gen max_td_cycle_`w' = (`w' * inc_qaly - inc_cost) / cycles
}

**********
* Summary Statistics
**********

// Store incremental outcome summaries
qui {
    sum inc_os
    scalar mean_inc_os = r(mean)
    _pctile inc_os, p(2.5 97.5)
    scalar inc_os_lo = r(r1)
    scalar inc_os_hi = r(r2)
    
    sum inc_qaly
    scalar mean_inc_qaly = r(mean)
    _pctile inc_qaly, p(2.5 97.5)
    scalar inc_qaly_lo = r(r1)
    scalar inc_qaly_hi = r(r2)
    
    sum inc_cTXd
    scalar mean_inc_cTXd = r(mean)
    _pctile inc_cTXd, p(2.5 97.5)
    scalar inc_cTXd_lo = r(r1)
    scalar inc_cTXd_hi = r(r2)
    
    sum inc_cNTd
    scalar mean_inc_cNTd = r(mean)
    _pctile inc_cNTd, p(2.5 97.5)
    scalar inc_cNTd_lo = r(r1)
    scalar inc_cNTd_hi = r(r2)
    
    sum inc_cost
    scalar mean_inc_cost = r(mean)
    _pctile inc_cost, p(2.5 97.5)
    scalar inc_cost_lo = r(r1)
    scalar inc_cost_hi = r(r2)
    
    sum cycles
    scalar mean_cycles = r(mean)
}

// Create matrix to store TPC results
local n_wtp = floor((`wtp_max' - `wtp_min') / `wtp_step') + 1
matrix tpc_results = J(`n_wtp', 4, .)
matrix colnames tpc_results = wtp mean lo hi

local row = 1
forval w = `wtp_min'(`wtp_step')`wtp_max' {
    qui sum max_td_cycle_`w'
    local mean = r(mean)
    qui _pctile max_td_cycle_`w', p(2.5 97.5)
    local lo = r(r1)
    local hi = r(r2)
    
    matrix tpc_results[`row', 1] = `w'
    matrix tpc_results[`row', 2] = `mean'
    matrix tpc_results[`row', 3] = `lo'
    matrix tpc_results[`row', 4] = `hi'
    
    local row = `row' + 1
}

**********
* Create TPC Dataset for Plotting
**********

// Convert matrix to dataset
clear
svmat tpc_results, names(col)

label var wtp "Willingness-to-Pay (A$/QALY)"
label var mean "Threshold Price per Cycle (Mean)"
label var lo "Lower 95% CI"
label var hi "Upper 95% CI"

**********
* Generate TPC Figure
**********

twoway (rarea lo hi wtp, color(navy%20) lwidth(none)) ///
       (line mean wtp, lcolor(navy) lwidth(medthick)), ///
    title("Threshold Price Curve") ///
	title("Tec-Dara vs SoC at Line ${line}", margin(small)) ///
    xtitle("Willingness-to-Pay Threshold (per QALY)", margin(medium)) ///
    ytitle("Tec-Dara Price (per cycle)", margin(medium)) ///
   	xscale(range(0 151000)) ///
	xlabel(0 "A$0" 50000 "A$50K" 100000 "A$100K" 150000 "A$150K") ///
	yscale(range(-5000 5000)) ///
	ylabel(-5000 "-A$5K" 0 "A$0" 5000 "A$5K") ///
	yline(0, lcolor(maroon) lpattern(solid)) ///
    legend(order(2 "Mean" 1 "95% CI") rows(1)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "analyses/td_pre/tpc_td_l${line}.png", replace width(1200)

**********
* Display and Export Results
**********

// Display incremental outcomes summary
di _n "{hline 70}"
di "INCREMENTAL OUTCOMES (TD vs SoC) - `n_bs' bootstrap samples"
di "{hline 70}"
di "Inc OS (years):     " %6.2f scalar(mean_inc_os) " (" %6.2f scalar(inc_os_lo) ", " %6.2f scalar(inc_os_hi) ")"
di "Inc QALYs:          " %6.3f scalar(mean_inc_qaly) " (" %6.3f scalar(inc_qaly_lo) ", " %6.3f scalar(inc_qaly_hi) ")"
di "Inc Treatment Cost: " %12.0fc scalar(mean_inc_cTXd) " (" %12.0fc scalar(inc_cTXd_lo) ", " %12.0fc scalar(inc_cTXd_hi) ")"
di "Inc Non-TX Cost:    " %12.0fc scalar(mean_inc_cNTd) " (" %12.0fc scalar(inc_cNTd_lo) ", " %12.0fc scalar(inc_cNTd_hi) ")"
di "Inc Total Cost:     " %12.0fc scalar(mean_inc_cost) " (" %12.0fc scalar(inc_cost_lo) ", " %12.0fc scalar(inc_cost_hi) ")"
di "Mean TD Cycles:     " %6.1f scalar(mean_cycles)
di "{hline 70}"

// Display TPC at key thresholds
di _n "{hline 70}"
di "THRESHOLD PRICE PER CYCLE (95% CI)"
di "{hline 70}"
di "WTP Threshold" _col(25) "Mean" _col(40) "95% CI"
di "{hline 70}"

qui foreach w in 25000 50000 75000 100000 150000 {
    qui sum mean if wtp == `w'
    local m = r(mean)
    qui sum lo if wtp == `w'
    local l = r(mean)
    qui sum hi if wtp == `w'
    local h = r(mean)
    
    local wtp_fmt = string(`w'/1000, "%3.0f") + "K"
    noi di "$`wtp_fmt'/QALY" _col(20) "$" %9.0fc `m' _col(35) "($" %9.0fc `l' ", $" %9.0fc `h' ")"
}
di "{hline 70}"
