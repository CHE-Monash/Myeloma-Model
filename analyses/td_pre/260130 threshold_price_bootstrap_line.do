**********
* EpiMAP Myeloma - Threshold Price Curve (Combined Lines)
*
* Calculates the maximum cost-effective total cost increment
* at various WTP thresholds, with 95% confidence intervals
* for Lines 2, 3, and 4
**********

clear all
set more off

cd "/Users/adami/Documents/Monash/Research/Blood Disorders/EpiMAP-Local/Myeloma/Simulation"

**********
* Configuration
**********

// WTP thresholds to evaluate (A$/QALY)
local wtp_min = 0
local wtp_max = 150000
local wtp_step = 5000

// Lines to analyse
local lines "2 3 4"

**********
* Process Each Line
**********

local first = 1
foreach line of local lines {
    
    use "analyses/td_pre/simulated/bootstrap_td_pre_`line'.dta", clear
    local n_bs = _N
    
    **********
    * Incremental Outcomes (for each iteration)
    **********
    
    gen inc_cost = cTotald1 - cTotald0
    gen inc_qaly = qTotald1 - qTotald0
    gen inc_os = (os1 - os0)/12
    gen inc_cTXd = cTXd1 - cTXd0
    gen inc_cNTd = cNTd1 - cNTd0
    gen cycles = txd_l`line'_1 * 30.4375 / 28
    
    **********
    * Threshold Price (for each iteration)
    **********
    
    forval w = `wtp_min'(`wtp_step')`wtp_max' {
        gen max_td_cycle_`w' = (`w' * inc_qaly - inc_cost) / cycles
    }
    
    **********
    * Summary Statistics
    **********
    
    // Generate arm-specific outcomes
    gen os_td = os1 / 12
    gen os_soc = os0 / 12
    
    // Store outcome summaries for each arm and incremental
    qui {
        // TD arm
        sum os_td
        scalar mean_os_td_`line' = r(mean)
        _pctile os_td, p(2.5 97.5)
        scalar os_td_lo_`line' = r(r1)
        scalar os_td_hi_`line' = r(r2)
        
        sum qTotald1
        scalar mean_qaly_td_`line' = r(mean)
        _pctile qTotald1, p(2.5 97.5)
        scalar qaly_td_lo_`line' = r(r1)
        scalar qaly_td_hi_`line' = r(r2)
        
        sum cTXd1
        scalar mean_cTXd_td_`line' = r(mean)
        _pctile cTXd1, p(2.5 97.5)
        scalar cTXd_td_lo_`line' = r(r1)
        scalar cTXd_td_hi_`line' = r(r2)
        
        sum cNTd1
        scalar mean_cNTd_td_`line' = r(mean)
        _pctile cNTd1, p(2.5 97.5)
        scalar cNTd_td_lo_`line' = r(r1)
        scalar cNTd_td_hi_`line' = r(r2)
        
        sum cTotald1
        scalar mean_cost_td_`line' = r(mean)
        _pctile cTotald1, p(2.5 97.5)
        scalar cost_td_lo_`line' = r(r1)
        scalar cost_td_hi_`line' = r(r2)
        
        // SoC arm
        sum os_soc
        scalar mean_os_soc_`line' = r(mean)
        _pctile os_soc, p(2.5 97.5)
        scalar os_soc_lo_`line' = r(r1)
        scalar os_soc_hi_`line' = r(r2)
        
        sum qTotald0
        scalar mean_qaly_soc_`line' = r(mean)
        _pctile qTotald0, p(2.5 97.5)
        scalar qaly_soc_lo_`line' = r(r1)
        scalar qaly_soc_hi_`line' = r(r2)
        
        sum cTXd0
        scalar mean_cTXd_soc_`line' = r(mean)
        _pctile cTXd0, p(2.5 97.5)
        scalar cTXd_soc_lo_`line' = r(r1)
        scalar cTXd_soc_hi_`line' = r(r2)
        
        sum cNTd0
        scalar mean_cNTd_soc_`line' = r(mean)
        _pctile cNTd0, p(2.5 97.5)
        scalar cNTd_soc_lo_`line' = r(r1)
        scalar cNTd_soc_hi_`line' = r(r2)
        
        sum cTotald0
        scalar mean_cost_soc_`line' = r(mean)
        _pctile cTotald0, p(2.5 97.5)
        scalar cost_soc_lo_`line' = r(r1)
        scalar cost_soc_hi_`line' = r(r2)
        
        // Incremental
        sum inc_os
        scalar mean_inc_os_`line' = r(mean)
        _pctile inc_os, p(2.5 97.5)
        scalar inc_os_lo_`line' = r(r1)
        scalar inc_os_hi_`line' = r(r2)
        
        sum inc_qaly
        scalar mean_inc_qaly_`line' = r(mean)
        _pctile inc_qaly, p(2.5 97.5)
        scalar inc_qaly_lo_`line' = r(r1)
        scalar inc_qaly_hi_`line' = r(r2)
        
        sum inc_cTXd
        scalar mean_inc_cTXd_`line' = r(mean)
        _pctile inc_cTXd, p(2.5 97.5)
        scalar inc_cTXd_lo_`line' = r(r1)
        scalar inc_cTXd_hi_`line' = r(r2)
        
        sum inc_cNTd
        scalar mean_inc_cNTd_`line' = r(mean)
        _pctile inc_cNTd, p(2.5 97.5)
        scalar inc_cNTd_lo_`line' = r(r1)
        scalar inc_cNTd_hi_`line' = r(r2)
        
        sum inc_cost
        scalar mean_inc_cost_`line' = r(mean)
        _pctile inc_cost, p(2.5 97.5)
        scalar inc_cost_lo_`line' = r(r1)
        scalar inc_cost_hi_`line' = r(r2)
        
        sum cycles
        scalar mean_cycles_`line' = r(mean)
    }
    
    // Create matrix to store TPC results for this line
    local n_wtp = floor((`wtp_max' - `wtp_min') / `wtp_step') + 1
    matrix tpc_`line' = J(`n_wtp', 4, .)
    matrix colnames tpc_`line' = wtp mean lo hi
    
    local row = 1
    forval w = `wtp_min'(`wtp_step')`wtp_max' {
        qui sum max_td_cycle_`w'
        local mean = r(mean)
        qui _pctile max_td_cycle_`w', p(2.5 97.5)
        local lo = r(r1)
        local hi = r(r2)
        
        matrix tpc_`line'[`row', 1] = `w'
        matrix tpc_`line'[`row', 2] = `mean'
        matrix tpc_`line'[`row', 3] = `lo'
        matrix tpc_`line'[`row', 4] = `hi'
        
        local row = `row' + 1
    }
    
    // Convert matrix to temporary dataset
    clear
    svmat tpc_`line', names(col)
    gen line = `line'
    
    if `first' {
        tempfile combined
        save `combined', replace
        local first = 0
    }
    else {
        append using `combined'
        save `combined', replace
    }
}

**********
* Process Subsequent LoT Receipt from Bootstrap Data
**********

foreach line of local lines {
    // Load bootstrap data for this line (contains both CE and pathways data)
    use "analyses/td_pre/simulated/bootstrap_td_pre_`line'.dta", clear
    
    // Calculate mean patient counts at subsequent lines (from line+1 to 9)
    local start_lot = `line' + 1
    forval lot = `start_lot'/9 {
        cap confirm variable int1_lot`lot'
        if !_rc {
            qui {
                sum int1_lot`lot'
                scalar lot`lot'_td_`line' = r(mean)
                
                sum int0_lot`lot'
                scalar lot`lot'_soc_`line' = r(mean)
                
                scalar lot`lot'_diff_`line' = scalar(lot`lot'_td_`line') - scalar(lot`lot'_soc_`line')
            }
        }
        else {
            scalar lot`lot'_td_`line' = .
            scalar lot`lot'_soc_`line' = .
            scalar lot`lot'_diff_`line' = .
        }
    }
}

**********
* Reshape for Plotting
**********

use `combined', clear

// Reshape wide by line
rename mean mean_
rename lo lo_
rename hi hi_
reshape wide mean_ lo_ hi_, i(wtp) j(line)

label var wtp "Willingness-to-Pay (A$/QALY)"

**********
* Generate Combined TPC Figure
**********

twoway (rarea lo_2 hi_2 wtp, color(navy%15) lwidth(none)) ///
       (rarea lo_3 hi_3 wtp, color(maroon%15) lwidth(none)) ///
       (rarea lo_4 hi_4 wtp, color(forest_green%15) lwidth(none)) ///
       (line mean_2 wtp, lcolor(navy) lwidth(medthick)) ///
       (line mean_3 wtp, lcolor(maroon) lwidth(medthick)) ///
       (line mean_4 wtp, lcolor(forest_green) lwidth(medthick)), ///
    title("Tec-Dara vs SoC by Line of Therapy", margin(small)) ///
    xtitle("Willingness-to-Pay Threshold (per QALY)", margin(medium)) ///
    ytitle("Tec-Dara Price (per cycle)", margin(medium)) ///
    xscale(range(0 151000)) ///
    xlabel(0 "A$0" 50000 "A$50K" 100000 "A$100K" 150000 "A$150K") ///
    yscale(range(-5000 5000)) ///
    ylabel(-5000 "-A$5K" 0 "A$0" 5000 "A$5K" 10000 "A$10K") ///
    yline(0, lcolor(gs8) lpattern(solid)) ///
    legend(order(4 "Line 2" 5 "Line 3" 6 "Line 4") rows(1)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "analyses/td_pre/tpc_td_combined.png", replace width(1200)

**********
* Display Results
**********

// Display outcomes summary for each line
foreach line of local lines {
    di _n "{hline 90}"
    di "OUTCOMES (TD vs SoC) - LINE `line'"
    di "{hline 90}"
    di _col(25) "Tec-Dara" _col(50) "SoC" _col(70) "Incremental"
    di "{hline 90}"
    
    di "OS (years)" ///
        _col(20) %6.2f scalar(mean_os_td_`line') " (" %5.2f scalar(os_td_lo_`line') ", " %5.2f scalar(os_td_hi_`line') ")" ///
        _col(45) %6.2f scalar(mean_os_soc_`line') " (" %5.2f scalar(os_soc_lo_`line') ", " %5.2f scalar(os_soc_hi_`line') ")" ///
        _col(70) %6.2f scalar(mean_inc_os_`line') " (" %5.2f scalar(inc_os_lo_`line') ", " %5.2f scalar(inc_os_hi_`line') ")"
    
    di "QALYs" ///
        _col(20) %6.3f scalar(mean_qaly_td_`line') " (" %5.3f scalar(qaly_td_lo_`line') ", " %5.3f scalar(qaly_td_hi_`line') ")" ///
        _col(45) %6.3f scalar(mean_qaly_soc_`line') " (" %5.3f scalar(qaly_soc_lo_`line') ", " %5.3f scalar(qaly_soc_hi_`line') ")" ///
        _col(70) %6.3f scalar(mean_inc_qaly_`line') " (" %5.3f scalar(inc_qaly_lo_`line') ", " %5.3f scalar(inc_qaly_hi_`line') ")"
    
    di "Treatment Cost" ///
        _col(20) %9.0fc scalar(mean_cTXd_td_`line') " (" %9.0fc scalar(cTXd_td_lo_`line') ", " %9.0fc scalar(cTXd_td_hi_`line') ")" ///
        _col(45) %9.0fc scalar(mean_cTXd_soc_`line') " (" %9.0fc scalar(cTXd_soc_lo_`line') ", " %9.0fc scalar(cTXd_soc_hi_`line') ")" ///
        _col(70) %9.0fc scalar(mean_inc_cTXd_`line') " (" %9.0fc scalar(inc_cTXd_lo_`line') ", " %9.0fc scalar(inc_cTXd_hi_`line') ")"
    
    di "Non-TX Cost" ///
        _col(20) %9.0fc scalar(mean_cNTd_td_`line') " (" %9.0fc scalar(cNTd_td_lo_`line') ", " %9.0fc scalar(cNTd_td_hi_`line') ")" ///
        _col(45) %9.0fc scalar(mean_cNTd_soc_`line') " (" %9.0fc scalar(cNTd_soc_lo_`line') ", " %9.0fc scalar(cNTd_soc_hi_`line') ")" ///
        _col(70) %9.0fc scalar(mean_inc_cNTd_`line') " (" %9.0fc scalar(inc_cNTd_lo_`line') ", " %9.0fc scalar(inc_cNTd_hi_`line') ")"
    
    di "Total Cost" ///
        _col(20) %9.0fc scalar(mean_cost_td_`line') " (" %9.0fc scalar(cost_td_lo_`line') ", " %9.0fc scalar(cost_td_hi_`line') ")" ///
        _col(45) %9.0fc scalar(mean_cost_soc_`line') " (" %9.0fc scalar(cost_soc_lo_`line') ", " %9.0fc scalar(cost_soc_hi_`line') ")" ///
        _col(70) %9.0fc scalar(mean_inc_cost_`line') " (" %9.0fc scalar(inc_cost_lo_`line') ", " %9.0fc scalar(inc_cost_hi_`line') ")"
    
    di "{hline 90}"
    di "Mean TD Cycles: " %6.1f scalar(mean_cycles_`line')
    di "{hline 90}"
}

// Display TPC at key thresholds for each line
di _n "{hline 90}"
di "THRESHOLD PRICE PER CYCLE (95% CI) BY LINE"
di "{hline 90}"
di "WTP Threshold" _col(20) "Line 2" _col(45) "Line 3" _col(70) "Line 4"
di "{hline 90}"

foreach w in 25000 50000 75000 100000 150000 {
    local wtp_fmt = string(`w'/1000, "%3.0f") + "K"
    
    qui sum mean_2 if wtp == `w'
    local m2 = r(mean)
    qui sum mean_3 if wtp == `w'
    local m3 = r(mean)
    qui sum mean_4 if wtp == `w'
    local m4 = r(mean)
    
    di "$`wtp_fmt'/QALY" _col(20) "$" %9.0fc `m2' _col(45) "$" %9.0fc `m3' _col(70) "$" %9.0fc `m4'
}
di "{hline 90}"

**********
* Display Subsequent LoT Receipt Table
**********

di _n "{hline 70}"
di "SUBSEQUENT LINE OF THERAPY RECEIPT (mean patient counts)"
di "{hline 70}"

foreach line of local lines {
    local start_lot = `line' + 1
    if `start_lot' <= 9 {
        di _n "Evaluation at Line `line' (Subsequent Lines `start_lot'-9)"
        di "{hline 70}"
        di "LoT" _col(20) "Tec-Dara" _col(35) "SoC" _col(50) "Difference"
        di "{hline 70}"
        
        forval lot = `start_lot'/9 {
            cap scalar dir lot`lot'_td_`line'
            if !_rc & !missing(scalar(lot`lot'_td_`line')) {
                local diff = scalar(lot`lot'_diff_`line')
                local sign = cond(`diff' >= 0, "+", "")
                
                di "Line `lot'" ///
                    _col(20) %8.1f scalar(lot`lot'_td_`line') ///
                    _col(35) %8.1f scalar(lot`lot'_soc_`line') ///
                    _col(50) "`sign'" %7.1f `diff'
            }
        }
        di "{hline 70}"
    }
}

