**********
* EpiMAP Myeloma - Threshold Price Curve (Combined Lines)
*
* Calculates the maximum cost-effective total cost increment
* at various WTP thresholds, with 95% confidence intervals
* for Lines 2, 3, and 4
*
* Naming Convention:
*   Uses pre-calculated variables from bootstrap dataset:
*   cost_{component}[_d]_{arm}, qaly_total[_d]_{arm}, inc_*
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
    * Derived Variables for Threshold Price
    **********
    
    // Treatment cycles (needed for per-cycle threshold price)
    gen cycles = txd_L`line'_1 * 30.4375 / 28
    
    // OS in years (for display)
    gen os_yrs_1 = os_1 / 12
    gen os_yrs_0 = os_0 / 12
    gen inc_os_yrs = inc_os / 12
    
    **********
    * Threshold Price (for each iteration)
    **********
    
    forval w = `wtp_min'(`wtp_step')`wtp_max' {
        gen max_td_cycle_`w' = (`w' * inc_qaly_d - inc_cost_d) / cycles
    }
    
    **********
    * Summary Statistics
    **********
    
    qui {
        // Intervention arm
        sum os_yrs_1
        scalar mean_os_1_`line' = r(mean)
        _pctile os_yrs_1, p(2.5 97.5)
        scalar os_1_lo_`line' = r(r1)
        scalar os_1_hi_`line' = r(r2)
        
        sum qaly_total_d_1
        scalar mean_qaly_1_`line' = r(mean)
        _pctile qaly_total_d_1, p(2.5 97.5)
        scalar qaly_1_lo_`line' = r(r1)
        scalar qaly_1_hi_`line' = r(r2)
        
        sum cost_tx_d_1
        scalar mean_cost_tx_d_1_`line' = r(mean)
        _pctile cost_tx_d_1, p(2.5 97.5)
        scalar cost_tx_d_1_lo_`line' = r(r1)
        scalar cost_tx_d_1_hi_`line' = r(r2)
        
        sum cost_nt_d_1
        scalar mean_cost_nt_d_1_`line' = r(mean)
        _pctile cost_nt_d_1, p(2.5 97.5)
        scalar cost_nt_d_1_lo_`line' = r(r1)
        scalar cost_nt_d_1_hi_`line' = r(r2)
        
        sum cost_total_d_1
        scalar mean_cost_1_`line' = r(mean)
        _pctile cost_total_d_1, p(2.5 97.5)
        scalar cost_1_lo_`line' = r(r1)
        scalar cost_1_hi_`line' = r(r2)
        
        // Control arm
        sum os_yrs_0
        scalar mean_os_0_`line' = r(mean)
        _pctile os_yrs_0, p(2.5 97.5)
        scalar os_0_lo_`line' = r(r1)
        scalar os_0_hi_`line' = r(r2)
        
        sum qaly_total_d_0
        scalar mean_qaly_0_`line' = r(mean)
        _pctile qaly_total_d_0, p(2.5 97.5)
        scalar qaly_0_lo_`line' = r(r1)
        scalar qaly_0_hi_`line' = r(r2)
        
        sum cost_tx_d_0
        scalar mean_cost_tx_d_0_`line' = r(mean)
        _pctile cost_tx_d_0, p(2.5 97.5)
        scalar cost_tx_d_0_lo_`line' = r(r1)
        scalar cost_tx_d_0_hi_`line' = r(r2)
        
        sum cost_nt_d_0
        scalar mean_cost_nt_d_0_`line' = r(mean)
        _pctile cost_nt_d_0, p(2.5 97.5)
        scalar cost_nt_d_0_lo_`line' = r(r1)
        scalar cost_nt_d_0_hi_`line' = r(r2)
        
        sum cost_total_d_0
        scalar mean_cost_0_`line' = r(mean)
        _pctile cost_total_d_0, p(2.5 97.5)
        scalar cost_0_lo_`line' = r(r1)
        scalar cost_0_hi_`line' = r(r2)
        
        // Incremental (pre-calculated in bootstrap dataset)
        sum inc_os_yrs
        scalar mean_inc_os_`line' = r(mean)
        _pctile inc_os_yrs, p(2.5 97.5)
        scalar inc_os_lo_`line' = r(r1)
        scalar inc_os_hi_`line' = r(r2)
        
        sum inc_qaly_d
        scalar mean_inc_qaly_`line' = r(mean)
        _pctile inc_qaly_d, p(2.5 97.5)
        scalar inc_qaly_lo_`line' = r(r1)
        scalar inc_qaly_hi_`line' = r(r2)
        
        sum inc_cost_tx_d
        scalar mean_inc_cost_tx_d_`line' = r(mean)
        _pctile inc_cost_tx_d, p(2.5 97.5)
        scalar inc_cost_tx_d_lo_`line' = r(r1)
        scalar inc_cost_tx_d_hi_`line' = r(r2)
        
        sum inc_cost_nt_d
        scalar mean_inc_cost_nt_d_`line' = r(mean)
        _pctile inc_cost_nt_d, p(2.5 97.5)
        scalar inc_cost_nt_d_lo_`line' = r(r1)
        scalar inc_cost_nt_d_hi_`line' = r(r2)
        
        sum inc_cost_d
        scalar mean_inc_cost_`line' = r(mean)
        _pctile inc_cost_d, p(2.5 97.5)
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
                scalar lot`lot'_1_`line' = r(mean)
                
                sum int0_lot`lot'
                scalar lot`lot'_0_`line' = r(mean)
                
                scalar lot`lot'_diff_`line' = scalar(lot`lot'_1_`line') - scalar(lot`lot'_0_`line')
            }
        }
        else {
            scalar lot`lot'_1_`line' = .
            scalar lot`lot'_0_`line' = .
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
    di "OUTCOMES (Intervention vs Control) - LINE `line'"
    di "{hline 90}"
    di _col(25) "Intervention" _col(50) "Control" _col(70) "Incremental"
    di "{hline 90}"
    
    di "OS (years)" ///
        _col(20) %6.2f scalar(mean_os_1_`line') " (" %5.2f scalar(os_1_lo_`line') ", " %5.2f scalar(os_1_hi_`line') ")" ///
        _col(45) %6.2f scalar(mean_os_0_`line') " (" %5.2f scalar(os_0_lo_`line') ", " %5.2f scalar(os_0_hi_`line') ")" ///
        _col(70) %6.2f scalar(mean_inc_os_`line') " (" %5.2f scalar(inc_os_lo_`line') ", " %5.2f scalar(inc_os_hi_`line') ")"
    
    di "QALYs" ///
        _col(20) %6.3f scalar(mean_qaly_1_`line') " (" %5.3f scalar(qaly_1_lo_`line') ", " %5.3f scalar(qaly_1_hi_`line') ")" ///
        _col(45) %6.3f scalar(mean_qaly_0_`line') " (" %5.3f scalar(qaly_0_lo_`line') ", " %5.3f scalar(qaly_0_hi_`line') ")" ///
        _col(70) %6.3f scalar(mean_inc_qaly_`line') " (" %5.3f scalar(inc_qaly_lo_`line') ", " %5.3f scalar(inc_qaly_hi_`line') ")"
    
    di "Treatment Cost" ///
        _col(20) %9.0fc scalar(mean_cost_tx_d_1_`line') " (" %9.0fc scalar(cost_tx_d_1_lo_`line') ", " %9.0fc scalar(cost_tx_d_1_hi_`line') ")" ///
        _col(45) %9.0fc scalar(mean_cost_tx_d_0_`line') " (" %9.0fc scalar(cost_tx_d_0_lo_`line') ", " %9.0fc scalar(cost_tx_d_0_hi_`line') ")" ///
        _col(70) %9.0fc scalar(mean_inc_cost_tx_d_`line') " (" %9.0fc scalar(inc_cost_tx_d_lo_`line') ", " %9.0fc scalar(inc_cost_tx_d_hi_`line') ")"
    
    di "Non-TX Cost" ///
        _col(20) %9.0fc scalar(mean_cost_nt_d_1_`line') " (" %9.0fc scalar(cost_nt_d_1_lo_`line') ", " %9.0fc scalar(cost_nt_d_1_hi_`line') ")" ///
        _col(45) %9.0fc scalar(mean_cost_nt_d_0_`line') " (" %9.0fc scalar(cost_nt_d_0_lo_`line') ", " %9.0fc scalar(cost_nt_d_0_hi_`line') ")" ///
        _col(70) %9.0fc scalar(mean_inc_cost_nt_d_`line') " (" %9.0fc scalar(inc_cost_nt_d_lo_`line') ", " %9.0fc scalar(inc_cost_nt_d_hi_`line') ")"
    
    di "Total Cost" ///
        _col(20) %9.0fc scalar(mean_cost_1_`line') " (" %9.0fc scalar(cost_1_lo_`line') ", " %9.0fc scalar(cost_1_hi_`line') ")" ///
        _col(45) %9.0fc scalar(mean_cost_0_`line') " (" %9.0fc scalar(cost_0_lo_`line') ", " %9.0fc scalar(cost_0_hi_`line') ")" ///
        _col(70) %9.0fc scalar(mean_inc_cost_`line') " (" %9.0fc scalar(inc_cost_lo_`line') ", " %9.0fc scalar(inc_cost_hi_`line') ")"
    
    di "{hline 90}"
    di "Mean Intervention Cycles: " %6.1f scalar(mean_cycles_`line')
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
        di "LoT" _col(20) "Intervention" _col(35) "Control" _col(50) "Difference"
        di "{hline 70}"
        
        forval lot = `start_lot'/9 {
            cap scalar dir lot`lot'_1_`line'
            if !_rc & !missing(scalar(lot`lot'_1_`line')) {
                local diff = scalar(lot`lot'_diff_`line')
                local sign = cond(`diff' >= 0, "+", "")
                
                di "Line `lot'" ///
                    _col(20) %8.1f scalar(lot`lot'_1_`line') ///
                    _col(35) %8.1f scalar(lot`lot'_0_`line') ///
                    _col(50) "`sign'" %7.1f `diff'
            }
        }
        di "{hline 70}"
    }
}
