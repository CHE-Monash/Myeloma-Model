**********
* Monash Myeloma Model - Process Bootstraps
*
* Purpose: Process bootstrap simulation results 
* Usage: Set configuration, then run: results, maxbs(500)
* Output: Summary dataset
*
* Naming Convention:
*   cost_{component}[_d]_{arm}  - Cost variables (tx, nt, total)
*   qaly_{component}[_d]_{arm}  - QALY variables
*   inc_{outcome}[_d]          - Incremental outcomes
*   _d suffix = discounted
*   _1 = intervention, _0 = control
**********

clear all

**********
// Configuration
**********
global analysis "td_pre"
global int1_name "Tec-Dara"
global int0_name "SoC"
global int1_file "td"
global int0_file "soc"
global line `1'
global input_path "analyses/$analysis/simulated/bootstrap"
global output_path "analyses/$analysis/simulated"
global file_name "predicted_B"
global drate = 0.05

// Output Flags
global OUT_SUMMARY = 1
global OUT_PATHWAYS = 1
global OUT_SURVIVAL = 1

**********
// Main Program
capture program drop results
program define results
    syntax [, maxbs(integer 500)]
    di _n "{hline 80}" _n "Processing Bootstrap Results: $int1_name vs $int0_name (Line $line)" _n "{hline 80}"
    global MaxBS = `maxbs'
    initialise_frames
    process_bootstrap_samples
    save_summary_dataset
    if $OUT_SUMMARY generate_summary_statistics
    if $OUT_PATHWAYS analyse_treatment_pathways
    di _n "{hline 80}" _n "Complete. Summary: $output_path/bootstrap_${analysis}_${line}.dta" _n "{hline 80}"
end

**********
// Initialise Frames
capture program drop initialise_frames
program define initialise_frames
    // Cost and QALY variables by arm (1=intervention, 0=control)
    local v "sample"
    local v "`v' cost_tx_1 cost_nt_1 cost_total_1"
    local v "`v' cost_tx_d_1 cost_nt_d_1 cost_total_d_1"
    local v "`v' cost_tx_0 cost_nt_0 cost_total_0"
    local v "`v' cost_tx_d_0 cost_nt_d_0 cost_total_d_0"
    local v "`v' qaly_total_1 qaly_total_d_1"
    local v "`v' qaly_total_0 qaly_total_d_0"
	local v "`v' os_1 os_0 txd_L${line}_1 txd_L${line}_0 txd_L${line}_med_1 txd_L${line}_med_0"
    
    cap frame drop bs_results
    frame create bs_results `v'
    
    // Build pathways frame variables dynamically based on evaluation line
    local pv "sample int1_n int0_n int1_cr int1_vg int1_pr int1_mr int1_sd int1_pd int0_cr int0_vg int0_pr int0_mr int0_sd int0_pd"
    local start_lot = $line + 1
    forval lot = `start_lot'/9 {
        local pv "`pv' int1_lot`lot' int0_lot`lot'"
    }
    cap frame drop bs_pathways
    frame create bs_pathways `pv'
end

**********
// Process Bootstrap Samples
capture program drop process_bootstrap_samples
program define process_bootstrap_samples
    forval b = 1/$MaxBS {
        if mod(`b',50)==0 | `b'==1 di "  Processing Sample `b' of $MaxBS..." 
        local f0 "$input_path/${int0_file}_${line}_${file_name}`b'.dta"
        cap use "`f0'", clear
        if _rc di "  NOT FOUND"
        gen Int = 0
        extract_arm_statistics `b' 0
        local ctrl "$input_path/temp_ctrl_`b'.dta"
        qui save `ctrl', replace
        local f1 "$input_path/${int1_file}_${line}_${file_name}`b'.dta"
        cap use "`f1'", clear
        if _rc di "  NOT FOUND"
        gen Int = 1
        extract_arm_statistics `b' 1
        qui append using `ctrl'
        erase "`ctrl'"
        extract_combined_metrics `b'
    }
end

**********
// Extract Arm Statistics
capture program drop extract_arm_statistics
program define extract_arm_statistics
    args b arm
    qui {
        local n = _N
        local L = $line
        count if BCR_L`L'==1
        local cr = r(N)/_N*100
        count if BCR_L`L'==2
        local vg = r(N)/_N*100
        count if BCR_L`L'==3
        local pr = r(N)/_N*100
        count if BCR_L`L'==4
        local mr = r(N)/_N*100
        count if BCR_L`L'==5
        local sd = r(N)/_N*100
        count if BCR_L`L'==6
        local pd = r(N)/_N*100
        
        // Count patients receiving subsequent lines (from line+1 to 9)
        local start_lot = $line + 1
        forval lot = `start_lot'/9 {
            cap count if TXD_L`lot'!=. & TXD_L`lot'>0
            local lot`lot' = cond(_rc==0, r(N), 0)
        }
    }
    
    if `arm'==0 {
        // Store control arm values in globals
        foreach x in n cr vg pr mr sd pd {
            global int0_`x' = ``x''
        }
        local start_lot = $line + 1
        forval lot = `start_lot'/9 {
            global int0_lot`lot' = `lot`lot''
        }
    }
    if `arm'==1 {
        // Build frame post command dynamically
        local post_vals "(`b') (`n') ($int0_n) (`cr') (`vg') (`pr') (`mr') (`sd') (`pd') ($int0_cr) ($int0_vg) ($int0_pr) ($int0_mr) ($int0_sd) ($int0_pd)"
        local start_lot = $line + 1
        forval lot = `start_lot'/9 {
            local post_vals "`post_vals' (`lot`lot'') (${int0_lot`lot'})"
        }
        frame post bs_pathways `post_vals'
    }
end

**********
// Extract Combined Metrics
capture program drop extract_combined_metrics
program define extract_combined_metrics
    args b
    qui {
        foreach i in 1 0 {
            // Costs (undiscounted)
            sum cost_tx if Int==`i'
            local cost_tx_`i' = r(mean)
            sum cost_nt if Int==`i'
            local cost_nt_`i' = r(mean)
            sum cost_total if Int==`i'
            local cost_total_`i' = r(mean)
            
            // Costs (discounted)
            sum cost_tx_d if Int==`i'
            local cost_tx_d_`i' = r(mean)
            sum cost_nt_d if Int==`i'
            local cost_nt_d_`i' = r(mean)
            sum cost_total_d if Int==`i'
            local cost_total_d_`i' = r(mean)
            
            // QALYs
            sum qaly_total if Int==`i'
            local qaly_total_`i' = r(mean)
            sum qaly_total_d if Int==`i'
            local qaly_total_d_`i' = r(mean)
            
            // Survival and treatment duration
            sum OC_TIME_L if Int==`i'
            local os_`i' = r(mean)
            sum TXD_L${line} if Int==`i', d
            local txd_L${line}_`i' = r(mean)
            local txd_L${line}_med_`i' = r(p50)
        }
    }
    
    local p "(`b')"
    // Intervention arm
    local p "`p' (`cost_tx_1') (`cost_nt_1') (`cost_total_1')"
    local p "`p' (`cost_tx_d_1') (`cost_nt_d_1') (`cost_total_d_1')"
    // Control arm
    local p "`p' (`cost_tx_0') (`cost_nt_0') (`cost_total_0')"
    local p "`p' (`cost_tx_d_0') (`cost_nt_d_0') (`cost_total_d_0')"
    // QALYs
    local p "`p' (`qaly_total_1') (`qaly_total_d_1')"
    local p "`p' (`qaly_total_0') (`qaly_total_d_0')"
    // OS and TXD
	local p "`p' (`os_1') (`os_0') (`txd_L${line}_1') (`txd_L${line}_0') (`txd_L${line}_med_1') (`txd_L${line}_med_0')"
    
    frame post bs_results `p'
end

**********
// Save Summary Dataset
capture program drop save_summary_dataset
program define save_summary_dataset
    // Merge pathways data into results
    frame change bs_results
    
    // Merge pathways into results using frlink
    frlink 1:1 sample, frame(bs_pathways)
    
    // Get all variables from bs_pathways except sample
    frame bs_pathways: ds
    local all_pathway_vars = r(varlist)
    local pathway_vars : list all_pathway_vars - sample
    
    frget `pathway_vars', from(bs_pathways)
    drop bs_pathways
    
    // Generate incremental outcomes (calculated once here)
    gen inc_cost_d = cost_total_d_1 - cost_total_d_0
    gen inc_qaly_d = qaly_total_d_1 - qaly_total_d_0
    gen icer = inc_cost_d / inc_qaly_d if abs(inc_qaly_d) > 1e-10
    gen inc_os = os_1 - os_0
    gen inc_cost = cost_total_1 - cost_total_0
    gen inc_qaly = qaly_total_1 - qaly_total_0
    gen inc_cost_tx_d = cost_tx_d_1 - cost_tx_d_0
    gen inc_cost_nt_d = cost_nt_d_1 - cost_nt_d_0
    
    // Formats
    format cost_* inc_cost* %12.0fc
    format qaly_* inc_qaly* %6.3f
    format icer %12.0fc
    format os_* inc_os %6.2f
    
    cap mkdir "$output_path"
    save "$output_path/bootstrap_${analysis}_${line}.dta", replace
    di "Summary saved: `c(N)' bootstrap iterations"
    frame change default
end

**********
// Generate Summary Statistics
capture program drop generate_summary_statistics
program define generate_summary_statistics
    frame change bs_results
    di _n "{hline 80}" _n "COST-EFFECTIVENESS RESULTS" _n "{hline 80}"
    display_group_stats
    
    if $OUT_SURVIVAL di _n "{hline 80}" _n "SURVIVAL OUTCOMES" _n "{hline 80}"
    if $OUT_SURVIVAL display_survival_stats
    
    frame change default
end

capture program drop display_group_stats
program define display_group_stats
    qui sum inc_cost_d
    local icm = r(mean)
    _pctile inc_cost_d, p(2.5 97.5)
    local icl = r(r1)
    local ich = r(r2)
    
    qui sum inc_qaly_d
    local iqm = r(mean)
    _pctile inc_qaly_d, p(2.5 97.5)
    local iql = r(r1)
    local iqh = r(r2)
    
    qui sum icer
    local icerm = r(mean)
    _pctile icer, p(2.5 97.5)
    local icerl = r(r1)
    local icerh = r(r2)
    
    di _n "All Patients" _n "{hline 60}"
    di "Inc Cost:  " %12.0fc `icm' " (" %12.0fc `icl' ", " %12.0fc `ich' ")"
    di "Inc QALYs: " %6.3f `iqm' " (" %6.3f `iql' ", " %6.3f `iqh' ")"
    di "ICER:      " %12.0fc `icerm' " (" %12.0fc `icerl' ", " %12.0fc `icerh' ")"
end

capture program drop display_survival_stats
program define display_survival_stats
    qui sum os_0
    local os0m = r(mean)
    _pctile os_0, p(2.5 97.5)
    local os0l = r(r1)
    local os0h = r(r2)
    
    qui sum os_1
    local os1m = r(mean)
    _pctile os_1, p(2.5 97.5)
    local os1l = r(r1)
    local os1h = r(r2)
    
    qui sum inc_os
    local iosm = r(mean)
    _pctile inc_os, p(2.5 97.5)
    local iosl = r(r1)
    local iosh = r(r2)
    
    di _n "All Patients" _n "{hline 60}"
    di "OS ($int0_name): " %6.2f `os0m' " mths (" %6.2f `os0l' ", " %6.2f `os0h' ")"
    di "OS ($int1_name): " %6.2f `os1m' " mths (" %6.2f `os1l' ", " %6.2f `os1h' ")"
    di "Inc OS:          " %6.2f `iosm' " mths (" %6.2f `iosl' ", " %6.2f `iosh' ")"
end

**********
// Analyse Treatment Pathways
capture program drop analyse_treatment_pathways
program define analyse_treatment_pathways
    frame change bs_results
    di _n "{hline 90}" _n "LINE $line BCR (95% CI)" _n "{hline 90}"
    foreach r in cr vg pr mr sd pd {
        qui sum int1_`r'
        local i1m = r(mean)
        _pctile int1_`r', p(2.5 97.5)
        local i1l = r(r1)
        local i1h = r(r2)
        qui sum int0_`r'
        local i0m = r(mean)
        _pctile int0_`r', p(2.5 97.5)
        local i0l = r(r1)
        local i0h = r(r2)
        cap drop diff_`r'
        gen diff_`r' = int1_`r' - int0_`r'
        qui sum diff_`r'
        local dm = r(mean)
        _pctile diff_`r', p(2.5 97.5)
        local dl = r(r1)
        local dh = r(r2)
        local nm ""
        if "`r'"=="cr" local nm "CR"
        if "`r'"=="vg" local nm "VG"
        if "`r'"=="pr" local nm "PR"
        if "`r'"=="mr" local nm "MR"
        if "`r'"=="sd" local nm "SD"
        if "`r'"=="pd" local nm "PD"
		di %-6s "`nm'" %5.1f `i1m' "% (" %4.1f `i1l' "-" %4.1f `i1h' ")  " %5.1f `i0m' "% (" %4.1f `i0l' "-" %4.1f `i0h' ")  " %5.1f `dm' "% (" %4.1f `dl' "," %4.1f `dh' ")"
    }
    
    // Display subsequent LOT receipt (from line+1 to 9)
    local start_lot = $line + 1
    if `start_lot' <= 9 {
        di _n "{hline 90}" _n "SUBSEQUENT LOT RECEIPT (Lines `start_lot'-9)" _n "{hline 90}"
        forval lot = `start_lot'/9 {
            qui sum int1_lot`lot'
            local i1m = r(mean)
            _pctile int1_lot`lot', p(2.5 97.5)
            local i1l = r(r1)
            local i1h = r(r2)
            qui sum int0_lot`lot'
            local i0m = r(mean)
            _pctile int0_lot`lot', p(2.5 97.5)
            local i0l = r(r1)
            local i0h = r(r2)
            cap drop diff_lot`lot'
            gen diff_lot`lot' = int1_lot`lot' - int0_lot`lot'
            qui sum diff_lot`lot'
            local dm = r(mean)
            _pctile diff_lot`lot', p(2.5 97.5)
            local dl = r(r1)
            local dh = r(r2)
			di "L`lot' " %7.0f `i1m' " (" %7.0f `i1l' "-" %7.0f `i1h' ")  " %7.0f `i0m' " (" %7.0f `i0l' "-" %7.0f `i0h' ")  " %7.0f `dm' " (" %7.0f `dl' "," %7.0f `dh' ")"
        }
    }
    frame change default
end

results, maxbs(500)
