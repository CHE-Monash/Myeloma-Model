**********
* EpiMAP Myeloma - Process Bootstraps
*
* Purpose: Process bootstrap simulation results 
* Usage: Set configuration, then run: results, maxbs(500)
* Output: Summary dataset
**********

clear all

**********
// Configuration
**********
global analysis "td_l4_pre"
global int1_name "Tec-Dara"
global int0_name "SoC"
global int1_file "td"
global int0_file "soc"
global start_line 4
global input_path "analyses/$analysis/simulated/bootstrap"
global output_path "analyses/$analysis/simulated"
global file_name "predicted_B"
global drate = 0.05

// Subgroups Flags (1=include, 0=exclude)
global SG_AGE70 = 0
global SG_AGE75 = 0
global SG_SEX = 0
global SG_ISS = 0
global SG_ECOG = 0
global SG_RENAL = 0

// Output Flags
global OUT_SUMMARY = 1
global OUT_PATHWAYS = 0
global OUT_SURVIVAL = 0

**********
// Main Program
capture program drop results
program define results
    syntax [, maxbs(integer 500)]
    di _n "{hline 80}" _n "Processing Bootstrap Results: $int1_name vs $int0_name (Line $start_line)" _n "{hline 80}"
    global MaxBS = `maxbs'
    global Line = $start_line
    initialise_frames
    process_bootstrap_samples
    save_summary_dataset
    if $OUT_SUMMARY generate_summary_statistics
    if $OUT_PATHWAYS analyse_treatment_pathways
    di _n "{hline 80}" _n "Complete. Summary: $output_path/summary_$analysis.dta" _n "{hline 80}"
end

**********
// Initialise Frames
capture program drop initialise_frames
program define initialise_frames
    local v "sample cTX1 cNT1 cTotal1 cTXd1 cNTd1 cTotald1 cTX0 cNT0 cTotal0 cTXd0 cNTd0 cTotald0"
	local v "`v' qTotal1 qTotald1 qTotal0 qTotald0 os1 os0 txd_l4_1 txd_l4_0"
    if $SG_AGE70 local v "`v' q1_a70 q0_a70 c1_a70 c0_a70 os1_a70 os0_a70 q1_a70n q0_a70n c1_a70n c0_a70n os1_a70n os0_a70n"
    if $SG_AGE75 local v "`v' q1_a75 q0_a75 c1_a75 c0_a75 os1_a75 os0_a75 q1_a75n q0_a75n c1_a75n c0_a75n os1_a75n os0_a75n"
    if $SG_SEX local v "`v' q1_m q0_m c1_m c0_m os1_m os0_m q1_f q0_f c1_f c0_f os1_f os0_f"
    if $SG_ISS local v "`v' q1_iss1 q0_iss1 c1_iss1 c0_iss1 os1_iss1 os0_iss1 q1_iss2 q0_iss2 c1_iss2 c0_iss2 os1_iss2 os0_iss2 q1_iss3 q0_iss3 c1_iss3 c0_iss3 os1_iss3 os0_iss3"
    if $SG_ECOG local v "`v' q1_ecog0 q0_ecog0 c1_ecog0 c0_ecog0 os1_ecog0 os0_ecog0 q1_ecog1 q0_ecog1 c1_ecog1 c0_ecog1 os1_ecog1 os0_ecog1 q1_ecog2 q0_ecog2 c1_ecog2 c0_ecog2 os1_ecog2 os0_ecog2"
    if $SG_RENAL local v "`v' q1_ckd q0_ckd c1_ckd c0_ckd os1_ckd os0_ckd q1_nockd q0_nockd c1_nockd c0_nockd os1_nockd os0_nockd"
    cap frame drop bs_results
    frame create bs_results `v'
    cap frame drop bs_pathways
    frame create bs_pathways sample int1_n int0_n int1_cr int1_vgpr int1_pr int1_mr int1_sd int1_pd int0_cr int0_vgpr int0_pr int0_mr int0_sd int0_pd int1_lot5 int0_lot5 int1_lot6 int0_lot6 int1_lot7 int0_lot7 int1_lot8 int0_lot8 int1_lot9 int0_lot9
end

**********
// Process Bootstrap Samples
capture program drop process_bootstrap_samples
program define process_bootstrap_samples
    forval b = 1/$MaxBS {
        if mod(`b',50)==0 | `b'==1 di "  Processing Sample `b' of $MaxBS..." 
        local f0 "$input_path/${int0_file}_${start_line}_${file_name}`b'.dta"
        cap use "`f0'", clear
        if _rc di "  NOT FOUND"
        gen Int = 0
        extract_arm_statistics `b' 0
		local ctrl "$input_path/temp_ctrl_`b'.dta"
        qui save `ctrl', replace
        local f1 "$input_path/${int1_file}_${start_line}_${file_name}`b'.dta"
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
        local L = $start_line
        count if BCR_L`L'==1
        local cr = r(N)/_N*100
        count if BCR_L`L'==2
        local vgpr = r(N)/_N*100
        count if BCR_L`L'==3
        local pr = r(N)/_N*100
        count if BCR_L`L'==4
        local mr = r(N)/_N*100
        count if BCR_L`L'==5
        local sd = r(N)/_N*100
        count if BCR_L`L'==6
        local pd = r(N)/_N*100
        forval lot=5/9 {
            cap count if TXD_L`lot'!=. & TXD_L`lot'>0
            local lot`lot' = cond(_rc==0, r(N)/_N*100, 0)
        }
    }
	if `arm'==0 foreach x in n cr vgpr pr mr sd pd lot5 lot6 lot7 lot8 lot9 {
			global int0_`x' = ``x''
		}
	if `arm'==1 frame post bs_pathways (`b') (`n') ($int0_n) (`cr') (`vgpr') (`pr') (`mr') (`sd') (`pd') ($int0_cr) ($int0_vgpr) ($int0_pr) ($int0_mr) ($int0_sd) ($int0_pd) (`lot5') ($int0_lot5) (`lot6') ($int0_lot6) (`lot7') ($int0_lot7) (`lot8') ($int0_lot8) (`lot9') ($int0_lot9)
end

**********
// Extract Combined Metrics
capture program drop extract_combined_metrics
program define extract_combined_metrics
    args b
    qui {
        foreach i in 1 0 {
            sum cTX if Int==`i'
            local cTX`i' = r(mean)
            sum cNT if Int==`i'
            local cNT`i' = r(mean)
            sum cTotal if Int==`i'
            local cTotal`i' = r(mean)
            sum cTXd if Int==`i'
            local cTXd`i' = r(mean)
            sum cNTd if Int==`i'
            local cNTd`i' = r(mean)
            sum cTotald if Int==`i'
            local cTotald`i' = r(mean)
            sum qTotal if Int==`i'
            local qTotal`i' = r(mean)
            sum qTotald if Int==`i'
            local qTotald`i' = r(mean)
            sum OC_TIME_L4 if Int==`i'
            local os`i' = r(mean)
			sum TXD_L4 if Int==`i'
			local txd_l4_`i' = r(mean)
        }
    }
    
    if $SG_AGE70 qui foreach a in 1 0 {
        local suf = cond(`a'==1,"_a70","_a70n")
        foreach i in 1 0 {
            sum qTotald if Int==`i' & Age70==`a'
            local q`i'`suf' = r(mean)
            sum cTotald if Int==`i' & Age70==`a'
            local c`i'`suf' = r(mean)
            sum OC_TIME if Int==`i' & Age70==`a'
            local os`i'`suf' = r(mean)/12
        }
    }
    
    if $SG_AGE75 qui foreach a in 1 0 {
        local suf = cond(`a'==1,"_a75","_a75n")
        foreach i in 1 0 {
            sum qTotald if Int==`i' & Age75==`a'
            local q`i'`suf' = r(mean)
            sum cTotald if Int==`i' & Age75==`a'
            local c`i'`suf' = r(mean)
            sum OC_TIME if Int==`i' & Age75==`a'
            local os`i'`suf' = r(mean)/12
        }
    }
    
    if $SG_SEX qui foreach s in 1 0 {
        local suf = cond(`s'==1,"_m","_f")
        foreach i in 1 0 {
            sum qTotald if Int==`i' & Male==`s'
            local q`i'`suf' = r(mean)
            sum cTotald if Int==`i' & Male==`s'
            local c`i'`suf' = r(mean)
            sum OC_TIME if Int==`i' & Male==`s'
            local os`i'`suf' = r(mean)/12
        }
    }
    
    if $SG_ISS qui foreach iss in 1 2 3 {
        foreach i in 1 0 {
            sum qTotald if Int==`i' & RISS==`iss'
            local q`i'_iss`iss' = r(mean)
            sum cTotald if Int==`i' & RISS==`iss'
            local c`i'_iss`iss' = r(mean)
            sum OC_TIME if Int==`i' & RISS==`iss'
            local os`i'_iss`iss' = r(mean)/12
        }
    }
    
    if $SG_ECOG qui foreach e in 0 1 {
        foreach i in 1 0 {
            sum qTotald if Int==`i' & ECOGcc==`e'
            local q`i'_ecog`e' = r(mean)
            sum cTotald if Int==`i' & ECOGcc==`e'
            local c`i'_ecog`e' = r(mean)
            sum OC_TIME if Int==`i' & ECOGcc==`e'
            local os`i'_ecog`e' = r(mean)/12
        }
    }
    if $SG_ECOG qui foreach i in 1 0 {
        sum qTotald if Int==`i' & ECOGcc>=2
        local q`i'_ecog2 = r(mean)
        sum cTotald if Int==`i' & ECOGcc>=2
        local c`i'_ecog2 = r(mean)
        sum OC_TIME if Int==`i' & ECOGcc>=2
        local os`i'_ecog2 = r(mean)/12
    }
    
    if $SG_RENAL qui foreach r in 1 0 {
        local suf = cond(`r'==1,"_ckd","_nockd")
        foreach i in 1 0 {
            sum qTotald if Int==`i' & CKD==`r'
            local q`i'`suf' = r(mean)
            sum cTotald if Int==`i' & CKD==`r'
            local c`i'`suf' = r(mean)
            sum OC_TIME if Int==`i' & CKD==`r'
            local os`i'`suf' = r(mean)/12
        }
    }
    
	local p "(`b') (`cTX1') (`cNT1') (`cTotal1') (`cTXd1') (`cNTd1') (`cTotald1') (`cTX0') (`cNT0') (`cTotal0') (`cTXd0') (`cNTd0') (`cTotald0') (`qTotal1') (`qTotald1') (`qTotal0') (`qTotald0') (`os1') (`os0') (`txd_l4_1') (`txd_l4_0')"
    if $SG_AGE70 local p "`p' (`q1_a70') (`q0_a70') (`c1_a70') (`c0_a70') (`os1_a70') (`os0_a70') (`q1_a70n') (`q0_a70n') (`c1_a70n') (`c0_a70n') (`os1_a70n') (`os0_a70n')"
    if $SG_AGE75 local p "`p' (`q1_a75') (`q0_a75') (`c1_a75') (`c0_a75') (`os1_a75') (`os0_a75') (`q1_a75n') (`q0_a75n') (`c1_a75n') (`c0_a75n') (`os1_a75n') (`os0_a75n')"
    if $SG_SEX local p "`p' (`q1_m') (`q0_m') (`c1_m') (`c0_m') (`os1_m') (`os0_m') (`q1_f') (`q0_f') (`c1_f') (`c0_f') (`os1_f') (`os0_f')"
    if $SG_ISS local p "`p' (`q1_iss1') (`q0_iss1') (`c1_iss1') (`c0_iss1') (`os1_iss1') (`os0_iss1') (`q1_iss2') (`q0_iss2') (`c1_iss2') (`c0_iss2') (`os1_iss2') (`os0_iss2') (`q1_iss3') (`q0_iss3') (`c1_iss3') (`c0_iss3') (`os1_iss3') (`os0_iss3')"
    if $SG_ECOG local p "`p' (`q1_ecog0') (`q0_ecog0') (`c1_ecog0') (`c0_ecog0') (`os1_ecog0') (`os0_ecog0') (`q1_ecog1') (`q0_ecog1') (`c1_ecog1') (`c0_ecog1') (`os1_ecog1') (`os0_ecog1') (`q1_ecog2') (`q0_ecog2') (`c1_ecog2') (`c0_ecog2') (`os1_ecog2') (`os0_ecog2')"
    if $SG_RENAL local p "`p' (`q1_ckd') (`q0_ckd') (`c1_ckd') (`c0_ckd') (`os1_ckd') (`os0_ckd') (`q1_nockd') (`q0_nockd') (`c1_nockd') (`c0_nockd') (`os1_nockd') (`os0_nockd')"
    frame post bs_results `p'
end

**********
// Save Summary Dataset
capture program drop save_summary_dataset
program define save_summary_dataset
    frame change bs_results 
    gen ic = cTotald1 - cTotald0
    gen iq = qTotald1 - qTotald0
    gen icer = ic/iq if abs(iq)>1e-10
    gen ios = os1 - os0
    gen ic_undis = cTotal1 - cTotal0
    gen iq_undis = qTotal1 - qTotal0
    
    if $SG_AGE70 foreach s in a70 a70n {
        gen ic_`s' = c1_`s' - c0_`s'
        gen iq_`s' = q1_`s' - q0_`s'
        gen icer_`s' = ic_`s'/iq_`s' if abs(iq_`s')>1e-10
        gen ios_`s' = os1_`s' - os0_`s'
    }
    
    if $SG_AGE75 foreach s in a75 a75n {
        gen ic_`s' = c1_`s' - c0_`s'
        gen iq_`s' = q1_`s' - q0_`s'
        gen icer_`s' = ic_`s'/iq_`s' if abs(iq_`s')>1e-10
        gen ios_`s' = os1_`s' - os0_`s'
    }
    
    if $SG_SEX foreach s in m f {
        gen ic_`s' = c1_`s' - c0_`s'
        gen iq_`s' = q1_`s' - q0_`s'
        gen icer_`s' = ic_`s'/iq_`s' if abs(iq_`s')>1e-10
        gen ios_`s' = os1_`s' - os0_`s'
    }
    
    if $SG_ISS foreach s in iss1 iss2 iss3 {
        gen ic_`s' = c1_`s' - c0_`s'
        gen iq_`s' = q1_`s' - q0_`s'
        gen icer_`s' = ic_`s'/iq_`s' if abs(iq_`s')>1e-10
        gen ios_`s' = os1_`s' - os0_`s'
    }
    
    if $SG_ECOG foreach s in ecog0 ecog1 ecog2 {
        gen ic_`s' = c1_`s' - c0_`s'
        gen iq_`s' = q1_`s' - q0_`s'
        gen icer_`s' = ic_`s'/iq_`s' if abs(iq_`s')>1e-10
        gen ios_`s' = os1_`s' - os0_`s'
    }
    
    if $SG_RENAL foreach s in ckd nockd {
        gen ic_`s' = c1_`s' - c0_`s'
        gen iq_`s' = q1_`s' - q0_`s'
        gen icer_`s' = ic_`s'/iq_`s' if abs(iq_`s')>1e-10
        gen ios_`s' = os1_`s' - os0_`s'
    }
    
    format c* ic* %12.0fc
    format q* iq* %6.3f
    format icer* %12.0fc
    format os* ios* %6.2f
    label var sample "Bootstrap sample"
    label var ic "Incremental cost (discounted)"
    label var iq "Incremental QALYs (discounted)"
    label var icer "ICER ($/QALY)"
    label var ios "Incremental OS (years)"
    cap mkdir "$output_path"
    save "$output_path/summary_$analysis.dta", replace
    di "Summary saved: `c(N)' bootstrap iterations"
    frame change default
end

**********
// Generate Summary Statistics
capture program drop generate_summary_statistics
program define generate_summary_statistics
    frame change bs_results
    di _n "{hline 80}" _n "COST-EFFECTIVENESS RESULTS" _n "{hline 80}"
    display_group_stats "" "All Patients"
    
    if $SG_ISS display_group_stats "_iss1" "ISS Stage I"
    if $SG_ISS display_group_stats "_iss2" "ISS Stage II"
    if $SG_ISS display_group_stats "_iss3" "ISS Stage III"
    
    if $SG_ECOG display_group_stats "_ecog0" "ECOG 0"
    if $SG_ECOG display_group_stats "_ecog1" "ECOG 1"
    if $SG_ECOG display_group_stats "_ecog2" "ECOG 2+"
    
    if $OUT_SURVIVAL di _n "{hline 80}" _n "SURVIVAL OUTCOMES" _n "{hline 80}"
    if $OUT_SURVIVAL display_survival_stats
    
    frame change default
end

capture program drop display_group_stats
program define display_group_stats
    args suf lbl
    local c0 = cond("`suf'"=="","cTotald0","c0`suf'")
    local c1 = cond("`suf'"=="","cTotald1","c1`suf'")
    local q0 = cond("`suf'"=="","qTotald0","q0`suf'")
    local q1 = cond("`suf'"=="","qTotald1","q1`suf'")
    qui sum ic`suf'
    local icm = r(mean)
    _pctile ic`suf', p(2.5 97.5)
    local icl = r(r1)
    local ich = r(r2)
    qui sum iq`suf'
    local iqm = r(mean)
    _pctile iq`suf', p(2.5 97.5)
    local iql = r(r1)
    local iqh = r(r2)
    qui sum icer`suf'
    local icerm = r(mean)
    _pctile icer`suf', p(2.5 97.5)
    local icerl = r(r1)
    local icerh = r(r2)
    di _n "`lbl'" _n "{hline 60}"
    di "Inc Cost:  " %12.0fc `icm' " (" %12.0fc `icl' ", " %12.0fc `ich' ")"
    di "Inc QALYs: " %6.3f `iqm' " (" %6.3f `iql' ", " %6.3f `iqh' ")"
    di "ICER:      " %12.0fc `icerm' " (" %12.0fc `icerl' ", " %12.0fc `icerh' ")"
end

capture program drop display_survival_stats
program define display_survival_stats
    qui sum os0
    local os0m = r(mean)
    _pctile os0, p(2.5 97.5)
    local os0l = r(r1)
    local os0h = r(r2)
    qui sum os1
    local os1m = r(mean)
    _pctile os1, p(2.5 97.5)
    local os1l = r(r1)
    local os1h = r(r2)
    qui sum ios
    local iosm = r(mean)
    _pctile ios, p(2.5 97.5)
    local iosl = r(r1)
    local iosh = r(r2)
    di _n "All Patients" _n "{hline 60}"
    di "OS ($int0_name): " %6.2f `os0m' " yrs (" %6.2f `os0l' ", " %6.2f `os0h' ")"
    di "OS ($int1_name): " %6.2f `os1m' " yrs (" %6.2f `os1l' ", " %6.2f `os1h' ")"
    di "Inc OS:          " %6.2f `iosm' " yrs (" %6.2f `iosl' ", " %6.2f `iosh' ")"
    if $SG_ISS foreach iss in 1 2 3 {
		di _n "ISS Stage `iss'" _n "{hline 60}"
        qui sum os0_iss`iss'
        local os0m = r(mean)
        _pctile os0_iss`iss', p(2.5 97.5)
        local os0l = r(r1)
        local os0h = r(r2)
        qui sum os1_iss`iss'
        local os1m = r(mean)
        _pctile os1_iss`iss', p(2.5 97.5)
        local os1l = r(r1)
        local os1h = r(r2)
        qui sum ios_iss`iss'
        local iosm = r(mean)
        _pctile ios_iss`iss', p(2.5 97.5)
        local iosl = r(r1)
        local iosh = r(r2)
        di "OS ($int0_name): " %6.2f `os0m' " yrs (" %6.2f `os0l' ", " %6.2f `os0h' ")"
        di "OS ($int1_name): " %6.2f `os1m' " yrs (" %6.2f `os1l' ", " %6.2f `os1h' ")"
        di "Inc OS:          " %6.2f `iosm' " yrs (" %6.2f `iosl' ", " %6.2f `iosh' ")"
    }
    if $SG_ECOG foreach e in 0 1 2 {
        local elbl = cond(`e'==2,"ECOG 2+","ECOG `e'")
        di _n "`elbl'" _n "{hline 60}"
        qui sum os0_ecog`e'
        local os0m = r(mean)
        _pctile os0_ecog`e', p(2.5 97.5)
        local os0l = r(r1)
        local os0h = r(r2)
        qui sum os1_ecog`e'
        local os1m = r(mean)
        _pctile os1_ecog`e', p(2.5 97.5)
        local os1l = r(r1)
        local os1h = r(r2)
        qui sum ios_ecog`e'
        local iosm = r(mean)
        _pctile ios_ecog`e', p(2.5 97.5)
        local iosl = r(r1)
        local iosh = r(r2)
        di "OS ($int0_name): " %6.2f `os0m' " yrs (" %6.2f `os0l' ", " %6.2f `os0h' ")"
        di "OS ($int1_name): " %6.2f `os1m' " yrs (" %6.2f `os1l' ", " %6.2f `os1h' ")"
        di "Inc OS:          " %6.2f `iosm' " yrs (" %6.2f `iosl' ", " %6.2f `iosh' ")"
    }
end

**********
// Analyse Treatment Pathways
capture program drop analyse_treatment_pathways
program define analyse_treatment_pathways
    frame change bs_pathways
		save "$output_path/Pathways_$analysis.dta", replace
        di _n "{hline 90}" _n "LINE $start_line BCR (95% CI)" _n "{hline 90}"
        foreach r in cr vgpr pr mr sd pd {
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
            if "`r'"=="vgpr" local nm "VGPR"
            if "`r'"=="pr" local nm "PR"
            if "`r'"=="mr" local nm "MR"
            if "`r'"=="sd" local nm "SD"
            if "`r'"=="pd" local nm "PD"
            di %-6s "`nm'" %5.1f `i1m' "% (" %4.1f `i1l' "-" %4.1f `i1h' ")  " %5.1f `i0m' "% (" %4.1f `i0l' "-" %4.1f `i0h' ")  " %+5.1f `dm' "% (" %+4.1f `dl' "," %+4.1f `dh' ")"
        }
        di _n "{hline 90}" _n "SUBSEQUENT LOT RECEIPT" _n "{hline 90}"
        forval lot=5/9 {
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
            di "L`lot' " %5.1f `i1m' "% (" %4.1f `i1l' "-" %4.1f `i1h' ")  " %5.1f `i0m' "% (" %4.1f `i0l' "-" %4.1f `i0h' ")  " %+5.1f `dm' "% (" %+4.1f `dl' "," %+4.1f `dh' ")"
        }
    frame change default
end

results, maxbs(500)
