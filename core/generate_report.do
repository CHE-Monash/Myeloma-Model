**********
* EpiMAP Myeloma - Generate Report
* 
* Purpose: Create PDF report using putpdf
*          Integrates with EpiMAP v2.0 dispatcher system
* 
* Output: ${int}_${data}_${scenario}.pdf
**********


**********
* Create Output Directory
capture mkdir "$simulated_path/report"
local report_dir "$simulated_path/report"

* Load Data
qui use "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}_${scenario}.dta", clear

**********
* Start PDF
**********

capture putpdf clear
set graphics off
putpdf begin

// Title page
putpdf paragraph, halign(center)
putpdf text ("EpiMAP Myeloma v2.0"), bold font(,18)

putpdf paragraph
putpdf text ("Simulation Report"), bold font(,16)

putpdf table settings = (8, 2), border(all)
putpdf table settings(1,1) = ("Setting"), bold
putpdf table settings(1,2) = ("Value"), bold
putpdf table settings(2,1) = ("Analysis")
putpdf table settings(2,2) = ("$analysis")
putpdf table settings(3,1) = ("Intervention")
putpdf table settings(3,2) = ("$int")
putpdf table settings(4,1) = ("Line")
putpdf table settings(4,2) = ("$line")
putpdf table settings(5,1) = ("Data")
putpdf table settings(5,2) = ("$data")
putpdf table settings(6,1) = ("Patient IDs")
putpdf table settings(6,2) = ("$min_id to $max_id")
putpdf table settings(7,1) = ("Scenario")
putpdf table settings(7,2) = ("$scenario")
putpdf table settings(8,1) = ("Report Date")
putpdf table settings(8,2) = ("`c(current_date)'")

**********
* Patients
**********

putpdf paragraph
putpdf text ("Patients"), bold font(,16)

// Sample size
quietly count
local total_n = string(r(N), "%9.0fc")

// Sex
quietly count if Male == 1
local male_n = string(r(N), "%9.0fc")
local male_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("Sample Size: `total_n'"), linebreak
putpdf text ("Male: `male_n' (`male_pct'%)")

// Age
putpdf paragraph
putpdf text ("Age at Diagnosis"), bold

quietly summarize Age_DN, detail
local mean_age = string(r(mean), "%4.1f")
local sd_age = string(r(sd), "%4.1f")
local median_age = string(r(p50), "%4.1f")
local p25_age = string(r(p25), "%4.1f")
local p75_age = string(r(p75), "%4.1f")

putpdf table age_sum = (3, 2), border(all)
putpdf table age_sum(1,1) = ("Statistic"), bold
putpdf table age_sum(1,2) = ("Years"), bold
putpdf table age_sum(2,1) = ("Mean (SD)")
putpdf table age_sum(2,2) = ("`mean_age' (`sd_age')")
putpdf table age_sum(3,1) = ("Median [IQR]")
putpdf table age_sum(3,2) = ("`median_age' [`p25_age' - `p75_age']")

putpdf paragraph

// ECOG
quietly count if ECOGcc == 0
local ecog0_n = string(r(N), "%9.0fc")
local ecog0_pct = string(100*r(N)/_N, "%4.1f")
quietly count if ECOGcc == 1
local ecog1_n = string(r(N), "%9.0fc")
local ecog1_pct = string(100*r(N)/_N, "%4.1f")
quietly count if ECOGcc == 2
local ecog2_n = string(r(N), "%9.0fc")
local ecog2_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("ECOG"), bold

putpdf table ecog_tbl = (4, 3), border(all)
putpdf table ecog_tbl(1,1) = ("ECOG Score"), bold
putpdf table ecog_tbl(1,2) = ("N"), bold
putpdf table ecog_tbl(1,3) = ("%"), bold
putpdf table ecog_tbl(2,1) = ("0")
putpdf table ecog_tbl(2,2) = ("`ecog0_n'")
putpdf table ecog_tbl(2,3) = ("`ecog0_pct'")
putpdf table ecog_tbl(3,1) = ("1")
putpdf table ecog_tbl(3,2) = ("`ecog1_n'")
putpdf table ecog_tbl(3,3) = ("`ecog1_pct'")
putpdf table ecog_tbl(4,1) = ("2+")
putpdf table ecog_tbl(4,2) = ("`ecog2_n'")
putpdf table ecog_tbl(4,3) = ("`ecog2_pct'")


// R-ISS
quietly count if RISS == 1
local riss1_n = string(r(N), "%9.0fc")
local riss1_pct = string(100*r(N)/_N, "%4.1f")
quietly count if RISS == 2
local riss2_n = string(r(N), "%9.0fc")
local riss2_pct = string(100*r(N)/_N, "%4.1f")
quietly count if RISS == 3
local riss3_n = string(r(N), "%9.0fc")
local riss3_pct = string(100*r(N)/_N, "%4.1f")

putpdf paragraph
putpdf text ("R-ISS"), bold

putpdf table riss_tbl = (4, 3), border(all)
putpdf table riss_tbl(1,1) = ("R-ISS Stage"), bold
putpdf table riss_tbl(1,2) = ("N"), bold
putpdf table riss_tbl(1,3) = ("%"), bold
putpdf table riss_tbl(2,1) = ("I")
putpdf table riss_tbl(2,2) = ("`riss1_n'")
putpdf table riss_tbl(2,3) = ("`riss1_pct'")
putpdf table riss_tbl(3,1) = ("II")
putpdf table riss_tbl(3,2) = ("`riss2_n'")
putpdf table riss_tbl(3,3) = ("`riss2_pct'")
putpdf table riss_tbl(4,1) = ("III")
putpdf table riss_tbl(4,2) = ("`riss3_n'")
putpdf table riss_tbl(4,3) = ("`riss3_pct'")

putpdf pagebreak

**********
* Treatment
**********

putpdf paragraph
putpdf text ("Treatments"), bold font(,16)

// Line 1
quietly count if TXR_L1 < .
local l1_total = r(N)
    
	// Get regimen codes from Mata
	mata: st_matrix("regimen_codes", oL1_TXR)
	mata: st_numscalar("n_regimens", cols(oL1_TXR))
	local n_regimens = n_regimens 

	// Extract codes from matrix to locals
	forval i = 1/`n_regimens' {
		local code_`i' = regimen_codes[1,`i']
	}

	// Count each specific regimen and store its name
	forval i = 1/`n_regimens' {
		local code =`code_`i''
		quietly count if TXR_L1 == `code'
		local reg_`i'_n = string(r(N), "%9.0fc")
		local reg_`i'_pct = string(100*r(N)/`l1_total', "%4.1f")
	  
		// Assign name based on code
		if `code' == 0  local reg_`i'_name "Other"
		if `code' == 2  local reg_`i'_name "TCd"
		if `code' == 4  local reg_`i'_name "VCd"
		if `code' == 7  local reg_`i'_name "Rd"
		if `code' == 9  local reg_`i'_name "VTd"
		if `code' == 31 local reg_`i'_name "VRd"
		if `code' == 49 local reg_`i'_name "Kd"
		if `code' == 56 local reg_`i'_name "Pd"
		if `code' == 80 local reg_`i'_name "DVd"
	}

	// Create table (n_regimens + 1 for header)
	local n_rows = `n_regimens' + 1
	putpdf table txr_l1_tbl = (`n_rows', 3), border(all)
	putpdf table txr_l1_tbl(1,1) = ("Line 1 Regimen"), bold
	putpdf table txr_l1_tbl(1,2) = ("N"), bold
	putpdf table txr_l1_tbl(1,3) = ("%"), bold

	// Fill in all regimens
	forval i = 1/`n_regimens' {
		local row = `i' + 1
		putpdf table txr_l1_tbl(`row',1) = ("`reg_`i'_name'")
		putpdf table txr_l1_tbl(`row',2) = ("`reg_`i'_n'")
		putpdf table txr_l1_tbl(`row',3) = ("`reg_`i'_pct'")
	}

// Autologous Stem Cell Transplant
qui count if SCT_L1 < .
local sct_total = r(N)

qui count if SCT_L1 == 1
local sct_n = string(r(N), "%9.0fc")
local sct_pct = string(100*r(N)/`sct_total', "%4.1f")
qui count if SCT_L1 == 0
local nosct_n = string(r(N), "%9.0fc")
local nosct_pct = string(100*r(N)/`sct_total', "%4.1f")
    
putpdf table sct_tbl = (3, 3), border(all)
putpdf table sct_tbl(1,1) = ("Received ASCT"), bold
putpdf table sct_tbl(1,2) = ("N"), bold
putpdf table sct_tbl(1,3) = ("%"), bold
putpdf table sct_tbl(2,1) = ("No")
putpdf table sct_tbl(2,2) = ("`nosct_n'")
putpdf table sct_tbl(2,3) = ("`nosct_pct'")
putpdf table sct_tbl(3,1) = ("Yes")
putpdf table sct_tbl(3,2) = ("`sct_n'")
putpdf table sct_tbl(3,3) = ("`sct_pct'")

// Line 2    
quietly count if TXR_L2 < .
local l2_total = r(N)
    
	// Get regimen codes from Mata
	mata: st_matrix("regimen_codes", oL2_TXR)
	mata: st_numscalar("n_regimens", cols(oL2_TXR))
	local n_regimens = n_regimens 

	// Extract codes from matrix to locals
	forval i = 1/`n_regimens' {
		local code_`i' = regimen_codes[1,`i']
	}

	// Count each specific regimen and store its name
	forval i = 1/`n_regimens' {
		local code =`code_`i''
		quietly count if TXR_L2 == `code'
		local reg_`i'_n = string(r(N), "%9.0fc")
		local reg_`i'_pct = string(100*r(N)/`l2_total', "%4.1f")
	  
		// Assign name based on code
		if `code' == 0  local reg_`i'_name "Other"
		if `code' == 2  local reg_`i'_name "TCd"
		if `code' == 4  local reg_`i'_name "VCd"
		if `code' == 5  local reg_`i'_name "Vd"
		if `code' == 7  local reg_`i'_name "Rd"
		if `code' == 9  local reg_`i'_name "VTd"
		if `code' == 31 local reg_`i'_name "VRd"
		if `code' == 49 local reg_`i'_name "Kd"
		if `code' == 56 local reg_`i'_name "Pd"
		if `code' == 80 local reg_`i'_name "DVd"
	}

	// Create table (n_regimens + 1 for header)
	local n_rows = `n_regimens' + 1
	putpdf table txr_l2_tbl = (`n_rows', 3), border(all)
	putpdf table txr_l2_tbl(1,1) = ("Line 2 Regimen"), bold
	putpdf table txr_l2_tbl(1,2) = ("N"), bold
	putpdf table txr_l2_tbl(1,3) = ("%"), bold

	// Fill in all regimens
	forval i = 1/`n_regimens' {
		local row = `i' + 1
		putpdf table txr_l2_tbl(`row',1) = ("`reg_`i'_name'")
		putpdf table txr_l2_tbl(`row',2) = ("`reg_`i'_n'")
		putpdf table txr_l2_tbl(`row',3) = ("`reg_`i'_pct'")
	}

// Best Clinical Response 
putpdf paragraph
putpdf text ("Best Clinical Response"), bold

	// Line 1
	forval b = 1/6 {
		qui count if BCR_L1 == `b'
		local `b'_n_l1 = string(r(N), "%9.0fc")
		local `b'_pct_l1 = string(100*r(N)/`l1_total', "%4.1f")
	}
  
	// ASCT
	qui count if SCT_L1 == 1
	local sct_total = r(N)
    forval b = 1/4 {
		qui count if BCR_SCT == `b'
		local `b'_n_asct = string(r(N), "%9.0fc")
		local `b'_pct_asct = string(100*r(N)/`sct_total', "%4.1f")
	}

    // Line 1
	forval b = 1/6 {
		qui count if BCR_L2 == `b'
		local `b'_n_l2 = string(r(N), "%9.0fc")
		local `b'_pct_l2 = string(100*r(N)/`l2_total', "%4.1f")
	}

    putpdf table bcr_tbl = (7, 4), border(all)
    
	putpdf table bcr_tbl(1,1) = ("Response"), bold
    putpdf table bcr_tbl(1,2) = ("Line 1"), bold
    putpdf table bcr_tbl(1,3) = ("ASCT"), bold
	putpdf table bcr_tbl(1,4) = ("Line 2"), bold
	
    putpdf table bcr_tbl(2,1) = ("Complete Response")
    putpdf table bcr_tbl(2,2) = ("`1_pct_l1'%")
	putpdf table bcr_tbl(2,3) = ("`1_pct_asct'%")
	putpdf table bcr_tbl(2,4) = ("`1_pct_l2'%")
	
    putpdf table bcr_tbl(3,1) = ("V.Good Partial Response")
    putpdf table bcr_tbl(3,2) = ("`2_pct_l1'%")
	putpdf table bcr_tbl(3,3) = ("`2_pct_asct'%")
	putpdf table bcr_tbl(3,4) = ("`2_pct_l2'%")
	
    putpdf table bcr_tbl(4,1) = ("Partial Response")
    putpdf table bcr_tbl(4,2) = ("`3_pct_l1'%")
	putpdf table bcr_tbl(4,3) = ("`3_pct_asct'%")
	putpdf table bcr_tbl(4,4) = ("`3_pct_l2'%")
	
    putpdf table bcr_tbl(5,1) = ("Minimal Response")
    putpdf table bcr_tbl(5,2) = ("`4_pct_l1'%")
	putpdf table bcr_tbl(5,3) = ("`4_pct_asct'%")
	putpdf table bcr_tbl(5,4) = ("`4_pct_l2'%")
	
    putpdf table bcr_tbl(6,1) = ("Stable Disease")
    putpdf table bcr_tbl(6,2) = ("`5_pct_l1'%")
	putpdf table bcr_tbl(6,4) = ("`5_pct_l2'%")
	
    putpdf table bcr_tbl(7,1) = ("Progressive Disease")
    putpdf table bcr_tbl(7,2) = ("`6_pct_l1'%")
	putpdf table bcr_tbl(7,4) = ("`6_pct_l2'%")

**********
* Overall Survival
**********

putpdf paragraph
putpdf text ("Overall Survival Results"), bold font(,16)

// Summary statistics
quietly summarize OC_TIME, detail
local mean = string(r(mean)/12, "%6.2f")
local sd = string(r(sd)/12, "%5.2f")
local median = string(r(p50)/12, "%6.2f")
local p25 = string(r(p25)/12, "%6.2f")
local p75 = string(r(p75)/12, "%6.2f")

putpdf paragraph
putpdf text ("Summary Statistics"), bold

putpdf table os_sum = (3, 2), border(all)
putpdf table os_sum(1,1) = ("Statistic"), bold
putpdf table os_sum(1,2) = ("Years"), bold
putpdf table os_sum(2,1) = ("Mean (SD)")
putpdf table os_sum(2,2) = ("`mean' (`sd')")
putpdf table os_sum(3,1) = ("Median [IQR]")
putpdf table os_sum(3,2) = ("`median' [`p25'-`p75']")

// Survival at key time points
putpdf paragraph
putpdf text ("Survival at Key Time Points"), bold

putpdf table surv_time = (6, 2), border(all)
putpdf table surv_time(1,1) = ("Time Point"), bold
putpdf table surv_time(1,2) = ("Survival % (95% CI)"), bold

local row = 2
foreach year in 1 2 3 5 10 {
    quietly count if OC_TIME/12 >= `year'
    local pct = (r(N) / _N) * 100    
    local pct_str = string(`pct', "%5.1f")
    
    putpdf table surv_time(`row',1) = ("`year'-year")
    putpdf table surv_time(`row',2) = ("`pct_str'%")
    local row = `row' + 1
}

// Generate and insert KM curves
capture mkdir "$simulated_path/report/figures"

// Overall KM
preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
sts graph, ///
    xtitle("Months") ytitle("Probability") title("") ///
    ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
    xlabel(0(24)240) ci risktable legend(off) ///
    graphregion(color(white)) name(os_overall, replace)
graph export "$simulated_path/report/figures/os_overall.png", replace width(1200)

putpdf paragraph
putpdf text ("Overall Survival"), bold linebreak(2)
putpdf image "$simulated_path/report/figures/os_overall.png", width(6)

// By ASCT
gen asct = SCT_L1
label define asct_lbl 0 "No ASCT" 1 "ASCT"
label values asct asct_lbl
    
sts graph, by(asct) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "No ASCT") label(2 "ASCT")) ///
	name(os_asct, replace)
graph export "$simulated_path/report/figures/os_asct.png", replace width(1200)
    
putpdf paragraph
putpdf text ("Overall Survival by ASCT Status"), bold
putpdf image "$simulated_path/report/figures/os_asct.png", width(6)
restore

// By BCR L1 / ASCT
preserve
gen bcr_group = .
replace bcr_group = 1 if SCT_L1 == 0 & BCR_L1 == 1 | SCT_L1 == 1 & BCR_SCT == 1
replace bcr_group = 2 if SCT_L1 == 0 & BCR_L1 == 2 | SCT_L1 == 1 & BCR_SCT == 2
replace bcr_group = 3 if SCT_L1 == 0 & BCR_L1 == 3 | SCT_L1 == 1 & BCR_SCT == 3
replace bcr_group = 4 if SCT_L1 == 0 & BCR_L1 == 4 | SCT_L1 == 1 & BCR_SCT == 4
replace bcr_group = 5 if BCR_L1 == 5
replace bcr_group = 6 if BCR_L1 == 6
    
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
sts graph, by(bcr_group) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "CR") label(2 "VGPR") label(3 "PR") label(4 "MR") label(5 "SD") label(6 "PD") rows(2)) ///
	name(os_bcr, replace)
graph export "$simulated_path/report/figures/os_bcr.png", replace width(1200)
   
putpdf pagebreak
putpdf paragraph
putpdf text ("Overall Survival by BCR LoT 1 / ASCT"), bold
putpdf image "$simulated_path/report/figures/os_bcr.png", width(6)
restore

// By Age
preserve
gen age_group = .
replace age_group = 1 if Age_DN < 65
replace age_group = 2 if Age_DN >= 65 & Age_DN < 75
replace age_group = 3 if Age_DN >= 75 & Age_DN < .

stset OC_TIME if OC_TIME < 240, failure(OC_MORT)

sts graph, by(age_group) ///
    xtitle("Months") ytitle("Probability") title("") ///
    ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
    graphregion(color(white)) ///
    legend(label(1 "<65") label(2 "65-74") label(3 "≥75") rows(1)) ///
    name(os_age, replace)
graph export "$simulated_path/report/figures/os_age.png", replace width(1200)

putpdf paragraph
putpdf text ("Overall Survival by Age Group"), bold
putpdf image "$simulated_path/report/figures/os_age.png", width(6)
restore

// By R-ISS
preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
sts graph, by(RISS) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "Stage I") label(2 "Stage II") label(3 "Stage III") rows(1)) ///
	name(os_riss, replace)
graph export "$simulated_path/report/figures/os_riss.png", replace width(1200)

putpdf pagebreak
putpdf paragraph
putpdf text ("Overall Survival by R-ISS Stage"), bold
putpdf image "$simulated_path/report/figures/os_riss.png", width(6)
restore

// By ECOG
preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
sts graph, by(ECOGcc) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "ECOG 0") label(2 "ECOG 1") label(3 "ECOG 2+") rows(1)) ///
	name(os_ecog, replace)
graph export "$simulated_path/report/figures/os_ecog.png", replace width(1200)
    
putpdf paragraph
putpdf text ("Overall Survival by ECOG Status"), bold
putpdf image "$simulated_path/report/figures/os_ecog.png", width(6)
restore

putpdf pagebreak

**********
* Economic Outcomes
**********

putpdf paragraph
putpdf text ("Economic Outcomes"), bold font(,16)

// Get discount rate for display
local drate_pct = string($drate * 100, "%3.1f")

**********
* Costs
**********

putpdf paragraph
putpdf text ("Costs (Discounted at `drate_pct'%)"), bold font(,14)

// Total Costs Summary
quietly summarize cTotald, detail
local n_cost = string(r(N), "%9.0fc")
local mean_cost = string(r(mean), "%12.0fc")
local sd_cost = string(r(sd), "%12.0fc")
local median_cost = string(r(p50), "%12.0fc")
local p25_cost = string(r(p25), "%12.0fc")
local p75_cost = string(r(p75), "%12.0fc")

putpdf paragraph
putpdf text ("Total Costs Summary"), bold

putpdf table cost_sum = (3, 2), border(all)
putpdf table cost_sum(1,1) = ("Statistic"), bold
putpdf table cost_sum(1,2) = ("AUD"), bold
putpdf table cost_sum(2,1) = ("Mean (SD)")
putpdf table cost_sum(2,2) = ("$`mean_cost' ($`sd_cost')")
putpdf table cost_sum(3,1) = ("Median [IQR]")
putpdf table cost_sum(3,2) = ("$`median_cost' [$`p25_cost' - $`p75_cost']")

// Cost Components
putpdf paragraph
putpdf text ("Cost Components (Mean)"), bold

quietly summarize cTXd
local mean_tx = string(r(mean), "%12.0fc")
quietly summarize cNTd
local mean_nt = string(r(mean), "%12.0fc")

// ASCT costs (recipients only)
quietly summarize cTX_ASCTd if SCT_L1 == 1
if r(N) > 0 {
	local mean_asct = string(r(mean), "%12.0fc")
	local n_asct = string(r(N), "%9.0fc")
}
else {
	local mean_asct = "N/A"
	local n_asct = "0"
}

// Maintenance costs (recipients only)
quietly summarize cTX_MNTd if MNT == 1
if r(N) > 0 {
	local mean_mnt = string(r(mean), "%12.0fc")
	local n_mnt = string(r(N), "%9.0fc")
}
else {
	local mean_mnt = "N/A"
	local n_mnt = "0"
}

// Non-treatment cost components
quietly summarize cNT_Hospd
local mean_hosp = string(r(mean), "%12.0fc")
quietly summarize cNT_Commd
local mean_comm = string(r(mean), "%12.0fc")
quietly summarize cNT_Emerd
local mean_emer = string(r(mean), "%12.0fc")

putpdf table cost_comp = (8, 2), border(all)
putpdf table cost_comp(1,1) = ("Component"), bold
putpdf table cost_comp(1,2) = ("Mean (AUD)"), bold
putpdf table cost_comp(2,1) = ("Treatment Costs (Total)")
putpdf table cost_comp(2,2) = ("$`mean_tx'")
putpdf table cost_comp(3,1) = ("  ASCT (n=`n_asct')")
putpdf table cost_comp(3,2) = ("$`mean_asct'")
putpdf table cost_comp(4,1) = ("  Maintenance (n=`n_mnt')")
putpdf table cost_comp(4,2) = ("$`mean_mnt'")
putpdf table cost_comp(5,1) = ("Non-Treatment Costs (Total)")
putpdf table cost_comp(5,2) = ("$`mean_nt'")
putpdf table cost_comp(6,1) = ("  Hospitalisation")
putpdf table cost_comp(6,2) = ("$`mean_hosp'")
putpdf table cost_comp(7,1) = ("  Community Care")
putpdf table cost_comp(7,2) = ("$`mean_comm'")
putpdf table cost_comp(8,1) = ("  Emergency")
putpdf table cost_comp(8,2) = ("$`mean_emer'")

// Treatment Costs by Line of Therapy
putpdf paragraph
putpdf text ("Treatment Costs by Line of Therapy"), bold

// Count rows needed (only include lines with patients)
local n_lines = 0
forval l = 1/9 {
	quietly count if cTX_L`l'd != . & cTX_L`l'd > 0
	if r(N) > 0 local n_lines = `l'
}

// Create table with header + lines
local n_rows = `n_lines' + 1
putpdf table cost_line = (`n_rows', 3), border(all)
putpdf table cost_line(1,1) = ("Line"), bold
putpdf table cost_line(1,2) = ("N Treated"), bold
putpdf table cost_line(1,3) = ("Mean Cost (AUD)"), bold

local row = 2
forval l = 1/`n_lines' {
	quietly count if cTX_L`l'd != . & cTX_L`l'd > 0
	local n_l = string(r(N), "%9.0fc")
	quietly summarize cTX_L`l'd if cTX_L`l'd > 0
	if r(N) > 0 {
		local mean_l = string(r(mean), "%12.0fc")
	}
	else {
		local mean_l = "0"
	}
	
	putpdf table cost_line(`row',1) = ("Line `l'")
	putpdf table cost_line(`row',2) = ("`n_l'")
	putpdf table cost_line(`row',3) = ("$`mean_l'")
	local row = `row' + 1
}

**********
* QALYs
**********

putpdf pagebreak
putpdf paragraph
putpdf text ("Quality-Adjusted Life Years (Discounted at `drate_pct'%)"), bold font(,14)

// Total QALYs Summary
quietly summarize qTotald, detail
local n_qaly = string(r(N), "%9.0fc")
local mean_qaly = string(r(mean), "%5.2f")
local sd_qaly = string(r(sd), "%5.2f")
local median_qaly = string(r(p50), "%5.2f")
local p25_qaly = string(r(p25), "%5.2f")
local p75_qaly = string(r(p75), "%5.2f")

putpdf paragraph
putpdf text ("Total QALYs Summary"), bold

putpdf table qaly_sum = (3, 2), border(all)
putpdf table qaly_sum(1,1) = ("Statistic"), bold
putpdf table qaly_sum(1,2) = ("QALYs"), bold
putpdf table qaly_sum(2,1) = ("Mean (SD)")
putpdf table qaly_sum(2,2) = ("`mean_qaly' (`sd_qaly')")
putpdf table qaly_sum(3,1) = ("Median [IQR]")
putpdf table qaly_sum(3,2) = ("`median_qaly' [`p25_qaly' - `p75_qaly']")

// QALY Components by Health State
putpdf paragraph
putpdf text ("QALYs by Health State (Mean)"), bold

quietly summarize qTFI_DNd
local q_tfi_dn = string(r(mean), "%5.3f")
quietly summarize qTXD_L1d
local q_txd_l1 = string(r(mean), "%5.3f")
quietly summarize qTFI_L1d
local q_tfi_l1 = string(r(mean), "%5.3f")
quietly summarize qTXD_L2d
local q_txd_l2 = string(r(mean), "%5.3f")
quietly summarize qPostL2d
local q_post = string(r(mean), "%5.3f")

putpdf table qaly_comp = (6, 3), border(all)
putpdf table qaly_comp(1,1) = ("Health State"), bold
putpdf table qaly_comp(1,2) = ("Utility Weight"), bold
putpdf table qaly_comp(1,3) = ("Mean QALYs"), bold
putpdf table qaly_comp(2,1) = ("TFI Pre-L1")
putpdf table qaly_comp(2,2) = ("0.72")
putpdf table qaly_comp(2,3) = ("`q_tfi_dn'")
putpdf table qaly_comp(3,1) = ("L1 Treatment")
putpdf table qaly_comp(3,2) = ("0.63")
putpdf table qaly_comp(3,3) = ("`q_txd_l1'")
putpdf table qaly_comp(4,1) = ("TFI Post-L1")
putpdf table qaly_comp(4,2) = ("0.72")
putpdf table qaly_comp(4,3) = ("`q_tfi_l1'")
putpdf table qaly_comp(5,1) = ("L2 Treatment")
putpdf table qaly_comp(5,2) = ("0.67")
putpdf table qaly_comp(5,3) = ("`q_txd_l2'")
putpdf table qaly_comp(6,1) = ("Post-L2")
putpdf table qaly_comp(6,2) = ("0.63")
putpdf table qaly_comp(6,3) = ("`q_post'")

**********
* Undiscounted vs Discounted Comparison
**********

putpdf paragraph
putpdf text ("Discounted vs Undiscounted Comparison"), bold font(,14)

// Get undiscounted values
quietly summarize cTotal
local mean_cost_undisc = string(r(mean), "%12.0fc")
quietly summarize cTotald
local mean_cost_disc = string(r(mean), "%12.0fc")

quietly summarize qTotal
local mean_qaly_undisc = string(r(mean), "%5.2f")
quietly summarize qTotald
local mean_qaly_disc = string(r(mean), "%5.2f")

putpdf table disc_comp = (3, 3), border(all)
putpdf table disc_comp(1,1) = ("Outcome"), bold
putpdf table disc_comp(1,2) = ("Undiscounted"), bold
putpdf table disc_comp(1,3) = ("Discounted (`drate_pct'%)"), bold
putpdf table disc_comp(2,1) = ("Mean Total Cost (AUD)")
putpdf table disc_comp(2,2) = ("$`mean_cost_undisc'")
putpdf table disc_comp(2,3) = ("$`mean_cost_disc'")
putpdf table disc_comp(3,1) = ("Mean QALYs")
putpdf table disc_comp(3,2) = ("`mean_qaly_undisc'")
putpdf table disc_comp(3,3) = ("`mean_qaly_disc'")

**********
* Costs and QALYs by Subgroup
**********

putpdf paragraph
putpdf text ("Economic Outcomes by Subgroup"), bold font(,14)

// By ASCT Status
putpdf paragraph
putpdf text ("By ASCT Status"), bold

quietly summarize cTotald if SCT_L1 == 0
local cost_noasct = string(r(mean), "%12.0fc")
quietly summarize qTotald if SCT_L1 == 0
local qaly_noasct = string(r(mean), "%5.2f")
quietly count if SCT_L1 == 0
local n_noasct = string(r(N), "%9.0fc")

quietly summarize cTotald if SCT_L1 == 1
local cost_asct = string(r(mean), "%12.0fc")
quietly summarize qTotald if SCT_L1 == 1
local qaly_asct = string(r(mean), "%5.2f")
quietly count if SCT_L1 == 1
local n_asct = string(r(N), "%9.0fc")

putpdf table asct_econ = (3, 4), border(all)
putpdf table asct_econ(1,1) = ("ASCT Status"), bold
putpdf table asct_econ(1,2) = ("N"), bold
putpdf table asct_econ(1,3) = ("Mean Cost (AUD)"), bold
putpdf table asct_econ(1,4) = ("Mean QALYs"), bold
putpdf table asct_econ(2,1) = ("No ASCT")
putpdf table asct_econ(2,2) = ("`n_noasct'")
putpdf table asct_econ(2,3) = ("$`cost_noasct'")
putpdf table asct_econ(2,4) = ("`qaly_noasct'")
putpdf table asct_econ(3,1) = ("ASCT")
putpdf table asct_econ(3,2) = ("`n_asct'")
putpdf table asct_econ(3,3) = ("$`cost_asct'")
putpdf table asct_econ(3,4) = ("`qaly_asct'")

// By Age Group
putpdf paragraph
putpdf text ("By Age Group"), bold

forval a = 1/3 {
	if `a' == 1 {
		local age_cond "Age_DN < 65"
		local age_lab "<65"
	}
	if `a' == 2 {
		local age_cond "Age_DN >= 65 & Age_DN < 75"
		local age_lab "65-74"
	}
	if `a' == 3 {
		local age_cond "Age_DN >= 75 & Age_DN < ."
		local age_lab "≥75"
	}
	
	quietly count if `age_cond'
	local n_age`a' = string(r(N), "%9.0fc")
	quietly summarize cTotald if `age_cond'
	local cost_age`a' = string(r(mean), "%12.0fc")
	quietly summarize qTotald if `age_cond'
	local qaly_age`a' = string(r(mean), "%5.2f")
}

putpdf table age_econ = (4, 4), border(all)
putpdf table age_econ(1,1) = ("Age Group"), bold
putpdf table age_econ(1,2) = ("N"), bold
putpdf table age_econ(1,3) = ("Mean Cost (AUD)"), bold
putpdf table age_econ(1,4) = ("Mean QALYs"), bold
putpdf table age_econ(2,1) = ("<65")
putpdf table age_econ(2,2) = ("`n_age1'")
putpdf table age_econ(2,3) = ("$`cost_age1'")
putpdf table age_econ(2,4) = ("`qaly_age1'")
putpdf table age_econ(3,1) = ("65-74")
putpdf table age_econ(3,2) = ("`n_age2'")
putpdf table age_econ(3,3) = ("$`cost_age2'")
putpdf table age_econ(3,4) = ("`qaly_age2'")
putpdf table age_econ(4,1) = ("≥75")
putpdf table age_econ(4,2) = ("`n_age3'")
putpdf table age_econ(4,3) = ("$`cost_age3'")
putpdf table age_econ(4,4) = ("`qaly_age3'")

// By R-ISS Stage
putpdf paragraph
putpdf text ("By R-ISS Stage"), bold

forval r = 1/3 {
	quietly count if RISS == `r'
	local n_riss`r' = string(r(N), "%9.0fc")
	quietly summarize cTotald if RISS == `r'
	local cost_riss`r' = string(r(mean), "%12.0fc")
	quietly summarize qTotald if RISS == `r'
	local qaly_riss`r' = string(r(mean), "%5.2f")
}

putpdf table riss_econ = (4, 4), border(all)
putpdf table riss_econ(1,1) = ("R-ISS Stage"), bold
putpdf table riss_econ(1,2) = ("N"), bold
putpdf table riss_econ(1,3) = ("Mean Cost (AUD)"), bold
putpdf table riss_econ(1,4) = ("Mean QALYs"), bold
putpdf table riss_econ(2,1) = ("Stage I")
putpdf table riss_econ(2,2) = ("`n_riss1'")
putpdf table riss_econ(2,3) = ("$`cost_riss1'")
putpdf table riss_econ(2,4) = ("`qaly_riss1'")
putpdf table riss_econ(3,1) = ("Stage II")
putpdf table riss_econ(3,2) = ("`n_riss2'")
putpdf table riss_econ(3,3) = ("$`cost_riss2'")
putpdf table riss_econ(3,4) = ("`qaly_riss2'")
putpdf table riss_econ(4,1) = ("Stage III")
putpdf table riss_econ(4,2) = ("`n_riss3'")
putpdf table riss_econ(4,3) = ("$`cost_riss3'")
putpdf table riss_econ(4,4) = ("`qaly_riss3'")

**********
* Save PDF
**********

set graphics on
local output_file "`report_dir'/${int}_${data}_${scenario}.pdf"
putpdf save "`output_file'", replace

di as result _n "✓ Report saved: `output_file'"
