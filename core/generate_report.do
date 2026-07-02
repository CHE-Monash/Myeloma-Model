**********
* Monash Myeloma Model - Generate Report
* 
* Purpose: Create PDF report using putpdf
*          Integrates with Monash Myeloma Model v2.0 dispatcher system
* 
* Output: ${int}_${data}_${scenario}.pdf
**********


**********
* Paths — simulated output lives under an optional scenario subfolder (empty
* $scenario => simulated_path itself), matching simulate.do / export_results.do.
local sim_out = cond("$scenario" == "", "$simulated_path", "$simulated_path/$scenario")
capture mkdir "`sim_out'"
capture mkdir "`sim_out'/report"
local report_dir "`sim_out'/report"

* Load Data
qui use "`sim_out'/${int}_${line}_${data}.dta", clear


**********
* Start PDF
**********

capture putpdf clear
set graphics off
putpdf begin

// Title page
putpdf paragraph, halign(center)
putpdf text ("Monash Myeloma Model v3.0"), bold font(,18)

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

// Line 1 Regimen × BCR Cross-tabulation
quietly count if TXR_L1 < .
local l1_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL1_TXR)
mata: st_numscalar("n_regimens", cols(oL1_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
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
    
    // Count N for this regimen
    quietly count if TXR_L1 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l1_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD)
    forval b = 1/6 {
        quietly count if TXR_L1 == `code' & BCR_L1 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 8 rows (header + N + 6 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l1 = (8, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l1(1,1) = ("Line 1"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l1(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l1(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l1(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR" "VGPR" "PR" "MR" "SD" "PD" "'
forval b = 1/6 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l1(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l1(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

// ASCT BCR Table with Patient Characteristics
quietly count if TXR_L1 < .
local l1_total = r(N)

quietly count if SCT_L1 == 1
local sct_n = r(N)
local sct_n_fmt = string(r(N), "%9.0fc")
local sct_pct = string(100*r(N)/`l1_total', "%4.1f")

quietly count if SCT_L1 == 0
local nosct_n = r(N)
local nosct_n_fmt = string(r(N), "%9.0fc")
local nosct_pct = string(100*r(N)/`l1_total', "%4.1f")

// BCR within ASCT patients (BCR_SCT: 1=CR, 2=VGPR, 3=PR, 4=MR)
forval b = 1/4 {
    quietly count if SCT_L1 == 1 & BCR_SCT == `b'
    if `sct_n' > 0 {
        local sct_bcr`b' = string(100*r(N)/`sct_n', "%4.1f")
    }
    else {
        local sct_bcr`b' = "—"
    }
}

// Age breakdown for ASCT patients
quietly count if SCT_L1 == 1 & Age_L1S < 65
local sct_age1 = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 1 & Age_L1S >= 65 & Age_L1S < 70
local sct_age2 = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 1 & Age_L1S >= 70 & Age_L1S < 75
local sct_age3 = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 1 & Age_L1S >= 75
local sct_age4 = string(100*r(N)/`sct_n', "%4.1f")

// Age breakdown for No ASCT patients
quietly count if SCT_L1 == 0 & Age_L1S < 65
local nosct_age1 = string(100*r(N)/`nosct_n', "%4.1f")
quietly count if SCT_L1 == 0 & Age_L1S >= 65 & Age_L1S < 70
local nosct_age2 = string(100*r(N)/`nosct_n', "%4.1f")
quietly count if SCT_L1 == 0 & Age_L1S >= 70 & Age_L1S < 75
local nosct_age3 = string(100*r(N)/`nosct_n', "%4.1f")
quietly count if SCT_L1 == 0 & Age_L1S >= 75
local nosct_age4 = string(100*r(N)/`nosct_n', "%4.1f")

// BCR L1 CR/VGPR for ASCT vs No ASCT
quietly count if SCT_L1 == 1 & (BCR_L1 == 1 | BCR_L1 == 2)
local sct_crvgpr = string(100*r(N)/`sct_n', "%4.1f")
quietly count if SCT_L1 == 0 & (BCR_L1 == 1 | BCR_L1 == 2)
local nosct_crvgpr = string(100*r(N)/`nosct_n', "%4.1f")

// Create table: 6 rows × 5 columns
putpdf table txr_bcr_sct = (6, 5), border(all)

// Header row
putpdf table txr_bcr_sct(1,1) = ("ASCT"), bold
putpdf table txr_bcr_sct(1,2) = (""), bold
putpdf table txr_bcr_sct(1,3) = (""), bold
putpdf table txr_bcr_sct(1,4) = ("ASCT"), bold
putpdf table txr_bcr_sct(1,5) = ("No ASCT"), bold

// Row 2: N / Age < 65
putpdf table txr_bcr_sct(2,1) = ("N"), bold
putpdf table txr_bcr_sct(2,2) = ("`sct_n_fmt' (`sct_pct'%)")
putpdf table txr_bcr_sct(2,3) = ("Age < 65")
putpdf table txr_bcr_sct(2,4) = ("`sct_age1'%")
putpdf table txr_bcr_sct(2,5) = ("`nosct_age1'%")

// Row 3: CR / Age 65-69
putpdf table txr_bcr_sct(3,1) = ("CR")
putpdf table txr_bcr_sct(3,2) = ("`sct_bcr1'%")
putpdf table txr_bcr_sct(3,3) = ("Age >= 65 & < 70")
putpdf table txr_bcr_sct(3,4) = ("`sct_age2'%")
putpdf table txr_bcr_sct(3,5) = ("`nosct_age2'%")

// Row 4: VGPR / Age 70-74
putpdf table txr_bcr_sct(4,1) = ("VGPR")
putpdf table txr_bcr_sct(4,2) = ("`sct_bcr2'%")
putpdf table txr_bcr_sct(4,3) = ("Age >= 70 & < 75")
putpdf table txr_bcr_sct(4,4) = ("`sct_age3'%")
putpdf table txr_bcr_sct(4,5) = ("`nosct_age3'%")

// Row 5: PR / Age 75+
putpdf table txr_bcr_sct(5,1) = ("PR")
putpdf table txr_bcr_sct(5,2) = ("`sct_bcr3'%")
putpdf table txr_bcr_sct(5,3) = ("Age >= 75")
putpdf table txr_bcr_sct(5,4) = ("`sct_age4'%")
putpdf table txr_bcr_sct(5,5) = ("`nosct_age4'%")

// Row 6: MR / BCR L1 CR/VGPR
putpdf table txr_bcr_sct(6,1) = ("MR")
putpdf table txr_bcr_sct(6,2) = ("`sct_bcr4'%")
putpdf table txr_bcr_sct(6,3) = ("BCR L1 CR/VGPR")
putpdf table txr_bcr_sct(6,4) = ("`sct_crvgpr'%")
putpdf table txr_bcr_sct(6,5) = ("`nosct_crvgpr'%")

// Line 2 Regimen × BCR Cross-tabulation
quietly count if TXR_L2 < .
local l2_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL2_TXR)
mata: st_numscalar("n_regimens", cols(oL2_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
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
    
    // Count N for this regimen
    quietly count if TXR_L2 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l2_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD)
    forval b = 1/6 {
        quietly count if TXR_L2 == `code' & BCR_L2 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 8 rows (header + N + 6 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l2 = (8, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l2(1,1) = ("Line 2"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l2(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l2(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l2(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR" "VGPR" "PR" "MR" "SD" "PD" "'
forval b = 1/6 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l2(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l2(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

// Line 3 Regimen × BCR Cross-tabulation
quietly count if TXR_L3 < .
local l3_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL3_TXR)
mata: st_numscalar("n_regimens", cols(oL3_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
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
    
    // Count N for this regimen
    quietly count if TXR_L3 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l3_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR/VGPR, 3=PR/MR, 5=SD/PD)
    forval b = 1(2)5 {
        quietly count if TXR_L3 == `code' & BCR_L3 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 5 rows (header + N + 3 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l3 = (5, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l3(1,1) = ("Line 3"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l3(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l3(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l3(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR/VGPR" "PR/MR" "SD/PD"'
forval b = 1/3 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l3(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l3(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

// Line 4 Regimen × BCR Cross-tabulation
quietly count if TXR_L4 < .
local l4_total = r(N)

// Get regimen codes from Mata
mata: st_matrix("regimen_codes", oL4_TXR)
mata: st_numscalar("n_regimens", cols(oL4_TXR))
local n_regimens = n_regimens 

// Extract codes from matrix to locals and build name mapping
forval i = 1/`n_regimens' {
    local code_`i' = regimen_codes[1,`i']
    local code = `code_`i''
    
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
    
    // Count N for this regimen
    quietly count if TXR_L4 == `code'
    local reg_`i'_n = r(N)
    local reg_`i'_n_fmt = string(r(N), "%9.0fc")
    local reg_`i'_pct = string(100*r(N)/`l4_total', "%4.1f")
    
    // Count BCR within this regimen (BCR 1=CR/VGPR, 3=PR/MR, 5=SD/PD)
    forval b = 1(2)5 {
        quietly count if TXR_L4 == `code' & BCR_L4 == `b'
        if `reg_`i'_n' > 0 {
            local reg_`i'_bcr`b' = string(100*r(N)/`reg_`i'_n', "%4.1f")
        }
        else {
            local reg_`i'_bcr`b' = "—"
        }
    }
}

// Create table: 5 rows (header + N + 3 BCR) × (1 + n_regimens) columns
local n_cols = `n_regimens' + 1
putpdf table txr_bcr_l4 = (5, `n_cols'), border(all)

// Header row
putpdf table txr_bcr_l4(1,1) = ("Line 4"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l4(1,`col') = ("`reg_`i'_name'"), bold
}

// N row
putpdf table txr_bcr_l4(2,1) = ("N"), bold
forval i = 1/`n_regimens' {
    local col = `i' + 1
    putpdf table txr_bcr_l4(2,`col') = ("`reg_`i'_n_fmt' (`reg_`i'_pct'%)")
}

// BCR rows
local bcr_names `" "CR/VGPR" "PR/MR" "SD/PD"'
forval b = 1/3 {
    local row = `b' + 2
    local bcr_label : word `b' of `bcr_names'
    putpdf table txr_bcr_l4(`row',1) = ("`bcr_label'")
    
    forval i = 1/`n_regimens' {
        local col = `i' + 1
        putpdf table txr_bcr_l4(`row',`col') = ("`reg_`i'_bcr`b''%")
    }
}

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
    graphregion(color(white)) name(os, replace)
graph export "$simulated_path/report/figures/os.png", replace width(1200)

putpdf paragraph
putpdf text ("Overall Survival"), bold linebreak(2)
putpdf image "$simulated_path/report/figures/os.png", width(7)

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
putpdf pagebreak

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
putpdf pagebreak

preserve
stset OC_TIME if OC_TIME < 240, failure(OC_MORT)
    
sts graph, by(RISS) ///
	xtitle("Months") ytitle("Probability") title("") ///
	ylabel(0(0.2)1, angle(0)) xlabel(0(24)240) ///
	graphregion(color(white)) ///
	legend(label(1 "Stage I") label(2 "Stage II") label(3 "Stage III") rows(1)) ///
	name(os_riss, replace)
graph export "$simulated_path/report/figures/os_riss.png", replace width(1200)

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

**********
* Lines of Therapy
**********

putpdf pagebreak

putpdf paragraph
putpdf text ("Lines of Therapy Distribution"), bold font(,16)

// Count patients receiving each line using TXR_L* variables
// If TXR_L# is not missing, patient received that line

// Calculate N receiving each line
local n_total = _N
forvalues l = 1/9 {
	quietly count if !missing(TXR_L`l')
	local n_l`l' = r(N)
	local pct_l`l' = (`n_l`l'' / `n_total') * 100
}

// Determine max line reached per patient
gen LOT_MAX = 0
forvalues l = 1/9 {
	replace LOT_MAX = `l' if !missing(TXR_L`l')
}

// Count patients by their maximum line reached
forvalues l = 1/9 {
	quietly count if LOT_MAX == `l'
	local n_max_l`l' = r(N)
	local pct_max_l`l' = (`n_max_l`l'' / `n_total') * 100
}

// Display to log
di _col(5) "Line" _col(15) "Received" _col(30) "%" _col(40) "Max Line" _col(55) "%"
di "{hline 60}"
forvalues l = 1/9 {
	di _col(5) "L`l'" _col(15) %9.0fc `n_l`l'' _col(30) %5.1f `pct_l`l'' _col(40) %9.0fc `n_max_l`l'' _col(55) %5.1f `pct_max_l`l''
}

// Create PDF table
putpdf paragraph

putpdf table lot_tbl = (11, 5), border(all)
putpdf table lot_tbl(1,1) = ("Line"), bold
putpdf table lot_tbl(1,2) = ("N Received"), bold
putpdf table lot_tbl(1,3) = ("% Received"), bold
putpdf table lot_tbl(1,4) = ("N Max Line"), bold
putpdf table lot_tbl(1,5) = ("% Max Line"), bold

forvalues l = 1/9 {
	local row = `l' + 1
	putpdf table lot_tbl(`row',1) = ("L`l'")
	putpdf table lot_tbl(`row',2) = ("`=string(`n_l`l'', "%9.0fc")'")
	putpdf table lot_tbl(`row',3) = ("`=string(`pct_l`l'', "%5.1f")'%")
	putpdf table lot_tbl(`row',4) = ("`=string(`n_max_l`l'', "%9.0fc")'")
	putpdf table lot_tbl(`row',5) = ("`=string(`pct_max_l`l'', "%5.1f")'%")
}

putpdf table lot_tbl(11,1) = ("Total"), bold
putpdf table lot_tbl(11,2) = ("`=string(`n_total', "%9.0fc")'"), bold
putpdf table lot_tbl(11,3) = ("—")
putpdf table lot_tbl(11,4) = ("`=string(`n_total', "%9.0fc")'"), bold
putpdf table lot_tbl(11,5) = ("100.0%"), bold

// Summary statistics
quietly summarize LOT_MAX, detail
local mean_lot = string(r(mean), "%4.2f")
local median_lot = string(r(p50), "%4.0f")

putpdf paragraph
putpdf text ("Mean lines received: `mean_lot'; Median: `median_lot'")

// Mortality by maximum line reached

putpdf paragraph
putpdf text ("Mortality by Maximum Line Reached"), bold

putpdf table mort_lot = (11, 4), border(all)
putpdf table mort_lot(1,1) = ("Line"), bold
putpdf table mort_lot(1,2) = ("Deaths"), bold
putpdf table mort_lot(1,3) = ("Deaths during TXD"), bold
putpdf table mort_lot(1,4) = ("Deaths during TFI"), bold

// Row for DN (patients who died before L1)
quietly count if MOR_DN == 1
local n_died_dn = r(N)
putpdf table mort_lot(2,1) = ("DN")
putpdf table mort_lot(2,2) = ("`=string(`n_died_dn', "%9.0fc")'")
putpdf table mort_lot(2,3) = ("—")
putpdf table mort_lot(2,4) = ("`=string(`n_died_dn', "%9.0fc")'")

// L1 to L9
forvalues l = 1/9 {
	local row = `l' + 2
	
	// N with max line = L`l'
	quietly count if LOT_MAX == `l'
	local n_lot = r(N)
	
	// Deaths on treatment (MOR_L#S)
	quietly count if LOT_MAX == `l' & MOR_L`l'S == 1
	local n_died_s = r(N)
	
	// Deaths at exit/TFI (MOR_L#E)
	quietly count if LOT_MAX == `l' & MOR_L`l'E == 1
	local n_died_e = r(N)
	
	putpdf table mort_lot(`row',1) = ("L`l'")
	putpdf table mort_lot(`row',2) = ("`=string(`n_lot', "%9.0fc")'")
	putpdf table mort_lot(`row',3) = ("`=string(`n_died_s', "%9.0fc")'")
	putpdf table mort_lot(`row',4) = ("`=string(`n_died_e', "%9.0fc")'")
}

**********
* Treatment Duration
**********

putpdf pagebreak
putpdf paragraph
putpdf text ("Treatment Duration (TXD)"), bold font(,16)

// Summary statistics for TXD by line
putpdf paragraph
putpdf text ("TXD Summary Statistics (months)"), bold

// Create summary table for TXD
local nrows = 10
putpdf table txd_sum = (`nrows', 6), border(all)
putpdf table txd_sum(1,1) = ("Line"), bold
putpdf table txd_sum(1,2) = ("N"), bold
putpdf table txd_sum(1,3) = ("Mean"), bold
putpdf table txd_sum(1,4) = ("SD"), bold
putpdf table txd_sum(1,5) = ("Median"), bold
putpdf table txd_sum(1,6) = ("IQR"), bold

forvalues l = 1/9 {
	local row = `l' + 1
	quietly summarize TXD_L`l' if TXD_L`l' > 0, detail
	local n = r(N)
	local mean = string(r(mean), "%5.1f")
	local sd = string(r(sd), "%5.1f")
	local median = string(r(p50), "%5.1f")
	local iqr = "[" + string(r(p25), "%4.1f") + "-" + string(r(p75), "%4.1f") + "]"
		
	putpdf table txd_sum(`row',1) = ("L`l'")
	putpdf table txd_sum(`row',2) = ("`=string(`n', "%9.0fc")'")
	putpdf table txd_sum(`row',3) = ("`mean'")
	putpdf table txd_sum(`row',4) = ("`sd'")
	putpdf table txd_sum(`row',5) = ("`median'")
	putpdf table txd_sum(`row',6) = ("`iqr'")
		
	di "TXD L`l': N=" %9.0fc `n' ", Mean=" %5.1f r(mean) ", Median=" %5.1f r(p50)
}

// Generate KM curves for TXD (Lines 1-9)
preserve

gen patient_id = _n
tempfile base_data
save `base_data'

// Stack all lines into long format
clear
local first = 1
forvalues l = 1/9 {
	use `base_data', clear
	keep patient_id TXD_L`l'
	rename TXD_L`l' txd_time
	gen line = `l'
	keep if !missing(txd_time) & txd_time > 0
	gen txd_event = 1
			
	if `first' {
		tempfile txd_stacked
		save `txd_stacked'
		local first = 0
	}
	else {
		append using `txd_stacked'
	save `txd_stacked', replace
	}
}

label define line_lbl 1 "L1" 2 "L2" 3 "L3" 4 "L4" 5 "L5" 6 "L6" 7 "L7" 8 "L8" 9 "L9"
label values line line_lbl

stset txd_time, failure(txd_event)

sts graph, by(line) ///
	ytitle("Proportion on treatment") xtitle("Months") title("") ///
	ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
	legend(order(1 "L1" 2 "L2" 3 "L3" 4 "L4" 5 "L5" ///
           6 "L6" 7 "L7" 8 "L8" 9 "L9") ///
			rows(2) size(small) pos(6)) ///
	scheme(s2color) ///
	graphregion(color(white)) ///
	name(txd, replace)

graph export "$simulated_path/report/figures/txd.png", replace width(1600)

restore

putpdf paragraph
putpdf text ("Treatment Duration by Line of Therapy"), bold
putpdf image "$simulated_path/report/figures/txd.png", width(8)


**********
* Treatment-free Intervals
**********

putpdf pagebreak
putpdf paragraph
putpdf text ("Treatment-Free Interval (TFI)"), bold font(,16)

// Summary statistics for TFI
putpdf paragraph
putpdf text ("TFI Summary Statistics (months)"), bold

// Create summary table for TFI (DN→L1 through L8→L9)
local nrows = 10  // Header + DN + L1-L8
putpdf table tfi_sum = (`nrows', 7), border(all)
putpdf table tfi_sum(1,1) = ("Interval"), bold
putpdf table tfi_sum(1,2) = ("N"), bold
putpdf table tfi_sum(1,3) = ("Mean"), bold
putpdf table tfi_sum(1,4) = ("SD"), bold
putpdf table tfi_sum(1,5) = ("Median"), bold
putpdf table tfi_sum(1,6) = ("IQR"), bold
putpdf table tfi_sum(1,7) = ("Zero TFI %"), bold

// Row 2: DN→L1
capture confirm variable TFI_DN
if !_rc {
	quietly summarize TFI_DN if !missing(TFI_DN), detail
	local n = r(N)
	local mean = string(r(mean), "%5.1f")
	local sd = string(r(sd), "%5.1f")
	local median = string(r(p50), "%5.1f")
	local iqr = "[" + string(r(p25), "%4.1f") + "-" + string(r(p75), "%4.1f") + "]"
	
	quietly count if TFI_DN == 0 & !missing(TFI_DN)
	local n_zero = r(N)
	local pct_zero = string((`n_zero' / `n') * 100, "%4.1f")
	
	putpdf table tfi_sum(2,1) = ("DN→L1")
	putpdf table tfi_sum(2,2) = ("`=string(`n', "%9.0fc")'")
	putpdf table tfi_sum(2,3) = ("`mean'")
	putpdf table tfi_sum(2,4) = ("`sd'")
	putpdf table tfi_sum(2,5) = ("`median'")
	putpdf table tfi_sum(2,6) = ("`iqr'")
	putpdf table tfi_sum(2,7) = ("`pct_zero'%")
	
	di "TFI DN→L1: N=" %9.0fc `n' ", Mean=" `mean' ", Median=" `median' ", Zero TFI=" `pct_zero' "%"
}

// Rows 3-10: L1→L2 through L8→L9
forvalues l = 1/8 {
	local row = `l' + 2
	local next_line = `l' + 1
	capture confirm variable TFI_L`l'
	if !_rc {
		quietly summarize TFI_L`l' if !missing(TFI_L`l'), detail
		local n = r(N)
		local mean = string(r(mean), "%5.1f")
		local sd = string(r(sd), "%5.1f")
		local median = string(r(p50), "%5.1f")
		local iqr = "[" + string(r(p25), "%4.1f") + "-" + string(r(p75), "%4.1f") + "]"
		
		quietly count if TFI_L`l' == 0 & !missing(TFI_L`l')
		local n_zero = r(N)
		local pct_zero = string((`n_zero' / `n') * 100, "%4.1f")
		
		putpdf table tfi_sum(`row',1) = ("L`l'→L`next_line'")
		putpdf table tfi_sum(`row',2) = ("`=string(`n', "%9.0fc")'")
		putpdf table tfi_sum(`row',3) = ("`mean'")
		putpdf table tfi_sum(`row',4) = ("`sd'")
		putpdf table tfi_sum(`row',5) = ("`median'")
		putpdf table tfi_sum(`row',6) = ("`iqr'")
		putpdf table tfi_sum(`row',7) = ("`pct_zero'%")
		
		di "TFI L`l'→L`next_line': N=" %9.0fc `n' ", Mean=" `mean' ", Median=" `median' ", Zero TFI=" `pct_zero' "%"
	}
}
	
// Generate KM curves for TFI (DN→L1 through L8→L9)
preserve

gen patient_id = _n
tempfile base_data
save `base_data'

// Stack all TFI intervals into long format
clear
local first = 1

// TFI DN
use `base_data', clear
keep patient_id TFI_DN
rename TFI_DN tfi_time
gen interval = 0
keep if !missing(tfi_time)
gen tfi_event = 1
			
tempfile tfi_stacked
save `tfi_stacked'
local first = 0

// TFI L1 to L8
forvalues l = 1/8 {
use `base_data', clear
	keep patient_id TFI_L`l'
	rename TFI_L`l' tfi_time
	gen interval = `l'
	keep if !missing(tfi_time)
	gen tfi_event = 1
				
	if `first' {
		tempfile tfi_stacked
		save `tfi_stacked'
		local first = 0
	}
	else {
	append using `tfi_stacked'
	save `tfi_stacked', replace
	}
}

label define int_lbl 0 "DN" 1 "L1" 2 "L2" 3 "L3" 4 "L4" ///
						 5 "L5" 6 "L6" 7 "L7" 8 "L8"
label values interval int_lbl

stset tfi_time, failure(tfi_event)

sts graph, by(interval) ///
	ytitle("Proportion not yet starting next line") xtitle("Months") title("") ///
	ylabel(0(0.2)1, angle(0) format(%3.1f)) ///
	legend(order(1 "DN" 2 "L1" 3 "L2" 4 "L3" 5 "L4" ///
           6 "L5" 7 "L6" 8 "L7" 9 "L8") ///
			rows(2) size(small) pos(6)) ///
	scheme(s2color) ///
	graphregion(color(white)) ///
	name(tfi, replace)

graph export "$simulated_path/report/figures/tfi.png", replace width(1600)

restore
	
putpdf paragraph
putpdf text ("Treatment-Free Interval by LoT"), bold
putpdf image "$simulated_path/report/figures/tfi.png", width(8)

**********
* Economic Outcomes
**********

putpdf pagebreak
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
quietly summarize cost_total_d, detail
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

quietly summarize cost_tx_d
local mean_tx = string(r(mean), "%12.0fc")
quietly summarize cost_nt_d
local mean_nt = string(r(mean), "%12.0fc")

// ASCT costs (recipients only)
quietly summarize cost_tx_asct_d if SCT_L1 == 1
if r(N) > 0 {
	local mean_asct = string(r(mean), "%12.0fc")
	local n_asct = string(r(N), "%9.0fc")
}
else {
	local mean_asct = "N/A"
	local n_asct = "0"
}

// Maintenance costs (recipients only)
quietly summarize cost_tx_mnt_d if MNT == 1
if r(N) > 0 {
	local mean_mnt = string(r(mean), "%12.0fc")
	local n_mnt = string(r(N), "%9.0fc")
}
else {
	local mean_mnt = "N/A"
	local n_mnt = "0"
}

// Non-treatment costs: process_data.do applies a single blended per-month rate
// (hospitalisation + community + emergency) and stores only the combined
// cost_nt_d, so report it as one line rather than an artificial component split.

putpdf table cost_comp = (5, 2), border(all)
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

// Treatment Costs by Line of Therapy
putpdf paragraph
putpdf text ("Treatment Costs by Line of Therapy"), bold

// Count rows needed (only include lines with patients)
local n_lines = 0
forval l = 1/9 {
	quietly count if cost_tx_L`l'_d != . & cost_tx_L`l'_d > 0
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
	quietly count if cost_tx_L`l'_d != . & cost_tx_L`l'_d > 0
	local n_l = string(r(N), "%9.0fc")
	quietly summarize cost_tx_L`l'_d if cost_tx_L`l'_d > 0
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
quietly summarize qaly_total_d, detail
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

quietly summarize qaly_tfi_DN_d
local q_tfi_dn = string(r(mean), "%5.3f")
quietly summarize qaly_txd_L1_d
local q_txd_l1 = string(r(mean), "%5.3f")
quietly summarize qaly_tfi_L1_d
local q_tfi_l1 = string(r(mean), "%5.3f")
quietly summarize qaly_txd_L2_d
local q_txd_l2 = string(r(mean), "%5.3f")
quietly summarize qaly_post_L2_d
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
quietly summarize cost_total
local mean_cost_undisc = string(r(mean), "%12.0fc")
quietly summarize cost_total_d
local mean_cost_disc = string(r(mean), "%12.0fc")

quietly summarize qaly_total
local mean_qaly_undisc = string(r(mean), "%5.2f")
quietly summarize qaly_total_d
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

quietly summarize cost_total_d if SCT_L1 == 0
local cost_noasct = string(r(mean), "%12.0fc")
quietly summarize qaly_total_d if SCT_L1 == 0
local qaly_noasct = string(r(mean), "%5.2f")
quietly count if SCT_L1 == 0
local n_noasct = string(r(N), "%9.0fc")

quietly summarize cost_total_d if SCT_L1 == 1
local cost_asct = string(r(mean), "%12.0fc")
quietly summarize qaly_total_d if SCT_L1 == 1
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
	quietly summarize cost_total_d if `age_cond'
	local cost_age`a' = string(r(mean), "%12.0fc")
	quietly summarize qaly_total_d if `age_cond'
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
	quietly summarize cost_total_d if RISS == `r'
	local cost_riss`r' = string(r(mean), "%12.0fc")
	quietly summarize qaly_total_d if RISS == `r'
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
local output_file "`report_dir'/${int}_${line}_${data}.pdf"
putpdf save "`output_file'", replace
n di as result _n "Report saved at `output_file'."
