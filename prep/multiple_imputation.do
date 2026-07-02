**********
* Monash Myeloma Model 
* Multiple Imputation
* 
* Purpose: Multiple Imputation of MRDR Long
*
* Author: Adam Irving
* Date: June 2026
**********

clear
clear mata

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $data_path (git-ignored)
	
local Data "$data_cut"
	
**********	
// Settings
/*
global imp 2
global boot 0
*/

global imp `1'
global boot `2'
global min_bs `3'
global max_bs `4'
global sample `5'   // "" = full cohort (main model); "train"/"test" = OOS fold (analyses/oos/)

* OOS routing: when $sample is set, restrict to that fold (split crosswalk written by
* analyses/oos/prep/oos_split.do) and write outputs under ${data_path}/oos/. Empty = main model.
if "$sample" == "" {
	global mi_outdir ""
	global mi_outtag ""
}
else {
	global mi_outdir "oos/"
	global mi_outtag "_$sample"
	capture mkdir "${data_path}/oos"
	capture mkdir "${data_path}/oos/bootstrap"
}

**********
// Define function
cap program drop multiple_imputation
program define multiple_imputation
	args RN1 RN2 RN3
	
	// Drop variables to be derived from imputed
	qui cap drop ISS
	qui cap drop RISS

	// MI settings
	mi set wide 
	mi register imputed Albumin AlkalinePhosphatase BMPlasmaCells LactateDehydrogenase SerumB2Microglobulin SerumCalcium SerumCreatinine eGFR EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore CM_CRD CM_PLM CM_DBT CM_MLG Male FISHRisk ExtraMedullaryD LyticLesion ECOGcc Para Lambda Kappa FLC dPara dLambda dKappa dFLC BCR
	mi register regular Age CLine
	mi describe
	
	// Diagnosis imputation	
	cap noi mi impute chained (regress) Albumin AlkalinePhosphatase BMPlasmaCells LactateDehydrogenase SerumB2Microglobulin SerumCalcium SerumCreatinine eGFR EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore Para Lambda Kappa FLC (logit, augment) Male CM_CRD CM_PLM CM_DBT CM_MLG FISHRisk ExtraMedullaryD LyticLesion (ologit, augment) ECOGcc = Age if Event0 == 3, add($imp) rseed(`RN1')
	if _rc {
		exit _rc
	}
		// ISS
		qui mi passive: gen ISS = 1 if SerumB2Microglobulin < 3.5 & Albumin >= 35 & Event0 == 3
		qui mi passive: replace ISS = 3 if SerumB2Microglobulin >= 5.5 & Event0 == 3
		qui mi passive: replace ISS = 2 if ISS == . & Event0 == 3
		label variable ISS "ISS (MI)"
			
		// LDHRisk
		qui mi passive: gen LDHRisk = 0 if LactateDehydrogenase <= LDHUpperLimit & Event0 == 3
		qui mi passive: replace LDHRisk = 1 if LactateDehydrogenase > LDHUpperLimit & Event0 == 3
		label variable LDHRisk "LDH Risk (MI)"
		
		// CKD
		qui mi passive: gen CM_CKD = 1 if eGFR <= 59 & Event0 == 3
		qui mi passive: replace CM_CKD = 0 if eGFR > 59 & Event0 == 3
		label variable CM_CKD "Chronic Kidney Disease"
		
		// CM Score
		qui mi passive: gen CM = CM_CKD + CM_CRD + CM_PLM + CM_DBT + CM_MLG
		qui mi passive: gen CMc = CM
		qui mi passive: replace CMc = 3 if CMc == 4 | CMc == 5
		
		// Carryforward
		local vars "Male ECOGcc ISS CMc CM_CKD CM_CRD CM_PLM CM_DBT CM_MLG"
		foreach v of local vars {
			qui mi xeq 0/$imp: bysort ID_BS (Date0): replace `v' = `v'[_n-1] if `v' == .
		}
		
		// RISS
		qui mi passive: gen RISS = 1 if ISS == 1 & LDHRisk == 0 & FISHRisk == 0 & Event0 == 3
		qui mi passive: replace RISS = 3 if ISS == 3 & (LDHRisk == 1 | FISHRisk == 1) & Event0 == 3
		qui mi passive: replace RISS = 2 if RISS == . & Event0 == 3
		qui mi xeq 0/$imp: bysort ID_BS (RISS): replace RISS = RISS[_n-1] if RISS == .
		label variable RISS "Revised ISS (MI)"

	// BCR imputation
	
		// TXR imputation
		cap noi mi impute chained (regress) dPara dLambda dKappa dFLC (ologit, augment) BCR = Age i.ECOGcc i.CLine if CStart == 1 & Duration != ., replace rseed(`RN2')
		if _rc {
			exit _rc
		}
		
		// ASCT imputation
		cap noi mi impute chained (regress) dPara dLambda dKappa dFLC (ologit, augment) BCR = Age i.ECOGcc if Event0 == 100, replace rseed(`RN3')
		if _rc {
			exit _rc
		}
		
		// Carryforward
		qui mi xeq 0/$imp: bysort ID_BS (Date0): replace BCR = BCR[_n-1] if BCR == . & Duration != .
		
		// Collapse BCR for ASCT - small n
		qui mi xeq 0/$imp: replace BCR = 4 if (BCR == 5 | BCR == 6) & Event0 == 100
			
		// Generate previous BCR (collapse ChemoS or SCT)
		preserve
		qui keep if CStart == 1 | Event0 == 100
		qui mi xeq 0/$imp: bysort ID: gen pBCR = BCR[_n-1]
		label values pBCR BCR_label
		qui keep ID_BS Event0 Date0 pBCR
		tempfile temp
		save `temp'	
		restore
			
		qui mmerge ID_BS Event0 Date0 using `temp', type(1:1) unmatched(master) nolabel
		qui mi update
		qui mi xeq 0/$imp: bysort ID: replace pBCR = pBCR[_n-1] if pBCR == . & Duration != .
		qui mi xeq 0/$imp: label values pBCR BCR_label
				
		// Generate BCR_L1 to BCR_L9
		forvalues l = 1/9 {
			qui mi passive: gen BCR_L`l' = BCR if Event0 == `l'0
			qui mi xeq 0/$imp: bysort ID_BS (BCR_L`l'): replace BCR_L`l' = BCR_L`l'[_n-1] if BCR_L`l' == .
			label values BCR_L`l' BCR_label
		}
				
		// Generate BCR_SCT
		qui mi passive: gen BCR_SCT = BCR if Event0 == 100
		qui mi xeq 0/$imp: bysort ID_BS (BCR_SCT): replace BCR_SCT = BCR_SCT[_n-1] if BCR_SCT == . 
		qui mi xeq 0/$imp: replace BCR_SCT = 0 if BCR_SCT == .  // Set to 0 for No SCT
		qui mi xeq 0/$imp: label values BCR_SCT BCR_label	
		
	// Unregister, keep, sort & order
	mi unregister AlkalinePhosphatase BMPlasmaCells SerumCalcium SerumCreatinine EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore ExtraMedullaryD LyticLesion Para Lambda Kappa FLC dPara dLambda dKappa dFLC
	keep ID ID_BS Event* Date* Age* Male ECOGcc ISS RISS SCT MNT CM* BCR* pBCR* Reg* OS Line Duration CID CLine CStart CEnd Country F_* CN_* Year Albumin SerumB2Microglobulin LactateDehydrogenase LDHUpperLimit LDHRisk FISHRisk eGFR _* Bortezomib Carfilzomib Cisplatin Cyclophosphamide Daratumamab Dexamethasone Doxorubicin Elotuzamab Etoposide Lenalidomide Melphalan Methylprednisolone Panobinostat Prednisolone Thalidomide Pomalidomide Ixazomib TXD* TFI*
	sort ID_BS Date0
	order $core Age Male ECOGcc RISS BCR Reg Regimen Line Duration

end

**********
// Execute based on arguments

if "$boot" == "0" {
	
	// Create temp folder on the shared drive (the locally-synced repo confuses Google Drive;
	// $data_path is mounted whenever this script runs).
	local repo = c(pwd)
	cap mkdir "${data_path}/temp"
	cd "${data_path}/temp"
		
	// Open MRDR Long Data
	use "${data_path}/MRDR Long.dta"
	gen ID_BS = ID

	// OOS: restrict to the requested fold (train/test) before imputing
	if "$sample" != "" {
		merge m:1 ID using "${data_path}/oos/oos_split.dta", keep(match) keepusing(fold) nogen
		keep if fold == "$sample"
		drop fold
	}

	// Draw random numbers
	local RN1 = 3949
	local RN2 = 6192
	local RN3 = 8273
		
	// Execute function
		multiple_imputation `RN1' `RN2' `RN3'
					
	// Save Long MI
	save "${data_path}/${mi_outdir}MRDR Long MI${mi_outtag}.dta", replace
			
	// Create Wide MI to create synthetic patient datasets
	keep if Event0 == 3
				
		// Convert _1_var to var_1
		foreach v in Male ECOGcc RISS ISS LDHRisk FISHRisk CMc CM_CKD CM_CRD CM_PLM CM_DBT CM_MLG {
			forvalues i = 1/$imp {
				gen `v'_`i' = _`i'_`v'
				drop _`i'_`v'
			}
		}
					
		// Drop non-MI variables
		keep ID Event0 Date0 Age Male* ECOGc* RISS* ISS* LDHRisk* FISHRisk* CMc* CM_* _mi_miss
		drop Male ECOGcc RISS ISS LDHRisk FISHRisk CMc CM_CKD CM_CRD CM_PLM CM_DBT CM_MLG
				
		// Turn MI imps into rows
		mi unset
		reshape long Male_ ECOGcc_ RISS_ ISS_ LDHRisk_ FISHRisk_ CMc_ CM_CKD_ CM_CRD_ CM_PLM_ CM_DBT_ CM_MLG_, i(ID) j(Imp)
		rename Male_ Male
		rename ECOGcc_ ECOGcc
		rename RISS_ RISS
		rename ISS_ ISS
		rename LDHRisk_ LDHRisk
		rename FISHRisk_ FISHRisk
		rename CMc_ CMc
		rename CM_CKD_ CM_CKD
		rename CM_CRD_ CM_CRD
		rename CM_PLM_ CM_PLM
		rename CM_DBT_ CM_DBT
		rename CM_MLG_ CM_MLG
		drop mi_miss
		
		// Save Wide MI
		gen MRDR = 1
		order MRDR ID Imp Event0 Date0 Age Male ECOGcc RISS ISS LDHRisk FISHRisk CMc CM_CKD CM_CRD CM_PLM CM_DBT CM_MLG
		save "${data_path}/${mi_outdir}MRDR Wide MI${mi_outtag}.dta", replace
			
	// Return to the repo root, then delete temp folder
	// (must cd out first, or Stata is left sitting in the just-removed directory)
	cd "`repo'"
	cap rmdir "${data_path}/temp"
}
else if "$boot" == "1" {
	forval b = $min_bs / $max_bs {
	
		// Open MRDR Long Data
		use "${data_path}/MRDR Long.dta"

		// OOS: restrict to the requested fold before resampling/imputing
		if "$sample" != "" {
			merge m:1 ID using "${data_path}/oos/oos_split.dta", keep(match) keepusing(fold) nogen
			keep if fold == "$sample"
			drop fold
		}

		// Create per-iteration temp folder on the shared drive - needed for M3 array jobs
		// (keeps scratch off the locally-synced repo; $data_path is mounted at run time)
		local curdir = c(pwd)
		capture mkdir "${data_path}/temp"
		capture mkdir "${data_path}/temp/job_`b'"
		cd "${data_path}/temp/job_`b'"
				
		// Base random numbers
		local RN1 = 7394
		local RN2 = 1392
		local RN3 = 4771

		// Retry settings
		local maxtries = 50
		local success = 0
		local try = 0

		// Snapshot the clean (pre-resample) data for restarts
		preserve

		while `success' == 0 & `try' < `maxtries' {
			local ++try

			// Restore clean data before each attempt
			restore, preserve

			// Attempt-specific, deterministic seeds
			local RS  = 5839 * `b' + 100003 * `try'
			local a1  = `RN1' + 911 * `try'
			local a2  = `RN2' + 911 * `try'
			local a3  = `RN3' + 911 * `try'

			// Resample with this attempt's seed
			set seed `RS'
			bsample, cluster(ID) idcluster(ID_BS)

			// Try imputation; trap perfect-predictor / convergence failures
			cap noi multiple_imputation `a1' `a2' `a3'

			if _rc == 0 {
				local success = 1
			}
			else {
				di as txt "Bootstrap `b': attempt `try' failed (rc=" _rc "), retrying with new seeds"
			}
		}

		// Drop the snapshot
		restore, not

		if `success' == 0 {
			di as error "Bootstrap `b' failed after `maxtries' attempts"
			exit 459
		}

		// Save Long MI
		cd "`curdir'"
		save "${data_path}/${mi_outdir}bootstrap/MRDR Long MI${mi_outtag} B`b'.dta", replace
			
		// Delete temporary folder
		cap rmdir "${data_path}/temp/job_`b'"
	}

	// Remove the now-empty parent temp dir
	cap rmdir "${data_path}/temp"
}
