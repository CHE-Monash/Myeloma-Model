**********
* Monash Myeloma Model
* Risk Equations
* 
* Purpose: Estimate risk equations for use in simulation
*
* Author: Adam Irving
* Date: December 2025
**********

clear
clear mata

local Data "$data_cut"
	
**********
// Capture arguments
global analysis `1'
global coeffs `2'
global min_year `3'
global max_year `4'
global boot `5'
global min_bs `6'
global max_bs `7'
global sample `8'   // "" = standard MI data; "train"/"test" = OOS fold under ${data_path}/oos/

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $data_path (git-ignored)

* OOS routing: read the OOS-fold imputed data when $sample is set (output is already keyed by
* $analysis, so OOS coefficients land in analyses/oos/coefficients/). Empty = main model.
if "$sample" == "" {
	global re_indir ""
	global re_intag ""
}
else {
	global re_indir "oos/"
	global re_intag "_$sample"
}
		
**********
// Define functions
cap program drop save_coefs
program define save_coefs
	args mat
	
	mata: b`mat' = st_matrix("e(b_mi)")
	mata: rb`mat' = st_matrixrowstripe("e(b_mi)")
	mata: cb`mat' = st_matrixcolstripe("e(b_mi)")
	
	scalar ncol = colsof(e(b_mi))
	forval i = 1/`=ncol'{
		mata: cb`mat'[`i',1] = "`i'"
	}
	scalar drop ncol
	global Coeffs $Coeffs b`mat' rb`mat' cb`mat'
	
	// For parametric survival analysis - save functional form
	if ("e(ecmd_mi)" != "") { 
		mata: fb`mat' = st_global("e(ecmd_mi)")
		global Coeffs $Coeffs fb`mat'
	}
	
	// For multinomial/ordered logt: save outcome levels
	capture mata: o`mat' = st_matrix("e(out)")
	if _rc == 0 {
		global Coeffs $Coeffs o`mat'
	}

	ereturn clear
end

cap program drop save_max
program define save_max
	args mat
	
	mata: max`mat' = st_numscalar("r(max)")
	global Coeffs $Coeffs max`mat'
end

cap program drop risk_equations
program define risk_equations
	
	// OS distribution
	global dOS "w"
		
	// TXD distribution
	global dTXD "w"
		
	// TFI distribution
	global dTFI = "w"
			
	// Pick up list of regimens by line
	qui do "analyses/$analysis/outcomes/txr_$coeffs.do"
		
	// Reset global
	global Coeffs

	// Risk equations
	
	di "Overall Survival"

	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 3) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS b0.OS#b5.BCR, d($dOS)
	save_coefs OS

	mata: _matrix_list(bOS, rbOS, cbOS)

	di "Diagnosis"
	
		// ASCT
	
		mi estimate: logit SCT Age Age2 Male i.ECOGcc i.RISS Age70 Age75 i.CMc if(Event0 == 3)
		save_coefs DN_SCT
		
		mata: _matrix_list(bDN_SCT, rbDN_SCT, cbDN_SCT)
			
		// Treatment-free Interval
		
		mi stset Date1, id(ID_BS) failure(Event1 == 10) origin(Event1 == 3) scale(30.4375)
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS SCT, d($dTFI)
		save_coefs DN_TFI
		
		mata: _matrix_list(bDN_TFI, rbDN_TFI, cbDN_TFI)
		
	di "Line 1" 
	
		// Treatment Regimen
	
		qui tab TXR_L1 // Check if modelling specific regimens
		if `r(r)' > 1 {
			mi estimate: mlogit TXR_L1 Age Age2 Male i.ECOGcc i.RISS SCT if(Event0 == 10 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
			*ASCT covariate, could separate by ASCT in the future
			save_coefs L1_TXR	
			
			mata: _matrix_list(bL1_TXR, rbL1_TXR, cbL1_TXR)	
		}
		else { // If all patients receive 'other'
			mata: oL1_TXR = 0
			global Coeffs $Coeffs oL1_TXR
		}

		// L1 Best Clinical Response
		
		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS SCT i.TXR_L1 if(Event0 == 10)
		save_coefs L1_BCR	
		
		mata: _matrix_list(bL1_BCR, rbL1_BCR, cbL1_BCR)

		// Treatment Duration

			// Fixed w/ ASCT - L1_TXD_ASCT
				
				// Cut off 1 - 60 for optimal fit with 3 splines, 90 for 2 splines
				scalar L1_TXD_ASCT_C1 = 60
				mata: L1_TXD_ASCT_C1 = st_numscalar("L1_TXD_ASCT_C1")
							
				// Cut off 2 - 120 for optimal fit with 3 splines
				scalar L1_TXD_ASCT_C2 = 120
				mata: L1_TXD_ASCT_C2 = st_numscalar("L1_TXD_ASCT_C2")
						
				global Coeffs $Coeffs L1_TXD_ASCT_C1 L1_TXD_ASCT_C2 
				
				// Spline 1
				gen S1_exit = Date1 + `=L1_TXD_ASCT_C1'
				mi stset Date1 if(F_L1 != 1 & CN_L1 != 1 & Duration < 730 & TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) exit(S1_exit) scale(30.4375)
				
				mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD) 
				save_coefs L1_TXD_ASCT_S1	
				
				mata: _matrix_list(bL1_TXD_ASCT_S1, rbL1_TXD_ASCT_S1, cbL1_TXD_ASCT_S1)
				
				// Spline 2
				gen S2_enter = S1_exit + 1
				gen S2_exit = Date1 + `=L1_TXD_ASCT_C2'
					
				mi stset Date1 if(F_L1 != 1 & CN_L1 != 1 & Duration < 730 & TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) enter(S2_enter) exit(S2_exit) scale(30.4375)
				
				mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
				save_coefs L1_TXD_ASCT_S2		
				
				mata: _matrix_list(bL1_TXD_ASCT_S2, rbL1_TXD_ASCT_S2, cbL1_TXD_ASCT_S2)
				
				// Spline 3
				gen S3_enter = S2_exit + 1
				mi stset Date1 if(F_L1 != 1 & CN_L1 != 1 & Duration < 730 & TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) enter(S3_enter) scale(30.4375)
				
				mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
				save_coefs L1_TXD_ASCT_S3	
				
				drop S1_exit S2_enter S2_exit S3_enter
				mata: _matrix_list(bL1_TXD_ASCT_S3, rbL1_TXD_ASCT_S3, cbL1_TXD_ASCT_S3)
						
			// Fixed w/o SCT  - L1_TXD_NoASCT 
			
			mi stset Date1 if(F_L1 != 1 & CN_L1 != 1 & Duration < 730 & TXR_L1 != 7 & SCT == 0), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) scale(30.4375)
			
			mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR i.TXR_L1, d($dTXD)
			save_coefs L1_TXD_NoASCT
			
			mata: _matrix_list(bL1_TXD_NoASCT, rbL1_TXD_NoASCT, cbL1_TXD_NoASCT)
				
			// Continuous - L1_TXD_Cont
			
			count if(TXR_L1 == 7)
			if r(N) > 0 { // If Rd is included
				mi stset Date1 if(F_L1 != 1 & CN_L1 != 1 & TXR_L1 == 7), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) scale(30.4375)
				save_max L1_TXD_Cont
				
				mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
				save_coefs L1_TXD_Cont
				}
				else { // If Rd is not included (need empty matrices)
					mata: bL1_TXD_Cont = .
					mata: rbL1_TXD_Cont = .
					mata: cbL1_TXD_Cont = .
					mata: fbL1_TXD_Cont = .
					mata: maxL1_TXD_Cont = .
					global Coeffs $Coeffs bL1_TXD_Cont rbL1_TXD_Cont cbL1_TXD_Cont fbL1_TXD_Cont maxL1_TXD_Cont
				}
				*mata: _matrix_list(bL1_TXD_Cont, rbL1_TXD_Cont, cbL1_TXD_Cont)

		// ASCT
		
		mi estimate: logit SCT Age Age2 Male i.ECOGcc i.RISS i.BCR Age70 Age75 i.CMc if(Event1 == 11 & TXR_L1 != 7 & BCR != 6)
		save_coefs L1_SCT	
		
		mata: _matrix_list(bL1_SCT, rbL1_SCT, cbL1_SCT)
				
		// ASCT Best Clinical Response
		
		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L1 if(Event0 == 100 & BCR_L1 != 6)	
		save_coefs SCT_BCR
		
		mata: _matrix_list(bSCT_BCR, rbSCT_BCR, cbSCT_BCR)
			
		// MNT

			// ASCT
			mi estimate: logit MNT Age Age2 Male i.ECOGcc i.RISS i.TXR_L1 i.BCR_SCT if(Event1 == 11 & SCT == 1)
			save_coefs MNT_ASCT
			
			mata: _matrix_list(bMNT_ASCT, rbMNT_ASCT, cbMNT_ASCT)
			
			// No ASCT
			mi estimate: logit MNT Age Age2 Male i.ECOGcc i.RISS i.TXR_L1 i.BCR_L1 if(Event1 == 11 & SCT == 0)
			save_coefs MNT_NoASCT	
			
			mata: _matrix_list(bMNT_NoASCT, rbMNT_NoASCT, cbMNT_NoASCT)	
		
		// Treatment-free Interval

			// ASCT 
			mi stset Date1 if(SCT == 1 & BCR_SCT != 0), id(ID_BS) failure(Event1 == 20) origin(Event1 == 11) scale(30.4375)
			save_max L1_TFI_ASCT
			
			mi estimate: streg Age Age2 Male i.ECOGcc i.RISS MNT i.BCR_SCT, d($dTFI)	
			save_coefs L1_TFI_ASCT	
			
			mata: _matrix_list(bL1_TFI_ASCT, rbL1_TFI_ASCT, cbL1_TFI_ASCT)
				
			// No ASCT
			mi stset Date1 if(SCT == 0), id(ID_BS) failure(Event1 == 20) origin(Event1 == 11) scale(30.4375)
			save_max L1_TFI_NoASCT
			
			mi estimate: streg Age Age2 Male i.ECOGcc i.RISS MNT i.BCR_L1, d($dTFI)	
			save_coefs L1_TFI_NoASCT
			
			mata: _matrix_list(bL1_TFI_NoASCT, rbL1_TFI_NoASCT, cbL1_TFI_NoASCT)

	di "Line 2" 
		
		// Treatment Regimen
	
		qui tab TXR_L2 // Check if modelling specific regimens
		if `r(r)' > 1 {
			mi estimate: mlogit TXR_L2 Age Age2 Male i.ECOGcc i.RISS i.pBCR if(Event0 == 20 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
			save_coefs L2_TXR
			
			mata: _matrix_list(bL2_TXR, rbL2_TXR, cbL2_TXR)
		}
		else { // If all patients receive 'other'
			mata: oL2_TXR = 0
			global Coeffs $Coeffs oL2_TXR
		}
		
		// Best Clinical Response

		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L1 i.BCR_SCT i.TXR_L2 if(Event0 == 20)
		save_coefs L2_BCR
		
		mata: _matrix_list(bL2_BCR, rbL2_BCR, cbL2_BCR)
			
		// Treatment Duration

		mi stset Date1 if(F_L2 != 1 & CN_L2 != 1), id(ID_BS) failure(Event1 == 21) origin(Event1 == 20) scale(30.4375)
		save_max L2_TXD
		
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L2 i.TXR_L2, d($dTXD)
		save_coefs L2_TXD	
		
		mata: _matrix_list(bL2_TXD, rbL2_TXD, cbL2_TXD)

		// Treatment-free Interval

		mi stset Date1, id(ID_BS) failure(Event1 == 30) origin(Event1 == 21) scale(30.4375)
		save_max L2_TFI	
		
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L2, d($dTFI)	
		save_coefs L2_TFI	
		
		mata: _matrix_list(bL2_TFI, rbL2_TFI, cbL2_TFI)

	di "Line 3"
	
		// Treatment Regimen

		qui tab TXR_L3 // Check if modelling specific regimens
		if `r(r)' > 1 {
			mi estimate: mlogit TXR_L3 Age Age2 Male i.ECOGcc i.RISS i.BCR_L2 if(Event0 == 30 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
			save_coefs L3_TXR
			
			mata: _matrix_list(bL3_TXR, rbL3_TXR, cbL3_TXR)
		}
		else { // If all patients receive 'other'
			mata: oL3_TXR = 0
			global Coeffs $Coeffs oL3_TXR
		}
		
		// Best Clinical Response
		
		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L2 i.TXR_L3 if(Event0 == 30)
		save_coefs L3_BCR
		
		mata: _matrix_list(bL3_BCR, rbL3_BCR, cbL3_BCR)
			
		// Treatment Duration

		mi stset Date1 if(F_L3 != 1 & CN_L3 != 1), id(ID_BS) failure(Event1 == 31) origin(Event1 == 30) scale(30.4375)
		save_max L3_TXD
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L3 i.TXR_L3, d($dTXD)
		save_coefs L3_TXD

		mata: _matrix_list(bL3_TXD, rbL3_TXD, cbL3_TXD)
		
		// Treatment-free Interval

		mi stset Date1, id(ID_BS) failure(Event1 == 40) origin(Event1 == 31) scale(30.4375)
		save_max L3_TFI	

		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L3, d($dTFI)	
		save_coefs L3_TFI

		mata: _matrix_list(bL3_TFI, rbL3_TFI, cbL3_TFI)
		
	di "Line 4"
	
		// Treatment Regimen

		qui tab TXR_L4 // Check if modelling specific regimens
		if `r(r)' > 1 {
			mi estimate: mlogit TXR_L4 Age Age2 Male i.ECOGcc i.RISS i.BCR_L3 if(Event0 == 40 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
			save_coefs L4_TXR
				
			mata: _matrix_list(bL4_TXR, rbL4_TXR, cbL4_TXR)
		}
		else { // If all patients receive 'other'
			mata: oL4_TXR = 0
			global Coeffs $Coeffs oL4_TXR
		}
			
		// Best Clinical Response
		
		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L3 i.TXR_L4 if(Event0 == 40)
		save_coefs L4_BCR		

		mata: _matrix_list(bL4_BCR, rbL4_BCR, cbL4_BCR)
		
		// Treatment Duration

		mi stset Date1 if(F_L4 != 1 & CN_L4 != 1), id(ID_BS) failure(Event1 == 41) origin(Event1 == 40) scale(30.4375)
		save_max L4_TXD	
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L4 i.TXR_L4, d($dTXD)
		save_coefs L4_TXD
			
		mata: _matrix_list(bL4_TXD, rbL4_TXD, cbL4_TXD)
		
		// Treatment-free Interval

		mi stset Date1, id(ID_BS) failure(Event1 == 50) origin(Event1 == 41) scale(30.4375)
		save_max L4_TFI
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L4, d($dTFI)	
		save_coefs L4_TFI	

		mata: _matrix_list(bL4_TFI, rbL4_TFI, cbL4_TFI)
		
	di "Line 5"
	
		// Best Clinical Response

		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L4 if(Event0 == 50)
		save_coefs L5_BCR		
		
		mata: _matrix_list(bL5_BCR, rbL5_BCR, cbL5_BCR)
		
		// Treatment Duration

		mi stset Date1 if(F_L5 != 1 & CN_L5 != 1), id(ID_BS) failure(Event1 == 51) origin(Event1 == 50) scale(30.4375)
		save_max L5_TXD	
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L5, d($dTXD)
		save_coefs L5_TXD
			
		mata: _matrix_list(bL5_TXD, rbL5_TXD, cbL5_TXD)
		
		// Treatment-free Interval

		mi stset Date1, id(ID_BS) failure(Event1 == 60) origin(Event1 == 51) scale(30.4375)
		save_max L5_TFI
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L5, d($dTFI)	
		save_coefs L5_TFI	

		mata: _matrix_list(bL5_TFI, rbL5_TFI, cbL5_TFI)
		
	di "Line 6"
/*	
		// Best Clinical Response

		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.pBCR if(Event0 == 60)
		save_coefs L6_BCR		
		
		mata: _matrix_list(bL6_BCR, rbL6_BCR, cbL6_BCR)
		
		// Treatment Duration

		mi stset Date1 if(F_L6 != 1 & CN_L6 != 1), id(ID_BS) failure(Event1 == 61) origin(Event1 == 60) scale(30.4375)
		save_max L6_TXD	
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
		save_coefs L6_TXD
			
		mata: _matrix_list(bL6_TXD, rbL6_TXD, cbL6_TXD)
		
		// Treatment-free Interval

		mi stset Date1, id(ID_BS) failure(Event1 == 70) origin(Event1 == 61) scale(30.4375)
		save_max L6_TFI
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTFI)	
		save_coefs L6_TFI	

		mata: _matrix_list(bL6_TFI, rbL6_TFI, cbL6_TFI)
*/	
	di "Line 6 to 9 (LX)"
	
		// Best Clinical Response

		mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.pBCR if(Event0 == 60 | Event0 == 70 | Event0 == 80 | Event0 == 90)
		save_coefs LX_BCR		
		
		mata: _matrix_list(bLX_BCR, rbLX_BCR, cbLX_BCR)	
		
		// Treatment Duration

		qui cap drop fail
		qui gen fail = 1 if(Event1 == 61 | Event1 == 71 | Event1 == 81 | Event1 == 91)
		qui replace fail = 0 if(Event1 == 60 | Event1 == 70 | Event1 == 80 | Event1 == 90)
		qui bysort ID_BS (Date0): replace fail = fail[_n-1] if(fail == .)
		
		mi stset Date1 if(F_L7 != 1 & F_L8 != 1 & F_L9 != 1), id(ID_BS) failure(fail) origin(Event1 == 60) exit(time .) scale(30.4375)
		save_max LX_TXD	
		
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
		save_coefs LX_TXD

		mata: _matrix_list(bLX_TXD, rbLX_TXD, cbLX_TXD)
		
		// Treatment-free Interval
		
		qui cap drop fail
		qui gen fail = 1 if(Event1 == 70 | Event1 == 80 | Event1 == 90)
		qui replace fail = 0 if(Event1 == 61 | Event1 == 71 | Event1 == 81)
		qui bysort ID_BS (Date0): replace fail = fail[_n-1] if(fail == .)
		
		mi stset Date1, id(ID_BS) failure(fail) origin(Event1 == 61) exit(time .) scale(30.4375)
		save_max LX_TFI
			
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTFI)	
		save_coefs LX_TFI	

		mata: _matrix_list(bLX_TFI, rbLX_TFI, cbLX_TFI)
	
end
	
**********
// Execute based on arguments

if("$boot" == "0") {
		
	// Load data
	use "${data_path}/${re_indir}MRDR Long MI${re_intag}.dta", replace
			
	// Execute function
	risk_equations

	// Save matrices
	mata: mata matsave "analyses/$analysis/coefficients/coefficients_$coeffs" $Coeffs, replace
}
else if("$boot" == "1") {
		
	forval b = $min_bs/$max_bs {
		
		// Load data
		use "${data_path}/${re_indir}bootstrap/MRDR Long MI${re_intag} B`b'.dta", replace
			
		// Shift dates so bootstrap samples of the same ID_BS do not have the same Date0/Date1
		qui bysort ID Event0 (ID_BS): gen tag0 = _n if Date0[_n-1] == Date0 & Date0 != .
		qui bysort ID Event1 (ID_BS): gen tag1 = _n if Date1[_n-1] == Date1 & Date1 != .
		qui replace Date0 = Date0 + tag0 - 1 if tag0 != .
		qui replace Date1 = Date1 + tag1 - 1 if tag1 != .
		qui drop tag0 tag1
			
		// Execute function
		risk_equations
				
		// Save matrices
		mata: mata matsave "analyses/$analysis/coefficients/bootstrap/coefficients_${coeffs}_B`b'" $Coeffs, replace
	}
}
