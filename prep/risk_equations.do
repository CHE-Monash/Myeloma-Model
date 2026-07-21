**********
* Monash Myeloma Model - Risk Equations
*
* Purpose: Fit every model risk equation (grouped by outcome type) from the MI data and matsave
*          the coefficients to analyses/$analysis/coefficients/. Called by the prep pipeline.
* Usage:   do risk_equations.do <analysis> <coeffs> <minyr> <maxyr> <boot> <minbs> <maxbs> [sample]
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
* $analysis, so OOS coefficients land in analyses/default/coefficients/). Empty = main model.
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

cap program drop gen_mnr
program define gen_mnr
	// MNR_L1 arrives from the extraction as a RAW 6-level maintenance drug code (see
	// docs/refractory.md 2). Keep that as MNR_drug and rebuild MNR_L1 as the levels this
	// analysis models, per $MNR_L1 from outcomes/mnr_<coeffs>.do; anything unlisted falls to
	// 0 = 'other'. Mirrors gen_txr, and is what lets one extraction serve both the historical
	// window the OOS validation scores against and a current-paradigm window: the maintenance
	// mix stepped in 2020, so no single fixed list serves both (docs/refractory.md 7.4).
	// Re-entrant: rebuilds from MNR_drug if already called.
	cap confirm variable MNR_drug
	if _rc {
		cap confirm variable MNR_L1
		if _rc {
			di as err "gen_mnr: MNR_L1 is not in the MI data."
			di as err "         Rebuild MRDR Long (prep/data_extraction.do), then re-run"
			di as err "         prep/multiple_imputation.do - its keep list must name MNR_L1."
			exit 111
		}
		rename MNR_L1 MNR_drug
	}
	cap drop MNR_L1
	gen byte MNR_L1 = 0 if !mi(MNR_drug)
	foreach r of global MNR_L1 {
		replace MNR_L1 = `r' if MNR_drug == `r'
	}
	cap label values MNR_L1 MNR_label
end

cap program drop gen_txr
program define gen_txr
	// Build TXR_L1..L9 from the per-line regimen code lists in $TXR_L1..$TXR_L9 (declared by the
	// analysis's outcomes/txr_<coeffs>.do). Any regimen not listed for a line falls into 0 ('other').
	forval l = 1/9 {
		cap drop TXR_L`l'
		gen TXR_L`l' = 0 if Event0 == `l'0
		foreach r of global TXR_L`l' {
			replace TXR_L`l' = `r' if Event0 == `l'0 & Regimen == `r'
		}
		bysort ID (TXR_L`l'): replace TXR_L`l' = TXR_L`l'[_n-1] if TXR_L`l' == .
	}
end

cap program drop risk_equations
program define risk_equations

	// OS distribution
	global dOS "w"

	// TXD distribution
	global dTXD "w"

	// TFI distribution
	global dTFI = "lognormal"

	// TXR variables
	forval l = 1/9 {
		global TXR_L`l' ""
	}
	qui do "analyses/$analysis/outcomes/txr_$coeffs.do"
	gen_txr

	// MNR variable. Not every analysis declares a maintenance regimen list; without one every
	// regimen falls to 'other' and L1_MNR is skipped below, which is right for the $line 2
	// analyses where maintenance is never costed (docs/refractory.md 7.3).
	global MNR_L1 ""
	cap qui do "analyses/$analysis/outcomes/mnr_$coeffs.do"
	gen_mnr

	// Reset global
	global Coeffs

	// Risk equations are grouped BY OUTCOME TYPE (not by line). Each block below is a
	// self-contained estimation (its own stset/if), so grouping is purely organisational -
	// $Coeffs/matsave key by name, so order does not affect the saved coefficients. Within a
	// group the equations run in pathway order (diagnosis -> L1 -> ... -> LX). See the engine's
	// core/outcomes/sim_*.do for how each is consumed.

	***** OVERALL SURVIVAL *****
	di "Overall Survival"
	// Per-line OS family: one Weibull per pathway stage, window-censored via origin()/exit().
	// Same covariates throughout (Age Age2 Male i.ECOGcc i.RISS CM_* + stage BCR); stages differ
	// only by the stset window and BCR source. See core/outcomes/sim_os.do for the firing map.

	// OS_DN: DN -> L1S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 3) exit(Event1 == 10) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT, d($dOS)
	save_coefs OS_DN
	mata: _matrix_list(bOS_DN, rbOS_DN, cbOS_DN)

	// OS_L1S: L1S -> L1E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 10) exit(Event1 == 11) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L1, d($dOS)
	save_coefs OS_L1S
	mata: _matrix_list(bOS_L1S, rbOS_L1S, cbOS_L1S)

	// OS_L1E_NoASCT: L1E -> L2S, No ASCT
	mi stset Date1 if(F_OS != 1 & SCT == 0), id(ID_BS) failure(Event1 == 104) origin(Event1 == 11) exit(Event1 == 20) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L1, d($dOS)
	save_coefs OS_L1E_NoASCT
	mata: _matrix_list(bOS_L1E_NoASCT, rbOS_L1E_NoASCT, cbOS_L1E_NoASCT)

	// OS_L1E_ASCT: L1E -> L2S, ASCT (BCR_SCT, 4 levels)
	mi stset Date1 if(F_OS != 1 & SCT == 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 11) exit(Event1 == 20) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_SCT, d($dOS)
	save_coefs OS_L1E_ASCT
	mata: _matrix_list(bOS_L1E_ASCT, rbOS_L1E_ASCT, cbOS_L1E_ASCT)

	// OS_L2S: L2S -> L2E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 20) exit(Event1 == 21) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L2 LenRefr_Tx_in, d($dOS)
	save_coefs OS_L2S
	mata: _matrix_list(bOS_L2S, rbOS_L2S, cbOS_L2S)

	// OS_L2E: L2E -> L3S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 21) exit(Event1 == 30) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L2 LenRefr_Tx_in, d($dOS)
	save_coefs OS_L2E
	mata: _matrix_list(bOS_L2E, rbOS_L2E, cbOS_L2E)

	// OS_L3S: L3S -> L3E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 30) exit(Event1 == 31) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L3 LenRefr_Tx_in, d($dOS)
	save_coefs OS_L3S
	mata: _matrix_list(bOS_L3S, rbOS_L3S, cbOS_L3S)

	// OS_L3E: L3E -> L4S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 31) exit(Event1 == 40) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L3 LenRefr_Tx_in, d($dOS)
	save_coefs OS_L3E
	mata: _matrix_list(bOS_L3E, rbOS_L3E, cbOS_L3E)

	// OS_L4S: L4S -> L4E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 40) exit(Event1 == 41) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L4 LenRefr_Tx_in, d($dOS)
	save_coefs OS_L4S
	mata: _matrix_list(bOS_L4S, rbOS_L4S, cbOS_L4S)

	// OS_L4E: L4E -> L5S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 41) exit(Event1 == 50) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L4 LenRefr_Tx_in, d($dOS)
	save_coefs OS_L4E
	mata: _matrix_list(bOS_L4E, rbOS_L4E, cbOS_L4E)

	// OS_L5S: L5S -> L5E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 50) exit(Event1 == 51) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L5, d($dOS)
	save_coefs OS_L5S
	mata: _matrix_list(bOS_L5S, rbOS_L5S, cbOS_L5S)

	// OS_L5E: L5E -> L6S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 51) exit(Event1 == 60) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L5, d($dOS)
	save_coefs OS_L5E
	mata: _matrix_list(bOS_L5E, rbOS_L5E, cbOS_L5E)

	// OS_L6plus: L6S onward (single conditional model for the sparse deep tail; running BCR)
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 60) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR, d($dOS)
	save_coefs OS_L6plus
	mata: _matrix_list(bOS_L6plus, rbOS_L6plus, cbOS_L6plus)

	***** ASCT (TRANSPLANT DECISION) *****
	di "ASCT"
	// Transplant-eligibility logits: at diagnosis (DN_SCT) and at L1 end (L1_SCT).
	// Comorbidities: 4 individual flags (CKD/cardiac/lung/diabetes).

	// DN
	mi estimate: logit SCT Age Age2 Male i.ECOGcc i.RISS Age70 Age75 CM_CKD CM_CRD CM_PLM CM_DBT if(Event0 == 3)
	save_coefs DN_SCT
	mata: _matrix_list(bDN_SCT, rbDN_SCT, cbDN_SCT)

	// L1
	mi estimate, esampvaryok: logit SCT Age Age2 Male i.ECOGcc i.RISS i.BCR Age70 Age75 CM_CKD CM_CRD CM_PLM CM_DBT if(Event1 == 11 & TXR_L1 != 7 & BCR != 6)
	save_coefs L1_SCT
	mata: _matrix_list(bL1_SCT, rbL1_SCT, cbL1_SCT)

	***** TREATMENT REGIMEN (TXR) *****
	di "Treatment Regimen"
	// Regimen choice is availability-driven; L1 carries Age Age2 SCT (transplant gates some
	// regimens), L2+ carry Age plus LenRefr_Tx_in - Male/ECOG/RISS/prior-BCR carried no signal (see
	// scratch/txr_predictor_check.do), but len-refractory status drives L2/L3 regimen choice almost
	// deterministically (docs/refractory.md 3.4), which Age alone cannot reproduce. L1 gets no flag:
	// LenRefr_Tx_in is 0 at L1 by construction. Engine design mPat in core/outcomes/sim_txr.do must
	// match. If a line has only 'other' (r(r) == 1) an empty o<L>_TXR = 0 is stored instead.

	// L1
	qui tab TXR_L1 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L1 Age Age2 SCT if(Event0 == 10 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L1_TXR
		mata: _matrix_list(bL1_TXR, rbL1_TXR, cbL1_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL1_TXR = 0
		global Coeffs $Coeffs oL1_TXR
	}

	// L2
	qui tab TXR_L2 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L2 Age Age2 LenRefr_Tx_in if(Event0 == 20 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L2_TXR
		mata: _matrix_list(bL2_TXR, rbL2_TXR, cbL2_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL2_TXR = 0
		global Coeffs $Coeffs oL2_TXR
	}

	// L3
	qui tab TXR_L3 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L3 Age Age2 LenRefr_Tx_in if(Event0 == 30 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L3_TXR
		mata: _matrix_list(bL3_TXR, rbL3_TXR, cbL3_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL3_TXR = 0
		global Coeffs $Coeffs oL3_TXR
	}

	// L4
	qui tab TXR_L4 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L4 Age Age2 LenRefr_Tx_in if(Event0 == 40 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L4_TXR
		mata: _matrix_list(bL4_TXR, rbL4_TXR, cbL4_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL4_TXR = 0
		global Coeffs $Coeffs oL4_TXR
	}

	// L5
	qui tab TXR_L5 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L5 Age Age2 if(Event0 == 50 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L5_TXR
		mata: _matrix_list(bL5_TXR, rbL5_TXR, cbL5_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL5_TXR = 0
		global Coeffs $Coeffs oL5_TXR
	}

	***** BEST CLINICAL RESPONSE (BCR) *****
	di "Best Clinical Response"
	// Ordered logit for response at each line; each line conditions on the prior line's response
	// (and regimen). SCT_BCR is the post-transplant response (conditional on BCR_L1).

	// L1
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS SCT i.TXR_L1 if(Event0 == 10)
	save_coefs L1_BCR
	mata: _matrix_list(bL1_BCR, rbL1_BCR, cbL1_BCR)

	// ASCT (post-transplant response)
	mi estimate, esampvaryok: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L1 if(Event0 == 100 & BCR_L1 != 6)
	save_coefs SCT_BCR
	mata: _matrix_list(bSCT_BCR, rbSCT_BCR, cbSCT_BCR)

	// L2
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L1 i.BCR_SCT i.TXR_L2 if(Event0 == 20)
	save_coefs L2_BCR
	mata: _matrix_list(bL2_BCR, rbL2_BCR, cbL2_BCR)

	// L3
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L2 i.TXR_L3 if(Event0 == 30)
	save_coefs L3_BCR
	mata: _matrix_list(bL3_BCR, rbL3_BCR, cbL3_BCR)

	// L4
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L3 i.TXR_L4 if(Event0 == 40)
	save_coefs L4_BCR
	mata: _matrix_list(bL4_BCR, rbL4_BCR, cbL4_BCR)

	// L5
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L4 if(Event0 == 50)
	save_coefs L5_BCR
	mata: _matrix_list(bL5_BCR, rbL5_BCR, cbL5_BCR)

	// LX (lines 6-9)
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.pBCR if(Event0 == 60 | Event0 == 70 | Event0 == 80 | Event0 == 90)
	save_coefs LX_BCR
	mata: _matrix_list(bLX_BCR, rbLX_BCR, cbLX_BCR)

	// L6 (superseded by LX; kept for reference)
/*
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.pBCR if(Event0 == 60)
	save_coefs L6_BCR
	mata: _matrix_list(bL6_BCR, rbL6_BCR, cbL6_BCR)
*/

	***** LENALIDOMIDE-REFRACTORY, TREATMENT LINES (LENREFR_TX) *****
	di "Lenalidomide-refractory (treatment lines)"
	// LenRefr_Tx is a LATCHED state: once refractory to a lenalidomide treatment line the patient
	// stays refractory (cumulative-any, and downstream TXR/OS read only that latched flag). So the
	// only event to model is the 0 -> 1 FLIP, which can only happen to a NOT-YET-REFRACTORY patient.
	// The engine (sim_lenrefr.do) therefore draws only where LenRefr_Tx_in == 0:
	//     BCR in {5,6}  -> refractory        (definitional, no equation)
	//     BCR in {1-4}  -> Bernoulli(this logit)
	// and leaves already-refractory patients latched at 1. This logit is fitted on exactly that
	// population - not-yet-refractory (LenRefr_Tx_in == 0), responding (BCR 1-4), len line-start
	// rows - so prior state is NOT a covariate (it is 0 by construction here). Fitting it on the
	// full sample with LenRefr_Tx_in as a covariate, as an earlier draft did, mismatches the engine's
	// application population (the 4.4 "test with what the engine will have" rule).
	//
	// POOLED across lines with the line collapsed to L1 / L2 / L3+ (LENREFR_line): docs/refractory.md
	// 3.5 settles that shape - full line resolution and the L4+ split each add nothing (LR p > 0.95),
	// and line is the dominant predictor. Carries the pre-specified baseline set REGARDLESS of
	// significance (selecting on in-sample p-values does not replicate out of sample; the OOS
	// validation is the arbiter), plus i.BCR.
	//
	// FIT vs ENGINE gate: the fit conditions on Lenalidomide == 1 (the true drug binary, the clean
	// residual-arm definition). The engine has no drug binary, so sim_lenrefr applies this where the
	// DRAWN regimen is a lenalidomide code (analysis-declared, 'other' treated as non-len). That
	// asymmetry is deliberate and noted in docs/refractory.md 3.5 / 4.
	cap drop LENREFR_line
	gen byte LENREFR_line = min(Line, 3)
	// esampvaryok: the residual arm (BCR 1-4) is defined on imputed BCR, so its membership varies
	// across imputations - legitimate, and pooled the same way the ASCT equations are (docs note).
	mi estimate, esampvaryok: logit LineRefr Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT ///
		i.LENREFR_line i.BCR ///
		if(CStart == 1 & inlist(Event0, 10, 20, 30, 40, 50, 60, 70, 80, 90) & Lenalidomide == 1 & inlist(BCR, 1, 2, 3, 4) & LenRefr_Tx_in == 0)
	save_coefs LENREFR_TX
	mata: _matrix_list(bLENREFR_TX, rbLENREFR_TX, cbLENREFR_TX)

	// Store the len-regimen codes alongside the coefficients so sim_lenrefr.do knows which drawn
	// regimens count as lenalidomide (the engine analogue of the fit's Lenalidomide == 1 gate).
	// Declared per analysis in outcomes/txr_$coeffs.do (default "7 31"); an analysis that declares
	// none saves nothing here, and sim_lenrefr.do becomes a no-op via lenrefr_model_exists().
	if ("$LENREFR_regimens" != "") {
		mata: LENREFR_regimens = strtoreal(tokens(st_global("LENREFR_regimens")))
		global Coeffs $Coeffs LENREFR_regimens
	}

	***** TREATMENT DURATION (TXD) *****
	di "Treatment Duration"
	// Time on treatment within a line (Weibull). L1 splits by transplant: ASCT fixed-duration is a
	// 3-spline fit (cut-offs 60/120 mo), NoASCT is a single fit, and Continuous (Rd) is separate.
	// L2+ are single fits conditional on that line's response and regimen.

	// L1 - Fixed w/ ASCT (3 splines)
	// Cut off 1 - 60 for optimal fit with 3 splines, 90 for 2 splines
	scalar L1_TXD_ASCT_C1 = 60
	mata: L1_TXD_ASCT_C1 = st_numscalar("L1_TXD_ASCT_C1")

	// Cut off 2 - 120 for optimal fit with 3 splines
	scalar L1_TXD_ASCT_C2 = 120
	mata: L1_TXD_ASCT_C2 = st_numscalar("L1_TXD_ASCT_C2")

	global Coeffs $Coeffs L1_TXD_ASCT_C1 L1_TXD_ASCT_C2

	// Spline 1
	gen S1_exit = Date1 + `=L1_TXD_ASCT_C1'
	mi stset Date1 if(TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) exit(S1_exit) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_ASCT_S1
	mata: _matrix_list(bL1_TXD_ASCT_S1, rbL1_TXD_ASCT_S1, cbL1_TXD_ASCT_S1)

	// Spline 2
	gen S2_enter = S1_exit + 1
	gen S2_exit = Date1 + `=L1_TXD_ASCT_C2'
	mi stset Date1 if(TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) enter(S2_enter) exit(S2_exit) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_ASCT_S2
	mata: _matrix_list(bL1_TXD_ASCT_S2, rbL1_TXD_ASCT_S2, cbL1_TXD_ASCT_S2)

	// Spline 3
	gen S3_enter = S2_exit + 1
	mi stset Date1 if(TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) enter(S3_enter) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_ASCT_S3
	drop S1_exit S2_enter S2_exit S3_enter
	mata: _matrix_list(bL1_TXD_ASCT_S3, rbL1_TXD_ASCT_S3, cbL1_TXD_ASCT_S3)

	// L1 - Fixed w/o ASCT
	mi stset Date1 if(TXR_L1 != 7 & SCT == 0), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_NoASCT
	mata: _matrix_list(bL1_TXD_NoASCT, rbL1_TXD_NoASCT, cbL1_TXD_NoASCT)

	// L1 - Continuous (Rd)
	count if(TXR_L1 == 7)
	if r(N) > 0 { // If Rd is included
		mi stset Date1 if(TXR_L1 == 7), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) scale(30.4375)
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

	// L2
	mi stset Date1, id(ID_BS) failure(Event1 == 21) origin(Event1 == 20) scale(30.4375)
	save_max L2_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L2 i.TXR_L2, d($dTXD)
	save_coefs L2_TXD
	mata: _matrix_list(bL2_TXD, rbL2_TXD, cbL2_TXD)

	// L3
	mi stset Date1, id(ID_BS) failure(Event1 == 31) origin(Event1 == 30) scale(30.4375)
	save_max L3_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L3 i.TXR_L3, d($dTXD)
	save_coefs L3_TXD
	mata: _matrix_list(bL3_TXD, rbL3_TXD, cbL3_TXD)

	// L4
	mi stset Date1, id(ID_BS) failure(Event1 == 41) origin(Event1 == 40) scale(30.4375)
	save_max L4_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L4 i.TXR_L4, d($dTXD)
	save_coefs L4_TXD
	mata: _matrix_list(bL4_TXD, rbL4_TXD, cbL4_TXD)

	// L5
	mi stset Date1, id(ID_BS) failure(Event1 == 51) origin(Event1 == 50) scale(30.4375)
	save_max L5_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L5, d($dTXD)
	save_coefs L5_TXD
	mata: _matrix_list(bL5_TXD, rbL5_TXD, cbL5_TXD)

	// LX (lines 6-9)
	qui cap drop fail
	qui gen fail = 1 if(Event1 == 61 | Event1 == 71 | Event1 == 81 | Event1 == 91)
	qui replace fail = 0 if(Event1 == 60 | Event1 == 70 | Event1 == 80 | Event1 == 90)
	qui bysort ID_BS (Date0): replace fail = fail[_n-1] if(fail == .)
	mi stset Date1, id(ID_BS) failure(fail) origin(Event1 == 60) exit(time .) scale(30.4375)
	save_max LX_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
	save_coefs LX_TXD
	mata: _matrix_list(bLX_TXD, rbLX_TXD, cbLX_TXD)

	// L6 (superseded by LX; kept for reference)
/*
	mi stset Date1, id(ID_BS) failure(Event1 == 61) origin(Event1 == 60) scale(30.4375)
	save_max L6_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
	save_coefs L6_TXD
	mata: _matrix_list(bL6_TXD, rbL6_TXD, cbL6_TXD)
*/

	***** MAINTENANCE (MNT) *****
	di "Maintenance"
	// Maintenance-therapy logit after L1 induction, split by transplant status.

	// ASCT
	mi estimate: logit MNT Age Age2 Male i.ECOGcc i.RISS i.TXR_L1 i.BCR_SCT if(Event1 == 11 & SCT == 1)
	save_coefs MNT_ASCT
	mata: _matrix_list(bMNT_ASCT, rbMNT_ASCT, cbMNT_ASCT)

	// No ASCT
	mi estimate: logit MNT Age Age2 Male i.ECOGcc i.RISS i.TXR_L1 i.BCR_L1 if(Event1 == 11 & SCT == 0)
	save_coefs MNT_NoASCT
	mata: _matrix_list(bMNT_NoASCT, rbMNT_NoASCT, cbMNT_NoASCT)

	***** MAINTENANCE REGIMEN AND DURATION (MNR / MND) *****
	di "Maintenance regimen and duration"
	// The maintenance analogue of TXR_L1 / TXD_L1, sat here because both condition on the MNT
	// logit above: MNT decides WHETHER, L1_MNR WHICH, L1_MND HOW LONG. Both fit at the same
	// Event1 == 11 row as MNT, and only among MNT == 1. Consumed only by cost_tx_mnt, so the
	// out-of-sample validation fits them but never uses them. See docs/refractory.md 7.4.

	// Which regimen. Lenalidomide and thalidomide only (docs/refractory.md 7.4) - lenalidomide
	// is the base, thalidomide the alternative, so this is effectively a logit and the engine
	// never produces an 'other' maintenance regimen. Year-windowed exactly as TXR_L1 is, so the
	// mix reflects the era; the list comes from outcomes/mnr_$coeffs.do ("1 5"). If a window
	// leaves only one regimen (r(r) == 1, e.g. thalidomide empty in a modern window) then
	// oL1_MNR = 1 is stored and sim_mnr assigns every maintenance patient lenalidomide.
	qui tab MNR_L1 if(Event1 == 11 & MNT == 1 & inlist(MNR_L1, 1, 5))
	if `r(r)' > 1 {
		mi estimate: mlogit MNR_L1 Age Age2 Male i.ECOGcc i.RISS SCT ///
			if(Event1 == 11 & MNT == 1 & inlist(MNR_L1, 1, 5) & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(1)
		save_coefs L1_MNR
		mata: _matrix_list(bL1_MNR, rbL1_MNR, cbL1_MNR)
	}
	else {
		mata: oL1_MNR = 1
		global Coeffs $Coeffs oL1_MNR
	}

	// How long: maintenance DURATION via parametric survival on the maintenance EVENTS in the
	// skeleton, exactly as L1_TFI is fitted (docs/refractory.md 4.4). The extraction keeps the L1
	// maintenance start (110) and end (111) events, so this is a normal stset:
	//     origin(Event1 == 110)          maintenance start
	//     failure(Event1 == 20 111)      maintenance end - the recorded end (111), or L2 start (20)
	//
	// GAP-DEPENDENT via ln(gap). Lenalidomide maintenance runs to progression, so its duration
	// scales with the gap; thalidomide is a fixed ~10-month course and does not. Without a gap term
	// the survival draw is gap-INDEPENDENT and the simulated share falls with gap length while the
	// registry share rises - so ln(gap) enters with a REGIMEN INTERACTION (MND_lntfi_thal), giving
	// lenalidomide a slope near 1 and thalidomide a slope near 0. In MONTHS, to match the engine's
	// mTFI (sim_mnd.do reads ln(drawn TFI_L1)).
	//
	// COMPLETE GAPS. ln(gap) uses TFI_L1 (L1E to L2), which is missing for patients still on
	// maintenance at the cut, so they drop from this fit. That is the price of a gap covariate that
	// is NOT derived from the censoring point: a censor-filled gap was tried and pinned the share at
	// 1.0, because for a censored patient the filled gap IS their own censoring time, and the
	// survival likelihood then inflates the predicted duration above the gap (docs/refractory.md 4.4).
	// The engine still uses the drawn (complete) gap, so the relationship transfers.
	//
	// SPLIT BY TRANSPLANT, exactly as L1_TFI is: the ASCT arm keys on the post-transplant response
	// (i.BCR_SCT), the no-ASCT arm on the post-induction response (i.BCR_L1). Otherwise the
	// covariates match. Lenalidomide and thalidomide only (inlist(MNR_L1, 1, 5)). Log-normal.
	//
	// The engine design matrices (sim_mnd.do) carry every factor level, so a level EMPTY in an arm
	// (e.g. no-ASCT SD/PD patients rarely get maintenance) drops from e(b) and trips the design
	// guard - collapse that level if so.
	qui gen double MND_lntfi      = ln(TFI_L1 / 30.4375)     // L1E-to-L2 gap, months; missing (dropped) where no L2
	qui gen double MND_lntfi_thal = MND_lntfi * (MNR_L1 == 5)

	// ASCT
	mi stset Date1 if(SCT == 1 & BCR_SCT != 0 & MNT == 1 & inlist(MNR_L1, 1, 5)), ///
		id(ID_BS) failure(Event1 == 20 111) origin(Event1 == 110) scale(30.4375)
	save_max L1_MND_ASCT
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.MNR_L1 i.BCR_SCT MND_lntfi MND_lntfi_thal, d($dTFI)
	save_coefs L1_MND_ASCT
	mata: _matrix_list(bL1_MND_ASCT, rbL1_MND_ASCT, cbL1_MND_ASCT)

	// No ASCT. Restrict to RESPONDERS (BCR_L1 in 1-4 = CR/VGPR/PR/MR): maintenance is a
	// post-response therapy, so SD/PD (5/6) do not get it by definition. The handful in the data
	// are noise, and on the small no-ASCT sample (smaller still on the OOS train fold) an SD or PD
	// cell empties and drops from i.BCR_L1, which trips the engine design guard. Dropping them
	// keeps i.BCR_L1 at a full 4 levels. sim_mnd.do excludes the same patients (docs 4.4). The
	// ASCT arm already needs no filter - BCR_SCT is collapsed to 1-4 at transplant.
	mi stset Date1 if(SCT == 0 & MNT == 1 & inlist(MNR_L1, 1, 5) & inlist(BCR_L1, 1, 2, 3, 4)), ///
		id(ID_BS) failure(Event1 == 20 111) origin(Event1 == 110) scale(30.4375)
	save_max L1_MND_NoASCT
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.MNR_L1 i.BCR_L1 MND_lntfi MND_lntfi_thal, d($dTFI)
	save_coefs L1_MND_NoASCT
	mata: _matrix_list(bL1_MND_NoASCT, rbL1_MND_NoASCT, cbL1_MND_NoASCT)

	***** TREATMENT-FREE INTERVAL (TFI) *****
	di "Treatment-free Interval"
	// Gap from the end of one line's treatment to the start of the next (log-normal). DN_TFI is the
	// diagnosis-to-L1 interval; L1 splits by transplant (with maintenance MNT); L2+ are single fits.

	// DN
	mi stset Date1, id(ID_BS) failure(Event1 == 10) origin(Event1 == 3) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS SCT, d($dTFI)
	save_coefs DN_TFI
	mata: _matrix_list(bDN_TFI, rbDN_TFI, cbDN_TFI)

	// L1 - ASCT
	mi stset Date1 if(SCT == 1 & BCR_SCT != 0), id(ID_BS) failure(Event1 == 20) origin(Event1 == 11) scale(30.4375)
	save_max L1_TFI_ASCT
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS MNT i.BCR_SCT, d($dTFI)
	save_coefs L1_TFI_ASCT
	mata: _matrix_list(bL1_TFI_ASCT, rbL1_TFI_ASCT, cbL1_TFI_ASCT)

	// L1 - No ASCT
	mi stset Date1 if(SCT == 0), id(ID_BS) failure(Event1 == 20) origin(Event1 == 11) scale(30.4375)
	save_max L1_TFI_NoASCT
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS MNT i.BCR_L1, d($dTFI)
	save_coefs L1_TFI_NoASCT
	mata: _matrix_list(bL1_TFI_NoASCT, rbL1_TFI_NoASCT, cbL1_TFI_NoASCT)

	// L2
	mi stset Date1, id(ID_BS) failure(Event1 == 30) origin(Event1 == 21) scale(30.4375)
	save_max L2_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L2, d($dTFI)
	save_coefs L2_TFI
	mata: _matrix_list(bL2_TFI, rbL2_TFI, cbL2_TFI)

	// L3
	mi stset Date1, id(ID_BS) failure(Event1 == 40) origin(Event1 == 31) scale(30.4375)
	save_max L3_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L3, d($dTFI)
	save_coefs L3_TFI
	mata: _matrix_list(bL3_TFI, rbL3_TFI, cbL3_TFI)

	// L4
	mi stset Date1, id(ID_BS) failure(Event1 == 50) origin(Event1 == 41) scale(30.4375)
	save_max L4_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L4, d($dTFI)
	save_coefs L4_TFI
	mata: _matrix_list(bL4_TFI, rbL4_TFI, cbL4_TFI)

	// L5
	mi stset Date1, id(ID_BS) failure(Event1 == 60) origin(Event1 == 51) scale(30.4375)
	save_max L5_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L5, d($dTFI)
	save_coefs L5_TFI
	mata: _matrix_list(bL5_TFI, rbL5_TFI, cbL5_TFI)

	// LX (lines 6-9)
	qui cap drop fail
	qui gen fail = 1 if(Event1 == 70 | Event1 == 80 | Event1 == 90)
	qui replace fail = 0 if(Event1 == 61 | Event1 == 71 | Event1 == 81)
	qui bysort ID_BS (Date0): replace fail = fail[_n-1] if(fail == .)
	mi stset Date1, id(ID_BS) failure(fail) origin(Event1 == 61) exit(time .) scale(30.4375)
	save_max LX_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTFI)
	save_coefs LX_TFI
	mata: _matrix_list(bLX_TFI, rbLX_TFI, cbLX_TFI)

	// L6 (superseded by LX; kept for reference)
/*
	mi stset Date1, id(ID_BS) failure(Event1 == 70) origin(Event1 == 61) scale(30.4375)
	save_max L6_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTFI)
	save_coefs L6_TFI
	mata: _matrix_list(bL6_TFI, rbL6_TFI, cbL6_TFI)
*/

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
