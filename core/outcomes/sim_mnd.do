**********
* Monash Myeloma Model - Sim MND (L1 maintenance duration)
*
* Purpose: Draw L1 maintenance DURATION by parametric survival, among patients on maintenance
*          (MNT == 1). Continuous time in months. The maintenance analogue of sim_txd_l1.do, and
*          split by transplant exactly as sim_tfi_l1.do is.
*
* Notes:   SIMPLE-FIRST design (docs/refractory.md 4.4). A survival fit on the duration uses every
*          maintenance patient - those who came off (failure) and those still on maintenance at
*          the data cut (censored) - so it sidesteps the complete-gap selection that a share of
*          TFI could not. process_data.do CAPS the drawn duration at the realised TFI_L1: a draw
*          that overshoots the gap means the patient stayed on maintenance until relapse, and
*          capping at the realised gap also inherits sim_mort's death curtailment. The cap is
*          therefore applied at billing time, not here.
*
*          Two arms, keyed on transplant like L1_TFI: the ASCT arm carries i.BCR_SCT, the no-ASCT
*          arm i.BCR_L1. Lenalidomide and thalidomide only; sim_mnr never draws 'other'.
*
* ORDER:   AFTER sim_mnr.do (needs vMNR), sim_bcr_asct.do (needs mBCR) AND sim_tfi_l1.do (needs the
*          drawn gap mTFI[.,2] for the ln(TFI) covariate). process_data.do also caps at the realised
*          TFI_L1 later.
*
*          Must match the fit in prep/risk_equations.do:
*              streg Age Age2 Male i.ECOGcc i.RISS i.MNR_L1 i.BCR_SCT|i.BCR_L1 MND_lntfi MND_lntfi_thal
*          e(b) order: Age Age2 Male, ECOG(0,1,2), RISS(1,2,3), MNR(1,5), BCR(...), ln(TFI),
*          ln(TFI)xthal, _cons, aux. ln(TFI) makes lenalidomide duration scale with the gap;
*          the base-level dummies are 0-coef in e(b), as sim_tfi_l1.do does.
**********

mata {
	vCoefA = get_mnd_coef_asct()
	vCoefN = get_mnd_coef_noasct()

	if (cols(vCoefA) > 0 | cols(vCoefN) > 0) {

		// Alive, eligible AND receiving maintenance. Same population sim_mnr drew for.
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idx) > 0) {

			// ---- ASCT arm: i.BCR_SCT (mBCR column 10, levels 1-4) ----
			if (cols(vCoefA) > 0) {
				iA = idx[selectindex(vSCT_L1[idx] :== 1)]
				if (rows(iA) > 0) {
					vB1 = (mBCR[iA, 10] :== 1)
					vB2 = (mBCR[iA, 10] :== 2)
					vB3 = (mBCR[iA, 10] :== 3)
					vB4 = (mBCR[iA, 10] :== 4)

					// ln(drawn TFI_L1, months) + thalidomide interaction: lenalidomide duration
					// scales with the gap, thalidomide does not (see risk_equations.do). Floor the
					// gap so a draw that rounded to 0 does not give ln(0).
					vLnTa = ln(rowmax((mTFI[iA, 2], J(rows(iA), 1, 0.01))))

					mPatA = (vAge[iA], vAge2[iA], vMale[iA],
							 vECOG0[iA], vECOG1[iA], vECOG2[iA],
							 vRISS1[iA], vRISS2[iA], vRISS3[iA],
							 (vMNR[iA] :== 1), (vMNR[iA] :== 5),
							 vB1, vB2, vB3, vB4,
							 vLnTa, (vLnTa :* (vMNR[iA] :== 5)),
							 vCons[iA])
					nPredA = cols(mPatA)

					if (cols(vCoefA) != nPredA + 1) {
						errprintf("sim_mnd (ASCT): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_ASCT has %g. A BCR/regimen level was likely empty in the fit.\n",
							nPredA, nPredA + 1, cols(vCoefA))
						exit(459)
					}

					vBetaA = vCoefA[1, 1..nPredA]'
					auxA   = vCoefA[1, cols(vCoefA)]
					vXBa   = mPatA * vBetaA
					vRNa   = rnDraw(iA, rn_mnd())
					vOCa   = calcSurvTime(vXBa, vRNa, fbL1_MND_ASCT, auxA)
					vMND[iA] = rowmin((vOCa, J(rows(iA), 1, maxL1_MND_ASCT)))
				}
			}

			// ---- No-ASCT arm: i.BCR_L1 (mBCR column 1, levels 1-6) ----
			if (cols(vCoefN) > 0) {
				// RESPONDERS only (BCR_L1 in 1-4): maintenance is a post-response therapy, so
				// SD/PD (5/6) do not get it. Same restriction as the fit (risk_equations.do); the
				// few SD/PD patients with vMNT == 1 simply get no duration and so no cost.
				iN = idx[selectindex((vSCT_L1[idx] :== 0) :& (mBCR[idx, 1] :<= 4))]
				if (rows(iN) > 0) {
					wB1 = (mBCR[iN, 1] :== 1)
					wB2 = (mBCR[iN, 1] :== 2)
					wB3 = (mBCR[iN, 1] :== 3)
					wB4 = (mBCR[iN, 1] :== 4)

					vLnTn = ln(rowmax((mTFI[iN, 2], J(rows(iN), 1, 0.01))))

					mPatN = (vAge[iN], vAge2[iN], vMale[iN],
							 vECOG0[iN], vECOG1[iN], vECOG2[iN],
							 vRISS1[iN], vRISS2[iN], vRISS3[iN],
							 (vMNR[iN] :== 1), (vMNR[iN] :== 5),
							 wB1, wB2, wB3, wB4,
							 vLnTn, (vLnTn :* (vMNR[iN] :== 5)),
							 vCons[iN])
					nPredN = cols(mPatN)

					if (cols(vCoefN) != nPredN + 1) {
						errprintf("sim_mnd (NoASCT): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_NoASCT has %g. A BCR/regimen level was likely empty in the fit (SD/PD patients rarely get maintenance).\n",
							nPredN, nPredN + 1, cols(vCoefN))
						exit(459)
					}

					vBetaN = vCoefN[1, 1..nPredN]'
					auxN   = vCoefN[1, cols(vCoefN)]
					vXBn   = mPatN * vBetaN
					vRNn   = rnDraw(iN, rn_mnd())
					vOCn   = calcSurvTime(vXBn, vRNn, fbL1_MND_NoASCT, auxN)
					vMND[iN] = rowmin((vOCn, J(rows(iN), 1, maxL1_MND_NoASCT)))
				}
			}
		}
	}
	else {
		// No model - leave vMND missing. process_data.do bills nothing where the duration is
		// missing, which is the safe direction: no maintenance cost beats silently reverting to
		// the old whole-gap bill.
		errprintf("sim_mnd: bL1_MND_ASCT / bL1_MND_NoASCT not found - maintenance duration will not be costed. Re-run prep/risk_equations.do.\n")
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnd_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
