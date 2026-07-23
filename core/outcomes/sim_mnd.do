**********
* Monash Myeloma Model - Sim MND (L1 maintenance duration)
*
* Purpose: Draw L1 maintenance DURATION by parametric survival, among patients on maintenance
*          (MNT == 1). Continuous time in months.
*
* Notes:   TWO ARMS, SPLIT BY REGIMEN AND POOLED ACROSS TRANSPLANT - the opposite way round from the
*          earlier version, and for a measured reason. Lenalidomide and thalidomide are different
*          processes (thalidomide a fixed course, 9% censored; lenalidomide running to progression,
*          47%) with non-overlapping ancillaries, and an AFT carries one ancillary per equation. The
*          two TRANSPLANT arms have near-identical KM curves, so SCT is a covariate instead. BCR
*          drops out as a consequence, since the arms key on different response variables. Full
*          reasoning in prep/risk_equations.do.
*
*          NO ln(TFI) COVARIATE. The old fit conditioned duration on the drawn gap, which restricted
*          the FIT to the quarter of patients with an observed L2 - and the short-maintenance
*          quarter at that. Removing it uses the whole maintenance population.
*
*          THE ORDERING IS NOW ENFORCED DOWNSTREAM, NOT HERE. Dropping the gap covariate alone would
*          leave MND and TFI independent, and ~40% of patients would draw maintenance longer than
*          their own gap - which process_data.do's clip would then pull back, destroying the median.
*          Instead sim_tfi_l1.do runs AFTER this file and draws the gap TRUNCATED below at the
*          maintenance already drawn. Nothing is clipped, so nothing is dragged down.
*
*          THALIDOMIDE IS CAPPED AT 18 MONTHS, matching the censoring in its fit. 21 of 289 complete
*          episodes have recorded ends beyond 18 months on a drug given as a ~12-month course and not
*          prescribed in Australia since 2020; they carry ~29% of all thalidomide maintenance months.
*          A documented judgement, not a measurement - see risk_equations.do.
*
* ORDER:   AFTER sim_mnr.do (needs vMNR) and sim_bcr_asct.do (needs mBCR).
*          BEFORE sim_tfi_l1.do, which now depends on vMND. This is a REVERSAL of the previous order
*          and the whole point of the design: maintenance first, then the gap that must contain it.
*
*          Must match the fits in prep/risk_equations.do:
*              len   streg Age Age2 Male i.ECOGcc i.RISS SCT
*              thal  streg Age Age2 Male i.ECOGcc i.RISS SCT   (exit at 18 months)
*          e(b) order: Age Age2 Male, ECOG(0,1,2), RISS(1,2,3), SCT, _cons, aux. The base-level
*          dummies are 0-coef in e(b), as sim_tfi_l1.do does.
**********

mata {
	vCoefL = get_mnd_coef_len()
	vCoefT = get_mnd_coef_thal()

	if (cols(vCoefL) > 0 | cols(vCoefT) > 0) {

		// Alive, eligible AND receiving maintenance. Same population sim_mnr drew for.
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idx) > 0) {

			// ---- Lenalidomide ----
			if (cols(vCoefL) > 0) {
				iL = idx[selectindex(vMNR[idx] :== 1)]
				if (rows(iL) > 0) {
					mPatL = (vAge[iL], vAge2[iL], vMale[iL],
							 vECOG0[iL], vECOG1[iL], vECOG2[iL],
							 vRISS1[iL], vRISS2[iL], vRISS3[iL],
							 vSCT_L1[iL],
							 vCons[iL])
					nPredL = cols(mPatL)

					if (cols(vCoefL) != nPredL + 1) {
						errprintf("sim_mnd (len): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_LEN has %g. An ECOG/RISS level was likely empty in the fit.\n",
							nPredL, nPredL + 1, cols(vCoefL))
						exit(459)
					}

					vBetaL = vCoefL[1, 1..nPredL]'
					auxL   = vCoefL[1, cols(vCoefL)]
					vXBl   = mPatL * vBetaL
					vRNl   = rnDraw(iL, rn_mnd())
					vOCl   = calcSurvTime(vXBl, vRNl, fbL1_MND_LEN, auxL)
					vMND[iL] = rowmin((vOCl, J(rows(iL), 1, maxL1_MND_LEN)))
				}
			}

			// ---- Thalidomide ----
			if (cols(vCoefT) > 0) {
				iT = idx[selectindex(vMNR[idx] :== 5)]
				if (rows(iT) > 0) {
					mPatT = (vAge[iT], vAge2[iT], vMale[iT],
							 vECOG0[iT], vECOG1[iT], vECOG2[iT],
							 vRISS1[iT], vRISS2[iT], vRISS3[iT],
							 vSCT_L1[iT],
							 vCons[iT])
					nPredT = cols(mPatT)

					if (cols(vCoefT) != nPredT + 1) {
						errprintf("sim_mnd (thal): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_THAL has %g. An ECOG/RISS level was likely empty in the fit.\n",
							nPredT, nPredT + 1, cols(vCoefT))
						exit(459)
					}

					vBetaT = vCoefT[1, 1..nPredT]'
					auxT   = vCoefT[1, cols(vCoefT)]
					vXBt   = mPatT * vBetaT
					vRNt   = rnDraw(iT, rn_mnd())
					vOCt   = calcSurvTime(vXBt, vRNt, fbL1_MND_THAL, auxT)
					// maxL1_MND_THAL is set to 18 in risk_equations.do, overriding the observed
					// maximum, because the records beyond that point are the ones the fit has just
					// declared untrustworthy.
					vMND[iT] = rowmin((vOCt, J(rows(iT), 1, maxL1_MND_THAL)))
				}
			}
		}
	}
	else {
		// No model - leave vMND missing. process_data.do bills nothing where the duration is
		// missing, which is the safe direction: no maintenance cost beats silently reverting to
		// the old whole-gap bill. sim_tfi_l1.do also falls back to an untruncated draw.
		errprintf("sim_mnd: bL1_MND_LEN / bL1_MND_THAL not found - maintenance duration will not be costed. Re-run prep/risk_equations.do.\n")
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnd_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
