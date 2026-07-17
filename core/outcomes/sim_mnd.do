**********
* Monash Myeloma Model - Sim MND
*
* Purpose: Draw the L1 maintenance duration among patients receiving maintenance
*          (MNT == 1), as a SHARE of TFI_L1 rather than as a duration in months.
*          The maintenance analogue of sim_txd_l1.do.
*
* Notes:   THIS MODULE DRAWS A SHARE, NOT A DURATION. TFI_L1 contains the maintenance
*          duration (TFI_L1 = TTM + MND_L1 + trueTFI), so a duration drawn independently
*          has to be truncated back to the gap - and a truncated patient is billed for the
*          WHOLE gap, which is exactly the defect being fixed. Measured, not assumed: an
*          independently drawn duration curtailed at the patient's own observed TFI_L1
*          overshoots 42.5% of the time against an observed 0.7%. The share is bounded by
*          construction and cannot do that. See docs/refractory.md 7.4.
*
*          The multiplication MND_L1 = MNS_L1 * TFI_L1 happens in core/process_data.do, NOT
*          here, and that is deliberate on two counts. sim_tfi_l1.do runs AFTER this module,
*          so TFI_L1 does not exist yet; and mTFI is later curtailed at death by sim_mort.do,
*          so deferring the multiply means a patient who dies mid-gap inherits the
*          curtailment for free rather than being billed for a gap they never lived.
*
*          Consumed only by cost_tx_mnt, so this runs but does nothing in the $line 2
*          analyses, where maintenance is never costed (7.3).
*
*          Must match the fit in prep/risk_equations.do:
*              betareg MNS_L1 Age Age2 Male i.ECOGcc i.RISS SCT i.MNR_L1
*          betareg returns TWO equations: the mean (logit link) and the scale (log link).
*          save_coefs stores both, so the coefficient vector is [mean equation | scale], and
*          the FINAL column is ln(phi) - the same "aux is the last column" idiom the
*          parametric survival models use (see calcSurvTime in core/mata_functions.do).
**********

mata {
	if (mnd_model_exists()) {

		vCoef = get_mnd_coef()

		// Filter for alive, eligible AND receiving maintenance. Must be the same population
		// sim_mnr.do drew for - MNR_L1 is a covariate below.
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idx) > 0) {

			// Maintenance regimen dummies. The modelled drug codes come from the analysis's
			// mnr_$coeffs.do, so build one dummy per level in the SAME order gen_mnr / mlogit
			// saw them (ascending code), base level included - i.MNR_L1 carries 0b.MNR_L1.
			vMNRout = get_mnr_outcome()
			nRegimens = cols(vMNRout)

			mMNRdum = J(rows(idx), 0, .)
			for (r = 1; r <= nRegimens; r++) {
				mMNRdum = mMNRdum, (vMNR[idx] :== vMNRout[1, r])
			}

			// Assemble patient matrix. ECOG/RISS dummies include the base level, as in
			// sim_mnt.do / sim_mnr.do.
			mPat = (vAge[idx], vAge2[idx], vMale[idx],
					vECOG0[idx], vECOG1[idx], vECOG2[idx],
					vRISS1[idx], vRISS2[idx], vRISS3[idx],
					vSCT_L1[idx])

			mPat = mPat, mMNRdum
			mPat = mPat, vCons[idx]

			nPredictors = cols(mPat)

			// Guard: mean equation + a single scale column. A mismatch means the design here
			// has drifted from the fit, or betareg returned an unexpected layout.
			if (cols(vCoef) != nPredictors + 1) {
				errprintf("sim_mnd: design/coefficient mismatch - mPat has %g columns so %g coefficients were expected (mean equation + ln(phi)), but the coefficient vector has %g\n",
					nPredictors, nPredictors + 1, cols(vCoef))
				exit(459)
			}

			// Mean equation -> mu via the logit link
			vBeta = vCoef[1, 1..nPredictors]'
			vXB = mPat * vBeta
			vMu = 1 :/ (1 :+ exp(-vXB))

			// Scale equation -> phi via the log link (final column is ln(phi))
			phi = exp(vCoef[1, cols(vCoef)])

			// Beta(a, b) draw by inverting the CDF with ONE uniform. Beta(mu*phi, (1-mu)*phi)
			// is betareg's parameterisation. invibeta() is the inverse incomplete beta, i.e.
			// the Beta quantile function, so invibeta(a, b, u) with u ~ U(0,1) IS a Beta draw
			// - exact, and one uniform per patient. The textbook alternative (X/(X+Y) from two
			// Gammas) would need two uniforms and break the one-slot-per-event CRN layout that
			// core/rng_slots.do exists to hold.
			vA = vMu :* phi
			vB = (1 :- vMu) :* phi
			vRN = rnDraw(idx, rn_mnd())
			vOC = invibeta(vA, vB, vRN)

			// Update vector. Patients with MNT == 0 keep the missing they were initialised
			// with in mata_setup.do - they have no maintenance and so no share.
			vMNS[idx] = vOC
		}
	}
	else {
		// No model - leave vMNS missing. process_data.do bills nothing where the share is
		// missing, which is the safe direction: it reverts to no maintenance cost rather
		// than silently reverting to the old whole-gap bill.
		errprintf("sim_mnd: bL1_MND not found - maintenance duration will not be costed. Re-run prep/risk_equations.do.\n")
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnd_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
