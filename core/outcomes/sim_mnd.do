**********
* Monash Myeloma Model - Sim MND
*
* Purpose: Draw the L1 maintenance duration among patients receiving maintenance (MNT == 1),
*          as a SHARE OF THE WINDOW the maintenance can occupy. The maintenance analogue of
*          sim_txd_l1.do.
*
* Notes:   THIS DRAWS A SHARE OF THE WINDOW, NOT A DURATION AND NOT A SHARE OF THE GAP.
*
*          W = TFI_L1 - TTM(sct)  is the window: the gap less the time from L1 ending to
*          maintenance starting. TTM is two constants keyed on transplant (L1_MND_TTM0 /
*          L1_MND_TTM1, fitted in prep/risk_equations.do and carried in the coefficient file),
*          because TTM is a function of transplant and essentially nothing else - ~0.5 months
*          without ASCT against ~4.8 with it.
*
*          Not a duration: TFI_L1 CONTAINS the duration, so a duration drawn independently
*          overshoots the gap 42.5% of the time against an observed 0.7%, and every overshoot
*          is truncated back to the whole gap - which is the defect being fixed.
*          Not a share of the gap: that rises with gap length for lenalidomide (0.564 -> 0.831)
*          purely as arithmetic, since maintenance runs to progression and the share is
*          1 - TTM/gap. Share of the WINDOW is flat (0.883 -> 0.906).
*          See docs/refractory.md 7.4.
*
* ORDER:   this runs AFTER sim_tfi_l1.do (it needs mTFI[., 2]) and BEFORE sim_mort.do. The
*          gap it reads is therefore the DRAWN one, not yet curtailed at death - deliberate,
*          because the share is a property of the patient's intended gap. core/process_data.do
*          forms MND_L1 = MNS_L1 * (realised TFI_L1 - TTM), so a patient who dies mid-gap has
*          their maintenance curtailed proportionally rather than being billed for a gap they
*          never lived.
*
*          Consumed by cost_tx_mnt. The TTM offset is what makes the maintenance END date
*          computable rather than just its length. It stands on the costing evidence (form C beat
*          billing a share of the whole gap), and the end date is also what any future
*          maintenance refractory rule has to key on (docs/refractory.md 4.4).
*
*          Must match the fit in prep/risk_equations.do:
*              betareg MNS_L1 Age Age2 Male i.ECOGcc i.RISS SCT c.MND_lnW##i.MNR_L1
*          betareg returns TWO equations, mean (logit link) then scale (log link). save_coefs
*          stores both, so the FINAL column of the coefficient vector is ln(phi) - the same
*          "aux is the last column" idiom the parametric survival models use.
**********

mata {
	if (mnd_model_exists()) {

		vCoef = get_mnd_coef()

		// Filter for alive, eligible AND receiving maintenance. Same population sim_mnr drew for.
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idx) > 0) {

			// The window. TFI_L1 is mTFI column 2, in months, drawn by sim_tfi_l1 just above.
			// TTM is keyed on transplant, so no-ASCT patients get ~0.5 months and ASCT ~4.8.
			vTTM = J(rows(idx), 1, L1_MND_TTM0)
			vTTM = vTTM :+ (vSCT_L1[idx] :* (L1_MND_TTM1 - L1_MND_TTM0))
			vW = mTFI[idx, 2] :- vTTM

			// A constant TTM can exceed a short gap, which would give a non-positive window and
			// an undefined ln(). Floor it: such a patient had no room for maintenance, and
			// process_data.do floors the billed duration at 0 for the same reason.
			vW = rowmax((vW, J(rows(idx), 1, 0.01)))
			vLnW = ln(vW)

			// Maintenance regimen dummies, in the SAME order gen_mnr / betareg saw them
			// (ascending code, base level included - i.MNR_L1 carries 0b.MNR_L1).
			vMNRout = get_mnr_outcome()
			nRegimens = cols(vMNRout)
			mMNRdum = J(rows(idx), 0, .)
			for (r = 1; r <= nRegimens; r++) {
				mMNRdum = mMNRdum, (vMNR[idx] :== vMNRout[1, r])
			}

			// Assemble in the fit's order: Age Age2 Male i.ECOGcc i.RISS SCT lnW i.MNR_L1
			// i.MNR_L1#c.lnW _cons. ECOG/RISS dummies include the base level, as in sim_mnt.do.
			mPat = (vAge[idx], vAge2[idx], vMale[idx],
					vECOG0[idx], vECOG1[idx], vECOG2[idx],
					vRISS1[idx], vRISS2[idx], vRISS3[idx],
					vSCT_L1[idx], vLnW)

			mPat = mPat, mMNRdum                       // i.MNR_L1
			mPat = mPat, (mMNRdum :* vLnW)             // i.MNR_L1#c.MND_lnW
			mPat = mPat, vCons[idx]

			nPredictors = cols(mPat)

			// Guard: mean equation + a single scale column.
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

			// Beta(mu*phi, (1-mu)*phi) by inverting the CDF with ONE uniform. invibeta() is the
			// Beta quantile function, so invibeta(a, b, u) with u ~ U(0,1) IS a Beta draw -
			// exact, and one uniform per patient, which is what the CRN slot layout requires.
			// THE DRAW MATTERS. phi is small (~1.4), i.e. the share is widely spread, so a
			// conditional MEAN would hand every patient much the same share and collapse that
			// spread. betareg is used rather than fracreg precisely because it returns the
			// precision as well as the mean, which is what makes the draw possible (4.4).
			vA = vMu :* phi
			vB = (1 :- vMu) :* phi
			vRN = rnDraw(idx, rn_mnd())
			vOC = invibeta(vA, vB, vRN)

			// Update. MNT == 0 patients keep the missing set in mata_setup.do.
			vMNS[idx] = vOC
		}
	}
	else {
		// No model - leave vMNS missing. process_data.do bills nothing where the share is
		// missing, which is the safe direction: no maintenance cost beats silently reverting to
		// the old whole-gap bill.
		errprintf("sim_mnd: bL1_MND not found - maintenance duration will not be costed. Re-run prep/risk_equations.do.\n")
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnd_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
