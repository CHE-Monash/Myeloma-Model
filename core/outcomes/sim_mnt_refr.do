**********
* Monash Myeloma Model - Sim MNT-refr (maintenance len-refractory)
*
* Purpose: set LenRefr_Mnt_in (vLenRefr_Mnt_in) - whether the patient becomes refractory to their L1
*          lenalidomide MAINTENANCE. A single per-patient event, drawn at L1E from the L1_MNTREFR
*          logit and constant from L2 on. Consumed by OS_L2..L4 (docs/refractory.md 4.4).
*
*          Predictor is the TAIL = TFI_L1 - TTM - MND_L1 (the gap between maintenance ending and L2
*          starting), the IMWG 60-day clock: a short tail is progression on/near lenalidomide
*          (refractory), a long tail a later relapse (sensitive). All three inputs exist at L1E:
*            mTFI[.,2]  drawn TFI_L1 (MONTHS)      - sim_tfi_l1
*            vMND       drawn maintenance duration (MONTHS) - sim_mnd (before process_data's cap)
*            TTM        0.49mo (no-ASCT) / 4.81mo (ASCT), the transplant-keyed constants
*          so tail = mTFI[.,2] - TTM - vMND, in months (no unit conversion - both already months).
*
* Fired by: core/simulation_engine.do at OMC 3 (L1E), AFTER sim_mnd. Only for lenalidomide
*           maintenance patients (vMNT == 1 & vMNR == 1); everyone else stays 0.
* Reads:    vLenRefr_Mnt_in, mTFI, vMND, vSCT_L1, vMNT, vMNR, bL1_MNTREFR, baseline vectors, rn_mntrefr().
* Writes:   vLenRefr_Mnt_in (0/1).
**********

mata {
	if (mntrefr_model_exists()) {

		bMR = get_mntrefr_coef()   // 1 x 13 logit coefficients (factor-expanded, base = 0)

		// Alive, reached L1E, and on lenalidomide maintenance
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1) :& (vMNR :== 1))
		if (rows(idx) > 0) {

			// Tail (months) and ln(TFI months). TTM: 0.49 no-ASCT, 4.81 ASCT.
			vSCTa = vSCT_L1[idx]
			vTTM  = 4.81 :* vSCTa :+ 0.49 :* (1 :- vSCTa)
			vTail = mTFI[idx, 2] :- vTTM :- vMND[idx]
			vLn   = ln(rowmax((mTFI[idx, 2], J(rows(idx), 1, 0.01))))   // guard ln(0)

			// Design: Age Age2 Male i.ECOGcc i.RISS SCT ln(TFI) tail _cons - order MUST match the fit
			mPat = (vAge[idx], vAge2[idx], vMale[idx],
			        vECOG0[idx], vECOG1[idx], vECOG2[idx],
			        vRISS1[idx], vRISS2[idx], vRISS3[idx],
			        vSCT_L1[idx], vLn, vTail,
			        vCons[idx])

			if (cols(mPat) != cols(bMR)) {
				errprintf("sim_mnt_refr: design/coefficient mismatch - mPat has %g columns but coefficient vector has %g\n", cols(mPat), cols(bMR))
				exit(459)
			}

			vXB = mPat * bMR'
			vPR = 1 :/ (1 :+ exp(-vXB))

			// One CRN draw per patient; refractory when p > u
			vRN = rnDraw(idx, rn_mntrefr())
			vLenRefr_Mnt_in[idx] = (vPR :> vRN)
		}
	}
}
