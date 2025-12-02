**********
* SIM TFI L1
*
* Purpose: Treatment-free Interval at Line 1 End (time from L1E to L2S)
* Method: Parametric survival analysis, split by ASCT status
* Outcome: Continuous time (months)
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
	idx = selectindex((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC))	
	if (rows(idx) > 0) {

		// GROUP 1: ASCT Patients`'
		idxASCT = idx[selectindex(vSCT_L1[idx] :== 1)]
		if (rows(idxASCT) > 0) {
			
			//Grab BCR to ASCT
			vBCR1 = (mBCR[idxASCT, 10] :== 1)
			vBCR2 = (mBCR[idxASCT, 10] :== 2)
			vBCR3 = (mBCR[idxASCT, 10] :== 3)
			vBCR4 = (mBCR[idxASCT, 10] :== 4)	
			
			// Assemble patient matrix (ASCT patients cannot have BCR = 5 or 6)
			mPat_ASCT = (vAge[idxASCT], vAge2[idxASCT], vMale[idxASCT], 
						 vECOG0[idxASCT], vECOG1[idxASCT], vECOG2[idxASCT], 
			             vRISS1[idxASCT], vRISS2[idxASCT], vRISS3[idxASCT], 
						 vMNT[idxASCT], 
			             vBCR1, vBCR2, vBCR3, vBCR4, 
						 vCons[idxASCT])
			
			// Extract coefficients for ASCT
			nPredictors = cols(mPat_ASCT)
			coef_ASCT = bL1_TFI_ASCT[1, 1..nPredictors]' 
			aux_ASCT = bL1_TFI_ASCT[1, cols(bL1_TFI_ASCT)]
			
			// Calculate XB
			vXB_ASCT = mPat_ASCT * coef_ASCT
			
			// Calculate outcome (survival time)
			vRN_ASCT = runiform(rows(idxASCT), 1)
			vOC[idxASCT] = calcSurvTime(vXB_ASCT, vRN_ASCT, fbL1_TFI_ASCT, aux_ASCT)
			
			// Curtail if beyond maximum observed
			vOC[idxASCT] = rowmin((vOC[idxASCT], J(rows(idxASCT), 1, maxL1_TFI_ASCT)))
		}
		
		// GROUP 2: NoASCT Patients
		idxNoASCT = idx[selectindex(vSCT_L1[idx] :== 0)]
		if (rows(idxNoASCT) > 0) {

			//Grab BCR to L1
			vBCR1 = (mBCR[idxNoASCT, 1] :== 1)
			vBCR2 = (mBCR[idxNoASCT, 1] :== 2)
			vBCR3 = (mBCR[idxNoASCT, 1] :== 3)
			vBCR4 = (mBCR[idxNoASCT, 1] :== 4)
			vBCR5 = (mBCR[idxNoASCT, 1] :== 5)	
			vBCR6 = (mBCR[idxNoASCT, 1] :== 6)	
		
			// Assemble patient matrix (NoASCT includes BCR 5 and 6)
			mPat_NoASCT = (vAge[idxNoASCT], vAge2[idxNoASCT], vMale[idxNoASCT], 
						   vECOG0[idxNoASCT], vECOG1[idxNoASCT], vECOG2[idxNoASCT],
			               vRISS1[idxNoASCT], vRISS2[idxNoASCT], vRISS3[idxNoASCT], 
						   vMNT[idxNoASCT],
			               vBCR1, vBCR2, vBCR3, vBCR4, vBCR5, vBCR6, 
						   vCons[idxNoASCT])
			
			// Extract coefficients for NoASCT
			nPredictors = cols(mPat_NoASCT)
			vCoef_NoASCT = bL1_TFI_NoASCT[1, 1..nPredictors]'
			aux_NoASCT = bL1_TFI_NoASCT[1, cols(bL1_TFI_NoASCT)]
			
			// Calculate XB
			vXB_NoASCT = mPat_NoASCT * vCoef_NoASCT
			
			// Calculate outcome (survival time)
			vRN_NoASCT = runiform(rows(idxNoASCT), 1)
			vOC[idxNoASCT] = calcSurvTime(vXB_NoASCT, vRN_NoASCT, fbL1_TFI_NoASCT, aux_NoASCT)
			
			// Curtail if beyond maximum observed
			vOC[idxNoASCT] = rowmin((vOC[idxNoASCT], J(rows(idxNoASCT), 1, maxL1_TFI_NoASCT)))

		}
		
		// Update matrices
		mTFI[idx, 2] = round(vOC[idx], 0.1)
		mTNE[idx, OMC] = round(vOC[idx], 0.1)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] + mTNE[idx, OMC]
	}
}
