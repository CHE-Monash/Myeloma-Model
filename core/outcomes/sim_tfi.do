**********
* SIM TFI (L2+) 
*
* Purpose: Treatment-free Interval (time from LXE to L(X+1)S)
* Method: Parametric survival analysis
* Outcome: Continous time (months)
**********

mata {
    // Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
    idx = selectindex((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC))
	if (rows(idx) > 0) {

		// Select coefficients
		if (Line == 2) {
			vCoef = bL2_TFI
			dist = fbL2_TFI
			maxTFI = maxL2_TFI
		}
		else if (Line == 3) {
			vCoef = bL3_TFI
			dist = fbL3_TFI
			maxTFI = maxL3_TFI
		}
		else if (Line == 4) {
			vCoef = bL4_TFI
			dist = fbL4_TFI
			maxTFI = maxL4_TFI
		}
		else if (Line == 5) {
			vCoef = bL5_TFI
			dist = fbL5_TFI
			maxTFI = maxL5_TFI
		}
		else if (Line >= 6) {
			vCoef = bLX_TFI
			dist = fbLX_TFI
			maxTFI = maxLX_TFI
		}
		
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx],
		        vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// Add BCR
		vBCR = mBCR[idx, Line]
		vBCR_1 = (vBCR :== 1)
		vBCR_2 = (vBCR :== 2)
		vBCR_3 = (vBCR :== 3)
		vBCR_4 = (vBCR :== 4)
		vBCR_5 = (vBCR :== 5)
		vBCR_6 = (vBCR :== 6)
		mPat = (mPat, vBCR_1, vBCR_2, vBCR_3, vBCR_4, vBCR_5, vBCR_6)
	
		// Add constant
		mPat = mPat, vCons[idx]
		
		// Extract coefficients
		aux = vCoef[1, cols(vCoef)]
		nPredictors = cols(mPat)
		vCoef = vCoef[1, 1..nPredictors]'
		
		// Calculate XB and OC
		vXB = mPat * vCoef
		vRN = runiform(rows(idx), 1)
		vOC = calcSurvTime(vXB, vRN, dist, aux)
		
		// Curtail if beyond maximum observed
		vOC = rowmin((vOC, J(rows(vOC), 1, maxTFI)))
		
		// Update matrices
		mTFI[idx, Line+1] = round(vOC, 0.01) // Col 1 of mTFI is DN, so TFI_L3 goes in column 4
		mTNE[idx, OMC] = round(vOC, 0.01)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] + mTNE[idx, OMC]
	}
}
