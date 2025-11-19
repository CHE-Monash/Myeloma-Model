**********
* SIM TXD (L2+)
*
* Purpose: Calculate treatment duration for L2 onwards
* Method: Parametric survival analysis
* Outcome: Continuous time (months)
**********
mata {
	// Filter for alive and eligible
	idx = selectindex((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1))
	if (rows(idx) > 0) {
		
		// Determine model stage and BCR structure
		if (Line == 1) {
			vCoef = bL2_TXD
			dist = fbL2_TXD
			vTXR = oL2_TXR
			BCR_cat = 6 
		}
		else if (Line == 2) {
			vCoef = bL3_TXD
			dist = fbL3_TXD
			vTXR = oL3_TXR
			BCR_cat = 6 
		}
		else if (Line == 3) {
			vCoef = bL4_TXD
			dist = fbL4_TXD
			vTXR = oL4_TXR
			BCR_cat = 3  
		}
		else if (Line >= 4) {
			vCoef = bLX_TXD
			dist = fbLX_TXD
			vTXR = J(1, 0, .)
			BCR_cat = 3  
		}
		
		nTXR = cols(vTXR)
		
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// Add treatment regimen dummies
		currentTX = mTXR[idx, Line]
		if (nTXR >= 1) {
			vTXR1 = (currentTX :== vTXR[1,1])
			mPat = (mPat, vTXR1)
		}
		if (nTXR >= 2) {
			vTXR2 = (currentTX :== vTXR[1,2])
			mPat = (mPat, vTXR2)
		}
		if (nTXR >= 3) {
			vTXR3 = (currentTX :== vTXR[1,3])
			mPat = (mPat, vTXR3)
		}
		
		// Add previous BCR
		if (BCR_cat == 6) {
			prevBCR = mBCR[idx, Line]
			vBCR_CR = (prevBCR :== 1)
			vBCR_VG = (prevBCR :== 2)
			vBCR_PR = (prevBCR :== 3)
			vBCR_MR = (prevBCR :== 4)
			vBCR_SD = (prevBCR :== 5)
			vBCR_PD = (prevBCR :== 6)
			mPat = (mPat, vBCR_CR, vBCR_VG, vBCR_PR, vBCR_MR, vBCR_SD, vBCR_PD)
		}
		else if (BCR_cat == 3) {
			prevBCR = mBCR[idx, Line]
			vBCR_CR = (prevBCR :== 1)
			vBCR_PR = (prevBCR :== 3)
			vBCR_SD = (prevBCR :== 5)
			mPat = (mPat, vBCR_CR, vBCR_PR, vBCR_SD)
		}
		
		// Add constant
		mPat = (mPat, vCons[idx])
		
		// Extract coefficients
		aux = vCoef[1, cols(vCoef)]
		nPredictors = cols(mPat)
		vCoef = vCoef[1, 1..nPredictors]'
		
		// Calculate XB and duration
		vXB = mPat * vCoef
		vRN = runiform(rows(idx), 1)
		vOC[idx] = calcSurvTime(vXB, vRN, dist, aux)
		
		// Store outcomes
		mTXD[idx, LX+1] = round(vOC[idx], 0.1)
		mTNE[idx, OMC] = round(vOC[idx], 0.1)
		mTSD[idx, OMC+1] = mTSD[idx, OMC-1] :+ mTNE[idx, OMC]
	}
}
