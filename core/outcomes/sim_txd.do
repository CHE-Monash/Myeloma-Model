**********
* SIM TXD
*
* Purpose: Calculate treatment duration for L2 onwards
* Method: Parametric survival analysis
* Outcome: Continuous time (months)
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
	idx = selectindex((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC))
	if (rows(idx) > 0) {
		
		// Determine model stage and BCR structure
		if (Line == 1) {
			vCoef = bL2_TXD
			dist = fbL2_TXD
			vTXR = oL2_TXR
			maxTXD = maxL2_TXD
		}
		else if (Line == 2) {
			vCoef = bL3_TXD
			dist = fbL3_TXD
			vTXR = oL3_TXR
			maxTXD = maxL3_TXD			
		}
		else if (Line == 3) {
			vCoef = bL4_TXD
			dist = fbL4_TXD
			vTXR = oL4_TXR
			maxTXD = maxL4_TXD			
		}
		else if (Line == 4) {
			vCoef = bL5_TXD
			dist = fbL4_TXD
			vTXR = J(1, 0, .)
			maxTXD = maxL5_TXD			
		}
		else if (Line == 5) {
			vCoef = bL6_TXD
			dist = fbL6_TXD
			vTXR = J(1, 0, .)
			maxTXD = maxL6_TXD			
		}
		else if (Line == 6) {
			vCoef = bL7_TXD
			dist = fbL7_TXD
			vTXR = J(1, 0, .)
			maxTXD = maxL7_TXD			
		}
		else if (Line >= 7) {
			vCoef = bLX_TXD
			dist = fbLX_TXD
			vTXR = J(1, 0, .)
			maxTXD = maxLX_TXD			
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
		mPat = (mPat, vCons[idx])
		
		// Extract coefficients
		aux = vCoef[1, cols(vCoef)]
		nPredictors = cols(mPat)
		vCoef = vCoef[1, 1..nPredictors]'
		
		// Calculate XB and OC
		vXB = mPat * vCoef
		vRN = runiform(rows(idx), 1)
		vOC = calcSurvTime(vXB, vRN, dist, aux)
		
		// Curtail if beyond maximum observed
		vOC = rowmin((vOC, J(rows(vOC), 1, maxTXD)))

		// Update matrices
		mTXD[idx, Line] = round(vOC, 0.01)
		mTNE[idx, OMC] = round(vOC, 0.01)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] :+ mTNE[idx, OMC]
	}
}
