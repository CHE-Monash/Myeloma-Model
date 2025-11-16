**********
* SIM TFI DN - Vectorised Implementation
* 
* Purpose: Calculate Treatment-free Interval at Diagnosis (time from diagnosis to Line 1 start)
* Method: Parametric survival analysis
* Outcome: Continuous survival time (in days)
**********
	
mata {
	// Initialize outcome
	vOC = J(st_nobs(), 1, .)
		
	// Filter for eligible
	idx = selectindex(mState[., 1] :<= OMC + 1)
	if (rows(idx) > 0) {
		
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx], 
				vSCT_DN[idx], vCons[idx])
		
		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef= bDN_TFI[1,1..nPredictors]'
		dist = fbDN_TFI
		aux = bDN_TFI[1, cols(bDN_TFI)]
			
		// Calculate XB
		vXB = mPat * vCoef
			
		// Generate random numbers
		vRN = runiform(rows(idx), 1)
			
		// Calculate outcome (survival time)
		vOC[idx] = calcSurvTime(vXB, vRN, dist, aux)
	}	
		
	// Update matrices	
	mTFI[., LX+1] = round(vOC, 0.1)                                
	mTNE[., OMC] = round(vOC, 0.1)
	mTSD[., OMC+1] = mTSD[., OMC] + mTNE[., OMC]
}
