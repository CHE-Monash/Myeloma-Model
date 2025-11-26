**********
* SIM TFI DN
* 
* Purpose: Treatment-free Interval at Diagnosis (Time from DN to L1S)
* Method: Parametric survival analysis
* Outcome: Continuous time (months)
**********
	
mata {
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
		aux = bDN_TFI[1, cols(bDN_TFI)]
			
		// Calculate XB
		vXB = mPat * vCoef
			
		// Calculate outcome 
		vRN = runiform(rows(idx), 1)
		vOC = calcSurvTime(vXB, vRN, fbDN_TFI, aux)
		
		// Update matrices	
		mTFI[idx, 1] = round(vOC, 0.1)  
		mTNE[idx, OMC] = round(vOC, 0.1)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] + mTNE[idx, OMC]
	}	
}
