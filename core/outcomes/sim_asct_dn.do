**********
* SIM SCT DN
* 
* Purpose: Determine ASCT eligibility at diagnosis
* Method: Logistic regression
* Outcome: Binary (0 = Not eligible, 1 = Eligible)
**********
	
mata {
	// Filter for eligible patients
	idx = selectindex(mState[., 1] :<= OMC + 1)
	if (rows(idx) > 0) {
	
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx], 
				vAge70[idx], vAge75[idx], 
				vCMc0[idx], vCMc1[idx], vCMc2[idx], vCMc3[idx],
				vCons[idx])

		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef = bDN_SCT[1, 1..nPredictors]'

		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		vPR = 1 :/ (1 :+ exp(-vXB))
		
		// Generate random numbers
		vRN = runiform(rows(idx), 1)
			
		// Determine outcome 
		vOC = (vPR :> vRN)
		
		// Update matrices 
		vSCT_DN[idx] = vOC
	}
}
