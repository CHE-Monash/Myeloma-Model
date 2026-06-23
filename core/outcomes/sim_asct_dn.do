**********
* SIM ASCT DN
* 
* Purpose: Determine ASCT eligibility at diagnosis
* Method: Logistic regression
* Outcome: Binary (0 = Not eligible, 1 = Eligible)
**********
	
mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for eligible patients
	idx = selectindex(mState[., 1] :<= OMC)
	if (rows(idx) > 0) {
	
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx], 
				vAge70[idx], vAge75[idx], 
				vCM0[idx], vCM1[idx], vCM2[idx], vCM3[idx],
				vCons[idx])

		// Extract coefficients
		nPredictors = cols(mPat)
		// Guard: design columns must equal the coefficient count (no cutpoints/ancillary here).
		if (nPredictors != cols(bDN_SCT)) {
			errprintf("sim_asct_dn: design/coefficient mismatch - mPat has %g columns but coefficient vector has %g\n", nPredictors, cols(bDN_SCT))
			exit(459)
		}
		vCoef = bDN_SCT[1, 1..nPredictors]'

		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		vPR = 1 :/ (1 :+ exp(-vXB))
		
		// Generate random numbers
		vRN = rnDraw(idx, rn_asct_dn())
			
		// Determine outcome 
		vOC = (vPR :> vRN)
		
		// Update matrix
		vSCT_DN[idx] = vOC
	}
}
