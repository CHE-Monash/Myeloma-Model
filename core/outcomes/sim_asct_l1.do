**********
* SIM ASCT L1
* 
* Purpose: Determine receipt of ASCT at Line 1 End
* Method: Logistic regression with outcome filters: TXR != 7 AND BCR != 6
* Outcome: Binary (0 = No ASCT, 1 = ASCT)
**********
	
mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
	if (rows(idx) > 0) {
		
		// Patient vectors
		vBCR_1 = (mBCR[idx, 1] :== 1)
		vBCR_2 = (mBCR[idx, 1] :== 2)
		vBCR_3 = (mBCR[idx, 1] :== 3)
		vBCR_4 = (mBCR[idx, 1] :== 4)
		vBCR_5 = (mBCR[idx, 1] :== 5)
		
		// Assemble patient matrix - no BCR == 6
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
		        vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx], 
				vBCR_1, vBCR_2, vBCR_3, vBCR_4, vBCR_5,
				vAge70[idx], vAge75[idx], 
				vCM0[idx], vCM1[idx], vCM2[idx], vCM3[idx], 
				vCons[idx])

		// Extract coefficients
		nPredictors = cols(mPat)
		// Guard: design columns must equal the coefficient count (no cutpoints/ancillary here).
		if (nPredictors != cols(bL1_SCT)) {
			errprintf("sim_asct_l1: design/coefficient mismatch - mPat has %g columns but coefficient vector has %g\n", nPredictors, cols(bL1_SCT))
			exit(459)
		}
		vCoef = bL1_SCT[1, 1..nPredictors]'

		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		vPR = 1 :/ (1 :+ exp(-vXB))
		
		// Generate random numbers
		vRN = rnDraw(idx, rn_asct_l1())
		
		// Extract TXR from matrix
		vTXR_L1_e = mTXR[idx, 1]
			
		// Determine outcome: PR > RN AND TXR != 7 AND BCR != 6
		vOC = (vPR :> vRN) :& (mTXR[idx, Line] :!= 7) :& (mBCR[idx, Line] :!= 6)		
		
		// Update matrix
		vSCT_L1[idx] = vOC
	}
}
