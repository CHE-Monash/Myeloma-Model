**********
* SIM BCR ASCT

* Purpose: Determine BCR to SCT
* Mehod: Multinomial logit
* Outcome: Categorical 1 to 4, no 5/6 
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive, eligible and valid SCT patients
	idx = selectindex((mMOR[., OMC-1] :== 0) :& 
	                  (mState[., 1] :<= OMC) :&
	                  (vSCT_L1 :== 1))
	if (rows(idx) > 0) {
		
		// Extract L1 BCR
		vBCR_L1 = mBCR[., 1]
		
		// Create BCR dummy variables
		pBCR_1 = (vBCR_L1 :== 1)
		pBCR_2 = (vBCR_L1 :== 2)
		pBCR_3 = (vBCR_L1 :== 3)
		pBCR_4 = (vBCR_L1 :== 4)
		pBCR_5 = (vBCR_L1 :== 5)
		
		// Assemble patient matrix - no BCR = 6
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
		        vECOG0[idx], vECOG1[idx], vECOG2[idx], 
		        vRISS1[idx], vRISS2[idx], vRISS3[idx],
		        pBCR_1[idx], pBCR_2[idx], pBCR_3[idx], pBCR_4[idx], pBCR_5[idx])
		
		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef = bSCT_BCR[1, 1..nPredictors]'
		
		// Extract cut points
		nCutPoints = 3
		cutPointIndices = (cols(bSCT_BCR) - nCutPoints + 1)..cols(bSCT_BCR)
		cutPoints = bSCT_BCR[1, cutPointIndices]
		
		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		cumProbs = calcOrdLogitProbs(vXB, cutPoints)
				
		// Assign outcomes
		vRN = runiform(rows(idx), 1)
		categoryValues = (1, 2, 3, 4)
		vOC = assignOrdOutcome(vRN, cumProbs, categoryValues)
		
		// Update matrix
		mBCR[idx, 10] = vOC
	}
}
