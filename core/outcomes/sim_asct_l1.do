**********
* SIM SCT L1
* 
* Purpose: Determine receipt of ASCT at Line 1 End
* Method: Logistic regression with outcome filters: TXR != 7 AND BCR != 6
* Outcome: Binary (0 = No ASCT, 1 = ASCT)
**********
	
mata {
	// Filter for alive and eligible
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
	if (rows(idx) > 0) {
		
		// Patient vectors
		vBCR1 = (mBCR[idx, Line] :== 1)
		vBCR2 = (mBCR[idx, Line] :== 2)
		vBCR3 = (mBCR[idx, Line] :== 3)
		vBCR4 = (mBCR[idx, Line] :== 4)
		vBCR5 = (mBCR[idx, Line] :== 5)
		
		// Assemble patient matrix - no BCR == 6
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
		        vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx], 
				vBCR1[idx], vBCR2[idx], vBCR3[idx], vBCR4[idx], vBCR5[idx],
				vAge[idx], vAge75[idx], 
				vCMc0[idx], vCMc1[idx], vCMc2[idx], vCMc3[idx], 
				vCons[idx])

		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef = bL1_SCT[1, 1..nPredictors]'

		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		vPR = 1 :/ (1 :+ exp(-vXB))
		
		// Generate random numbers
		vRN = runiform(rows(idx), 1)
		
		// Extract TXR from matrix
		vTXR_L1_e = mTXR[idx, 1]
			
		// Determine outcome: PR > RN AND TXR != 7 AND BCR != 6
		vOC = (vPR :> vRN) :& (mTXR[idx, Line] :!= 7) :& (mBCR[idx, Line] :!= 6)		
		
		// Update matrix
		vSCT_L1[idx] = vOC
	}
}
