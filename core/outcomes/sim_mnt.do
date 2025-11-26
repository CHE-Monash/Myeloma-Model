**********
* SIM MNT

* Purpose: Determine receipt of Maintenance Therapy
* Method: Logistic regression
* Outcome: Binary (1 = MNT / 0 = No MNT)
**********

mata {
	// Filter for alive and eligible patients
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
	if (rows(idx) > 0) {
		
		// Extract current treatment regimen
		currentTX = mTXR[., Line]
		
		// Determine number of treatment regimen dummies
		if (Line == 1) vTXR = oL1_TXR
		else if (Line == 2) vTXR = oL2_TXR
		else if (Line == 3) vTXR = oL3_TXR
		else if (Line == 4) vTXR = oL4_TXR
		
		nRegimens = cols(vTXR)
		
		// Create treatment regimen dummies
		if (nRegimens >= 2) TXR_is_R2 = (currentTX :== vTXR[1, 2])
		if (nRegimens >= 3) TXR_is_R3 = (currentTX :== vTXR[1, 3])
		if (nRegimens >= 4) TXR_is_R4 = (currentTX :== vTXR[1, 4])
		
		// Extract L1 BCR
		vBCR_L1 = mBCR[., 1]
		
		// Create BCR dummy variables (reference = 1)
		pBCR_1 = (vBCR_L1 :== 1)
		pBCR_2 = (vBCR_L1 :== 2)
		pBCR_3 = (vBCR_L1 :== 3)
		pBCR_4 = (vBCR_L1 :== 4)
		pBCR_5 = (vBCR_L1 :== 5)
		pBCR_6 = (vBCR_L1 :== 6)
		
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
		        vECOG0[idx], vECOG1[idx], vECOG2[idx], 
		        vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// Add treatment regimen dummies
		if (nRegimens >= 2) mPat = mPat, TXR_is_R2[idx]
		if (nRegimens >= 3) mPat = mPat, TXR_is_R3[idx]
		if (nRegimens >= 4) mPat = mPat, TXR_is_R4[idx]
		
		// Add SCT
		mPat = mPat, vSCT_L1[idx]
		
		// Add BCR dummies
		mPat = mPat, (pBCR_1[idx], pBCR_2[idx], pBCR_3[idx], 
		              pBCR_4[idx], pBCR_5[idx], pBCR_6[idx])
		
		// Add constant
		mPat = mPat, vCons[idx]
		
		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef = bMNT[1, 1..nPredictors]'
		
		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probability
		vPR = 1 :/ (1 :+ exp(-vXB))
				
		// Determine outcome
		vRN = runiform(rows(idx), 1)
		vOC = (vPR :> vRN) :* 1
		
		// Update vector
		vMNT[idx] = vOC
	}
}
