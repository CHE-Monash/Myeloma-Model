**********
* SIM BCR 
*
* Purpose: Determine Best Clinical Response (BCR)
* Method: Ordered logit
* Outcome: 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Determine structure
	nCategories = 6
	categoryValues = (1, 2, 3, 4, 5, 6)
	nCutPoints = 5
	
	// Extract previous BCR
	if (Line >= 2) {
		vBCR = mBCR[., Line]
		vBCR_1 = (vBCR :== 1)
		vBCR_2 = (vBCR :== 2)
		vBCR_3 = (vBCR :== 3)
		vBCR_4 = (vBCR :== 4)
		vBCR_5 = (vBCR :== 5)
		vBCR_6 = (vBCR :== 6)
	}
	if (Line == 2) { // Need BCR_SCT
		vBCR_SCT = mBCR[., 10]
		vBCR_SCT_0 = (vBCR_SCT :== 0) // No ASCT patients
		vBCR_SCT_1 = (vBCR_SCT :== 1)
		vBCR_SCT_2 = (vBCR_SCT :== 2)
		vBCR_SCT_3 = (vBCR_SCT :== 3)
		vBCR_SCT_4 = (vBCR_SCT :== 4)
	}
	
	// Filter for alive and eligble
	idx = selectindex((mMOR[., OMC - 1] :== 0) :& (mState[., 1] :<= OMC))
	if (rows(idx) > 0) {
	
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx],
				vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// BCR_L1 - Add SCT
		if (Line == 1) mPat = mPat, vSCT_DN[idx]
		
		// BCR_L2 onwards - Add previous BCR
		if (Line >= 2 ) {
			mPat = mPat, (vBCR_1[idx], vBCR_2[idx], vBCR_3[idx], 
						  vBCR_4[idx], vBCR_5[idx], vBCR_6[idx])
		}
		// BCR_L2 only - Add BCR_SCT
		if (Line == 2) {
			mPat = mPat, (vBCR_SCT_0[idx], vBCR_SCT_1[idx], vBCR_SCT_2[idx], 
						  vBCR_SCT_3[idx], vBCR_SCT_4[idx])
		}
		
		// Add TXR dummies
		if (Line <= 4) {
			if (Line == 1) vTXR = oL1_TXR
			if (Line == 2) vTXR = oL2_TXR
			if (Line == 3) vTXR = oL3_TXR
			if (Line == 4) vTXR = oL4_TXR
			
			currentTX = mTXR[idx, Line]
			
			if (cols(vTXR) >= 1) mPat = mPat, (currentTX :== vTXR[1, 1])
			if (cols(vTXR) >= 2) mPat = mPat, (currentTX :== vTXR[1, 2])
			if (cols(vTXR) >= 3) mPat = mPat, (currentTX :== vTXR[1, 3])
			if (cols(vTXR) >= 4) mPat = mPat, (currentTX :== vTXR[1, 4])
		}
		
		nPredictors = cols(mPat)
		
		// Extract coefficients
		if (Line == 1) vCoef_full = bL1_BCR
		if (Line == 2) vCoef_full = bL2_BCR
		if (Line == 3) vCoef_full = bL3_BCR
		if (Line == 4) vCoef_full = bL4_BCR
		if (Line == 5) vCoef_full = bL5_BCR
		if (Line == 6) vCoef_full = bL6_BCR
		if (Line == 7) vCoef_full = bL7_BCR
		if (Line >= 8) vCoef_full = bLX_BCR
		
		vCoef = vCoef_full[1, 1..nPredictors]'
		cutPointIndices = (cols(vCoef_full) - nCutPoints + 1)..cols(vCoef_full)
		cutPoints = vCoef_full[1, cutPointIndices]
				
		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		cumProbs = calcOrdLogitProbs(vXB, cutPoints)
		
		// Assign outcomes
		vRN = runiform(rows(idx), 1)
		vOC = assignOrdOutcome(vRN, cumProbs, categoryValues)
	
		// Update matrix
		mBCR[idx, Line] = vOC
	}	
}

// Check for override file, execute if it exists
mata: st_local("current_line", strofreal(Line))
if `current_line' == ${line} {
	local override_file "${outcomes_path}/sim_bcr_override.do"
	capture confirm file "`override_file'"
	if _rc == 0 {
		do `override_file'
	}
}
