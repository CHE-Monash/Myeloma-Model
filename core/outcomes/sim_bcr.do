**********
* SIM BCR 
*
* Purpose: Determine Best Clinical Response (BCR)
* Method: Ordered logit
* Outcome:
* 6-category (L1,L2): 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD
* 3-category (L3+): 1=CR/VGPR, 3=PR/MR, 5=SD/PD
**********

mata {
	// Determine structure
	if (Line <= 2) {
		nCategories = 6
		categoryValues = (1, 2, 3, 4, 5, 6)
		nCutPoints = 5
	}
	else {
		nCategories = 3
		categoryValues = (1, 3, 5)
		nCutPoints = 2
	}
	
	// Extract previous BCR
	if (Line == 2 | Line == 3) { // L1, L2 6-category
		pBCR = mBCR[., Line-1]
		pBCR_CR = (pBCR :== 1)
		pBCR_VG = (pBCR :== 2)
		pBCR_PR = (pBCR :== 3)
		pBCR_MR = (pBCR :== 4)
		pBCR_SD = (pBCR :== 5)
		pBCR_PD = (pBCR :== 6)
	}
	if (Line == 2) { // Also need BCR_SCT
		BCR_SCT = mBCR[., 10]
		pBCR_SCT_0 = (BCR_SCT :== 0) // No ASCT patients
		pBCR_SCT_1 = (BCR_SCT :== 1)
		pBCR_SCT_2 = (BCR_SCT :== 2)
		pBCR_SCT_3 = (BCR_SCT :== 3)
		pBCR_SCT_4 = (BCR_SCT :== 4)
	}
	else if (Line >= 4) { // L3+ 3-category
		pBCR = mBCR[., Line-1]
		pBCR_CR = (pBCR :== 1)
		pBCR_PR = (pBCR :== 3)
		pBCR_SD = (pBCR :== 5)
	}
	
	// Filter for alive and eligble
	idx = selectindex((mMOR[., OMC - 1] :== 0) :& (mState[., 1] :<= OMC + 1))
	if (rows(idx) > 0) {
	
		// Build patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx],
				vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// Add SCT
		if (Line == 1) mPat = mPat, vSCT_DN[idx]
		
		// Add previous BCR
		if (Line == 2 | Line == 3) {
			mPat = mPat, (pBCR_CR[idx], pBCR_VG[idx], pBCR_PR[idx], 
						  pBCR_MR[idx], pBCR_SD[idx], pBCR_PD[idx])
		}	
		if (Line == 2) {
			mPat = mPat, (pBCR_SCT_0[idx], pBCR_SCT_1[idx], pBCR_SCT_2[idx], 
						  pBCR_SCT_3[idx], pBCR_SCT_4[idx])
		}
		else if (Line >= 4) {
			mPat = mPat, (pBCR_CR[idx], pBCR_PR[idx], pBCR_SD[idx])
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
		if (Line >= 5) vCoef_full = bLX_BCR
		
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
	local override_file "${analysis_path}/outcomes/sim_bcr_override.do"
	capture confirm file "`override_file'"
	if _rc == 0 {
		qui do `override_file'
	}
}
