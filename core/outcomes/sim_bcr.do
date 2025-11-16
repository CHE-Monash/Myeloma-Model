**********
*SIM BCR Vector

* Purpose: Vectorised Best Clinical Response (BCR) calculation for all lines
*          of therapy (L1-L9)
*
* Key Features:
*   - Single implementation works across all therapy lines
*   - Handles 6 categories (L1-L2): CR, nCR, VGPR, MR, PR, SD/PD
*   - Handles 3 categories (L3-L9): CR, VGPR, PR
*   - Dynamically adjusts predictors based on line number
*
* Categories:
* 6-category (L1-L2): 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD
* 3-category (L3+): 1=CR/VGPR, 3=PR/MR, 5=SD/PD
**********

mata:
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
	if (Line >= 2) {
		if (Line == 2 | Line == 3) { // 6-category
			pBCR = mBCR[., Line]
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
		else {
			// L4+ 3-category
			prevLineIdx = Line - 1
			prevBCR = mBCR[., prevLineIdx]
			pBCR_CR = (prevBCR :== 1)
			pBCR_PR = (prevBCR :== 3)
			pBCR_SD = (prevBCR :== 5)
		}
	}
	
	// Extract TXR dummies
	if (Line <= 4) {
		if (Line == 1) oVector = oL1_TXR
		if (Line == 2) oVector = oL2_TXR
		if (Line == 3) oVector = oL3_TXR
		if (Line == 4) oVector = oL4_TXR
		
		currentTX = mTXR[., Line]
		
		if (cols(oVector) >= 2) TXR_is_R2 = (currentTX :== oVector[1, 2])
		if (cols(oVector) >= 3) TXR_is_R3 = (currentTX :== oVector[1, 3])
		if (cols(oVector) >= 4) TXR_is_R4 = (currentTX :== oVector[1, 4])
		
		nRegimens = cols(oVector)
	}
	else {
		nRegimens = 0
	}
	
	// Build patient matrix
	mPat = (vAge, vAge2, vMale, vECOG0, vECOG1, vECOG2, vRISS1, vRISS2, vRISS3)
	
	// Add SCT
	if (Line == 1) mPat = mPat, vSCT_DN
	
	// Add previous BCR
	if (Line == 2 | Line == 3) {
		mPat = mPat, (pBCR_CR, pBCR_VG, pBCR_PR, pBCR_MR, pBCR_SD, pBCR_PD)
	}	
	if (Line == 2) {
		mPat = mPat, (pBCR_SCT_0, pBCR_SCT_1, pBCR_SCT_2, pBCR_SCT_3, pBCR_SCT_4)
	}
	else if (Line >= 4) {
		mPat = mPat, (pBCR_CR, pBCR_PR, pBCR_SD)
	}
		
	// Add TXR dummies
	if (Line <= 4) {
		if (Line == 1) oVector = oL1_TXR
		if (Line == 2) oVector = oL2_TXR
		if (Line == 3) oVector = oL3_TXR
		if (Line == 4) oVector = oL4_TXR
		
		currentTX = mTXR[., Line]
		
		if (cols(oVector) >= 1) mPat = mPat, (currentTX :== oVector[1, 1])
		if (cols(oVector) >= 2) mPat = mPat, (currentTX :== oVector[1, 2])
		if (cols(oVector) >= 3) mPat = mPat, (currentTX :== oVector[1, 3])
		if (cols(oVector) >= 4) mPat = mPat, (currentTX :== oVector[1, 4])
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
	XB = mPat * vCoef
	
	// Calculate probabilities
	cumProbs = calcOrdLogitProbs(XB, cutPoints)
	
	// Assign outcomes
	RN = runiform(nObs, 1)
	vOutcome = assignOrdOutcome(RN, cumProbs, categoryValues)
	
	// Update matrices
	for (i = 1; i <= nObs; i++) {
		if (mMOR[i, OMC - 1] == 0 & mState[i, 1] <= OMC + 1) {
			finalOutcome = vOutcome[i]
		}
		else if (mState[i, 1] > OMC + 1) {
			finalOutcome = mBCR[i, LX]
		}
		else {
			finalOutcome = .
		}
			
		if (mMOR[i, OMC - 1] == 0) {
			mBCR[i, LX] = finalOutcome
		}
	}
end

// Check for override file, execute if it exists
mata: st_local("current_line", strofreal(Line))
if "${Line}" == "`current_line'" {
	local override_file "${analysis_path}/outcomes/sim_bcr_override_${int}_l${line}.do"
	capture confirm file "`override_file'"
	if _rc == 0 {
		di "Exists"
		quietly do "`override_file'"
	}
}
