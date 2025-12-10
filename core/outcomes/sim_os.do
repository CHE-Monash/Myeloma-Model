**********
* SIM OS 
* 
* Purpose: Calculate Overall Survival for any line of therapy
* Method: Parametric survival with BCR interacted with Line
* Outcome: Continuous time (months)
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Determine coefficient segments based on OMC
	if (OMC <= 2) { 
		segments = 0
	}
	else if (OMC == 3 | OMC == 4) {
		// L1E/L2S: need to loop through both ASCT and non-ASCT
		segments = (1, 2)
	}
	else {
		// OMC 5,6 → segment 3; OMC 7,8 → segment 4; OMC 9,10 → segment 5; etc.
		segments = floor((OMC - 1) / 2)
		if (segments > 7) segments = 7
	}

	// Loop through segments
	for (s = 1; s <= cols(segments); s++) { // cols(segments) is 1 except at L1E/L2S
		segment = segments[s]
    
		// Determine BCR column based on segment
		if (segment == 0) {
			vBCR = J(Obs, 1, 5)  // No BCR yet, use SD as placeholder
		}
		else if (segment == 1) {  // L1E/L2S Non-ASCT
			vBCR = mBCR[., 1]
		}
		else if (segment == 2) {  // L1E/L2S ASCT
			vBCR = mBCR[., 10]
		}
		else {
			vBCR = mBCR[., Line]
		}
		
		// Create BCR dummy variables
		vBCR1 = (vBCR :== 1)
		vBCR2 = (vBCR :== 2)
		vBCR3 = (vBCR :== 3)
		vBCR4 = (vBCR :== 4)
		vBCR5 = (vBCR :== 5)
		vBCR6 = (vBCR :== 6)
		
		// Calculate BCR coefficient start position
		bcrStart = 10 + segment * 6
		
		// Build coefficient column vector
		coefCols = (1, 2, 3, 5, 6, 8, 9,                           // Base effects
					bcrStart, bcrStart+1, bcrStart+2,              // BCR 1-3
					bcrStart+3, bcrStart+4, bcrStart+5,            // BCR 4-6
					58)                                            // Constant
		
		// Determine patients for this segment
		if (segment == 0) { // DN or L1S
			if (OMC == 1) { // DN
				idx = selectindex(mState[., 1] :<= OMC)
			}
			else { // L1S
				idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
			}
		}
		else if (segment == 1) { // L1E / L2S No ASCT
			idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vSCT_L1 :== 0))
		}
		else if (segment == 2) { // L1E / L2S ASCT
			idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vSCT_L1 :== 1))
		}
		else {
			// L2E+
			idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
		}
			
		// Calculate for this segment if patients exist
		if (rows(idx) > 0) {
			
			// Assemble patient matrix - reference ECOG/RISS categories NOT required as using coefCols
			pMat = (vAge[idx], vAge2[idx], vMale[idx], 
					vECOG1[idx], vECOG2[idx], 
					vRISS2[idx], vRISS3[idx],
					vBCR1[idx], vBCR2[idx], vBCR3[idx], 
					vBCR4[idx], vBCR5[idx], vBCR6[idx],
					vCons[idx])
				
			// Extract coefficients
			vCoef = bOS[1, coefCols]'
			aux =  bOS[1, cols(bOS)]
				
			// Calculate XB
			vXB = pMat * vCoef
				
			// Calculate probability of survival to current time point
			if (OMC >= 2) {
				vPR = calcSurvProb(vXB, mTSD[idx, OMC], fbOS, aux)
				vRN = runiform(rows(idx), 1) :* vPR // Conditional on survival to mTSD
			}
			else {
				vRN = runiform(rows(idx), 1)  // At diagnosis, no conditioning needed
			}
				
			// Calculate survival time
			vOC = calcSurvTime(vXB, vRN, fbOS, aux)
				
			// Update matrix
			mOS[idx, OMC] = round(vOC, 0.1)
		}
	}
}

// 	Segment	Stage			OMC		bcrCol
// 	0		DN / L1S		1 / 2	placeholder
//	1		L1E - No ASCT	3		1
//	2		L1E - ASCT		3		10
//	1		L2S - No ASCT	4		1
//	2		L2S - ASCT		4		10
//	3		L2E				5		2
//	3		L3S				6		2
//	4		L3E				7		3
//	4		L4S				8		3
//	5		L4E				9		4
//	5		L5S				10		4
//	6		L5E				11		5
//	6		L6S				12		5
//	7		L6E				13		6
//	7		L7S				14		6
//	8		L7E				15		7
//	8		L8S				16		7
//	9		L8E				17		8
//	9		L9S				18		8
//	10		L9E				19		9
