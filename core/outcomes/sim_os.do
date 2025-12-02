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
	
	// Determine which column of mBCR to use based on OMC
	if (OMC == 1 | OMC == 2) {
		currentBCR = J(Obs, 1, 5)
		bcrCol = 0
	}
	else if (mod(OMC, 2) == 0) {
		bcrCol = floor((OMC - 2) / 2)
		currentBCR = mBCR[., bcrCol]
	}
	else {
		bcrCol = floor((OMC - 1) / 2)
		currentBCR = mBCR[., bcrCol]
	}
    
	// Create BCR dummy variables
	vBCR1 = (currentBCR :== 1)
	vBCR2 = (currentBCR :== 2)
	vBCR3 = (currentBCR :== 3)
	vBCR4 = (currentBCR :== 4)
	vBCR5 = (currentBCR :== 5)
	vBCR6 = (currentBCR :== 6)
	
	// Determine coefficient segments based on Line
	// Segment maps Line×BCR coefficient position: BCR_start = 10 + segment*6
	if (Line == 0) { 
		// DN or L1S: No BCR yet
		segments = 0
	}
	else if (Line == 1) { 
		// L1E: Split by ASCT status
		segments = (1, 2)
	}
	else {
		// L2+ End phases
		segments = Line + 1
		
		// Cap at segment 7 (L6) - reuse L6 coefficients for L7-L9
		if (segments > 7) {
			segments = 7
		}	
	}
 
    // Loop through segments
    for (s = 1; s <= cols(segments); s++) {
        segment = segments[s]
        
        // Calculate BCR coefficient start position
        bcrStart = 10 + segment * 6
        
        // Build coefficient column vector
        coefCols = (1, 2, 3, 5, 6, 8, 9,                           // Base effects
                    bcrStart, bcrStart+1, bcrStart+2,              // BCR 1-3
                    bcrStart+3, bcrStart+4, bcrStart+5,            // BCR 4-6
                    58)                                            // Constant
        
		// Determine patients for this segment
		if (segment == 0) { // DN
			if (OMC == 1) {
				idx = selectindex(mState[., 1] :<= OMC)
			}
			else {
				idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
			}
		}
		else if (segment == 1) { // L1 No ASCT
			idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vSCT_DN :== 0))
		}
		else if (segment == 2) { // L1 ASCT
			idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vSCT_DN :== 1))
		}
		else {
			// L2+ segments: all alive patients
			idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
		}
        
        // Calculate for this segment if patients exist
        if (rows(idx) > 0) {
            // Build patient matrix - reference ECOG/RISS categories NOT required as using coefCols
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
