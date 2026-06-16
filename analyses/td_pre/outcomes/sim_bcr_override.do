**********
* SIM BCR OVERRIDE - TecDara Method
*
* Purpose: Override BCR at Line 2, 3 or 4 
* Method: Called after sim_bcr.do if Line = $line
*
* Author: Adam Irving
* Date: December 2025
**********

capture program drop bcr_override_xb
program define bcr_override_xb
	args coef_matrix
	mata {
		// Get alive, non-prevalent patients
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))		
		if (rows(idx) > 0) {
			
			nPatients = rows(idx)
			
			// Extract previous BCR
			BCR = mBCR[., Line-1]
			BCR_1 = (BCR :== 1)
			BCR_2 = (BCR :== 2)
			BCR_3 = (BCR :== 3)
			BCR_4 = (BCR :== 4)
			BCR_5 = (BCR :== 5)
			BCR_6 = (BCR :== 6)
			
			// Assemble patient matrix
			mPat = (vAge[idx], vAge2[idx], vMale[idx], 
					vECOG0[idx], vECOG1[idx], vECOG2[idx],
					BCR_1[idx], BCR_2[idx], BCR_3[idx],
					BCR_4[idx], BCR_5[idx], BCR_6[idx])
			
			// Add treatment indicators
			MRDR = J(nPatients, 1, 1)  // Always 1
			DVd = J(nPatients, 1, 0) // Always 0
                
			if (st_global("int") == "td") {
				TD = J(nPatients, 1, 1)
			}
			else {  // soc
				TD = J(nPatients, 1, 0)
			}
                
			mPat = mPat, (MRDR, TD, DVd)
			
			// Add line indicators
			if (Line == 2) {
				L2 = J(nPatients, 1, 1)
				L3 = J(nPatients, 1, 0)
				L4 = J(nPatients, 1, 0)
			}
			else if (Line == 3) {
				L2 = J(nPatients, 1, 0)
				L3 = J(nPatients, 1, 1)
				L4 = J(nPatients, 1, 0)
			}
			else if (Line == 4) {
				L2 = J(nPatients, 1, 0)
				L3 = J(nPatients, 1, 0)
				L4 = J(nPatients, 1, 1)
			}
			
			mPat = mPat, (L2, L3, L4)
		
			// Get coefficients (excluding TXR)
			vCoef_full = *findexternal(st_local("coef_matrix"))
						
			nPredictors = cols(mPat)
			vCoef = vCoef_full[1, 1..nPredictors]'
			
			// Calculate XB
			vXB = mPat * vCoef
            
            // Extract cutpoints from final 5 cells
            nCutPoints = 5
            cutPointIndices = (cols(vCoef_full) - nCutPoints + 1)..cols(vCoef_full)
            cutPoints = vCoef_full[1, cutPointIndices]
			            
            // Calculate probabilities from ordered logit
            cumProbs = calcOrdLogitProbs(vXB, cutPoints)

            // Stochastic assignment (CRN: same column as core BCR -> aligned across arms)
            vRN = rnDraw(idx, rn_bcr(Line))
			categoryValues = (1, 2, 3, 4, 5, 6)
            vOC = assignOrdOutcome(vRN, cumProbs, categoryValues)
            
			// Override BCR - convert bcr_idx to column vector for indexing
            mBCR[idx, Line] = vOC
        }
    }
end

// Call functions based on $boot 

if ($boot == 0) {
	if ("$int" == "td") {
		mata: mata matuse "$outcomes_path/transport/transport_td.mmat", replace
		bcr_override_xb bL4_BCR_T
	}
}
else if ($boot == 1) {
	if ("$int" == "td") {
		mata: mata matuse "$outcomes_path/transport/bootstrap/transport_td_B${b}.mmat"
		bcr_override_xb bL4_BCR_T
	}
}
