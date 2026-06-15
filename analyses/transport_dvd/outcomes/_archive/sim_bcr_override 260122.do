**********
* SIM BCR OVERRIDE - DVd L2 Method
*
* Purpose: Override BCR for DVd/Vd at Line 2 
* Method: Called after sim_bcr.do
*
* Author: Adam Irving
* Date: December 2025
**********

capture program drop bcr_override_xb
program define bcr_override_xb
    * Common setup: patient selection, covariate assembly, XB calculation
    * Returns: idx, vXB, nPatients stored via st_matrix/st_numscalar
	args coef_matrix
	mata {
		// Get alive, non-prevalent patients
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))		
		if (rows(idx) > 0) {
			
			nPatients = rows(idx)

			// Extract previous BCR conditional on SCT status
			pBCR = (vSCT_L1 :== 0) :* mBCR[., 1] + (vSCT_L1 :== 1) :* mBCR[., 10]

			// Create dummy indicators
			pBCR_1 = (pBCR :== 1)
			pBCR_2 = (pBCR :== 2)
			pBCR_3 = (pBCR :== 3)
			pBCR_4 = (pBCR :== 4)
			pBCR_5 = (pBCR :== 5)
			pBCR_6 = (pBCR :== 6)
			
			// Assemble patient matrix (without TXR)
			mPat = (vAge[idx], vAge2[idx], vMale[idx])
			
			mPat = mPat, (pBCR_1[idx], pBCR_2[idx], pBCR_3[idx],
						  pBCR_4[idx], pBCR_5[idx], pBCR_6[idx])
			
			// Add treatment indicators for scenario B
			if (st_global("scenario") == "B_transport") {
				MRDR = J(nPatients, 1, 1)  // Always 1					
				if (st_global("int") == "dvd") {
					DVd = J(nPatients, 1, 1)
					Vd = J(nPatients, 1, 0)
				}
				else {
					DVd = J(nPatients, 1, 0)
					Vd = J(nPatients, 1, 1)
				}	
				mPat = mPat, (MRDR, DVd, Vd)
			}			  
		
			// Get coefficients (excluding TXR)
			vCoef_full = *findexternal(st_local("coef_matrix"))
						
			nPredictors = cols(mPat)
			vCoef = vCoef_full[1, 1..nPredictors]'
			
			// Calculate XB
			vXB = mPat * vCoef
		
			// Store for allocation functions via Stata
			st_matrix("bcr_idx", idx)
			st_matrix("bcr_vXB", vXB)
			st_numscalar("bcr_nPatients", nPatients)
		}
		else {
			// No patients - store zero
			st_numscalar("bcr_nPatients", 0)
		}
	}
end		

capture program drop A_trial
program define A_trial
    * Rank-based allocation using target distribution
    args bcr_matrix
	mata {	
		// Retrieve from Stata
		bcr_idx = st_matrix("bcr_idx")
		bcr_vXB = st_matrix("bcr_vXB")
		bcr_nPatients = st_numscalar("bcr_nPatients")
		
		stata(`"noi di as text "  bcr_nPatients: "' + strofreal(bcr_nPatients))
		
		if (bcr_nPatients > 0) {		
			// Get BCR probability matrix by name
			mBCR_probs = *findexternal(st_local("bcr_matrix"))
						
			// Sort patients by XB (higher XB = worse prognosis)
			sortorder = order(bcr_vXB, -1)
			
			// Calculate cumulative probabilities 
			cumProbs = runningsum(mBCR_probs)
						
			// Assign outcomes deterministically based on prognosis ranking
			vOC = J(bcr_nPatients, 1, .)
			
			for (i = 1; i <= bcr_nPatients; i++) {
				// Find this patient's rank (1 = worst prognosis, nPatients = best)
				rank = sortorder[i]
				
				// What percentile is this rank?
				percentile = rank / bcr_nPatients
				
				// Assign BCR category based on cumulative distribution
				if (percentile <= cumProbs[1]) {
					vOC[i] = 1  // CR (best)
				}
				else if (percentile <= cumProbs[2]) {
					vOC[i] = 2  // VGPR
				}
				else if (percentile <= cumProbs[3]) {
					vOC[i] = 3  // PR
				}
				else if (percentile <= cumProbs[4]) {
					vOC[i] = 4  // MR 
				}
				else if (percentile <= cumProbs[5]) {
					vOC[i] = 5  // SD
				}
				else {
					vOC[i] = 6  // PD (worst) 
				}
			}
			
			// Print BCR distribution
			for (j = 1; j <= 6; j++) {
				count = sum(vOC :== j)
				pct = 100 * count / bcr_nPatients
			}
			
			// Override BCR - convert bcr_idx to column vector for indexing
			bcr_idx = bcr_idx'
			mBCR[bcr_idx, Line] = vOC
		}
	}
end

capture program drop B_transport
program define B_transport
    * Cutpoint-based allocation from transported regression
    args coef_matrix
    mata {	
		// Retrieve from Stata
		bcr_idx = st_matrix("bcr_idx")
		bcr_vXB = st_matrix("bcr_vXB")
		bcr_nPatients = st_numscalar("bcr_nPatients")
		        
        if (bcr_nPatients > 0) {
            vCoef_full = *findexternal(st_local("coef_matrix"))
            
            // Extract cutpoints from final 5 cells
            nCutPoints = 5
            cutPointIndices = (cols(vCoef_full) - nCutPoints + 1)..cols(vCoef_full)
            cutPoints = vCoef_full[1, cutPointIndices]
			           
            // Calculate probabilities from ordered logit
            cumProbs = calcOrdLogitProbs(bcr_vXB, cutPoints)
            
            // Stochastic assignment
            vRN = runiform(bcr_nPatients, 1)
			categoryValues = (1, 2, 3, 4, 5, 6)
            vOC = assignOrdOutcome(vRN, cumProbs, categoryValues)
            
			// Override BCR - convert bcr_idx to column vector for indexing
			bcr_idx = bcr_idx'
            mBCR[bcr_idx, Line] = vOC
        }
    }
end

// Call functions based on $boot, $int & $scenario
if ($boot == 0) {
	if ("$scenario" == "A_trial") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/A_trial/bcr_dvd_l2.mmat", replace
			bcr_override_xb bL2_BCR
			A_trial mBCR_DVd_L2
		}
		if ("$int" == "vd") {
			mata: mata matuse "$outcomes_path/A_trial/bcr_vd_l2.mmat", replace
			bcr_override_xb bL2_BCR
			A_trial mBCR_Vd_L2
		}
	}
	else if ("$scenario" == "B_transport") {
		mata: mata matuse "$outcomes_path/B_transport/transport_dvd.mmat", replace
		bcr_override_xb bL2_BCR_T
		B_transport bL2_BCR_T
	}
}
else if ($boot == 1) {
	if ("$scenario" == "A_trial") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/A_trial/bootstrap/bcr_dvd_l2_B{$b}.mmat", replace
			bcr_override_xb bL2_BCR
			A_trial mBCR_DVd_L2
		}
		if ("$int" == "vd") {
			mata: mata matuse "$outcomes_path/A_trial/bootstrap/bcr_vd_l2_B{$b}.mmat", replace
			bcr_override_xb bL2_BCR
			A_trial mBCR_Vd_L2
		}
	}
	else if ("$scenario" == "B_transport") {
		mata: mata matuse "$outcomes_path/B_transport/bootstrap/transport_dvd_B{$b}.mmat", replace
		bcr_override_xb bL2_BCR_T
		B_transport bL2_BCR_T
	}
}
