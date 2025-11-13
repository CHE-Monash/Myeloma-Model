**********
* SIM BCR OVERRIDE - DVd L2 Method
*
* Purpose: Override BCR for DVd at Line 2 using Common Comparator Method
*
* Context: Called AFTER standard sim_bcr.do has run
*          Applies to entire cohort (all patients receive DVd at L2)
*          Only active when $Intervention = "DVd" AND Line = 2
*
* Author: Adam Irving
* Date: November 2025
**********

// Load Common Comparator Predictions		
local bsIteration = "$BSIteration"
local bcrDir = "$analysis_path/data/bcr"
		
if ("`bsIteration'" != "" & "`bsIteration'" != "0") {
    // Bootstrap mode
    capture mata: mata matuse "`bcrDir'/bootstrap/bcr_dvd_l2_B`bsIteration'.mmat"
}
else {
    // Point estimate mode
    capture mata: mata matuse "`bcrDir'/bcr_dvd_l2.mmat"
}

mata {
    // Get alive, non-prevalent patients
    idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
    if (rows(idx) > 0) {
		
	nPatients = rows(idx)
		
	// Calculate XB
		
		// Extract previous BCR for L2
		BCR_L1 = mBCR[., 1]
		pBCR_CR = (BCR_L1 :== 1)
		pBCR_VG = (BCR_L1 :== 2)
		pBCR_PR = (BCR_L1 :== 3)
		pBCR_MR = (BCR_L1 :== 4)
		pBCR_SD = (BCR_L1 :== 5)
		pBCR_PD = (BCR_L1 :== 6)
		
		BCR_SCT = mBCR[., 10]
		pBCR_SCT_0 = (BCR_SCT :== 0)
		pBCR_SCT_1 = (BCR_SCT :== 1)
		pBCR_SCT_2 = (BCR_SCT :== 2)
		pBCR_SCT_3 = (BCR_SCT :== 3)
		pBCR_SCT_4 = (BCR_SCT :== 4)
		
		// Build patient matrix 
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx],
				vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		mPat = mPat, (pBCR_CR[idx], pBCR_VG[idx], pBCR_PR[idx],
					  pBCR_MR[idx], pBCR_SD[idx], pBCR_PD[idx])
						
		mPat = mPat, (pBCR_SCT_0[idx], pBCR_SCT_1[idx], pBCR_SCT_2[idx], 
					  pBCR_SCT_3[idx], pBCR_SCT_4[idx])
		
		// Get coefficients (excluding TXR coefficients)
		vCoef = bL2_BCR
		nPredictors = cols(mPat)
		vCoef = vCoef[1, 1..nPredictors]'
		
		// Calculate prognostic index (XB)
		XB = mPat * vCoef

	// Assign BCR based on prognosis ranking
	
		// Sort patients by XB (higher XB = worse prognosis)
		sortorder = order(XB, -1)
		
		// Calculate cumulative probabilities
		cumProbs = runningsum(mBCR_DVd_L2')
		
		// Assign outcomes deterministically based on prognosis ranking
		vOut = J(nPatients, 1, .)
		
		for (i = 1; i <= nPatients; i++) {
			// Find this patient's rank (1 = worst prognosis, nPatients = best)
			rank = sortorder[i]
			
			// What percentile is this rank?
			percentile = rank / nPatients
			
			// Assign BCR category based on cumulative distribution
			if (percentile <= cumProbs[1]) {
				vOut[i] = 1  // CR (best)
			}
			else if (percentile <= cumProbs[2]) {
				vOut[i] = 2  // VGPR
			}
			else if (percentile <= cumProbs[3]) {
				vOut[i] = 3  // PR
			}
			else if (percentile <= cumProbs[4]) {
				vOut[i] = 4  // MR 
			}
			else if (percentile <= cumProbs[5]) {
				vOut[i] = 5  // SD
			}
			else {
				vOut[i] = 6  // PD (worst) 
			}
		}
		
	// Replace BCR outcome with Common Comparator predictions
	mBCR[idx, LX] = vOut
	
	}
}
