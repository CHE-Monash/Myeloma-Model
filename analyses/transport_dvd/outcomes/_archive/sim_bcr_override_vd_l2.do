**********
* SIM BCR OVERRIDE - DVd L2 Method
*
* Purpose: Override BCR for Vd at Line 2 
*
* Context: Called AFTER standard sim_bcr.do has run
*          Applies to entire cohort (all patients receive DVd at L2)
*          Only active $int = "Vd" AND $line = 2
*
* Author: Adam Irving
* Date: November 2025
**********

// Load BCR distribution based on bootstrap & scenario
if (${boot} == 0) {
	mata: mata matuse "${analysis_path}/data/bcr/A_trial/bcr_vd_l2.mmat"
}
else if (${boot} == 1) {
	mata: mata matuse "${analysis_path}/data/bcr/$A_trial/bootstrap/bcr_vd_l2_B`b'.mmat"
}

mata {
    // Get alive, non-prevalent patients
    idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
    if (rows(idx) > 0) {
		
	nPatients = rows(idx)
		
	// Extract previous BCR
	BCR_L1 = mBCR[., 1]
	pBCR_L1_CR = (BCR_L1 :== 1)
	pBCR_L1_VG = (BCR_L1 :== 2)
	pBCR_L1_PR = (BCR_L1 :== 3)
	pBCR_L1_MR = (BCR_L1 :== 4)
	pBCR_L1_SD = (BCR_L1 :== 5)
	pBCR_L1_PD = (BCR_L1 :== 6)
		
	BCR_SCT = mBCR[., 10]
	pBCR_SCT_0 = (BCR_SCT :== 0)
	pBCR_SCT_CR = (BCR_SCT :== 1)
	pBCR_SCT_VG = (BCR_SCT :== 2)
	pBCR_SCT_PR = (BCR_SCT :== 3)
	pBCR_SCT_MR = (BCR_SCT :== 4)
		
	// Assemble patient matrix (without TXR)
	mPat = (vAge[idx], vAge2[idx], vMale[idx], 
			vECOG0[idx], vECOG1[idx], vECOG2[idx],
			vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
	mPat = mPat, (pBCR_L1_CR[idx], pBCR_L1_VG[idx], pBCR_L1_PR[idx],
				  pBCR_L1_MR[idx], pBCR_L1_SD[idx], pBCR_L1_PD[idx])
						
	mPat = mPat, (pBCR_SCT_0[idx], pBCR_SCT_CR[idx], pBCR_SCT_VG[idx], 
				  pBCR_SCT_PR[idx], pBCR_SCT_MR[idx])
		
	// Get coefficients (excluding TXR)
	vCoef = bL2_BCR
	nPredictors = cols(mPat)
	vCoef = vCoef[1, 1..nPredictors]'
		
	// Calculate prognostic index (vXB)
	vXB = mPat * vCoef

	// Assign BCR based on prognosis ranking
	
		// Sort patients by vXB (higher vXB = worse prognosis)
		sortorder = order(vXB, -1)
		
		// Calculate cumulative probabilities
		cumProbs = runningsum(mBCR_Vd_L2')
		
		// Assign outcomes deterministically based on prognosis ranking
		vOC = J(nPatients, 1, .)
		
		for (i = 1; i <= nPatients; i++) {
			// Find this patient's rank (1 = worst prognosis, nPatients = best)
			rank = sortorder[i]
			
			// What percentile is this rank?
			percentile = rank / nPatients
			
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
		
	// Override BCR
	mBCR[idx, Line] = vOC
	
	}
}
