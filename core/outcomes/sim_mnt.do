**********
* Monash Myeloma Model - Sim MNT
*
* Purpose: Determine receipt of Maintenance Therapy via logistic regression. Binary outcome
*          (0 = no MNT, 1 = MNT).
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
					
	// Extract L1 treatment regimen
	L1_TXR = mTXR[., 1]
	
	// Determine number of treatment regimen dummies
	nRegimens = cols(oL1_TXR)
	
	// Filter for alive and eligible patients
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
	if (rows(idx) > 0) {

		// Create treatment regimen dummies (incl. the base regimen: coefficient vectors carry 0b.TXR_L1)
		TXR_is_R1 = (L1_TXR :== oL1_TXR[1, 1])
		if (nRegimens >= 2) TXR_is_R2 = (L1_TXR :== oL1_TXR[1, 2])
		if (nRegimens >= 3) TXR_is_R3 = (L1_TXR :== oL1_TXR[1, 3])
		if (nRegimens >= 4) TXR_is_R4 = (L1_TXR :== oL1_TXR[1, 4])
		
		// ASCT Group
		idxASCT = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vSCT_L1 :== 1))
		if (rows(idxASCT) > 0) {
			
			// Extract ASCT BCR
			vBCR_ASCT = mBCR[idxASCT, 10]
			
			// Create BCR dummy variables 
			vBCR_ASCT_1 = (vBCR_ASCT :== 1)
			vBCR_ASCT_2 = (vBCR_ASCT :== 2)
			vBCR_ASCT_3 = (vBCR_ASCT :== 3)
			vBCR_ASCT_4 = (vBCR_ASCT :== 4)		
			
			// Assemble patient matrix
			mPat = (vAge[idxASCT], vAge2[idxASCT], vMale[idxASCT], 
					vECOG0[idxASCT], vECOG1[idxASCT], vECOG2[idxASCT], 
					vRISS1[idxASCT], vRISS2[idxASCT], vRISS3[idxASCT])
			
			// Add treatment regimen dummies (base regimen first, to match the coefficient layout)
			mPat = mPat, TXR_is_R1[idxASCT]
			if (nRegimens >= 2) mPat = mPat, TXR_is_R2[idxASCT]
			if (nRegimens >= 3) mPat = mPat, TXR_is_R3[idxASCT]
			if (nRegimens >= 4) mPat = mPat, TXR_is_R4[idxASCT]
			
			// Add BCR dummies
			mPat = mPat, (vBCR_ASCT_1, vBCR_ASCT_2, 
						  vBCR_ASCT_3, vBCR_ASCT_4)
			
			// Add constant
			mPat = mPat, vCons[idxASCT]
			
			// Extract coefficients for ASCT group
			nPredictors = cols(mPat)
			// Guard: design columns must equal the coefficient count (no cutpoints/ancillary here).
			if (nPredictors != cols(bMNT_ASCT)) {
				errprintf("sim_mnt (ASCT): design/coefficient mismatch - mPat has %g columns but coefficient vector has %g\n", nPredictors, cols(bMNT_ASCT))
				exit(459)
			}
			vCoef = bMNT_ASCT[1, 1..nPredictors]'
			
			// Calculate XB
			vXB = mPat * vCoef
			
			// Calculate probability
			vPR = 1 :/ (1 :+ exp(-vXB))
					
			// Determine outcome
			vRN = rnDraw(idxASCT, rn_mnt(1))
			vOC = (vPR :> vRN) :* 1
			
			// Update vector
			vMNT[idxASCT] = vOC
		}

		// No ASCT Group
		idxNoASCT = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vSCT_L1 :== 0))
		if (rows(idxNoASCT) > 0) {
			
			// Extract L1 BCR
			vBCR_L1 = mBCR[idxNoASCT, 1]
			
			// Create BCR dummy variables 
			vBCR_L1_1 = (vBCR_L1 :== 1)
			vBCR_L1_2 = (vBCR_L1 :== 2)
			vBCR_L1_3 = (vBCR_L1 :== 3)
			vBCR_L1_4 = (vBCR_L1 :== 4)
			vBCR_L1_5 = (vBCR_L1 :== 5)
			vBCR_L1_6 = (vBCR_L1 :== 6)
			
			// Assemble patient matrix
			mPat = (vAge[idxNoASCT], vAge2[idxNoASCT], vMale[idxNoASCT], 
					vECOG0[idxNoASCT], vECOG1[idxNoASCT], vECOG2[idxNoASCT], 
					vRISS1[idxNoASCT], vRISS2[idxNoASCT], vRISS3[idxNoASCT])
			
			// Add treatment regimen dummies (base regimen first, to match the coefficient layout)
			mPat = mPat, TXR_is_R1[idxNoASCT]
			if (nRegimens >= 2) mPat = mPat, TXR_is_R2[idxNoASCT]
			if (nRegimens >= 3) mPat = mPat, TXR_is_R3[idxNoASCT]
			if (nRegimens >= 4) mPat = mPat, TXR_is_R4[idxNoASCT]
			
			// Add BCR dummies
			mPat = mPat, (vBCR_L1_1, vBCR_L1_2, vBCR_L1_3, 
						  vBCR_L1_4, vBCR_L1_5, vBCR_L1_6)
			
			// Add constant
			mPat = mPat, vCons[idxNoASCT]
			
			// Extract coefficients for No ASCT group
			nPredictors = cols(mPat)
			// Guard: design columns must equal the coefficient count (no cutpoints/ancillary here).
			if (nPredictors != cols(bMNT_NoASCT)) {
				errprintf("sim_mnt (NoASCT): design/coefficient mismatch - mPat has %g columns but coefficient vector has %g\n", nPredictors, cols(bMNT_NoASCT))
				exit(459)
			}
			vCoef = bMNT_NoASCT[1, 1..nPredictors]'
			
			// Calculate XB
			vXB = mPat * vCoef
			
			// Calculate probability
			vPR = 1 :/ (1 :+ exp(-vXB))
					
			// Determine outcome
			vRN = rnDraw(idxNoASCT, rn_mnt(2))
			vOC = (vPR :> vRN) :* 1
			
			// Update vector
			vMNT[idxNoASCT] = vOC
		}
	}
}
