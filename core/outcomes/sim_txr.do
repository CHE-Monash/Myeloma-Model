**********
*SIM TXR

*Purpose: Select treatment regimen
*Method: Multinomial logit
*Outcome: Regimen number
**********

mata {
	
	if (txr_model_exists(Line)) {
	// Model exists - run multinomial logit
			
	vCoef = get_txr_coef(Line)
	vTXR = get_txr_outcome(Line)
		
		// Select patients
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
		if (rows(idx) > 0) {
		
			nRegimens = cols(vTXR)
			
			vRN = runiform(rows(mMOR), 1)
			
			vXB1 = J(rows(idx), 1, 1)
			vXB2 = J(rows(idx), 1, 0)
			vXB3 = J(rows(idx), 1, 0)
			vXB4 = J(rows(idx), 1, 0)
			
			// Patient characteristics
			vAge_e = mAge[idx, OMC]
			vAge2_e = vAge_e :^ 2
			vMale_e = vMale[idx]
			vECOG0_e = (vECOG[idx] :== 0)
			vECOG1_e = (vECOG[idx] :== 1)
			vECOG2_e = (vECOG[idx] :== 2)
			vRISS1_e = (vRISS[idx] :== 1)
			vRISS2_e = (vRISS[idx] :== 2)
			vRISS3_e = (vRISS[idx] :== 3)
			vCons_e = vCons[idx]
			
			// Build patient matrix (L1)
			if (Line == 0) {
				vSCT_e = vSCT_DN[idx]
				mPat = (vAge_e, vAge2_e, vMale_e, 
						   vECOG0_e, vECOG1_e, vECOG2_e, 
						   vRISS1_e, vRISS2_e, vRISS3_e, 
						   vSCT_e, 
						   vCons_e)
				
				nPredictors = cols(mPat)

				if (nRegimens >= 2) {
					coef_XB2 = vCoef[1, (nPredictors+1)..(2*nPredictors)]'
					vXB2 = exp(mPat * coef_XB2)
				}

				if (nRegimens >= 3) {
					coef_XB3 = vCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
					vXB3 = exp(mPat * coef_XB3)
				}

				if (nRegimens >= 4) {
					coef_XB4 = vCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
					vXB4 = exp(mPat * coef_XB4)
				}
			}
			
			// Build patient matrix (L2 or L3)
			else if (Line == 1 | Line == 2) {
				prevBCR = mBCR[idx, Line]
				pBCR_1 = (prevBCR :== 1)
				pBCR_2 = (prevBCR :== 2)
				pBCR_3 = (prevBCR :== 3)
				pBCR_4 = (prevBCR :== 4)
				pBCR_5 = (prevBCR :== 5)
				pBCR_6 = (prevBCR :== 6)
				
				mPat = (vAge_e, vAge2_e, vMale_e, 
						vECOG0_e, vECOG1_e, vECOG2_e, 
						vRISS1_e, vRISS2_e, vRISS3_e, 
						pBCR_1, pBCR_2, pBCR_3, pBCR_4, pBCR_5, pBCR_6, 
						vCons_e)
				
				nPredictors = cols(mPat)

				if (nRegimens >= 2) {
					coef_XB2 = vCoef[1, (nPredictors+1)..(2*nPredictors)]'
					vXB2 = exp(mPat * coef_XB2)
				}

				if (nRegimens >= 3) {
					coef_XB3 = vCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
					vXB3 = exp(mPat * coef_XB3)
				}

				if (nRegimens >= 4) {
					coef_XB4 = vCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
					vXB4 = exp(mPat * coef_XB4)
				}
			}
			// Build patient matrix (L4+)
			else {
				prevBCR = mBCR[idx, Line]
				pBCR_1 = (prevBCR :== 1)
				pBCR_3 = (prevBCR :== 3)
				pBCR_5 = (prevBCR :== 5)
				
				mPat = (vAge_e, vAge2_e, vMale_e, 
						vECOG0_e, vECOG1_e, vECOG2_e, 
						vRISS1_e, vRISS2_e, vRISS3_e, 
						pBCR_1, pBCR_3, pBCR_5, 
						vCons_e)
				
				if (nRegimens >= 2) {
					coef_XB2 = vCoef[1, (nPredictors+1)..(2*nPredictors)]'
					vXB2 = exp(mPat * coef_XB2)
				}
				
				if (nRegimens >= 3) {
					coef_XB3 = vCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
					vXB3 = exp(mPat * coef_XB3)
				}
				
				if (nRegimens >= 4) {
					coef_XB4 = vCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
					vXB4 = exp(mPat * coef_XB4)
				}
			}
			
			// Calculate probabilities
			vPR1 = vXB1 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4)
			vPR2 = vPR1 :+ (vXB2 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR3 = vPR2 :+ (vXB3 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR4 = vPR3 :+ (vXB4 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			
			vRN_e = vRN[idx]
			vOut = J(rows(idx), 1, .)
			
			// Determine Outcome
			vOut = (vRN_e :< vPR1) :* vTXR[1,1]

			if (nRegimens >= 2) {
				vOut = vOut :+ ((vRN_e :>= vPR1) :& (vRN_e :< vPR2)) :* vTXR[1,2]
			}
			
			if (nRegimens >= 3) {
				vOut = vOut :+ ((vRN_e :>= vPR2) :& (vRN_e :< vPR3)) :* vTXR[1,3]
			}
			
			if (nRegimens >= 4) {
				vOut = vOut :+ ((vRN_e :>= vPR3) :& (vRN_e :< vPR4)) :* vTXR[1,4]
			}
			
			// Update matrix
			mTXR[idx, Line+1] = vOut
		}	
	}
	else {
		// No model - assign pooled "Other" treatment
			
		idxAlive = selectindex(mMOR[., OMC-1] :== 0)
		idxDead = selectindex(mMOR[., OMC-1] :!= 0)
			
		// Assign "Other" (code 0) to alive patients
		if (rows(idxAlive) > 0) {
			mTXR[idxAlive, LX+1] = J(rows(idxAlive), 1, 0)
		}
			
		// Set missing for dead patients
		if (rows(idxDead) > 0) {
			mTXR[idxDead, LX+1] = J(rows(idxDead), 1, .)
		}
			
	}
}

// Check for override file, execute if it exists
mata: st_local("current_line", strofreal(Line+1))
if "${Line}" == "`current_line'" {
	local override_file "${analysis_path}/outcomes/sim_txr_override_${Int}_l${Line}.do"
	capture confirm file "`override_file'"
	if _rc == 0 {
		di "Overriding TXR"
		quietly do "`override_file'"
	}
}

