**********
* SIM TXR
*
* Purpose: Select treatment regimen
* Method: Multinomial logit
* Outcome: Regimen number
**********

mata {
	if (txr_model_exists(Line)) {
	// Model exists - run multinomial logit
			
	vCoef = get_txr_coef(Line)
	vTXR = get_txr_outcome(Line)
		
		// Filter for alive and eligible
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
		if (rows(idx) > 0) {
		
			nRegimens = cols(vTXR)
			
			vXB1 = J(rows(idx), 1, 1)
			vXB2 = J(rows(idx), 1, 0)
			vXB3 = J(rows(idx), 1, 0)
			vXB4 = J(rows(idx), 1, 0)
			
			// Assemble patient matrix (L1)
			if (Line == 0) {
				mPat = (vAge[idx], vAge2[idx], vMale[idx], 
						vECOG0[idx], vECOG1[idx], vECOG2[idx], 
						vRISS1[idx], vRISS2[idx], vRISS3[idx], 
						vSCT_DN[idx], vCons[idx])
				
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
				vBCR_1 = (prevBCR :== 1)
				vBCR_2 = (prevBCR :== 2)
				vBCR_3 = (prevBCR :== 3)
				vBCR_4 = (prevBCR :== 4)
				vBCR_5 = (prevBCR :== 5)
				vBCR_6 = (prevBCR :== 6)
				
				mPat = (vAge[idx], vAge2[idx], vMale[idx], 
						vECOG0[idx], vECOG1[idx], vECOG2[idx], 
						vRISS1[idx], vRISS2[idx], vRISS3[idx],  
						vBCR_1, vBCR_2, vBCR_3, vBCR_4, vBCR_5, vBCR_6, 
						vCons[idx])
				
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
				vBCR_1 = (prevBCR :== 1)
				vBCR_3 = (prevBCR :== 3)
				vBCR_5 = (prevBCR :== 5)
				
				mPat = (vAge[idx], vAge2[idx], vMale[idx], 
						vECOG0[idx], vECOG1[idx], vECOG2[idx], 
						vRISS1[idx], vRISS2[idx], vRISS3[idx],  
						vBCR_1, vBCR_3, vBCR_5, 
						vCons[idx])
						
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
			
			// Calculate probabilities
			vPR1 = vXB1 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4)
			vPR2 = vPR1 :+ (vXB2 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR3 = vPR2 :+ (vXB3 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR4 = vPR3 :+ (vXB4 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			
			// Determine Outcome
			vRN = runiform(rows(idx), 1)
			vOC = (vRN :< vPR1) :* vTXR[1,1]

			if (nRegimens >= 2) {
				vOC = vOC :+ ((vRN :>= vPR1) :& (vRN :< vPR2)) :* vTXR[1,2]
			}
			
			if (nRegimens >= 3) {
				vOC = vOC :+ ((vRN :>= vPR2) :& (vRN :< vPR3)) :* vTXR[1,3]
			}
			
			if (nRegimens >= 4) {
				vOC = vOC :+ ((vRN :>= vPR3) :& (vRN :< vPR4)) :* vTXR[1,4]
			}
			
			// Update matrix
			mTXR[idx, Line+1] = vOC
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
if "${line}" == "`current_line'" {
	local override_file "${analysis_path}/outcomes/sim_txr_override_${int}_l${line}.do"
	capture confirm file "`override_file'"
	if _rc == 0 {
		di "Overriding TXR"
		quietly do "`override_file'"
	}
}

