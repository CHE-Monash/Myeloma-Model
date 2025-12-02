**********
* SIM TXD L1
*
* Purpose: Treatment Duration for Line 1
* Method: Parametric survival analysis
* Outcome: Continuous time (months)
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC)
	idx = selectindex(vEligible)
	if (rows(idx) > 0) {
		
		// Group 1: Fixed + ASCT (with splines)
		vASCT = (mTXR[.,1] :!= 7) :& (vSCT_DN :== 1) :& vEligible
		idxASCT = selectindex(vASCT)
		if (rows(idxASCT) > 0) {
			
			// Assemble patient matrix
			mPat_ASCT = (vAge[idxASCT], vAge2[idxASCT], vMale[idxASCT], 
					   vECOG0[idxASCT], vECOG1[idxASCT], vECOG2[idxASCT], 
					   vRISS1[idxASCT], vRISS2[idxASCT], vRISS3[idxASCT])
			nCR = cols(oL1_TXR)
			if (nCR >= 1) {
				vCR1_ASCT = (mTXR[idxASCT, 1] :== oL1_TXR[1,1])
				mPat_ASCT = (mPat_ASCT , vCR1_ASCT)
			}
			if (nCR >= 2) {
				vCR2_ASCT = (mTXR[idxASCT, 1] :== oL1_TXR[1,2])
				mPat_ASCT = (mPat_ASCT , vCR2_ASCT)
			}
			if (nCR >= 3) {
				vCR3_ASCT = (mTXR[idxASCT, 1] :== oL1_TXR[1,3])
				mPat_ASCT = (mPat_ASCT , vCR3_ASCT)	
			}
			mPat_ASCT = (mPat_ASCT, vCons[idxASCT])

			// --- SPLINE 1 ---
				// Extract coefficients
				nPredictors = cols(mPat_ASCT)
				coef_S1 = bL1_TXD_ASCT_S1[1, 1..nPredictors]'
				aux_S1 = bL1_TXD_ASCT_S1[1, cols(bL1_TXD_ASCT_S1)]
				
				// Calculate XB
				vXB_S1 = mPat_ASCT * coef_S1
				
				// Calculate survival time
				vRN_ASCT = runiform(rows(idxASCT), 1)
				vDur_S1 = calcSurvTime(vXB_S1, vRN_ASCT, fbL1_TXD_ASCT_S1, aux_S1)
			
				// --- SPLINE 2 (for those beyond Cutoff 1) ---
					vBeyondC1 = (vDur_S1 :> L1_TXD_ASCT_C1)
					idxBeyondC1 = selectindex(vBeyondC1)

					if (rows(idxBeyondC1) > 0) {
						// Patient matrix
						mPat_ASCT_C1 = mPat_ASCT[idxBeyondC1, .]
						
						// Extract coefficients
						nPredictors = cols(mPat_ASCT_C1)
						coef_S2 = bL1_TXD_ASCT_S2[1, 1..nPredictors]'
						aux_S2 = bL1_TXD_ASCT_S2[1, cols(bL1_TXD_ASCT_S2)]
						
						// Calculate XB
						vXB_S2 = mPat_ASCT_C1 * coef_S2
						
						// Conditional RN: survive past C1
						vSurvC1 = exp(-exp(vXB_S2) :* (L1_TXD_ASCT_C1 :^ exp(aux_S2)))
						vRN_S2 = runiform(rows(idxBeyondC1), 1) :* vSurvC1
						
						// Recalculate duration with Spline 2
						vDur_S2 = calcSurvTime(vXB_S2, vRN_S2, fbL1_TXD_ASCT_S2, aux_S2)
						
						// Update durations for those beyond C1
						vDur_S1[idxBeyondC1] = vDur_S2
				
					// --- SPLINE 3 (for those beyond Cutoff 2) ---
						vBeyondC2 = (vDur_S2 :> L1_TXD_ASCT_C2)
						idxBeyondC2_local = selectindex(vBeyondC2)
						
						if (rows(idxBeyondC2_local) > 0) {
							// Patient matrix
							idxBeyondC2 = idxBeyondC1[idxBeyondC2_local]
							mPat_ASCT_C2 = mPat_ASCT[idxBeyondC2, .]
							
							// Extract coefficients
							nPredictors = cols(mPat_ASCT_C2)
							coef_S3 = bL1_TXD_ASCT_S3[1, 1..nPredictors]'
							aux_S3 = bL1_TXD_ASCT_S3[1, cols(bL1_TXD_ASCT_S3)]
							
							// Calculate XB
							vXB_S3 = mPat_ASCT_C2 * coef_S3
							
							// Conditional RN: survive past C2
							vSurvC2 = exp(-exp(vXB_S3) :* (L1_TXD_ASCT_C2 :^ exp(aux_S3)))
							vRN_S3 = runiform(rows(idxBeyondC2), 1) :* vSurvC2
							
							// Recalculate duration with Spline 3
							vDur_S3 = calcSurvTime(vXB_S3, vRN_S3, fbL1_TXD_ASCT_S3, aux_S3)
							
							// Update durations for those beyond C2
							vDur_S1[idxBeyondC2] = vDur_S3
						}
					}
			
			// Store final durations for ASCT patients
			vOC[idxASCT] = vDur_S1
		}
	
	// Group 2: Fixed + No ASCT
		vNoASCT = (mTXR[.,1] :!= 7) :& (vSCT_DN :== 0) :& vEligible
		idxNoASCT = selectindex(vNoASCT)
		if (rows(idxNoASCT) > 0) {

			// Patient matrix
			mPat_NoASCT = (vAge[idxNoASCT], vAge2[idxNoASCT], vMale[idxNoASCT], 
					       vECOG0[idxNoASCT], vECOG1[idxNoASCT], vECOG2[idxNoASCT], 
					       vRISS1[idxNoASCT], vRISS2[idxNoASCT], vRISS3[idxNoASCT])

			if (nCR >= 1) {
				vCR1_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_TXR[1,1])
				mPat_NoASCT = (mPat_NoASCT , vCR1_NoASCT)
			}
			if (nCR >= 2) {
				vCR2_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_TXR[1,2])
				mPat_NoASCT = (mPat_NoASCT , vCR2_NoASCT)
			}
			if (nCR >= 3) {
				vCR3_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_TXR[1,3])
				mPat_NoASCT = (mPat_NoASCT , vCR3_NoASCT)	
			}	
			mPat_NoASCT = (mPat_NoASCT , vCons[idxNoASCT])
						
			// Extract coefficients 
			nPredictors = cols(mPat_NoASCT)
			coef_NoASCT = bL1_TXD_NoASCT[1, 1..nPredictors]'
			aux_NoASCT = bL1_TXD_NoASCT[1, cols(bL1_TXD_NoASCT)]

			// Calculate XB
			vXB_NoASCT = mPat_NoASCT * coef_NoASCT
			
			// Calculate survival time
			vRN_NoASCT = runiform(rows(idxNoASCT), 1)
			vOC[idxNoASCT] = calcSurvTime(vXB_NoASCT, vRN_NoASCT, fbL1_TXD_NoASCT, aux_NoASCT)
		}
				
		// Group 3: Continuous Therapy (CR == 7)
		if (anyof(oL1_TXR, 7)) {
			vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
			vCont = (mTXR[.,1] :== 7) :& vEligible
			idxCont = selectindex(vCont)
			if (rows(idxCont) > 0) {
					
				// Assemble patient matrix	
				mPat_Cont = (vAge[idxCont], vAge2[idxCont], vMale[idxCont],
							 vECOG0[idxCont], vECOG1[idxCont], vECOG2[idxCont], 
							 vRISS1[idxCont], vRISS2[idxCont], vRISS3[idxCont], 
							 vCons[idxCont])
				
				// Extract coefficients 
				nPredictors = cols(mPat_Cont)
				coef_Cont = bL1_TXD_Cont[1, 1..nPredictors]'
				dist_Cont = fbL1_TXD_Cont
				aux_Cont = bL1_TXD_Cont[1, cols(bL1_TXD_Cont)]
					
				// Calculate XB
				vXB_Cont = mPat_Cont * coef_Cont
					
				// Calculate survival time
				vRN_Cont = runiform(rows(idxCont), 1)
				vOC[idxCont] = calcSurvTime(vXB_Cont, vRN_Cont, dist_Cont, aux_Cont)
					
				// Curtail if beyond observed maximum
				vOC[idxCont] = (vOC[idxCont] :> maxL1_TXD_Cont) :* maxL1_TXD_Cont :+ (vOC[idxCont] :<= maxL1_TXD_Cont) :* vOC[idxCont]
			}
		}
		
		// Update outcome matrices
		mTXD[idx, 1] = round(vOC[idx], 0.1)
		mTNE[idx, OMC] = round(vOC[idx], 0.1)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] :+ mTNE[idx, OMC]	
	}
}
