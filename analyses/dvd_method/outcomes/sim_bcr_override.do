**********
* SIM BCR OVERRIDE - DVd Method
*
* Purpose: Override BCR for DVd/Vd at Line 2 
* Method: Called after sim_bcr.do
*
* Author: Adam Irving
* Date: January 2026
**********

capture program drop A_trial
program define A_trial
    args bcr_matrix
	mata {	
		// Get alive, non-prevalent patients
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))		
		if (rows(idx) > 0) {
			
			nPatients = rows(idx)

            // mPat already assembled by sim_bcr.do		
			// Get BCR probability matrix by name
			mBCR_probs = *findexternal(st_local("bcr_matrix"))
						
			// Get ranks (1 = worst prognosis/highest XB, nPatients = best)
			ranks = invorder(order(vXB, -1))
			percentiles = ranks / nPatients
			
			// Calculate cumulative probabilities 
			cumProbs = runningsum(mBCR_probs)
			
			// Assign BCR based on rank percentile (vectorised)
			vOC = 1 :* (percentiles :<= cumProbs[1]) + 
				  2 :* (percentiles :> cumProbs[1] :& percentiles :<= cumProbs[2]) +
				  3 :* (percentiles :> cumProbs[2] :& percentiles :<= cumProbs[3]) +
				  4 :* (percentiles :> cumProbs[3] :& percentiles :<= cumProbs[4]) +
				  5 :* (percentiles :> cumProbs[4] :& percentiles :<= cumProbs[5]) +
				  6 :* (percentiles :> cumProbs[5])
			
			// Override BCR
			mBCR[idx, Line] = vOC
		}
	}
end

capture program drop B_transport
program define B_transport
    mata {	
		// Get alive, non-prevalent patients
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))		
		if (rows(idx) > 0) {
			
			nPatients = rows(idx)

            // Grab mPat assembled by sim_bcr.do and remove TXR columns
            //mPat_T = mPat[., 1..(cols(mPat) - cols(oL2_TXR))]
			mPat_T = J(nPatients, 2, 1)

            // Get Stage 2 coefficients
            nPredictors = cols(mPat_T)
            vCoef = bL2_BCR_T1[1, 1..nPredictors]'
            
            // Calculate XB with Stage 1 coefficients
            vXB_T = mPat_T * vCoef
            
            // Extract Stage 1 cutpoints
            nCutPoints = 5
            cutPointIndices = (cols(bL2_BCR_T1) - nCutPoints + 1)..cols(bL2_BCR_T1)
			cutPoints = bL2_BCR_T1[1, cutPointIndices]
			           
            // Calculate probabilities from ordered logit
            cumProbs = calcOrdLogitProbs(vXB_T, cutPoints)
            
            // Stochastic assignment
            vRN = runiform(nPatients, 1)
			categoryValues = (1, 2, 3, 4, 5, 6)
            vOC = assignOrdOutcome(vRN, cumProbs, categoryValues)
            
			// Override BCR
            mBCR[idx, Line] = vOC
        }
    }
end

// Call functions based on $boot, $int & $scenario
if ($boot == 0) {
	if ("$scenario" == "A_trial") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/A_trial/bcr_dvd_l2.mmat", replace
			A_trial mBCR_DVd_L2
		}
		if ("$int" == "vd") {
			mata: mata matuse "$outcomes_path/A_trial/bcr_vd_l2.mmat", replace
			A_trial mBCR_Vd_L2
		}
	}
	else if ("$scenario" == "B_transport") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/B_transport/transport_dvd.mmat", replace
			B_transport
		}
	}
}
else if ($boot == 1) {
	if ("$scenario" == "A_trial") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/A_trial/bootstrap/bcr_dvd_l2_B{$b}.mmat", replace
			A_trial mBCR_DVd_L2
		}
		if ("$int" == "vd") {
			mata: mata matuse "$outcomes_path/A_trial/bootstrap/bcr_vd_l2_B{$b}.mmat", replace
			A_trial mBCR_Vd_L2
		}
	}
	else if ("$scenario" == "B_transport") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/B_transport/bootstrap/transport_dvd_B{$b}.mmat", replace
			B_transport
		}
	}
}
