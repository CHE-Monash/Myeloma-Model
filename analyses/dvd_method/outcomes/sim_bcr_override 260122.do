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

            // mPat already assembled by sim_bcr.do
            
            // Get Stage 2 coefficients
            vCoef_stage2 = *findexternal(st_local("coef_matrix"))
            nPredictors = cols(mPat)
            vCoef = vCoef_stage2[1, 1..nPredictors]'
            
            // Calculate XB with Stage 2 coefficients
            vXB = mPat * vCoef
            
            // For DVd arm: add treatment offset
            if (st_global("int") == "dvd") {
                offset = strtoreal(st_local("offset_DVd"))
                vXB = vXB :+ offset
            }
            
            // Store for allocation functions
            st_matrix("bcr_idx", idx)
            st_matrix("bcr_vXB", vXB)
            st_numscalar("bcr_nPatients", nPatients)
		}
	}
end		

capture program drop allocate_rank
program define allocate_rank
    * Rank-based allocation using target distribution
    args bcr_matrix
	mata {	
		// Retrieve from Stata
		bcr_idx = st_matrix("bcr_idx")
		bcr_vXB = st_matrix("bcr_vXB")
		bcr_nPatients = st_numscalar("bcr_nPatients")
		
		if (bcr_nPatients > 0) {		
			// Get BCR probability matrix by name
			mBCR_probs = *findexternal(st_local("bcr_matrix"))
						
			// Get ranks (1 = worst prognosis/highest XB, nPatients = best)
			ranks = invorder(order(bcr_vXB, -1))
			percentiles = ranks / bcr_nPatients
			
			// Calculate cumulative probabilities 
			cumProbs = runningsum(mBCR_probs)
			
			// Assign BCR based on rank percentile (vectorised)
			vOC = 1 :* (percentiles :<= cumProbs[1]) + 
				  2 :* (percentiles :> cumProbs[1] :& percentiles :<= cumProbs[2]) +
				  3 :* (percentiles :> cumProbs[2] :& percentiles :<= cumProbs[3]) +
				  4 :* (percentiles :> cumProbs[3] :& percentiles :<= cumProbs[4]) +
				  5 :* (percentiles :> cumProbs[4] :& percentiles :<= cumProbs[5]) +
				  6 :* (percentiles :> cumProbs[5])
			
			// Override BCR - convert bcr_idx to column vector for indexing
			bcr_idx = bcr_idx'
			mBCR[bcr_idx, Line] = vOC
		}
	}
end

capture program drop allocate_cutpoint
program define allocate_cutpoint
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
			allocate_rank mBCR_DVd_L2
		}
		if ("$int" == "vd") {
			mata: mata matuse "$outcomes_path/A_trial/bcr_vd_l2.mmat", replace
			bcr_override_xb bL2_BCR
			allocate_rank mBCR_Vd_L2
		}
	}
	else if ("$scenario" == "B_transport") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/B_transport/transport_dvd.mmat", replace
			bcr_override_xb bL2_BCR_T
			allocate_cutpoint bL2_BCR_T
		}
	}
}
else if ($boot == 1) {
	if ("$scenario" == "A_trial") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/A_trial/bootstrap/bcr_dvd_l2_B{$b}.mmat", replace
			bcr_override_xb bL2_BCR
			allocate_rank mBCR_DVd_L2
		}
		if ("$int" == "vd") {
			mata: mata matuse "$outcomes_path/A_trial/bootstrap/bcr_vd_l2_B{$b}.mmat", replace
			bcr_override_xb bL2_BCR
			allocate_rank mBCR_Vd_L2
		}
	}
	else if ("$scenario" == "B_transport") {
		if ("$int" == "dvd") {
			mata: mata matuse "$outcomes_path/B_transport/bootstrap/transport_dvd_B{$b}.mmat", replace
			bcr_override_xb bL2_BCR_T
			allocate_cutpoint bL2_BCR_T
		}
	}
}
