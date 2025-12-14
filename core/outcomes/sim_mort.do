**********
* SIM MORT
* 
* Purpose: Determine which patients have died
* Outcome: Binary (0 = Alive, 1 = Dead)
*
* Design principle: 
*   - mTNE/mTSD are NEVER modified at death (preserves audit trail of simulated times)
*   - mTFI/mTXD ARE curtailed at death (reflects actual time experienced)
*   - Use mTFI/mTXD (not mTNE/mTSD) when calculating results post-simulation
**********
mata {
	// Filters
	vStateValid = (mState[.,1] :<= OMC)
	vWasAlive = (OMC == 1 ? J(rows(mMOR), 1, 1) : (mMOR[.,OMC-1] :== 0))
	vEligible = vStateValid :& vWasAlive
	
	// Mortality determination
	// Patient dies if cumulative time at end of this OMC >= their OS
	// mTSD[.,OMC+1] should already be set (= mTSD[.,OMC] + mTNE[.,OMC])
	vDies = (mTSD[.,OMC+1] :>= mOS[.,OMC]) :& vEligible
	vLives = (mTSD[.,OMC+1] :< mOS[.,OMC]) :& vEligible
	
	// Get indices
	idxEligible = selectindex(vEligible)
	idxDies = selectindex(vDies)
	
	// Age limit: cap survival time if age at death exceeds maximum
	vExceedsLimit = ((mAge[.,1] :+ (mOS[.,OMC] / 12)) :> Limit) :& vDies
	mOS[.,OMC] = vExceedsLimit :* ((Limit :- mAge[.,1]) * 12) :+ (!vExceedsLimit) :* mOS[.,OMC]
	
	// Update mMOR - only update eligible patients
	if (rows(idxEligible) > 0) {
		vDiesElig = vDies[idxEligible]
		mMOR[idxEligible, OMC] = vDiesElig
	}
	
	// Update mOC - only for deaths
	if (rows(idxDies) > 0) {
		vOSdies = mOS[idxDies, OMC]
		mOC[idxDies, 1] = vOSdies
		mOC[idxDies, 2] = J(rows(idxDies), 1, 1)
	}
	
	// Curtail mTFI or mTXD - only for deaths
	if (rows(idxDies) > 0) {
		vTSDsafe = editmissing(mTSD[idxDies, OMC], 0)
		vOSdies = mOS[idxDies, OMC]
		vTimeMonths = rowmax((vOSdies :- vTSDsafe, J(rows(idxDies), 1, 0)))
		
		if (mod(OMC, 2) == 0) {
			// Even OMC: treatment duration (TXD)
			mTXD[idxDies, OMC/2] = rowmin((mTXD[idxDies, OMC/2], vTimeMonths))
		}
		else {
			// Odd OMC: treatment-free interval (TFI)
			mTFI[idxDies, (OMC+1)/2] = rowmin((mTFI[idxDies, (OMC+1)/2], vTimeMonths))
		}
	}
	
	// NOTE: mTNE and mTSD are NOT modified at death
	// This preserves the full audit trail of simulated times
	// For results, use mTFI/mTXD (actual experienced durations)
}
