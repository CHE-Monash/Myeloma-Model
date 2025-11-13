**********
* SIM TXR OVERRIDE - DVd L2 Method
*
* Purpose: Override Line 2 treatment to DVd (TXR=80) for all patients
*
* Context: Called AFTER standard sim_txr.do has run
*          Applies to entire DVd cohort
*          Only active when $Intervention = "DVd" AND Line = 2
*
* Author: Adam Irving
* Date: November 2025
**********

mata {
 
    // Get alive, non-prevalent patients
    idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
    if (rows(idx) > 0) {
		
		// Override Line 2 treatment to DVD
		mTXR[., 2] = J(rows(mTXR), 1, 80)
    }
}	
    