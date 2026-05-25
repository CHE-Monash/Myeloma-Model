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
    // Filter for alive and eligible
    idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
    if (rows(idx) > 0) {
		
		// Override Line 2 treatment to DVd
		mTXR[idx, 2] = J(rows(idx), 1, 80)
    }
}
