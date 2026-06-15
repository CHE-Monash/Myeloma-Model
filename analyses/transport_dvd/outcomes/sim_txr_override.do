**********
* SIM TXR Override - DVd L2 Method
*
* Purpose: Override Line 2 treatment to DVd (TXR=80) or Vd (TXR=5) for all patients
*
* Context: Called AFTER standard sim_txr.do has run
*          Applies to entire DVd cohort
*          Only active when $analysis = "transport_dvd" AND Line = 2
*
* Author: Adam Irving
* Date: November 2025
**********

capture program drop txr_override 
program define txr_override
	args tx_code
	mata {
		// Filter for alive and eligible
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
		if (rows(idx) > 0) {
			
			// Override Line 2 treatment 
			mTXR[idx, 2] = J(rows(idx), 1, strtoreal(st_local("tx_code")))
		}
	}

end

// Call function based on $int
if ("$int" == "dvd") {
	txr_override 80
}
if ("$int" == "vd") {
	txr_override 5
}

