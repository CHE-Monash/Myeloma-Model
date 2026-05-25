**********
* SIM TXR Override - TD L4 Pre
*
* Purpose: Override Line 2 treatment to Td (TXR=999)
*
* Context: Called AFTER standard sim_txr.do has run
*          Applies to entire cohort
*          Only active when $analysis = "dvd-l2-method" AND Line = 2
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
			mTXR[idx, $line] = J(rows(idx), 1, strtoreal(st_local("tx_code")))
		}
	}
end

// Call function based on $int
if ("${int}" == "td") {
	txr_override 999
}
