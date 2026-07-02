**********
* Monash Myeloma Model - Process Data
*
* Generalised version that processes simulations starting at line $line
* - If $line == 1: Full pathway from diagnosis, calendar-date discounting
* - If $line > 1: From line $line onwards, relative-time discounting (L$line start = time 0)
*
* Naming Convention:
*   cost_{component}[_L#][_d]  - Cost variables (tx, nt, total)
*   qaly_{component}[_d]       - QALY variables
*   _d suffix = discounted
**********

capture program drop process_data
program process_data

local L = $line
if `L' == 0 local L = 1   // Line 0 = full pathway from diagnosis
local maxL = 9

di as text "Processing Simulated Data (Starting Line: `L')"

* Free the CRN matrix - all draws are done by end of simulation; not needed here,
*   and process_data is the peak-memory stage (mSum assembly + getmata).
cap mata: mata drop mRN

* Create mSum in Mata
	mata: mSum = vID , vMale , vECOG , vRISS , vISS , vCM , vCKD , vAge70 , vAge75 , vSCT_DN , vSCT_L1 , vMNT , /// 
			mAge , mOS , mTNE , mTSD , mMOR , mOC , mTXR , mTXD , mBCR , mTFI , mState
	
* Column names for mSum, in assembly order below.
* (getmata errors on a name/column count mismatch, which guards this alignment.)
	local varnames ID Male ECOGcc RISS ISS CMc CM_CKD Age70 Age75 SCT_DN SCT_L1 MNT ///
		Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E ///
		OS_DN OS_L1S OS_L1E OS_L2S OS_L2E OS_L3S OS_L3E OS_L4S OS_L4E OS_L5S OS_L5E OS_L6S OS_L6E OS_L7S OS_L7E OS_L8S OS_L8E OS_L9S OS_L9E ///
		TNE_DN TNE_L1S TNE_L1E TNE_L2S TNE_L2E TNE_L3S TNE_L3E TNE_L4S TNE_L4E TNE_L5S TNE_L5E TNE_L6S TNE_L6E TNE_L7S TNE_L7E TNE_L8S TNE_L8E TNE_L9S TNE_L9E ///
		TSD_DN TSD_L1S TSD_L1E TSD_L2S TSD_L2E TSD_L3S TSD_L3E TSD_L4S TSD_L4E TSD_L5S TSD_L5E TSD_L6S TSD_L6E TSD_L7S TSD_L7E TSD_L8S TSD_L8E TSD_L9S TSD_L9E ///
		MOR_DN MOR_L1S MOR_L1E MOR_L2S MOR_L2E MOR_L3S MOR_L3E MOR_L4S MOR_L4E MOR_L5S MOR_L5E MOR_L6S MOR_L6E MOR_L7S MOR_L7E MOR_L8S MOR_L8E MOR_L9S MOR_L9E ///
		OC_TIME OC_MORT ///
		TXR_L1 TXR_L2 TXR_L3 TXR_L4 TXR_L5 TXR_L6 TXR_L7 TXR_L8 TXR_L9 ///
		TXD_L1 TXD_L2 TXD_L3 TXD_L4 TXD_L5 TXD_L6 TXD_L7 TXD_L8 TXD_L9 ///
		BCR_L1 BCR_L2 BCR_L3 BCR_L4 BCR_L5 BCR_L6 BCR_L7 BCR_L8 BCR_L9 BCR_SCT ///
		TFI_DN TFI_L1 TFI_L2 TFI_L3 TFI_L4 TFI_L5 TFI_L6 TFI_L7 TFI_L8 /// 
		State DateDN ///
				

* Write the Mata matrix straight to named variables. getmata reads mSum
* directly, bypassing the st_matrix()/svmat round-trip that dominated runtime
* (st_matrix on the full ~101k x 148 matrix was ~60s; getmata is < 1s).
	drop _all
	getmata (`varnames') = mSum, double
	cap mata: mata drop mSum
	
	format DateDN %td
	order ID Male ECOGcc RISS CMc SCT_L1 MNT CM_CKD
		
* Label
	label values State State_lbl

**********
* Time Variables
**********

* Survival time from starting line to death (months)
	if `L' == 1 {
		qui gen OC_TIME_L = OC_TIME
	}
	else {
		qui gen OC_TIME_L = OC_TIME - TSD_L`L'S
	}

* Relative time markers (starting line = time 0)
	qui {
		if `L' == 1 {
			gen TSD_L1S_ref = TFI_DN if TFI_DN != .
			replace TSD_L1S_ref = 0 if TFI_DN == .
			gen TSD_L1E_ref = TSD_L1S_ref + TXD_L1 if TXD_L1 != .
			local first_l = 2
		}
		else if `L' > 1 {
			gen TSD_L`L'S_ref = 0
			gen TSD_L`L'E_ref = TXD_L`L' if TXD_L`L' != .
			local first_l = `=`L'+1'
		}
		local prev = `L'
		forval l = `first_l'/`maxL' {
			local tfi_idx = `=`l'-1'
			gen TSD_L`l'S_ref = TSD_L`prev'E_ref + TFI_L`tfi_idx' if TFI_L`tfi_idx' != .
			gen TSD_L`l'E_ref = TSD_L`l'S_ref + TXD_L`l' if TXD_L`l' != .
			local prev = `l'
		}
	}

* Generate Dates & Years from L onwards
	qui {
		forval l = `L'/`maxL' {
			gen DateL`l'S = DateDN + (TSD_L`l'S * 30.4375)
			gen DateL`l'E = DateDN + (TSD_L`l'E * 30.4375)
		}
		
		if `L' == 1 {
			gen DateSCT = DateL1E + 1 if SCT_L1 == 1
			gen YearSCT = yofd(DateSCT)
		}
		gen DateMOR = DateDN + (OC_TIME * 30.4375)
		format Date* %td
		
		gen YearDN = yofd(DateDN)
		forval l = `L'/`maxL' {
			gen YearL`l' = yofd(DateL`l'S)
		}
		gen YearMOR = yofd(DateMOR)
	}

**********
* Costs
**********

	local cVCd = 902
	local cVRd = 1776
	local cRd = 1608
	local cKd = 15025
	local cDVd = 12110
	local cPd = 2291
	local cVd = 724
	local cOther = 1612 // VTd / TCd / Td / Vd
	
	local cASCT = 41723
	local cMNT = 1329
	
	local cHosp = 38743
	local cComm = 10928
	local cEmer = 2476
	
	local ln_r = ln(1 + $drate)

* Treatment costs by line (undiscounted)
	forval l = `L'/`maxL' {
		qui gen cost_tx_L`l' = 0
		
		* L1-specific regimens with cycle caps
		if `l' == 1 {
			qui replace cost_tx_L1 = `cVCd' * min(4, TXD_L1 * 30.4375 / 21) if TXR_L1 == 4
			qui replace cost_tx_L1 = `cVRd' * min(5, TXD_L1 * 30.4375 / 21) if TXR_L1 == 31
		}
		
		* Standard regimens (all lines)
		qui replace cost_tx_L`l' = `cRd' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 7
		qui replace cost_tx_L`l' = `cKd' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 49
		qui replace cost_tx_L`l' = `cDVd' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 80
		qui replace cost_tx_L`l' = `cPd' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 56
		qui replace cost_tx_L`l' = `cVd' * min(8, TXD_L2 * 30.4375 / 21) if TXR_L`l' == 5
		qui replace cost_tx_L`l' = `cOther' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 0
	}

* L1-specific costs (ASCT, MNT)
	qui gen cost_tx_asct = 0
	qui gen cost_tx_mnt = 0
	if `L' == 1 {
		qui replace cost_tx_asct = `cASCT' if SCT_L1 == 1
		qui replace cost_tx_mnt = `cMNT' * (TFI_L1 * 30.4375 / 28) if MNT == 1
	}

* Total undiscounted treatment cost
	qui gen cost_tx = cost_tx_asct + cost_tx_mnt
	forval l = `L'/`maxL' {
		qui replace cost_tx = cost_tx + cost_tx_L`l'
	}

* Non-treatment costs (undiscounted)
	qui gen cost_nt = (`cHosp' + `cComm' + `cEmer') * (OC_TIME_L / 12)

* Total undiscounted cost
	qui gen cost_total = cost_tx + cost_nt

* Discounted treatment costs
	qui gen cost_tx_d = 0
	
	* First line (starts at time 0 or after TFI_DN for L1)
	if `L' == 1 {
		qui gen cost_tx_L1_d = cost_tx_L1 * ((1 + $drate)^(-TSD_L1S_ref/12) - (1 + $drate)^(-TSD_L1E_ref/12)) / (`ln_r' * TXD_L1/12) if cost_tx_L1 > 0 & TXD_L1 > 0
		local first_l = 2
	}
	else {
		qui gen cost_tx_L`L'_d = cost_tx_L`L' * (1 - (1 + $drate)^(-TSD_L`L'E_ref/12)) / (`ln_r' * TXD_L`L'/12) if cost_tx_L`L' > 0 & TXD_L`L' > 0
		local first_l = `=`L'+1'
	}
	qui replace cost_tx_d = cost_tx_d + cost_tx_L`L'_d if cost_tx_L`L'_d != .
	
	* Subsequent lines
	forval l = `first_l'/`maxL' {
		qui gen cost_tx_L`l'_d = cost_tx_L`l' * ((1 + $drate)^(-TSD_L`l'S_ref/12) - (1 + $drate)^(-TSD_L`l'E_ref/12)) / (`ln_r' * TXD_L`l'/12) if cost_tx_L`l' > 0 & TXD_L`l' > 0
		qui replace cost_tx_d = cost_tx_d + cost_tx_L`l'_d if cost_tx_L`l'_d != .
	}
	
	* ASCT and MNT (L1 only)
	if `L' == 1 {
		qui gen cost_tx_asct_d = cost_tx_asct * (1 + $drate)^(-TSD_L1E_ref/12) if SCT_L1 == 1
		qui replace cost_tx_d = cost_tx_d + cost_tx_asct_d if cost_tx_asct_d != .
		
		qui gen cost_tx_mnt_d = cost_tx_mnt * ((1 + $drate)^(-TSD_L1E_ref/12) - (1 + $drate)^(-TSD_L2S_ref/12)) / (`ln_r' * TFI_L1/12) if cost_tx_mnt > 0 & TFI_L1 > 0
		qui replace cost_tx_d = cost_tx_d + cost_tx_mnt_d if cost_tx_mnt_d != .
	}
	
	* Non-treatment costs (discounted)
	qui gen cost_nt_d = (`cHosp' + `cComm' + `cEmer') * (1 - (1 + $drate)^(-OC_TIME_L/12)) / `ln_r'

* Total discounted costs
	qui gen cost_total_d = cost_tx_d + cost_nt_d

**********
* QALYs
**********

* L1 analysis
	if `L' == 1 {
		
		local uTFI = 0.72
		local uTXD_L1 = 0.63
		local uTXD_L2 = 0.67
		local uPostL2 = 0.63
		
		qui gen PreL2 = 0
		
		* TFI_DN
		qui gen qaly_tfi_DN = 0
		qui replace qaly_tfi_DN = (TFI_DN / 12) * `uTFI' if TFI_DN != .
		qui replace PreL2 = PreL2 + TFI_DN if TFI_DN != .
		
		* TXD_L1
		qui gen qaly_txd_L1 = 0
		qui replace qaly_txd_L1 = (TXD_L1 / 12) * `uTXD_L1' if TXD_L1 != .
		qui replace PreL2 = PreL2 + TXD_L1 if TXD_L1 != .
		
		* TFI_L1
		qui gen qaly_tfi_L1 = 0
		qui replace qaly_tfi_L1 = (TFI_L1 / 12) * `uTFI' if TFI_L1 != .
		qui replace PreL2 = PreL2 + TFI_L1 if TFI_L1 != .
		
		* TXD_L2
		qui gen qaly_txd_L2 = 0
		qui replace qaly_txd_L2 = (TXD_L2 / 12) * `uTXD_L2' if TXD_L2 != .
		qui replace PreL2 = PreL2 + TXD_L2 if TXD_L2 != .
		
		* PostL2
		qui gen qaly_post_L2 = 0
		qui replace qaly_post_L2 = ((OC_TIME_L - PreL2) / 12) * `uPostL2' if OC_TIME_L > PreL2 & PreL2 < .
		drop PreL2
		
		* Total undiscounted QALYs
		qui gen qaly_total = qaly_tfi_DN + qaly_txd_L1 + qaly_tfi_L1 + qaly_txd_L2 + qaly_post_L2
		
		* Discounting (using _ref time markers)
		
		* TFI_DN: from time 0 to TSD_L1S_ref
		qui gen qaly_tfi_DN_d = 0
		qui replace qaly_tfi_DN_d = qaly_tfi_DN * (1 - (1 + $drate)^(-TSD_L1S_ref/12)) / (`ln_r' * TSD_L1S_ref/12) if TSD_L1S_ref > 0 & qaly_tfi_DN > 0
		
		* TXD_L1: from TSD_L1S_ref to TSD_L1E_ref
		qui gen qaly_txd_L1_d = 0
		qui replace qaly_txd_L1_d = qaly_txd_L1 * ((1 + $drate)^(-TSD_L1S_ref/12) - (1 + $drate)^(-TSD_L1E_ref/12)) / (`ln_r' * TXD_L1/12) if TXD_L1 > 0 & qaly_txd_L1 > 0
		
		* TFI_L1: from TSD_L1E_ref to TSD_L2S_ref
		qui gen qaly_tfi_L1_d = 0
		qui replace qaly_tfi_L1_d = qaly_tfi_L1 * ((1 + $drate)^(-TSD_L1E_ref/12) - (1 + $drate)^(-TSD_L2S_ref/12)) / (`ln_r' * TFI_L1/12) if TFI_L1 > 0 & qaly_tfi_L1 > 0
		
		* TXD_L2: from TSD_L2S_ref to TSD_L2E_ref
		qui gen qaly_txd_L2_d = 0
		qui replace qaly_txd_L2_d = qaly_txd_L2 * ((1 + $drate)^(-TSD_L2S_ref/12) - (1 + $drate)^(-TSD_L2E_ref/12)) / (`ln_r' * TXD_L2/12) if TXD_L2 > 0 & qaly_txd_L2 > 0
		
		* PostL2: from TSD_L2E_ref to OC_TIME_L
		qui gen qaly_post_L2_d = 0
		qui replace qaly_post_L2_d = qaly_post_L2 * ((1 + $drate)^(-TSD_L2E_ref/12) - (1 + $drate)^(-OC_TIME_L/12)) / (`ln_r' * (OC_TIME_L - TSD_L2E_ref)/12) if qaly_post_L2 > 0 & OC_TIME_L > TSD_L2E_ref
		
		* Total discounted QALYs
		qui gen qaly_total_d = qaly_tfi_DN_d + qaly_txd_L1_d + qaly_tfi_L1_d + qaly_txd_L2_d + qaly_post_L2_d
	}
* L2 analysis
	else if `L' == 2 {
		
		local uTXD_L2 = 0.67
		local uPostL2 = 0.63
		
		* TXD_L2
		qui gen qaly_txd_L2 = (TXD_L2 / 12) * `uTXD_L2' if TXD_L2 != .
		qui replace qaly_txd_L2 = 0 if qaly_txd_L2 == .
		
		* Post L2
		qui gen qaly_post_L2 = ((OC_TIME_L - TSD_L2E_ref) / 12) * `uPostL2' if OC_TIME_L > TSD_L2E_ref
		qui replace qaly_post_L2 = 0 if qaly_post_L2 == .
		
		* Total undiscounted QALYs
		qui gen qaly_total = qaly_txd_L2 + qaly_post_L2
		
		* Discounting (TSD_L2S_ref = 0 by construction)
		
		* TXD_L2: from TSD_L2S_ref (=0) to TSD_L2E_ref
		qui gen qaly_txd_L2_d = qaly_txd_L2 * ((1 + $drate)^(-TSD_L2S_ref/12) - (1 + $drate)^(-TSD_L2E_ref/12)) / (`ln_r' * TXD_L2/12) if TXD_L2 > 0 & qaly_txd_L2 > 0
		qui replace qaly_txd_L2_d = 0 if qaly_txd_L2_d == .
		
		* Post: from TSD_L2E_ref to OC_TIME_L
		qui gen qaly_post_L2_d = qaly_post_L2 * ((1 + $drate)^(-TSD_L2E_ref/12) - (1 + $drate)^(-OC_TIME_L/12)) / (`ln_r' * (OC_TIME_L - TSD_L2E_ref)/12) if qaly_post_L2 > 0 & OC_TIME_L > TSD_L2E_ref
		qui replace qaly_post_L2_d = 0 if qaly_post_L2_d == .
		
		* Total discounted QALYs
		qui gen qaly_total_d = qaly_txd_L2_d + qaly_post_L2_d
	}
* L3+ analysis
	else {
		* L >= 3: Everything from line L start uses uPostL2
		local uPostL2 = 0.63
		
		* Total time from start of line L (TSD_L`L'S_ref = 0) to outcome
		qui gen qaly_total = (OC_TIME_L / 12) * `uPostL2' if OC_TIME_L > 0 & OC_TIME_L != .
		qui replace qaly_total = 0 if qaly_total == .
		
		* Discounting: from TSD_L`L'S_ref (=0) to OC_TIME_L
		qui gen qaly_total_d = qaly_total * ((1 + $drate)^(-TSD_L`L'S_ref/12) - (1 + $drate)^(-OC_TIME_L/12)) / (`ln_r' * OC_TIME_L/12) if qaly_total > 0 & OC_TIME_L > 0
		qui replace qaly_total_d = 0 if qaly_total_d == .
	}

end
