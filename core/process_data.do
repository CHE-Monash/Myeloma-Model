**********
	*EpiMAP Myeloma - Process Data
**********

capture program drop process_data
program process_data

di as text "Processing Simulated Data"

*Create mSum in Mata 
	mata: mSum = vID , vMale , vECOG , vRISS , vCM , vCKD , vAge70 , vAge75 , vSCT_DN , vSCT_L1 , vMNT , /// 
			mAge , mOS , mTNE , mTSD , mMOR , mOC , mTXR , mTXD , mBCR , mTFI , mState
	
*Convert mSum to stSum
	mata: st_matrix("stSum", mSum)
	drop _all
		
*Convert stSum to variables
	svmat double stSum
	
*Name variables
	local varnames ID Male ECOGcc RISS CMc CKD Age70 Age75 SCT_DN SCT_L1 MNT ///
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
				
	local varlength : word count `varnames'
		
	forvalues i = 1/`varlength'{
		local currentvar : word `i' of `varnames'
		rename stSum`i' `currentvar'
	}		
	
	format DateDN %td
	order ID Male ECOGcc RISS CMc SCT_L1 MNT CKD
		
*Label
	label values State State_lbl
	
*Generate Dates & Years
	qui {
		gen DateL1S = DateDN + (TNE_DN * 30.4375)
		gen DateL1E = DateL1S + (TNE_L1S * 30.4375)
		gen DateL2S = DateL1E + (TNE_L1E * 30.4375)
		gen DateL2E = DateL2S + (TNE_L2S * 30.4375)
		gen DateL3S = DateL2E + (TNE_L2E * 30.4375)
		gen DateL3E = DateL3S + (TNE_L3S * 30.4375)
		gen DateL4S = DateL3E + (TNE_L3E * 30.4375)
		gen DateL4E = DateL4S + (TNE_L4S * 30.4375)
		gen DateL5S = DateL4E + (TNE_L4E * 30.4375)
		gen DateL5E = DateL5S + (TNE_L5S * 30.4375)
		gen DateL6S = DateL5E + (TNE_L5E * 30.4375)
		gen DateL6E = DateL6S + (TNE_L6S * 30.4375)
		gen DateL7S = DateL6E + (TNE_L6E * 30.4375)
		gen DateL7E = DateL7S + (TNE_L7S * 30.4375)
		gen DateL8S = DateL7E + (TNE_L7E * 30.4375)
		gen DateL8E = DateL8S + (TNE_L8S * 30.4375)
		gen DateL9S = DateL8E + (TNE_L8E * 30.4375)
		gen DateL9E = DateL9S + (TNE_L9S * 30.4375)
		gen DateSCT = DateL1E + 1 if(SCT_L1 == 1) // Fix DateSCT 1 day after DateL1E
		gen DateMOR = DateDN + (OC_TIME * 30.4375)
		format Date* %td

		gen YearDN = yofd(DateDN)
		gen YearL1 = yofd(DateL1S)
		gen YearL2 = yofd(DateL2S)
		gen YearL3 = yofd(DateL3S)
		gen YearL4 = yofd(DateL4S)
		gen YearL5 = yofd(DateL5S)
		gen YearL6 = yofd(DateL6S)
		gen YearL7 = yofd(DateL7S)
		gen YearL8 = yofd(DateL8S)
		gen YearL9 = yofd(DateL9S)
		gen YearSCT = yofd(DateSCT)
		gen YearMOR = yofd(DateMOR)
	}

* Costs

	* Define costs (AUD)
	local cVCd = 902 	// VCd: 4 21-day cycles maximum
	local cVRd = 1776 	// VRd: 5 21-day cycles maximum
	local cRd = 1608 	// Rd: 28-day cycles until progression
	local cKd = 15025 	// Kd: 28-day cycles until progression
	local cDVd = 12110 	// DVd: 28-day cycles until progression
	local cPd = 2291 	// Pd: 28-day cycles until progression
	local cOther = 4016 // Other: 28-day cycles until progression
	local cASCT = 41723 // One-time cost for stem cell transplant
	local cMNT = 1329 	// Rd / Td: 28-day cycles until progression
	local cHosp = 38743 // Hospitalisation costs per year
	local cComm = 10928 // Community care costs per year
	local cEmer = 2476 	// Emergency admission costs per year
	
	*Treatment costs

		* L1
		qui gen cTX_L1 = 0
		qui replace cTX_L1 = `cVCd' * min(4, TXD_L1 * 30.4375 / 21) if TXR_L1 == 4 // VCd (TXR = 4): Fixed-duration, max 4 × 21-day cycles
		qui replace cTX_L1 = `cVRd' * min(5, TXD_L1 * 30.4375 / 21) if TXR_L1 == 31 // VRd (TXR = 31): Fixed-duration, max 5 × 21-day cycles
		qui replace cTX_L1 = `cRd' * (TXD_L1 * 30.4375 / 28) if TXR_L1 == 7 // Rd (CR = 7): Continuous until progression, 28-day cycles
		qui replace cTX_L1 = `cOther' * (TXD_L1 * 30.4375 / 28) if TXR_L1 == 0 // Other (TXR = 0): Continuous until progression, 28-day cycles
		
		* ASCT
		qui gen cTX_ASCT = 0
		qui replace cTX_ASCT = `cASCT' if SCT_L1 == 1 

		* MNT
		qui gen cTX_MNT = 0
		qui replace cTX_MNT = `cMNT' * (TFI_L1 * 30.4375 / 28) if MNT == 1`'

		* L2
		qui gen cTX_L2 = 0
		qui replace cTX_L2 = `cDVd' * (TXD_L2 * 30.4375 / 28) if TXR_L2 == 80 // DVd (TXR = 80): Continuous until progression
		qui replace cTX_L2 = `cRd' * (TXD_L2 * 30.4375 / 28) if TXR_L2 == 7 // Rd (TXR = 7): Continuous until progression
		qui replace cTX_L2 = `cOther' * (TXD_L2 * 30.4375 / 28) if TXR_L2 == 0 // Other (TXR = 0): Continuous until progression

		* L3
		qui gen cTX_L3 = 0
		qui replace cTX_L3 = `cKd' * (TXD_L3 * 30.4375 / 28) if TXR_L3 == 49 // Kd (TXR = 49): Continuous until progression
		qui replace cTX_L3 = `cRd' * (TXD_L3 * 30.4375 / 28) if TXR_L3 == 7 // Rd (TXR = 7): Continuous until progression
		qui replace cTX_L3 = `cOther' * (TXD_L3 * 30.4375 / 28) if TXR_L3 == 0 //Other (TXR = 0): Continuous until progression

		* L4
		qui gen cTX_L4 = 0
		qui replace cTX_L4 = `cKd' * (TXD_L4 * 30.4375 / 28) if TXR_L4 == 49 // Kd (TXR = 49): Continuous until progression
		qui replace cTX_L4 = `cPd' * (TXD_L4 * 30.4375 / 28) if TXR_L4 == 56 // Pd (TXR = 56): Continuous until progression
		qui replace cTX_L4 = `cOther' * (TXD_L4 * 30.4375 / 28) if TXR_L4 == 0 // Other (TXR = 0): Continuous until progression

		* L5 - L9
		forval l = 5/9 {
			qui gen cTX_L`l' = 0
			qui replace cTX_L`l' = `cOther' * (TXD_L`l' * 30.4375 / 28) ///
				if TXD_L`l' != .
		}

		* Total undiscounted treatment cost
		qui gen cTX = cTX_L1 + cTX_ASCT + cTX_MNT + cTX_L2 + cTX_L3 + cTX_L4 + cTX_L5 + cTX_L6 + cTX_L7 + cTX_L8 + cTX_L9
	
	* Non-treatment costs (Hospital, Community and Emergency)
	qui gen cNT_Hosp = `cHosp' * (OC_TIME / 12)
	qui gen cNT_Comm = `cComm' * (OC_TIME / 12)
	qui gen cNT_Emer = `cEmer' * (OC_TIME / 12)
	qui gen cNT = cNT_Hosp + cNT_Comm + cNT_Emer
	
	*Total undiscounted cost
	qui gen cTotal = cTX + cNT
	
	* Discounting

		* Pre-calculate ln(1+r) for efficiency
		local ln_r = ln(1 + $drate)

		* Treatment Costs
		qui gen cTXd = 0
	
			* L1 to L9
			forval l = 1/9 {
				* Apply continuous discounting for uniform accrual
				qui gen cTX_L`l'd = cTX_L`l' * ((1 + $drate)^(-TSD_L`l'S/12) - (1 + $drate)^(-TSD_L`l'E/12)) / (`ln_r' * (TSD_L`l'E - TSD_L`l'S) / 12) if cTX_L`l' != . & cTX_L`l' > 0
				* Accumulate discounted costs
				qui replace cTXd = cTXd + cTX_L`l'd if cTX_L`l'd != .
			}
			
			* ASCT
			qui gen cTX_ASCTd = cTX_ASCT/(1+$drate)^((DateSCT - DateDN)/365.25) if SCT_L1 == 1
			qui replace cTXd = cTXd + cTX_ASCTd if cTX_ASCTd != .
			
			* MNT
			qui gen cTX_MNTd = cTX_MNT * ((1 + $drate)^(-TSD_L1E/12) - (1 + $drate)^(-TSD_L2S/12)) / (`ln_r' * (TSD_L2S - TSD_L1E) / 12) if cTX_MNT != . & cTX_MNT > 0
			qui replace cTXd = cTXd + cTX_MNTd if cTX_MNTd != .
		
		* Non-treatment costs
		qui gen cNT_Hospd = `cHosp' * (1 - (1 + $drate)^(-OC_TIME/12)) / `ln_r'
		qui gen cNT_Commd = `cComm' * (1 - (1 + $drate)^(-OC_TIME/12)) / `ln_r'
		qui gen cNT_Emerd = `cEmer' * (1 - (1 + $drate)^(-OC_TIME/12)) / `ln_r'
		qui gen cNTd = cNT_Hospd + cNT_Commd + cNT_Emerd
	
	* Total discounted costs
	qui gen cTotald = cTXd + cNTd	

* QALYs

	* Utility weights (from Acaster et al.)
	local uTFI = 0.72     	  // Treatment-free interval
	local uTXD_L1  = 0.63     // L1 treatment 
	local uTXD_L2  = 0.67     // L2 treatment
	local uPostL2 = 0.63	  // Post L2

	* Keep track of PreL2 time
	qui gen PreL2 = 0

	* TFI_DN
	qui gen qTFI_DN = 0
	qui replace qTFI_DN = (TFI_DN / 12) * `uTFI' if TFI_DN != .
	qui replace PreL2 = PreL2 + TFI_DN if TFI_DN != .

	* TXD_L1
	qui gen qTXD_L1 = 0
	qui replace qTXD_L1 = (TXD_L1 / 12) * `uTXD_L1' if TXD_L1 != .
	qui replace PreL2 = PreL2 + TXD_L1 if TXD_L1 != .

	* TFI_L1
	qui gen qTFI_L1 = 0
	qui replace qTFI_L1 = (TFI_L1 / 12) * `uTFI' if TFI_L1 != .
	qui replace PreL2 = PreL2 + TFI_L1 if TFI_L1 != .

	* TXD_L2 
	qui gen qTXD_L2 = 0
	qui replace qTXD_L2 = (TXD_L2 / 12) * `uTXD_L2' if TXD_L2 != .
	qui replace PreL2 = PreL2 + TXD_L2 if TXD_L2 != .

	* PostL2
	qui gen qPostL2 = 0
	qui replace qPostL2 = ((OC_TIME - PreL2) / 12) * `uPostL2' if OC_TIME > PreL2 & PreL2 < .
	drop PreL2
	
	* Calculate total undiscounted QALYs
	qui gen qTotal = qTFI_DN + qTXD_L1 + qTFI_L1 + qTXD_L2 + qPostL2

	* Discounting

		* TFI_DN
		qui gen qTFI_DNd = qTFI_DN * (1 - (1 + $drate)^(-TSD_L1S/12)) / (`ln_r' * TSD_L1S / 12) if qTFI_DN != .

		* TXD_L1
		qui gen qTXD_L1d = qTXD_L1 * ((1 + $drate)^(-TSD_L1S/12) - (1 + $drate)^(-TSD_L1E/12)) / (`ln_r' * (TSD_L1E - TSD_L1S) / 12) if qTXD_L1 != .
			
		* TFI_L1
		qui gen qTFI_L1d = qTFI_L1 * ((1 + $drate)^(-TSD_L1E/12) - (1 + $drate)^(-TSD_L2S/12)) / (`ln_r' * (TSD_L2S - TSD_L1E) / 12) if qTFI_L1 != .

		* TXD_L2
		qui gen qTXD_L2d = qTXD_L2 * ((1 + $drate)^(-TSD_L2S/12) - (1 + $drate)^(-TSD_L2E/12)) / (`ln_r' * (TSD_L2E - TSD_L2S) / 12) if qTXD_L2 != .
			
		* PostL2
		qui gen qPostL2d = qPostL2 * ((1 + $drate)^(-TSD_L2E/12) - (1 + $drate)^(-OC_TIME/12)) / (`ln_r' * (OC_TIME - TSD_L2E) / 12) if qPostL2 != .

		* Calculate total discounted QALYs
		qui gen qTotald = qTFI_DNd + qTXD_L1d + qTFI_L1d + qTXD_L2d + qPostL2d
		
end
