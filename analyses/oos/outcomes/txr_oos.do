**********
* TXR
*
* Out-of-sample (70/30) validation -- regimen definitions identical to the base model
* (the OOS analysis validates the base-model structure on held-out patients).
**********

	*2 - Thal/Cycl/Dexa
	*4 - Bort/Cycl/Dexa
	*7 - Lena/Dexa
	*9 - Bort/Thal/Dexa
	*31 - Bort/Lena/Dexa
	*49 - Carf/Dexa
	*56 - Poma/Dexa
	*80 - Dara/Bort/Dexa
/*
	cap program drop regimens
	program define regimens
		foreach Event of local 0 { // Arguments automaticaly stored in local 0
			di _n "Event = `Event'"
			qui preserve
			qui contract Regimen if Event0 == `Event' & yofd(Date0) >= 2020 & yofd(Date0) <= 2025, freq(n) percent(pct)
			qui gsort -n
			qui format pct %4.1f
			list in 1/5, clean noobs
			qui restore
		}
	end

	regimens 10 20 30 40
*/

	*L1
		gen TXR_L1 = 0 if Event0 == 10
		local TXR_L1 4 31
		foreach l of local TXR_L1 {
			replace TXR_L1 = `l' if Event0 == 10 & Regimen == `l'
		}
		bysort ID (TXR_L1): replace TXR_L1 = TXR_L1[_n-1] if TXR_L1 == .

	*L2
		gen TXR_L2 = 0 if Event0 == 20
		local TXR_L2 7 80
		foreach l of local TXR_L2 {
			replace TXR_L2 = `l' if Event0 == 20 & Regimen == `l'
		}
		bysort ID (TXR_L2): replace TXR_L2 = TXR_L2[_n-1] if TXR_L2 == .

	*L3
		gen TXR_L3 = 0 if Event0 == 30
		local TXR_L3 49 7
		foreach l of local TXR_L3 {
			replace TXR_L3 = `l' if Event0 == 30 & Regimen == `l'
		}
		bysort ID (TXR_L3): replace TXR_L3 = TXR_L3[_n-1] if TXR_L3 == .

	*L4
		gen TXR_L4 = 0 if Event0 == 40
		local TXR_L4 49 56
		foreach l of local TXR_L4 {
			replace TXR_L4 = `l' if Event0 == 40 & Regimen == `l'
		}
		bysort ID (TXR_L4): replace TXR_L4 = TXR_L4[_n-1] if TXR_L4 == .
