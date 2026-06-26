**********
	* SIM TXR - VRd Post (regimen definitions)
**********

	*2 - Thal/Cycl/Dexa
	*4 - Bort/Cycl/Dexa
	*7 - Lena/Dexa
	*9 - Bort/Thal/Dexa
	*31 - Bort/Lena/Dexa
	*49 - Carf/Dexa
	*56 - Poma/Dexa
	*80 - Dara/Bort/Dexa
	
	*tab Regimen if Event0 == 20 & yofd(Date0) >= 2019 & yofd(Date0) <= 2024, sort
	
	*L1 Chemo Regimens
		gen TXR_L1 = 0 if(Event0 == 10)
		local TXR_L1 4 7 31
		foreach l of local TXR_L1 {		
			replace TXR_L1 = `l' if(Event0 == 10 & Regimen == `l')
		}
		bysort ID (TXR_L1): replace TXR_L1 = TXR_L1[_n-1] if(TXR_L1 == .)
		
	*L2 Chemo Regimens
		gen TXR_L2 = 0 if(Event0 == 20)
		local TXR_L2 7 80 
		foreach l of local TXR_L2 {
			replace TXR_L2 = `l' if(Event0 == 20 & Regimen == `l')
		}
		bysort ID (TXR_L2): replace TXR_L2 = TXR_L2[_n-1] if(TXR_L2 == .)
		
	*L3 Chemo Regimens
		gen TXR_L3 = 0 if(Event0 == 30)
		local TXR_L3 7 49 
		foreach l of local TXR_L3 {
			replace TXR_L3 = `l' if(Event0 == 30 & Regimen == `l')
		}
		bysort ID (TXR_L3): replace TXR_L3 = TXR_L3[_n-1] if(TXR_L3 == .)
		
	*L4 Chemo Regimens
		gen TXR_L4 = 0 if(Event0 == 40)
		local TXR_L4 7 49
		foreach l of local TXR_L4 {
			replace TXR_L4 = `l' if(Event0 == 40 & Regimen == `l')
		}
		bysort ID (TXR_L4): replace TXR_L4 = TXR_L4[_n-1] if(TXR_L4 == .)

