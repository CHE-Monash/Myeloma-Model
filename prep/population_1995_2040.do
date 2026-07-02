**********
* Monash Myeloma Model - Population 1995 to 2040
*
* Purpose: This code generates a synthetic cohort of MM patients to be simultaed
* Timeframes:
*		- 1995 to 2020 from AIHW data 
*		- 2021 to 2040 from Daffodil Centre data
* Outcome: population_1995_2040_(1 to 10).dta
**********

clear
set more off
	
if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $data_path (git-ignored)

**********	
*Process function
capture program drop population
program define population

	*Open Daffodil Centre data
		import delimited "patients/population_forecast.csv", case(preserve) clear
		
	*Globals
		global AgeGroups 1 2 3 4 5 6 7 8 9
		global Years 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030 2031 2032 2033 2034 2035 2036 2037 2038 2039 2040
		
	*IncidentMales1 (2021 to 2040)
		matrix def IncidentMales1 = J(9,20,.)
		matrix colnames IncidentMales1 = $Years
		matrix rownames IncidentMales1 = $AgeGroups
		foreach i of global AgeGroups {
			foreach j of global Years {
				qui gen temp = Incidence if Sex == "Male" & AgeGroup == `i' & Year == `j'
				qui sum temp
				matrix IncidentMales1[`i',`j'-2020] = `r(mean)'
				qui drop temp
			}
		}
		matlist IncidentMales1
		
	*IncidentFemales (2021 to 2040)
		matrix def IncidentFemales1 = IncidentMales1
		foreach i of global AgeGroups {
			foreach j of global Years {
				qui gen temp = Incidence if Sex == "Female" & AgeGroup == `i' & Year == `j'
				qui sum temp
				matrix IncidentFemales1[`i',`j'-2020] = `r(mean)'
				qui drop temp
			}
		}
		matlist IncidentFemales1

	*Open AIHW historical data
		import delimited "patients/population_historical.csv", case(preserve) clear
		
	*Globals
		global AgeGroups 1 2 3 4 5 6 7 8 9
		global Years 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020

	*IncidentMales2 (1995 to 2020)
		matrix def IncidentMales2 = J(9,26,.)
		matrix colnames IncidentMales2 = $Years
		matrix rownames IncidentMales2 = $AgeGroups
		foreach i of global AgeGroups {
			foreach j of global Years {
				qui gen temp = Incidence if Sex == "Male" & AgeGroup == `i' & Year == `j'
				qui sum temp
				matrix IncidentMales2[`i',`j'-1994] = `r(mean)'
				qui drop temp
			}
		}
		matlist IncidentMales2
		
	*IncidentFemales2 (1995 to 2020)
		matrix def IncidentFemales2 = IncidentMales2
		foreach i of global AgeGroups {
			foreach j of global Years {
				qui gen temp = Incidence if Sex == "Female" & AgeGroup == `i' & Year == `j'
				qui sum temp
				matrix IncidentMales2[`i',`j'-1994] = `r(mean)'
				qui drop temp
			}
		}
		matlist IncidentFemales2
	
	*Combine Incident1 and Incident2
		matrix def IncidentMales = IncidentMales2 , IncidentMales1
		matlist IncidentMales
		
		matrix def IncidentFemales = IncidentFemales2 , IncidentFemales1
		matlist IncidentFemales

	*Create IncidentMales and IncidentFemales in each AgeGroup and Year
		global Years 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026 2027 2028 2029 2030 2031 2032 2033 2034 2035 2036 2037 2038 2039 2040
		drop _all
		gen Male = .
		gen AgeGroup = .
		gen Year = .
		foreach a of global AgeGroups {
			foreach y of global Years {
				qui sum Year
				scalar obs = r(N) + IncidentMales[`a', `y' - 1994]
				set obs `=obs' 
				replace Male = 1 if(Male == .)
				replace AgeGroup = `a' if(AgeGroup == .)
				replace Year = `y' if(Year == .)
				qui sum Year
				scalar obs = r(N) + IncidentFemales[`a', `y' - 1994]
				set obs `=obs' 
				replace Male = 0 if(Male == .)
				replace AgeGroup = `a' if(AgeGroup == .)
				replace Year = `y' if(Year == .)
			}
		}
		
		*Append and impute remaining variables
			append using "${data_path}/MRDR Wide MI.dta"
			order MRDR ID Year AgeGroup Age Male ECOGcc RISS ISS LDHRisk FISHRisk CMc
			replace ID = _n if ID == .
			
		*AgeGroup
			qui replace AgeGroup = 1 if Age >= 20 & Age < 50
			qui replace AgeGroup = 2 if Age >= 50 & Age < 55
			qui replace AgeGroup = 3 if Age >= 55 & Age < 60
			qui replace AgeGroup = 4 if Age >= 60 & Age < 65
			qui replace AgeGroup = 5 if Age >= 65 & Age < 70
			qui replace AgeGroup = 6 if Age >= 70 & Age < 75
			qui replace AgeGroup = 7 if Age >= 75 & Age < 80
			qui replace AgeGroup = 8 if Age >= 80 & Age < 85
			qui replace AgeGroup = 9 if Age >= 85 & Age < 100
			
		*Age bounds
			gen Age_lower = 21 if AgeGroup == 1
			gen Age_upper = 49.9 if AgeGroup == 1
			qui replace Age_lower = 50 if AgeGroup == 2
			qui replace Age_upper = 54.9 if AgeGroup == 2
			qui replace Age_lower = 55 if AgeGroup == 3
			qui replace Age_upper = 59.9 if AgeGroup == 3
			qui replace Age_lower = 60 if AgeGroup == 4
			qui replace Age_upper = 64.9 if AgeGroup == 4
			qui replace Age_lower = 65 if AgeGroup == 5
			qui replace Age_upper = 69.9 if AgeGroup == 5
			qui replace Age_lower = 70 if AgeGroup == 6
			qui replace Age_upper = 74.9 if AgeGroup == 6
			qui replace Age_lower = 75 if AgeGroup == 7
			qui replace Age_upper = 79.9 if AgeGroup == 7
			qui replace Age_lower = 80 if AgeGroup == 8
			qui replace Age_upper = 84.9 if AgeGroup == 8
			qui replace Age_lower = 85 if AgeGroup == 9
			qui replace Age_upper = 100 if AgeGroup == 9	
			
		*MI set
			mi set wide
			*mi register imputed Age ECOGcc RISS CMc CM_CKD
			mi register imputed Age ECOGcc ISS LDHRisk FISHRisk CM_CKD CM_CRD CM_PLM CM_DBT CM_LVR CM_PNR CM_MLG
			mi register regular AgeGroup Male Age_lower Age_upper
			mi describe
			
		*Imputation
			*local seed = 9854 * `Sample'
			mi impute chained (truncreg, ll(Age_lower) ul(Age_upper)) Age (ologit) ECOGcc ISS (logit) LDHRisk FISHRisk CM_CKD CM_CRD CM_PLM CM_DBT CM_LVR CM_PNR CM_MLG = AgeGroup Male, add(1) rseed(9854)
			mi unset
			
		*Drop MRDR and reshape Synthetic
			drop if MRDR != .
			drop mi_miss
			foreach var of varlist _all {
				qui count if !missing(`var')
				if r(N) == 0 {
					drop `var'
				}
			}
		
		*Remove underscores in varnames	
			foreach var of varlist _all {
				local newname = subinstr("`var'", "_1_", "", .)
				if "`var'" != "`newname'" {
					rename `var' `newname'
				}
			}
			
		*Create RISS
			gen RISS = 1 if ISS == 1 & LDHRisk == 0 & FISHRisk == 0
			replace RISS = 3 if ISS == 3 & (LDHRisk == 1 | FISHRisk == 1)
			replace RISS = 2 if RISS == .
			
		*Create CM Score
			gen CM = CM_CKD + CM_CRD + CM_PLM + CM_DBT + CM_LVR + CM_PNR + CM_MLG
			gen CMc = CM
			replace CMc = 3 if CMc == 4 | CMc == 5 | CMc == 6
			
		*Allocate DateDN 
			sort Year AgeGroup Male
			replace ID = _n	// Create ID sorted on YearN
			gen DateDN = .
			foreach y of global Years {
				di `y'
				
			*Scalars for how many patients to allocate per day
				qui sum ID if(Year == `y')
				scalar D_`y' = r(N)/364 // Using 364 so no patients slip into the following year
				scalar n_`y' = `=r(min)'
				scalar N_`y' = `=r(max)'
				
			*Locals for 1 Jan in each Year
				local Jan`y' = "1jan" + "`y'"
				local st_JanYear`y' = date("`Jan`y''", "DMY")

			*Loop over scalars
				forvalues i = `=scalar(n_`y')'/`=scalar(N_`y')' {
					qui replace DateDN = `st_JanYear`y'' + round((`i'-`=scalar(n_`y')')/`=scalar(D_`y')') if(ID == `i' & Year == `y')
				}
			}
		
		*Clean
			drop Age_lower Age_upper AgeGroup
			format DateDN %td
			replace Age = round(Age, 0.1)
			rename Year YearDN
			
		*Core
			gen State = 1 // Diagnosis
			gen SCT_DN = .
			gen SCT_L1 = .
			gen MNT = .
			gen Age70 = Age >= 70
			gen Age75 = Age >= 75
			
		forval l = 1/9 {
			gen TXR_L`l' = .
			gen TXD_L`l' = .
			gen TFI_L`l' = .
			gen BCR_L`l' = .
		}
			gen BCR_SCT = .
			rename TFI_L9 TFI_DN
			
		local State "DN L1S L1E L2S L2E L3S L3E L4S L4E L5S L5E L6S L6E L7S L7E L8S L8E L9S L9E"
		foreach s of local State {
			gen Age_`s' = .
			gen TNE_`s' = .
			gen TSD_`s' = .
			gen MOR_`s' = .
		}	
			replace Age_DN = Age
			drop Age
			
		order ID YearDN DateDN State Male ECOGcc RISS ISS CMc CM_CKD Age70 Age75 SCT_DN SCT_L1 MNT Age* TSD* TNE* TXR* TXD* TFI_DN TFI* BCR* MOR*

end

**********
*Execute based on arguments	

forval s = 1/10 {
				
	*Execute function	
		population `s'
			
	*Sample
		gen Sample = `s'	
		
	*Save
		save "patients/population_1995_2040_`s'.dta", replace
}
