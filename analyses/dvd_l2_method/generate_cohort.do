**********
* EpiMAP Myeloma - DVd L2 Predicted Cohort Generation
*
* Purpose: Generate predicted patient cohort for DVd L2 pre-market analysis
*          Identifies all patients reaching L2 during 2021-2025
*
* Inputs:  patients/population_1995_2040_[1-10].dta
* Outputs: analyses/dvd_l2_method/data/patients/predicted_[1-10].dta
*          analyses/dvd_l2_method/data/patients/predicted_cohort.dta
*
* Author: Adam Irving
* Date: November 2025
**********

clear all
set more off
macro drop _all

**********
* Configuration
**********

// Set working directory
cd "/Users/adami/Documents/Monash/Research/Blood Disorders/Myeloma/EpiMAP/Simulation"

// Parameters
local n_samples = 10
local diagnosis_year = 2020
local l2_start_year = 2021
local l2_end_year = 2025

**********
* STEP 1: Generate predicted patients from each population sample
**********

forval s = 1/`n_samples' {
	
	// Define settings
	global analysis		"dvd_l2_method"    	// Analysis name
	global int			"dvd"            	// Intervention
	global line         "2"             	// Line being assessed
	global coeffs		"dvd_l2_pre"    	// Which coefficients
	global data		    "population"    	// Which patients
	global min_year		"1995"				// Patients diagnosed from
	global max_year		"2025"				// Patients diagnosed until
	global min_id       "1"             	// First patient ID
	global max_id       "101212"            // Last patient ID (<=101212) 
	global boot		    "0"             	// Bootstrap flag
	global min_bs 		"1"             	// First bootstrap
	global max_bs 		"10"             	// Last bootstrap
	global report       "0"             	// Report flag
	global scenario		""					// A_trial / B_ccbm / C_mrdr

	// Execute simulation
	do "EpiMAP_Myeloma.do" ///
		`analysis' `int' `line' `coeffs' `data' `min_year' `max_year' ///
		`min_id' `max_id' `boot' `min_bs' `max_bs' `report'
    
    // Filters 
	keep if MOR_L1E == 0 											// Alive at L2S
    keep if YearL2 >= `l2_start_year' & YearL2 <= `l2_end_year'	// L2S between dates
      
    // Clean for simulation 
    replace State = 4
	replace Age_L2S = .
	replace TXR_L2 = .
	gen Age70 = Age_DN >= 70
	gen Age75 = Age_DN >= 75
    
    // Reset ID
    replace ID = _n
    
    // Add sample identifier
    gen Sample = `s'
    
    // Save predicted population
    save "analyses/dvd_l2_method/data/patients/predicted_`s'.dta", replace
}

**********
* STEP 2: Combine all predicted populations
**********

clear
local files_found = 0

// Load first available file
use "analyses/dvd_l2_method/data/patients/predicted_1.dta", clear
erase "analyses/dvd_l2_method/data/patients/predicted_1.dta"

// Append remaining files
forval s = 2/`n_samples' {
	append using "analyses/dvd_l2_method/data/patients/predicted_`s'.dta", nolabel
	erase "analyses/dvd_l2_method/data/patients/predicted_`s'.dta" 
}

**********
* STEP 3: Summary and validation
**********

// Overall summary
di as text "Total patients across all samples: " as result _N
di as text "Population samples found: " as result `files_found' as text " of `n_samples'"

// Sample distribution
di as text _newline "Distribution by sample:"
tab Sample, missing

// Year distribution
di as text _newline "Distribution by L2 start year:"
tab YearL2, missing

// Treatment regimen distribution at L2
di as text _newline "Treatment regimen distribution at L2:"
tab TXR_L2, missing

// Age distribution
di as text _newline "Age distribution at L2:"
summarize Age_L2S, detail

// Months to L2 distribution
di as text _newline "Time from diagnosis to L2:"
summarize TSD_L2S, detail

// Previous line response distribution
di as text _newline "Line 1 best clinical response distribution:"
tab BCR_L1, missing

// Save combined cohort
local outfile "analyses/dvd_l2_method/data/patients/patients_$analysis.dta"
save "`outfile'", replace

