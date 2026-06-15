**********
* EpiMAP Myeloma - DVd L2 Predicted Cohort Generation
*
* Purpose: Generate predicted patient cohort for DVd L2 analysis
*          Identifies all patients reaching L2 during 2020-2025
*
* Inputs:  patients/population_1995_2040_[1-10].dta
* Outputs: analyses/td_l4_pre/patients/predicted_cohort.dta
*
* Author: Adam Irving
* Date: June 2026
**********

clear all
set more off
macro drop _all

cd "/Users/adami/Documents/Monash/Vault/research/models/myeloma model/repo"

**********
* Configuration
**********

// Parameters
local n_samples = 10
local start_year = 2020
local end_year = 2025

**********
* Generate predicted patients from each population sample
**********

forval s = 1/`n_samples' {
	
	mata: mata clear
	
	// Define settings
	global analysis		"transport_dvd"    	// Analysis name
	global int			"dvd"            	// Intervention
	global line         "2"             	// Line being assessed
	global coeffs		"dvd_pre"    		// Which coefficients
	global data		    "population_`s'"	// Which patients (independent sample s: population_1995_2040_`s')
	global min_year		"1995"				// Patients diagnosed from
	global max_year		"2020"				// Patients diagnosed until
	global min_id       "1"             	// First patient ID
	global max_id       "101212"            // Last patient ID (<=101212) 
	global boot		    "0"             	// Bootstrap flag
	global cost_year	"2020"				// Cost year
	global drate		"0.05"				// Annual discount rate (PBAC = 5%)

	// Set Paths
	global coefficients_path    "analyses/$analysis/coefficients"
	global outcomes_path		"analyses/$analysis/outcomes"
	global patients_path        "analyses/$analysis/patients"
	global simulated_path       "analyses/$analysis/simulated"
	global populations_path     "data/populations"

	// Load programs
	run "core/load_patients.do"
	run "core/mata_setup.do"
	run "core/simulation_engine.do"
	run "core/process_data.do"
	
	// Load coefficients
    qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
    
    // Load utility functions
    run "core/mata_functions.do"
    
    // Execute simulation pipeline
    load_patients
    mata_setup
    simulation
    process_data
    
    // Filters 
	keep if MOR_L`= ${line} - 1'E == 0
    keep if YearL${line} >= `start_year' & YearL${line} <= `end_year'
	
	// Clean for simulation 
    replace State = ${line} * 2
	replace Age_L${line}S = .
	replace TXR_L${line} = .
    replace ID = _n
	replace DateDN = td(1jan2020) - (TSD_L${line}S * 12) // All patients must start L4 on 1 Jan 2020
	replace YearDN = yofd(DateDN)
	drop DateL* DateMOR YearL* YearMOR c* q* OC_TIME_L TSD_*_ref 
	cap drop DateSCT YearSCT
    
    // Add sample identifier
    gen Sample = `s'
    
    // Save predicted population
    save "analyses/$analysis/patients/predicted_`s'.dta", replace
}

**********
* Combine all predicted populations
**********

clear

// Load first file
use "analyses/$analysis/patients/predicted_1.dta", clear
erase "analyses/$analysis/patients/predicted_1.dta"

// Append remaining files
forval s = 2/`n_samples' {
	append using "analyses/$analysis/patients/predicted_`s'.dta", nolabel
	erase "analyses/$analysis/patients/predicted_`s'.dta" 
}

// Overall summary
di as text "Total patients across all samples: " as result _N

// Sample distribution
tab Sample, missing

// Time Since Diagnosis (months)
summarize TSD_L${line}S, detail

// Save combined cohort
save "analyses/$analysis/patients/patients_${analysis}_${line}.dta", replace

