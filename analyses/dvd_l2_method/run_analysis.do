****************************************
* EpiMAP Myeloma v2.0 - Execution Script
* 
* Purpose: Execute simulation based on settings
*
* Author: EpiMAP Research Team
* Date: October 2025
*****************************************

clear all
macro drop _all

cd "/Users/adami/Documents/Monash/Research/Blood Disorders/Myeloma/EpiMAP/Simulation"

// Define settings
global analysis		"dvd_l2_method"    	// Analysis name
global int			"dvd"            	// Intervention
global line         "2"             	// Line being assessed (1/9)
global coeffs		"dvd_l2_pre"    	// Which coefficients
global data		    "predicted"   	 	// Which patients
global min_year		"2000"				// Patients diagnosed from (>= 1995)
global max_year		"2025"				// Patients diagnosed until (<= 2040)
global min_id       "1"             	// First patient ID (>= 1)
global max_id       "10"     	      	// Last patient ID (Population <= 101,212, Prediction <= 79991)
global boot		    "0"             	// Bootstrap flag (0/1)
global min_bs 		""	             	// First bootstrap
global max_bs 		"" 	            	// Last bootstrap
global report       "0"             	// Report flag (0/1)
global scenario		"A_trial"			// A_trial / B_ccbm / C_mrdr

// Execute simulation
do "EpiMAP_Myeloma.do" ///
		`analysis' `int' `line' `coeffs' `data' `min_year' `max_year' ///
		`min_id' `max_id' `boot' `min_bs' `max_bs' `report' `scenario'
