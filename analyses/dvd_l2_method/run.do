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
global line         "2"             	// Line being assessed
global coeffs		"dvd_l2_method"    	// Which coefficients
global data		    "population"    	// Which patients
global min_id       "1"             	// First patient ID
global max_id       "10000"            	// Last patient ID  
global boot		    "1"             	// Bootstrap flag
global min_bs 		"1"             	// First bootstrap
global max_bs 		"10"             	// Last bootstrap
global report       "0"             	// Report flag

// Execute simulation
do "EpiMAP_Myeloma.do" ///
    `analysis' `int' `line' `coeffs' `data' `min_id' `max_id' /// 
	`boot' `min_bs' `max_bs' `report'
