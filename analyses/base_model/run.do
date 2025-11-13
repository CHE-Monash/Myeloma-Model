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
global analysis		"base_model"    	// Analysis name
global int			"all"            	// Intervention
global line         "1"             	// Line being assessed
global coeffs		"base_model"    	// Which coefficients
global data		    "population"    	// Which patients
global min_id       "1"             	// First patient ID
global max_id       "10000"            	// Last patient ID  
global boot		    "0"             	// Bootstrap flag
global min_bs 		""			     	// First bootstrap
global max_bs 		""            	 	// Last bootstrap
global report       "0"             	// Report flag

// Execute simulation
do "EpiMAP_Myeloma.do" ///
    `analysis' `int' `line' `coeffs' `data' `min_id' `max_id' /// 
	`boot' `min_bs' `max_bs' `report'
