**********
* EpiMAP Myeloma - Base Model Dispatcher
* 
* Purpose: Execute analysis
**********

clear all
set more off

**********
* Configuration
**********

// Set working directory
cap cd "/Users/adami/Documents/Monash/Research/Blood Disorders/EpiMAP-Local/Myeloma/Simulation"

// Analysis settings
global analysis     "base_model"    	// Analysis name
global int          "all"               // Intervention
global line         "0"                 // Line being assessed (1-9)
global coeffs       "base_model"        // Coefficient set (dvd_l2_pre / dvd_l2_post)
global data         "population"        // Patient data (predicted / population)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2040"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID (>= 1)
global max_id       "10"            // Last patient ID (<= 101,212)
global boot         "0"                 // Bootstrap flag (0/1)
global min_bs       ""                  // First bootstrap iteration
global max_bs       ""                  // Last bootstrap iteration
global cost_year	"2025"				// Price year for all costs (AUD)
global drate		"0.05"				// Annual discount rate (PBAC = 5%)
global report       "0"                 // Generate report (0/1)
global scenario     ""           		// Scenario

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path		"analyses/$analysis/outcomes"
global patients_path        "analyses/$analysis/patients"
global simulated_path       "analyses/$analysis/simulated"

**********
* Load Programs
**********

run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"

**********
* Execute Simulation
**********

if ("$boot" == "0") {
	
// No Bootstrapping
    
    // Load coefficients
    qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
    
    // Load utility functions
    run "core/mata_functions.do"
    
    // Execute simulation pipeline
    load_patients
    mata_setup
    simulation
    // process_data
    
    // Save results
    save "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}_${scenario}.dta", replace
	
	// Validate results
	run "core/validation.do"
	
	// Generate report
	if ("$report" == "1") qui do "core/generate_report.do"

}
else {
	
	// Bootstrapping
	forvalues b = $min_bs/$max_bs {
	global BSIteration "`b'"
	mata: mata clear 
      
		// Load coefficients
		qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_$coeffs_B`b'"
			
		// Load utility functions
		run "core/mata_functions.do"
			
		// Execute simulation pipeline
		load_patients
		mata_setup
		simulation
        process_data
        
        // Save results
        save "$simulated_path/bootstrap/${int}_${line}_${data}_${min_id}_${max_id}_${scenario}_B`b'.dta", replace
        
        di as text "Iteration `b' completed"
    }
}
