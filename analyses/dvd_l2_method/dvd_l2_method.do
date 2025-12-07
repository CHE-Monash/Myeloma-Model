**********
* EpiMAP Myeloma - DVd L2 Method Dispatcher
* 
* Purpose: Execute analysis
**********

clear all
macro drop _all
set more off

**********
* Configuration
**********

// Set working directory
cd "/Users/adami/Documents/Monash/Research/Blood Disorders/EpiMAP-Local/Myeloma/Simulation"

// Analysis settings
global analysis     "base_model"     // Analysis name
global int          "all"               // Intervention
global line         "0"                 // Line being assessed (1-9)
global coeffs       "base_model"       // Coefficient set (dvd_l2_pre / dvd_l2_post)
global data         "population"         // Patient data (predicted / population)
global min_year     "2000"              // Patients diagnosed from (>= 1995)
global max_year     "2025"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID (>= 1)
global max_id       "100000"            // Last patient ID
global boot         "0"                 // Bootstrap flag (0/1)
global min_bs       ""                  // First bootstrap iteration
global max_bs       ""                  // Last bootstrap iteration
global cost_year	"2025"				// Price year for all costs (AUD)
global drate		"0.05"				// Annual discount rate (PBAC = 5%)
global report       "1"                 // Generate report (0/1)
global scenario     ""            // Scenario (1_trial / 2_ccbm / 3_mrdr)

**********
* Set Paths
**********

global analysis_path        "analyses/$analysis"
global coefficients_path    "$analysis_path/data/coefficients"
global patients_path        "$analysis_path/data/patients"
global simulated_path       "$analysis_path/data/simulated"
global populations_path     "data/populations"

// Create output directories if needed
capture mkdir "$simulated_path"
capture mkdir "$simulated_path/bootstrap"
capture mkdir "$simulated_path/report"

**********
* Load Programs
**********

qui do "core/load_patients.do"
qui do "core/mata_setup.do"
qui do "core/simulation_engine.do"
qui do "core/process_data.do"

**********
* Execute Simulation
**********

if ("$boot" == "0") {
	
// No Bootstrapping
    
    // Load coefficients
    qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
    
    // Load utility functions
    qui do "core/mata_functions.do"
    
    // Execute simulation pipeline
    load_patients
    mata_setup
    simulation
    process_data
    
    // Save results
    save "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}_${scenario}.dta", replace
	
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
		qui do "core/mata_functions.do"
			
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
