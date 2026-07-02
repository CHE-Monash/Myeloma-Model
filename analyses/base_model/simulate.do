**********
* Monash Myeloma Model - base_model: simulate.do (simulation dispatcher)
*
* Current-practice projection over the synthetic population (all regimens). Orchestrated by run.do;
* on the HPC it is sbatch'd directly (it never sources run.do).
*
* Point estimate: $boot 0. Bootstrap: $boot 1 with $min_bs/$max_bs over the coefficient resamples.
**********

* optional positional args for the bootstrap run, read into locals BEFORE clear all:
*   do simulate.do            -> point estimate ($boot 0)
*   do simulate.do 1 1 500    -> bootstrap iterations 1-500 (HPC: pass an array chunk, e.g. 1 101 200)
local a_boot  `"`1'"'
local a_minbs `"`2'"'
local a_maxbs `"`3'"'

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

**********
* Configuration
**********


// Analysis settings
global analysis     "base_model"    	// Analysis name
global int          "all"               // Intervention
global line         "0"                 // Line being assessed (1-9)
global coeffs       "base_model"        // Coefficient set (dvd_l2_pre / dvd_l2_post)
global data         "population"        // Patient data (predicted / population)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2040"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID (>= 1)
global max_id       "101212"            // Last patient ID (<= 101,212)
global boot         "0"                 // Bootstrap flag (0/1)
global min_bs       ""                  // First bootstrap iteration
global max_bs       ""                  // Last bootstrap iteration
if `"`a_boot'"'  != "" global boot   `"`a_boot'"'           // positional args override (for the bootstrap run)
if `"`a_minbs'"' != "" global min_bs `"`a_minbs'"'
if `"`a_maxbs'"' != "" global max_bs `"`a_maxbs'"'
global cost_year	"2025"				// Price year for all costs (AUD)
global drate		"0.05"				// Annual discount rate (PBAC = 5%)
global report       "1"                 // Generate report (0/1)
global scenario     ""           		// Scenario

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path		"analyses/$analysis/outcomes"
global patients_path        "analyses/$analysis/patients"
global simulated_path       "analyses/$analysis/simulated"

// Output partition for the simulated .dta: scenario is an optional subfolder
// (empty $scenario => simulated_path itself), matching export_results.do.
// Filename convention, reused across analyses:
//   <int>_<line>_<data>.dta                  point estimate
//   bootstrap/<int>_<line>_<data>_B<b>.dta   bootstrap iteration b
// The patient-ID range is a run setting (recorded in the report), not in the name.
global sim_out = cond("$scenario" == "", "$simulated_path", "$simulated_path/$scenario")

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
    process_data
    
    // Save results
    capture mkdir "$sim_out"
    save "$sim_out/${int}_${line}_${data}.dta", replace
	
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
		qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"
			
		// Load utility functions
		run "core/mata_functions.do"
			
		// Execute simulation pipeline
		load_patients
		mata_setup
		simulation
        process_data
        
        // Save results
        capture mkdir "$sim_out"
        capture mkdir "$sim_out/bootstrap"
        save "$sim_out/bootstrap/${int}_${line}_${data}_B`b'.dta", replace
        
        di as text "Iteration `b' completed"
    }
}
