**********
* Monash Myeloma Model - Calibrated Transport Dispatcher
*
* Purpose: Simulation dispatcher for the DVd L2 calibrated transport analysis; configures the run
*          and executes the engine. Optional positional args: <boot> <min_bs> <max_bs> <scenario>.
**********

* Optional positional args (bootstrap run + scenario select), read into locals BEFORE clear all:
local a_boot  `1'
local a_minbs `2'
local a_maxbs `3'
local a_scen  `4'   // scenario override (A_trial / B_transport / C_mrdr); default set below

clear all
set more off
global cohort_file  ""

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $scratch_path (git-ignored; see config.example.do)

**********
* Configuration
**********

// Analysis settings
global analysis     "transport_dvd"    	// Analysis name
global int1         "dvd"               // Intervention
global int0			"vd"				// Comparator ("" for single arm)
global line         "2"                 // Line being assessed (1-9)
global coeffs       "transport_dvd"     // Coefficient set (transport_dvd)
global data         "predicted"         // Patient data (predicted / synthetic)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2025"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID (>= 1)
global max_id       "50000"             // Last patient ID (Prediction 105955)
global cost_year	"2025"				// Cost year
global drate		"0.05"				// Annual discount rate (PBAC = 5%)
global report       "0"                 // Generate report (0/1)
global scenario     "B_transport"     	// Scenario (A_trial / B_transport / C_mrdr)

local single_arm = ("$int0" == "")

// Bootstrap settings
global boot         "0"                 // Bootstrap flag (0/1)
if "`a_boot'"  != "" global boot   "`a_boot'"
if "`a_minbs'" != "" global min_bs "`a_minbs'"
if "`a_maxbs'" != "" global max_bs "`a_maxbs'"
if "`a_scen'"  != "" global scenario "`a_scen'"   // overrides the default above (run.do loops scenarios)

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path		"analyses/$analysis/outcomes"
global patients_path        "analyses/$analysis/patients"
global simulated_path       "analyses/$analysis/simulated"
global scratch_path			"${scratch_path}/$analysis/simulated"

**********
* Load Programs
**********

run "core/export_results.do"
run "core/run_pipeline.do"

**********
* Execute Simulation
**********

// No Bootstrap
if ($boot == 0) {
	    
    // Intervention arm
	global int "$int1"
	
		// Load coefficients
		qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
	
		// Execute pipeline
		run_pipeline
		export_results
		
		// Save results
		save "$simulated_path/$scenario/${int}_${line}_${data}.dta", replace
    
		if (`single_arm') {
		
			// Validate results
			run "core/validation.do"
			
			// Generate report
			if ("$report" == "1") run "core/generate_report.do"
		}

    // Comparator arm
	if (!`single_arm') {
		
		global int "$int0"
		qui mata: mata clear
		
		// Load coefficients
		qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
	
		// Execute pipeline
		run_pipeline
		export_results
		
		// Save results
		save "$simulated_path/$scenario/${int}_${line}_${data}.dta", replace
		
		// Calculate ICER
        gen arm = 0
		
			// Append intervention arm
			preserve
			use "$simulated_path/$scenario/${int1}_${line}_${data}.dta", clear
			gen arm = 1
			tempfile intervention
			qui save `intervention'
			restore
			append using `intervention'
			
			tab BCR_L2 TXR_L2, col
			
			// Calculate mean outcomes by arm
			collapse (mean) cost=cost_total_d qaly=qaly_total_d ly=OC_TIME, by(arm)
			local cost0 = cost[1]
			local cost1 = cost[2]
			local qaly0 = qaly[1]
			local qaly1 = qaly[2]
			
			di "Inc Cost: $" %7.0fc (`cost1' - `cost0')
			di "Inc QALY: " %6.2f (`qaly1' - `qaly0')
			di "ICER: $" %8.0fc ((`cost1' - `cost0') / (`qaly1' - `qaly0'))
	}
}
// Bootstrap
else if ($boot == 1) {

	// Ensure bootstrap output folder exists
	cap mkdir "$scratch_path/$scenario/bootstrap"

	// Iteration
	forvalues b = $min_bs / $max_bs {
	global b "`b'"
	mata: mata clear 
	
		// Intervention arm
		global int "$int1"
		  
			// Load coefficients
			qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"
			
			// Execute pipeline
			run_pipeline
			
			// Save results
			save "$scratch_path/$scenario/bootstrap/${int}_${line}_${data}_B`b'.dta", replace
			
		// Comparator arm
		global int "$int0"
		qui mata: mata clear
		
			// Load coefficients
			qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"
			
			// Execute pipeline
			run_pipeline
			
			// Save results
			save "$scratch_path/$scenario/bootstrap/${int}_${line}_${data}_B`b'.dta", replace
    }
}
