**********
* EpiMAP Myeloma - td_l4_pre
* 
* Purpose: Execute analysis
**********

clear all
set more off

cap cd "/Users/adami/Documents/Monash/Research/Blood Disorders/EpiMAP-Local/Myeloma/Simulation"

**********
* Configuration
**********

global analysis     "td_pre"			// Analysis name
global int1         "td"				// Intervention (td / soc)
global int0			"soc"				// Comparator ("" for single arm)
global line         "2"					// Line being assessed (1-9)
global coeffs       "td_pre"			// Coefficient set (dvd_l2_pre / dvd_l2_post)
global data         "predicted"			// Patient data (predicted / population)
global min_year     "1995"				// Patients diagnosed from (>= 1995)
global max_year     "2040"				// Patients diagnosed until (<= 2040)
global min_id       "1"					// First patient ID (>= 1)
global max_id       "999999"			// Last patient ID (Pop <= 101,212, Pred <= 46,628)
global boot         "0"					// Bootstrap flag (0/1)
global min_bs       `1'					// First bootstrap iteration
global max_bs       `2'					// Last bootstrap iteration
global cost_year	"2025"				// Price year for all costs (AUD)
global drate		"0.05"				// Annual discount rate (PBAC = 5%)
global report       ""					// Generate report (0/1)
global scenario     ""					// Scenario

local single_arm = ("$int0" == "")

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
run "core/run_pipeline.do"

cap program drop simulation_pipeline
program define simulation_pipeline
		
	// Shared engine pass: entry -> processed per-patient outcomes
	//   (run_pipeline loads mata_functions + rng_slots, then runs the pass + builds mRN)
	run_pipeline
	
end

**********
* Execute Simulation
**********

// No Bootstrap
if ("$boot" == "0") {
	    
    // Intervention arm
	global int "$int1"
	
		// Load coefficients
		qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
	
		// Simulate
		simulation_pipeline
		
		// Save results
		save "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}.dta", replace
    
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
	
		// Simulate
		simulation_pipeline
		
		// Save results
		save "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}.dta", replace
		
		// Calculate ICER
        gen arm = 0
		
			// Append intervention arm
			append using "$simulated_path/${int1}_${line}_${data}_${min_id}_${max_id}.dta"
			replace arm = 1 if arm == .
			tab BCR_L$line arm, col
			
			// Calculate mean outcomes by arm
			local outcomes cost_total_d qaly_total_d
			foreach o of local outcomes {
				forval i = 0/1 {
					qui sum `o' if arm == `i'
					local `o'_`i' = r(mean)
				}
			}

			di "Inc Cost: $" %7.0fc (`cost_total_d_1' - `cost_total_d_0')
			di "Inc QALY: " %6.2f (`qaly_total_d_1' - `qaly_total_d_0')
			di "ICER: $" %8.0fc ((`cost_total_d_1' - `cost_total_d_0') / (`qaly_total_d_1' - `qaly_total_d_0'))
	}
}
// Bootstrap
qui if ("$boot" == "1") {
	
	// Iteration
	forvalues b = $min_bs / $max_bs {
	global b "`b'"

		// Loop over intervantions
		foreach int in $int1 $int0 {
			
			mata: mata clear
			global int `int'
			
			// Load coefficients
			qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"
			
			// Simulate
			simulation_pipeline
			
			// Save results
			save "$simulated_path/bootstrap/${int}_${line}_${data}_B`b'.dta", replace
		}
    }
}
