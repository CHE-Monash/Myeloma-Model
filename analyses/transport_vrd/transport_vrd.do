**********
* Monash Myeloma Model
* Tansport VRd - Analysis Dispatcher
* 
* Purpose: Execute analysis
**********

clear all
set more off

cap cd "/Users/adami/Documents/Monash/Vault/research/models/myeloma model/repo"
*cap cd "~/em76/adam"

**********
* Configuration
**********

global analysis     "transport_vrd"    	// Analysis name
global int1         "vrd"               // Intervention
global int0			"rd"				// Comparator ("" for single arm)
global line         "1"                 // Line being assessed (1-9)
global coeffs       "vrd_post"       	// Coefficient set (dvd_pre / dvd_post)
global data         "predicted"         // Patient data (predicted / population)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2040"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID (>= 1)
global max_id       "105955"            // Last patient ID (Prediction 105955)
global boot         "1"                 // Bootstrap flag (0/1)
global min_bs       `1'                 // First bootstrap iteration
global max_bs       `2'                 // Last bootstrap iteration
global cost_year	"2020"				// Cost year
global drate		"0.05"				// Annual discount rate (PBAC = 5%)
global report       "1"                 // Generate report (0/1)
global scenario     "A_trial"     		// Scenario (A_trial / B_transport / C_mrdr)

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
run "core/export_results.do"
run "core/run_pipeline.do"

cap program drop simulation_pipeline
program define simulation_pipeline
		
	// Shared engine pass: entry -> processed per-patient outcomes
	//   (run_pipeline loads mata_functions + rng_slots, then runs the pass + builds mRN)
	run_pipeline

	// Export CSV results (runs by default; skips bootstrap internally)
	export_results
	
end

**********
* Execute Simulation
**********

// No Bootstrap
if ($boot == 0) {
	    
    // Intervention arm
	global int "$int1"
	
		// Load coefficients
		qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
	
		// Simulate
		simulation_pipeline
		
		// Save results
		save "$simulated_path/$scenario/${int}_${line}_${data}_${min_id}_${max_id}.dta", replace
    
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
		save "$simulated_path/$scenario/${int}_${line}_${data}_${min_id}_${max_id}.dta", replace
		
		// Calculate ICER
        gen arm = 0
		
			// Append intervention arm
			preserve
			use "$simulated_path/$scenario/${int1}_${line}_${data}_${min_id}_${max_id}.dta", clear
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
	cap mkdir "$simulated_path/$scenario/bootstrap"

	// Iteration
	forvalues b = $min_bs / $max_bs {
	global b "`b'"
	mata: mata clear 
	
		// Intervention arm
		global int "$int1"
		  
			// Load coefficients
			qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"
			
			// Simulate
			simulation_pipeline
			
			// Save results
			save "$simulated_path/$scenario/bootstrap/${int}_${line}_${data}_B`b'.dta", replace
			
		// Comparator arm
		global int "$int0"
		qui mata: mata clear
		
			// Load coefficients
			qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"
			
			// Simulate
			simulation_pipeline
			
			// Save results
			save "$simulated_path/$scenario/bootstrap/${int}_${line}_${data}_B`b'.dta", replace
    }
}
