**********
* Monash Myeloma Model - VRd Line 1 Post-Market Dispatcher
*
* Purpose: VRd at line 1, post-market impact. Uses a coefficient set in which
*          VRd is excluded from the risk equations; $int toggles the comparison
*          (SoC = VRd-eligible patients receive historical alternatives;
*           VRd = VRd available).
**********

clear all
set more off

if "$repo_root" != "" cd "$repo_root"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

**********
* Configuration
**********

// Analysis settings
global analysis     "vrd_post"          // Analysis name
global int          "VRd"               // Intervention scenario (SoC / VRd)
global line         "1"                 // Line being assessed (1-9)
global coeffs       "vrd_post"          // Coefficient set (VRd excluded)
global data         "predicted"         // Patient data (predicted / population)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2040"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID (>= 1)
global max_id       "999999"            // Last patient ID (high cap = whole cohort)
global boot         "0"                 // Bootstrap flag (0/1)
global min_bs       ""                  // First bootstrap iteration
global max_bs       ""                  // Last bootstrap iteration
global cost_year    "2025"              // Price year for all costs (AUD)
global drate        "0.05"              // Annual discount rate (PBAC = 5%)
global report       "0"                 // Generate report (0/1)
global scenario     ""                  // Scenario

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path        "analyses/$analysis/outcomes"
global patients_path        "analyses/$analysis/patients"
global simulated_path       "analyses/$analysis/simulated"

// Predicted cohort uses a legacy filename; point load_patients at it directly.
// (Rename to patients_vrd_post_1.dta to drop this override.)
global cohort_file          "$patients_path/patients_vrd_l1_post.dta"

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
    save "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}.dta", replace

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
        save "$simulated_path/bootstrap/${int}_${line}_${data}_${min_id}_${max_id}_B`b'.dta", replace

        di as text "Iteration `b' completed"
    }
}
