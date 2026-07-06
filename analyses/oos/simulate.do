**********
* Monash Myeloma Model - OOS (70/30): simulate.do (simulation dispatcher)
*
* Simulates the held-out 30% patients using coefficients trained on the 70%, then compared to their
* observed outcomes with validate_oos.do. Orchestrated by run.do; on the HPC it is sbatch'd directly
* (it never sources run.do). Mirrors analyses/base_model/simulate.do, but:
*   - the input cohort is the real held-out patients (analyses/oos/patients/oos_cohort.dta) via
*     $cohort_file, not the synthetic population; and
*   - coefficients come from analyses/oos/coefficients/ (trained on the 70%).
*
* Point estimate: $boot 0. Prediction intervals: $boot 1 with $min_bs/$max_bs over the 70% bootstrap
* coefficient sets (one simulated dataset per resample; validate_oos/bootstrap_validation aggregate).
**********

* Optional positional args for the bootstrap run, read into locals BEFORE clear all:
*   do simulate.do            -> point estimate ($boot 0)
*   do simulate.do 1 1 500    -> bootstrap sims 1-500 (HPC: pass an array chunk, e.g. 1 101 200)
local a_boot  `"`1'"'
local a_minbs `"`2'"'
local a_maxbs `"`3'"'

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

**********
* Configuration
**********

// Analysis settings
global analysis     "oos"                                   // Analysis name
global int          "all"                                   // Intervention
global line         "0"                                     // Line being assessed (0-9)
global coeffs       "oos"                                   // Coefficient set (-> coefficients_oos)
global data         "oos"                                   // not "population_N" -> use $cohort_file
global cohort_file  "analyses/oos/patients/oos_cohort.dta"  // real held-out 30% patients
global min_year     "1995"                                  // Patients diagnosed from (>= 1995)
global max_year     "2040"                                  // Patients diagnosed until (<= 2040)
global min_id       "1"                                     // First patient ID (>= 1)
global max_id       "101212"                                // Last patient ID (<= 101,212)
global cost_year    "2025"                                  // Cost year (AUD)
global drate        "0.05"                                  // Annual discount rate (PBAC = 5%)
global report       "0"                                     // Generate report (0/1)
global scenario     ""                                      // Scenario

// Bootstrap settings
global boot         "0"                                     // Bootstrap flag (0/1)
if `"`a_boot'"'  != "" global boot   `"`a_boot'"'
if `"`a_minbs'"' != "" global min_bs `"`a_minbs'"'
if `"`a_maxbs'"' != "" global max_bs `"`a_maxbs'"'

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path        "analyses/$analysis/outcomes"
global patients_path        "analyses/$analysis/patients"
global simulated_path       "analyses/$analysis/simulated"

// Output partition for the simulated .dta: scenario is an optional subfolder
global sim_out = cond("$scenario" == "", "$simulated_path", "$simulated_path/$scenario")

**********
* Load Programs
**********

run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"
run "core/export_results.do"
run "core/run_pipeline.do"

**********
* Execute Simulation
**********

if ("$boot" == "0") {

// No Bootstrapping

    // Load coefficients (70%-trained)
    qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"

    // Execute pipeline
    run_pipeline

    // Save results
    capture mkdir "$sim_out"
    save "$sim_out/${int}_${line}_${data}.dta", replace

    // Validate results
    run "core/validation.do"

    // Generate report
    if ("$report" == "1") qui do "core/generate_report.do"

}
else {

    // Bootstrapping: one simulated dataset per 70% resample
    forvalues b = $min_bs/$max_bs {
    global BSIteration "`b'"
    mata: mata clear

        // Load coefficients
        qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"

        // Execute pipeline
        run_pipeline

        // Save results
        capture mkdir "$sim_out"
        capture mkdir "$sim_out/bootstrap"
        save "$sim_out/bootstrap/${int}_${line}_${data}_B`b'.dta", replace

        di as text "Iteration `b' completed"
    }
}
