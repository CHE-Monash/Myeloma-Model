**********
* Monash Myeloma Model - Simulate (template dispatcher)
*
* Purpose: the standard simulation dispatcher -- a configuration block of globals, then the shared engine
*          pass (run_pipeline). COPY THIS FOLDER to analyses/<your_analysis>/, adapt, and fill in the
*          <...> placeholders in the Configuration block.
* Usage:   orchestrated by run.do; on the HPC it is sbatch'd directly (never sources run.do). Point
*          estimate: $boot 0. Bootstrap: $boot 1 with $min_bs/$max_bs over coefficient resamples.
*          Optional positional args (run.do / HPC arrays): boot min_bs max_bs [scenario].
* Notes:   Worked examples -- analyses/base_model/ (simplest single-run dispatcher; start here);
*          analyses/transport_dvd/ (scenarios A/B/C + outcome overrides, outcomes/sim_bcr_override.do).
**********

* Optional positional args, read into locals BEFORE clear all:
local a_boot  `1'
local a_minbs `2'
local a_maxbs `3'
local a_scen  `4'   // scenario override; delete if your analysis has no scenarios

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

**********
* Configuration  -- EDIT THESE
**********

// Analysis settings
global analysis     "template"          // <- your analysis name (this folder under analyses/)
global int          "all"               // Intervention label (drives regimen logic; e.g. all / VRd / dvd)
global line         "0"                 // Line assessed (0 = all lines; 1-9 = a single line, e.g. for an override)
global coeffs       "template"          // Coefficient set -> coefficients_<coeffs>.mmat + outcomes/txr_<coeffs>.do
global data         "population"        // Patient data (population / predicted)
global min_year     "1995"              // Patients diagnosed from (>= 1995)
global max_year     "2040"              // Patients diagnosed until (<= 2040)
global min_id       "1"                 // First patient ID
global max_id       "999999"            // Last patient ID
global cost_year    "2025"              // Cost year (AUD)
global drate        "0.05"              // Annual discount rate (PBAC = 5%)
global report       "0"                 // Generate PDF report (0/1)
global scenario     ""                  // Scenario label (optional; partitions outputs + selects overrides)

// Bootstrap settings
global boot         "0"                 // Bootstrap flag (0/1)
if "`a_boot'"  != "" global boot     "`a_boot'"
if "`a_minbs'" != "" global min_bs   "`a_minbs'"
if "`a_maxbs'" != "" global max_bs   "`a_maxbs'"
if "`a_scen'"  != "" global scenario "`a_scen'"   // run.do / HPC arrays loop scenarios via this arg

**********
* Set Paths
**********

global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path		"analyses/$analysis/outcomes"   // core/outcomes/sim_bcr.do & sim_txr.do
global patients_path        "analyses/$analysis/patients"   //   auto-run sim_bcr_override.do /
global simulated_path       "analyses/$analysis/simulated"  //   sim_txr_override.do from here (if present)

// Output partition for the simulated .dta: scenario is an optional subfolder
global sim_out = cond("$scenario" == "", "$simulated_path", "$simulated_path/$scenario")

// If your predicted cohort has a non-standard filename, point load_patients at it directly:
// global cohort_file       "$patients_path/<your_cohort>.dta"

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

    // Load coefficients
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

    // Bootstrapping
    forvalues b = $min_bs/$max_bs {
    global b "`b'"
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
