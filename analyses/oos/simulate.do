**********
* Monash Myeloma Model - OOS (70/30): simulate.do (simulation dispatcher)
*
* Simulates the held-out 30% patients using coefficients trained on the 70%, then compare to their
* observed outcomes with validate_oos.do. Orchestrated by run.do; on the HPC it is sbatch'd directly
* (it never sources run.do). Mirrors analyses/base_model/simulate.do, but:
*   - the input cohort is the real held-out patients (analyses/oos/patients/oos_cohort.dta) via
*     $cohort_file, not the synthetic population; and
*   - coefficients come from analyses/oos/coefficients/ (trained on the 70%).
*
* Point estimate: $boot 0. Prediction intervals: $boot 1 with $min_bs/$max_bs over the 70%
* bootstrap coefficient sets (one simulated dataset per resample; validate_oos.do aggregates).
**********

* optional positional args for the bootstrap simulation, read into locals BEFORE clear all:
*   do simulate.do            -> point estimate ($boot 0)
*   do simulate.do 1 1 500    -> bootstrap sims 1-500 (HPC: pass an array chunk, e.g. 1 101 200)
local a_boot  `"`1'"'
local a_minbs `"`2'"'
local a_maxbs `"`3'"'

clear all
set more off

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

**********
* Configuration
**********
global analysis     "oos"
global int          "all"
global line         "0"
global coeffs       "oos"                                   // -> coefficients_oos[.mmat / _B*]
global data         "oos"                                   // not "population_N" -> use $cohort_file
global cohort_file  "analyses/oos/patients/oos_cohort.dta"  // real held-out 30% patients
global min_year     "1995"
global max_year     "2040"
global min_id       "1"
global max_id       "101212"
global boot         "0"                                     // 0 = point estimate; 1 = bootstrap PIs
global min_bs       ""
global max_bs       ""
if `"`a_boot'"'  != "" global boot   `"`a_boot'"'           // args override (for the bootstrap run)
if `"`a_minbs'"' != "" global min_bs `"`a_minbs'"'
if `"`a_maxbs'"' != "" global max_bs `"`a_maxbs'"'
global cost_year    "2025"
global drate        "0.05"
global report       "0"
global scenario     ""

**********
* Set Paths
**********
global coefficients_path    "analyses/$analysis/coefficients"
global outcomes_path        "analyses/$analysis/outcomes"
global patients_path        "analyses/$analysis/patients"
global simulated_path       "analyses/$analysis/simulated"
capture mkdir "$simulated_path"
capture mkdir "$simulated_path/bootstrap"

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

    // Point estimate (70%-trained coefficients)
    qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"
    run "core/mata_functions.do"

    load_patients
    mata_setup
    simulation
    process_data

    save "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}.dta", replace
    run "core/validation.do"
}
else {

    // Bootstrap prediction intervals: one simulated dataset per 70% resample
    forvalues b = $min_bs/$max_bs {
        global BSIteration "`b'"
        mata: mata clear

        qui mata: mata matuse "$coefficients_path/bootstrap/coefficients_${coeffs}_B`b'"
        run "core/mata_functions.do"

        load_patients
        mata_setup
        simulation
        process_data

        save "$simulated_path/bootstrap/${int}_${line}_${data}_${min_id}_${max_id}_B`b'.dta", replace
    }
}
