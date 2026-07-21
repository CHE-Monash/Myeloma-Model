**********
* Monash Myeloma Model - Validate in-sample (face validity)
*
* Purpose: Compare the full-fit PROJECTION simulation (coeffs=full on the synthetic incidence
*          population, analyses/default/simulated/all_0_synthetic.dta) against the in-sample registry
*          benchmarks (scratch/benchmarks/, built by `prep/generate_benchmarks.do` with no args).
*          Reuses the shared engine validate_outcomes.do via its $val_targets / $val_simfile globals.
* Notes:   Run after the projection point estimate `simulate.do 0` ("" scenario) has produced
*          simulated/all_0_synthetic.dta, and after generate_benchmarks.do has (re)built
*          scratch/benchmarks.
*
*          THIS IS A FACE-VALIDITY CHECK, NOT THE MAINSTAY VALIDATION. It scores the SYNTHETIC
*          projection cohort against the WHOLE-registry benchmarks -- so a miss can reflect the
*          synthetic incidence cohort differing from the registry (era mix, age/stage structure) as
*          much as a model defect, and it is NOT held-out. The mainstay calibration test is the
*          out-of-sample 70/30 track (validate_outsample.do + bootstrap_validation.do); see
*          docs/validation.md. Use this only for a quick in-sample sanity read against the registry.
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"

* Log so the run can be reviewed after the fact (scratch/ is git-ignored).
capture log close _all
log using "scratch/validate_insample.log", replace text

* Point the shared validator at the in-sample benchmarks + the full-fit synthetic projection dataset
global val_targets "scratch/benchmarks"
global val_simfile "analyses/default/simulated/all_0_synthetic.dta"

do "analyses/default/validate_outcomes.do"

* Reset the globals after the run
global val_targets ""
global val_simfile ""

capture log close _all
