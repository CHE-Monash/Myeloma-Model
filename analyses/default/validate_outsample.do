**********
* Monash Myeloma Model - Validate out-of-sample
*
* Purpose: Compare the simulated held-out 30% (analyses/default/simulated/outsample/) against those
*          patients' OBSERVED outcomes (analyses/default/targets/) -- OS & TXD/TFI by horizon survival,
*          BCR by category, pathways by competing-risks CIF, ASCT among L1-end reachers. Reuses the
*          shared engine validate_outcomes.do via its $val_targets / $val_simfile globals.
* Notes:   Run after `simulate.do 0 . . outsample` (the point-estimate out-of-sample run, $boot 0).
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

* Log so the run can be reviewed after the fact (scratch/ is git-ignored). Mirrors validate_insample.do.
capture log close _all
log using "scratch/validate_outsample.log", replace text

* Point the shared validator at the out-of-sample targets + the out-of-sample simulated dataset
global val_targets "analyses/default/targets"
global val_simfile "analyses/default/simulated/outsample/all_0_test.dta"

do "analyses/default/validate_outcomes.do"

* Reset the globals after the run
global val_targets ""
global val_simfile ""

capture log close _all

* ----------------------------------------------------------------------------------------------
* PREDICTION-INTERVAL CALIBRATION (the headline out-of-sample metric) is in bootstrap_validation.do.
*   The block above validates the POINT estimate (one 70%-trained coefficient set). For calibration:
*     do "analyses/default/simulate.do" 1 1 500 outsample   // 500 bootstrap sims of the held-out 30%
*     do "analyses/default/bootstrap_validation.do"          // 95% percentile PI per outcome; coverage
*   It forms [p2.5, p97.5] across the 500 resamples per target and checks whether each held-out
*   OBSERVED value falls inside (percentile method, as in the 2024 paper). Designed to run on the HPC;
*   see run.do bootstrap section.
* ----------------------------------------------------------------------------------------------
