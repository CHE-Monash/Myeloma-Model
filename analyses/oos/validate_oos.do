**********
* Monash Myeloma Model - Validate OOS
*
* Purpose: Compare the simulated held-out 30% (analyses/oos/simulated/) against those patients'
*          OBSERVED outcomes (analyses/oos/targets/) -- OS & TXD/TFI by horizon survival, BCR by
*          category, pathways by competing-risks CIF, ASCT among L1-end reachers. Reuses the shared
*          engine validate_outcomes.do via its $val_targets / $val_simfile globals.
* Notes:   Run after analyses/oos/simulate.do (point estimate, $boot 0).
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

* Point the shared validator at the OOS targets + the OOS simulated dataset
global val_targets "analyses/oos/targets"
global val_simfile "analyses/oos/simulated/all_0_oos.dta"

do "analyses/oos/validate_outcomes.do"

* Reset so a later in-sample validation run in the same session uses its defaults
global val_targets ""
global val_simfile ""

* ----------------------------------------------------------------------------------------------
* PREDICTION-INTERVAL CALIBRATION (the headline OOS metric) is implemented in bootstrap_validation.do.
*   The block above validates the POINT estimate (one 70%-trained coefficient set). For calibration:
*     do "analyses/oos/simulate.do" 1 1 500     // 500 bootstrap sims of the held-out 30%
*     do "analyses/oos/bootstrap_validation.do" // 95% percentile PI per outcome; coverage vs observed
*   It forms [p2.5, p97.5] across the 500 resamples per target and checks whether each held-out
*   OBSERVED value falls inside (percentile method, as in the 2024 paper). Designed to run on the HPC;
*   see run.do bootstrap section (a)-(c).
* ----------------------------------------------------------------------------------------------
