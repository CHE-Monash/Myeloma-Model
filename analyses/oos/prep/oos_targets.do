**********
* Monash Myeloma Model - OOS targets
*
* Purpose: Build the OBSERVED outcome targets for the held-out 30% (OS, BCR, TXD, TFI, pathways) --
*          what the model's predictions for those patients are compared against. Reuses
*          prep/generate_benchmarks.do (the estimator source of truth) with the test-fold data + an
*          OOS output dir as ARGUMENTS; identical estimators keep the comparison fair.
* Notes:   Invoke with do, not run (run suppresses output; args used because generate_benchmarks runs
*          clear all). Requires ${data_path}/oos/MRDR Long MI_test.dta (run.do step 1 first) ->
*          analyses/oos/targets/ (13 csv files).
**********

if "$repo_path" != "" cd "$repo_path"   // cd using the session's repo_path BEFORE clear all wipes it
clear all
set more off
capture run "config.do"

di as text "Working directory : " c(pwd)
di as text "data_path          : [$data_path]"

capture mkdir "analyses/oos/targets"

local testdata "${data_path}/oos/MRDR Long MI_test.dta"
capture confirm file "`testdata'"
if _rc {
	di as error "Test-fold imputed data not found:"
	di as error "  `testdata'"
	di as error "Run analyses/oos/run.do step 1 (multiple_imputation ... test) first."
	exit 601
}

di as text "Building OOS targets from `testdata' ..."
do "prep/generate_benchmarks.do" "`testdata'" "analyses/oos/targets"
di as text "OOS targets written to analyses/oos/targets/"
