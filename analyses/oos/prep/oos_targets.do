**********
* Monash Myeloma Model - OOS (70/30): held-out validation targets
*
* Purpose: Build the OBSERVED outcome targets for the held-out 30% (OS, BCR, TXD, TFI, pathways) --
*          what the model's predictions for those patients are compared against.
*
* This is the dedicated OOS targets script, but it deliberately REUSES prep/generate_benchmarks.do
* (the single source of truth for the estimators: KM, the M12/M24 horizon-survival columns, the
* competing-risks pathways CIF, the L1-end-reacher ASCT denominator) by passing it the test-fold
* imputed data and an OOS output directory as ARGUMENTS. Identical estimators are what make the OOS
* comparison fair; arguments (not globals) are used because generate_benchmarks.do runs clear all.
* (Invoke with do, not run -- run suppresses all screen output.)
*
* Requires (run analyses/oos/run.do step 1, the test fold, first):
*   ${data_path}/oos/MRDR Long MI_test.dta     (test fold, imputed separately)
* Output:
*   analyses/oos/targets/   -- the 13 csv files (same schema the model validates against)
*
* Run from the repository root.
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
