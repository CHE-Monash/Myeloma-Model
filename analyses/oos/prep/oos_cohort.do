**********
* Monash Myeloma Model - OOS (70/30): held-out simulation cohort
*
* Purpose: Turn the held-out 30% patients' BASELINE (diagnosis) covariates into a simulation input
*          cohort, in the same schema as patients/population_1995_2040_*.dta, so the engine can
*          simulate the real test patients (not a synthetic incident population). The simulated
*          trajectories are then compared to those same patients' observed outcomes (oos_targets.do).
*
* Input:   ${data_path}/oos/MRDR Wide MI_test.dta   (test fold, imputed; one row per patient x imp,
*          at diagnosis; from prep/multiple_imputation.do with $sample=="test")
* Output:  analyses/oos/patients/oos_cohort.dta
*
* The finalisation block below mirrors prep/population_1995_2040.do (cohort schema must match what
* core/load_patients.do / core/mata_setup.do read). Run from the repository root.
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

capture mkdir "analyses/oos/patients"

* ---- Which imputation to use for the input covariates ----
* One completed imputation is enough for the cohort's baseline covariates; coefficient/bootstrap
* uncertainty is carried on the model side. (Could loop Imp to add input uncertainty if desired.)
local use_imp = 1

use "${data_path}/oos/MRDR Wide MI_test.dta", clear
keep if Imp == `use_imp'

* ---- Baseline fields (mirror population_1995_2040.do finalisation) ----
gen DateDN = Date0
format DateDN %td
gen YearDN = yofd(Date0)
replace Age = round(Age, 0.1)

* Core covariates / flags
gen State  = 1            // Diagnosis
gen SCT_DN = .
gen SCT_L1 = .
gen MNT    = .
gen Age70  = Age >= 70
gen Age75  = Age >= 75

* Per-line outcome placeholders (engine fills these)
forval l = 1/9 {
	gen TXR_L`l' = .
	gen TXD_L`l' = .
	gen TFI_L`l' = .
	gen BCR_L`l' = .
}
gen BCR_SCT = .
gen TFI_DN  = .

* Per-state placeholders
local State "DN L1S L1E L2S L2E L3S L3E L4S L4E L5S L5E L6S L6E L7S L7E L8S L8E L9S L9E"
foreach s of local State {
	gen Age_`s' = .
	gen TNE_`s' = .
	gen TSD_`s' = .
	gen MOR_`s' = .
}
replace Age_DN = Age
drop Age

order ID YearDN DateDN State Male ECOGcc RISS ISS CM_CKD Age70 Age75 SCT_DN SCT_L1 MNT Age* TSD* TNE* TXR* TXD* TFI_DN TFI* BCR* MOR*

gen Sample = 1
label data "OOS held-out 30% simulation cohort (imputation `use_imp')"
save "analyses/oos/patients/oos_cohort.dta", replace

di _n "OOS cohort: " _N " held-out patients written to analyses/oos/patients/oos_cohort.dta"
