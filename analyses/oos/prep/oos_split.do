**********
* Monash Myeloma Model - OOS (70/30): patient split
*
* Purpose: Assign each patient (ID) to TRAIN (70%) or TEST (30%), once, with a fixed seed, on the
*          PRE-IMPUTATION data -- so imputation and risk-equation fitting never see the held-out
*          30%. Everything downstream (multiple_imputation.do, risk_equations.do via $sample, and
*          the OOS prep/validation scripts) keys off this crosswalk.
*
* Output:  ${data_path}/oos/oos_split.dta   (ID, fold)
*
* Run once, from the repository root. Restricted MRDR data ($data_path via config.do).
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

* ---- Tunables ----
local train_frac = 0.70
local seed       = 20260630      // fixed for reproducibility -- do not change once results are out

capture mkdir "${data_path}/oos"

use "${data_path}/MRDR Long.dta", clear

* One row per patient (the split is at patient level so all of a patient's records share a fold)
bysort ID: keep if _n == 1
keep ID

* NOTE: a simple random split is used here. To match / improve on the 2024 PLOS ONE analysis,
* consider STRATIFYING by diagnosis era (and/or key prognostic covariates) so train and test are
* balanced -- e.g. sort by the stratifier, then assign within strata. Confirm against the 2024 methods.

set seed `seed'
gen double _u = runiform()
gen str5 fold = cond(_u <= `train_frac', "train", "test")
drop _u

label data "OOS 70/30 split crosswalk (seed `seed')"
di _n "OOS split:"
tab fold

save "${data_path}/oos/oos_split.dta", replace
