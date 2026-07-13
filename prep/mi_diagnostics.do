**********
* Monash Myeloma Model - MI diagnostics
*
* Purpose: Missing-data diagnostics for the imputed BASELINE COVARIATES that enter the risk
*          equations -- demographics (Male), performance status (ECOGcc), R-ISS (and its inputs),
*          and comorbidities (CM_CKD/CRD/PLM/DBT) -- plus the first-line response (BCR_L1), read off
*          the saved MRDR Long MI.dta at each variable's own record. Reports per variable: the
*          pre-imputation percent missing, then the relative increase in variance (RVI), the
*          fraction of missing information (FMI), the relative efficiency of the current number of
*          imputations (RE) and the imputation-adjusted degrees of freedom (DF) -- the figures that
*          justify the number of imputations and go in a manuscript's missing-data methods.
* Usage:   do "prep/mi_diagnostics.do"   (needs the MRDR drive; reads $data_path/MRDR Long MI.dta,
*          the main-model output of multiple_imputation.do run with $imp set to the reporting M).
* Notes:   Diagnostics only -- reads the imputed data, never rebuilds it. Deliberately kept out of
*          multiple_imputation.do so the production build (including the 500-iteration bootstrap on
*          the HPC) stays lean and its logs uncluttered. FMI and RE only settle as M grows, so read
*          them off the main M=10 model, not a 2-imputation run. These covariates also enter the
*          OS (streg) and ASCT (logit) equations; for the FMI/DF of those COEFFICIENTS as fitted,
*          add `vartable` / `dftable` to the matching `mi estimate` calls in prep/risk_equations.do.
**********

clear all
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"
set linesize 200                  // keep the wide mi vartable/dftable on one line in the log

* Write a plain-text log to scratch/ (local + git-ignored via *.log; the MRDR drive is not readable
* off-Stata) so the results can be reviewed after a drive run. Named log -> coexists with any log
* already open in an interactive session. Convention documented in CLAUDE.md.
cap mkdir "scratch"
cap log close mi_diag
log using "scratch/mi_diagnostics.log", replace text name(mi_diag)

* Baseline covariates entering the risk equations, at the diagnosis record (Event0 == 3, one row
* per patient), grouped by how the FMI is best summarised:
*   binary / continuous -> mean ;  ordinal (>2 levels) -> proportion (FMI per level).
local mean_vars  Male FISHRisk CM_CRD CM_PLM CM_DBT eGFR CM_CKD   // 0/1 + eGFR; CM_CKD passive (from eGFR)
local prop_vars  ECOGcc RISS                                      // ordinal; ECOGcc imputed, RISS passive
* Directly-imputed inputs listed for the missingness table (labs feed ISS/R-ISS):
local miss_vars  Male FISHRisk ECOGcc eGFR CM_CRD CM_PLM CM_DBT Albumin SerumB2Microglobulin LactateDehydrogenase

* -- Pre-imputation missingness at diagnosis (context) --
use "${data_path}/MRDR Long.dta", clear
di as text _n(2) "{hline 74}"
di as text "Percent missing at diagnosis (Event0 == 3), before imputation"
di as text "{hline 74}"
misstable summarize `miss_vars' if Event0 == 3

* -- MI diagnostics at the current M --
use "${data_path}/MRDR Long MI.dta", clear
mi describe

di as text _n(2) "{hline 74}"
di as text "Binary / continuous covariates -- variance information (RVI, FMI, rel. efficiency)"
di as text "{hline 74}"
mi estimate, vartable: mean `mean_vars' if Event0 == 3

di as text _n(2) "{hline 74}"
di as text "Binary / continuous covariates -- degrees of freedom (DF)"
di as text "{hline 74}"
mi estimate, dftable: mean `mean_vars' if Event0 == 3

di as text _n(2) "{hline 74}"
di as text "Ordinal covariates (ECOGcc, R-ISS) -- variance information, FMI per level"
di as text "{hline 74}"
mi estimate, vartable: proportion `prop_vars' if Event0 == 3

* First-line best clinical response (BCR_L1), among patients reaching L1 (Event0 == 10). Reported
* separately: BCR is time-varying (imputed at line starts), not a diagnosis-record covariate.
di as text _n(2) "{hline 74}"
di as text "First-line response (BCR_L1) -- variance information, FMI per level (at L1 start)"
di as text "{hline 74}"
mi estimate, vartable: proportion BCR_L1 if Event0 == 10

di as text _n "Diagnostics log saved to scratch/mi_diagnostics.log"
cap log close mi_diag
