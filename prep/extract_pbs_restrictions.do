**********
* Monash Myeloma Model - PBS multiple-myeloma restriction map (cost-engine Stage C)
*
* Purpose: Map each modelled drug to its multiple-myeloma PBS restrictions (eligibility), so the
*          model's regimen/line assignments can be checked against actual PBS listing rules and the
*          per-line item/max-amount can be pinned. Companion to prep/extract_pbs_costs.do (prices);
*          this extracts the *eligibility* side (condition, treatment phase, combination, line signal).
* Usage:   do "prep/extract_pbs_restrictions.do"
* Source:  $pbs_src  = tables_as_csv/ dir of a dated PBS API CSV download (set in config.do).
* Inputs:  <pbs_src>/items.csv                        - drug <-> pbs_code, max amount
*          <pbs_src>/item-restriction-relationships.csv - pbs_code <-> res_code
*          <pbs_src>/restrictions.csv                  - res_code -> treatment phase, authority, criteria text
* Output:  prep/inputs/pbs_restrictions.csv - one row per (drug, restriction): condition, phase, combination
*          partners, line signal, authority method, a likely modelled regimen, and a trimmed criteria summary.
* Method:  Keep restrictions whose text mentions "myeloma". Parse the (HTML-stripped) restriction text for
*          condition (relapsed/refractory -> RRMM; previously untreated/newly diagnosed -> newly-dx),
*          combination partners ("in combination with ..."), and a line signal (prior-therapy wording).
*          `likely_regimen` is a best-effort label from drug + combination — verify before relying on it.
* Notes:   Eligibility reference only; the priced items + max amounts live in pbs_prices_<year>.csv. A restriction
*          typically spans several items (strengths/programs) — `n_items` counts how many of the drug's
*          priced codes it governs.
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

local SRC "$pbs_src"
if "`SRC'" == "" local SRC "../data/2026-07-01-PBS-API-CSV-files/tables_as_csv"
local SCHED "2026-07-01"
local OUT   "prep/inputs"

di as text _n "=== Extracting PBS multiple-myeloma restrictions from `SRC' (schedule `SCHED') ===" _n

* ---- 1. Our drugs' items: pbs_code -> drug ----
import delimited "`SRC'/items.csv", varnames(1) case(preserve) stringcols(_all) clear
gen DRUG = upper(drug_name)
keep if inlist(DRUG,"BORTEZOMIB","CARFILZOMIB","CYCLOPHOSPHAMIDE","DARATUMUMAB") ///
      | inlist(DRUG,"LENALIDOMIDE","POMALIDOMIDE","THALIDOMIDE","DEXAMETHASONE")
gen drug = proper(lower(DRUG))
keep pbs_code drug
duplicates drop
tempfile items
save `items'

* ---- 2. pbs_code -> res_code (restrict to our codes) ----
import delimited "`SRC'/item-restriction-relationships.csv", varnames(1) case(preserve) stringcols(_all) clear
keep pbs_code res_code
merge m:1 pbs_code using `items', keep(match) nogen      // keep only our drugs' codes
tempfile itemres
save `itemres'

* ---- 3. Restrictions -> keep multiple myeloma, parse ----
import delimited "`SRC'/restrictions.csv", varnames(1) case(preserve) stringcols(_all) clear
keep res_code treatment_phase authority_method li_html_text
merge 1:m res_code using `itemres', keep(match) nogen     // one row per (res_code, pbs_code, drug)

* keep MM restrictions
gen lc = lower(li_html_text)
keep if strpos(lc, "myeloma")

* plain-text criteria (strip HTML tags, collapse whitespace)
gen plain = ustrregexra(li_html_text, "<[^>]+>", " ")
replace plain = ustrregexra(plain, "[ ]+", " ")
replace plain = strtrim(plain)
gen plc = lower(plain)

* condition
gen condition = cond(strpos(plc,"relapsed")|strpos(plc,"refractory"), "RRMM", ///
                cond(strpos(plc,"previously untreated")|strpos(plc,"newly diagnosed"), "newly-dx", "?"))

* combination partners ("in combination with <...>" up to a stop)
gen combination = ""
replace combination = strtrim(ustrregexs(1)) if ustrregexm(plc, "combination with ([a-z ]+?)(,| for | in | as |;|\.| the | patient| this )")
replace combination = subinstr(combination, " and", " +", .)

* line signal: prefer the (authoritative) treatment-phase wording, then fall back to prior-therapy text
gen lph = lower(treatment_phase)
gen line_signal = ""
replace line_signal = "L1 (untreated/1st-line)" if strpos(lph,"first line")|strpos(lph,"first-line")
replace line_signal = "L2 (2nd-line)"           if line_signal=="" & (strpos(lph,"second line")|strpos(lph,"second-line"))
replace line_signal = "L3 (3rd-line)"           if line_signal=="" & (strpos(lph,"third line")|strpos(lph,"third-line"))
replace line_signal = "L1 (untreated/1st-line)" if line_signal=="" & (strpos(plc,"previously untreated")|strpos(plc,"newly diagnosed"))
replace line_signal = "L3+ (>=2 prior)"          if line_signal=="" & (strpos(plc,"at least 2 prior")|strpos(plc,"at least two prior"))
replace line_signal = "L2+ (>=1 prior)"          if line_signal=="" & strpos(plc,"at least one prior")

* likely modelled regimen (best-effort — verify)
gen likely_regimen = ""
replace likely_regimen = "Kd"          if drug=="Carfilzomib"  & strpos(combination,"dexamethasone") & !strpos(combination,"lenalidomide")
replace likely_regimen = "KRd (n/m)"   if drug=="Carfilzomib"  & strpos(combination,"lenalidomide")
replace likely_regimen = "DVd"         if drug=="Daratumumab"  & strpos(combination,"bortezomib")
replace likely_regimen = "Dara 1L (n/m)" if drug=="Daratumumab" & likely_regimen=="" & condition=="newly-dx"
replace likely_regimen = "Pd"          if drug=="Pomalidomide" & strpos(combination,"dexamethasone")
replace likely_regimen = "Pom-triple (n/m)" if drug=="Pomalidomide" & likely_regimen==""
replace likely_regimen = "Rd / R-maint" if drug=="Lenalidomide"
replace likely_regimen = "V-based (VCd/VRd/Vd/VTd)" if drug=="Bortezomib"
replace likely_regimen = "T-based (VTd/TCd/Td)"     if drug=="Thalidomide"

* trimmed criteria summary (first 220 chars)
gen criteria_summary = substr(plain, 1, 220)

* ---- 4. Collapse to one row per (drug, res_code); count items governed ----
bysort drug res_code: gen n_items = _N
bysort drug res_code: keep if _n==1
gen sched = "`SCHED'"
keep  drug condition line_signal treatment_phase likely_regimen combination authority_method n_items res_code criteria_summary sched
order drug condition line_signal treatment_phase likely_regimen combination authority_method n_items res_code criteria_summary sched
gsort drug condition treatment_phase res_code  // res_code tie-break = reproducible row order
export delimited "`OUT'/pbs_restrictions.csv", replace
di as result _n "Wrote `OUT'/pbs_restrictions.csv (" _N " restrictions)"
