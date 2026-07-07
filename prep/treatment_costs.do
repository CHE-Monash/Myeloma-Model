**********
* Monash Myeloma Model - Treatment cost calculator
*
* Purpose: Build the per-cycle regimen drug costs and the (inflated) non-drug costs the economic
*          model uses, for a chosen price year, from tidy input tables. Replaces the manual
*          "treatment costs" spreadsheet: the derivation is now reproducible and versioned.
* Usage:   do "prep/treatment_costs.do" [target_year]      (default 2025)
* Inputs:  prep/inputs/treatment_regimens.csv  - dosing spec (STABLE: dose, basis, schedule, cycles)
*          prep/inputs/drug_prices.csv         - PBS DPMA per (drug, strength) x YEAR (add a year = append rows)
*          prep/inputs/other_costs.csv         - non-drug costs (ASCT, hospital, community, emergency) x source year
*          prep/inputs/cost_index.csv          - year x ABS price index (for inflating the non-drug costs)
* Output:  prep/inputs/treatment_costs_<year>.csv  - the c* per-cycle drug costs + inflated non-drug costs
* Method:  drug cost = actual PBS price for the target year (carry the latest price <= target forward if a
*          year is missing); non-drug costs = source-year value x ABS index[target]/index[source].
*          Per-drug per-cycle = (DPMA/units_per_pack) * vials_per_admin * admins / cycles, where
*          vials_per_admin = round(dose*bodysize/strength) [Excel-style half-up], bodysize = BSA (/m2),
*          weight (/kg) or 1 (flat). Regimen cost = sum over its drugs. "Other" = mean(VTd,TCd,Td,Vd);
*          maintenance = usage-weighted blend of lenalidomide (R) and thalidomide (T).
* Notes:   Does NOT modify core/process_data.do (which still carries the hardcoded locals). This script
*          regenerates and VALIDATES those values; wiring process_data to read the output CSV is a
*          separate change. Body-size is a population average (see below); per-patient dosing is a
*          possible future extension (matters mainly for carfilzomib /m2 and daratumumab /kg).
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

* ---- Config ----
local target = "`1'"
if "`target'" == "" local target 2025
local BSA    = 1.93        // population-average body surface area (m^2), MRDR
local WEIGHT = 81.10597    // population-average weight (kg), MRDR
local IN     "prep/inputs"

di as text _n "=== Treatment costs for price year `target' ===" _n

* ---- 1. Drug prices for the target year (actual PBS; carry latest <= target forward) ----
import delimited "`IN'/drug_prices.csv", varnames(1) case(preserve) clear
keep if year <= `target'
gsort drug strength_mg -year
by drug strength_mg: keep if _n == 1
keep drug strength_mg units_per_pack year dpma
rename (year dpma) (price_year dpma_used)
tempfile prices
save `prices'

* ---- 2. Dosing spec + join price by (drug, strength) ----
import delimited "`IN'/treatment_regimens.csv", varnames(1) case(preserve) clear
capture destring admins_override units_override, replace force
merge m:1 drug strength_mg using `prices', keep(master match)
quietly count if _merge == 1
if r(N) > 0 {
    di as error "treatment_costs: `r(N)' drug/strength rows have no price at or before `target':"
    list regimen drug strength_mg if _merge == 1, noobs
    exit 459
}
drop _merge

* ---- 3. Per-drug per-cycle cost ----
gen double bodysize = cond(basis == "m2", `BSA', cond(basis == "kg", `WEIGHT', 1))
gen double vials    = cond(!missing(units_override), units_override, floor(dose*bodysize/strength_mg + 0.5))
gen double admins   = cond(!missing(admins_override), admins_override, freq_per_cycle*cycles)
gen double per_cycle = (dpma_used/units_per_pack)*vials*admins/cycles

* ---- 4. Collapse to per-regimen, then the derived "Other" and maintenance blends ----
collapse (sum) per_cycle, by(regimen)
foreach rg in VCd VRd Rd Kd DVd Pd Vd VTd TCd Td R T {
    quietly summarize per_cycle if regimen == "`rg'", meanonly
    scalar c_`rg' = r(mean)
}
scalar c_Other = (c_VTd + c_TCd + c_Td + c_Vd)/4          // pooled "other" regimen (code 0)
scalar c_MNT   = (1002/1504)*c_R + (502/1504)*c_T          // maintenance: MRDR usage-weighted R/T

* ---- 5. Non-drug costs, inflated source-year -> target via the ABS index ----
import delimited "`IN'/cost_index.csv", varnames(1) case(preserve) clear
quietly summarize index if year == `target', meanonly
scalar idx_t = r(mean)
if missing(idx_t) {
    di as text "  (no ABS index for `target'; non-drug costs left un-inflated)"
    scalar idx_t = .
}
tempfile idx
keep year index
rename (year index) (source_year idx_s)
save `idx'

import delimited "`IN'/other_costs.csv", varnames(1) case(preserve) clear
merge m:1 source_year using `idx', keep(master match) nogen
gen double factor = cond(missing(idx_t) | missing(idx_s), 1, idx_t/idx_s)
gen double cost_infl = cost*factor
levelsof component, local(comps) clean
foreach cp of local comps {
    quietly summarize cost_infl if component == "`cp'", meanonly
    scalar `cp' = r(mean)
}

* ---- 6. Validation against the 2025 reference (proves the port reproduces the spreadsheet) ----
if `target' == 2025 {
    di as text _n "Validation vs 2025 reference (Dexa=17.39, Kd 2-vial fix):"
    * name  reference  (tolerance 0.5)
    local refs "VCd 898.02 VRd 1774.96 Rd 1606.55 Kd 25028.49 DVd 12109.93 Pd 2289.08 Vd 724.03 Other 1610.67 MNT 1328.62"
    local nfail 0
    tokenize "`refs'"
    while "`1'" != "" {
        local got = c_`1'
        local d = abs(`got' - `2')
        local ok = cond(`d' < 0.5, "OK", "FAIL")
        if "`ok'" == "FAIL" local ++nfail
        di as text "  c`1' " %10.2f `got' "  ref " %10.2f `2' "   `ok'"
        macro shift 2
    }
    if `nfail' == 0 di as result "  All regimen costs reproduce the reference."
    else            di as error  "  `nfail' mismatch(es) - investigate before use."
}

* ---- 7. Write the output cost set for this year ----
clear
set obs 13
gen str12 parameter = ""
gen double value = .
gen str8 unit = ""
gen str40 note = ""
local i 0
foreach p in cVCd cVRd cRd cKd cDVd cPd cVd cOther cMNT {
    local ++i
    local nm = subinstr("`p'", "c", "", 1)
    replace parameter = "`p'" in `i'
    replace value = c_`nm' in `i'
    replace unit = "AUD/cyc" in `i'
}
foreach p in cASCT cHosp cComm cEmer {
    local ++i
    replace parameter = "`p'" in `i'
    replace value = `p' in `i'
    replace unit = cond("`p'"=="cASCT","AUD","AUD/yr") in `i'
    replace note = "inflated to `target'" in `i'
}
format value %12.2f
export delimited "`IN'/treatment_costs_`target'.csv", replace
di as result _n "Wrote `IN'/treatment_costs_`target'.csv"
list parameter value unit, noobs sep(0)
