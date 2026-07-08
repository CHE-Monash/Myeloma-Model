**********
* Monash Myeloma Model - Cost index (ABS CPI) builder
*
* Purpose: Write prep/inputs/cost_index.csv from the ABS Consumer Price Index, hardcoded here so the
*          non-drug cost deflator is preserved in code and easy to extend. treatment_costs.do inflates
*          the non-treatment (health-state) costs from their source year to the target cost year by
*          index[target]/index[source] (carrying the latest year <= target forward). Drug prices are
*          NOT inflated - they use actual dated PBS values.
* Usage:   do "prep/build_cost_index.do"
* Series:  ABS 6401.0 Consumer Price Index, Australia (table 610101) - "Index Numbers ; All groups CPI ;
*          Australia ;", series ID A2325846C, ORIGINAL, index reference base 2011-12 = 100. Values below
*          are the JUNE-quarter index for each year (Yap 2025 non-treatment costs are in June-2019 AUD).
*          To add a future year: append one row to the input block and re-run.
* Output:  prep/inputs/cost_index.csv  (year, index, note)
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
local OUT "prep/inputs"

* ABS CPI All groups Australia (A2325846C), June quarter, index (2011-12 = 100)
input int year double index
2015 107.5
2016 108.6
2017 110.7
2018 113.0
2019 114.8
2020 114.4
2021 118.8
2022 126.1
2023 133.7
2024 138.8
2025 141.7
end

gen str90 note = "ABS CPI All groups Australia (series A2325846C, 6401.0/610101), June quarter, base 2011-12=100"
replace note = note + "; Yap 2025 non-treatment source year" if year == 2019
replace note = note + "; latest target (2019->2025 factor 141.7/114.8 = 1.234)" if year == 2025

format index %6.1f
export delimited "`OUT'/cost_index.csv", replace datafmt
di as result _n "Wrote `OUT'/cost_index.csv (" _N " years, ABS A2325846C)"
list year index, noobs sep(0)
