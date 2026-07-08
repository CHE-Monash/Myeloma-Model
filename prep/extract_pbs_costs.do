**********
* Monash Myeloma Model - PBS cost extraction
*
* Purpose: Build the small, committed PBS cost inputs the treatment-cost calculator needs, by
*          filtering a dated full PBS Schedule extract down to just the modelled myeloma drugs.
*          The raw extract (~40 MB, git-ignored, in the sibling data/ folder) stays out of the repo;
*          this script records exactly how the committed subset was produced, so the whole cost
*          chain is reproducible from a re-downloadable, dated source.
* Usage:   do "prep/extract_pbs_costs.do"
* Source:  $pbs_src  = tables_as_csv/ dir of a dated PBS API CSV download (set in config.do).
*          Swappable later for the PBS REST API (data.pbs.gov.au) — keep the source isolated here.
* Inputs:  <pbs_src>/items.csv        - per-item detail (AEMP, form, program, pack/vial content, max amount)
*          <pbs_src>/fees.csv         - per-program dispensing / EFC preparation / container fees
*          <pbs_src>/copayments.csv   - general & concessional patient co-payments
*          <pbs_src>/markup-bands.csv  - wholesale + pharmacy mark-up bands (General Schedule path)
* Output:  prep/inputs/pbs_prices.csv     - AEMP per (drug, strength, program, pack), cheapest brand
*          prep/inputs/pbs_fees.csv       - the fees the DPMQ build uses, per program
*          prep/inputs/pbs_copayments.csv - general / concessional co-payment
*          prep/inputs/pbs_markups.csv    - General-Schedule mark-up bands (dexamethasone, oral cyclo)
* Method:  DPMQ (built downstream in treatment_costs.do) = AEMP + mark-ups + fees, by program:
*            EFC injectables (IN public / IP private) : sum(vial AEMP) + EFC preparation fee
*            Section 100 HSD orals (HB public / HS private): pack AEMP + dispensing fee
*            General Schedule orals (GE)              : pack AEMP + wholesale + pharmacy mark-up + fee
*          Prices are AEMP (ex-manufacturer = claimed_price, else determined_price). s19A temporary
*          shortage-import brands are dropped; within each (drug, strength, program, pack) the
*          cheapest brand is kept (captures generic price disclosure).
* Notes:   Prices are a dated snapshot — record the schedule date with any analysis. Perspective is
*          the Australian health system (full DPMQ, co-payment NOT netted); co-payments are extracted
*          only so a cost-to-government view can be produced on request.
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

* ---- Config ----
local SRC "$pbs_src"
if "`SRC'" == "" local SRC "../data/2026-07-01-PBS-API-CSV-files/tables_as_csv"
local SCHED "2026-07-01"          // schedule date of this extract (provenance stamp)
local OUT   "prep/inputs"

di as text _n "=== Extracting PBS cost inputs from `SRC' (schedule `SCHED') ===" _n

* ===========================================================================
* 1. PRICES  -  AEMP per modelled (drug, strength, program, pack), cheapest brand
* ===========================================================================
import delimited "`SRC'/items.csv", varnames(1) case(preserve) stringcols(_all) clear
gen DRUG = upper(drug_name)
keep if inlist(DRUG,"BORTEZOMIB","CARFILZOMIB","CYCLOPHOSPHAMIDE","DARATUMUMAB") ///
      | inlist(DRUG,"LENALIDOMIDE","POMALIDOMIDE","THALIDOMIDE","DEXAMETHASONE")

* AEMP = ex-manufacturer price: claimed_price, else determined_price
destring claimed_price determined_price, gen(cp dp) force
gen double aemp = cond(!missing(cp), cp, dp)

* drop s19A temporary shortage-import brands (inflated prices, not the standard listing)
drop if ustrregexm(lower(li_form), "s19a")

* strength (mg per vial / capsule / tablet): injectables from vial_content (normalise mcg->mg via
* unit_of_measure, e.g. bortezomib is listed in mcg); orals parsed from the form text
destring vial_content pack_size maximum_amount, gen(vc psize maxamt) force
gen byte is_mcg = inlist(lower(unit_of_measure),"mcg","microgram","micrograms")
gen double strength_mg = cond(is_mcg, vc/1000, vc)
replace maxamt = maxamt/1000 if is_mcg & !missing(maxamt)   // keep max_amount in mg too
replace strength_mg = real(ustrregexs(1))      if missing(strength_mg) & ustrregexm(li_form,"([0-9.]+) mg")
replace strength_mg = real(ustrregexs(1))/1000 if missing(strength_mg) & ustrregexm(li_form,"([0-9.]+) microgram")

* keep only the modelled program x form per drug (public + private paths retained)
gen byte keep_it = 0
* EFC injectables (public IN / private IP)
replace keep_it = 1 if inlist(DRUG,"BORTEZOMIB","CARFILZOMIB") & inlist(program_code,"IN","IP")
replace keep_it = 1 if DRUG=="DARATUMUMAB" & inlist(program_code,"IN","IP") & ustrregexm(li_form,"I\.V\.")
* Section 100 HSD oral immunomodulators (public HB / private HS)
replace keep_it = 1 if inlist(DRUG,"LENALIDOMIDE","POMALIDOMIDE","THALIDOMIDE") & inlist(program_code,"HB","HS")
* General Schedule oral tablets (dexamethasone, cyclophosphamide)
replace keep_it = 1 if inlist(DRUG,"DEXAMETHASONE","CYCLOPHOSPHAMIDE") & program_code=="GE" & ustrregexm(li_form,"Tablet")
keep if keep_it
drop if missing(strength_mg) | missing(aemp)

* cheapest brand within (drug, strength, program, pack)
bysort DRUG strength_mg program_code psize (aemp): keep if _n==1

* tidy output columns
gen drug     = proper(lower(DRUG))
gen setting  = cond(inlist(program_code,"IN","HB"), "public", ///
               cond(inlist(program_code,"IP","HS"), "private", "community"))
gen kind     = cond(inlist(program_code,"IN","IP"), "injectable", "oral")
gen sched    = "`SCHED'"
drop pack_size                       // original string column; use destring'd psize instead
rename (program_code li_form pbs_code brand_name maxamt psize) ///
       (program form code brand max_amount_mg pack_size)
keep  drug strength_mg pack_size aemp program setting kind form max_amount_mg code brand sched
order drug strength_mg pack_size aemp program setting kind form max_amount_mg code brand sched
gsort drug kind program strength_mg pack_size
format aemp %9.2f
export delimited "`OUT'/pbs_prices.csv", replace datafmt
di as result "  wrote `OUT'/pbs_prices.csv (" _N " rows)"

* ===========================================================================
* 2. FEES  -  the fees the DPMQ build uses, per modelled program
* ===========================================================================
import delimited "`SRC'/fees.csv", varnames(1) case(preserve) stringcols(_all) clear
keep if inlist(program_code,"IN","IP","HB","HS","GE")
destring dispensing_fee_ready_prepared efc_preparation_fee efc_diluent_fee efc_distribution_fee, ///
         gen(disp_fee prep_fee dil_fee dist_fee) force
gen setting = cond(inlist(program_code,"IN","HB"), "public", ///
              cond(inlist(program_code,"IP","HS"), "private", "community"))
gen note = cond(inlist(program_code,"IN","IP"), "EFC: sum(vial AEMP) + prep_fee", ///
           cond(inlist(program_code,"HB","HS"), "S100 HSD: pack AEMP + disp_fee", ///
           "General Schedule: pack AEMP + mark-up + disp_fee"))
rename program_code program
keep  program setting disp_fee prep_fee dil_fee dist_fee note
order program setting disp_fee prep_fee dil_fee dist_fee note
gsort program
export delimited "`OUT'/pbs_fees.csv", replace
di as result "  wrote `OUT'/pbs_fees.csv (" _N " rows)"

* ===========================================================================
* 3. CO-PAYMENTS  -  general & concessional (2026)
* ===========================================================================
import delimited "`SRC'/copayments.csv", varnames(1) case(preserve) stringcols(_all) clear
destring general concessional, gen(g c) force
clear
set obs 2
gen str12 type   = cond(_n==1, "general", "concessional")
gen double copay = .
* pull the two scalars back in via a tiny re-read (avoids a reshape on a 1-row table)
preserve
import delimited "`SRC'/copayments.csv", varnames(1) case(preserve) stringcols(_all) clear
destring general concessional, gen(g c) force
scalar gg = g[1]
scalar cc = c[1]
restore
replace copay = gg if type=="general"
replace copay = cc if type=="concessional"
gen sched = "`SCHED'"
format copay %6.2f
export delimited "`OUT'/pbs_copayments.csv", replace datafmt
di as result "  wrote `OUT'/pbs_copayments.csv (" _N " rows)"

* ===========================================================================
* 4. MARK-UP BANDS  -  General Schedule community pharmacy (s90-cp): wholesale + pharmacy
* ===========================================================================
import delimited "`SRC'/markup-bands.csv", varnames(1) case(preserve) stringcols(_all) clear
keep if program_code=="GE" & dispensing_rule_mnem=="s90-cp"
keep if inlist(markup_band_code,"W","C")     // W = wholesale mark-up, C = pharmacy (community) mark-up
destring limit variable offset fixed, gen(lo_limit pct off fixed_amt) force
gen markup = cond(markup_band_code=="W","wholesale","pharmacy")
keep  markup lo_limit pct off fixed_amt
order markup lo_limit pct off fixed_amt
gsort markup lo_limit
di as text "  (DPMQ mark-up per band = pct% * (AEMP + off) + fixed_amt, using the band whose lo_limit <= AEMP)"
export delimited "`OUT'/pbs_markups.csv", replace
di as result "  wrote `OUT'/pbs_markups.csv (" _N " rows)"

di as result _n "=== PBS cost inputs extracted (schedule `SCHED') ==="
