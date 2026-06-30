**********
* Monash Myeloma Model - OOS (70/30): bootstrap prediction-interval validation
*
* Purpose: The HEADLINE out-of-sample metric. validate_oos.do checks the POINT estimate (one
*          70%-trained coefficient set) against the held-out 30%'s observed outcomes with fixed
*          tolerances. This script instead builds a 95% PREDICTION INTERVAL for each outcome from the
*          500 bootstrap simulations (each a patient-cluster resample of the 70%, re-imputed + re-fit
*          + re-simulated on the held-out 30%) and asks, by the PERCENTILE METHOD: does the held-out
*          OBSERVED value fall inside the bootstrap 95% interval [p2.5, p97.5]? -- the calibration check
*          used in the 2024 PLOS ONE paper. Coverage (% of targets inside) should approach ~95% if the
*          model is well-calibrated.
*
* Outcomes (identical estimators to validate_outcomes.do, computed per bootstrap):
*   OS  - KM survival at 3 & 5 yr, by line (L1 no-ASCT, ASCT, L2, L3) and BCR
*   BCR - response distribution (%) by line
*   TXD - % still on treatment at 12 & 24 mo, by line and BCR
*   TFI - % still treatment-free at 12 & 24 mo, by line and BCR
*   PATH- ASCT rate (among L1-end reachers) and L2-L5 reach rates
*   (the TFI-L2 response-ORDERING check is a pass/fail, not a target value, so it is not covered here)
*
* ---------------------------------------------------------------------------------------------------
* DESIGNED TO RUN ON THE HPC (the 500 bootstrap .dta live there). To run it there, copy:
*     analyses/oos/bootstrap_validation.do          (this file)
*     analyses/oos/targets/                          (the 13 observed-target CSVs -- movable)
*     analyses/oos/simulated/bootstrap         (the folder of 500 bootstrap _B<n>.dta simulations)
*   then, from the repository root:  stata -b do analyses/oos/bootstrap_validation.do
*   Positional args override the paths (see usage below); edit the defaults if you prefer. SANITY-CHECK
*   FIRST with a small run:  stata -b do analyses/oos/bootstrap_validation.do 5
*   Pull back:                       analyses/oos/results/oos_bootstrap_validation.md   (review)
*                                    analyses/oos/results/oos_bootstrap_coverage.csv    (machine-readable)
* NB: not run locally at authoring time (no bootstrap files present) -- smoke-test with bv_maxbs(5).
* ---------------------------------------------------------------------------------------------------
**********

* ---- optional positional args, read into LOCALS before clear all (which would wipe globals) ----
* usage:  do bootstrap_validation.do [maxbs] [targets_dir] [simdir] [out_dir] [minbs] [stub]
*   e.g.  do analyses/oos/bootstrap_validation.do 5                 // 5-resample sanity check
*         do analyses/oos/bootstrap_validation.do 500 mytargets mysims myout   // HPC paths
local a_maxbs  `"`1'"'
local a_targ   `"`2'"'
local a_simdir `"`3'"'
local a_out    `"`4'"'
local a_minbs  `"`5'"'
local a_stub   `"`6'"'

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes home
capture run "config.do"                  // machine-specific paths (git-ignored; harmless if absent on HPC)

**********
* Configuration -- positional args override the defaults
**********
global bv_targets "analyses/oos/targets"                      // observed-target CSVs (movable to HPC)
global bv_simpath  "analyses/oos/simulated/bootstrap"          // the 500 bootstrap sims (live on HPC)
global bv_stub    "all_0_oos_1_101212"                        // file = <stub>_B<b>.dta
global bv_minbs   "1"
global bv_maxbs   "500"
global bv_out     "analyses/oos/results"
if `"`a_maxbs'"'  != "" global bv_maxbs  `"`a_maxbs'"'
if `"`a_targ'"'   != "" global bv_targets `"`a_targ'"'
if `"`a_simdir'"' != "" global bv_simpath `"`a_simdir'"'
if `"`a_out'"'    != "" global bv_out    `"`a_out'"'
if `"`a_minbs'"'  != "" global bv_minbs  `"`a_minbs'"'
if `"`a_stub'"'   != "" global bv_stub   `"`a_stub'"'
capture mkdir "$bv_out"

di "{hline 90}"
di "OOS BOOTSTRAP PREDICTION-INTERVAL VALIDATION  (percentile method, 95% interval)"
di "Resamples: $bv_minbs .. $bv_maxbs   |  targets: $bv_targets   |  sims: $bv_simpath"
di "{hline 90}"

**********
* Load observed targets (held-out 30%) into matrices -- same import block as validate_outcomes.do
**********
import delimited "${bv_targets}/os_l1_noasct.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L1_NoASCT_bench)
import delimited "${bv_targets}/os_asct.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_ASCT_bench)
import delimited "${bv_targets}/os_l2.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L2_bench)
import delimited "${bv_targets}/os_l3.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L3_bench)
import delimited "${bv_targets}/bcr.csv", clear case(preserve)
mkmat N CR VG PR MR SD PD, matrix(BCR_bench)
import delimited "${bv_targets}/txd_l1_noasct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TXD_L1_NoASCT_bench)
import delimited "${bv_targets}/txd_l1_asct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TXD_L1_ASCT_bench)
import delimited "${bv_targets}/tfi_l1_noasct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L1_NoASCT_bench)
import delimited "${bv_targets}/tfi_l1_asct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L1_ASCT_bench)
import delimited "${bv_targets}/tfi_l2.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L2_bench)
import delimited "${bv_targets}/pathways.csv", clear case(preserve)
mkmat N ASCT L2 L3 L4 L5 L6 L7 L8, matrix(Pathways_bench)
di "Targets loaded."

**********
* Helpers (matrices are global, so they are visible inside these programs)
**********

* OS: KM survival at 3 & 5 yr for each BCR group, clocked from a line's start (OC_TIME - <tsd>).
cap program drop bv_os
program define bv_os
    args h b tsd gvar nbcr bench lab sct
    if "`sct'" == "0"      local cond "& SCT_L1 == 0"
    else if "`sct'" == "1" local cond "& SCT_L1 == 1"
    else                   local cond ""
    capture drop _ostime
    qui gen double _ostime = OC_TIME - `tsd'
    qui stset _ostime, failure(OC_MORT==1) id(ID)
    forvalues bcr = 1/`nbcr' {
        foreach yr in 3 5 {
            local col = cond(`yr'==3, 5, 7)
            local mo  = `yr' * 12
            local obs = `bench'[`bcr', `col'] * 100
            capture drop _osv
            qui sts generate _osv = s if `gvar' == `bcr' `cond'
            qui summarize _osv if `gvar' == `bcr' `cond' & _t <= `mo'
            if r(N) > 0 & !missing(`obs') {
                post `h' (`b') ("OS") ("OS/`lab'/BCR`bcr'/`yr'yr") (`obs') (`=r(min)*100')
            }
            capture drop _osv
        }
    }
end

* Horizon survival (TXD/TFI): % still on-treatment / treatment-free at 12 & 24 mo, by BCR group.
cap program drop bv_hz
program define bv_hz
    args h b fam tvar fail gvar nbcr bench lab sct
    if "`sct'" == "0"      local cond "& SCT_L1 == 0"
    else if "`sct'" == "1" local cond "& SCT_L1 == 1"
    else                   local cond ""
    qui stset `tvar', failure(`fail') id(ID)
    forvalues bcr = 1/`nbcr' {
        local c = 7
        foreach mo in 12 24 {
            local obs = `bench'[`bcr', `c']
            capture drop _hv
            qui sts generate _hv = s if `gvar' == `bcr' `cond'
            qui summarize _hv if `gvar' == `bcr' `cond' & _t <= `mo'
            if r(N) > 0 & !missing(`obs') {
                post `h' (`b') ("`fam'") ("`fam'/`lab'/BCR`bcr'/`mo'mo") (`obs') (`=r(min)*100')
            }
            capture drop _hv
            local ++c
        }
    }
end

* BCR distribution (%) for a line.
cap program drop bv_bcr
program define bv_bcr
    args h b gvar nbcr brow lab
    qui count if !missing(`gvar')
    local n = r(N)
    if `n' == 0 exit
    forvalues bcr = 1/`nbcr' {
        local obs = BCR_bench[`brow', `bcr'+1]
        qui count if `gvar' == `bcr'
        local sim = r(N) / `n' * 100
        if !missing(`obs') post `h' (`b') ("BCR") ("BCR/`lab'/`bcr'") (`obs') (`sim')
    }
end

**********
* 1. Per-bootstrap: load the resample's simulated 30% and post every outcome  (long: b, family, metric, obs, sim)
**********
tempname pf
postfile `pf' int b str8 family str44 metric double obs double sim using "$bv_out/bootstrap_iterations.dta", replace

local nfiles = 0
forvalues b = $bv_minbs/$bv_maxbs {
    capture use "$bv_simpath/${bv_stub}_B`b'.dta", clear
    if _rc {
        di as error "  missing bootstrap file B`b' -- skipped"
        continue
    }
    local ++nfiles

    // OS (4 lines)
    bv_os `pf' `b' TSD_L1S BCR_L1  6 OS_L1_NoASCT_bench L1NA 0
    bv_os `pf' `b' TSD_L1E BCR_SCT 4 OS_ASCT_bench      ASCT none
    bv_os `pf' `b' TSD_L2S BCR_L2  6 OS_L2_bench        L2   none
    bv_os `pf' `b' TSD_L3S BCR_L3  6 OS_L3_bench        L3   none

    // BCR distribution (4 lines)
    bv_bcr `pf' `b' BCR_L1  6 1 L1
    bv_bcr `pf' `b' BCR_SCT 4 2 ASCT
    bv_bcr `pf' `b' BCR_L2  6 3 L2
    bv_bcr `pf' `b' BCR_L3  6 4 L3

    // TXD (% on treatment) -- L1 NoASCT & ASCT share stset TXD_L1, failure(MOR_L1S==0)
    bv_hz `pf' `b' TXD TXD_L1 MOR_L1S==0 BCR_L1  6 TXD_L1_NoASCT_bench L1NA 0
    bv_hz `pf' `b' TXD TXD_L1 MOR_L1S==0 BCR_SCT 4 TXD_L1_ASCT_bench   L1AS 1

    // TFI (% treatment-free)
    bv_hz `pf' `b' TFI TFI_L1 MOR_L1E==0 BCR_L1  6 TFI_L1_NoASCT_bench L1NA 0
    bv_hz `pf' `b' TFI TFI_L1 MOR_L1E==0 BCR_SCT 4 TFI_L1_ASCT_bench   L1AS 1
    bv_hz `pf' `b' TFI TFI_L2 MOR_L2E==0 BCR_L2  6 TFI_L2_bench        L2   none

    // PATHWAYS -- ASCT among L1-end reachers; L2-L5 reach among all
    qui count if OC_TIME > TSD_L1E & !missing(TSD_L1E)
    local nl1e = r(N)
    qui count if SCT_L1 == 1
    local sct_n = r(N)
    if `nl1e' > 0 {
        local obs = Pathways_bench[1, 2]
        post `pf' (`b') ("PATH") ("PATH/ASCT") (`obs') (`=`sct_n'/`nl1e'*100')
    }
    qui count
    local nt = r(N)
    forvalues line = 2/5 {
        qui count if !missing(BCR_L`line')
        local sim = r(N) / `nt' * 100
        local obs = Pathways_bench[1, `line'+1]
        post `pf' (`b') ("PATH") ("PATH/L`line'reach") (`obs') (`sim')
    }

    if (`b' == $bv_minbs | mod(`b', 50) == 0) di as text "  ... processed resample `b'"
}
postclose `pf'
di as text _n "Processed `nfiles' bootstrap file(s)."

**********
* 2. Per-target 95% percentile interval + coverage check (observed inside [p2.5, p97.5]?)
**********
use "$bv_out/bootstrap_iterations.dta", clear

* Only assess targets populated in ALL processed resamples (`nfiles'); cells empty/degenerate in
* some bootstraps are flagged included==0 and excluded from coverage (listed in the report).
tempname cf
postfile `cf' str8 family str44 metric double observed double lo95 double med double hi95 long nsim byte inside byte included ///
    using "$bv_out/oos_bootstrap_coverage.dta", replace

quietly levelsof metric, local(mets)
foreach mm of local mets {
    qui _pctile sim if metric == "`mm'", p(2.5 50 97.5)
    local lo  = r(r1)
    local med = r(r2)
    local hi  = r(r3)
    qui summarize obs if metric == "`mm'", meanonly
    local ob  = r(mean)
    qui count if metric == "`mm'" & !missing(sim)
    local ns  = r(N)
    qui levelsof family if metric == "`mm'", local(fm) clean
    local ins = (`ob' >= `lo' & `ob' <= `hi')
    local inc = (`ns' == `nfiles')
    post `cf' ("`fm'") ("`mm'") (`ob') (`lo') (`med') (`hi') (`ns') (`ins') (`inc')
}
postclose `cf'

**********
* 3. Report -- write a readable .md (the pull-back artifact) + machine-readable .csv
**********
use "$bv_out/oos_bootstrap_coverage.dta", clear
gsort family metric

qui count if included
local ntot = r(N)
qui count if included & inside
local ncov = r(N)
local opct = cond(`ntot' > 0, 100*`ncov'/`ntot', .)
qui count if !included
local nexcl = r(N)

tempname rf
file open `rf' using "$bv_out/oos_bootstrap_validation.md", write replace
file write `rf' "# OOS bootstrap prediction-interval validation" _n _n
file write `rf' "Percentile method: 95% interval [p2.5, p97.5] of the simulated outcome across `nfiles' bootstrap resample(s) (B$bv_minbs..B$bv_maxbs). Each target asks whether the **held-out observed** value falls inside that interval; overall coverage near ~95% indicates good calibration." _n _n
file write `rf' "Only targets populated in **all `nfiles'** resamples are assessed (`nexcl' sparse cell(s) excluded -- listed at the end). Outcomes run to L3 (OS/BCR/TFI) and L5 (pathways reach); nothing beyond line 6." _n _n
file write `rf' "Targets: \`$bv_targets'  |  Sims: \`$bv_simpath'" _n _n
file write `rf' "## Overall coverage" _n _n
file write `rf' "**`ncov' / `ntot' assessed targets inside the bootstrap 95% interval (`=string(`opct',"%4.1f")'%)**" _n _n
file write `rf' "_~95% = well calibrated; well below = intervals too narrow or systematic bias; well above = intervals too wide / over-dispersed._" _n _n
file write `rf' "## Coverage by family (assessed targets)" _n _n
file write `rf' "| Family | Inside | Assessed | Coverage |" _n
file write `rf' "|---|---|---|---|" _n
foreach fm in OS BCR TXD TFI PATH {
    qui count if family == "`fm'" & included
    local ft = r(N)
    if `ft' > 0 {
        qui count if family == "`fm'" & included & inside
        local fc = r(N)
        file write `rf' "| `fm' | `fc' | `ft' | `=string(100*`fc'/`ft',"%4.1f")'% |" _n
    }
}
file write `rf' _n "## Per-target intervals (assessed)" _n _n
file write `rf' "| Family | Target | Observed | lo95 | hi95 | Median | n | Inside? |" _n
file write `rf' "|---|---|---|---|---|---|---|---|" _n
forvalues i = 1/`=_N' {
    if included[`i'] {
        local ins = cond(inside[`i'], "yes", "**NO**")
        file write `rf' "| `=family[`i']' | `=metric[`i']' | `=string(observed[`i'],"%5.1f")' | `=string(lo95[`i'],"%5.1f")' | `=string(hi95[`i'],"%5.1f")' | `=string(med[`i'],"%5.1f")' | `=nsim[`i']' | `ins' |" _n
    }
}
if `nexcl' > 0 {
    file write `rf' _n "## Excluded (sparse -- not populated in all `nfiles' resamples)" _n _n
    file write `rf' "| Family | Target | n |" _n
    file write `rf' "|---|---|---|" _n
    forvalues i = 1/`=_N' {
        if !included[`i'] file write `rf' "| `=family[`i']' | `=metric[`i']' | `=nsim[`i']' |" _n
    }
}
file close `rf'

* Machine-readable copy for figures ('included' flags the assessed targets)
export delimited family metric observed lo95 med hi95 nsim inside included ///
    using "$bv_out/oos_bootstrap_coverage.csv", replace

di _n "{hline 90}"
di "OVERALL COVERAGE: " %3.0f `ncov' " / " %3.0f `ntot' " assessed inside the bootstrap 95% interval  (" %5.1f `opct' "%)   [`nexcl' sparse excluded]"
di "{hline 90}"
di "Outputs in $bv_out :"
di "  oos_bootstrap_validation.md    (readable report -- pull this back)"
di "  oos_bootstrap_coverage.csv     (per-target intervals + coverage)"
di "  bootstrap_iterations.dta       (per-resample long data, retained)"
