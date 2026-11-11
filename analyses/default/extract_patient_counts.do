**********
* Monash Myeloma Model - patient counts for the Ig demand interface
*
* Purpose: Collapse each simulation to the small table the Ig demand model consumes -- patients
*          STARTING each line, patients PREVALENT on each line at 1 January, and three diagnostic
*          "on the line at some point in the year" measures, by calendar year -- and write ONE
*          tidy CSV. The 500 bootstrap .dta stay on the HPC; only the CSV is pulled back.
*
* Usage:   stata -b do analyses/default/extract_patient_counts.do [maxbs] [simdir] [outdir] [minbs] [stub] [y0] [y1]
*   e.g.   stata -b do analyses/default/extract_patient_counts.do 0      // point estimate only (b = 0)
*          stata -b do analyses/default/extract_patient_counts.do 5      // smoke test, B1-B5 first
*          stata -b do analyses/default/extract_patient_counts.do 500    // the full run
*          Run from the REPO ROOT (the script refuses to run from anywhere else).
*
* Output:  $pc_out/patient_counts_bootstrap.csv   -- one row per (b, line, year), when resamples
*          are requested. With maxbs = 0 (point estimate only) the file is
*          $pc_out/patient_counts_point.csv instead, so a point-only run never overwrites a
*          bootstrap extraction that took 500 HPC simulations to produce.
*          b = 0 is the point-estimate run, processed through the IDENTICAL code so the extraction
*          can be reconciled against the committed patients_myeloma_model.csv before any bootstrap
*          number is believed.
*
* Measures (per line k and calendar year y):
*   starting   YearLkS == y: patients starting line k in y. Verified to reproduce the committed
*              interface file exactly (11 Sep 2026, all 99 cells).
*   prevalent  ON line k ON 1 JANUARY of y: DateLkS <= 1 Jan y and end of the TREATMENT spell
*              >= 1 Jan y, where the spell end is DateLkE, or DateMOR if DateLkE is missing, or
*              unbounded if both are missing. This is the interface measure the Ig demand model
*              reads as the prevalent stock at listing. It is a point-in-time count, deliberately
*              on the treatment spell rather than the line state: a patient in the treatment-free
*              interval after line k is counted again when they start line k+1 in the incident
*              flow, so including them here would double-count. A patient on treatment at 1 Jan
*              who dies later that year was genuinely prevalent at 1 Jan, which is why a missing
*              treatment end is filled by the date of death rather than dropped.
*   on_tx      on treatment at any point in the year: DateLkS <= 31 Dec & DateLkE >= 1 Jan
*   on_tx_mor  as on_tx, with a missing treatment end filled by the date of death
*   on_state   in the line-k state at year granularity: YearLk <= y, next line not yet
*              started, and not dead before y
*   The last three are diagnostics kept from the search for the definition behind the retired
*   hand-built `on_line` rows, which none of them reproduced; nothing downstream reads them.
*
* Note:    only the ~29 variables needed are read from each file, not all 260, so each load is fast.
**********

* ---- positional args into LOCALS before clear all (which would wipe globals) ----
local a_maxbs  `"`1'"'
local a_simdir `"`2'"'
local a_out    `"`3'"'
local a_minbs  `"`4'"'
local a_stub   `"`5'"'
local a_y0     `"`6'"'
local a_y1     `"`7'"'

clear all
set more off
set linesize 200

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"                  // machine-specific paths; harmless if absent on the HPC

**********
* Configuration -- positional args override the defaults
**********
global pc_simpath "analyses/default/simulated/bootstrap"        // the 500 bootstrap sims
global pc_point   "analyses/default/simulated/all_0_synthetic.dta"   // point estimate -> b = 0
global pc_stub    "all_0_synthetic"                             // file = <stub>_B<b>.dta
global pc_minbs   "1"
global pc_maxbs   "500"
global pc_out     "analyses/default/results"
global pc_y0      "2026"
global pc_y1      "2036"
if `"`a_maxbs'"'  != "" global pc_maxbs   `"`a_maxbs'"'
if `"`a_simdir'"' != "" global pc_simpath `"`a_simdir'"'
if `"`a_out'"'    != "" global pc_out     `"`a_out'"'
if `"`a_minbs'"'  != "" global pc_minbs   `"`a_minbs'"'
if `"`a_stub'"'   != "" global pc_stub    `"`a_stub'"'
if `"`a_y0'"'     != "" global pc_y0      `"`a_y0'"'
if `"`a_y1'"'     != "" global pc_y1      `"`a_y1'"'

* Output file: the bootstrap table when resamples are requested, a separate point-estimate table
* when they are not, so `... 0` cannot clobber a 500-resample extraction.
if $pc_maxbs < $pc_minbs global pc_file "patient_counts_point"
else                     global pc_file "patient_counts_bootstrap"

* ---- refuse to run from anywhere but the repo root -------------------------
* Every path below is relative to the repo root, and `capture run "config.do"`
* above silently does nothing if you are not there. Without this check, running
* from analyses/default/ resolves $pc_out to analyses/default/analyses/default/
* results, and the level-by-level mkdir below would cheerfully create it.
capture confirm file "analyses/default/simulate.do"
if _rc {
    di as error "{hline 90}"
    di as error "Not at the repo root: analyses/default/simulate.do was not found."
    di as error "Working directory is: `c(pwd)'"
    di as error ""
    di as error "  cd <repo root>"
    di as error "  stata -b do analyses/default/extract_patient_counts.do 5"
    di as error "{hline 90}"
    exit 601
}

* ---- named plain-text log so the run can be read back (repo convention) ----
cap mkdir "scratch"
cap log close pclog
capture log using "scratch/extract_patient_counts.log", replace text name(pclog)

* ---- create the output directory, one level at a time ----------------------
* Stata's mkdir does NOT create intermediate directories, and `capture` hides the
* failure, so a missing parent surfaces later as an opaque r(603) from postfile.
local _p "$pc_out"
local _abs = substr("`_p'", 1, 1) == "/"
local _dirs = subinstr("`_p'", "/", " ", .)
local _acc ""
foreach _part of local _dirs {
    if "`_acc'" == "" {
        if `_abs' local _acc "/`_part'"
        else      local _acc "`_part'"
    }
    else local _acc "`_acc'/`_part'"
    capture mkdir "`_acc'"
}

* ---- prove it is writable before doing any work ----------------------------
tempname _wt
capture file open `_wt' using "$pc_out/_pc_writetest.tmp", write replace
if _rc {
    di as error "{hline 90}"
    di as error "Cannot write to the output directory: $pc_out"
    di as error "Working directory is: `c(pwd)'"
    di as error "Run this do-file from the repo root, or pass an output directory"
    di as error "as the third positional argument."
    di as error "{hline 90}"
    exit 693
}
file close `_wt'
capture erase "$pc_out/_pc_writetest.tmp"

di "{hline 90}"
di "PATIENT-COUNT EXTRACTION FOR THE Ig DEMAND INTERFACE"
di "Working dir: `c(pwd)'"
di "Resamples: $pc_minbs .. $pc_maxbs   |  sims: $pc_simpath   |  years: $pc_y0 - $pc_y1"
di "Output:     $pc_out/$pc_file.csv"
di "{hline 90}"

* Variables actually needed -- keeps each load small.
local keepvars "YearMOR DateMOR"
forvalues k = 1/9 {
    local keepvars "`keepvars' YearL`k' DateL`k'S DateL`k'E"
}

**********
* Counter: posts one row per line x year for the dataset currently in memory
**********
cap program drop pc_counts
program define pc_counts
    args h b
    forvalues k = 1/9 {
        capture confirm variable YearL`k'
        if _rc continue

        capture drop _s _e_tx _e_txm _nxt
        qui gen double _s     = DateL`k'S
        qui gen double _e_tx  = DateL`k'E
        qui gen double _e_txm = DateL`k'E
        qui replace   _e_txm  = DateMOR if missing(_e_txm)

        local nxt = `k' + 1
        qui gen double _nxt = .
        capture confirm variable YearL`nxt'
        if _rc == 0 qui replace _nxt = YearL`nxt'

        forvalues y = $pc_y0/$pc_y1 {
            local jan = mdy(1,1,`y')
            local dec = mdy(12,31,`y')

            qui count if YearL`k' == `y'
            local n_start = r(N)

            * prevalent: on the treatment spell at 1 January. Spell end = DateLkE, else DateMOR,
            * else unbounded (a missing _e_txm means neither is recorded, so the spell is open).
            qui count if !missing(_s) & _s <= `jan' & (missing(_e_txm) | _e_txm >= `jan')
            local n_prev = r(N)

            qui count if !missing(_s) & _s <= `dec' & !missing(_e_tx) & _e_tx >= `jan'
            local n_tx = r(N)

            qui count if !missing(_s) & _s <= `dec' & !missing(_e_txm) & _e_txm >= `jan'
            local n_txm = r(N)

            qui count if !missing(YearL`k') & YearL`k' <= `y' ///
                       & (missing(_nxt) | _nxt > `y') ///
                       & (missing(YearMOR) | YearMOR >= `y')
            local n_st = r(N)

            post `h' (`b') (`k') (`y') (`n_start') (`n_prev') (`n_tx') (`n_txm') (`n_st')
        }
        capture drop _s _e_tx _e_txm _nxt
    }
end

**********
* 1. Point estimate (b = 0) then every bootstrap
**********
tempname pf
postfile `pf' int b byte line int year long starting long prevalent long on_tx long on_tx_mor long on_state ///
    using "$pc_out/$pc_file.dta", replace

local nfiles = 0
local nmiss  = 0

capture use `keepvars' using "$pc_point", clear
if _rc {
    capture use "$pc_point", clear
}
if _rc {
    di as error "  point-estimate file not found at $pc_point -- b = 0 skipped"
}
else {
    pc_counts `pf' 0
    di as text "  ... processed point estimate (b = 0), N = " _N
}

forvalues b = $pc_minbs/$pc_maxbs {
    local f "$pc_simpath/${pc_stub}_B`b'.dta"
    capture use `keepvars' using "`f'", clear
    if _rc {
        capture use "`f'", clear
    }
    if _rc {
        di as error "  missing bootstrap file B`b' -- skipped"
        local ++nmiss
        continue
    }
    local ++nfiles
    pc_counts `pf' `b'
    if (`b' == $pc_minbs | mod(`b', 25) == 0) di as text "  ... processed resample `b'"
}
postclose `pf'

di as text _n "Processed `nfiles' bootstrap file(s); `nmiss' missing."

**********
* 2. Export the tidy CSV -- this is the only file to pull back
**********
use "$pc_out/$pc_file.dta", clear
label variable b          "0 = point estimate, 1..B = bootstrap resample"
label variable line       "treatment line (1-9)"
label variable year       "calendar year"
label variable starting   "patients starting this line in this year"
label variable prevalent  "on the line's treatment spell at 1 January (end filled by date of death)"
label variable on_tx      "on treatment at any point in the year (diagnostic)"
label variable on_tx_mor  "on_tx, missing treatment end filled by date of death (diagnostic)"
label variable on_state   "in the line state at year granularity (diagnostic)"
order b line year starting prevalent on_tx on_tx_mor on_state
sort b line year
export delimited using "$pc_out/$pc_file.csv", replace

qui levelsof b if b > 0, local(bs)
local nb : word count `bs'
qui count
di _n "{hline 90}"
di "Wrote $pc_out/$pc_file.csv"
di "  rows: " r(N) "   |  bootstraps: `nb'   |  lines: 1-9   |  years: $pc_y0 - $pc_y1"
di "  b = 0 is the point estimate: its `starting' rows must reproduce patients_myeloma_model.csv"
di "  exactly, and its `prevalent' rows are what that file's prevalent rows are built from."
di "{hline 90}"
di "Pull back this one CSV only (a few MB). Leave the .dta on the HPC."
di "{hline 90}"

* point-estimate sanity print for lines 4 and 5, the two the Ig demand model reads
list line year starting prevalent on_tx on_tx_mor if b == 0 & inlist(line, 3, 4, 5) & year == $pc_y0, noobs abbrev(10)

cap log close pclog
