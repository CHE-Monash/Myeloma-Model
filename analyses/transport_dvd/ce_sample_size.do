**********
* EpiMAP Myeloma - Cost-Effectiveness Sample-Size / Monte Carlo Precision
*
* Purpose: From a scenario's PSA (bootstrap) output, isolate the PARAMETER
*          uncertainty SD of the incremental QALY/cost and find the simulated
*          cohort size N at which Monte Carlo (first-order) error is negligible
*          relative to it. Produces the supplementary convergence figure that
*          defends the chosen N.
*
* Reads (as written by the bootstrap dispatcher; same paths as bootstrap_summary.do):
*   simulated/<scenario>/bootstrap/{dvd,vd}_<line>_<data>_B<b>.dta
*   each per-patient with cost_total_d, qaly_total_d.
*
* Needs sigma_pp (per-patient SD of the PAIRED increment) from ce_precision.do,
*   measured on the SAME CRN engine. It is N-independent, so one value serves
*   every N (config below).
*
* Method (O'Hagan, Stevenson & Madan 2007 ANOVA decomposition):
*   Across PSA draws b, with independent first-order draws per draw,
*       Var(dQ_b) = Var_param + sigma_pp^2 / N_boot
*   => Var_param = Var(dQ_b) - sigma_pp^2 / N_boot      (the MC component removed)
*   MC SD at any N is sigma_pp/sqrt(N); choose N so it is <= k * SD_param.
*   At N_boot = 50k the subtracted term is tiny, so the corrected and raw PSA SDs
*   nearly coincide - the correction matters only at small N.
*
* Reusable: edit the Configuration block for any analysis/scenario with a bootstrap.
*
* Output: results/ce_sample_size.csv   (decomposition + required N)
*         results/ce_sample_size.png/.pdf
*
* Author: Adam Irving
* Date: June 2026
**********

clear all
set more off

**********
* Configuration
**********
cap cd "~/em76/adam/analyses/transport_dvd"   // <-- where the bootstrap output lives (matches bootstrap_summary.do)

local scenario  "B_transport"
local line      2
local data      "predicted"
local maxbs     500
local bootN     50000          // patients per PSA iteration (the N the bootstrap was run at)

* sigma_pp: per-patient SD of the paired increment, from ce_precision.csv (CRN engine)
local spp_dQ    0.5538         // sd_pp_dQ   (UPDATE from results/ce_precision.csv)
local spp_dC    119585         // sd_pp_dC

local wtp       50000          // WTP (A$/QALY) for INMB; AUS has no single value - vary as needed
local kfrac     0.10           // MC target: MC SD <= kfrac * parameter SD
local chosen_n  50000          // cohort size to mark/defend on the figure

local sim       "simulated"
local results   "results"
cap mkdir "`results'"

**********
* 1. Collect per-iteration incremental means: dQ_b, dC_b, dNMB_b
**********
tempname pf
postfile `pf' int b double dQ dC dNMB ///
    using "`results'/_ss_iters.dta", replace

local nfound = 0
forval i = 1/`maxbs' {
    local fdvd "`sim'/`scenario'/bootstrap/dvd_`line'_`data'_B`i'.dta"
    local fvd  "`sim'/`scenario'/bootstrap/vd_`line'_`data'_B`i'.dta"
    capture confirm file "`fdvd'"
    local ok = (_rc == 0)
    capture confirm file "`fvd'"
    local ok = `ok' & (_rc == 0)
    if (!`ok') {
        di as error "  missing pair at b=`i' - skipped"
        continue
    }
    qui use "`fdvd'", clear
    qui summarize qaly_total_d, meanonly
    local q1 = r(mean)
    qui summarize cost_total_d, meanonly
    local c1 = r(mean)
    qui use "`fvd'", clear
    qui summarize qaly_total_d, meanonly
    local q0 = r(mean)
    qui summarize cost_total_d, meanonly
    local c0 = r(mean)

    local dq   = `q1' - `q0'
    local dc   = `c1' - `c0'
    local dnmb = `wtp'*`dq' - `dc'
    post `pf' (`i') (`dq') (`dc') (`dnmb')
    local ++nfound
    if (`i'==1 | mod(`i',50)==0) di as text "  iteration `i' of `maxbs'"
}
postclose `pf'
di as text "Iterations collected: `nfound'"
if (`nfound' == 0) {
    di as error "No bootstrap files found - check the cd path and simulated/<scenario>/bootstrap/ names."
    exit 601
}

**********
* 2. Variance decomposition: total (PSA) = parameter + Monte Carlo
**********
use "`results'/_ss_iters.dta", clear

* --- Incremental QALY ---
qui summarize dQ
scalar mQ    = r(mean)
scalar vtotQ = r(Var)
qui centile dQ, centile(2.5 97.5)
scalar loQ = r(c_1)
scalar hiQ = r(c_2)
scalar vmcQ   = (`spp_dQ'^2)/`bootN'
scalar vparQ  = max(0, vtotQ - vmcQ)
scalar sdtotQ = sqrt(vtotQ)
scalar sdparQ = sqrt(vparQ)
scalar mcsdQ_chosen = `spp_dQ'/sqrt(`chosen_n')
scalar mcfracQ = mcsdQ_chosen / sdparQ
scalar reqNQ   = (`spp_dQ'/(`kfrac'*sdparQ))^2

* --- Incremental cost ---
qui summarize dC
scalar mC     = r(mean)
scalar vtotC  = r(Var)
scalar vparC  = max(0, vtotC - (`spp_dC'^2)/`bootN')
scalar sdparC = sqrt(vparC)
scalar mcfracC = (`spp_dC'/sqrt(`chosen_n')) / sdparC

* --- INMB (report PSA SD; MC negligible at this N) ---
qui summarize dNMB
scalar mNMB     = r(mean)
scalar sdtotNMB = r(sd)
qui centile dNMB, centile(2.5 97.5)
scalar loNMB = r(c_1)
scalar hiNMB = r(c_2)

**********
* 3. Console summary
**********
di _n "{hline 74}"
di as text "Cost-effectiveness sample-size analysis - `scenario' (PSA N=`bootN', `nfound' iters)"
di "{hline 74}"
di as text "Incremental QALY (dQ):"
di as text "   PSA mean dQ              = " as result %8.4f mQ as text "  (95% CI " %7.4f loQ " to " %7.4f hiQ ")"
di as text "   Total PSA SD            = " as result %8.4f sdtotQ
di as text "   - MC component (N=`bootN')   = " as result %8.5f mcsdQ_chosen
di as text "   = Parameter SD           = " as result %8.4f sdparQ
di as text "   MC / param at N=`chosen_n'   = " as result %5.1f (100*mcfracQ) as text " %"
di as text "   N s.t. MC <= " as result %2.0f (100*`kfrac') as text "% of param SD = " as result %9.0fc reqNQ
di _n as text "Incremental cost (dC):  parameter SD = " as result %9.0fc sdparC ///
    as text "   (MC/param at N=`chosen_n' = " as result %4.1f (100*mcfracC) as text "%)"
di as text "INMB @ " as result %6.0fc `wtp' as text "/QALY: mean = " as result %9.0fc mNMB ///
    as text "   PSA SD = " as result %9.0fc sdtotNMB
di "{hline 74}"

**********
* 4. Summary CSV
**********
clear
set obs 1
gen scenario          = "`scenario'"
gen psa_N             = `bootN'
gen iters             = `nfound'
gen wtp               = `wtp'
gen k_target          = `kfrac'
gen chosen_N          = `chosen_n'
gen dQ_mean           = mQ
gen dQ_lo95           = loQ
gen dQ_hi95           = hiQ
gen sd_total_dQ       = sdtotQ
gen mc_sd_dQ_chosenN  = mcsdQ_chosen
gen sd_param_dQ       = sdparQ
gen sigma_pp_dQ       = `spp_dQ'
gen mcfrac_dQ_chosenN = mcfracQ
gen reqN_dQ           = reqNQ
gen sd_param_dC       = sdparC
gen NMB_mean          = mNMB
gen NMB_sd            = sdtotNMB
export delimited using "`results'/ce_sample_size.csv", replace
di as text "Saved `results'/ce_sample_size.csv"

**********
* 5. Convergence figure: MC SD of dQ vs N, against the parameter-SD reference
**********
clear
local npts = 300
set obs `npts'
local nmin = 1000
local nmax = 300000
gen double N = round(exp(ln(`nmin') + (_n-1)/(`npts'-1)*(ln(`nmax')-ln(`nmin'))))
gen double mcsd = `spp_dQ'/sqrt(N)

local par   = sdparQ
local tgt   = `kfrac'*sdparQ
local ychos = `spp_dQ'/sqrt(`chosen_n')

twoway ///
    (line mcsd N, lcolor(navy) lwidth(medthick)) ///
    (scatteri `ychos' `chosen_n', msymbol(O) mcolor(black) msize(medium)) ///
    , ///
    yline(`par', lcolor(red) lpattern(dash)) ///
    yline(`tgt', lcolor(orange) lpattern(shortdash)) ///
    xline(`chosen_n', lcolor(gs11) lpattern(dot)) ///
    xscale(log) ///
    xlabel(1000 2000 5000 10000 20000 50000 100000 300000, angle(45) labsize(small)) ///
    ylabel(, angle(0) format(%5.3f) labsize(small)) ///
    xtitle("Simulated cohort size, N") ///
    ytitle("Monte Carlo SD of incremental QALY") ///
    text(`par' `nmin' "Parameter SD = `=string(`par',"%5.3f")'", place(ne) size(small) color(red)) ///
    text(`tgt' `nmin' "`=100*`kfrac''% of parameter SD", place(se) size(small) color(orange)) ///
    text(`ychos' `chosen_n' "  N = `chosen_n'", place(e) size(small)) ///
    legend(off) ///
    title("Monte Carlo error vs parameter uncertainty (`scenario')", size(medium)) ///
    note("MC SD = sigma_pp/sqrt(N), sigma_pp = `=string(`spp_dQ',"%4.2f")'. Required N for MC <= `=100*`kfrac''% of parameter SD: `=string(reqNQ,"%9.0fc")'.") ///
    graphregion(color(white)) plotregion(margin(medium))

graph export "`results'/ce_sample_size.png", replace width(2200)
cap graph export "`results'/ce_sample_size.pdf", replace
di as text "Saved `results'/ce_sample_size.png (.pdf)"
