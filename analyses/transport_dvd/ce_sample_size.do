**********
* EpiMAP Myeloma - Cost-Effectiveness Sample-Size / Monte Carlo Precision
*                  (multi-scenario)
*
* Purpose: For EACH scenario's PSA (bootstrap) output, isolate the PARAMETER
*          uncertainty SD of the incremental QALY and express the Monte Carlo
*          (first-order) error as a fraction of it at the chosen cohort size N.
*          All scenarios share ONE cohort, so N must satisfy the most demanding
*          scenario - reported here as the "binding" scenario.
*
* Method (per scenario; O'Hagan, Stevenson & Madan 2007 ANOVA decomposition):
*   Across PSA draws b, with independent first-order draws per draw,
*       Var(dQ_b) = Var_param + sigma_pp^2 / N_boot
*   => SD_param = sqrt(Var(dQ_b) - sigma_pp^2 / N_boot)
*      MC SD at N = sigma_pp/sqrt(N);  ratio = MC SD / SD_param;  pick N so ratio <= kfrac.
*
* sigma_pp per scenario comes from each scenario's ce_precision run (Config block).
*
* Reads:  simulated/<scenario>/bootstrap/{dvd,vd}_<line>_<data>_B<b>.dta
* Output: results/ce_sample_size.csv        (one row per scenario + binding flag)
*         results/ce_sample_size_ratio.png  (MC/parameter ratio vs N, all scenarios)
*
* Reusable: edit the Configuration block (scenarios, sigma_pp, N, kfrac).
*
* Author: Adam Irving
* Date: June 2026
**********

clear all
set more off
cap cd "~/em76_scratch2/adam/transport_dvd"   // where the bootstrap output lives

**********
* Configuration
**********
local scenarios "A_trial B_transport C_mrdr"
local line      2
local data      "predicted"
local maxbs     500
local bootN     50000          // patients per PSA iteration (the N the bootstrap was run at)
local kfrac     0.05           // MC target: MC SD <= kfrac * parameter SD
local chosen_n  50000          // cohort size to evaluate / defend

local sim     "simulated"
local results "results"
cap mkdir "`results'"

* sigma_pp per scenario, from each scenario's ce_precision run. ce_precision is
* run on a DIFFERENT machine (Mac) to this one (HPC), so supply it either way:
*   1. paste the values below (from ce_precision's console output / ce_precision_sigma.csv), OR
*   2. copy ce_precision_sigma.csv into results/ - if present it OVERRIDES the values.
local spp_dQ_A_trial      0.5137
local spp_dC_A_trial      110907
local spp_dQ_B_transport  0.5137
local spp_dC_B_transport  110907
local spp_dQ_C_mrdr       0.5137
local spp_dC_C_mrdr       110907

local sigfile "`results'/ce_precision_sigma.csv"
local have_sig = 0
capture confirm file "`sigfile'"
if (_rc == 0) {
    local have_sig = 1
    tempfile sigtab
    import delimited "`sigfile'", clear case(lower) varnames(1)
    keep scenario sd_pp_dq sd_pp_dc
    qui save "`sigtab'"
    di as text "sigma_pp: using `sigfile'"
}
else {
    di as text "sigma_pp: table not found - using the config values"
}

**********
* 1. Per-scenario decomposition
**********
tempname S
postfile `S' str16 scenario double dQ dQlo dQhi double sd_param_dQ sigpp_dQ ///
    double mcsd_dQ mcfrac_dQ reqN_dQ sd_param_dC int iters ///
    using "`results'/_ss_byscenario.dta", replace

foreach s of local scenarios {

    if (`have_sig') {
        use "`sigtab'", clear
        qui summarize sd_pp_dq if scenario == "`s'", meanonly
        local spp_dQ = r(mean)
        qui summarize sd_pp_dc if scenario == "`s'", meanonly
        local spp_dC = r(mean)
    }
    else {
        local spp_dQ "`spp_dQ_`s''"
        local spp_dC "`spp_dC_`s''"
    }
    if ("`spp_dQ'" == "") {
        di as error "  `s': no sigma_pp (table or config) - skipped"
        continue
    }
    if (`spp_dQ' == . | `spp_dQ' == 0) {
        di as error "  `s': sigma_pp missing/zero - skipped"
        continue
    }

    * collect per-iteration increments for this scenario
    tempfile iters_s
    tempname P
    postfile `P' int b double dQ dC using "`iters_s'", replace
    local nfound = 0
    forval i = 1/`maxbs' {
        local fdvd "`sim'/`s'/bootstrap/dvd_`line'_`data'_B`i'.dta"
        local fvd  "`sim'/`s'/bootstrap/vd_`line'_`data'_B`i'.dta"
        capture confirm file "`fdvd'"
        local ok = (_rc == 0)
        capture confirm file "`fvd'"
        local ok = `ok' & (_rc == 0)
        if (!`ok') continue
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
        post `P' (`i') (`q1' - `q0') (`c1' - `c0')
        local ++nfound
    }
    postclose `P'
    if (`nfound' == 0) {
        di as error "  `s': no bootstrap files found under `sim'/`s'/bootstrap/ - skipped"
        continue
    }

    use "`iters_s'", clear
    qui summarize dQ
    local mQ    = r(mean)
    local vtotQ = r(Var)
    qui centile dQ, centile(2.5 97.5)
    local loQ = r(c_1)
    local hiQ = r(c_2)
    local sdparQ  = sqrt(max(0, `vtotQ' - (`spp_dQ'^2)/`bootN'))
    local mcsdQ   = `spp_dQ'/sqrt(`chosen_n')
    local mcfracQ = `mcsdQ'/`sdparQ'
    local reqNQ   = (`spp_dQ'/(`kfrac'*`sdparQ'))^2

    qui summarize dC
    local sdparC = sqrt(max(0, r(Var) - (`spp_dC'^2)/`bootN'))

    post `S' ("`s'") (`mQ') (`loQ') (`hiQ') (`sdparQ') (`spp_dQ') ///
        (`mcsdQ') (`mcfracQ') (`reqNQ') (`sdparC') (`nfound')

    di as text "  `s': dQ=" as result %7.4f `mQ' as text "  SD_param=" %7.4f `sdparQ' ///
        as text "  MC/param@`chosen_n'=" %4.1f (100*`mcfracQ') "%" ///
        as text "  reqN(`=100*`kfrac''%)=" %9.0fc `reqNQ'
}
postclose `S'

**********
* 2. Summary table + binding scenario
**********
use "`results'/_ss_byscenario.dta", clear
count
if (r(N) == 0) {
    di as error "No scenarios processed - check the cd path / bootstrap folders."
    exit 601
}
qui summarize reqN_dQ
local reqN_binding = r(max)
gen byte binding = (reqN_dQ == `reqN_binding')
gsort -reqN_dQ

di _n "{hline 78}"
di as text "Sample-size by scenario (PSA N=`bootN', chosen N=`chosen_n', target k=`=100*`kfrac''%)"
di "{hline 78}"
list scenario dQ sd_param_dQ mcfrac_dQ reqN_dQ binding, clean noobs
export delimited using "`results'/ce_sample_size.csv", replace

local adeq = cond(`chosen_n' >= `reqN_binding', "ADEQUATE for all scenarios", "NOT adequate - raise N")
di _n as text "Binding scenario requires N >= " as result %9.0fc `reqN_binding' ///
    as text "  ->  chosen N=`chosen_n' is " as result "`adeq'"
di "{hline 78}"

* stash per-scenario sigma_pp and parameter SD for the figure
foreach s of local scenarios {
    qui summarize sigpp_dQ if scenario == "`s'", meanonly
    local sig_`s' = r(mean)
    qui summarize sd_param_dQ if scenario == "`s'", meanonly
    local sdp_`s' = r(mean)
}

**********
* 3. Figure: MC SD as a fraction of parameter SD vs N, one line per scenario
*    (each scenario meets the target where its curve drops below kfrac)
**********
clear
local npts = 300
set obs `npts'
gen double N = round(exp(ln(1000) + (_n-1)/(`npts'-1)*(ln(300000)-ln(1000))))

local plot ""
local k = 0
foreach s of local scenarios {
    if ("`sig_`s''" == "" | "`sdp_`s''" == ".") continue
    local ++k
    gen double ratio_`k' = (`sig_`s''/sqrt(N)) / `sdp_`s''
    label variable ratio_`k' "`s'"
    local plot "`plot' (line ratio_`k' N)"
}

twoway `plot' ///
    , xline(`chosen_n', lcolor(gs10) lpattern(dot)) ///
    xscale(log) ///
    xlabel(1000 5000 10000 20000 50000 100000 300000, angle(45) labsize(small)) ///
    ylabel(0(0.05)0.3, angle(0) format(%4.2f) labsize(small)) ///
    xtitle("Simulated cohort size, N") ///
    ytitle("Monte Carlo SD as a fraction of parameter SD") ///
    text(0.28 `chosen_n' "N = `chosen_n'", place(e) size(small)) ///
    title("MC error vs parameter uncertainty, by scenario", size(medium)) ///
    legend(pos(2) ring(0) cols(1) size(small) region(lstyle(none))) ///
    graphregion(color(white)) plotregion(margin(medium))

graph export "`results'/ce_sample_size_ratio.png", replace width(2400)
cap graph export "`results'/ce_sample_size_ratio.pdf", replace
di as text "Saved -> results/ce_sample_size.csv and results/ce_sample_size_ratio.png"
