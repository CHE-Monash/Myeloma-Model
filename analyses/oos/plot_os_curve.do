**********
* Monash Myeloma Model - OOS: redraw the OS validation figure (2024 PLOS ONE Fig 2 style) locally
*
* Purpose: Draw the whole-population OS validation figure -- observed (validation-cohort) KM 95% CI vs
*          the simulated cohort's 95% CI, plus the monthly p-value testing the difference -- from
*          os_wholepop_curve_validation.csv, which bootstrap_validation.do writes on the HPC. Use this
*          when HPC batch graphics are disabled, or to restyle the figure without re-running the
*          validation. Reads only the pulled-back CSV (no drive, no bootstrap sims).
*
* Run from the repository root:  do "analyses/oos/plot_os_curve.do"
*   optional arg 1 = path to the validation CSV
*                    (default: analyses/oos/results/os_wholepop_curve_validation.csv)
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

local csv `"`1'"'
if `"`csv'"' == "" local csv "analyses/oos/results/os_wholepop_curve_validation.csv"
capture confirm file "`csv'"
if _rc {
    di as error "Validation CSV not found: `csv'"
    di as error "Run analyses/oos/bootstrap_validation.do first (it writes os_wholepop_curve_validation.csv), then pull it back."
    exit 601
}

import delimited "`csv'", clear case(preserve)

* Headline (as in the 2024 paper): fraction of months 1..120 with no significant difference (p >= 0.05).
qui count if month >= 1 & month <= 120
local nmo = r(N)
qui count if month >= 1 & month <= 120 & pvalue >= 0.05
local nns = r(N)
di as text _n "OS validation: no significant difference in `nns'/`nmo' months (" ///
    %4.1f cond(`nmo' > 0, 100*`nns'/`nmo', .) "%); p>=0.05" _n

twoway ///
    (rarea obs_lo obs_hi month, color(navy%22) lwidth(none)) ///
    (rarea sim_lo sim_hi month, color(cranberry%22) lwidth(none)) ///
    (line  s_obs    month, lcolor(navy)      lwidth(medthin)) ///
    (line  sim_mean month, lcolor(cranberry) lwidth(medthin)) ///
    (line  pvalue   month, yaxis(2) lcolor(gs7) lpattern(dash)) ///
  , yline(0.05, axis(2) lcolor(gs11) lpattern(shortdash)) ///
    ytitle("Overall survival (%)") yscale(range(0 100)) ylabel(0(20)100, angle(0)) ///
    ytitle("Monthly p-value", axis(2)) yscale(range(0 1) axis(2)) ylabel(0(.2)1, axis(2) angle(0)) ///
    xtitle("Months since diagnosis") xlabel(0(12)120) ///
    legend(order(1 "Observed (validation) 95% CI" 2 "Simulated 95% CI" 5 "Monthly p-value") ///
           rows(1) size(vsmall) region(lstyle(none))) ///
    title("Out-of-sample overall survival: observed vs simulated", size(medsmall)) ///
    name(os_curve, replace)

graph export "analyses/oos/results/os_wholepop_curve.png", replace width(1800)
graph save   "analyses/oos/results/os_wholepop_curve.gph", replace
di as text "Wrote analyses/oos/results/os_wholepop_curve.png (+ .gph)"
