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

* Confidence-interval BOUNDS as lines, no shading: observed (validation) bounds solid black, simulated
* bounds dotted -- so the overlap is legible. The monthly p-value and the p=0.05 threshold are light
* blue (2nd axis) to set them apart from the black KM CI bounds.
gen double year = month/12
local pcol "74 143 199"     // light blue for the p-value lines
* Bounds as lines only -> the simulated dotted line is one object shared by plot + legend, so they
* match exactly (no scatter/proxy). Validation solid black; simulated dashed.
twoway ///
    (line obs_lo year, lcolor(black) lpattern(solid) lwidth(medthin)) ///
    (line obs_hi year, lcolor(black) lpattern(solid) lwidth(medthin)) ///
    (line sim_lo year, lcolor(black) lpattern(shortdash) lwidth(medthin)) ///
    (line sim_hi year, lcolor(black) lpattern(shortdash) lwidth(medthin)) ///
    (line pvalue year, yaxis(2) lcolor("`pcol'") lwidth(vthin)  lpattern(shortdash)) ///
    (function y = 0.05, range(0 10) yaxis(2) lcolor("`pcol'") lwidth(vthin) lpattern(solid)) ///
  , ytitle("Overall survival (%)") yscale(range(0 100)) ylabel(0(20)100, angle(0)) ///
    ytitle("p-value", axis(2) color("`pcol'")) yscale(range(0 1) axis(2) lcolor("`pcol'") lwidth(vthin)) ///
    ylabel(0(.2)1, axis(2) angle(0) format(%2.1f) labcolor("`pcol'") tlcolor("`pcol'") tlwidth(vthin)) ///
    xtitle("Years since diagnosis", margin(t=3)) xlabel(0(2)10) ///
    legend(order(1 "Validation 95% CI" 3 "Simulated 95% CI" 5 "Monthly p-value" 6 "p = 0.05") ///
           rows(1) size(vsmall) symxsize(8) region(lstyle(none))) ///
    graphregion(color(white)) plotregion(color(white)) ///
    name(os_curve, replace)

graph export "analyses/oos/results/os_wholepop_curve.png", replace width(1800)
graph save   "analyses/oos/results/os_wholepop_curve.gph", replace
di as text "Wrote analyses/oos/results/os_wholepop_curve.png (+ .gph)"
