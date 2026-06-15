**********
* EpiMAP Myeloma - DVd Calibrated Transport
* compare_scenarios.do  (Tier 3: cross-scenario bootstrap aggregation)
*
* Purpose: Collate the three scenarios' 500-iteration BOOTSTRAP results into the
*          cross-scenario comparison used in the manuscript. Bootstrap-only: there
*          are no deterministic files; the point estimate is the bootstrap mean,
*          with percentile 95% CIs.
*            (1) DVd BCR distribution per scenario (mean + 95% CI)
*            (2) MAE of predicted (A, B) vs observed (C), PAIRED per iteration
*                against the same C_b, with the A-vs-B reduction (abs + %) and CI
*            (3) DVd-vs-Vd ICER / incremental cost / incremental QALY per scenario
*                (point + 95% CI), side-by-side
*
* Design: all three scenarios are bootstrapped the same way (A = resampled CASTOR
*         BCR; B = calibrated-transport bootstrap; C = resampled observed cohort),
*         so every CI reflects the same uncertainty sources and is comparable.
*         The MAE reduction pairs MAE_A_b and MAE_B_b against a COMMON C_b each
*         iteration, so C's uncertainty largely cancels in the A-vs-B difference
*         while still widening the absolute MAE and ICER CIs. A/B ("pre") and C
*         ("post") resamples are independent, aligned by index b.
*
* Run once, AFTER the 6x500 bootstrap simulations exist. Reads:
*   simulated/<scenario>/bootstrap/{dvd,vd}_2_predicted_B<b>.dta
* Writes: results/{bcr_distributions,mae_comparison,icer_comparison}.csv, results.md
*         results/bootstrap_iterations.dta  (per-iteration data, retained)
*
* NB: untested here (no Stata); sanity-check on first run with maxbs(5).
**********

clear all
set more off

* Run from the analysis folder (analyses/transport_dvd).
cap cd "~/em76/adam/analyses/transport_dvd"

**********
* Configuration
**********
global analysis "transport_dvd"
global line     2
global data     "predicted"
global maxbs    500
global sim      "simulated"
global results  "results"
cap mkdir "$results"

local bcr_labs "CR VGPR PR MR SD PD"

**********
* Helper: BCR distribution (%) + mean discounted cost & QALY from one saved dta
**********
cap program drop summ_dta
program define summ_dta, rclass
    args file
    capture use "`file'", clear
    if _rc {
        di as error "  MISSING: `file'"
        forval k = 1/6 {
            return scalar b`k' = .
        }
        return scalar cost = .
        return scalar qaly = .
        exit
    }
    local L = $line
    qui count
    local N = r(N)
    forval k = 1/6 {
        qui count if BCR_L`L' == `k'
        return scalar b`k' = 100*r(N)/`N'
    }
    qui summarize cost_total_d, meanonly
    return scalar cost = r(mean)
    qui summarize qaly_total_d, meanonly
    return scalar qaly = r(mean)
end

**********
* Helper: stash mean + 95% percentile CI of a variable into named scalars
**********
cap program drop stash
program define stash
    args v nm
    qui summarize `v', meanonly
    scalar `nm'_m = r(mean)
    qui centile `v', centile(2.5 97.5)
    scalar `nm'_lo = r(c_1)
    scalar `nm'_hi = r(c_2)
end

**********
* 1. Build per-iteration bootstrap dataset (DVd BCR + DVd-vs-Vd economics by scenario)
**********
tempname pf
postfile `pf' int b ///
    double a1 a2 a3 a4 a5 a6 aic adc adq ///
    double bb1 bb2 bb3 bb4 bb5 bb6 bic bdc bdq ///
    double c1 c2 c3 c4 c5 c6 cic cdc cdq ///
    using "$results/bootstrap_iterations.dta", replace

forval i = 1/$maxbs {
    if (`i' == 1 | mod(`i', 50) == 0) {
        di as text "  ... bootstrap iteration `i' of $maxbs"
    }
    foreach s in A_trial B_transport C_mrdr {
        local p = lower(substr("`s'",1,1))         // a / b / c
        local pp = cond("`p'"=="b","bb","`p'")     // avoid clash with postfile var "b"
        // DVd arm: BCR distribution + cost/QALY
        summ_dta "$sim/`s'/bootstrap/dvd_${line}_${data}_B`i'.dta"
        forval k = 1/6 {
            local `pp'`k' = r(b`k')
        }
        local `p'c1 = r(cost)
        local `p'q1 = r(qaly)
        // Vd arm: cost/QALY
        summ_dta "$sim/`s'/bootstrap/vd_${line}_${data}_B`i'.dta"
        local `p'c0 = r(cost)
        local `p'q0 = r(qaly)
        // incremental + ICER (DVd vs Vd)
        local `p'dc = ``p'c1' - ``p'c0'
        local `p'dq = ``p'q1' - ``p'q0'
        local `p'ic = cond(abs(``p'dq')>1e-9, ``p'dc'/``p'dq', .)
    }
    post `pf' (`i') ///
        (`a1') (`a2') (`a3') (`a4') (`a5') (`a6') (`aic') (`adc') (`adq') ///
        (`bb1') (`bb2') (`bb3') (`bb4') (`bb5') (`bb6') (`bic') (`bdc') (`bdq') ///
        (`c1') (`c2') (`c3') (`c4') (`c5') (`c6') (`cic') (`cdc') (`cdq')
}
postclose `pf'

**********
* 2. MAE vs the RAW OBSERVED DVd distribution (fixed benchmark) + reduction, CIs
*    Benchmark = the actual L2 DVd best-response marginal in the MRDR (n=533
*    evaluable; 565 starts - 32 missing BCR). Held fixed (not resampled): it is
*    the observed reality we validate against. C_mrdr is retained for the BCR
*    table and the economics, but NOT as the accuracy benchmark.
**********
use "$results/bootstrap_iterations.dta", clear

* Observed DVd BCR % (CR VGPR PR MR SD PD), n=533
local o1 = 11.44
local o2 = 23.08
local o3 = 33.40
local o4 = 7.69
local o5 = 15.38
local o6 = 9.01

gen MAE_A = (abs(a1-`o1')+abs(a2-`o2')+abs(a3-`o3')+abs(a4-`o4')+abs(a5-`o5')+abs(a6-`o6'))/6
gen MAE_B = (abs(bb1-`o1')+abs(bb2-`o2')+abs(bb3-`o3')+abs(bb4-`o4')+abs(bb5-`o5')+abs(bb6-`o6'))/6
gen red    = MAE_A - MAE_B
gen pctred = 100*red/MAE_A

stash MAE_A  maeA
stash MAE_B  maeB
stash red    red
stash pctred pctred

foreach p in a b c {
    local pp = cond("`p'"=="b","bb","`p'")
    forval k = 1/6 {
        stash `pp'`k' bcr_`p'`k'
    }
    stash `p'ic icer_`p'
    stash `p'dc dcost_`p'
    stash `p'dq dqaly_`p'
    // ICER point estimate = ratio of means (more stable than mean of per-iteration ratios)
    scalar icpt_`p' = dcost_`p'_m / dqaly_`p'_m
}

**********
* 3. Write CSVs  (point = bootstrap mean; ICER point = ratio of means)
**********
* --- BCR distributions: scenario x category, bootstrap mean + 95% CI ---
tempname fb
file open `fb' using "$results/bcr_distributions.csv", write replace
file write `fb' "scenario,bcr,mean,lo95,hi95" _n
foreach p in a b c {
    local sc = cond("`p'"=="a","A_trial",cond("`p'"=="b","B_transport","C_mrdr"))
    forval k = 1/6 {
        local lab : word `k' of `bcr_labs'
        file write `fb' "`sc',`lab',`=string(bcr_`p'`k'_m,"%6.2f")',`=string(bcr_`p'`k'_lo,"%6.2f")',`=string(bcr_`p'`k'_hi,"%6.2f")'" _n
    }
}
file close `fb'

* --- MAE comparison: bootstrap mean + 95% CI ---
tempname fm
file open `fm' using "$results/mae_comparison.csv", write replace
file write `fm' "metric,mean,lo95,hi95" _n
file write `fm' "MAE_A_traditional,`=string(maeA_m,"%6.3f")',`=string(maeA_lo,"%6.3f")',`=string(maeA_hi,"%6.3f")'" _n
file write `fm' "MAE_B_calibrated_transport,`=string(maeB_m,"%6.3f")',`=string(maeB_lo,"%6.3f")',`=string(maeB_hi,"%6.3f")'" _n
file write `fm' "MAE_reduction_abs,`=string(red_m,"%6.3f")',`=string(red_lo,"%6.3f")',`=string(red_hi,"%6.3f")'" _n
file write `fm' "MAE_reduction_pct,`=string(pctred_m,"%6.2f")',`=string(pctred_lo,"%6.2f")',`=string(pctred_hi,"%6.2f")'" _n
file close `fm'

* --- ICER comparison: scenario, ICER (ratio of means) + inc cost + inc QALY, with CIs ---
tempname fi
file open `fi' using "$results/icer_comparison.csv", write replace
file write `fi' "scenario,icer_point,icer_lo95,icer_hi95,inc_cost_mean,inc_cost_lo95,inc_cost_hi95,inc_qaly_mean,inc_qaly_lo95,inc_qaly_hi95" _n
foreach p in a b c {
    local sc = cond("`p'"=="a","A_trial",cond("`p'"=="b","B_transport","C_mrdr"))
    file write `fi' "`sc'," ///
        "`=string(icpt_`p',"%12.0f")',`=string(icer_`p'_lo,"%12.0f")',`=string(icer_`p'_hi,"%12.0f")'," ///
        "`=string(dcost_`p'_m,"%12.0f")',`=string(dcost_`p'_lo,"%12.0f")',`=string(dcost_`p'_hi,"%12.0f")'," ///
        "`=string(dqaly_`p'_m,"%7.4f")',`=string(dqaly_`p'_lo,"%7.4f")',`=string(dqaly_`p'_hi,"%7.4f")'" _n
}
file close `fi'

**********
* 4. Short results.md (headline numbers, bootstrap mean + 95% CI)
**********
tempname fr
file open `fr' using "$results/results.md", write replace
file write `fr' "# DVd Calibrated Transport - cross-scenario results" _n _n
file write `fr' "Bootstrap mean with 95% percentile CI (500 iterations; A and B bootstrapped). MAE is vs the raw observed DVd distribution (n=533, fixed); C_mrdr is retained for the BCR table and economics. Generated by compare_scenarios.do." _n _n
file write `fr' "## BCR prediction accuracy (DVd, vs observed C)" _n
file write `fr' "- MAE, Traditional (A): `=string(maeA_m,"%4.1f")' pp (95% CI `=string(maeA_lo,"%4.1f")'-`=string(maeA_hi,"%4.1f")')" _n
file write `fr' "- MAE, Calibrated Transport (B): `=string(maeB_m,"%4.1f")' pp (95% CI `=string(maeB_lo,"%4.1f")'-`=string(maeB_hi,"%4.1f")')" _n
file write `fr' "- Reduction: `=string(red_m,"%4.1f")' pp; `=string(pctred_m,"%2.0f")'% (95% CI `=string(pctred_lo,"%2.0f")'-`=string(pctred_hi,"%2.0f")'%)" _n _n
file write `fr' "## ICER, DVd vs Vd (discounted; point = ratio of mean inc cost / mean inc QALY)" _n
file write `fr' "- A Traditional: `=string(icpt_a,"%9.0fc")' (95% CI `=string(icer_a_lo,"%9.0fc")' to `=string(icer_a_hi,"%9.0fc")')" _n
file write `fr' "- B Calibrated Transport: `=string(icpt_b,"%9.0fc")' (95% CI `=string(icer_b_lo,"%9.0fc")' to `=string(icer_b_hi,"%9.0fc")')" _n
file write `fr' "- C Observed: `=string(icpt_c,"%9.0fc")' (95% CI `=string(icer_c_lo,"%9.0fc")' to `=string(icer_c_hi,"%9.0fc")')" _n _n
file write `fr' "Note: ICER percentile CIs can be unstable where the incremental QALY is near zero (notably scenario C); the inc-cost and inc-QALY CIs in icer_comparison.csv are more directly interpretable there." _n
file close `fr'

di _n "{hline 70}"
di "compare_scenarios.do complete. Outputs in: $results"
di "  bcr_distributions.csv / mae_comparison.csv / icer_comparison.csv / results.md"
di "  (per-iteration data retained in bootstrap_iterations.dta)"
di "{hline 70}"
