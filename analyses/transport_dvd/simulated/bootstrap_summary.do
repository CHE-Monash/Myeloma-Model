**********
* Monash Myeloma Model - DVd Calibrated Transport - OOS bootstrap aggregation
*
* Purpose: Collate the three scenarios' 500-iteration bootstrap results into the cross-scenario
*          comparison used in the manuscript (Tier 3). Bootstrap-only: no deterministic files, the
*          point estimate is the bootstrap mean with percentile 95% CIs. Produces (1) DVd AND Vd BCR
*          distribution per scenario (bcr_distributions.csv, with an `intervention` column dvd/vd);
*          (2) MAE of predicted (A, B) vs the RAW OBSERVED DVd distribution (n=399, held fixed;
*          see section 2), with the A-vs-B reduction (abs + %), CI and exceedance; (3) DVd-vs-Vd
*          ICER / incremental cost / QALY
*          per scenario (point + 95% CI), side-by-side.
* Notes:   Run once, after the 6x500 bootstrap sims exist. Reads
*          simulated/<scenario>/bootstrap/{dvd,vd}_2_predicted_B<b>.dta; writes
*          results/{bcr_distributions,mae_comparison,icer_comparison}.csv, results.md and
*          results/bootstrap_iterations.dta. All scenarios are bootstrapped the same way (A =
*          resampled CASTOR BCR; B = calibrated-transport bootstrap; C = fitted registry-model reference),
*          so CIs are comparable. MAE_A_b and MAE_B_b are scored against the same fixed observed
*          benchmark and the reduction is paired within replicate; the MAE intervals and exceedance
*          are therefore conditional on the observed distribution and carry no sampling uncertainty
*          in it (stated in the manuscript's Methods). C_mrdr, the fitted registry-model economic reference, feeds
*          the BCR table and the economics, not the benchmark. (Header corrected 12 Oct 2026; it had
*          described an earlier paired-C_b design with an n=533 benchmark.)
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

* Run from the analysis folder (analyses/transport_dvd).

**********
* Configuration
**********
global analysis "transport_dvd"
global line     2
global data     "predicted"
// Gaps are filled in place by hpc/fill_gaps.sh, which renames spare replicates generated above 500
// down into the failed indices - MI file, coefficients and all three transport artefacts together,
// then verifies 1..500 is complete. So this stays 500 and the set is contiguous.
global maxbs    500
* Resolve paths whether this is run from the REPO ROOT (as run.do does) or from the
* analysis folder. Previously these were bare relative paths, so running from the root
* silently created an empty <repo>/results/ and then found no inputs: the same class of
* failure already documented in outcomes/calibrated_transport.do.
if (fileexists("analyses/$analysis/simulated/bootstrap_summary.do")) {
    global base "analyses/$analysis"      // repo root
}
else {
    global base "."                        // already in analyses/transport_dvd
}
global sim      "$base/simulated"
global results  "$base/results"
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
*
*    Stage 1 reads the 3,000 per-replicate .dta files under each scenario's simulated bootstrap folder, which
*    live on the HPC and are not pulled back. Stage 2 needs only bootstrap_iterations.dta
*    and the benchmark hard-coded below, and the per-replicate distributions do not depend
*    on that benchmark. So if the benchmark changes, the summary can be refreshed LOCALLY
*    without re-bootstrapping:  do bootstrap_summary.do 0
**********
*    Decide automatically: rebuild when the per-replicate files are actually present (the
*    HPC), otherwise reuse bootstrap_iterations.dta (the usual local case). Pass an argument
*    to override -  do bootstrap_summary.do 0  - rather than setting a global: `clear all'
*    does NOT drop globals, so a $rebuild left over from an earlier run in the same session
*    used to survive into the next one and abort it at the probe with nothing written.
args rebuild_arg
local probe "$sim/A_trial/bootstrap/dvd_${line}_${data}_B1.dta"
if "`rebuild_arg'" != "" {
    global rebuild `rebuild_arg'
}
else {
    if (fileexists("`probe'"))                              global rebuild 1
    else if (fileexists("$results/bootstrap_iterations.dta")) global rebuild 0
    else                                                     global rebuild 1
}
if $rebuild == 0 {
    di as text "Stage 1 skipped: reusing $results/bootstrap_iterations.dta"
    confirm file "$results/bootstrap_iterations.dta"
}
else if (!fileexists("`probe'")) {
    di as error "Per-replicate bootstrap files not found, e.g. `probe'"
    di as error "They are produced on the HPC and are not pulled back. To refresh the"
    di as error "summaries from the existing iteration dataset, run: do bootstrap_summary.do 0"
    error 601
}
if $rebuild {
tempname pf
postfile `pf' int b ///
    double a1 a2 a3 a4 a5 a6 av1 av2 av3 av4 av5 av6 aic adc adq ///
    double bb1 bb2 bb3 bb4 bb5 bb6 bbv1 bbv2 bbv3 bbv4 bbv5 bbv6 bic bdc bdq ///
    double c1 c2 c3 c4 c5 c6 cv1 cv2 cv3 cv4 cv5 cv6 cic cdc cdq ///
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
        // Vd arm: BCR distribution + cost/QALY
        summ_dta "$sim/`s'/bootstrap/vd_${line}_${data}_B`i'.dta"
        forval k = 1/6 {
            local `pp'v`k' = r(b`k')
        }
        local `p'c0 = r(cost)
        local `p'q0 = r(qaly)
        // incremental + ICER (DVd vs Vd)
        local `p'dc = ``p'c1' - ``p'c0'
        local `p'dq = ``p'q1' - ``p'q0'
        local `p'ic = cond(abs(``p'dq')>1e-9, ``p'dc'/``p'dq', .)
    }
    post `pf' (`i') ///
        (`a1') (`a2') (`a3') (`a4') (`a5') (`a6') (`av1') (`av2') (`av3') (`av4') (`av5') (`av6') (`aic') (`adc') (`adq') ///
        (`bb1') (`bb2') (`bb3') (`bb4') (`bb5') (`bb6') (`bbv1') (`bbv2') (`bbv3') (`bbv4') (`bbv5') (`bbv6') (`bic') (`bdc') (`bdq') ///
        (`c1') (`c2') (`c3') (`c4') (`c5') (`c6') (`cv1') (`cv2') (`cv3') (`cv4') (`cv5') (`cv6') (`cic') (`cdc') (`cdq')
}
postclose `pf'
}   // end stage 1 (skipped when $rebuild == 0)

**********
* 2. MAE vs the RAW OBSERVED DVd distribution (fixed benchmark) + reduction, CIs
*    Benchmark = the second-line DVd best-response marginal in the MRDR (held-out
*    cohort: Event0==20, Regimen==80, BCR-evaluable; same cohort as
*    mrdr_dvd_baseline.csv). HARD-CODED because this script runs on the HPC, which
*    cannot reach the raw registry. The values below are produced by the LOCAL run
*    of calibrated_transport.do, which writes results/observed_dvd_bcr.csv and prints
*    these `local' lines ready to paste here -- update them if the held-out cohort
*    changes. Held fixed (not resampled); C is the fitted registry-model economic reference; accuracy uses the separate fixed recorded DVd benchmark. C_mrdr is
*    retained for the BCR table/economics but NOT this benchmark.
**********
use "$results/bootstrap_iterations.dta", clear

* Observed DVd BCR % (CR VGPR PR MR SD PD) + n  -- from results/observed_dvd_bcr.csv
* 10 Aug 2026: was n=529 (11.34/23.06/33.46/7.94/15.12/9.07). That cohort predated the 13 Jul
* nomaster fix, when the m=0 master was BCR carry-forward filled and 130 patients with no recorded
* L2 response passed the drop-if-missing filter carrying their L1 response. See
* scratch/bcr/benchmark_n_check.do, which reproduces the old figures exactly by re-applying the LOCF.
local oN = 399
local o1 = 8.77
local o2 = 22.81
local o3 = 34.84
local o4 = 9.27
local o5 = 13.53
local o6 = 10.78

gen MAE_A = (abs(a1-`o1')+abs(a2-`o2')+abs(a3-`o3')+abs(a4-`o4')+abs(a5-`o5')+abs(a6-`o6'))/6
gen MAE_B = (abs(bb1-`o1')+abs(bb2-`o2')+abs(bb3-`o3')+abs(bb4-`o4')+abs(bb5-`o5')+abs(bb6-`o6'))/6
gen red    = MAE_A - MAE_B
gen pctred = 100*red/MAE_A

stash MAE_A  maeA
stash MAE_B  maeB
stash red    red
stash pctred pctred

* Exceedance probability: the share of replicates in which Calibrated Transport is
* closer to the observed benchmark than the trial-based prediction. This is the paired
* summary quoted in the manuscript, and it was previously computed ad hoc and never
* written out, so the reported figure had no results file behind it.
qui count if MAE_B < MAE_A
scalar exceed_pct = 100*r(N)/_N

foreach p in a b c {
    local pp = cond("`p'"=="b","bb","`p'")
    forval k = 1/6 {
        stash `pp'`k'  bcr_`p'`k'
        stash `pp'v`k' bcr_`p'v`k'
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
file write `fb' "scenario,intervention,bcr,mean,lo95,hi95" _n
foreach p in a b c {
    local sc = cond("`p'"=="a","A_trial",cond("`p'"=="b","B_transport","C_mrdr"))
    forval k = 1/6 {
        local lab : word `k' of `bcr_labs'
        file write `fb' "`sc',dvd,`lab',`=string(bcr_`p'`k'_m,"%6.2f")',`=string(bcr_`p'`k'_lo,"%6.2f")',`=string(bcr_`p'`k'_hi,"%6.2f")'" _n
        file write `fb' "`sc',vd,`lab',`=string(bcr_`p'v`k'_m,"%6.2f")',`=string(bcr_`p'v`k'_lo,"%6.2f")',`=string(bcr_`p'v`k'_hi,"%6.2f")'" _n
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
* paired summary: no interval, it is a proportion over the replicates, not a resampled statistic
file write `fm' "exceedance_pct,`=string(exceed_pct,"%6.2f")',," _n
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
file write `fr' "Bootstrap mean with 95% percentile CI (500 iterations; A and B bootstrapped). MAE is vs the raw observed DVd distribution (n=`oN', fixed); C_mrdr is retained for the BCR table and economics. Generated by bootstrap_summary.do." _n _n
file write `fr' "## BCR prediction accuracy (DVd, vs recorded validation benchmark)" _n
file write `fr' "- Observed benchmark (n=`oN'): CR `=string(`o1',"%4.1f")', VGPR `=string(`o2',"%4.1f")', PR `=string(`o3',"%4.1f")', MR `=string(`o4',"%4.1f")', SD `=string(`o5',"%4.1f")', PD `=string(`o6',"%4.1f")'" _n
file write `fr' "- MAE, Trial-based (A): `=string(maeA_m,"%4.1f")' pp (95% CI `=string(maeA_lo,"%4.1f")'-`=string(maeA_hi,"%4.1f")')" _n
file write `fr' "- MAE, Calibrated Transport (B): `=string(maeB_m,"%4.1f")' pp (95% CI `=string(maeB_lo,"%4.1f")'-`=string(maeB_hi,"%4.1f")')" _n
file write `fr' "- Reduction: `=string(red_m,"%4.1f")' pp; `=string(pctred_m,"%2.0f")'% (95% CI `=string(pctred_lo,"%2.0f")'-`=string(pctred_hi,"%2.0f")'%)" _n
file write `fr' "- Exceedance: Calibrated Transport closer to observed in `=string(exceed_pct,"%4.1f")'% of replicates" _n _n
file write `fr' "## ICER, DVd vs Vd (discounted; point = ratio of mean inc cost / mean inc QALY)" _n
file write `fr' "- A Trial-based: `=string(icpt_a,"%9.0fc")' (95% CI `=string(icer_a_lo,"%9.0fc")' to `=string(icer_a_hi,"%9.0fc")')" _n
file write `fr' "- B Calibrated Transport: `=string(icpt_b,"%9.0fc")' (95% CI `=string(icer_b_lo,"%9.0fc")' to `=string(icer_b_hi,"%9.0fc")')" _n
file write `fr' "- C Registry-based: `=string(icpt_c,"%9.0fc")' (95% CI `=string(icer_c_lo,"%9.0fc")' to `=string(icer_c_hi,"%9.0fc")')" _n _n
file write `fr' "Note: ICER percentile CIs can be unstable where the incremental QALY is near zero (notably scenario C); the inc-cost and inc-QALY CIs in icer_comparison.csv are more directly interpretable there." _n
file close `fr'

di _n "{hline 70}"
di "bootstrap_summary.do complete. Outputs in: $results"
di "  bcr_distributions.csv / mae_comparison.csv / icer_comparison.csv / results.md"
di "  (per-iteration data retained in bootstrap_iterations.dta)"
di "{hline 70}"
