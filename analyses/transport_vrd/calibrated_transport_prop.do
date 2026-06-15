**********
* EpiMAP Myeloma - VRd L1 Calibrated Transport  (PARTIAL PROPORTIONAL ODDS variant)
*
* Companion to calibrated_transport.do. Identical data and design (real MRDR records
* for the registry cohorts; synthetic SWOG arms); the only change is the
* Calibrated-Transport model. Instead of a single proportional-odds shift
*     ologit BCR MRDR VRd
* this fits a PARTIAL proportional-odds (generalised ordered logit) model
*     gologit2 BCR MRDR VRd, npl(MRDR)
* in which the cross-setting setting shift (MRDR) is allowed to differ across the
* response-category cut-points, while the treatment effect (VRd) is kept parallel.
*
* MOTIVATION: in the proportional-odds fit the VRd efficacy-effectiveness gap is
* concentrated at the PR->SD boundary (PR collapses, SD grows) while the TOP of the
* distribution (CR/VGPR) holds up in routine care. A single shift cannot represent
* that - it drags CR/VGPR down too (over-correction; >=VGPR 43.5 -> 33.5 vs 43.6
* observed). Relaxing proportional odds for MRDR lets the shift differ by cut-point.
* NOTE the residual risk: the Rd anchor's OWN gap also spans the top (Rd >=VGPR
* 31.8 -> 24.1 in the registry) whereas VRd's does not, so even a threshold-specific
* shift transported from Rd may over-correct VRd's >=VGPR - a treatment-by-setting
* interaction a flexible link cannot fix. This script quantifies how much it helps.
*
* gologit2 can yield non-monotone (occasionally negative) predicted probabilities
* when cut-point lines cross; such replicates are guarded in the bootstrap. Requires
* the user-written gologit2 (auto-installed below).
*
*   Registry (MRDR, first line) - ACTUAL patient records
*       Event0==10 & Regimen==7   : Rd  anchor  (MRDR=1, VRd=0)
*       Event0==10 & Regimen==31  : VRd target  (held-out observed RW VRd)
*       6-level BCR collapsed to 5: recode BCR (1=1)(2=2)(3=3)(4 5=4)(6=5)  (MR->SD)
*   Trial (SWOG S0777) - SYNTHETIC counts
*       VRd n=216: 34 60 82 34 6 ;  Rd n=214: 18 50 85 52 9
*
* USAGE:
*   do calibrated_transport_prop.do 0            -> deterministic point estimate + table
*   do calibrated_transport_prop.do 1 500 71523  -> bootstrap (reps, seedbase)
**********

clear all
global boot     `1'
global nreps    `2'
global seedbase `3'
if "$boot"     == "" global boot     0
if "$nreps"    == "" global nreps    500
if "$seedbase" == "" global seedbase 71523

* ---- registry source (same file transport_dvd uses) ----
global REG "/Volumes/shared/R-MNHS-SPHPM-EPM-TRU/EpiMAP/Myeloma/Data/251128/MRDR Long MI.dta"

* ---- ensure gologit2 is available ----
cap which gologit2
if _rc ssc install gologit2, replace

**********
* Helpers
**********

* fill_bcr5 : build pseudo-patients with BCR from 5 ordered category counts (CR..PD)
*             (used for the SYNTHETIC SWOG arms only)
cap program drop fill_bcr5
program define fill_bcr5
    args c1 c2 c3 c4 c5
    clear
    local n = `c1' + `c2' + `c3' + `c4' + `c5'
    set obs `n'
    gen BCR = .
    local s = 1
    forvalues k = 1/5 {
        local cnt = `c`k''
        if `cnt' > 0 {
            local e = `s' + `cnt' - 1
            replace BCR = `k' in `s'/`e'
            local s = `e' + 1
        }
    }
end

* props5 : write the 5-category proportions of BCR (in memory) into a 1x5 matrix
cap program drop props5
program define props5
    args matname
    quietly count
    local n = r(N)
    matrix `matname' = J(1,5,0)
    forvalues k = 1/5 {
        quietly count if BCR == `k'
        matrix `matname'[1,`k'] = r(N) / `n'
    }
end

* run_vrd_prop : two synthetic SWOG arms + two REAL registry cohorts (optionally
*   within-cohort bootstrapped); partial proportional-odds fit (gologit2, npl(MRDR));
*   predicts the Calibrated-Transport cell (MRDR=1,VRd=1). Returns r(maeTr) r(maeCT),
*   r(bad)=1 if the prediction was invalid, and globals pVRd (trial), pCT, pObs.
cap program drop run_vrd_prop
program define run_vrd_prop, rclass
    args doboot rdfile vrdfile
    tempfile g

    * ---- SWOG VRd (novel) : MRDR=0, VRd=1   [synthetic] ----
    fill_bcr5 34 60 82 34 6
    if `doboot' bsample
    props5 pVRd                                   // trial-based prediction (unchanged)
    gen MRDR = 0
    gen VRd  = 1
    save `g', replace

    * ---- SWOG Rd (comparator, trial) : MRDR=0, VRd=0  (reference)  [synthetic] ----
    fill_bcr5 18 50 85 52 9
    if `doboot' bsample
    gen MRDR = 0
    gen VRd  = 0
    append using `g'
    save `g', replace

    * ---- MRDR Rd (registry anchor) : MRDR=1, VRd=0   [REAL records] ----
    use "`rdfile'", clear
    if `doboot' bsample
    gen MRDR = 1
    gen VRd  = 0
    append using `g'
    save `g', replace

    * ---- Observed target: MRDR VRd (held out, NOT in the fit)   [REAL records] ----
    use "`vrdfile'", clear
    if `doboot' bsample
    props5 pObs

    * ---- Fit partial proportional-odds model on the four groups ----
    use `g', clear
    label define BCR5 1 "CR" 2 "VGPR" 3 "PR" 4 "SD" 5 "PD", replace
    label values BCR BCR5
    *   npl(MRDR) : relax proportional odds for the setting shift only;
    *   VRd stays parallel. (Alternative: , autofit  to let Wald tests decide.)
    gologit2 BCR MRDR VRd, npl(MRDR)

    * ---- Predict the Calibrated-Transport cell (MRDR=1, VRd=1) ----
    set obs `=_N + 1'
    replace MRDR = 1 in L
    replace VRd  = 1 in L
    predict double _ctp1 _ctp2 _ctp3 _ctp4 _ctp5, pr
    matrix pCT = J(1,5,0)
    local bad = 0
    forvalues k = 1/5 {
        local v = _ctp`k'[_N]
        if (`v' < 0 | `v' > 1 | missing(`v')) local bad = 1
        matrix pCT[1,`k'] = `v'
    }
    return scalar bad = `bad'

    * ---- MAE (percentage points, mean over the 5 categories) ----
    scalar _maeTr = 0
    scalar _maeCT = 0
    forvalues k = 1/5 {
        scalar _maeTr = _maeTr + abs(pVRd[1,`k'] - pObs[1,`k'])
        scalar _maeCT = _maeCT + abs(pCT[1,`k']  - pObs[1,`k'])
    }
    return scalar maeTr = 100 * _maeTr / 5
    return scalar maeCT = 100 * _maeCT / 5
end

**********
* Build the two real registry cohorts once (5-level BCR), into tempfiles
**********
tempfile rdAnchor vrdTarget
use "$REG", clear
cap mi extract 0, clear                 // observed data (no covariates are imputed in this model)
keep if Event0 == 10                    // first line
recode BCR (1=1)(2=2)(3=3)(4 5=4)(6=5), gen(BCR5)   // fold MR(4) into SD(5)->4
drop BCR
rename BCR5 BCR
drop if missing(BCR)
preserve
    keep if Regimen == 7                // Rd anchor
    keep BCR
    save `rdAnchor', replace
restore
keep if Regimen == 31                   // VRd target
keep BCR
save `vrdTarget', replace

**********
* Paths
**********
local Local "/Users/adami/Documents/Monash/Vault/research/models/myeloma model/repo"
local Out   "`Local'/analyses/transport_vrd/results"
cap mkdir "`Out'"

**********
* Execute
**********
if ($boot == 0) {
    * ---------------- Deterministic point estimate ----------------
    run_vrd_prop 0 "`rdAnchor'" "`vrdTarget'"
    scalar maeTr = r(maeTr)
    scalar maeCT = r(maeCT)
    if r(bad) di as error "WARNING: partial-PO prediction contains an out-of-range probability (lines crossed)."

    local labs CR VGPR PR SD PD
    di as txt _n "{hline 70}"
    di as txt "First-line VRd: trial vs Calibrated Transport (PARTIAL PO) vs observed"
    di as txt "{hline 70}"
    di as txt %-8s "category" %12s "trial(VRd)" %16s "CT partial-PO" %12s "observed"
    forvalues k = 1/5 {
        local nm : word `k' of `labs'
        di as txt %-8s "`nm'" ///
           as res %11.1f 100*pVRd[1,`k'] %15.1f 100*pCT[1,`k'] %12.1f 100*pObs[1,`k']
    }
    scalar vgTr = 100*(pVRd[1,1]+pVRd[1,2])
    scalar vgCT = 100*(pCT[1,1]+pCT[1,2])
    scalar vgOb = 100*(pObs[1,1]+pObs[1,2])
    scalar orTr = 100*(pVRd[1,1]+pVRd[1,2]+pVRd[1,3])
    scalar orCT = 100*(pCT[1,1]+pCT[1,2]+pCT[1,3])
    scalar orOb = 100*(pObs[1,1]+pObs[1,2]+pObs[1,3])
    di as txt %-8s ">=VGPR" as res %11.1f vgTr %15.1f vgCT %12.1f vgOb
    di as txt %-8s "ORR"    as res %11.1f orTr %15.1f orCT %12.1f orOb
    di as txt "{hline 70}"
    di as txt "MAE vs observed (pp):  trial = " as res %4.2f maeTr ///
       as txt "   CT partial-PO = " as res %4.2f maeCT
    di as txt "MAE reduction (pp) = " as res %4.2f (maeTr - maeCT)
    di as txt "(compare proportional-odds CT from calibrated_transport.do)"
    di as txt "{hline 70}"

    file open fh using "`Out'/transport_vrd_prop.csv", write replace
    file write fh "category,trial_VRd,calibrated_transport_partialPO,observed_VRd" _n
    forvalues k = 1/5 {
        local nm : word `k' of `labs'
        local tv = 100*pVRd[1,`k']
        local cv = 100*pCT[1,`k']
        local ov = 100*pObs[1,`k']
        file write fh "`nm',`tv',`cv',`ov'" _n
    }
    file write fh ">=VGPR,`=vgTr',`=vgCT',`=vgOb'" _n
    file write fh "ORR,`=orTr',`=orCT',`=orOb'" _n
    file write fh "MAE_pp,`=maeTr',`=maeCT'," _n
    file close fh
    di as txt "Written: `Out'/transport_vrd_prop.csv"
    type "`Out'/transport_vrd_prop.csv"
}
else if ($boot == 1) {
    * ---------------- Bootstrap: CIs + exceedance ----------------
    set seed $seedbase
    tempname pf
    postfile `pf' rep maeTr maeCT using "`Out'/vrd_prop_boot.dta", replace
    local nbad = 0
    forvalues b = 1/$nreps {
        cap run_vrd_prop 1 "`rdAnchor'" "`vrdTarget'"
        if _rc == 0 & r(bad) == 0 post `pf' (`b') (r(maeTr)) (r(maeCT))
        else local nbad = `nbad' + 1
    }
    postclose `pf'

    use "`Out'/vrd_prop_boot.dta", clear
    gen redux = maeTr - maeCT
    gen byte ctwin = maeCT < maeTr
    di as txt _n "{hline 60}"
    di as txt "VRd Calibrated Transport (PARTIAL PO) - bootstrap ($nreps reps)"
    di as txt "  (`nbad' replicate(s) dropped for non-convergence/invalid prediction)"
    di as txt "{hline 60}"
    quietly summarize maeTr
    di as txt "trial MAE      mean: " as res %4.2f r(mean)
    quietly summarize maeCT
    di as txt "CT partial-PO  mean: " as res %4.2f r(mean)
    quietly _pctile maeCT, p(2.5 97.5)
    di as txt "CT partial-PO  95% CI: " as res %4.2f r(r1) " - " %4.2f r(r2)
    quietly _pctile redux, p(2.5 97.5)
    di as txt "MAE reduction 95% CI: " as res %4.2f r(r1) " - " %4.2f r(r2)
    quietly summarize ctwin
    di as txt "Partial-PO CT more accurate (lower MAE) in " ///
       as res %4.1f 100*r(mean) as txt "% of replicates"
    di as txt "{hline 60}"
}
