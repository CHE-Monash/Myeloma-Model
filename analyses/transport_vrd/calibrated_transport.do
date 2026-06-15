**********
* EpiMAP Myeloma - VRd L1 Calibrated Transport  (SECOND EXEMPLAR)
* Trial-based vs Calibrated Transport prediction of real-world first-line VRd,
* validated against observed MRDR VRd. The registry cohorts are read from the ACTUAL
* MRDR patient records (as in transport_dvd's calibrated_transport.do); only the
* trial arms are synthetic, because SWOG is available only as aggregate counts.
*
*   Trial (SWOG S0777, Durie 2017, Lancet; Table 3, assessable) - SYNTHETIC counts
*       VRd  n=216 : CR34 VG60 PR82 SD34 PD6        (MRDR=0, VRd=1)
*       Rd   n=214 : CR18 VG50 PR85 SD52 PD9        (MRDR=0, VRd=0, reference)
*   Registry (MRDR, first line) - ACTUAL patient records
*       Event0==10 & Regimen==7   : Rd  anchor  (MRDR=1, VRd=0) -> identifies beta_MRDR
*       Event0==10 & Regimen==31  : VRd target  (held-out observed RW VRd)
*
* RESPONSE SCALE - 5 ORDERED CATEGORIES: CR > VGPR > PR > SD > PD.
*   SWOG does not report minimal response (MR), so the MRDR 6-category BCR
*   (1CR 2VG 3PR 4MR 5SD 6PD) is collapsed to 5 by folding MR into SD:
*       recode BCR (1=1)(2=2)(3=3)(4 5=4)(6=5)   ->  1CR 2VGPR 3PR 4SD(incl MR) 5PD
*   (This reproduces the tabulated anchor/target counts exactly: Rd 38/79/155/206/8,
*    n=486; VRd 239/422/352/489/14, n=1516.)
*
* MODEL:  ologit BCR MRDR VRd   (first line -> no line indicators)
*   beta_MRDR = trial->registry shift carried by the common comparator (Rd)
*   beta_VRd  = VRd-vs-Rd treatment effect (from the SWOG arms)
*   Calibrated-Transport prediction: real-world VRd at MRDR=1, VRd=1.
*   Trial-based prediction: the raw (resampled) SWOG VRd distribution.
*   Both scored against the observed MRDR VRd distribution.
*
* BOOTSTRAP: within-cohort resample of the two real registry cohorts AND the two
*   synthetic trial arms (n preserved); both predictions scored against the same
*   resampled observed target (paired) -> CIs + exceedance probability.
*
* USAGE:
*   do calibrated_transport.do 0            -> deterministic: point estimate, table,
*                                              AND the MRDR Rd 1L baseline for Table 1
*   do calibrated_transport.do 1 500 71523  -> bootstrap (reps, seedbase)
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

* run_vrd : one configuration. Builds the two synthetic SWOG arms + the two REAL
*   registry cohorts (optionally within-cohort bootstrapped), fits the ologit,
*   predicts the Calibrated-Transport cell (MRDR=1,VRd=1), returns r(maeTr) r(maeCT)
*   and the global matrices pVRd (trial), pCT (CT), pObs.
*   doboot - 1 = within-cohort bsample (n preserved); 0 = full cohorts
*   rdfile - .dta of the real MRDR Rd anchor (variable BCR, 5-level)
*   vrdfile- .dta of the real MRDR VRd target (variable BCR, 5-level)
cap program drop run_vrd
program define run_vrd, rclass
    args doboot rdfile vrdfile
    tempfile g

    * ---- SWOG VRd (novel) : MRDR=0, VRd=1   [synthetic] ----
    fill_bcr5 34 60 82 34 6
    if `doboot' bsample
    props5 pVRd                                   // trial-based prediction
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

    * ---- Fit on the four groups ----
    use `g', clear
    label define BCR5 1 "CR" 2 "VGPR" 3 "PR" 4 "SD" 5 "PD", replace
    label values BCR BCR5
    ologit BCR MRDR VRd

    * ---- Predict the Calibrated-Transport cell (MRDR=1, VRd=1) ----
    set obs `=_N + 1'
    replace MRDR = 1 in L
    replace VRd  = 1 in L
    predict double _ctp1 _ctp2 _ctp3 _ctp4 _ctp5, pr
    matrix pCT = J(1,5,0)
    forvalues k = 1/5 {
        matrix pCT[1,`k'] = _ctp`k'[_N]
    }

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

* mrdr_baseline_vrd : baseline characteristics for the MRDR Rd first-line anchor
*   (Event0==10 & Regimen==7) - the MRDR Rd column of the manuscript Table 1
*   (SWOG VRd | SWOG Rd | MRDR Rd). Computed on observed (m=0) data. Produces the four
*   characteristics SWOG reports and the MRDR can match: Age >=65, Male, ISS stage III,
*   ECOG >1 (each on its own non-missing denominator, the SWOG n/N convention). Mirrors
*   transport_dvd's mrdr_baseline.
*   regsrc - registry .dta ; regimen - MRDR regimen code (7 = Rd, 31 = VRd) ;
*   lab - short label for the screen banner ; outcsv - CSV path.
*   Restricted to BCR-evaluable patients so n matches the analysis cohort.
cap program drop mrdr_baseline_vrd
program define mrdr_baseline_vrd
    args regsrc regimen lab outcsv
    use "`regsrc'", clear
    cap mi extract 0, clear                 // observed data (baseline convention)
    keep if Event0 == 10                     // first line
    keep if Regimen == `regimen'             // 7 = Rd anchor ; 31 = VRd target
    drop if missing(BCR)                      // BCR-evaluable only (matches the analysis cohort)

    **********************************************************************
    * >>> CONFIRM VARIABLE NAMES (Age/Male/ISS match transport_dvd's baseline;
    *     confirm the ECOG variable name). Leave a name as "??" to skip that row.
    **********************************************************************
    local vAge    Age      // age at the start of this line (years)
    local vSex    Male     // sex variable
    local maleval 1        // value of `vSex' that denotes MALE
    local vISS    ISS      // ISS stage coded 1 / 2 / 3   (NOT R-ISS)
    local vECOG   ECOG     // ECOG performance status (>1 = grade 2+)   <-- CONFIRM NAME
    **********************************************************************

    cap file close fh
    file open fh using "`outcsv'", write replace
    file write fh "characteristic,value" _n
    quietly count
    file write fh "N (cohort),`r(N)'" _n

    * helper: write "lab %" as count(cond)/count(non-missing var)
    * ---- Age >=65 % ----
    if "`vAge'" != "??" {
        quietly count if !missing(`vAge')
        local d = r(N)
        quietly count if `vAge' >= 65 & !missing(`vAge')
        file write fh "Age >=65 years %,`=round(100*r(N)/`d',0.1)'" _n
    }
    else file write fh "Age >=65 years,NEEDS VARIABLE (vAge)" _n

    * ---- Male % ----
    if "`vSex'" != "??" & "`maleval'" != "??" {
        quietly count if !missing(`vSex')
        local d = r(N)
        quietly count if `vSex' == `maleval'
        file write fh "Male %,`=round(100*r(N)/`d',0.1)'" _n
    }
    else file write fh "Male,NEEDS VARIABLE (vSex/maleval)" _n

    * ---- ISS stage III % ----
    if "`vISS'" != "??" {
        quietly count if !missing(`vISS')
        local d = r(N)
        quietly count if `vISS' == 3
        file write fh "ISS stage III %,`=round(100*r(N)/`d',0.1)'" _n
    }
    else file write fh "ISS stage III,NEEDS VARIABLE (vISS)" _n

    * ---- ECOG >1 % ----
    if "`vECOG'" != "??" {
        capture confirm variable `vECOG'
        if !_rc {
            quietly count if !missing(`vECOG')
            local d = r(N)
            quietly count if `vECOG' > 1 & !missing(`vECOG')
            file write fh "ECOG >1 %,`=round(100*r(N)/`d',0.1)'" _n
        }
        else file write fh "ECOG >1,VARIABLE `vECOG' NOT FOUND - confirm name" _n
    }
    else file write fh "ECOG >1,NEEDS VARIABLE (vECOG)" _n

    file close fh
    di as txt _n "{hline 60}"
    di as txt "MRDR `lab' 1L baseline (BCR-evaluable) -> Table 1 (VRd panel):"
    di as res "  `outcsv'"
    di as txt "{hline 60}"
    type "`outcsv'"
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
    quietly count
    di as txt "MRDR Rd anchor (Event0==10 & Regimen==7): n = " r(N)
    save `rdAnchor', replace
restore
keep if Regimen == 31                   // VRd target
keep BCR
quietly count
di as txt "MRDR VRd target (Event0==10 & Regimen==31): n = " r(N)
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
    run_vrd 0 "`rdAnchor'" "`vrdTarget'"
    scalar maeTr = r(maeTr)
    scalar maeCT = r(maeCT)

    local labs CR VGPR PR SD PD
    di as txt _n "{hline 64}"
    di as txt "First-line VRd: trial-based vs Calibrated Transport vs observed"
    di as txt "{hline 64}"
    di as txt %-8s "category" %12s "trial(VRd)" %14s "CalibTransp" %12s "observed"
    forvalues k = 1/5 {
        local nm : word `k' of `labs'
        di as txt %-8s "`nm'" ///
           as res %11.1f 100*pVRd[1,`k'] %13.1f 100*pCT[1,`k'] %12.1f 100*pObs[1,`k']
    }
    scalar vgTr = 100*(pVRd[1,1]+pVRd[1,2])
    scalar vgCT = 100*(pCT[1,1]+pCT[1,2])
    scalar vgOb = 100*(pObs[1,1]+pObs[1,2])
    scalar orTr = 100*(pVRd[1,1]+pVRd[1,2]+pVRd[1,3])
    scalar orCT = 100*(pCT[1,1]+pCT[1,2]+pCT[1,3])
    scalar orOb = 100*(pObs[1,1]+pObs[1,2]+pObs[1,3])
    di as txt %-8s ">=VGPR" as res %11.1f vgTr %13.1f vgCT %12.1f vgOb
    di as txt %-8s "ORR"    as res %11.1f orTr %13.1f orCT %12.1f orOb
    di as txt "{hline 64}"
    di as txt "MAE vs observed (pp):  trial = " as res %4.2f maeTr ///
       as txt "   Calibrated Transport = " as res %4.2f maeCT
    di as txt "MAE reduction (pp) = " as res %4.2f (maeTr - maeCT)
    di as txt "{hline 64}"

    * ---- write a tidy CSV ----
    file open fh using "`Out'/transport_vrd.csv", write replace
    file write fh "category,trial_VRd,calibrated_transport,observed_VRd" _n
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
    di as txt "Written: `Out'/transport_vrd.csv"
    type "`Out'/transport_vrd.csv"

    * ---- Table 1 baselines (VRd panel): only on the deterministic run ----
    *      comparator anchor (Rd) and the held-out novel-treatment target (VRd)
    mrdr_baseline_vrd "$REG" 7  "Rd"  "`Out'/mrdr_rd_baseline.csv"
    mrdr_baseline_vrd "$REG" 31 "VRd" "`Out'/mrdr_vrd_baseline.csv"
}
else if ($boot == 1) {
    * ---------------- Bootstrap: CIs + exceedance ----------------
    set seed $seedbase
    tempname pf
    postfile `pf' rep maeTr maeCT using "`Out'/vrd_boot.dta", replace
    forvalues b = 1/$nreps {
        cap run_vrd 1 "`rdAnchor'" "`vrdTarget'"
        if _rc == 0 post `pf' (`b') (r(maeTr)) (r(maeCT))
    }
    postclose `pf'

    use "`Out'/vrd_boot.dta", clear
    gen redux = maeTr - maeCT
    gen byte ctwin = maeCT < maeTr
    di as txt _n "{hline 56}"
    di as txt "VRd Calibrated Transport - bootstrap ($nreps reps)"
    di as txt "{hline 56}"
    quietly summarize maeTr
    di as txt "trial MAE  mean: " as res %4.2f r(mean)
    quietly summarize maeCT
    di as txt "CT    MAE  mean: " as res %4.2f r(mean)
    quietly _pctile maeTr, p(2.5 97.5)
    di as txt "trial MAE  95% CI: " as res %4.2f r(r1) " - " %4.2f r(r2)
    quietly _pctile maeCT, p(2.5 97.5)
    di as txt "CT    MAE  95% CI: " as res %4.2f r(r1) " - " %4.2f r(r2)
    quietly _pctile redux, p(2.5 97.5)
    di as txt "MAE reduction 95% CI: " as res %4.2f r(r1) " - " %4.2f r(r2)
    quietly summarize ctwin
    di as txt "Calibrated Transport more accurate (lower MAE) in " ///
       as res %4.1f 100*r(mean) as txt "% of replicates"
    di as txt "{hline 56}"
}
