**********
* Monash Myeloma Model - DVd L2 Calibrated Transport - proportional-odds diagnostics
*
* Purpose: Companion to calibrated_transport.do (same data and design) that interrogates the
*          proportional-odds (parallel-lines) assumption the single-coefficient correction rests
*          on, focusing on whether the SETTING shift (MRDR) is proportional across the five BCR
*          cut-points. Three deterministic diagnostics: (1) Brant test (needs SPost13 brant;
*          skipped with a note if absent); (2) gologit2 npl() LR tests of the parallel model vs a
*          generalised ordered logit with MRDR / DVd / both freed; (3) predicted-distribution
*          comparison at the CT cell (MRDR=1, DVd=1, L2) - trial-based vs CT under proportional
*          odds vs CT with the setting shift freed vs both freed, against the observed DVd marginal.
* Usage:   do calibrated_transport_prop.do 0            -> deterministic tests + table
*          do calibrated_transport_prop.do 1 500 71523  -> bootstrap (reps, seedbase)
* Notes:   Bootstrap mode ($boot==1) within-cohort resamples the partial-PO (npl MRDR) vs
*          proportional-odds CT MAE for a CI on the MAE change and the share of replicates it helps;
*          self-contained (resamples the registry anchor in memory, no HPC tree). gologit2 can give
*          out-of-range predicted probabilities when cut-point lines cross; such replicates are
*          flagged (r(bad)) and dropped. n=135 anchor -> low power; read a non-significant LR test
*          as weak reassurance, not proof.
**********

clear all
global boot     `1'
global nreps    `2'
global seedbase `3'
if "$boot"     == "" global boot     0
if "$nreps"    == "" global nreps    500
if "$seedbase" == "" global seedbase 71523

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"                 // machine-specific paths: $data_path

* ---- registry source (same file calibrated_transport.do uses) ----
local REG "${data_path}/MRDR Long MI.dta"
local DVDCODE 80                        // DVd's MRDR Regimen code (Vd is 5)
local Out "analyses/transport_dvd/results"
cap mkdir "`Out'"

* ---- ensure gologit2 is available (Brant is optional; see note above) ----
cap which gologit2
if _rc ssc install gologit2, replace

**********
* Helpers (identical to calibrated_transport.do: 6-level BCR)
**********

* fill_bcr : pseudo-patients with BCR filled from 6 category counts (CR..PD)
cap program drop fill_bcr
program define fill_bcr
    args c1 c2 c3 c4 c5 c6
    clear
    local n = `c1' + `c2' + `c3' + `c4' + `c5' + `c6'
    set obs `n'
    gen BCR = .
    local start = 1
    forvalues k = 1/6 {
        local cnt = `c`k''
        if `cnt' > 0 {
            local end = `start' + `cnt' - 1
            replace BCR = `k' in `start'/`end'
            local start = `end' + 1
        }
    }
end

* bcr_props : 6-category proportions of BCR (in memory) into a 1x6 matrix
cap program drop bcr_props
program define bcr_props
    args matname
    quietly count
    local n = r(N)
    matrix `matname' = J(1,6,0)
    forvalues k = 1/6 {
        quietly count if BCR == `k'
        matrix `matname'[1,`k'] = r(N) / `n'
    }
end

* gen_line : draw a line of therapy (2..5; 5 = L5+) from cumulative thresholds
cap program drop gen_line
program define gen_line
    args t3 t4 t5
    gen RN = runiform()
    gen Line = 2
    replace Line = 3 if RN > `t3'
    replace Line = 4 if RN > `t4'
    replace Line = 5 if RN > `t5'
    drop RN
end

* run_dvd_prop : build the four cohorts (optionally within-cohort bootstrapped),
*   fit the proportional-odds CT model (ologit) and the partial-PO variant
*   (gologit2, npl(MRDR)), predict both at the CT cell (MRDR=1,DVd=1,L2), and
*   return MAEs vs the held-out observed DVd marginal. With `dotests'==1 it also
*   runs the Brant test and the gologit2 LR tests on the point estimate and leaves
*   the predicted-distribution matrices (pDVd_, pObs_, pCTpl, pCTnp, pCTboth) in
*   memory for the caller's table.
cap program drop run_dvd_prop
program define run_dvd_prop, rclass
    args doboot anchorfile targetfile dotests
    tempfile g

    * ---- CASTOR Vd arm (n=231) : MRDR=0, DVd=0 (reference)  [synthetic] ----
    fill_bcr 21 47 80 20 47 16
    if `doboot' bsample
    gen MRDR = 0
    gen DVd  = 0
    gen_line 0.457 0.757 0.887
    save `g'

    * ---- CASTOR DVd arm (n=238) : MRDR=0, DVd=1            [synthetic] ----
    fill_bcr 46 96 57 10 24 5
    if `doboot' bsample
    bcr_props pDVd_                                 // trial-based DVd prediction
    gen MRDR = 0
    gen DVd  = 1
    gen_line 0.486 0.765 0.912
    append using `g'
    save `g', replace

    * ---- MRDR Vd anchor (relapsed, lines 2-6) : MRDR=1, DVd=0   [REAL records] ----
    use "`anchorfile'", clear
    if `doboot' bsample
    gen MRDR = 1
    gen DVd  = 0
    append using `g'
    save `g', replace

    * ---- Observed target: MRDR DVd L2 (held out, NOT in the fit)  [REAL records] ----
    use "`targetfile'", clear
    if `doboot' bsample
    bcr_props pObs_

    * ---- assemble the fit data + the CT-cell prediction row ----
    use `g', clear
    gen line3 = (Line == 3)
    gen line4 = (Line == 4)
    gen line5 = (Line >= 5)
    label define BCR_lbl 1 "CR" 2 "VG" 3 "PR" 4 "MR" 5 "SD" 6 "PD", replace
    label values BCR BCR_lbl
    set obs `=_N + 1'                                // CT cell: MRDR=1, DVd=1, L2
    replace MRDR  = 1 in L
    replace DVd   = 1 in L
    replace line3 = 0 in L
    replace line4 = 0 in L
    replace line5 = 0 in L

    * ---- (a) proportional-odds CT == the main model ----
    ologit BCR MRDR DVd line3 line4 line5
    predict double _pl1 _pl2 _pl3 _pl4 _pl5 _pl6, pr
    matrix pCTpl = J(1,6,0)
    forvalues k = 1/6 {
        matrix pCTpl[1,`k'] = _pl`k'[_N]
    }
    drop _pl1 _pl2 _pl3 _pl4 _pl5 _pl6

    * ---- (b) partial-PO CT: setting shift freed ----
    gologit2 BCR MRDR DVd line3 line4 line5, npl(MRDR)
    predict double _np1 _np2 _np3 _np4 _np5 _np6, pr
    matrix pCTnp = J(1,6,0)
    local bad = 0
    forvalues k = 1/6 {
        local v = _np`k'[_N]
        if (`v' < 0 | `v' > 1 | missing(`v')) local bad = 1
        matrix pCTnp[1,`k'] = `v'
    }
    drop _np1 _np2 _np3 _np4 _np5 _np6
    return scalar bad = `bad'

    * ---- MAEs (pp, mean over 6 categories) vs observed ----
    scalar _mTr = 0
    scalar _mPl = 0
    scalar _mNp = 0
    forvalues k = 1/6 {
        scalar _mTr = _mTr + abs(pDVd_[1,`k'] - pObs_[1,`k'])
        scalar _mPl = _mPl + abs(pCTpl[1,`k']  - pObs_[1,`k'])
        scalar _mNp = _mNp + abs(pCTnp[1,`k']  - pObs_[1,`k'])
    }
    return scalar maeTr = 100 * _mTr / 6
    return scalar maePl = 100 * _mPl / 6
    return scalar maeNp = 100 * _mNp / 6

    * ---- (c) point-estimate tests (deterministic run only) ----
    if `dotests' {
        * parallel baseline for the LR tests (gologit2 ..., pl == ologit)
        gologit2 BCR MRDR DVd line3 line4 line5, pl
        estimates store _gpl

        gologit2 BCR MRDR DVd line3 line4 line5, npl(MRDR)
        lrtest _gpl ., force
        return scalar lr_mrdr_chi2 = r(chi2)
        return scalar lr_mrdr_df   = r(df)
        return scalar lr_mrdr_p    = r(p)

        gologit2 BCR MRDR DVd line3 line4 line5, npl(DVd)
        lrtest _gpl ., force
        return scalar lr_dvd_chi2 = r(chi2)
        return scalar lr_dvd_df   = r(df)
        return scalar lr_dvd_p    = r(p)

        gologit2 BCR MRDR DVd line3 line4 line5, npl(MRDR DVd)
        predict double _bo1 _bo2 _bo3 _bo4 _bo5 _bo6, pr
        matrix pCTboth = J(1,6,0)
        forvalues k = 1/6 {
            matrix pCTboth[1,`k'] = _bo`k'[_N]
        }
        drop _bo1 _bo2 _bo3 _bo4 _bo5 _bo6
        lrtest _gpl ., force
        return scalar lr_both_chi2 = r(chi2)
        return scalar lr_both_df   = r(df)
        return scalar lr_both_p    = r(p)

        * Brant test (optional; needs SPost13's brant). Overall + per-variable.
        ologit BCR MRDR DVd line3 line4 line5
        cap brant, detail
        if _rc == 0 {
            return scalar brant_chi2 = r(chi2)
            return scalar brant_df   = r(df)
            return scalar brant_p    = r(p)
            return scalar brant_ok   = 1
        }
        else return scalar brant_ok = 0
    }
end

**********
* Build the two real registry cohorts once (6-level BCR), into tempfiles
*   anchor : relapsed Vd, lines 2-6  (Event0 20..60 & Regimen==5), keeps BCR + Line
*   target : DVd L2                   (Event0==20    & Regimen==`DVDCODE'), keeps BCR
* Line is the registry's own variable (as in calibrated_transport.do); L6+ lumped into 5.
**********
tempfile anchorVd targetDVd
use "`REG'", clear
cap mi extract 0, clear                            // observed data (no covariate imputed here)
preserve
    keep if inlist(Event0, 20, 30, 40, 50, 60)
    keep if Regimen == 5
    drop if missing(BCR)
    replace Line = 5 if Line > 5                    // use the registry Line (as calibrated_transport.do); lump L6+ into 5
    keep BCR Line
    save `anchorVd', replace
restore
keep if Event0 == 20
keep if Regimen == `DVDCODE'
drop if missing(BCR)
keep BCR
save `targetDVd', replace

**********
* Execute
**********
if ($boot == 0) {
    * ---------------- Deterministic: tests + predicted-distribution table -------
    set seed $seedbase
    run_dvd_prop 0 "`anchorVd'" "`targetDVd'" 1
    * ---- stash all returns into scalars immediately (guard against r() being cleared) ----
    scalar maeTr = r(maeTr)
    scalar maePl = r(maePl)
    scalar maeNp = r(maeNp)
    scalar badf  = r(bad)
    scalar lrMc = r(lr_mrdr_chi2)
    scalar lrMd = r(lr_mrdr_df)
    scalar lrMp = r(lr_mrdr_p)
    scalar lrDc = r(lr_dvd_chi2)
    scalar lrDd = r(lr_dvd_df)
    scalar lrDp = r(lr_dvd_p)
    scalar lrBc = r(lr_both_chi2)
    scalar lrBd = r(lr_both_df)
    scalar lrBp = r(lr_both_p)
    scalar brOK = r(brant_ok)
    scalar brc  = r(brant_chi2)
    scalar brd  = r(brant_df)
    scalar brp  = r(brant_p)
    if badf di as error "WARNING: partial-PO prediction out of range (cut-point lines crossed)."

    * ---- proportional-odds tests ----
    di as txt _n "{hline 74}"
    di as txt "Proportional-odds tests for the DVd Calibrated-Transport model"
    di as txt "  ologit BCR MRDR DVd line3 line4 line5   (n anchor=135; total fit ~604)"
    di as txt "{hline 74}"
    di as txt "LR test (parallel vs term freed across cut-points):"
    di as txt "  setting shift MRDR : chi2(" as res lrMd as txt ") = " ///
       as res %6.2f lrMc as txt ",  p = " as res %5.3f lrMp
    di as txt "  treatment DVd      : chi2(" as res lrDd as txt ") = " ///
       as res %6.2f lrDc as txt ",  p = " as res %5.3f lrDp
    di as txt "  both (MRDR + DVd)  : chi2(" as res lrBd as txt ") = " ///
       as res %6.2f lrBc as txt ",  p = " as res %5.3f lrBp
    if brOK {
        di as txt "Brant test (overall) : chi2(" as res brd as txt ") = " ///
           as res %6.2f brc as txt ",  p = " as res %5.3f brp
    }
    else di as txt "Brant test           : skipped (SPost13 'brant' not installed; LR tests above suffice)"
    di as txt "  (significant p rejects proportional odds for that term; n=135 -> low power)"
    di as txt "{hline 74}"

    * ---- predicted-distribution comparison at the CT cell (MRDR=1,DVd=1,L2) ----
    local labs CR VGPR PR MR SD PD
    di as txt _n "{hline 86}"
    di as txt "Predicted BCR at DVd L2: trial vs CT (proportional odds) vs CT (npl) vs observed"
    di as txt "{hline 86}"
    di as txt %-8s "category" %10s "trial" %14s "CT par-PO" %14s "CT npl(MRDR)" %14s "CT npl(both)" %11s "observed"
    forvalues k = 1/6 {
        local nm : word `k' of `labs'
        di as txt %-8s "`nm'" ///
           as res %10.1f 100*pDVd_[1,`k'] %14.1f 100*pCTpl[1,`k'] ///
           %14.1f 100*pCTnp[1,`k'] %14.1f 100*pCTboth[1,`k'] %11.1f 100*pObs_[1,`k']
    }
    * >=VGPR summary row
    scalar vgTr = 100*(pDVd_[1,1]+pDVd_[1,2])
    scalar vgPl = 100*(pCTpl[1,1]+pCTpl[1,2])
    scalar vgNp = 100*(pCTnp[1,1]+pCTnp[1,2])
    scalar vgBo = 100*(pCTboth[1,1]+pCTboth[1,2])
    scalar vgOb = 100*(pObs_[1,1]+pObs_[1,2])
    di as txt %-8s ">=VGPR" as res %10.1f vgTr %14.1f vgPl %14.1f vgNp %14.1f vgBo %11.1f vgOb
    di as txt "{hline 86}"
    di as txt "MAE vs observed (pp):  trial = " as res %4.2f maeTr ///
       as txt " | CT par-PO = " as res %4.2f maePl ///
       as txt " | CT npl(MRDR) = " as res %4.2f maeNp
    di as txt "Change from relaxing PO (par-PO - npl) = " as res %4.2f (maePl - maeNp) as txt " pp"
    di as txt "(compare the proportional-odds CT in calibrated_transport.do / Figure 1)"
    di as txt "{hline 86}"

    * ---- write CSVs ----
    file open fd using "`Out'/po_diagnostics_dist.csv", write replace
    file write fd "category,trial_DVd,CT_parallelPO,CT_npl_MRDR,CT_npl_both,observed_DVd" _n
    forvalues k = 1/6 {
        local nm : word `k' of `labs'
        file write fd "`nm',`=100*pDVd_[1,`k']',`=100*pCTpl[1,`k']',`=100*pCTnp[1,`k']',`=100*pCTboth[1,`k']',`=100*pObs_[1,`k']'" _n
    }
    file write fd ">=VGPR,`=vgTr',`=vgPl',`=vgNp',`=vgBo',`=vgOb'" _n
    file write fd "MAE_pp,`=maeTr',`=maePl',`=maeNp',," _n
    file close fd

    file open ft using "`Out'/po_diagnostics_tests.csv", write replace
    file write ft "test,term,chi2,df,p" _n
    file write ft "LR,MRDR,`=lrMc',`=lrMd',`=lrMp'" _n
    file write ft "LR,DVd,`=lrDc',`=lrDd',`=lrDp'" _n
    file write ft "LR,MRDR+DVd,`=lrBc',`=lrBd',`=lrBp'" _n
    if brOK file write ft "Brant,overall,`=brc',`=brd',`=brp'" _n
    else file write ft "Brant,overall,NA,NA,NA (SPost13 not installed)" _n
    file close ft

    di as txt "Written: `Out'/po_diagnostics_dist.csv  and  `Out'/po_diagnostics_tests.csv"
    type "`Out'/po_diagnostics_dist.csv"
    type "`Out'/po_diagnostics_tests.csv"
}
else if ($boot == 1) {
    * ---------------- Bootstrap: CI on the MAE change from relaxing PO ----------
    set seed $seedbase
    tempname pf
    postfile `pf' rep maeTr maePl maeNp using "`Out'/po_prop_boot.dta", replace
    local nbad = 0
    forvalues b = 1/$nreps {
        cap run_dvd_prop 1 "`anchorVd'" "`targetDVd'" 0
        if _rc == 0 & r(bad) == 0 post `pf' (`b') (r(maeTr)) (r(maePl)) (r(maeNp))
        else local nbad = `nbad' + 1
    }
    postclose `pf'

    use "`Out'/po_prop_boot.dta", clear
    gen redux  = maePl - maeNp                       // PO MAE minus partial-PO MAE
    gen byte npwin = maeNp < maePl                    // partial-PO closer than PO
    di as txt _n "{hline 64}"
    di as txt "DVd partial-PO (npl MRDR) vs proportional-odds CT - bootstrap ($nreps reps)"
    di as txt "  (`nbad' replicate(s) dropped for non-convergence/invalid prediction)"
    di as txt "{hline 64}"
    quietly summarize maePl
    di as txt "proportional-odds CT MAE  mean: " as res %4.2f r(mean)
    quietly summarize maeNp
    di as txt "partial-PO (npl MRDR) MAE mean: " as res %4.2f r(mean)
    quietly _pctile maeNp, p(2.5 97.5)
    di as txt "partial-PO MAE 95% CI: " as res %4.2f r(r1) " - " %4.2f r(r2)
    quietly _pctile redux, p(2.5 97.5)
    di as txt "MAE change (PO - partial) 95% CI: " as res %4.2f r(r1) " - " %4.2f r(r2)
    quietly summarize npwin
    di as txt "Partial-PO closer than PO in " as res %4.1f 100*r(mean) as txt "% of replicates"
    di as txt "{hline 64}"
}
