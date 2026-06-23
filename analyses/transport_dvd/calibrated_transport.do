**********
* Monash Myeloma Model - DVd L2 Calibrated Transport
* COMBINED trial-resample generator for scenarios A_trial and B_transport
* B model = relapsed-Vd anchor + LINE indicator (matches the CASTOR population)
*
* Same structure as transport.do (one shared trial resample per replicate feeds
* both A and B), but the B_transport fit pools ALL relapsed-setting Vd and adds a
* line indicator instead of using L2 Vd only.
*
*   A_trial      : resampled CASTOR BCR distributions (the trial-based prediction)
*                  -> mBCR_Vd_L2, mBCR_DVd_L2   (6-category proportions)
*   B_transport  : ologit BCR MRDR DVd line3 line4 line5  (observed data, m=0;
*                  no covariate is imputed in this model, so MI is not used)
*                  fit on {registry relapsed Vd (Regimen==5, Event0 20..60): MRDR=1, DVd=0}
*                        + {resampled CASTOR Vd  (n=231): MRDR=0, DVd=0}  <- ref cell (L2)
*                        + {resampled CASTOR DVd (n=238): MRDR=0, DVd=1}
*                  -> bL2_BCR_T = (beta_MRDR, beta_DVd, line3, line4, line5 | 5 cutpoints)
*
* WHY RELAPSED + LINE: CASTOR enrolled relapsed/refractory patients with >=1 prior
* line (median 2 prior, range 1-9; Palumbo 2016) - it is NOT L2-only. Restricting
* the anchor to L2 Vd (n=71) both under-powers it AND mismatches CASTOR's line mix,
* baking a line difference into beta_MRDR. Pooling all relapsed Vd (n=135: L2 71,
* L3 28, L4 14, L5 17, L6 5) with a line indicator fixes the mismatch and ~doubles
* the anchor. First-line Vd (n=240) is excluded (not CASTOR-eligible).
*
* CASTOR line of therapy is GENERATED from Palumbo 2016 Table 1 prior-line shares
* (1 prior -> L2, 2 -> L3, 3 -> L4, >3 -> L5+), drawn per pseudo-patient like a
* covariate: DVd arm 48.6/27.9/14.7/8.8% ; Vd arm 45.7/30.0/13.0/11.3%. CASTOR
* reports arm-level BCR only, so BCR is assigned independently of line; the line
* effect is identified from the registry, the CASTOR line draw balances the line
* mix for the MRDR comparison. Reference line = L2; L5+ lumped into line5.
*
* WHY COMBINED / SHARED RANDOMNESS: beta_DVd and the cutpoints come from the CASTOR
* arms, which in B drive the DVd-vs-Vd gap (DVd - Vd in scenario B = beta_DVd). If
* the trial arms were held fixed, B's bootstrap would omit trial sampling
* uncertainty and its CIs would be artificially narrow. Generating A and B from the
* SAME in-memory trial draw at iteration b (set seed = SEEDBASE + b once per
* replicate) locks their trial uncertainty together by construction. The registry
* anchor (MRDR Long MI B`b') is pre-resampled upstream and consumes no RNG here, so
* the trial draw is fully determined by b and identical for A and B.
*
* PREDICTION (sim_bcr_override.do, B_transport): predict DVd at L2 - MRDR=1, DVd=1,
* line3=line4=line5=0. The existing override is UNCHANGED: it builds a 2-column
* mPat (MRDR, DVd) and reads the first two coefficients + the final five cutpoints,
* so the line coefficients (cells 3-5) are skipped - exactly correct at L2.
*
* DETERMINISTIC PATH ($boot==0): no resample; exact CASTOR Table 2 counts and the
* full relapsed registry anchor. Deterministic outputs go to the EpiMAP-Local /
* repo outcomes tree (as transport.do); bootstrap outputs go to the repo outcomes
* tree (where the dispatcher's sim_bcr_override.do reads them).
*
* Usage:  do B_transport_dvd_new.do <boot> <min_bs> <max_bs>
**********

clear
clear mata
global boot   `1'
global min_bs `2'
global max_bs `3'

**********
* Helpers
**********

* fill_bcr : build pseudo-patients with BCR filled from 6 category counts (CR..PD)
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

* bcr_props : write the 6-category proportions of BCR (in memory) into a 1x6 matrix
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
    args t3 t4 t5                            // cumulative P(line<=2), <=3, <=4
    gen RN = runiform()
    gen Line = 2
    replace Line = 3 if RN > `t3'
    replace Line = 4 if RN > `t4'
    replace Line = 5 if RN > `t5'
    drop RN
end

* run_one : one configuration. Draws the two trial arms (optionally bootstrapped),
*           writes A's resampled distributions, then fits B's line-adjusted ologit
*           on the relapsed registry anchor + the SAME trial arms.
*   regsrc  - registry .dta (anchor source)
*   doboot  - 1 = within-arm bootstrap the trial arms; 0 = deterministic counts
*   a_vd    - A_trial output path for the Vd BCR distribution  (mBCR_Vd_L2)
*   a_dvd   - A_trial output path for the DVd BCR distribution (mBCR_DVd_L2)
*   b_out   - B_transport output path for the coefficients (bL2_BCR_T)
cap program drop run_one
program define run_one
    args regsrc doboot a_vd a_dvd b_out
    tempfile tvd tdvd

    * ---- CASTOR Vd arm (n=231, Table 2: CR 21, VG 47, PR 80, MR 20, SD 47, PD 16) ----
    fill_bcr 21 47 80 20 47 16
    if `doboot' bsample                     // within-arm resample (n preserved)
    bcr_props pVd                           // A: resampled Vd proportions
    gen MRDR = 0
    gen DVd  = 0
    gen_line 0.457 0.757 0.887              // Vd arm prior lines 45.7/30.0/13.0/11.3
    save `tvd'

    * ---- CASTOR DVd arm (n=238, Table 2: CR 46, VG 96, PR 57, MR 10, SD 24, PD 5) ----
    fill_bcr 46 96 57 10 24 5
    if `doboot' bsample                     // within-arm resample (n preserved)
    bcr_props pDVd                          // A: resampled DVd proportions
    gen MRDR = 0
    gen DVd  = 1
    gen_line 0.486 0.765 0.912              // DVd arm prior lines 48.6/27.9/14.7/8.8
    save `tdvd'

    * ---- A_trial outputs : the resampled trial BCR distributions ----
    mata: mBCR_Vd_L2  = st_matrix("pVd")
    mata: mBCR_DVd_L2 = st_matrix("pDVd")
    mata: mata matsave "`a_vd'"  mBCR_Vd_L2,  replace
    mata: mata matsave "`a_dvd'" mBCR_DVd_L2, replace

    * ---- B_transport : line-adjusted ologit on relapsed anchor + the SAME trial draw ----
    use "`regsrc'", clear
    cap mi extract 0, clear                       // observed data only — no covariate is imputed in this model
    keep if inlist(Event0, 20, 30, 40, 50, 60)   // L2..L6 starts
    keep if Regimen == 5                          // Vd
    gen MRDR = 1
    gen DVd  = 0
    replace Line = 5 if Line > 5                   // Line already defined in the registry; lump L6+ into L5
    append using "`tvd'"
    append using "`tdvd'"
    label define BCR_lbl 1 "CR" 2 "VG" 3 "PR" 4 "MR" 5 "SD" 6 "PD", replace
    label values BCR BCR_lbl
    gen line3 = (Line == 3)
    gen line4 = (Line == 4)
    gen line5 = (Line >= 5)                        // L5+

    ologit BCR MRDR DVd line3 line4 line5

    mata: bL2_BCR_T  = st_matrix("e(b)")
    mata: rbL2_BCR_T = st_matrixrowstripe("e(b)")
    mata: cbL2_BCR_T = st_matrixcolstripe("e(b)")
    forvalues i = 1/`=colsof(e(b))' {
        mata: cbL2_BCR_T[`i',1] = "`i'"
    }
    mata: mata matsave "`b_out'" bL2_BCR_T rbL2_BCR_T cbL2_BCR_T, replace

    matrix drop pVd pDVd
end

* mrdr_baseline : baseline-characteristics table for the MRDR Vd anchor cohort
*   (the SAME cohort the B_transport ologit is fitted on: relapsed-setting Vd,
*    Event0 in 20..60 & Regimen==5, n=135). Produces the MRDR Vd column of the
*    manuscript Table 1 (CASTOR DVd | CASTOR Vd | MRDR Vd). Computed on the
*    observed (non-imputed, m=0) data, the convention for a baseline table.
*
*   regsrc  - registry .dta (same deterministic source as the B fit)
*   outcsv  - CSV path for the baseline table
*
*   >>> ACTION REQUIRED: confirm the registry variable names in the CONFIRM
*       block below. Rows whose variable is left as "??" are written to the CSV
*       as "NEEDS VARIABLE" and skipped, so the script always runs end-to-end.
cap program drop mrdr_baseline
program define mrdr_baseline
    args regsrc outcsv regimen lines
    use "`regsrc'", clear

    * ---- baseline on observed data (m=0) before subsetting. `cap' makes this a
    *      no-op if the source is not mi-set. Change 0 to another m if preferred. ----
    cap mi extract 0, clear
	
	* --- prior treatment
	forval l = 10(10)50 {
		gen PI_`l' = 1 if (Bortezomib == 1 | Carfilzomib == 1) & Event0 == `l'
		gen IMID_`l'= 1 if (Lenalidomide == 1 | Thalidomide == 1 | Pomalidomide == 1) & Event0 == `l'
		gen Alk_`l' = 1 if (Cyclophosphamide == 1 | Melphalan == 1) & Event0 == `l'
	}
	foreach c in PI IMID Alk {
		forval l = 10(10)50 {
			bysort ID (`c'_`l'): replace `c'_`l' = `c'_`l'[_n-1] if `c'_`l' == .
		}
	}

    * ---- restrict to the requested cohort (BCR-evaluable) ----
    keep if inlist(Event0, `lines')               // line-of-therapy starts (Event0 = 10*line)
    keep if Regimen == `regimen'                  // 5 = Vd anchor ; DVd target = its own code
    drop if missing(BCR)                           // BCR-evaluable only (matches the analysis cohort)

    * ---- derived line variables (no confirmation needed; from Event0) ----
    gen Prior = Line - 1                           // prior lines of therapy: 1..5
	
	* --- prior exposure = class received at any line STRICTLY BEFORE the Vd line ---
	*     Event0 on the kept row IS the line at which this patient received Vd, so
	*     `l' < Event0 counts only earlier lines and excludes the Vd line itself
	*     (important: Vd's own bortezomib must not be counted as "prior PI").
	foreach c in PI IMID Alk {
		gen byte prior`c' = 0
		forval l = 10(10)50 {
			replace prior`c' = 1 if `c'_`l' == 1 & `l' < Event0
		}
	}
	gen byte priorPI_IMID = (priorPI == 1 & priorIMID == 1)

    **********************************************************************
    * >>>>>>>>>>>>>>>>>>>>>  CONFIRM VARIABLE NAMES  <<<<<<<<<<<<<<<<<<<<<<
    * Replace each "??" with the actual MRDR variable name. Leave as "??" to
    * skip that row. (Line / prior-lines are already derived from Event0.)
    **********************************************************************
    local vAge    Age     // age at the start of this line of therapy (years)
    local vSex    Male     // sex variable
    local maleval 1     // value of `vSex' that denotes MALE (e.g. 1, or "Male")
    local vISS    ISS     // ISS stage, expected coded 1 / 2 / 3  (NOT R-ISS)
    local vDxDate DateDN     // date of diagnosis            (Stata date, days)
    local vLoTDate Date0    // date of start of THIS line   (Stata date, days)
    * prior-therapy exposure flags, coded 1 = yes / 0 = no (leave "??" if absent):
    local vASCT   SCT        // prior ASCT  (CONFIRM: is SCT line-prior, or the Vd-line value? see note)
    local vPI     priorPI    // prior proteasome inhibitor   (built above from PI_`l')
    local vIMiD   priorIMID  // prior immunomodulatory drug  (built above from IMID_`l')
    local vAlk    priorAlk   // prior alkylating agent        (built above from Alk_`l')
    * Note: "Prior PI + IMiD" is derived below from `vPI' & `vIMiD' if both given.
    **********************************************************************

    * ---- open CSV (close any handle left open by a previous aborted run) ----
    cap file close fh
    file open fh using "`outcsv'", write replace
    file write fh "characteristic,value" _n

    quietly count
    local N = r(N)
    file write fh "N,`N'" _n

    * ---------- Age ----------
    if "`vAge'" != "??" {
        quietly summarize `vAge', detail
        local med = r(p50)
        local lo  = r(min)
        local hi  = r(max)
        file write fh "Age median (range) y,`med' (`lo'-`hi')" _n
        foreach band in "lt65 . 65" "65to74 65 75" "ge75 75 ." {
            gettoken nm rest : band
            gettoken a  b    : rest
            if "`a'"=="." quietly count if `vAge' <  `b'
            else if "`b'"=="." quietly count if `vAge' >= `a'
            else quietly count if `vAge' >= `a' & `vAge' < `b'
            file write fh "Age `nm' %,`=round(100*r(N)/`N',0.1)'" _n
        }
    }
    else file write fh "Age,NEEDS VARIABLE (vAge)" _n

    * ---------- Sex ----------
    if "`vSex'" != "??" & "`maleval'" != "??" {
        quietly count if `vSex' == `maleval'
        file write fh "Male %,`=round(100*r(N)/`N',0.1)'" _n
    }
    else file write fh "Male,NEEDS VARIABLE (vSex/maleval)" _n

    * ---------- ISS stage ----------
    if "`vISS'" != "??" {
        forvalues s = 1/3 {
            quietly count if `vISS' == `s'
            file write fh "ISS `s' %,`=round(100*r(N)/`N',0.1)'" _n
        }
        quietly count if missing(`vISS')
        file write fh "ISS missing %,`=round(100*r(N)/`N',0.1)'" _n
    }
    else file write fh "ISS,NEEDS VARIABLE (vISS)" _n

    * ---------- Time since diagnosis (years) ----------
    if "`vDxDate'" != "??" & "`vLoTDate'" != "??" {
        gen double _tsd = (`vLoTDate' - `vDxDate') / 365.25
        quietly summarize _tsd, detail
        file write fh "Time since diagnosis median (range) y,`=round(r(p50),0.1)' (`=round(r(min),0.1)'-`=round(r(max),0.1)')" _n
        drop _tsd
    }
    else file write fh "Time since diagnosis,NEEDS VARIABLE (vDxDate/vLoTDate)" _n

    * ---------- Prior lines of therapy (derived; no confirmation needed) ----------
    quietly summarize Prior, detail
    file write fh "Prior lines median (range),`=r(p50)' (`=r(min)'-`=r(max)')" _n
    foreach p in 1 2 3 {
        quietly count if Prior == `p'
        file write fh "Prior lines = `p' %,`=round(100*r(N)/`N',0.1)'" _n
    }
    quietly count if Prior >= 4
    file write fh "Prior lines >=4 %,`=round(100*r(N)/`N',0.1)'" _n

    * ---------- Prior therapy exposures (flags coded 1 = yes) ----------
    if "`vASCT'" != "??" {
        quietly count if `vASCT' == 1
        file write fh "Prior ASCT %,`=round(100*r(N)/`N',0.1)'" _n
    }
    else file write fh "Prior ASCT,NEEDS VARIABLE (vASCT)" _n

    if "`vPI'" != "??" {
        quietly count if `vPI' == 1
        file write fh "Prior PI %,`=round(100*r(N)/`N',0.1)'" _n
    }
    else file write fh "Prior PI,NEEDS VARIABLE (vPI)" _n

    if "`vIMiD'" != "??" {
        quietly count if `vIMiD' == 1
        file write fh "Prior IMiD %,`=round(100*r(N)/`N',0.1)'" _n
    }
    else file write fh "Prior IMiD,NEEDS VARIABLE (vIMiD)" _n

    if "`vAlk'" != "??" {
        quietly count if `vAlk' == 1
        file write fh "Prior Alkylator %,`=round(100*r(N)/`N',0.1)'" _n
    }
    else file write fh "Prior Alkylator,NEEDS VARIABLE (vAlk)" _n
    * Prior PI + IMiD (both), only if both variables confirmed
    if "`vPI'" != "??" & "`vIMiD'" != "??" {
        quietly count if `vPI' == 1 & `vIMiD' == 1
        file write fh "Prior PI + IMiD %,`=round(100*r(N)/`N',0.1)'" _n
    }
    else file write fh "Prior PI + IMiD,NEEDS VARIABLE (vPI & vIMiD)" _n

    file close fh

    * ---- also echo to the Results window so you can read/paste directly ----
    di as txt _n "{hline 60}"
    di as txt "MRDR baseline (BCR-evaluable, n=`N') written to:"
    di as res "  `outcsv'"
    di as txt "{hline 60}"
    type "`outcsv'"
end

**********
* Execute
**********
if ($boot == 0) {
    * ---- Deterministic point estimate: no trial resample (fixed seed for the
    *      CASTOR line draw, so the point estimate is reproducible) ----
    set seed 71523
    local Local "/Users/adami/Documents/Monash/Vault/research/models/myeloma model/repo"
    local Out   "`Local'/analyses/transport_dvd/outcomes"
    cap mkdir "`Out'/A_trial"

    run_one ///
        "/Volumes/shared/R-MNHS-SPHPM-EPM-TRU/EpiMAP/Myeloma/Data/251128/MRDR Long MI.dta" ///
        0 ///
        "`Out'/A_trial/bcr_vd_l2.mmat" ///
        "`Out'/A_trial/bcr_dvd_l2.mmat" ///
        "`Out'/B_transport/transport_dvd.mmat"

    * ---- Table 1 baselines (MRDR columns): only on the deterministic run ----
    *      comparator anchor (Vd, relapsed lines 2-6) + held-out novel target (DVd, L2).
    *      (fill the CONFIRM block inside mrdr_baseline; set DVDCODE below.)
    cap mkdir "`Local'/analyses/transport_dvd/results"
    local DVDCODE 80    // DVd's MRDR Regimen code (Vd is 5)
    mrdr_baseline ///
        "/Volumes/shared/R-MNHS-SPHPM-EPM-TRU/EpiMAP/Myeloma/Data/251128/MRDR Long MI.dta" ///
        "`Local'/analyses/transport_dvd/results/mrdr_vd_baseline.csv" ///
        5 20,30,40,50,60
    mrdr_baseline ///
        "/Volumes/shared/R-MNHS-SPHPM-EPM-TRU/EpiMAP/Myeloma/Data/251128/MRDR Long MI.dta" ///
        "`Local'/analyses/transport_dvd/results/mrdr_dvd_baseline.csv" ///
        `DVDCODE' 20
}
else if ($boot == 1) {
    * ---- Bootstrap: shared trial draw per replicate (seed = SEEDBASE + b) ----
    local SEEDBASE 71523
    forval b = $min_bs / $max_bs {
        set seed `=`SEEDBASE' + `b''
        run_one ///
            "~/em76/adam/data/251128/bootstrap/MRDR Long MI B`b'.dta" ///
            1 ///
            "outcomes/A_trial/bootstrap/bcr_vd_l2_B`b'.mmat" ///
            "outcomes/A_trial/bootstrap/bcr_dvd_l2_B`b'.mmat" ///
            "outcomes/B_transport/bootstrap/transport_dvd_B`b'.mmat"
    }
}
