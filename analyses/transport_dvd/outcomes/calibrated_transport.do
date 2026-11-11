**********
* Monash Myeloma Model - DVd L2 Calibrated Transport
*
* Purpose: Combined trial-resample generator for scenarios A_trial and B_transport (one shared
*          trial resample per replicate feeds both). A_trial -> resampled CASTOR BCR distributions
*          (mBCR_Vd_L2, mBCR_DVd_L2). B_transport -> ologit BCR MRDR DVd line3 line4 line5, fit on
*          registry relapsed Vd (MRDR=1) + resampled CASTOR Vd (ref cell, L2) + resampled CASTOR
*          DVd, giving bL2_BCR_T = (beta_MRDR, beta_DVd, line3-5 | 5 cutpoints).
* Usage:   do calibrated_transport.do <boot> <min_bs> <max_bs>
* Notes:   B pools ALL relapsed Vd (n=135) with a line indicator instead of L2 Vd only (n=71) to
*          match CASTOR's line mix (relapsed, >=1 prior; Palumbo 2016) and avoid baking a line
*          difference into beta_MRDR; first-line Vd excluded. CASTOR line is generated from Palumbo
*          Table 1 prior-line shares, drawn per pseudo-patient; ref line = L2, L5+ -> line5. A and B
*          share the SAME in-memory trial draw per iteration (seed = SEEDBASE + b) so their trial
*          uncertainty is locked together; the registry anchor is pre-resampled upstream and uses no
*          RNG here. Prediction (sim_bcr_override.do, B_transport): DVd at L2 = MRDR=1, DVd=1, lines=0;
*          the override reads the first two coefficients + five cutpoints, skipping the line cells
*          (correct at L2). $boot==0 deterministic path uses exact CASTOR Table 2 counts and the full
*          registry anchor; deterministic outputs -> repo/EpiMAP-Local outcomes tree, bootstrap
*          outputs -> repo outcomes tree.
**********

clear
clear mata
global boot   `1'
global min_bs `2'
global max_bs `3'

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $data_path (git-ignored; see config.example.do)

**********
* Helpers
**********

* export_coefs : write the full fitted regression table to a tidy CSV — every
*   parameter (MRDR, DVd, line3-5 and the five cut-points) with coefficient,
*   standard error, z, p and 95% CI. Reads r(table) (captured into `rtname')
*   and the e(b) column names. Feeds Supplementary Table S3.
mata:
void export_coefs(string scalar rtname, string scalar path)
{
    real matrix    rt
    string matrix  cs
    real scalar    j, fh
    rt = st_matrix(rtname)                         // rows: 1 b, 2 se, 3 z, 4 p, 5 ll, 6 ul
    cs = st_matrixcolstripe("e(b)")                // col 2 = parameter name
    unlink(path)                                    // fopen("w") errors if the file exists
    fh = fopen(path, "w")
    fput(fh, "term,coef,se,z,p,ci_lower,ci_upper")
    for (j=1; j<=cols(rt); j++) {
        fput(fh, cs[j,2] + "," +
            strtrim(strofreal(rt[1,j], "%12.6f")) + "," +
            strtrim(strofreal(rt[2,j], "%12.6f")) + "," +
            strtrim(strofreal(rt[3,j], "%12.4f")) + "," +
            strtrim(strofreal(rt[4,j], "%12.3g")) + "," +
            strtrim(strofreal(rt[5,j], "%12.6f")) + "," +
            strtrim(strofreal(rt[6,j], "%12.6f")))
    }
    fclose(fh)
}
end

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
    args regsrc doboot a_vd a_dvd b_out coef_csv
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
    matrix rtab = r(table)                          // full coef table (b/se/z/p/CI), all terms

    mata: bL2_BCR_T  = st_matrix("e(b)")
    mata: rbL2_BCR_T = st_matrixrowstripe("e(b)")
    mata: cbL2_BCR_T = st_matrixcolstripe("e(b)")
    forvalues i = 1/`=colsof(e(b))' {
        mata: cbL2_BCR_T[`i',1] = "`i'"
    }
    mata: mata matsave "`b_out'" bL2_BCR_T rbL2_BCR_T cbL2_BCR_T, replace

    * ---- full regression table -> results CSV (point-estimate run only; the
    *      bootstrap calls pass no `coef_csv', so this is skipped there) ----
    if "`coef_csv'" != "" {
        mata: export_coefs("rtab", "`coef_csv'")
    }

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

    * ---------- Period: treatment-start dates of this cohort (for the Methods, which
    *            state the calibration and validation periods; co-author comment C009) ----------
    if "`vLoTDate'" != "??" {
        quietly summarize `vLoTDate'
        local d0 : display %tdDD_Mon_CCYY r(min)
        local d1 : display %tdDD_Mon_CCYY r(max)
        file write fh "Treatment start first,`=trim("`d0'")'" _n
        file write fh "Treatment start last,`=trim("`d1'")'" _n
        gen int _yr = year(`vLoTDate')
        quietly levelsof _yr, local(yrs)
        foreach y of local yrs {
            quietly count if _yr == `y'
            file write fh "Treatment start in `y' n,`=r(N)'" _n
        }
        drop _yr
    }
    else file write fh "Treatment start period,NEEDS VARIABLE (vLoTDate)" _n

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
    local Out   "analyses/transport_dvd/outcomes"
    cap mkdir "`Out'/A_trial"
    cap mkdir "analyses/transport_dvd/results"

    run_one ///
        "${data_path}/MRDR Long MI.dta" ///
        0 ///
        "`Out'/A_trial/bcr_vd_l2.mmat" ///
        "`Out'/A_trial/bcr_dvd_l2.mmat" ///
        "`Out'/B_transport/transport_dvd.mmat" ///
        "analyses/transport_dvd/results/transport_dvd_coefs.csv"

    * ---- Table 1 baselines (MRDR columns): only on the deterministic run ----
    *      comparator anchor (Vd, relapsed lines 2-6) + held-out novel target (DVd, L2).
    *      (fill the CONFIRM block inside mrdr_baseline; set DVDCODE below.)
    cap mkdir "analyses/transport_dvd/results"
    local DVDCODE 80    // DVd's MRDR Regimen code (Vd is 5)
    mrdr_baseline ///
        "${data_path}/MRDR Long MI.dta" ///
        "analyses/transport_dvd/results/mrdr_vd_baseline.csv" ///
        5 20,30,40,50,60
    mrdr_baseline ///
        "${data_path}/MRDR Long MI.dta" ///
        "analyses/transport_dvd/results/mrdr_dvd_baseline.csv" ///
        `DVDCODE' 20

    * ---- Observed DVd BCR marginal for the held-out validation cohort (n + 6-cat
    *      counts/%). The HPC bootstrap_summary.do must HARD-CODE this benchmark (no
    *      registry access on the cluster); this writes results/observed_dvd_bcr.csv
    *      and prints the `local' lines ready to paste there. Same cohort as
    *      mrdr_dvd_baseline.csv (Event0==20, Regimen==DVDCODE, BCR-evaluable). ----
    use "${data_path}/MRDR Long MI.dta", clear
    cap mi extract 0, clear
    keep if Event0 == 20
    keep if Regimen == `DVDCODE'
    drop if missing(BCR)
    qui count
    local oN = r(N)
    cap file close obf
    file open obf using "analyses/transport_dvd/results/observed_dvd_bcr.csv", write replace
    file write obf "bcr,count,pct" _n
    local oblabs CR VGPR PR MR SD PD
    di as txt _n "{hline 64}"
    di as txt "Observed DVd BCR benchmark (n=`oN') -- paste into simulated/bootstrap_summary.do:"
    di as res "    local oN = `oN'"
    forval k = 1/6 {
        qui count if BCR == `k'
        local cnt = r(N)
        local pct = 100*`cnt'/`oN'
        local oblab : word `k' of `oblabs'
        file write obf "`oblab',`cnt',`=string(`pct',"%6.2f")'" _n
        di as res "    local o`k' = " %6.2f `pct'
    }
    file write obf "N,`oN'," _n
    file close obf
    di as txt "{hline 64}"

    * ---- Observed Vd BCR marginal for the calibration cohort (n + 6-cat counts/%).
    *      Same cohort as mrdr_vd_baseline.csv and as the anchor the model is fitted on:
    *      relapsed-setting Vd (Event0 20-60, Regimen==5), BCR-evaluable -> n=93.
    *      Written for Figure 1, which shows the efficacy-effectiveness gap measured on
    *      the comparator (trial Vd vs registry Vd) beside the gap predicted for DVd.
    *      NOTE: this is the OBSERVED registry marginal. It is not the "vd" rows of
    *      results/bcr_distributions.csv, which are the simulated Vd arm of the CE
    *      cohort and do not reproduce the fitted setting shift. ----
    use "${data_path}/MRDR Long MI.dta", clear
    cap mi extract 0, clear
    keep if inlist(Event0, 20, 30, 40, 50, 60)
    keep if Regimen == 5
    drop if missing(BCR)
    qui count
    local vN = r(N)
    cap file close vbf
    file open vbf using "analyses/transport_dvd/results/observed_vd_bcr.csv", write replace
    file write vbf "bcr,count,pct" _n
    local vblabs CR VGPR PR MR SD PD
    di as txt _n "{hline 64}"
    di as txt "Observed Vd BCR marginal, calibration cohort (n=`vN') -- for Figure 1:"
    forval k = 1/6 {
        qui count if BCR == `k'
        local cnt = r(N)
        local pct = 100*`cnt'/`vN'
        local vblab : word `k' of `vblabs'
        file write vbf "`vblab',`cnt',`=string(`pct',"%6.2f")'" _n
        di as res "    `vblab' = " %6.2f `pct' "%  (`cnt')"
    }
    file write vbf "N,`vN'," _n
    file close vbf
    di as txt "{hline 64}"
}
else if ($boot == 1) {
    * ---- Bootstrap: shared trial draw per replicate (seed = SEEDBASE + b) ----
    local SEEDBASE 71523
    // Paths are REPO-ROOT relative, as the deterministic branch's `Out' is. They were previously
    // "outcomes/..." with no analyses/transport_dvd/ prefix, which resolves to <repo>/outcomes -
    // a directory that does not exist, so every write failed with "could not be opened for output".
    // This is the first-run check flagged in run.do; it is now closed.
    local Out "analyses/transport_dvd/outcomes"
    cap mkdir "`Out'/A_trial/bootstrap"
    cap mkdir "`Out'/B_transport/bootstrap"
    forval b = $min_bs / $max_bs {
        set seed `=`SEEDBASE' + `b''
        run_one ///
            "${data_path}/bootstrap/MRDR Long MI B`b'.dta" ///
            1 ///
            "`Out'/A_trial/bootstrap/bcr_vd_l2_B`b'.mmat" ///
            "`Out'/A_trial/bootstrap/bcr_dvd_l2_B`b'.mmat" ///
            "`Out'/B_transport/bootstrap/transport_dvd_B`b'.mmat"
    }
}
