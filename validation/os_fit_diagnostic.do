**********
* Monash Myeloma Model - OS parametric fit diagnostic
*
* Purpose: Investigate why simulated overall survival sits below the MRDR benchmark
*          (worse at 5y than 3y, worse in poorer responders -- the classic signature
*          of a parametric distribution under-fitting the survival tail).
*
* Method:  IN-SAMPLE fit-quality check, independent of the simulation. For each line/BCR
*          stratum it compares the registry Kaplan-Meier survival against several parametric
*          families fit to the SAME data, at 1/3/5/10-year horizons, plus AIC/BIC. The engine
*          currently samples OS from a Weibull ($dOS = "w" in prep/risk_equations.do); if a
*          more flexible family (log-normal / log-logistic / generalised gamma) reproduces the
*          KM tail markedly better, the OS validation gap is a distributional-form issue and the
*          remedy is to change the OS family rather than to debug the engine.
*
* Data:    Restricted MRDR (uses $data_dir from config.do) -- run from the repository root.
* Output:  A per-stratum, per-BCR table (KM vs each distribution at each horizon) + AIC.
*
* Notes:   Fits are intercept-only WITHIN each BCR group, so the parametric survival curve is a
*          single curve directly comparable to that group's KM (this matches how validation
*          aggregates by BCR). A per-group Weibull here is actually MORE flexible than the engine
*          (which shares one shape parameter across BCR), so any tail under-fit seen here is a
*          lower bound on the engine's.
**********

clear all
set more off

if "$repo_root" != "" cd "$repo_root"   // cd to repo root only if config.do set it
capture run "config.do"                 // machine-specific paths ($data_dir; git-ignored)

* Save a plain-text copy of the output so it can be reviewed after the run
* (interactive runs otherwise leave output only in the Results window).
capture log close _all
log using "validation/os_fit_diagnostic_output.txt", text replace

use "${data_dir}/MRDR Long MI.dta", clear
mi extract 1, clear

* Horizons (months) and their year labels -- kept in lockstep
global HMON  12 36 60 120
global HYR   1 3 5 10

* Distributions to compare (Stata streg names). First is the engine's current family.
global DISTS weibull lognormal loglogistic ggamma gompertz exponential


**********
* Reusable routine: one stratum (line), comparing KM vs each parametric family by BCR
**********

capture program drop os_diag
program define os_diag
	* args: origin_event  bcr_var  extra_if   label
	args origin bcrvar extra label

	di _n(2) "{hline 92}"
	di "STRATUM: `label'   (origin = Event1==`origin'; BCR = `bcrvar'; filter: `extra')"
	di "{hline 92}"

	forvalues g = 1/6 {
		preserve

		* Restrict to this BCR group within the stratum
		quietly keep if `bcrvar' == `g' `extra'
		quietly count
		if r(N) == 0 {
			restore
			continue
		}

		* Time origin = line start; failure = death (matches generate_benchmarks.do OS setup)
		quietly stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == `origin') scale(30.4375)
		quietly count if _d == 1
		local nevt = r(N)
		quietly count if !missing(_t)
		local nobs = r(N)
		if `nevt' < 10 {
			di _n "  BCR=`g': only `nevt' events (N=`nobs') -- skipped (too few to fit)."
			restore
			continue
		}

		* ---- Kaplan-Meier survival at each horizon (= last KM step at or before h) ----
		quietly sts generate km_s = s
		foreach h of global HMON {
			quietly summarize km_s if _t <= `h', meanonly
			if r(N) > 0  local km_`h' = r(min)
			else         local km_`h' = .
		}
		drop km_s

		* ---- Parametric fits: predicted survival at each horizon + AIC ----
		* Predicted S(h) is read off by setting _t=h. We also set the entry time _t0=0 so that
		* predict,surv returns the UNCONDITIONAL S(h); the long format has delayed entry, and
		* without this predict returns the late-entry-conditional S(h)/S(_t0), which exceeds 1
		* for records that entered after h. (_t / _t0 are stset system vars -- stash & restore.)
		quietly gen double tt_orig  = _t
		quietly gen double tt0_orig = _t0
		foreach d of global DISTS {
			capture quietly streg, distribution(`d')
			if _rc {
				foreach h of global HMON {
					local p_`d'_`h' = .
				}
				local aic_`d' = .
			}
			else {
				* Intercept-only -> S(h) is constant across obs; read it off by setting _t=h, _t0=0
				foreach h of global HMON {
					quietly replace _t  = `h'
					quietly replace _t0 = 0
					capture quietly predict double pred_s, surv
					if _rc == 0 {
						quietly summarize pred_s, meanonly
						local p_`d'_`h' = r(mean)
						drop pred_s
					}
					else {
						local p_`d'_`h' = .
					}
				}
				quietly replace _t  = tt_orig      // restore stset times before next fit
				quietly replace _t0 = tt0_orig
				quietly estat ic
				matrix ic_S = r(S)
				local aic_`d' = ic_S[1,5]
			}
		}
		drop tt_orig tt0_orig

		* ---- Print table for this BCR group ----
		di _n "  BCR=`g'   (N=`nobs', events=`nevt')"
		di "  " %-7s "Horizon" "  " %7s "KM" _continue
		foreach d of global DISTS {
			di "  " %10s abbrev("`d'",10) _continue
		}
		di ""
		local i = 1
		foreach h of global HMON {
			local yr : word `i' of $HYR
			di "  " %-7s "`yr'-yr" "  " %7.3f `km_`h'' _continue
			foreach d of global DISTS {
				di "  " %10.3f `p_`d'_`h'' _continue
			}
			di ""
			local ++i
		}
		di "  " %-7s "AIC" "  " %7s "-" _continue
		foreach d of global DISTS {
			di "  " %10.0f `aic_`d'' _continue
		}
		di ""

		restore
	}
end


**********
* Run the diagnostic on the strata that drove the OS validation failures
**********

* L1, NoASCT  (origin = L1 start = Event1==10; BCR_L1; transplant-ineligible arm)
os_diag 10 BCR_L1 "& SCT == 0" "L1 (No ASCT)"

* L2  (origin = L2 start = Event1==20; BCR_L2)
os_diag 20 BCR_L2 "" "L2"

* L3  (origin = L3 start = Event1==30; BCR_L3)
os_diag 30 BCR_L3 "" "L3"

di _n(2) "{hline 92}"
di "Reading the table: where the Weibull column sits BELOW the KM column at 5-yr/10-yr while"
di "log-normal / log-logistic / generalised-gamma track KM, the OS shortfall is a tail-fit"
di "problem in the Weibull family. Compare AIC across families (lower = better) to pick the"
di "replacement for \$dOS in prep/risk_equations.do."
di "{hline 92}"


**********
* PART B: Is the OS shortfall better resolved by a FLEXIBLE family or by SEPARATE-BY-LINE fits?
*
* On the engine's own OS stset (origin = diagnosis, every line/state spell stacked on one shared
* time scale), compare the TOTAL AIC of:
*   - POOLED, one shared shape across states (engine analogue: covariates + i.OS##i.BCR), and
*   - SEPARATE-BY-LINE, the model refit within each OS state (own shape + own covariate effects),
*     the per-state AICs summed,
* each for Weibull, generalised gamma (flexible, official Stata) and Royston-Parmar splines
* (stpm2, if installed). Lower total AIC is better (AIC penalises the extra parameters that
* separate fits and flexible families spend; it is additive across the disjoint per-line samples).
* BIC is omitted here because its sample-size penalty does not sum cleanly across separate fits.
**********

di _n(2) "{hline 92}"
di "PART B: POOLED vs SEPARATE-BY-LINE, and Weibull vs flexible families (AIC tournament)"
di "{hline 92}"

* Engine OS stset (origin = diagnosis; id() handles the stacked per-state spells / delayed entry)
quietly stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 3) scale(30.4375)

* Common estimation sample so POOLED and SEPARATE use exactly the same records: complete
* covariates/response, a valid spell, and only OS states with enough deaths to fit (>=20).
capture drop esamp inlvl
quietly gen byte esamp = !missing(Age, Age2, Male, ECOGcc, RISS, BCR, OS) & _st == 1
quietly levelsof OS if esamp, local(allos)
local uselevels ""
foreach L of local allos {
	quietly count if esamp & OS == `L' & _d == 1
	if r(N) >= 20 local uselevels "`uselevels' `L'"
}
quietly gen byte inlvl = 0
foreach L of local uselevels {
	quietly replace inlvl = 1 if OS == `L'
}
quietly replace esamp = esamp & inlvl
quietly count if esamp
di "Estimation sample: " r(N) " state-spells over OS levels:`uselevels'"

tempname R

* ---- Pooled (one model, shared shape) ----
* Families here are the fast, stable ones: Weibull (the engine's choice) and log-normal (a
* heavier-tailed flexible alternative -- Part A flagged it as the main competitor in later lines).
* Generalised gamma is NOT fit on the pooled saturated i.OS##i.BCR data: with 3 ancillary
* parameters it is slow/unstable to converge there (this is what hangs). It is already assessed
* per-cell in Part A. iterate() caps run time so a non-converging fit cannot hang.
foreach spec in weibull lognormal {
	capture quietly streg Age Age2 Male i.ECOGcc i.RISS i.OS##i.BCR if esamp, distribution(`spec') iterate(50)
	if _rc == 0 {
		quietly estat ic
		matrix `R' = r(S)
		local aic_pool_`spec' = `R'[1,5]
		local df_pool_`spec'  = `R'[1,4]
	}
	else {
		local aic_pool_`spec' = .
		local df_pool_`spec'  = .
	}
}

* ---- Separate by line (own shape + own covariate effects per state) ----
foreach spec in weibull lognormal {
	local saic = 0
	local sdf  = 0
	local nmod = 0
	local ok   = 1
	foreach L of local uselevels {
		capture quietly streg Age Age2 Male i.ECOGcc i.RISS i.BCR if esamp & OS == `L', distribution(`spec') iterate(50)
		if _rc == 0 {
			quietly estat ic
			matrix `R' = r(S)
			local saic = `saic' + `R'[1,5]
			local sdf  = `sdf'  + `R'[1,4]
			local ++nmod
		}
		else local ok = 0
	}
	local aic_sep_`spec' = cond(`ok', `saic', .)
	local df_sep_`spec'  = cond(`ok', `sdf', .)
	local n_sep_`spec'   = `nmod'
}

* ---- Royston-Parmar flexible parametric (stpm2), if available ----
local have_stpm2 = 0
capture which stpm2
if _rc == 0 local have_stpm2 = 1

if `have_stpm2' {
	* pooled
	capture quietly stpm2 Age Age2 Male i.ECOGcc i.RISS i.OS##i.BCR if esamp, df(3) scale(hazard) iterate(50)
	if _rc == 0 {
		quietly estat ic
		matrix `R' = r(S)
		local aic_pool_rp = `R'[1,5]
		local df_pool_rp  = `R'[1,4]
	}
	else {
		local aic_pool_rp = .
		local df_pool_rp  = .
	}
	* separate
	local saic = 0
	local sdf  = 0
	local ok   = 1
	foreach L of local uselevels {
		capture quietly stpm2 Age Age2 Male i.ECOGcc i.RISS i.BCR if esamp & OS == `L', df(3) scale(hazard) iterate(50)
		if _rc == 0 {
			quietly estat ic
			matrix `R' = r(S)
			local saic = `saic' + `R'[1,5]
			local sdf  = `sdf'  + `R'[1,4]
		}
		else local ok = 0
	}
	local aic_sep_rp = cond(`ok', `saic', .)
	local df_sep_rp  = cond(`ok', `sdf', .)
}

* ---- SEPARATE BY LINE, each fitted from its OWN line start (line-start clock) ----
* This is the deployable per-line model and matches the scale OS is benchmarked/validated on:
* origin = each line's start (cf. the blocks above, which use the engine's diagnosis clock).
* Because the time origin differs, these AICs are comparable ACROSS these families but NOT to the
* diagnosis-clock rows -- cross-check the two clocks via Part A's predicted-vs-KM (line-start scale).
* Lines and their start events (Event1): L1 NoASCT=10, L1 ASCT=11, L2=20, L3=30, L4=40; each line
* uses its own response variable.
* NB: one local per line -- ';' is not a command separator in default Stata (#delimit cr).
local nlines = 5
local o1 "10"
local f1 "& SCT == 0"
local b1 "BCR_L1"
local o2 "11"
local f2 "& SCT == 1"
local b2 "BCR_SCT"
local o3 "20"
local f3 ""
local b3 "BCR_L2"
local o4 "30"
local f4 ""
local b4 "BCR_L3"
local o5 "40"
local f5 ""
local b5 "BCR_L4"

foreach spec in weibull lognormal {
	local saic = 0
	local sdf  = 0
	local nmod = 0
	forvalues k = 1/`nlines' {
		capture quietly stset Date1 if(F_OS != 1 `f`k''), id(ID_BS) failure(Event1 == 104) origin(Event1 == `o`k'') scale(30.4375)
		if _rc continue
		quietly count if _d == 1
		if r(N) < 20 continue
		capture quietly streg Age Age2 Male i.ECOGcc i.RISS i.`b`k'', distribution(`spec') iterate(50)
		if _rc continue
		quietly estat ic
		matrix `R' = r(S)
		local saic = `saic' + `R'[1,5]
		local sdf  = `sdf'  + `R'[1,4]
		local ++nmod
	}
	local aic_ls_`spec' = `saic'
	local df_ls_`spec'  = `sdf'
	local n_ls_`spec'   = `nmod'
}

if `have_stpm2' {
	local saic = 0
	local sdf  = 0
	local nmod = 0
	forvalues k = 1/`nlines' {
		capture quietly stset Date1 if(F_OS != 1 `f`k''), id(ID_BS) failure(Event1 == 104) origin(Event1 == `o`k'') scale(30.4375)
		if _rc continue
		quietly count if _d == 1
		if r(N) < 20 continue
		capture quietly stpm2 Age Age2 Male i.ECOGcc i.RISS i.`b`k'', df(3) scale(hazard) iterate(50)
		if _rc continue
		quietly estat ic
		matrix `R' = r(S)
		local saic = `saic' + `R'[1,5]
		local sdf  = `sdf'  + `R'[1,4]
		local ++nmod
	}
	local aic_ls_rp = `saic'
	local df_ls_rp  = `sdf'
	local n_ls_rp   = `nmod'
}

* ---- Report ----
di _n "  GROUP 1 -- DIAGNOSIS clock (the engine's current time scale). AIC comparable WITHIN this group;"
di    "  this group answers: within the engine's clock, do per-line shapes beat one pooled model?"
di    "  Specification                          df    total AIC"
di    "  ------------------------------------  -----  ----------"
foreach spec in weibull lognormal {
	di "  Pooled (engine)       -- " %-11s "`spec'" "  " %5.0f `df_pool_`spec'' "  " %10.0f `aic_pool_`spec''
}
if `have_stpm2' di "  Pooled (engine)       -- RP splines    " %5.0f `df_pool_rp' "  " %10.0f `aic_pool_rp'
foreach spec in weibull lognormal {
	di "  Separate by line      -- " %-11s "`spec'" "  " %5.0f `df_sep_`spec'' "  " %10.0f `aic_sep_`spec'' "   (`n_sep_`spec'' models)"
}
if `have_stpm2' di "  Separate by line      -- RP splines    " %5.0f `df_sep_rp' "  " %10.0f `aic_sep_rp'

di _n "  GROUP 2 -- LINE-START clock (each line clocked from its OWN start = the benchmark / validation"
di    "  scale, and the deployable per-line model). AIC comparable WITHIN this group; this group"
di    "  answers: which family for the per-line, line-start OS model?"
di    "  Specification                          df    total AIC"
di    "  ------------------------------------  -----  ----------"
foreach spec in weibull lognormal {
	di "  Separate by line      -- " %-11s "`spec'" "  " %5.0f `df_ls_`spec'' "  " %10.0f `aic_ls_`spec'' "   (`n_ls_`spec'' models)"
}
if `have_stpm2' di "  Separate by line      -- RP splines    " %5.0f `df_ls_rp' "  " %10.0f `aic_ls_rp' "   (`n_ls_rp' models)"
if !`have_stpm2' di _n "  (Royston-Parmar rows skipped: stpm2 not installed -- 'ssc install stpm2' to enable.)"

di _n "  How to read it:"
di    "  - GROUP 1, within the diagnosis clock: Pooled vs Separate shows whether per-line shapes help."
di    "  - GROUP 2: the family choice for the per-line model fitted on the validation (line-start) scale."
di    "  - GROUP 1 and GROUP 2 use DIFFERENT time origins, so their AICs are NOT comparable to each"
di    "    other. To judge the line-start model against the registry, use Part A's predicted-vs-KM"
di    "    (also line-start): does the chosen family track KM at 5-10y by line x BCR?"
di "{hline 92}"

capture drop esamp inlvl

capture log close
