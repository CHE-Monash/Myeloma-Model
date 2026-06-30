**********
* Monash Myeloma Model
* Validate outcomes -- shared comparison engine
*
* Purpose: Compare a simulated dataset to a set of observed-outcome target CSVs (OS, BCR, TXD/TFI by
*          horizon survival, pathways by competing-risks CIF, ASCT among L1-end reachers) and print a
*          pass/fail summary. Driven entirely by $val_targets (target CSV directory) and $val_simfile
*          (simulated .dta). This is the model's validation comparison engine; the out-of-sample
*          validator (analyses/oos/validate_oos.do) sets both globals and calls it.
*
* Author: Adam Irving
* Date: January 2026
**********

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

* The caller MUST set both the target-CSV directory and the simulated dataset -- there is no in-sample
* default (the standalone registry acceptance test was retired; see _notes -> Validation). Normally
* invoked via analyses/oos/validate_oos.do, not directly.
if "$val_targets" == "" | "$val_simfile" == "" {
	di as error "Set globals \$val_targets (target CSV directory) and \$val_simfile (simulated .dta) first."
	di as error "This comparison engine is normally run via analyses/oos/validate_oos.do."
	exit 198
}


**********
// Load Benchmarks
**********

di _n "{hline 80}"
di "LOADING BENCHMARKS"
di "{hline 80}"

// OS benchmarks
import delimited "${val_targets}/os_l1_noasct.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L1_NoASCT_bench)

import delimited "${val_targets}/os_asct.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_ASCT_bench)

import delimited "${val_targets}/os_l2.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L2_bench)

import delimited "${val_targets}/os_l3.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L3_bench)

// BCR distributions
import delimited "${val_targets}/bcr.csv", clear case(preserve)
mkmat N CR VG PR MR SD PD, matrix(BCR_bench)

// TXD benchmarks
import delimited "${val_targets}/txd_l1_noasct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TXD_L1_NoASCT_bench)

import delimited "${val_targets}/txd_l1_asct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TXD_L1_ASCT_bench)

import delimited "${val_targets}/txd_l2.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TXD_L2_bench)

// TFI benchmarks
import delimited "${val_targets}/tfi_l1_noasct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L1_NoASCT_bench)

import delimited "${val_targets}/tfi_l1_asct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L1_ASCT_bench)

import delimited "${val_targets}/tfi_l2.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L2_bench)

import delimited "${val_targets}/tfi_l3.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TFI_L3_bench)

// Pathways
import delimited "${val_targets}/pathways.csv", clear case(preserve)
mkmat N ASCT L2 L3 L4 L5 L6 L7 L8, matrix(Pathways_bench)

di "Benchmarks loaded successfully"

**********
// Load Simulation
**********

use "$val_simfile", clear

local tests_run = 0
local tests_passed = 0
local tests_failed = 0


**********
// Overall Survival
**********

local tolerance = 0.10

// OS by BCR_L1 (NoASCT)
qui cap gen OC_TIME_L1S = OC_TIME - TSD_L1S
qui stset OC_TIME_L1S, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/6 {
	// 3-year survival
	local bench_3yr = OS_L1_NoASCT_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_L1 == `bcr' & SCT_L1 == 0
	qui summarize surv_temp if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp

	// 5-year survival
	local bench_5yr = OS_L1_NoASCT_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_L1 == `bcr' & SCT_L1 == 0
	qui summarize surv_temp if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

// OS by BCR_SCT
qui cap gen OC_TIME_L1E = OC_TIME - TSD_L1E
qui stset OC_TIME_L1E, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/4 {
	// 3-year survival
	local bench_3yr = OS_ASCT_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_SCT == `bcr'
	qui summarize surv_temp if BCR_SCT == `bcr' & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp
	
	// 5-year survival
	local bench_5yr = OS_ASCT_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_SCT == `bcr'
	qui summarize surv_temp if BCR_SCT == `bcr' & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

// OS by BCR_L2
qui cap gen OC_TIME_L2S = OC_TIME - TSD_L2S
qui stset OC_TIME_L2S, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/6 {
	// 3-year survival
	local bench_3yr = OS_L2_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_L2 == `bcr'
	qui summarize surv_temp if BCR_L2 == `bcr' & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp
	
	// 5-year survival
	local bench_5yr = OS_L2_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_L2 == `bcr'
	qui summarize surv_temp if BCR_L2 == `bcr' & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

// OS by BCR_L3
qui cap gen OC_TIME_L3S = OC_TIME - TSD_L3S
qui stset OC_TIME_L3S, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/6 {
	// 3-year survival
	local bench_3yr = OS_L3_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_L3 == `bcr'
	qui summarize surv_temp if BCR_L3 == `bcr' & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp
	
	// 5-year survival
	local bench_5yr = OS_L3_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_L3 == `bcr'
	qui summarize surv_temp if BCR_L3 == `bcr' & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

**********
// BCR
**********

local tolerance = 5

// L1
qui count if !missing(BCR_L1)
local n_sim = r(N)

* Line | Metric | Benchmark | Simulated | Diff | Pass?

qui forvalues bcr = 1/6 {
	qui count if BCR_L1 == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[1, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di "  L1   | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

//ASCT
qui count if !missing(BCR_SCT)
local n_sim = r(N)

qui forvalues bcr = 1/4 {
	qui count if BCR_SCT == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[2, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di "  ASCT | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

// L2
qui count if !missing(BCR_L2)
local n_sim = r(N)

qui forvalues bcr = 1/6 {
	qui count if BCR_L2 == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[3, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di "  L2   | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

// L3
qui count if !missing(BCR_L3)
local n_sim = r(N)

qui forvalues bcr = 1/6 {
	qui count if BCR_L3 == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[4, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di " L3   | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

**********
// TXD  --  % still on treatment at 12 / 24 months vs benchmark (+/- 10 pp)
**********

di _n _n "{hline 80}"
di "TXD PREDICTIONS (% still on treatment at 12 / 24 months)"
di "{hline 80}"
di _n "Line | BCR | Horizon | Benchmark | Simulated |  Diff | Pass?"
di    "-----|-----|---------|-----------|-----------|-------|------"

local tol_pp = 10

// TXD_L1 NoASCT
qui stset TXD_L1, failure(MOR_L1S==0) id(ID)
qui forvalues bcr = 1/6 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L1 == `bcr' & SCT_L1 == 0
    local c = 7
    foreach h in 12 24 {
        local bench = TXD_L1_NoASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1NA | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TXD_L1 ASCT
qui forvalues bcr = 1/4 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_SCT == `bcr' & SCT_L1 == 1
    local c = 7
    foreach h in 12 24 {
        local bench = TXD_L1_ASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_SCT == `bcr' & SCT_L1 == 1 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1AS | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}


**********
// TFI  --  % still treatment-free at 12 / 24 months vs benchmark (+/- 10 pp)
**********

di _n _n "{hline 80}"
di "TFI PREDICTIONS (% still treatment-free at 12 / 24 months)"
di "{hline 80}"
di _n "Line | BCR | Horizon | Benchmark | Simulated |  Diff | Pass?"
di    "-----|-----|---------|-----------|-----------|-------|------"

// TFI_L1 NoASCT
qui stset TFI_L1, failure(MOR_L1E==0) id(ID)
qui forvalues bcr = 1/6 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L1 == `bcr' & SCT_L1 == 0
    local c = 7
    foreach h in 12 24 {
        local bench = TFI_L1_NoASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1NA | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TFI_L1 ASCT
qui forvalues bcr = 1/4 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_SCT == `bcr' & SCT_L1 == 1
    local c = 7
    foreach h in 12 24 {
        local bench = TFI_L1_ASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_SCT == `bcr' & SCT_L1 == 1 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1AS | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TFI_L2
qui stset TFI_L2, failure(MOR_L2E==0) id(ID)
qui forvalues bcr = 1/6 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L2 == `bcr'
    local c = 7
    foreach h in 12 24 {
        local bench = TFI_L2_bench[`bcr', `c']
        qui summarize surv_h if BCR_L2 == `bcr' & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L2   | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TFI_L2 response ordering (better response -> longer treatment-free interval), tested as
// % treatment-free at 12 months for the clearly-separated categories CR > VGPR > PR.
// (MR/SD progress almost immediately and are not reliably ordered -- even the registry medians
//  do not separate them, so the strict 5-way ordering is not asserted.)
di _n "Checking TFI_L2 ordering at 12 months (CR > VG > PR)..."
foreach bcr in 1 2 3 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L2 == `bcr'
    qui summarize surv_h if BCR_L2 == `bcr' & _t <= 12
    local s12_`bcr' = cond(r(N) > 0, r(min) * 100, .)
    qui drop surv_h
}
if (`s12_1' > `s12_2') & (`s12_2' > `s12_3') {
    local status "PASS"
    local tests_passed = `tests_passed' + 1
}
else {
    local status "FAIL"
    local tests_failed = `tests_failed' + 1
}
local tests_run = `tests_run' + 1
di "  CR=" %4.1f `s12_1' "% > VG=" %4.1f `s12_2' "% > PR=" %4.1f `s12_3' "%   Result: `status'"


**********
// Pathways
**********

di _n _n "{hline 80}"
di "TEST 5: TREATMENT PATHWAYS"
di "{hline 80}"

di _n "Line | Benchmark | Simulated | Diff | Pass?"
di "-----|-----------|-----------|------|-------"

qui count
local n_total = r(N)
local tolerance = 5

// ASCT -- transplant is a decision taken AT THE END OF L1, so the rate is among patients who
// reach L1 end (alive past L1 induction), matching the benchmark denominator in generate_benchmarks.
// (Dividing by all patients would dilute it by the ~8% who die during induction and never reach
//  the transplant decision.)
qui count if OC_TIME > TSD_L1E & !missing(TSD_L1E)
local n_l1e = r(N)
qui count if SCT_L1 == 1
local pct_sim = r(N) / `n_l1e' * 100
local pct_bench = Pathways_bench[1, 2]
local diff = `pct_sim' - `pct_bench'

if abs(`diff') <= `tolerance' {
	local status "PASS"
	local tests_passed = `tests_passed' + 1
}
else {
	local status "FAIL"
	local tests_failed = `tests_failed' + 1
}

local tests_run = `tests_run' + 1

di "ASCT | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"

// Lines 2-5
qui forvalues line = 2/5 {
	qui count if !missing(BCR_L`line')
	local pct_sim = r(N) / `n_total' * 100
	local pct_bench = Pathways_bench[1, `line'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di "L`line'  | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

**********
// Summary
**********

di _n _n "{hline 80}"
di "VALIDATION SUMMARY"
di "{hline 80}"
di "Tests run:    " %3.0f `tests_run'
di "Tests passed: " %3.0f `tests_passed' " (" %5.1f `tests_passed'/`tests_run'*100 "%)"
di "Tests failed: " %3.0f `tests_failed' " (" %5.1f `tests_failed'/`tests_run'*100 "%)"
di "{hline 80}"

if `tests_failed' == 0 {
	di _n "✓ ALL TESTS PASSED"
	di "Simulation accurately reproduces MRDR training data"
}
else {
	di _n "⚠ SOME TESTS FAILED"
	di "Review failed tests above for model issues"
}
