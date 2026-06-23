**********
* Monash Myeloma Model - Validate Simulation
* 
* Purpose: Compare simulation results to MRDR benchmarks
*
* Author: Adam Irving
* Date: January 2026
**********

cd "~/Documents/Monash/Research/Blood Disorders/EpiMAP-Local/Myeloma/Simulation"

**********
// Load Benchmarks
**********

di _n "{hline 80}"
di "LOADING BENCHMARKS"
di "{hline 80}"

// OS benchmarks
import delimited "validation/benchmarks/OS_L1_NoASCT.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L1_NoASCT_bench)

import delimited "validation/benchmarks/OS_ASCT.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_ASCT_bench)

import delimited "validation/benchmarks/OS_L2.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L2_bench)

import delimited "validation/benchmarks/OS_L3.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L3_bench)

// BCR distributions
import delimited "validation/benchmarks/BCR.csv", clear case(preserve)
mkmat N CR VG PR MR SD PD, matrix(BCR_bench)

// TXD benchmarks
import delimited "validation/benchmarks/TXD_L1_NoASCT.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TXD_L1_NoASCT_bench)

import delimited "validation/benchmarks/TXD_L1_ASCT.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TXD_L1_ASCT_bench)

import delimited "validation/benchmarks/TXD_L2.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TXD_L2_bench)

// TFI benchmarks
import delimited "validation/benchmarks/TFI_L1_NoASCT.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TFI_L1_NoASCT_bench)

import delimited "validation/benchmarks/TFI_L1_ASCT.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TFI_L1_ASCT_bench)

import delimited "validation/benchmarks/TFI_L2.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TFI_L2_bench)

import delimited "validation/benchmarks/TFI_L3.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TFI_L3_bench)

// Pathways
import delimited "validation/benchmarks/Pathways.csv", clear case(preserve)
mkmat N ASCT L2 L3 L4 L5 L6 L7 L8, matrix(Pathways_bench)

di "Benchmarks loaded successfully"

**********
// Load Simulation
**********

use "analyses/base_model/simulated/all_0_population_1_101212.dta", clear

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
	qui summarize surv_temp if BCR_L1 == `bcr' & SCT_L1 == 0 & abs(_t - 36) < 0.5
	if r(N) > 0 {
		local sim_3yr = r(mean) * 100
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
	qui summarize surv_temp if BCR_L1 == `bcr' & SCT_L1 == 0 & abs(_t - 60) < 0.5
	qui if r(N) > 0 {
		local sim_5yr = r(mean) * 100
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
	qui summarize surv_temp if BCR_SCT == `bcr' & abs(_t - 36) < 0.5
	if r(N) > 0 {
		local sim_3yr = r(mean) * 100
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
	qui summarize surv_temp if BCR_SCT == `bcr' & abs(_t - 60) < 0.5
	qui if r(N) > 0 {
		local sim_5yr = r(mean) * 100
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
	qui summarize surv_temp if BCR_L2 == `bcr' & abs(_t - 36) < 0.5
	if r(N) > 0 {
		local sim_3yr = r(mean) * 100
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
	qui summarize surv_temp if BCR_L2 == `bcr' & abs(_t - 60) < 0.5
	qui if r(N) > 0 {
		local sim_5yr = r(mean) * 100
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
	qui summarize surv_temp if BCR_L3 == `bcr' & abs(_t - 36) < 0.5
	if r(N) > 0 {
		local sim_3yr = r(mean) * 100
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
	qui summarize surv_temp if BCR_L3 == `bcr' & abs(_t - 60) < 0.5
	qui if r(N) > 0 {
		local sim_5yr = r(mean) * 100
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
// TXD
**********

local tolerance = 0.20

// TXD_L1 NoASCT
qui stset TXD_L1, failure(MOR_L1S==0) id(ID)

* BCR | Benchmark | Simulated | Ratio | Pass?

qui forvalues bcr = 1/6 {
	local bench_median = TXD_L1_NoASCT_bench[`bcr', 3]
	
	capture qui stsum if BCR_L1 == `bcr' & SCT_L1 == 0
	if _rc == 0 & r(N_sub) > 0 {
		local sim_median = r(p50)
		local ratio = `sim_median' / `bench_median'
		
		if abs(`ratio' - 1) <= `tolerance' {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | " %9.1f `bench_median' " | " %9.1f `sim_median' " | " %5.2f `ratio' " | `status'"
	}
}

// TXD_L1 ASCT
* BCR | Benchmark | Simulated | Ratio | Pass?

qui forvalues bcr = 1/4 {
	local bench_median = TXD_L1_ASCT_bench[`bcr', 3]
	
	capture qui stsum if BCR_SCT == `bcr' & SCT_L1 == 1
	if _rc == 0 & r(N_sub) > 0 {
		local sim_median = r(p50)
		local ratio = `sim_median' / `bench_median'
		
		if abs(`ratio' - 1) <= `tolerance' {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | " %9.1f `bench_median' " | " %9.1f `sim_median' " | " %5.2f `ratio' " | `status'"
	}
}

**********
// TFI
**********

di _n _n "{hline 80}"
di "TEST 2: TFI PREDICTIONS"
di "{hline 80}"

local tolerance = 0.20

// TFI_L1_NoASCT
qui stset TFI_L1, failure(MOR_L1E==0) id(ID)

di _n "TFI_L1 by BCR_L1 (NoASCT):"
di "BCR | Benchmark | Simulated | Ratio | Pass?"
di "----|-----------|-----------|-------|-------"

qui forvalues bcr = 1/6 {
	local bench_median = TFI_L1_NoASCT_bench[`bcr', 3]
	
	capture qui stsum if BCR_L1 == `bcr' & SCT_L1 == 0
	if _rc == 0 & r(N_sub) > 0 {
		local sim_median = r(p50)
		local ratio = `sim_median' / `bench_median'
	
		if abs(`ratio' - 1) <= `tolerance' {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "BCR=`bcr' | " %9.1f `bench_median' " | " %9.1f `sim_median' " | " %5.2f `ratio' " | `status'"
	}
}

// TFI_L1_ASCT
di _n "TFI_L1 by BCR_SCT (ASCT):"
di "BCR | Benchmark | Simulated | Ratio | Pass?"
di "----|-----------|-----------|-------|-------"

qui forvalues bcr = 1/4 {
	local bench_median = TFI_L1_ASCT_bench[`bcr', 3]
	
	capture qui stsum if BCR_SCT == `bcr' & SCT_L1 == 1
	if _rc == 0 & r(N_sub) > 0 {
		local sim_median = r(p50)
		local ratio = `sim_median' / `bench_median'
		
		if abs(`ratio' - 1) <= `tolerance' {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "BCR=`bcr' | " %9.1f `bench_median' " | " %9.1f `sim_median' " | " %5.2f `ratio' " | `status'"
	}
}

// TFI_L2
qui stset TFI_L2, failure(MOR_L2E==0) id(ID)

di _n "TFI_L2 by BCR_L2:"
di "BCR | Benchmark | Simulated | Ratio | Pass?"
di "----|-----------|-----------|-------|-------"

qui forvalues bcr = 1/6 {
	local bench_median = TFI_L2_bench[`bcr', 3]
	
	capture qui stsum if BCR_L2 == `bcr'
	if _rc == 0 & r(N_sub) > 0 {
		local sim_median = r(p50)
		local ratio = `sim_median' / `bench_median'
		
		if abs(`ratio' - 1) <= `tolerance' {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "BCR=`bcr' | " %9.1f `bench_median' " | " %9.1f `sim_median' " | " %5.2f `ratio' " | `status'"
	}
}

// Test ordering preserved
di _n "Checking TFI_L2 ordering (CR > VG > PR > MR > SD)..."
qui stsum if BCR_L2 == 1
local med_cr = r(p50)
qui stsum if BCR_L2 == 2
local med_vg = r(p50)
qui stsum if BCR_L2 == 3
local med_pr = r(p50)
qui stsum if BCR_L2 == 4
local med_mr = r(p50)
qui stsum if BCR_L2 == 5
local med_sd = r(p50)

if (`med_cr' > `med_vg') & (`med_vg' > `med_pr') & (`med_pr' > `med_mr') & (`med_mr' > `med_sd') {
	local status "PASS"
	local tests_passed = `tests_passed' + 1
}
else {
	local status "FAIL"
	local tests_failed = `tests_failed' + 1
}

local tests_run = `tests_run' + 1

di "  CR=" %5.1f `med_cr' " > VG=" %5.1f `med_vg' " > PR=" %5.1f `med_pr' " > MR=" %5.1f `med_mr' " > SD=" %5.1f `med_sd'
di "  Result: `status'"


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

// ASCT
qui count if SCT_L1 == 1
local pct_sim = r(N) / `n_total' * 100
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
