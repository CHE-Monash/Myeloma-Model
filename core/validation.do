**********
* EpiMAP Myeloma - Validation
*
* Purpose: Run systematic validation checks after simulation completes
* Usage: Run after simulation_engine.do
**********

qui {
	n di as text _n "====================="
	n di as text "Simulation Validation"
	n di as text "====================="

	scalar validation_errors = 0
	scalar validation_warnings = 0

	// 1. Mortality

	n di as text _n "--- 1. Mortality ---"

	// 1.1 Deaths should be recorded only once

		mata {
			death_inconsistencies = 0
			for (i = 1; i <= rows(mMOR); i++) {
				found_death = 0
				for (j = 1; j <= cols(mMOR); j++) {
					if (mMOR[i,j] == 1) found_death = 1
					if (found_death == 1 & mMOR[i,j] == 0) {
						death_inconsistencies++
						break
					}
				}
			}
			st_numscalar("temp_result", death_inconsistencies)
		}
		if (temp_result > 0) {
			n di as error "ERROR: `=temp_result' patients have death flag reversal"
			scalar validation_errors = validation_errors + 1
		}
		else {
			n di as text "PASS: No death flag reversals"
		}


	// 1.2 Every death should have a recorded time
	mata {
		deaths_total = sum(mOC[.,2] :== 1)
		deaths_with_time = sum((mOC[.,2] :== 1) :& (mOC[.,1] :< .))
		st_numscalar("deaths_total", deaths_total)
		st_numscalar("deaths_without_time", deaths_total - deaths_with_time)
	}
	if (deaths_without_time > 0) {
		n di as error "ERROR: `=deaths_without_time' deaths have missing OC_TIME"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All `=deaths_total' deaths have recorded OC_TIME"
	}

	// 1.3 Death time should be positive
	mata {
		death_times = select(mOC[.,1], mOC[.,2] :== 1)
		if (rows(death_times) > 0) {
			st_numscalar("negative_deaths", sum(death_times :< 0))
		}
		else {
			st_numscalar("negative_deaths", 0)
		}
	}
	if (negative_deaths > 0) {
		n di as error "ERROR: `=negative_deaths' deaths have negative OC_TIME"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All death times are non-negative"
	}

	// 2. Time
	n di as text _n "--- 2. Time ---"

	// 2.1 mTSD monotonicity
	mata {
		tsd_decreases = 0
		for (omc = 2; omc <= cols(mTSD); omc++) {
			if (omc > 1) {
				alive_at_omc = (mMOR[., omc-1] :== 0)
			}
			else {
				alive_at_omc = J(rows(mTSD), 1, 1)
			}
			tsd_diff = mTSD[., omc] - mTSD[., omc-1]
			tsd_decreases = tsd_decreases + sum(alive_at_omc :& (tsd_diff :< -0.001) :& (mTSD[.,omc] :< .))
		}
		st_numscalar("tsd_decreases", tsd_decreases)
	}
	if (tsd_decreases > 0) {
		n di as error "ERROR: mTSD decreased for `=tsd_decreases' patient-OMC combinations"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: mTSD is monotonically increasing"
	}

	// 2.2-2.5 Negative values checks
	mata {
		st_numscalar("negative_tne", sum(mTNE :< 0))
		st_numscalar("negative_tsd", sum(mTSD :< 0))
		st_numscalar("negative_tfi", sum(mTFI :< 0))
		st_numscalar("negative_txd", sum(mTXD :< 0))
	}
	
	if (negative_tne > 0) {
		n di as error "ERROR: `=negative_tne' negative values in mTNE"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All mTNE values are non-negative"
	}
	
	if (negative_tne > 0) {
		n di as error "ERROR: `=negative_tsd' negative values in mTSD"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All mTSD values are non-negative"
	}
	
	if (negative_tfi > 0) {
		n di as error "ERROR: `=negative_tfi' negative values in mTFI"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All mTFI values are non-negative"
	}

	if (negative_txd > 0) {
		n di as error "ERROR: `=negative_txd' negative values in mTXD"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All mTXD values are non-negative"
	}

	// 2.6 mTNE missingness by OMC (excluding death state OMC 19)
	mata {
		max_treatment_omc = min((cols(mTNE), 18))  // Exclude OMC 19 (death state)
		total_missing = 0
		
		for (omc = 1; omc <= max_treatment_omc; omc++) {
			// Patients who should have mTNE at this OMC:
			// - Alive (not dead at prior OMC)
			// - Started treatment at or before this OMC
			// - Not yet in death state (OMC 19)
			
			if (omc == 1) {
				should_have = (mState[.,1] :== 1)  // Only those starting at OMC 1
			}
			else {
				// Alive at prior OMC AND reached this OMC
				should_have = (mMOR[., omc-1] :== 0) :& (mState[.,1] :<= omc)
			}
			
			expected = sum(should_have)
			actual = sum(should_have :& (mTNE[., omc] :< .))
			missing = expected - actual
			
			if (missing > 0 & expected > 0) {
				pct = 100 * missing / expected
				printf("{txt}  OMC %2.0f: %6.0f/%6.0f eligible patients missing mTNE (%5.2f%%)\n", omc, missing, expected, pct)
				if (pct > 1) {
					st_numscalar("validation_warnings", st_numscalar("validation_warnings") + 1)
				}
				total_missing = total_missing + missing
			}
		}
		st_numscalar("total_tne_missing", total_missing)
	}
	if (total_tne_missing > 0) {
		n di as text "Total missing mTNE (OMC 1-18): `=total_tne_missing'"
	}
	else {
		n di as text "PASS: No missing mTNE"
	}

	// 2.7 Check odd OMCs: mTFI → mTNE
	mata {
		tfi_tne_errors = 0
		
		for (omc = 1; omc <= cols(mTNE); omc = omc + 2) {
			tfi_col = (omc + 1) / 2
			
			if (tfi_col <= cols(mTFI)) {
				tfi_exists_tne_missing = sum((mTFI[., tfi_col] :< .) :& (mTNE[., omc] :>= .))
				
				if (tfi_exists_tne_missing > 0) {
					printf("{error}  OMC %2.0f: %6.0f patients have mTFI[.,%g] but MISSING mTNE[.,%g]\n", omc, tfi_exists_tne_missing, tfi_col, omc)
					tfi_tne_errors = tfi_tne_errors + tfi_exists_tne_missing
				}
				else {
					printf("{txt}  OMC %2.0f: mTFI[.,%g] -> mTNE[.,%g] OK\n", omc, tfi_col, omc)
				}
			}
		}
		st_numscalar("tfi_tne_errors", tfi_tne_errors)
	}

	if (tfi_tne_errors > 0) {
		n di as error "ERROR: `=tfi_tne_errors' patients have mTFI calculated but not copied to mTNE"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: No patients have mTFI without mTNE"
	}

	// 2.8 Check even OMCs: mTXD → mTNE
	mata {
		txd_tne_errors = 0
		
		for (omc = 2; omc <= cols(mTNE); omc = omc + 2) {
			txd_col = omc / 2
			
			if (txd_col <= cols(mTXD)) {
				txd_exists_tne_missing = sum((mTXD[., txd_col] :< .) :& (mTNE[., omc] :>= .))
				
				if (txd_exists_tne_missing > 0) {
					printf("{error}  OMC %2.0f: %6.0f patients have mTXD[.,%g] but MISSING mTNE[.,%g]\n", omc, txd_exists_tne_missing, txd_col, omc)
					txd_tne_errors = txd_tne_errors + txd_exists_tne_missing
				}
				else {
					printf("{txt}  OMC %2.0f: mTXD[.,%g] -> mTNE[.,%g] OK\n", omc, txd_col, omc)
				}
			}
		}
		st_numscalar("txd_tne_errors", txd_tne_errors)
	}

	if (txd_tne_errors > 0) {
		n di as error "ERROR: `=txd_tne_errors' patients have mTXD calculated but not copied to mTNE"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: No patients have mTXD without mTNE"
	}
	
	// 3. Age
	n di as text _n "--- 3. Age ---"

	// 3.1 Age monotonicity
	mata {
		age_decreases = 0
		for (omc = 2; omc <= cols(mAge); omc++) {
			if (omc > 1) {
				alive_at_omc = (mMOR[., omc-1] :== 0)
			}
			else {
				alive_at_omc = J(rows(mAge), 1, 1)
			}
			age_diff = mAge[., omc] - mAge[., omc-1]
			age_decreases = age_decreases + sum(alive_at_omc :& (age_diff :< -0.001) :& (mAge[.,omc] :< .) :& (mAge[.,omc-1] :< .))
		}
		st_numscalar("age_decreases", age_decreases)
	}
	if (age_decreases > 0) {
		n di as error "ERROR: Age decreased for `=age_decreases' patient-OMC combinations"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: Age is monotonically increasing"
	}

	// 3.2 Age bounds
	mata {
		st_numscalar("age_out_of_range", sum(((mAge :< 18) :| (mAge :> 100)) :& (mAge :< .)))
	}
	if (age_out_of_range > 0) {
		n di as error "WARNING: `=age_out_of_range' age values outside 18-100 range"
		scalar validation_warnings = validation_warnings + 1
	}
	else {
		n di as text "PASS: All ages within reasonable bounds (18-100)"
	}

	// 3.3 Age missingness (abbreviated - just summary)
	mata {
		max_treatment_omc = min((cols(mAge), 18))
		total_age_missing = 0
		
		for (omc = 1; omc <= max_treatment_omc; omc++) {
			if (omc == 1) {
				should_have = (mState[.,1] :== 1)
			}
			else {
				should_have = (mMOR[., omc-1] :== 0) :& (mState[.,1] :<= omc)
			}
			expected = sum(should_have)
			actual = sum(should_have :& (mAge[., omc] :< .))
			total_age_missing = total_age_missing + (expected - actual)
		}
		st_numscalar("total_age_missing", total_age_missing)
	}
	if (total_age_missing > 0) {
		n di as text "INFO: `=total_age_missing' total missing age values (OMC 1-18)"
	}
	else {
		n di as text "PASS: No missing age values"
	}

	// 4. BCR
	n di as text _n "--- 4. BCR ----"
	
	// 4.1: L1 BCR (valid range: 1-6)
	mata {
		vBCR_L1 = mBCR[., 1]
		vL1_invalid = (vBCR_L1 :< 1) :| (vBCR_L1 :> 6)
		vL1_invalid = vL1_invalid :& (vBCR_L1 :< .)
		st_numscalar("bcr_l1_invalid", sum(vL1_invalid))
	}

	if (bcr_l1_invalid > 0) {
		n di as error "  FAIL: `=bcr_l1_invalid' L1 BCR values outside valid range 1-6"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All L1 BCR values in valid range 1-6"
	}

	// 4.2: ASCT BCR (valid range: 1-4)
	mata {
		vBCR_ASCT = mBCR[., 10]
		vASCT_invalid = (vBCR_ASCT :< 1) :| (vBCR_ASCT :> 4)
		vASCT_invalid = vASCT_invalid :& (vBCR_ASCT :< .)
		st_numscalar("bcr_asct_invalid", sum(vASCT_invalid))
	}

	if (bcr_asct_invalid > 0) {
		n di as error "  FAIL: `=bcr_asct_invalid' ASCT BCR values outside valid range 1-4"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All ASCT BCR values in valid range 1-4"
	}

	// 4.3: L2+ BCR (valid range: 1-6)
	mata {
		nL2plus_invalid = 0
		
		for (col = 2; col <= 9; col++) {
			vBCR_Ln = mBCR[., col]
			
			vInvalid = (vBCR_Ln :< 1) :| (vBCR_Ln :> 6)
			vInvalid = vInvalid :& (vBCR_Ln :< .)
					
			if (sum(vInvalid) > 0) {
				line_num = col + 1
				printf("  FAIL: %f invalid values in L%f BCR outside valid range 1-6\n", 
					   sum(vInvalid), line_num)
				nL2plus_invalid = nL2plus_invalid + sum(vInvalid)
			}
		}
		
		st_numscalar("bcr_l2plus_invalid", nL2plus_invalid)
	}

	if (bcr_l2plus_invalid > 0) {
		n di as error "  FAIL: Total `=bcr_l2plus_invalid' L2+ BCR values outside valid range 1-6"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All L2+ BCR values in valid range 1-6"
	}

	// 5 BCR Distribution
	n di as text _n "--- 5. BCR Distribution ---"

	// Line 1
	capture mata: st_numscalar("n_bcr", sum(mBCR[., 1] :< .))
	if _rc == 0 & n_bcr > 0 {
		mata {
			bcr_col = mBCR[., 1]
			n = sum(bcr_col :< .)
			st_numscalar("n_bcr", n)
			st_numscalar("pct1", 100 * sum(bcr_col :== 1) / n)
			st_numscalar("pct2", 100 * sum(bcr_col :== 2) / n)
			st_numscalar("pct3", 100 * sum(bcr_col :== 3) / n)
			st_numscalar("pct4", 100 * sum(bcr_col :== 4) / n)
			st_numscalar("pct5", 100 * sum(bcr_col :== 5) / n)
			st_numscalar("pct6", 100 * sum(bcr_col :== 6) / n)
		}
		n di as text "Line 1 (n=" %7.0fc n_bcr "): " ///
			"CR=" %4.1f pct1 "% " ///
			"VG=" %4.1f pct2 "% " ///
			"PR=" %4.1f pct3 "% " ///
			"MR=" %4.1f pct4 "% " ///
			"SD=" %4.1f pct5 "% " ///
			"PD=" %4.1f pct6 "%"
	}


	// ASCT (column 10, categories 1-4)
	capture mata: st_numscalar("n_bcr", sum(mBCR[., 10] :< .))
	if _rc == 0 & n_bcr > 0 {
		mata {
			bcr_col = mBCR[., 10]
			n = sum(bcr_col :< .)
			st_numscalar("n_bcr", n)
			st_numscalar("pct1", 100 * sum(bcr_col :== 1) / n)
			st_numscalar("pct2", 100 * sum(bcr_col :== 2) / n)
			st_numscalar("pct3", 100 * sum(bcr_col :== 3) / n)
			st_numscalar("pct4", 100 * sum(bcr_col :== 4) / n)
		}
		n di as text "ASCT   (n=" %7.0fc n_bcr "): " ///
			"CR=" %4.1f pct1 "% " ///
			"VG=" %4.1f pct2 "% " ///
			"PR=" %4.1f pct3 "% " ///
			"MR=" %4.1f pct4 "%" 
	}
	
	// Lines 2-9
	forvalues line = 2/9 {
		capture mata: st_numscalar("n_bcr", sum(mBCR[., `line'] :< .))
		if _rc == 0 & n_bcr > 0 {
			mata {
				bcr_col = mBCR[., `line']
				n = sum(bcr_col :< .)
				st_numscalar("n_bcr", n)
				st_numscalar("pct1", 100 * sum(bcr_col :== 1) / n)
				st_numscalar("pct2", 100 * sum(bcr_col :== 2) / n)
				st_numscalar("pct3", 100 * sum(bcr_col :== 3) / n)
				st_numscalar("pct4", 100 * sum(bcr_col :== 4) / n)
				st_numscalar("pct5", 100 * sum(bcr_col :== 5) / n)
				st_numscalar("pct6", 100 * sum(bcr_col :== 6) / n)
			}
			n di as text "Line `line' (n=" %7.0fc n_bcr "): " ///
				"CR=" %4.1f pct1 "% " ///
				"VG=" %4.1f pct2 "% " ///
				"PR=" %4.1f pct3 "% " ///
				"MR=" %4.1f pct4 "% " ///
				"SD=" %4.1f pct5 "% " ///
				"PD=" %4.1f pct6 "%"
		}
	}

	// 6. Overall Survival
	n di as text _n "--- 6. Overall Survival ---"

	mata {
		st_numscalar("negative_os", sum(mOS :< 0))
	}
	if (negative_os > 0) {
		n di as error "ERROR: `=negative_os' negative OS values"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All OS values are non-negative"
	}

	// OS Statistics - use Stata's summarize
	mata {
		// Get final OS for each patient
		final_os = J(rows(mOS), 1, .)
		for (i = 1; i <= rows(mOS); i++) {
			for (j = cols(mOS); j >= 1; j--) {
				if (mOS[i, j] < .) {
					final_os[i] = mOS[i, j]
					break
				}
			}
		}
		// Store to Stata
		st_store(., st_addvar("double", "temp_os"), final_os)
	}
	n di as text "OS Statistics (months):"
	quietly summarize temp_os, detail
	n di as text "  N: " %12.0fc r(N)
	n di as text "  Mean: " %12.1f r(mean)
	n di as text "  Median: " %12.1f r(p50)
	n di as text "  P10: " %12.1f r(p10)
	n di as text "  P25: " %12.1f r(p25)
	n di as text "  P75: " %12.1f r(p75)
	n di as text "  P90: " %12.1f r(p90)
	drop temp_os

	// 7. ASCT
	n di as text _n "--- 7. ASCT ---"

	mata {
		st_numscalar("asct_invalid", sum((vSCT_L1 :!= 0) :& (vSCT_L1 :!= 1) :& (vSCT_L1 :< .)))
		st_numscalar("asct_rate", 100 * sum(vSCT_L1 :== 1) / sum(vSCT_L1 :< .))
	}
	if (asct_invalid > 0) {
		n di as error "ERROR: `=asct_invalid' non-binary ASCT values"
		scalar validation_errors = validation_errors + 1
	}
	else {
		n di as text "PASS: All ASCT values are binary"
	}
	n di as text "ASCT rate: " %5.1f asct_rate "%"

	// 8. Cross-Matrix
	n di as text _n "--- 8. Cross-Matrix ---"

	mata {
		outcomes_after_death = 0
		for (i = 1; i <= rows(mMOR); i++) {
			death_omc = .
			for (j = 1; j <= cols(mMOR); j++) {
				if (mMOR[i, j] == 1) {
					death_omc = j
					break
				}
			}
			if (death_omc < .) {
				for (j = death_omc + 1; j <= cols(mTNE); j++) {
					if (mTNE[i, j] < .) {
						outcomes_after_death++
						break
					}
				}
			}
		}
		st_numscalar("outcomes_after_death", outcomes_after_death)
	}
	if (outcomes_after_death > 0) {
		n di as error "WARNING: `=outcomes_after_death' patients have mTNE values after death"
		scalar validation_warnings = validation_warnings + 1
	}
	else {
		n di as text "PASS: No outcomes recorded after death"
	}

	// 9. Patient Flow
	n di as text _n "--- 9. Patient Flow ---"
	
	n di as text "  " %5s "OMC" %12s "Stage" %11s "Alive" %11s "Dead" %12s "Cumulative" %11s "Started"
	n di as text "  " %5s "-----" %12s "------------" %11s "-----------" %11s "-----------" %12s "------------" %11s "-----------"

	mata {
		stage_labels = ("DN", "L1S", "L1E", "L2S", "L2E", "L3S", "L3E", "L4S", "L4E", "L5S", "L5E", "L6S", "L6E", "L7S", "L7E", "L8S", "L8E", "L9S", "L9E")'
		cum_dead = 0
		for (omc = 1; omc <= min((cols(mMOR), 19)); omc++) {
			st_numscalar("temp_omc", omc)
			st_numscalar("temp_alive", sum(mMOR[., omc] :== 0))
			st_numscalar("temp_dead", sum(mMOR[., omc] :== 1))
			cum_dead = sum(mMOR[., omc] :== 1) + sum(mMOR[., omc] :== .)
			st_numscalar("temp_cumulative", cum_dead)
			st_numscalar("temp_started", sum(mState[.,1] :== omc))
			
			stata(`"n di as text "  " %5.0f scalar(temp_omc) %12s ""' + stage_labels[omc] + `"" %11.0fc scalar(temp_alive) %11.0fc scalar(temp_dead) %12.0fc scalar(temp_cumulative) %11.0fc scalar(temp_started)"')
		}
	}

/*
	// 10. Patient-level Diagnositcs (uncomment if errors found)
	if (validation_errors > 0) {
		n di as text _n "--- 10. PATIENT-LEVEL DIAGNOSTICS ---"
		n di as text "Showing detailed traces for patients with issues..."
		
		mata {
			real scalar i, omc, last_alive_omc, sample_count, n_problems
			real colvector missing_os_alive
			
			missing_os_alive = J(rows(mOS), 1, 0)
			
			for (i = 1; i <= rows(mOS); i++) {
				last_alive_omc = 0
				for (omc = 1; omc <= cols(mMOR); omc++) {
					if (mMOR[i, omc] == 0) {
						last_alive_omc = omc
					}
				}
				
				if (last_alive_omc > 0 & last_alive_omc <= cols(mOS)) {
					if (mOS[i, last_alive_omc] >= .) {
						missing_os_alive[i] = 1
					}
				}
			}
			
			n_problems = sum(missing_os_alive)
			st_numscalar("n_problems", n_problems)
			
			if (n_problems > 0) {
				printf("\nFound %f patients with missing OS while alive\n", n_problems)
				printf("Showing first 3 patient traces:\n\n")
				
				sample_count = 0
				
				for (i = 1; i <= rows(mOS) & sample_count < 3; i++) {
					if (missing_os_alive[i] == 1) {
						sample_count++
						printf("Patient %f (State=%f, ASCT=%f, BCR_L1=%f):\n", 
							   i, mState[i,1], vSCT_L1[i], mBCR[i,1])
						
						printf("  OMC:  ")
						for (omc = 1; omc <= 10 & omc <= cols(mTNE); omc++) printf("%8.0f ", omc)
						printf("\n  mTNE: ")
						for (omc = 1; omc <= 10 & omc <= cols(mTNE); omc++) {
							if (mTNE[i, omc] < .) printf("%8.1f ", mTNE[i, omc])
							else printf("%8s ", ".")
						}
						printf("\n  mOS:  ")
						for (omc = 1; omc <= 10 & omc <= cols(mOS); omc++) {
							if (mOS[i, omc] < .) printf("%8.1f ", mOS[i, omc])
							else printf("%8s ", ".")
						}
						printf("\n  mMOR: ")
						for (omc = 1; omc <= 10 & omc <= cols(mMOR); omc++) {
							printf("%8.0f ", mMOR[i, omc])
						}
						printf("\n\n")
					}
				}
			}
		}
	}
*/

	// Validation Summary
	n di as text _n "--- Validation Outcome ---"

	if (validation_errors > 0) {
		n di as error "Validation Failed: `=validation_errors' error(s) found"
	}
	else if (validation_warnings > 0) {
		n di as text "Validation Passed with `=validation_warnings' warning(s)"
	}
	else {
		n di as result "Validation Passed: All checks passed"
	}
}
