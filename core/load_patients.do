**********
* Monash Myeloma Model - Load Patients
*
* Purpose: Load the cohort (synthetic, $cohort_file override, or production file) and
*          filter by year of diagnosis, disease stage and ID range; reset ID = _n.
**********

capture program drop load_patients
program define load_patients

	// Synthetic incidence cohort (default projection): "synthetic" or "synthetic_N" -> synthetic_1995_2040_N.
	if regexm("$data", "synthetic") {
		global data_type = "synthetic"
		if regexm("$data", "synthetic_([0-9]+)") global pop_number = regexs(1)
		else global pop_number = 1
	}
	else {
		global data_type = "$data" // Predicted
		global pop_number = 1  // Default cohort
	}

	// Determine data source
	if ("$data_type" == "synthetic") {
		use "patients/synthetic_1995_2040_${pop_number}.dta", clear
	}
	else if ("$cohort_file" != "") {
		// $cohort_file lets a caller (e.g. the default outsample scenario -> patients_test.dta, or
		// ce_precision) read a different cohort without overwriting the production file -- honoured for
		// any non-synthetic $data.
		use "$cohort_file", clear
	}
	else {
		use "$patients_path/patients_${analysis}_${line}.dta", clear
	}
	
	// Filters
	qui keep if YearDN >= $min_year & YearDN <= $max_year 	// Based on year of diagnosis
	qui keep if State <= ($line * 2) + 2 					// Based on disease stage
	qui replace ID = _n										// Reset ID
	qui keep if ID >= $min_id & ID <= $max_id				// Based on ID
	
	qui sum ID
	scalar Obs = r(N)
	di as text "Final sample size: `=Obs' patients"
	
	mata: Limit = 100 // Age limit
	
end
