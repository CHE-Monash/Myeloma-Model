**********
* Monash Myeloma Model - Load Patients
* 
* Purpose: Load and filter patient data file
**********

capture program drop load_patients
program define load_patients

	// Handle population-specific requests (e.g., population_1, population_2)
	if regexm("$data", "population_([0-9]+)") {
		global pop_number = regexs(1)
		global data_type = "population"
		di as text "Using specific population: ${pop_number}"
	}
	else {
		global data_type = "$data" // Predicted
		global pop_number = 1  // Default population
	}

	// Determine data source
	if ("$data_type" == "population") {
		use "patients/population_1995_2040_${pop_number}.dta", clear
	}
	else if ("$cohort_file" != "") {
		// $cohort_file lets a caller (e.g. ce_precision, analyses/oos) read a different cohort
		// without overwriting the production file -- honoured for any non-population $data.
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
