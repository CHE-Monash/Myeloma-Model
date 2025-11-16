**********
*EpiMAP Myeloma - Load Patients
**********

capture program drop load_patients
program define load_patients

	// Determine data source based on $data_type (not $data)
	if ("$data_type" == "population") {
		// Use population number extracted by main dispatcher
		use "patients/population_1995_2040_${pop_number}.dta", clear
		di as text "Loaded population ${pop_number}"
	}
	else if ("$data_type" == "predicted") {
		use "$patients_path/patients_${analysis}.dta", clear
		di as text "Loaded ${analysis} specific patient data"
	}
	else {
		di as error "Unknown data type: $data_type"
		di as error "Valid options: Population_#, Predicted"
		exit 198
	}
	
	// Filters
	qui keep if YearDN >= $min_year & YearDN <= $max_year 	// Based on year of diagnosis
	qui keep if State <= ($line * 2) + 2 					// Based on disease stage
	qui replace ID = _n
	qui keep if ID >= $min_id & ID <= $max_id				// Based on ID
	
	qui sum ID
	scalar Obs = r(N)
	di as text "Final sample size: `=Obs' patients"
	
	mata: Limit = 100 // Age limit
	
end
