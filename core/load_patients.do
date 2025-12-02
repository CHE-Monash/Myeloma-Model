**********
*EpiMAP Myeloma - Load Patients
**********

capture program drop load_patients
program define load_patients

	// Determine data source based on $data_type
	if ("$data_type" == "population") {
		use "patients/population_1995_2040_${pop_number}.dta", clear
	}
	else if ("$data_type" == "predicted") {
		use "$patients_path/patients_${analysis}.dta", clear
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
