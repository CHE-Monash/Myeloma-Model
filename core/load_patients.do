**********
*EpiMAP Myeloma - Load Patients
**********

capture program drop load_patients
program define load_patients

	// Determine data source based on $data_type (not $data)
	if ("$data_type" == "population") {
		// Use population number extracted by main dispatcher
		use "patients/population/1995-2040/population${pop_number}.dta", clear
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

	// Standard data preparation and filtering
	keep if State <= ($line * 2) + 2
	replace ID = _n
	keep if ID >= $min_id & ID <= $max_id
	
	quietly sum ID
	scalar Obs = r(N)
	di as text "Final sample size: `=Obs' patients (IDs $min_ID to $max_ID)"
	
	mata: Limit = 100 // Age limit
	
end
