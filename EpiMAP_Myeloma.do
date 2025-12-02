**********
*EpiMAP Myeloma v2.0 - Main Dispatcher
**********
	
	clear
	clear mata
	clear matrix
	set more off
	
**********
* Validate required globals are set
foreach req in analysis int line coeffs data min_year max_year min_id max_id boot {
    capture confirm existence $`req'
    if _rc {
        di as error "Error: Global `req' not set"
        di as error "This script must be called from run.do"
        exit 198
    }
}

**********
*Parse Data parameter for population-specific requests
	if regexm("$data", "population_([0-9]+)") {
		global pop_number = regexs(1)
		global data_type = "population"
		di "Using specific population: ${pop_number}"
	}
	else {
		global data_type = "$data"
		global pop_number = 1  // Default population
	}

**********
*Validate Analysis exists
	capture confirm file "analyses/$analysis/${analysis}.do"
	if _rc != 0 {
		di as error "Error: Analysis '$analysis' not found."
		di as error "Available analyses:"
		
		local analyses : dir "analyses" dirs "*"
		foreach analysis of local analyses {
			di "  - `analysis'"
		}
		exit 601
	}

**********
*Validate Population data if needed
	if ("$data_type" == "population") {
		capture confirm file "patients/population_1995_2040_${pop_number}.dta"
		if _rc != 0 {
			di as error "Error: Population ${pop_number} not found."
			di as error "Available populations:"
			local populations : dir "data/populations" files "population_*.dta"
			foreach pop of local populations {
				di "  - `pop'"
			}
			exit 601
		}
	}

**********
*Set paths based on analysis
	global analysis_path "analyses/$analysis"
	global coefficients_path "$analysis_path/data/coefficients"
	global patients_path "$analysis_path/data/patients" 
	global simulated_path "$analysis_path/data/simulated"
	global populations_path "data/populations"
	
	// Create results directory if it doesn't exist
	capture mkdir "$results_path"

**********
*Run analysis
	di as text "Running EpiMAP Myeloma Analysis: $analysis"
	di as text "Intervention: $int, Line: $line, Data: $data_type"
	if ("$data_type" == "population") di as text "Population: ${pop_number}"
	di as text "Patient range: $min_ID to $max_ID"
	
	// Run the specific analysis
	do "$analysis_path/${analysis}.do"
	
	di "Analysis completed successfully."
	
**********
*Generate report
	if ("$report" == "1" & "$boot" == "0") {
		qui do "core/generate_report.do"
	}

**********
*Final summary
	di as result "{bf:EpiMAP Myeloma v2.0 - Execution Complete}"
	di as text "Simulated data saved to:"
	if ("$boot" == "0") {
		di as result "  $simulated_path/${int}_${line}_${data}_${max_id}_${max_id}_${scenario}.dta"
	}
	else {
		di as result "  $simulated_path/bootstrap/ (samples ${min_BS} to ${max_BS})"
	}
	if ("$report" == "1") {
		di as text "Report saved to:"
		di as result "  $simulated_path/reports/"
	}
	