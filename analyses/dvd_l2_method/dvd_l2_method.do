**********
*EpiMAP Myeloma 
*
*Simulation for DVd L2 Method analysis
**********

**********
*Analysis Configuration
	local analysis_name "dvd_l2_method"
	local coefficient_file "coefficients_$coeffs"

**********
*Load Programs
	
	quietly do "core/load_patients.do"
	quietly do "core/mata_setup.do"
	quietly do "core/simulation_engine.do"
	quietly do "core/process_data.do"
	
**********
*Determine processing approach
	if("$boot" == "0") {
		// No Bootstrapping
		qui mata: mata matuse "$coefficients_path/`coefficient_file'"
				
		quietly do "core/mata_functions.do"
		load_patients
		mata_setup
		simulation
		process_data
		
		save "$simulated_path/${int}_${line}_${data}_${min_id}_${max_id}_${scenario}.dta", replace
		di as text "Results saved to: $simulated_path/${int}_${line}_${data}_${min_id}_${max_id}.dta"
	}
	else {
		// Bootstrapping
		di as text "Running bootstrap analysis with ${max_BS} iterations..."
		
		forvalues b = $min_bs / $max_bs {
			di as text _newline "========================================" 
			di as text "Processing bootstrap iteration `b' of ${max_BS}..."
			di as text "========================================"
			
			global BSIteration "`b'"
			
			mata: mata clear
			mata: mata matuse "$coefficients_path/bootstrap/`coefficient_file'_B`b'"
			
			quietly do "core/mata_functions.do"
			load_patients
			mata_setup
			simulation
			process_data
			
			save "$simulated_path/bootstrap/${int}_${line}_${data}_${min_id}_${max_id}_${scenario}_B`b'.dta", replace
			di as text "Bootstrap iteration `b' completed"
		}
		
		di as result _newline "All bootstrap iterations completed successfully!"
		di as text "Bootstrap results saved in: $simulated_path/bootstrap/"
		di as text "Files saved as: ${int}_${line}_${data}_${min_id}_${max_id}_B[$MinBS-$MaxBS].dta"
	}
