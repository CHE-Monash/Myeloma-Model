**********
* Monash Myeloma Model - Run Pipeline (shared lean engine pass)
*
* Purpose: The single definition of the core simulation pass (load_patients -> mata_setup
*          -> simulation -> process_data) used by every orchestrator, WITHOUT CSV export
*          (callers run core/export_results.do separately).
* Notes:   Caller must first run core/load_patients.do, core/mata_setup.do,
*          core/simulation_engine.do, core/process_data.do, load coefficients (mata matuse)
*          and set the usual globals ($data, $int, $line, $scenario, $coeffs, $boot, paths).
**********

capture program drop run_pipeline
program define run_pipeline
	// Compile the persistent Mata utility functions once per (cleared) Mata
	// state. The definition files error ("... already exists") if re-run while
	// their functions are still compiled, so load them only when absent: true on
	// a fresh session and again after the bootstrap loop's `mata clear`, but a
	// no-op if run_pipeline is called repeatedly within one live Mata state.
	capture mata: _rp_probe = &get_txr_coef()
	if (_rc) {
		run "core/mata_functions.do"
		run "core/rng_slots.do"
	}
	load_patients
	mata_setup
	simulation
	process_data
end
