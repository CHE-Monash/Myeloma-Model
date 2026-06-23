**********
* Monash Myeloma Model - Run Pipeline (shared lean engine pass)
*
* Purpose: One definition of the core simulation pass used by every orchestrator
*          (build_cohort_pool.do, transport_dvd.do, ce_convergence.do, ...).
*          Runs entry/diagnosis -> processed per-patient outcomes, WITHOUT the
*          CSV export (callers that want CSVs run core/export_results.do
*          separately). Replaces the per-script "simulation_pipeline" copies.
*
* Requires (run once by the caller, before calling run_pipeline):
*   core/load_patients.do, core/mata_setup.do, core/simulation_engine.do,
*   core/process_data.do
* Uses the usual globals the caller sets: $data, $int, $line, $scenario,
*   $coeffs, $boot, paths, etc. Coefficients must be loaded (mata matuse) before
*   the call.
**********

capture program drop run_pipeline
program define run_pipeline
	run "core/mata_functions.do"
	run "core/rng_slots.do"
	load_patients
	mata_setup
	simulation
	process_data
end
