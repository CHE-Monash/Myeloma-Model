**********
* Monash Myeloma Model - Run (template analysis runbook)
*
* Purpose: the analysis runbook -- the full ordered pipeline (fit -> simulate -> validate) as runnable
*          `do` lines, plus the bootstrap HPC plumbing. COPY THIS FOLDER to analyses/<your_analysis>/,
*          adapt, and fill in every <...>.
* Usage:   run from the repository root; MRDR data via $data_path (config.do); simulate.do is the
*          dispatcher. Needs shared prep outputs (built once by prep/): ${data_path}/MRDR Long MI.dta
*          (prep/multiple_imputation.do 10 0) and analyses/template/patients/<your_cohort>.dta.
* Notes:   Worked runbooks to copy from -- analyses/base_model/run.do (simplest); analyses/oos/run.do
*          (full HPC shell block: rsync + sbatch arrays); analyses/transport_dvd/run.do (scenarios + an
*          extra coefficient-generation step).
**********

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"


**********
* DETERMINISTIC (point estimate) -- runs locally, top to bottom
**********

* 1. Risk equations -> analyses/template/coefficients/coefficients_template.mmat
*    args: analysis coeffs min_year max_year boot [min_bs max_bs]   (loads outcomes/txr_template.do)
do "prep/risk_equations.do" template template 1995 2040 0

* 1b. (LINE-SPECIFIC analyses only) Build the line-L simulation cohort -> patients/patients_template_<line>.dta.
*     Skip for whole-population analyses (point $data at a population cohort instead). Needs step 1's coefficients.
* do "analyses/template/patients/cohort_pool.do"     // build the line-entry pool once (expensive)
* do "analyses/template/patients/draw_cohort.do"     // draw a fixed-size, seeded cohort from the pool

* 2. Simulate (point estimate, $boot 0) -> analyses/template/simulated/...
*    simulate.do also runs core/validation.do (and core/generate_report.do when $report 1).
*    If your analysis has scenarios, loop them via the 4th positional arg (see transport_dvd/run.do):
*    foreach scen in <A> <B> { do "analyses/template/simulate.do" 0 1 1 "`scen'" }
do "analyses/template/simulate.do"


**********
* BOOTSTRAP (prediction intervals) -- HEAVY; runs on the HPC.
* The `do ...` lines below are the LOCAL (serial, slow) equivalent. For the full HPC plumbing -- rsync
* code+data up, sbatch the generic array jobs (hpc/{multiple_imputation,risk_equations,simulate}.script),
* pull results back -- copy the /* ... */ shell block from analyses/oos/run.do and swap in `template`.
**********

*  (a) SHARED bootstrap MI (FULL data) -> ${data_path}/bootstrap/MRDR Long MI B1..500.dta
*      Reuse the main-model bootstrap MI if base_model already produced it.
*      LOCAL:  do "prep/multiple_imputation.do" 10 1 1 500
*
*  (b) Bootstrap risk equations -> analyses/template/coefficients/bootstrap/coefficients_template_B1..500
*      LOCAL:  do "prep/risk_equations.do" template template 1995 2040 1 1 500
*      HPC:    sbatch --export=ALL,ANALYSIS=template,COEFFS=template,MINYR=1995,MAXYR=2040 hpc/risk_equations.script
*
*  (c) 500 bootstrap simulations -> analyses/template/simulated/bootstrap/
*      LOCAL:  do "analyses/template/simulate.do" 1 1 500
*      HPC:    sbatch --export=ALL,ANALYSIS=template hpc/simulate.script
*
*  (d) Aggregate the 500 sims into prediction intervals -- add an aggregation script for your analysis
*      (cf. analyses/oos/bootstrap_validation.do, analyses/transport_dvd/simulated/bootstrap_summary.do).
