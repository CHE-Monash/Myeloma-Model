**********
* Monash Myeloma Model - base_model: run.do (analysis runbook)
*
* Current-practice projection over the synthetic population (all regimens), in order:
* risk equations -> simulate. Run from the repository root. MRDR data via $data_path (config.do).
* simulate.do is the simulation dispatcher.
*
* Depends on the SHARED prep outputs (built once by prep/, not per-analysis):
*   - ${data_path}/MRDR Long MI.dta            (prep/multiple_imputation.do 10 0)  -- risk-equation input
*   - patients/population_1995_2040_*.dta     (prep/population_1995_2040.do)       -- the simulation cohort
*
* Two tracks:
*   DETERMINISTIC (point estimate)  -- light; runs locally, top to bottom (the live `do` lines).
*   BOOTSTRAP (prediction intervals) -- HEAVY. The risk-equation and 500-replicate simulation steps run
*       on the HPC. This file is NOT run on the HPC -- the bootstrap section below is the canonical
*       record of the sequence (left commented), with the equivalent sbatch commands.
**********

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"


**********
* DETERMINISTIC (point estimate) -- runs locally, top to bottom
**********

* 1. Risk equations -> analyses/base_model/coefficients/coefficients_base_model.mmat
*    args: analysis coeffs min_year max_year boot [min_bs max_bs]   (loads outcomes/txr_base_model.do)
do "prep/risk_equations.do" base_model base_model 1995 2040 0

* 2. Simulate the population (point estimate, $boot 0) -> analyses/base_model/simulated/...
*    simulate.do also runs core/validation.do (and core/generate_report.do when $report 1).
do "analyses/base_model/simulate.do"


**********
* BOOTSTRAP (prediction intervals) -- HEAVY; runs on the HPC.
* Left commented: record of the sequence, not executed here. (When this analysis is first run on the
* cluster, add a /* ... */ HPC shell block like analyses/oos/run.do, submitting the generic
* hpc/risk_equations.script with `sbatch --export=ALL,ANALYSIS=base_model,COEFFS=base_model,MINYR=1995,MAXYR=2040`.)
**********

*  (a) Bootstrap risk equations -> analyses/base_model/coefficients/bootstrap/coefficients_base_model_B1..500
*      Requires the SHARED bootstrap MI ${data_path}/bootstrap/MRDR Long MI B<b>.dta
*      (prep/multiple_imputation.do 10 1 1 500 -- the main-model bootstrap MI; itself an HPC step).
*      LOCAL (slow, serial):
* do "prep/risk_equations.do" base_model base_model 1995 2040 1 1 500
*      HPC (array 1-500; from ~/em76/adam):
*        sbatch <base_model risk_equations.script>   // -> stata-mp prep/risk_equations.do base_model base_model 1995 2040 1 <task> <task>
*
*  (b) 500 bootstrap simulations of the population -> analyses/base_model/simulated/bootstrap/
*      LOCAL:  do "analyses/base_model/simulate.do" 1 1 500
*      HPC (array; chunk the range, e.g. 1 101 200):
*        sbatch <base_model simulate.script>          // -> stata-mp analyses/base_model/simulate.do 1 <task> <task>
