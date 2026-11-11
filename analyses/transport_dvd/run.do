**********
* Monash Myeloma Model - transport_dvd runbook
*
* Purpose: Held-out calibrated transport of the CASTOR DVd effect to a second-line MRDR cohort, in
*          order: risk equations -> transport coefficients -> simulate (x3 scenarios) -> aggregate.
*          Three scenarios, all simulated the same way (dispatcher's $scenario, 4th arg): A_trial
*          (resampled CASTOR BCR - trial-based prediction), B_transport (ologit BCR = MRDR + DVd +
*          line indicators - calibrated transport, the headline), C_mrdr (native MRDR BCR benchmark).
* Usage:   Run from the repository root. MRDR data via $data_path (config.do); simulate.do is the
*          simulation dispatcher.
* Notes:   Depends on shared prep output ${data_path}/MRDR Long MI.dta, built by prep/ (see
*          PREREQUISITES below). Three tracks: DETERMINISTIC (light, runs locally top to bottom via
*          the live `do` lines), COHORT (build-once, HEAVY, deliberately commented) and BOOTSTRAP
*          (the headline metric; HEAVY, runs on the HPC). This file is NOT run on the HPC; the
*          bootstrap section is the canonical record of that sequence, with sbatch and rsync plumbing
*          in a /* ... */ shell block at the end.
**********

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

**********
* PREREQUISITES -- shared prep, not part of this analysis. Rebuild only when the extraction or the
* imputation changes; both write to $data_path and feed step 1. Commented: they are slow and shared.
**********

do "prep/data_extraction.do"              // -> MRDR Long.dta
do "prep/multiple_imputation.do" 10 0     // -> MRDR Long MI.dta   (args: imp boot [min max sample])

**********
* DETERMINISTIC (point estimate) -- runs locally, top to bottom
**********

* 1. Risk equations on the FULL data -> analyses/transport_dvd/coefficients/coefficients_transport_dvd.mmat
*    args: analysis coeffs min_year max_year boot [min_bs max_bs]   (loads outcomes/txr_transport_dvd.do)
do "prep/risk_equations.do" transport_dvd transport_dvd 2020 2025 0

* 2. Transport coefficients (beta_MRDR / beta_DVd) + Table 1 baselines + observed DVd BCR benchmark
*    -> outcomes/{A_trial/bcr_*_l2.mmat, B_transport/transport_dvd.mmat}, results/*.csv
*    args: boot [min_bs max_bs]
do "analyses/transport_dvd/outcomes/calibrated_transport.do" 0

**********
* 3. COHORT (build once) -- HEAVY, and deliberately NOT part of a top-to-bottom run. Uncomment and
*    run by hand only when the coefficients, the engine or the synthetic cohorts have changed.
*
*    IT RUNS HERE, AFTER STEPS 1 AND 2, NOT BEFORE THEM: cohort_pool.do loads
*    coefficients_transport_dvd.mmat and simulates the full pathway with it, and ce_precision.do
*    reads step 2's BCR .mmat artefacts. Changing the cohort moves the ICER.
**********

*  (a) Line-entry pool. Simulates patients/synthetic_1995_2040_1..10.dta from diagnosis, keeps
*      everyone alive at end of L1 reaching L2 in the 2020-2025 case-mix window, blanks the L2
*      regimen and rebases the diagnosis clock.  -> patients/cohort_pool_2.dta
*      Do NOT widen the window to change the pool size: that shifts the case-mix. Size is (b)'s job.
do "analyses/transport_dvd/patients/sample_size/cohort_pool.do"
*
*  (b) Production CE cohort. Fixed-size, fixed-seed draw from the pool (N = 50,000, seed 20260616;
*      N is the Monte Carlo precision choice, see results/ce_sample_size_appendix.md).
*      -> patients/patients_transport_dvd_2.dta, which load_patients.do reads as $data "predicted"
do "analyses/transport_dvd/patients/sample_size/ce_cohort.do"
*
*  (c) (optional) Re-justify N. Only needed when the engine or cohort moves sigma_pp materially.
*      -> results/ce_precision_sigma.csv
do "analyses/transport_dvd/patients/sample_size/ce_precision.do"

**********
* 4. Simulate each scenario (point estimate, $boot 0). The dispatcher runs BOTH arms (DVd vs Vd) and
*    prints the ICER; the 4th positional arg selects the scenario (min_bs/max_bs unused at boot 0).
*    -> analyses/transport_dvd/simulated/<scenario>/{dvd,vd}_2_predicted.dta (+ econ/patient CSVs)
**********

foreach scen in A_trial B_transport C_mrdr {
    do "analyses/transport_dvd/simulate.do" 0 1 1 "`scen'"
}

* Figure 1 (requires both deterministic and bootstrap results) -> results/
do "analyses/transport_dvd/draw_figure1.do"

**********
* BOOTSTRAP (prediction intervals) -- the headline transport metric; HEAVY; runs on the HPC.
* The `do ...` lines below are the LOCAL (serial, slow) equivalent. The actual HPC plumbing -- rsync
* code+data up, sbatch the array jobs, pull results back -- is the /* ... */ shell block at the END of
* this file: it is a Stata block comment (inert here), so SELECT-AND-RUN those lines in a VS Code
* TERMINAL (bash). The slurm scripts (hpc/multiple_imputation.script, hpc/risk_equations.script,
* hpc/simulate.script) are generic array jobs 1-500 -- the per-analysis args ride on `sbatch --export=`.
**********

*  (a) SHARED bootstrap MI (FULL data, no fold) -> ${data_path}/bootstrap/MRDR Long MI B1..500.dta
*      This is the main-model bootstrap MI -- if base_model already produced it, REUSE it (skip this).
*      LOCAL equivalent (slow, serial):
* do "prep/multiple_imputation.do" 10 1 1 500
*
*  (b) Bootstrap risk equations -> analyses/transport_dvd/coefficients/bootstrap/coefficients_transport_dvd_B1..500
*      LOCAL:  do "prep/risk_equations.do" transport_dvd transport_dvd 1995 2025 1 1 500
*
*  (c) Bootstrap transport coefficients (shared trial draw per replicate; reads MRDR Long MI B<b>)
*      -> outcomes/{A_trial,B_transport}/bootstrap/*_B1..500.mmat
*      LOCAL:  do "analyses/transport_dvd/outcomes/calibrated_transport.do" 1 1 500
*
*  (d) 3x500 bootstrap simulations (one array per scenario) -> simulated/<scenario>/bootstrap/{dvd,vd}_2_predicted_B<b>.dta
*      LOCAL:  foreach scen in A_trial B_transport C_mrdr { do "analyses/transport_dvd/simulate.do" 1 1 500 "`scen'" }
*
*  (e) Cross-scenario aggregation (bootstrap mean + percentile 95% CI) -> results/  (the headline metric)
* do "analyses/transport_dvd/simulated/bootstrap_summary.do"
*
*  (f) (optional) Sample-size decomposition. Needs the bootstrap sims, so it sits here rather than
*      with the cohort build.  -> results/ce_sample_size.csv
* do "analyses/transport_dvd/patients/sample_size/ce_sample_size.do"


/* ================================================================================================
   HPC plumbing for steps (a)-(e) -- RUN BY HAND IN A VS CODE TERMINAL (bash), NOT in Stata.
   (This is a Stata block comment, so it is ignored when you `do` this file; select the lines you
   need and run them in a terminal.) FULL-data MI -> risk equations + calibrated transport -> 3x500
   scenario sims -> cross-scenario aggregation.

   # FIRST, in this terminal, load the machine-specific paths (git-ignored; sets user/repo_path/hpc/drive_path):
   source env.sh
   data=251128
   analysis=transport_dvd
   scratch=~/em76_scratch2/$user      # simulate.do writes bootstrap sims here ($scratch_path in hpc/config.do)

   # Directories on hpc (mirrors the repo tree + the data tree)
   ssh $hpc "mkdir -p em76/$user/data/$data/bootstrap"
   ssh $hpc "mkdir -p em76/$user/prep"
   ssh $hpc "mkdir -p em76/$user/prep/inputs"
   ssh $hpc "mkdir -p em76/$user/hpc"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/outcomes/A_trial/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/outcomes/B_transport/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/coefficients/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/patients"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/simulated/A_trial/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/simulated/B_transport/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/simulated/C_mrdr/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/results"

   # Config (HPC paths: sets repo_path, data_path, scratch_path) -> ~/em76/adam/config.do
   rsync -auvzce ssh "$repo_path"/hpc/config.do $hpc:~/em76/$user/config.do

   # Prep + engine code the bootstrap steps load
   rsync -auvzce ssh "$repo_path"/prep/multiple_imputation.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/prep/risk_equations.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/prep/inputs/ $hpc:~/em76/$user/prep/inputs/
   rsync -auvzce ssh "$repo_path"/core/ $hpc:~/em76/$user/core/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/simulate.do $hpc:~/em76/$user/analyses/$analysis/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/outcomes/ $hpc:~/em76/$user/analyses/$analysis/outcomes/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/patients/patients_${analysis}_2.dta $hpc:~/em76/$user/analyses/$analysis/patients/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/simulated/bootstrap_summary.do $hpc:~/em76_scratch2/$user/$analysis/simulated/

   # Data: raw MRDR Long (MI input). The bootstrap MI datasets are PRODUCED on hpc by step (a).
   rsync -auvzce ssh "$drive_path/data/$data/MRDR Long.dta" $hpc:~/em76/$user/data/$data/

   # Slurm scripts -> hpc/ (generic array jobs 1-500; per-analysis args via --export at submit)
   rsync -auvzce ssh "$repo_path"/hpc/multiple_imputation.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/risk_equations.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/simulate.script $hpc:~/em76/$user/hpc/

   # == Steps (a)-(c): submit as SEPARATE batches so you never queue more than one 500-array at a time
   #    (the HPC caps how many jobs you can submit at once). Send one, wait for it to finish
   #    (ssh $hpc "squeue -u $user"), then send the next. Order matters: MI FIRST -- both risk equations
   #    and calibrated transport read its B<b> datasets; those two are independent of each other.
   #    If the main-model bootstrap MI already exists on hpc, SKIP (a) and go straight to (b)/(c).

   # (a) Bootstrap MI -- 500-job array (IMP=10, full model):
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,IMP=10 hpc/multiple_imputation.script"
   #     ... wait until this array has finished before submitting (b) ...

   # (b) Bootstrap risk equations -- 500-job array (run after MI has finished). --parsable prints just the
   #     job id: note it as REJOB if you want to chain step (d) with a dependency.
   ssh $hpc "cd em76/$user ; sbatch --parsable --mail-user=$hpc_email --export=ALL,ANALYSIS=$analysis,COEFFS=$analysis,MINYR=2020,MAXYR=2025 hpc/risk_equations.script"

   # (c) Bootstrap calibrated transport -- 1 serial job, loops min_bs..max_bs internally (run after MI
   #     has finished; independent of (b)). Note the printed id as CTJOB if chaining step (d).
   #
   #     THE cd IS NOT OPTIONAL. calibrated_transport.do does `cd "$repo_path"' BEFORE it runs
   #     config.do, so on a fresh Stata session repo_path is still empty and the cd is skipped -
   #     every relative path then resolves from wherever the job happened to start. A --wrap job
   #     does not inherit the submit directory the way the hpc/*.script jobs do, so it started in
   #     $HOME and the writes failed with r(603). \$HOME survives both quoting layers.
   #
   #     RANGE: cover EVERY replicate, not just gap-fills. Replicate b pairs its BCR prediction with
   #     the risk-equation parameters from the SAME resample, so a rebuilt MI invalidates the old
   #     transport artefacts for all of them.
   #
   #     TIME: 500 iterations fitted inside 20 minutes in the July run; 30 is a small margin for the
   #     extra 12 and the changed extraction. Do not inflate it - a longer request queues longer.
   ssh $hpc "cd em76/$user ; sbatch --parsable --mail-user=$hpc_email --time=0-00:30:00 --wrap=\"cd \\\$HOME/em76/$user && module load stata/16u2021 && stata-mp analyses/$analysis/outcomes/calibrated_transport.do 1 1 512\""

   # == Step (d): 3x500 scenario sims -- ALSO one 500-array PER SCENARIO, so submit them ONE AT A TIME
   #    (run a single line below, wait, then the next). Each needs (b) AND (c) finished (their outputs
   #    must exist). If you'd rather queue them behind (b)/(c) instead of waiting, add
   #    --dependency=afterok:<REJOB>:<CTJOB> using the ids you noted above. ==
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,ANALYSIS=$analysis,SCENARIO=A_trial     hpc/simulate.script"
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,ANALYSIS=$analysis,SCENARIO=B_transport hpc/simulate.script"
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,ANALYSIS=$analysis,SCENARIO=C_mrdr      hpc/simulate.script"

   # == Step (e): aggregate. The sims land in SCRATCH; stage them into the analysis tree, then run
   #    bootstrap_summary on hpc (heavy sims stay on the cluster; pull only results/). ==
   for scen in A_trial B_transport C_mrdr ; dols

     ssh $hpc "rsync -a $scratch/$analysis/simulated/$scen/bootstrap/ em76/$user/analyses/$analysis/simulated/$scen/bootstrap/"
   done
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --time=0-00:30:00 --wrap=\"module load stata/16u2021 && stata-mp analyses/$analysis/simulated/bootstrap_summary.do\""
   rsync -auvzce ssh $hpc:~/em76_scratch2/$user/$analysis/results/ "$repo_path"/analyses/$analysis/results/

   # ---- FIRST-RUN CHECKS (two known path quirks; verify before trusting the aggregate) ----
   #  * calibrated_transport.do (boot branch) writes to RELATIVE  outcomes/{A_trial,B_transport}/bootstrap/,
   #    but the dispatcher's sim_bcr_override reads analyses/$analysis/outcomes/... -- confirm the sim step
   #    finds the B<b> transport coefficients (symlink/adjust the write path if not).
   #  * bootstrap_summary.do cd's to $repo_path then reads relative simulated/<scenario>/bootstrap/ -- confirm
   #    $sim resolves to analyses/$analysis/simulated (its header says "run from the analysis folder").
   ================================================================================================ */
