**********
* Monash Myeloma Model - OOS (70/30): run.do (analysis runbook)
*
* The full OOS pipeline in order -- prep -> simulate -> validate. Run from the repository root.
* MRDR data is read via $data_path (config.do). simulate.do is the simulation dispatcher.
*
* Two tracks:
*   DETERMINISTIC (point estimate)  -- light; runs locally, top to bottom (the live `do` lines).
*   BOOTSTRAP (prediction intervals) -- the headline OOS metric; HEAVY. The MI, risk-equation and
*       500-replicate simulation steps run on the HPC. This file is NOT run on the HPC -- the
*       bootstrap section below is the canonical record of the sequence (left commented), with the
*       equivalent sbatch commands and the rsync/submit/pull plumbing in a /* ... */ shell block at the end.
*
* The split (step 0) and per-fold imputation need the restricted MRDR drive; everything from the
* targets/cohort down operates on the imputed outputs.
**********

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

* 0. Split patients 70/30 -> ${data_path}/oos/oos_split.dta   (run ONCE, fixed seed; needs MRDR drive)
do "analyses/oos/prep/oos_split.do"


**********
* DETERMINISTIC (point estimate) -- runs locally, top to bottom
**********

* 1. Multiple imputation, each fold imputed SEPARATELY (split BEFORE imputation = no leakage).
*    args: imp boot min_bs max_bs sample
do "prep/multiple_imputation.do" 2 0 . . train         // train (70%) -> coefficients
do "prep/multiple_imputation.do" 2 0 . . test          // test  (30%) -> targets + cohort

* 2. Risk equations on the TRAIN fold -> analyses/oos/coefficients/coefficients_oos.mmat
*    args: analysis coeffs min_year max_year boot min_bs max_bs sample   (loads outcomes/txr_oos.do)
do "prep/risk_equations.do" oos oos 1995 2040 0 . . train

* 3. Held-out 30% validation targets (observed outcomes) -> analyses/oos/targets/ (13 csv)
do "analyses/oos/prep/oos_targets.do"

* 4. Held-out 30% simulation cohort -> analyses/oos/patients/oos_cohort.dta
do "analyses/oos/prep/oos_cohort.do"

* 5. Simulate the 30% with the 70%-trained coefficients ($boot 0)
*    -> analyses/oos/simulated/all_0_oos_1_101212.dta
do "analyses/oos/simulate.do"

* 6. Validate the point estimate vs the observed targets (fixed tolerances)
do "analyses/oos/validate_oos.do"


**********
* BOOTSTRAP (prediction intervals) -- the headline OOS metric; HEAVY; runs on the HPC.
* The `do ...` lines below are the LOCAL (serial, slow) equivalent. The actual HPC plumbing -- rsync
* code+data up, sbatch the array jobs, pull results back -- is the /* ... */ shell block at the END of
* this file: it is a Stata block comment (inert here), so SELECT-AND-RUN those lines in a VS Code
* TERMINAL (bash). The slurm scripts (hpc/oos_mi.script, hpc/oos_risk_equations.script) are array jobs 1-500.
**********

*  (a) TRAIN bootstrap: MI + risk equations -> analyses/oos/coefficients/bootstrap/coefficients_oos_B1..500
*      LOCAL equivalent (slow, serial):
* do "prep/multiple_imputation.do" 2 1 1 500 train
* do "prep/risk_equations.do" oos oos 1995 2040 1 1 500 train
*      HPC: the two array jobs (oos_mi.script then oos_risk_equations.script) -- see the shell block below.
*
*  (b) 500 bootstrap simulations of the held-out 30% -> analyses/oos/simulated/bootstrap/
*      LOCAL:  do "analyses/oos/simulate.do" 1 1 500
*      HPC: an array job over simulate.do 1 <task> <task> (slurm wrapper mirrors oos_mi.script; not yet created).
*
*  (c) Percentile-method PI coverage vs the held-out observed -> analyses/oos/results/  (the headline metric)
* do "analyses/oos/bootstrap_validation.do" 500


/* ================================================================================================
   HPC plumbing for steps (a)/(b) -- RUN BY HAND IN A VS CODE TERMINAL (bash), NOT in Stata.
   (This is a Stata block comment, so it is ignored when you `do` this file; select the lines you
   need and run them in a terminal.) Train-fold MI -> train-fold risk equations on MASSIVE M3 ->
   analyses/oos/coefficients/bootstrap/coefficients_oos_B1..500.
   Prereq (on the Mac, needs the drive): do analyses/oos/prep/oos_split.do -> ${data_path}/oos/oos_split.dta

   # FIRST, in this terminal, load the machine-specific paths (git-ignored; sets user/repo_path/hpc/drive_path):
   source hpc/env.sh
   data=251128
   analysis=oos

   # Directories on M3 (data tree + repo tree)
   ssh $hpc "mkdir -p em76/$user/data/$data/oos/bootstrap"
   ssh $hpc "mkdir -p em76/$user/prep"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/outcomes"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/coefficients/bootstrap"

   # Config (HPC paths: sets repo_path, data_path) -> ~/em76/adam/config.do
   rsync -auvzce ssh "$repo_path"/hpc/config.do $hpc:~/em76/$user/config.do

   # Prep code + the OOS regimen definition risk_equations loads
   rsync -auvzce ssh "$repo_path"/prep/multiple_imputation.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/prep/risk_equations.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/outcomes/txr_oos.do $hpc:~/em76/$user/analyses/$analysis/outcomes/

   # Data: MRDR Long + the 70/30 split crosswalk
   rsync -auvzce ssh "$drive_path/data/$data/MRDR Long.dta" $hpc:~/em76/$user/data/$data/
   rsync -auvzce ssh "$drive_path/data/$data/oos/oos_split.dta" $hpc:~/em76/$user/data/$data/oos/

   # Slurm scripts (array jobs 1-500)
   rsync -auvzce ssh "$repo_path"/hpc/oos_mi.script $hpc:~/em76/$user/
   rsync -auvzce ssh "$repo_path"/hpc/oos_risk_equations.script $hpc:~/em76/$user/

   # Submit: risk equations starts only after the whole MI array finishes OK
   ssh $hpc "cd em76/$user ; mi=\$(sbatch --parsable oos_mi.script) ; echo MI job \$mi ; sbatch --dependency=afterok:\$mi oos_risk_equations.script"

   # Monitor
   # ssh $hpc "squeue -u $user"

   # Pull the bootstrap coefficients back into the repo
   rsync -auvzce ssh $hpc:~/em76/$user/analyses/$analysis/coefficients/bootstrap/ "$repo_path"/analyses/$analysis/coefficients/bootstrap/
   ================================================================================================ */
