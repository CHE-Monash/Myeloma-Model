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
*    -> analyses/oos/simulated/all_0_oos.dta
do "analyses/oos/simulate.do"

* 6. Validate the point estimate vs the observed targets (fixed tolerances)
do "analyses/oos/validate_oos.do"


**********
* BOOTSTRAP (prediction intervals) -- the headline OOS metric; HEAVY; runs on the HPC.
* The `do ...` lines below are the LOCAL (serial, slow) equivalent. The actual HPC plumbing -- rsync
* code+data up, sbatch the array jobs, pull results back -- is the /* ... */ shell block at the END of
* this file: it is a Stata block comment (inert here), so SELECT-AND-RUN those lines in a VS Code
* TERMINAL (bash). The slurm scripts (hpc/multiple_imputation.script, hpc/risk_equations.script) are
* generic array jobs 1-500 -- the per-analysis args ride on the `sbatch --export=` lines below.
**********

*  (a) TRAIN bootstrap: MI + risk equations -> analyses/oos/coefficients/bootstrap/coefficients_oos_B1..500
*      LOCAL equivalent (slow, serial):
* do "prep/multiple_imputation.do" 2 1 1 500 train
* do "prep/risk_equations.do" oos oos 1995 2040 1 1 500 train
*      HPC: the two generic array jobs (multiple_imputation.script then risk_equations.script) -- see the shell block below.
*
*  (b) 500 bootstrap simulations of the held-out 30% -> analyses/oos/simulated/bootstrap/
*      LOCAL:  do "analyses/oos/simulate.do" 1 1 500
*      HPC: an array job over simulate.do 1 <task> <task> via hpc/simulate.script -- see the shell block below.
*
*  (c) Percentile-method PI coverage vs the held-out observed -> analyses/oos/results/  (the headline metric)
* do "analyses/oos/bootstrap_validation.do" 500


/* ================================================================================================
   HPC plumbing for steps (a)-(c) -- RUN BY HAND IN A VS CODE TERMINAL (bash), NOT in Stata.
   (This is a Stata block comment, so it is ignored when you `do` this file; select the lines you
   need and run them in a terminal.) Train-fold MI -> train-fold risk equations on MASSIVE hpc ->
   analyses/oos/coefficients/bootstrap/coefficients_oos_B1..500.
   Prereq (on the Mac, needs the drive): do analyses/oos/prep/oos_split.do -> ${data_path}/oos/oos_split.dta

   # FIRST, in this terminal, load the machine-specific paths (git-ignored; sets user/repo_path/hpc/drive_path):
   source env.sh
   data=251128
   analysis=oos

   # Directories on hpc (mirrors the repo tree + the data tree)
   ssh $hpc "mkdir -p em76/$user/data/$data/oos/bootstrap"
   ssh $hpc "mkdir -p em76/$user/prep"
   ssh $hpc "mkdir -p em76/$user/hpc"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/outcomes"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/coefficients/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/patients"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/simulated/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/results"

   # Config (HPC paths: sets repo_path, data_path) -> ~/em76/adam/config.do
   rsync -auvzce ssh "$repo_path"/hpc/config.do $hpc:~/em76/$user/config.do

   # Prep code + the OOS regimen definition risk_equations loads
   rsync -auvzce ssh "$repo_path"/prep/multiple_imputation.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/prep/risk_equations.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/outcomes/txr_oos.do $hpc:~/em76/$user/analyses/$analysis/outcomes/

   # Data: MRDR Long + the 70/30 split crosswalk
   rsync -auvzce ssh "$drive_path/data/$data/MRDR Long.dta" $hpc:~/em76/$user/data/$data/
   rsync -auvzce ssh "$drive_path/data/$data/oos/oos_split.dta" $hpc:~/em76/$user/data/$data/oos/

   # Slurm scripts -> hpc/ (mirrors the repo; generic array jobs 1-500, per-analysis args via --export at submit)
   rsync -auvzce ssh "$repo_path"/hpc/multiple_imputation.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/risk_equations.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/simulate.script $hpc:~/em76/$user/hpc/

   # == Step (a): submit MI + risk equations. Run ONE of [1]/[2], or [3] to chain both. ==
   #    --export args:  MI needs IMP, SAMPLE  |  risk equations needs ANALYSIS, COEFFS, MINYR, MAXYR, SAMPLE.

   # [1] MI only (the 500-resample bootstrap MI array):
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,IMP=2,SAMPLE=train hpc/multiple_imputation.script"

   # [2] Risk equations only (run after the MI array has finished):
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,ANALYSIS=oos,COEFFS=oos,MINYR=1995,MAXYR=2040,SAMPLE=train hpc/risk_equations.script"

   # [3] Both chained -- risk equations waits on the whole MI array (afterok):
   ssh $hpc "cd em76/$user ; mi=\$(sbatch --parsable --mail-user=$hpc_email --export=ALL,IMP=2,SAMPLE=train hpc/multiple_imputation.script) ; echo MI job \$mi ; sbatch --mail-user=$hpc_email --dependency=afterok:\$mi --export=ALL,ANALYSIS=oos,COEFFS=oos,MINYR=1995,MAXYR=2040,SAMPLE=train hpc/risk_equations.script"

   # Pull the bootstrap coefficients back into the repo (optional; they are already on hpc for step b)
   rsync -auvzce ssh $hpc:~/em76/$user/analyses/$analysis/coefficients/bootstrap/ "$repo_path"/analyses/$analysis/coefficients/bootstrap/
   
   # == Step (b): 500 bootstrap simulations of the held-out 30% (run after step a's coefficients are on hpc) ==
   # Reads coefficients_oos_B1..500 (already on hpc) + the held-out cohort; writes analyses/oos/simulated/bootstrap/.
   # Needs core/ + simulate.do + the cohort, which the MI/risk phase did not send:
   rsync -auvzce ssh "$repo_path"/core/ $hpc:~/em76/$user/core/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/simulate.do $hpc:~/em76/$user/analyses/$analysis/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/patients/oos_cohort.dta $hpc:~/em76/$user/analyses/$analysis/patients/

   # Submit the simulation array (simulate.do 1 <task> <task>; --export picks the analysis):
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --export=ALL,ANALYSIS=oos hpc/simulate.script"

   # == Step (c): PI coverage. Aggregate the 500 sims on hpc (heavy -- don't pull them); pull only results/. ==
   # Needs the target CSVs + the validator on hpc (verify bootstrap_validation.do's own file deps on first run):
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/targets/ $hpc:~/em76/$user/analyses/$analysis/targets/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/bootstrap_validation.do $hpc:~/em76/$user/analyses/$analysis/
   ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --time=0-00:30:00 --wrap=\"module load stata/16u2021 && stata-mp analyses/$analysis/bootstrap_validation.do 500\""
   rsync -auvzce ssh $hpc:~/em76/$user/analyses/$analysis/results/ "$repo_path"/analyses/$analysis/results/
   ================================================================================================ */
