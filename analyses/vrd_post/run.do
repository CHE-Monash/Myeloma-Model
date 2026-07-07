**********
* Monash Myeloma Model - vrd_post (VRd LoT 1 post-market): run.do (analysis runbook)
*
* Impact of VRd at line 1 vs historical practice, over the VRd-eligible predicted cohort, in order:
* risk equations -> simulate (x2 arms) -> validate. Run from the repository root. MRDR data via
* $data_path (config.do). simulate.do is the simulation dispatcher.
*
* Two arms, both simulated the same way (the dispatcher's $int, overridable as its 4th positional arg):
*   SoC   VRd-eligible patients receive their historical alternatives (VRd excluded)
*   VRd   VRd is available          -> the SoC-vs-VRd contrast is the post-market impact
*
* Depends on the SHARED prep outputs (built once by prep/, not per-analysis):
*   - ${data_path}/MRDR Long MI.dta            (prep/multiple_imputation.do 10 0)  -- risk-equation input
*   - analyses/vrd_post/patients/patients_vrd_l1_post.dta  -- the VRd-eligible L1 simulation cohort
*
* Two tracks:
*   DETERMINISTIC (point estimate)  -- light; runs locally, top to bottom (the live `do` lines).
*   BOOTSTRAP (prediction intervals) -- HEAVY. The risk-equation and 2x500-replicate simulation steps
*       run on the HPC. This file is NOT run on the HPC -- the bootstrap section below is the canonical
*       record of the sequence (left commented), with the equivalent sbatch commands and the
*       rsync/submit/pull plumbing in a /* ... */ shell block at the end.
**********

if "$repo_path" != "" cd "$repo_path"
capture run "config.do"


**********
* DETERMINISTIC (point estimate) -- runs locally, top to bottom
**********

* 1. Risk equations (VRd-excluded set) -> analyses/vrd_post/coefficients/coefficients_vrd_post.mmat
*    args: analysis coeffs min_year max_year boot [min_bs max_bs]   (loads outcomes/txr_vrd_post.do)
do "prep/risk_equations.do" vrd_post vrd_post 1995 2040 0

* 2. Simulate each arm (point estimate, $boot 0). The 4th positional selects the arm ($int);
*    min_bs/max_bs are unused at boot 0.
*    -> analyses/vrd_post/simulated/<arm>_1_predicted_1_999999.dta (+ validation)
foreach arm in SoC VRd {
    do "analyses/vrd_post/simulate.do" 0 1 1 "`arm'"
}


**********
* BOOTSTRAP (prediction intervals) -- HEAVY; runs on the HPC.
* The `do ...` lines below are the LOCAL (serial, slow) equivalent. The actual HPC plumbing -- rsync
* code+data up, sbatch the array jobs, pull results back -- is the /* ... */ shell block at the END of
* this file: it is a Stata block comment (inert here), so SELECT-AND-RUN those lines in a VS Code
* TERMINAL (bash). The slurm scripts (hpc/multiple_imputation.script, hpc/risk_equations.script,
* hpc/simulate.script) are generic array jobs 1-500 -- the per-analysis args ride on `sbatch --export=`.
**********

*  (a) SHARED bootstrap MI (FULL data) -> ${data_path}/bootstrap/MRDR Long MI B1..500.dta
*      This is the main-model bootstrap MI -- if base_model/transport_dvd already produced it, REUSE it.
*      LOCAL equivalent (slow, serial):
* do "prep/multiple_imputation.do" 10 1 1 500
*
*  (b) Bootstrap risk equations -> analyses/vrd_post/coefficients/bootstrap/coefficients_vrd_post_B1..500
*      LOCAL:  do "prep/risk_equations.do" vrd_post vrd_post 1995 2040 1 1 500
*
*  (c) 2x500 bootstrap simulations (one array per arm) -> simulated/bootstrap/<arm>_1_predicted_1_999999_B<b>.dta
*      LOCAL:  foreach arm in SoC VRd { do "analyses/vrd_post/simulate.do" 1 1 500 "`arm'" }
*
*  (d) Aggregate the SoC-vs-VRd bootstrap distribution into prediction intervals -- no aggregation
*      script exists for this analysis yet (cf. oos/bootstrap_validation.do, transport_dvd/
*      simulated/bootstrap_summary.do); add one when the PI summary is defined.


/* ================================================================================================
   HPC plumbing for steps (a)-(c) -- RUN BY HAND IN A VS CODE TERMINAL (bash), NOT in Stata.
   (This is a Stata block comment, so it is ignored when you `do` this file; select the lines you
   need and run them in a terminal.) FULL-data MI -> risk equations -> 2x500 arm sims.

   # FIRST, in this terminal, load the machine-specific paths (git-ignored; sets user/repo_path/hpc/drive_path):
   source env.sh
   data=251128
   analysis=vrd_post

   # Directories on hpc (mirrors the repo tree + the data tree)
   ssh $hpc "mkdir -p em76/$user/data/$data/bootstrap"
   ssh $hpc "mkdir -p em76/$user/prep"
   ssh $hpc "mkdir -p em76/$user/hpc"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/outcomes"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/coefficients/bootstrap"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/patients"
   ssh $hpc "mkdir -p em76/$user/analyses/$analysis/simulated/bootstrap"

   # Config (HPC paths: sets repo_path, data_path) -> ~/em76/adam/config.do
   rsync -auvzce ssh "$repo_path"/hpc/config.do $hpc:~/em76/$user/config.do

   # Prep + engine code the bootstrap steps load
   rsync -auvzce ssh "$repo_path"/prep/multiple_imputation.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/prep/risk_equations.do $hpc:~/em76/$user/prep/
   rsync -auvzce ssh "$repo_path"/core/ $hpc:~/em76/$user/core/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/simulate.do $hpc:~/em76/$user/analyses/$analysis/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/outcomes/txr_vrd_post.do $hpc:~/em76/$user/analyses/$analysis/outcomes/
   rsync -auvzce ssh "$repo_path"/analyses/$analysis/patients/patients_vrd_l1_post.dta $hpc:~/em76/$user/analyses/$analysis/patients/

   # Data: raw MRDR Long (MI input). The bootstrap MI datasets are PRODUCED on hpc by step (a).
   rsync -auvzce ssh "$drive_path/data/$data/MRDR Long.dta" $hpc:~/em76/$user/data/$data/

   # Slurm scripts -> hpc/ (generic array jobs 1-500; per-analysis args via --export at submit)
   rsync -auvzce ssh "$repo_path"/hpc/multiple_imputation.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/risk_equations.script $hpc:~/em76/$user/hpc/
   rsync -auvzce ssh "$repo_path"/hpc/simulate.script $hpc:~/em76/$user/hpc/

   # == Steps (a)-(b): MI -> risk equations, chained on the MI array (afterok). MI: IMP=10, SAMPLE empty. ==
   #    If the main-model bootstrap MI already exists on hpc, SKIP the MI submit and run only risk equations.
   ssh $hpc "cd em76/$user ; \
     mi=\$(sbatch --parsable --mail-user=$hpc_email --export=ALL,IMP=10 hpc/multiple_imputation.script) ; echo MI job \$mi ; \
     re=\$(sbatch --parsable --mail-user=$hpc_email --dependency=afterok:\$mi --export=ALL,ANALYSIS=$analysis,COEFFS=$analysis,MINYR=1995,MAXYR=2040 hpc/risk_equations.script) ; echo risk-eq job \$re ; \
     echo \"export REJOB=\$re\"    # note this job id for step (c)'s dependency"

   # == Step (c): 2x500 arm sims. One simulate array per arm, each waiting on the risk-equation array. ==
   #    The generic simulate.script passes SCENARIO as the dispatcher's 4th positional -> vrd_post maps it to $int.
   #    Fill REJOB from step (b)'s echo, then:
   for arm in SoC VRd ; do
     ssh $hpc "cd em76/$user ; sbatch --mail-user=$hpc_email --dependency=afterok:$REJOB \
       --export=ALL,ANALYSIS=$analysis,SCENARIO=$arm hpc/simulate.script"
   done

   # Pull the bootstrap sims back into the repo (or aggregate on hpc once a PI summary script exists)
   rsync -auvzce ssh $hpc:~/em76/$user/analyses/$analysis/simulated/bootstrap/ "$repo_path"/analyses/$analysis/simulated/bootstrap/
   ================================================================================================ */
