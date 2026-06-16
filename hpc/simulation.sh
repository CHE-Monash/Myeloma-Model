[M3 commands]
user=adam
analysis=transport_dvd
LocalPath=/Users/adami/Documents/Monash/Vault/research/models/"myeloma model"/repo/
M3Path=adami@m3-dtn.massive.org.au

# Make directories
ssh $M3Path "cd em76/$user/analyses/$analysis ; mkdir simulated"
ssh $M3Path "cd em76/$user/analyses/$analysis/simulated ; mkdir bootstrap"

# Send core
rsync -auvzce ssh $LocalPath/core $M3Path:~/em76/$user

# Send outcomes
rsync -auvzce ssh $LocalPath/analyses/$analysis/{outcomes,patients/patients_transport_dvd_2.dta} $M3Path:~/em76/$user/analyses/$analysis/

# Send analysis dispatcher
rsync -auvzce ssh $LocalPath/analyses/$analysis/$analysis.do $M3Path:~/em76/$user/analyses/$analysis

# Send slurm script
rsync -auvzce ssh $LocalPath/hpc/simulation.script $M3Path:~/em76/$user

# Run script
ssh $M3Path "cd em76/$user ; sbatch simulation.script"

# Send sample size 
rsync -auvzce ssh $LocalPath/analyses/$analysis/ce_sample_size.do $M3Path:~/em76/$user/analyses/$analysis/ce_sample_size.do

# Send bootstrap summary
rsync -auvzce ssh $LocalPath/analyses/$analysis/compare_scenarios.do $M3Path:~/em76/$user/analyses/$analysis

# Pull results
rsync -auvzce ssh $M3Path:~/em76/$user/analyses/$analysis/results/ $LocalPath/analyses/$analysis/results/





