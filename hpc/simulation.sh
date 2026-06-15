[M3 commands]
user=adam
analysis=dvd_method
LocalPath=/Users/adami/Documents/Monash/Vault/research/models/"myeloma model"/repo/
M3Path=adami@m3-dtn.massive.org.au

# Make directories
ssh $M3Path "cd em76/$user/analyses/$analysis ; mkdir simulated"
ssh $M3Path "cd em76/$user/analyses/$analysis/simulated ; mkdir bootstrap"

# Send core
rsync -auvzce ssh $LocalPath/core $M3Path:~/em76/$user

# Send outcomes
rsync -auvzce ssh $LocalPath/analyses/$analysis/{outcomes,patients} $M3Path:~/em76/$user/analyses/$analysis/

# Send analysis dispatcher
rsync -auvzce ssh $LocalPath/analyses/$analysis/$analysis.do $M3Path:~/em76/$user/analyses/$analysis

# Send slurm script
rsync -auvzce ssh $LocalPath/hpc/simulation.script $M3Path:~/em76/$user

# Run script
ssh $M3Path "cd em76/$user ; sbatch simulation.script"

# Send bootstrap summary
rsync -auvzce ssh $LocalPath/analyses/$analysis/compare_scenarios.do $M3Path:~/em76/$user/analyses/$analysis

# Pull bootstrap summary
rsync -auvzce ssh $M3Path:~/em76/$user/analyses/$analysis/results/ $LocalPath/analyses/$analysis/results/