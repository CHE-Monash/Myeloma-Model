[M3 commands]
user=adam
analysis=dvd_method
DrivePath=/Volumes/shared/R-MNHS-SPHPM-EPM-TRU/EpiMAP/Myeloma/
LocalPath=/Users/adami/Documents/Monash/Vault/research/models/"myeloma model"/repo/
M3Path=adami@m3-dtn.massive.org.au

# Directories
ssh $M3Path "mkdir em76/$user"
ssh $M3Path "cd em76/$user ; mkdir analyses"
ssh $M3Path "cd em76/$user/analyses ; mkdir $analysis"
ssh $M3Path "cd em76/$user/analyses/$analysis ; mkdir coefficients"
ssh $M3Path "cd em76/$user/analyses/$analysis/coefficients ; mkdir bootstrap"

#ssh $M3Path "cd em76/$user/analyses ; mkdir transport"
#ssh $M3Path "cd em76/$user/analyses/$analysis/transport ; mkdir bootstrap"

# Send files
rsync -auvzce ssh $DrivePath/risk_equations.do $M3Path:~/em76/$user
rsync -auvzce ssh $DrivePath/analyses/$analysis/txr_$analysis.do $M3Path:~/em76/$user/analyses/$analysis

rsync -auvzce ssh $DrivePath/analyses/$analysis/txr_dvd_pre.do $M3Path:~/em76/$user/analyses/$analysis
rsync -auvzce ssh $DrivePath/analyses/$analysis/txr_dvd_post.do $M3Path:~/em76/$user/analyses/$analysis

rsync -auvzce ssh $DrivePath/analyses/$analysis/transport.do $M3Path:~/em76/$user/analyses/$analysis

rsync -auvzce ssh $LocalPath/hpc/risk_equations.script $M3Path:~/em76/$user

rsync -auvzce ssh $DrivePath/analyses/$analysis/$analysis.do $M3Path:~/em76/$user/

# Run script
ssh $M3Path "cd em76/$user ; sbatch risk_equations.script"

# Pull files
coefficient_path=analyses/$analysis/coefficients
rsync -auvzce ssh $M3Path:~/em76/$user/$coefficient_path/bootstrap $LocalPath/Simulation/$coefficient_path