[M3 commands]
user=adam
data=251128
DrivePath=/Volumes/shared/R-MNHS-SPHPM-EPM-TRU/EpiMAP/Myeloma
LocalPath=/Users/adami/Documents/Monash/Vault/research/models/"myeloma model"/repo/
M3Path=adami@m3-dtn.massive.org.au

# Directories
ssh $M3Path "mkdir em76/$user"
ssh $M3Path "cd em76/$user ; mkdir data"
ssh $M3Path "cd em76/$user/data ; mkdir $data"
ssh $M3Path "cd em76/$user/data/$data ; mkdir bootstrap"

# Send to M3
rsync -auvzce ssh $DrivePath/multiple_imputation.do $M3Path:~/em76/$user
rsync -auvzce ssh $DrivePath/data/$data/"MRDR Long.dta" $M3Path:~/em76/$user/data/$data
rsync -auvzce ssh $LocalPath/hpc/multiple_imputation.script $M3Path:~/em76/$user

# Run script
ssh $M3Path "cd em76/$user ; sbatch multiple_imputation.script"

# Pull from M3
rsync -auvzce ssh $M3Path:~/em76/$user/data/$data/bootstrap $DrivePath/data/$data/