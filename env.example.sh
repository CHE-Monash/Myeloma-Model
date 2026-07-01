# env.example.sh — copy to env.sh and set your own paths.
#
# Bash counterpart of config.example.do. The real file (env.sh) is GIT-IGNORED (/env.sh, alongside
# /config.do), so your machine-specific paths never reach GitHub, while the committed run.do files
# reference only the variable names below. Names match config.do where the concept is shared
# ($repo_path, $drive_path).
#
# Usage: copy this file to env.sh, fill in your values, then in a VS Code terminal run
#     source env.sh
# before the rsync/ssh/sbatch lines in an analysis's run.do bootstrap block.

user=your_hpc_username                 # HPC account username
drive_path=/path/to/EpiMAP/Myeloma     # EpiMAP project base on your data drive (= $drive_path in config.do)
repo_path=/path/to/repo                # this repo on your machine        (= $repo_path  in config.do)
hpc=user@hpc-host.example.org          # ssh target for the cluster, so `ssh $hpc` works
hpc_email=you@example.org              # address for SLURM job notifications (passed as --mail-user at submit)
