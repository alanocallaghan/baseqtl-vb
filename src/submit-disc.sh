#!/bin/bash 
#SBATCH --partition=skylake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --time=10:00:00

. ~/miniconda3/etc/profile.d/conda.sh
conda activate r-41

echo "running src/analyse-discrepant.R -m GT" && Rscript src/analyse-discrepant.R -m GT
echo "running src/analyse-discrepant.R -m noGT" && Rscript src/analyse-discrepant.R -m noGT
echo "done"
