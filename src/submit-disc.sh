#!/bin/bash 
#SBATCH --partition=skylake
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --time=10:00:00

module load R/4.0.3
Rscript src/analyse-discrepant.R -m GT
Rscript src/analyse-discrepant.R -m noGT
