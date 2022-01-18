#!/usr/bin/env bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p skylake
#SBATCH -t 100:00:00

module load R/4.0.3
Rscript run-GT.R -m sampling
