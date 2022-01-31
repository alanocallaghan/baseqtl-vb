#!/usr/bin/env bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p skylake
#SBATCH -t 10:00:00

module load R/4.0.3
# Rscript run-GT.R -m sampling
Rscript run-noGT.R -m optimizing
Rscript run-noGT.R -m vb
Rscript run-noGT.R -m sampling
