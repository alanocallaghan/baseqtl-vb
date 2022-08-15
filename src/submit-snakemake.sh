#!/bin/bash 
#SBATCH --partition=skylake
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --time=10:00:00

# . ~/miniconda3/etc/profile.d/conda.sh
# conda activate
# snakemake -j 100 --profile slurm

# module activate R/4.2.0-icelake

echo "running Rscript src/plot-lm.R" && Rscript src/plot-lm.R
echo "running Rscript src/plot-times.R -m GT" && Rscript src/plot-times.R -m GT
echo "running Rscript src/plot-times.R -m noGT" && Rscript src/plot-times.R -m noGT
echo "running Rscript src/analyse-discrepant.R -m GT" && Rscript src/analyse-discrepant.R -m GT
echo "running Rscript src/analyse-discrepant.R -m noGT" && Rscript src/analyse-discrepant.R -m noGT

