#!/bin/bash 
#SBATCH --partition=skylake
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --time=50:00:00

# executable 
conda activate snakemake
snakemake -j 3 --profile slurm
