#!/usr/bin/env bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p skylake
#SBATCH -t 10:00:00
#SBATCH --cpus-per-task 16
#SBATCH --ntasks 1

. ~/miniconda3/etc/profile.d/conda.sh
conda activate pyro

# python3 ./gt_nb_pyro.py
python3 ./gt_nb_pymc3.py
