#!/bin/bash 
#SBATCH --partition=skylake
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --time=10:00:00

# . ~/miniconda3/etc/profile.d/conda.sh
# conda activate r-41
module activate R/4.2.0-icelake

Rscript src/plot-log-prob.R -m GT -i 1
Rscript src/plot-log-prob.R -m GT -i 2
Rscript src/plot-log-prob.R -m GT -i 3
Rscript src/plot-log-prob.R -m GT -i 4
Rscript src/plot-log-prob.R -m GT -i 5

Rscript src/plot-log-prob.R -m noGT -i 1
Rscript src/plot-log-prob.R -m noGT -i 2
Rscript src/plot-log-prob.R -m noGT -i 3
Rscript src/plot-log-prob.R -m noGT -i 4
Rscript src/plot-log-prob.R -m noGT -i 5
