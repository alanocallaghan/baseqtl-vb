#!/usr/bin/env bash

module load R/4.0.3

# Rscript src/plot-times.R -m GT -t 0.01
Rscript src/plot-times.R -m GT -t 0.001

# Rscript src/plot-times.R -m noGT -t 0.01
Rscript src/plot-times.R -m noGT -t 0.001

# Rscript src/plot.R -m vb -t 0.01
Rscript src/plot.R -m vb -t 0.001
Rscript src/plot.R -m sampling

# Rscript src/plot-noGT.R -m vb -t 0.01
Rscript src/plot-noGT.R -m vb -t 0.001
Rscript src/plot-noGT.R -m sampling
