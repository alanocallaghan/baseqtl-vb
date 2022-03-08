#!/usr/bin/env bash

sbatch src/submit-gtvb.sh
# sbatch src/submit-gthmc.sh

sbatch src/submit-nogtvb.sh
# sbatch src/submit-nogthmc.sh
