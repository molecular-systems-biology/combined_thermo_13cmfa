#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --cpus-per-task=CPUS
#SBATCH --time=HOURS:02:00
#SBATCH --partition=regular
#SBATCH --mem-per-cpu=MEMORY
#SBATCH --signal=B:TERM@120

## Load Python and activate environment
module purge
source ~/venvs/13C_INCA/bin/activate

##line to be appended
