#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --cpus-per-task=CPUS
#SBATCH --time=HOURS:MINS:00
#SBATCH --partition=regular
#SBATCH --mem-per-cpu=MEMORY

## Load Python and activate environment
module purge
module load Gurobi/12.0.0-GCCcore-13.2.0
module load R/4.4.1-gfbf-2023b
source ~/venvs/pta/bin/activate


##lines to be appended

