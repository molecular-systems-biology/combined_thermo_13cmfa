#!/bin/bash
#
# Script to setup pmo init pointssampling on the HPC cluster
#
# Nathan Fargier (template J.Losa)
# 2025.30.09


CPU=$1
HOURS=$2
MEMORY=$3
CONDITION=$4
# BNAME=$5
N_DRGSAMPLES=$5
N_INITPOINTS=$6
CHOSENSAMP_FILE=$7
FLUX_DF_FILE=$8

OUTPUT_FOLDER=$(pwd)

MINS=0

# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f todolist.sh
touch todolist.sh

# Read the chosen samples from CHOSENSAMP_FILE. This grabs the last
# numeric field on each line, so it works whether the file is a bare
# list of values, or a pandas Series dumped with an index/header
# (e.g. "0,12\n1,47\n...").
mapfile -t SAMPLES < <(awk -F',' 'NR>1 {print $2}' "$CHOSENSAMP_FILE")

if [ ${#SAMPLES[@]} -eq 0 ]; then
  echo "No numeric sample values found in $CHOSENSAMP_FILE, aborting." >&2
  exit 1
fi

for SAMPLE in "${SAMPLES[@]}"
do
  NAME="metsample_${CONDITION}_sample${SAMPLE}"

  # Change the settings for the job submission:
  cat ~/202601_13C_files/PTA_SCRIPTS/13C_pta_job_template.sh | sed "s/JOBNAME/$NAME/g" | sed "s/CPUS/$CPU/g" | sed "s/HOURS/$HOURS/g" | sed "s/MINS/$MINS/g" | sed "s/MEMORY/$MEMORY/g" > sample_script_$NAME.slurm

  # Append line to run the Python script with arguments:
  echo "python -u ~/202601_13C_files/METSAMPLE_PTA_SCRIPTS/run_pta_metsample_bootstrap.py $CONDITION $FLUX_DF_FILE $OUTPUT_FOLDER  \
    --n_drgsamples $N_DRGSAMPLES \
    --n_initpoint $N_INITPOINTS \
    --sample_n0 $SAMPLE \
    --max_threads $CPU \
    > ${OUTPUT_FOLDER}/log_${NAME}.log 2>&1" >> sample_script_$NAME.slurm

  echo "deactivate" >> sample_script_$NAME.slurm

  echo "sbatch  sample_script_${NAME}.slurm" >> todolist.sh

done

chmod +x todolist.sh