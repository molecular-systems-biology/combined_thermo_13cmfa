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
BNAME=$5
NUMSAMPLES=$6
NUMDIRSAMPLES=$7
MAXSTEPS=$8
NUM_INIT_STEPS=$9
WARMUP=${10}

MINS=0

##old version no warmup here
#get output folder - current folder
OUTPUT_FOLDER=$(pwd)

# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f todolist.sh
touch todolist.sh

NAME="${CONDITION}_${BNAME}"

# Change the settings for the job submission:
cat ~/202601_13C_files/PTA_SCRIPTS/13C_pta_job_template.sh | sed "s/JOBNAME/$NAME/g" | sed "s/CPUS/$CPU/g" | sed "s/HOURS/$HOURS/g" | sed "s/MINS/$MINS/g" | sed "s/MEMORY/$MEMORY/g" > sample_script_${NAME}.slurm

# Append line to run the Python script with arguments:
echo "python -u ~/202601_13C_files/PTA_SCRIPTS/run_pta.py sample $CONDITION $BNAME $OUTPUT_FOLDER  \
  --findpoints \
--num_samples $NUMSAMPLES --num_dir_samples $NUMDIRSAMPLES --max_steps $MAXSTEPS --num_initial_steps $NUM_INIT_STEPS --max_threads $CPU --warmup $WARMUP \
  > ${OUTPUT_FOLDER}/log_${NAME}.log 2>&1" >> sample_script_$NAME.slurm 

echo "deactivate" >> sample_script_$NAME.slurm

echo "sbatch  sample_script_${NAME}.slurm" >> todolist.sh


# Turn the 'todolist' into an executable:
chmod +x todolist.sh