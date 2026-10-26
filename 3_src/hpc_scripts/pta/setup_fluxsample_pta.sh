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
NUMSAMPLES=$5
NUMINPUTDRGSAMP=$6
PTA_FOLDER=$7
MINS=0

##old version no warmup here
#get output folder - current folder
OUTPUT_FOLDER=$(pwd)

NAME="fluxsample_${CONDITION}"

# Change the settings for the job submission:
cat ~/202601_13C_files/PTA_SCRIPTS/13C_pta_job_template.sh | sed "s/JOBNAME/$NAME/g" | sed "s/CPUS/$CPU/g" | sed "s/HOURS/$HOURS/g" | sed "s/MINS/$MINS/g" | sed "s/MEMORY/$MEMORY/g" > sample_script_$NAME.slurm

# Append line to run the Python script with arguments:
echo "python -u ~/202601_13C_files/PTA_SCRIPTS/pta_fluxsample.py $CONDITION $PTA_FOLDER $OUTPUT_FOLDER  \
  --num_samples $NUMSAMPLES \
  --num_inputdrG_samples $NUMINPUTDRGSAMP \
  > ${OUTPUT_FOLDER}/log_${NAME}.log 2>&1" >> sample_script_$NAME.slurm 
 
echo "deactivate" >> sample_script_$NAME.slurm

##then will be run without todolist