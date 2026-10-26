#!/bin/bash
#
# Script to setup simulations on the HPC cluster
### specific script to start classical 13C fitting 
# Nathan Fargier (template J.Losa)
# 2024.10.11

# Store the input arguments in variables with understandable names:
CPU=$1
HOURS=$2
MEMORY=$3
CONDITION=$4
LABELING_STRAT="combined"
MAX_ITERATIONS=$5
FLUX_FOLDER="~/202601_13C_files/INPUT_FILES/1_data/c_synthetic_data/2_fluxsamples_for13c/"
FLUX_FILE="${FLUX_FOLDER}/${CONDITION}_drawn_samples.csv"
##specific args for classic NLPCI
OPT_U_FILE=$6
#get output folder - current folder
OUTPUT_FOLDER=$(pwd)
# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f todolist.sh
touch todolist.sh

SLURM_MINS=$((MINS+1)) #add 1 min to job to allow process to save output 

# Setup batch scripts for each simulation:
for i in $(seq 0 203) #hardcoded n_reactions
do
    for dir in $(seq 0 1)
    do
  	NAME="${CONDITION}_flux${i}_direction${dir}_classicNLPCI"
        # Change the settings for the job submission:
        cat ~/202601_13C_files/FITTING_SCRIPTS/13C_fit_job_template.sh | sed "s/JOBNAME/$NAME/g" | sed "s/CPUS/$CPU/g" | sed "s/HOURS/$HOURS/g" | sed "s/MEMORY/$MEMORY/g" > script_$NAME.slurm

        # Append line to run the Python script with arguments:
        echo "python3 -u ~/202601_13C_files/FITTING_SCRIPTS/run_classic_NLPCI.py $CONDITION $LABELING_STRAT $OPT_U_FILE $i $dir 1e-6 $MAX_ITERATIONS $OUTPUT_FOLDER $NAME $FLUX_FILE  > ${OUTPUT_FOLDER}/log_${NAME}.log 2>&1" >> script_$NAME.slurm

        echo "deactivate" >> script_$NAME.slurm

        # Append job to the 'todolist':
        echo "sbatch  script_$NAME.slurm" >> todolist.sh
    done
done

# Turn the 'todolist' into an executable:
chmod +x todolist.sh




