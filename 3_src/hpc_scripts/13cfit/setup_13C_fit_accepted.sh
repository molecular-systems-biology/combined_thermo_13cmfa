#!/bin/bash
#
# Script to setup simulations on the HPC cluster
#
# Nathan Fargier (template J.Losa)
# 2026.01.09

# Store the input arguments in variables with understandable names:
CPU=$1
HOURS=$2
MEMORY=$3
CONDITION=$4
DATASET="mixed"
SAMPLE=$5
N_RESTARTS=$6
###solvers should be a list like "mma auglag"
MAX_ITERS=$7
LABELING_STRAT=$8 ##"combined" for both 

#forceflux arg true or false 
FORCE_FLUX=$9
FORCE_FLUX_ARG=""
if [ "$FORCE_FLUX" = "True" ]; then
    FORCE_FLUX_ARG="--force_flux"
fi
ACC_SAMPLES_FILE=${10}

##seed for bootstrapping of 13C data
SEED=${11}

BASE_FOLDER="~/202601_13C_files/INPUT_FILES/"
FLUX_FOLDER="~/202601_13C_files/INPUT_FILES/1_data/c_synthetic_data/2_fluxsamples_for13c/"
FLUX_FILE="${FLUX_FOLDER}/${CONDITION}_drawn_samples.csv"
#get output folder - current folder
OUTPUT_FOLDER=$(pwd)
# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f todolist.sh
touch todolist.sh

SLURM_MINS=$((MINS+1)) #add 1 min to job to allow process to save output 

# Setup batch scripts for each simulation:

while read -r i
do
if [ ! -f "${OUTPUT_FOLDER}/optu_${CONDITION}_combined_sample${i}.csv" ]
then
        echo "File does NOT exist → creating job for sample $i"

        NAME="${CONDITION}_sample${i}"
        # Change the settings for the job submission:
        cat ~/202601_13C_files/FITTING_SCRIPTS/13C_fit_job_template.sh | sed "s/JOBNAME/$NAME/g" | sed "s/CPUS/$CPU/g" | sed "s/HOURS/$HOURS/g" | sed "s/MEMORY/$MEMORY/g" > script_$NAME.slurm

        # Append line to run the Python script with arguments:
        echo "exec python -u $HOME/202601_13C_files/FITTING_SCRIPTS/run_fitting26.py $CONDITION  $i $BASE_FOLDER $OUTPUT_FOLDER $FLUX_FILE \
        --max_iters $MAX_ITERS \
        --n_restart $N_RESTARTS \
        --labeling_strategy $LABELING_STRAT \
        --flux_bounds \
        --fit_mdv \
        $FORCE_FLUX_ARG \
        --seed $SEED \
        > ${OUTPUT_FOLDER}/log_${NAME}.log 2>&1" >> script_$NAME.slurm 
        # Append job to the 'todolist':
        echo "sbatch  script_$NAME.slurm" >> todolist.sh
fi
done <  "$ACC_SAMPLES_FILE"
# Turn the 'todolist' into an executable:
chmod +x todolist.sh
