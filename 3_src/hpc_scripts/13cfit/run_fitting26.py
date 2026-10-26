#-------------------------------------------------------------------------------
# Name:        run_fitting.py

# Purpose:     build c13problem class and run script on HPC (folders defined here)
#              with arg parse the script is flexible so we can use this script to run either of the labeling strategies
# 
# 
# Author:      Nathan Fargier
#
# Created:     22/01/2025

#--------
"""Inputs :
-condition=string
labeling_strategy=string
sample_no=int
tol=1e-2
max_iters=1000
n_jobs = int 
dataframe flux sample file
folder MS data labeling 

Outputs : 
-optimal objective
-optimal freeflux vector (u)
-optimization successful or not
"""

import sys
import os 

sys.path.append(os.path.expanduser('~/'))
sys.path.append('/home4/p313427/')
sys.path.append('/home2/p313417/')
sys.path.append('/home6/p320093/')

import os
import numpy as np
from HELPER_FREEFLUX_NF import c13problem_class as c13
from datetime import datetime
import signal
import pandas as pd

### deactivate warnings 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse

parser = argparse.ArgumentParser(description="Run flux fitting with specified parameters.")

# required positional arguments
parser.add_argument("condition", type=str)
parser.add_argument("sample_no", type=int)
parser.add_argument("BASE_FOLDER", type=str) ##folder with data
parser.add_argument("OUTPUT_FOLDER", type=str)
parser.add_argument("FLUX_FILE", type=str)

# optional arguments
parser.add_argument("--solver", type=str,default ='gsubg' )
parser.add_argument("--tol", type=float, default=1e-6, help='tolerance for solution acceptance')
parser.add_argument("--max_iters", type=int, default=200 ,help='n iteration in a local solve')
parser.add_argument("--n_restart", type=int, default=3 ,help='number of local solve restarts')
parser.add_argument("--n_jobs", type=int, default=8,help='number of processes for the initial calculations on the network')
parser.add_argument('--labeling_strategy',type=str,default='combined',help='should be 1-13C, U-13C or combined')
parser.add_argument('--dataset',type=str,default='mixed',help='should be MAREIKE or MHLAB or mixed')

# boolean flags (store True if provided) ## shuold make these default as true
parser.add_argument("--flux_bounds", action="store_true",help='use flux bounds')
parser.add_argument("--fit_mdv", action="store_true",help='fit 13c mdv if true, else fit ex flux values')
parser.add_argument("--force_flux", action="store_true",help='add FF constraints if true')

#seed for bootstrapping 13c data
def seed_type(s):
    if s.lower() == 'none':
        return None
    return int(s)

parser.add_argument("--seed", type=seed_type, default=None, help='seed for bootstrapping 13c data')

args = parser.parse_args()

print(vars(args))##print args

INPUT_FOLDER = f"{args.BASE_FOLDER}1_data/b_processed_files/"
SYNDATA_FOLDER= f"{args.BASE_FOLDER}1_data/c_synthetic_data/"
MODEL_FILE = f"{args.BASE_FOLDER}1_data/c_synthetic_data/reactions_processed_ATP_2601.xlsx"

OUTPUT_FILENAME=f'{args.condition}_{args.labeling_strategy}_sample{args.sample_no}'
##init log 
today = datetime.now().strftime("%Y%m%d")
os.makedirs(args.OUTPUT_FOLDER, exist_ok = True)

import os, uuid

# Prefer a scheduler-provided unique ID if you have one (SLURM, PBS, etc.)
JOB_ID = (
    os.environ.get('SLURM_JOB_ID')
    or os.environ.get('SLURM_ARRAY_JOB_ID', '') + '_' + os.environ.get('SLURM_ARRAY_TASK_ID', '')
    or str(os.getpid())
)

FLAG_DIR = os.path.expanduser('~/HELPER_FREEFLUX_NF/flags')
os.makedirs(FLAG_DIR, exist_ok=True)
FLAG_FILE = os.path.join(FLAG_DIR, f'exit_flag_{JOB_ID}.txt')

def signal_handler(signum, frame):
    print(f"[Signal Handler] Caught signal {signum}, setting flag ({FLAG_FILE})")
    with open(FLAG_FILE, 'w') as f:
        f.write("TRUE")

if os.path.exists(FLAG_FILE):
    try:
        os.remove(FLAG_FILE)
    except FileNotFoundError:
        pass

def main():
    c13_problem=c13.c13_problem(args.condition,args.labeling_strategy,args.sample_no,MODEL_FILE,args.FLUX_FILE,INPUT_FOLDER,args.OUTPUT_FOLDER,OUTPUT_FILENAME,
                                tol=args.tol,max_iters=args.max_iters,solver=args.solver,force_flux=args.force_flux,
                                n_restart=args.n_restart,n_jobs=args.n_jobs,flux_bounds=args.flux_bounds,fit_mdv=args.fit_mdv,dataset=args.dataset,
                                FLAG_FILE=FLAG_FILE,seed=args.seed)

    signal.signal(signal.SIGTERM, signal_handler) ## can catch signal to stop the process in case of timelimit


    if args.force_flux:
        c13_problem.drg_df=pd.read_csv(f'{SYNDATA_FOLDER}3_drGVA_indivfluxbounds/{args.condition}_indivfluxbounds_drGbounds.csv',index_col=[0,1])
        
    if c13_problem.labeling_strategy=='combined':
        c13_problem.labeling_strategy='1-13C'
        c13_problem._set_fit()
        c13_problem.optModelA=c13_problem._build_optmodel()
        #then second model
        c13_problem.modelB=c13_problem._set_model()
        c13_problem.labeling_strategy='U-13C'
        c13_problem._set_fit(alternative_model=True)
        c13_problem.optModelB=c13_problem._build_optmodel(alternative_model=True)
        ##if we want to save the nullspace matrix on cluster : 
        # pd.DataFrame(c13_problem.modelA.null_space).to_csv(f'{c13_problem.OUTPUT_FOLDER}/nullspace_{c13_problem.OUTPUT_FILENAME}.csv')
        # now c13problem has the 2 model, we solve it using combined f and df in optim
        opt_obj, opt_u, is_success= c13_problem()

        print(f'>>>>>> Optimization successful : {is_success} with objective {opt_obj}')

        ## save the array in a csv file
        #os.makedirs(OUTPUT_FOLDER, exist_ok = True)
        np.savetxt(f'{c13_problem.OUTPUT_FOLDER}/optu_{c13_problem.OUTPUT_FILENAME}.csv',opt_u,delimiter=',',header=str(opt_obj))
    else:
        c13_problem._set_fit()
        c13_problem.optModelA=c13_problem._build_optmodel()
        c13_problem.optModelB=None
        opt_obj, opt_u, is_success= c13_problem()
        print(f'>>>>>> Optimization successful : {is_success} with objective {opt_obj}')

        ## save the array in a csv file
        #os.makedirs(OUTPUT_FOLDER, exist_ok = True)
        np.savetxt(f'{c13_problem.OUTPUT_FOLDER}/optu_{c13_problem.OUTPUT_FILENAME}.csv',opt_u,delimiter=',',header=str(opt_obj))

main()
