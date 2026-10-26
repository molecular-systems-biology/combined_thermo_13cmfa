#-------------------------------------------------------------------------------
# Name:        run_fitting.py

# Purpose:     build c13problem class and run script on HPC (folders defined here)
# 
# 
# Author:      Nathan Fargier
#
# Created:     22/01/2025

#--------
"""Inputs (in this order):
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
sys.path.append('/home4/p313427/freeflux/')
sys.path.append('/home4/p313427/')
#sys.path.append(os.path.abspath('../'))
import signal 
import os
import numpy as np
from HELPER_FREEFLUX_NF import classic_c13problem as classic_c13
from datetime import datetime

### deactivate warnings 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
###
BASE_FOLDER=os.path.expanduser('~/202601_13C_files/INPUT_FILES/1_data/')
INPUT_FOLDER=f'{BASE_FOLDER}b_processed_files/'
MODEL_FILE=f'{BASE_FOLDER}c_synthetic_data/reactions_processed_ATP.xlsx'

arguments=sys.argv[1:]
##check that there are 10 arguments
if len(arguments)!=6:
    raise ValueError('6 arguments are needed') 
condition,labeling_strategy,dataset,n_restart,OUTPUT_FOLDER,FLUX_FILE=arguments
n_jobs=6
#sample_no=int(sample_no)
n_restart=int(n_restart)
tol=1e-6
max_iters=200
flux_bounds=True # if flux_bounds=='True' else flux_bounds
fit_mdv=True # if fit_mdv=='True' else False
fit_fluxes=True # if fit_fluxes=='True' else False
#print(condition,labeling_strategy,dataset,sample_no,tol,max_iters,n_restart,n_jobs,flux_bounds,fit_mdv,fit_fluxes,OUTPUT_FOLDER,FLUX_FILE,EX_FLUX_FILE)
OUTPUT_FILENAME=f'classic_fitflux{fit_fluxes}_{condition}'
##init log 
 
today = datetime.now().strftime("%Y%m%d")
os.makedirs(OUTPUT_FOLDER, exist_ok = True)

def signal_handler(signum, frame):
    print(f"[Signal Handler] Caught signal {signum}, setting flag")
    with open(FLAG_FILE, 'w') as f:
        f.write("TRUE")

FLAG_FILE=os.path.expanduser('~/HELPER_FREEFLUX_NF/exit_flag.txt')
if os.path.exists(FLAG_FILE):
    try:
        os.remove(FLAG_FILE)
    except FileNotFoundError:
        pass

def main():
    c13_problem=classic_c13.classic_c13_problem(
        condition=condition,labeling_strategy=labeling_strategy,dataset=dataset,sample_no=0,MODEL_FILE=MODEL_FILE,
        FLUX_FILE=FLUX_FILE,MS_FOLDER=INPUT_FOLDER,OUTPUT_FOLDER=OUTPUT_FOLDER,OUTPUT_FILENAME=OUTPUT_FILENAME
                                                ,tol=tol,max_iters=max_iters,
                                n_restart=n_restart,n_jobs=n_jobs,flux_bounds=flux_bounds,fit_fluxes=fit_fluxes,df_ex_flux=None)
    
    signal.signal(signal.SIGTERM, signal_handler) ## can catch signal to stop the process in case of timelimit#dont think its implemented yet here

    print(f'{c13_problem.fit_fluxes} fit fluxes')
    if c13_problem.labeling_strategy=='combined':
        c13_problem.labeling_strategy='1-13C'
        c13_problem._set_fit()
#        np.savetxt(f'{c13_problem.OUTPUT_FOLDER}/nullspace_matrix_{c13_problem.OUTPUT_FILENAME}.csv',c13_problem.fitA.model.null_space,delimiter=',')
        c13_problem.optModelA=c13_problem._build_optmodel()
        print('net_flux_range ASTPAm : see if bounds are properly set' , c13_problem.fitA.model.net_fluxes_range['ASPTAm']) ##test to check VA bounds or sample bounds
        #then second model
        c13_problem.modelB=c13_problem._set_model()
        c13_problem.labeling_strategy='U-13C'
        c13_problem._set_fit(alternative_model=True)
        c13_problem.optModelB=c13_problem._build_optmodel(alternative_model=True)
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
