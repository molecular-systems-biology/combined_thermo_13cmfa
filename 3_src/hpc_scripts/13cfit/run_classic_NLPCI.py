#-------------------------------------------------------------------------------
# Name:        run_classic_NLPCI.py

# Purpose:     build classic c13problem class and run NLP CI search script on HPC (folders defined there)
# 
# 
# Author:      Nathan Fargier
#
# Created:     13/03/2025

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

sys.path.append(os.path.expanduser('~/freeflux/'))
sys.path.append(os.path.expanduser('~/'))
import pandas as pd
import os
import numpy as np
from HELPER_FREEFLUX_NF import classic_c13problem as classic_c13
from datetime import datetime
import signal
### deactivate warnings 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
###
BASE_FOLDER=os.path.expanduser('~/202601_13C_files/INPUT_FILES/1_data/')
# BASE_FOLDER=os.path.expanduser('~/thermo_flux/Nathan/glucose/202409_13C/1_data/')

INPUT_FOLDER=f'{BASE_FOLDER}b_processed_files/'
MODEL_FILE=f'{BASE_FOLDER}c_synthetic_data/reactions_processed_ATP.xlsx'

arguments=sys.argv[1:]
###here args should be  self.optu = optu ## has to be set after the first fit (call of classic_c13problem)
        # self.confidence_level = confidence_level
        # self.max_iterations = max_iterations
        #flux_indx
        ##direction
if len(arguments)!=10:
    raise ValueError('10 arguments are needed') 
condition,labeling_strategy,opt_u,flux_idx,direction,tol,max_iterations,OUTPUT_FOLDER,OUTPUT_FILENAME,FLUX_FILE=arguments
n_jobs=8
flux_idx=int(flux_idx)
direction=int(direction)
tol=float(tol)
max_iterations=int(max_iterations)
##opt_u is an text file that needs to be converted to a numpy array
opt_u=pd.read_csv(opt_u,header=0).values.T[0]
print(opt_u.shape)


print(condition,labeling_strategy,'mixed',flux_idx,direction,OUTPUT_FOLDER,FLUX_FILE)

# =f'classic_{condition}_flux{flux_idx}_dir{direction}'
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
    c13_problem=classic_c13.classic_c13_problem(optu=opt_u, flux_idx=flux_idx, direction=direction,
        condition=condition,labeling_strategy=labeling_strategy,dataset='mixed',sample_no=0,MODEL_FILE=MODEL_FILE,
        FLUX_FILE=FLUX_FILE,MS_FOLDER=INPUT_FOLDER,OUTPUT_FOLDER=OUTPUT_FOLDER,OUTPUT_FILENAME=OUTPUT_FILENAME
                                                ,tol=tol,max_iters=max_iterations,n_jobs=1)
    
    signal.signal(signal.SIGTERM, signal_handler) ## can catch signal to stop the process in case of timelimit ##dont think its implemented yet here

    c13_problem.labeling_strategy='1-13C'
    c13_problem._set_fit()
    c13_problem.optModelA=c13_problem._build_optmodel()
    print('net_flux_range ASTPAm : see if bounds are properly set' , c13_problem.fitA.model.net_fluxes_range['ASPTAm']) ##test to check VA bounds or sample bounds
    #then second model
    c13_problem.modelB=c13_problem._set_model()
    c13_problem.labeling_strategy='U-13C'
    c13_problem._set_fit(alternative_model=True)
    c13_problem.optModelB=c13_problem._build_optmodel(alternative_model=True)
    
    ## run search 
    obj_val,vi,is_success,all_cstr=c13_problem.nlp_opt_singleflux()


    print(f'>>>>>> Optimization successful : {is_success} with objective {obj_val} and max cstr {all_cstr}')

    ## save the array in a csv file
    #os.makedirs(OUTPUT_FOLDER, exist_ok = True)
    np.savetxt(f'{c13_problem.OUTPUT_FOLDER}/optu_{c13_problem.OUTPUT_FILENAME}.csv',vi,delimiter=',',header=str(obj_val))
    
main()
