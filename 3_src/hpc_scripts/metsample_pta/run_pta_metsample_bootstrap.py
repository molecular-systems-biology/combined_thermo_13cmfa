'''metabolite sampling on a flux sample using PTA 
first drG sample with fluxes fixed to the flux sample, then metabolite sample from the drG samples obtained'''

import json 
import os
import sys
with open(os.path.expanduser('~/202601_13C_files/PTA_SCRIPTS/sys_path.json')) as f: ## i managed to import pta when just doing python3.9 import pta on the shell so i dumped path in this file
    paths = json.load(f)

sys.path.append(os.path.expanduser('~/'))

import TFpta
import pta 
from pta.sampling import tfs


import argparse
from datetime import datetime
import pandas as pd 
import numpy as np
import cvxpy as cp

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
##argparse
def parse_args(args): 
    def int_from_exp(arg):
        """: The type argument of add_argument accept a callable that check and convert"""
        return int(float(arg))
        
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("condition", help="Condition name")
    # common_parser.add_argument("BNAME", help="Batch name")
    # parser.add_argument("chosensamp_file", help="csv of objective values")
    parser.add_argument('flux_df_file', help="csv of sampled fluxes")
    parser.add_argument("output_folder", help="Folder for saving output")
    parser.add_argument("--sample_n0", type=int, default=None, help="sector index")
    parser.add_argument('--n_drgsamples',default=6,type=int,help='integer setting for drG sampling')
    parser.add_argument("--n_initpoint", type=int, default=12, help="Number of initpoint for drg sampling")
    parser.add_argument("--max_threads", type=int, default=4)


    return parser.parse_args(args)

args=parse_args(sys.argv[1:])

##folders
BASE_FOLDERNAME='202601_13C'

base_folder=f'{os.path.expanduser("~/")}/{BASE_FOLDERNAME}_files/INPUT_FILES/'
INPUT_FOLDER=f'{base_folder}/1_data/c_synthetic_data/'
MODEL_FOLDER=f"{base_folder}/1_data/b_processed_files/"
MODEL_BOUNDS_FOLDER=f'{base_folder}/4_results/a_fluxspace/'
met_data_file=f'{MODEL_FOLDER}/yeast_met_data_formatted.csv'

FILENAME=f'{BASE_FOLDERNAME}_tf01'
cobrafile=f'{MODEL_FOLDER}/cobra_roundp3_ymodel.xml'
drGcovsqrtfile =f'{INPUT_FOLDER}/202601_13C_tf01_drG0_cov_sqrt.csv'
drG0file=f'{INPUT_FOLDER}/202601_13C_tf01_drG0mean.csv'
rign_file=f'{INPUT_FOLDER}/tmodel_rxn_ignore_snd.csv'
vbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_flux_bounds_{args.condition}.csv'
lncbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_lnconc_bounds_{args.condition}.csv'

def main():
    prep=TFpta.preprocess.Preprocess(cobrafile,rign_file,vbounds_file,met_data_file,args.condition,lncbounds_file,drG0file,drGcovsqrtfile)

    ##flux sample  

    df_samples=pd.read_csv(args.flux_df_file,index_col=0)
    df_samples.columns=[int(i) for i in range(len(df_samples.columns))] ## make the cols be  1-1000 instead of actual drawn samples number 

    v_sample=df_samples.iloc[:,args.sample_n0] 

    ##set model bounds as vsample 
    for rid,val in v_sample.items():
        prep.model.reactions.get_by_id(rid).bounds=(val,val)

    thermodynamic_space = TFpta.ThermodynamicSpaceMod(
        prep.Sconstrained,
        prep.rid_constrained,
        prep.model.metabolites,
        prep.drg0_prime_mean,
        prep.drg0_prime_cov_sqrt,
        prep.log_conc_cov,
        prep.log_conc_mean,
    )

    tbasis = TFpta.ThermodynamicSpaceBasismod(
        thermodynamic_space,
        explicit_log_conc=False,
        explicit_drg0=False,
        explicit_drg=True,
        min_eigenvalue=1e-10,
        ranktol=1e-5,
        ignoreconc_met_idx=prep.ignoreconc_met_idx,
    )

    problem = pta.PmoProblem(prep.model, thermodynamic_space, tbasis, solver="GUROBI")

    tfs_model = TFpta.SupTFSmodel(prep.model, thermodynamic_space, tbasis, problem, solver=cp.GUROBI)

    points_arr=tfs_model.get_initial_points(args.n_initpoint) ###

    n=args.n_drgsamples##
    tfs_model.num_samples=int(10**(n-1))
    tfs_model.num_dir_samples=int(10**(n+1))
    tfs_model.max_steps=int(10**(n+2))
    tfs_model.num_initial_steps=int(10**(n+1))
    tfs_model.max_threads=args.max_threads
    tfs_model.warmup=-1

    print(f'sampling {tfs_model.num_samples} points with {tfs_model.num_dir_samples} direction samples at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')
    result = tfs_model.sample(points_arr)


    min_eigenval=1e-5
    basis_conc= TFpta.ThermodynamicSpaceBasismod(
        thermodynamic_space,
        explicit_log_conc=True,
        explicit_drg0=False,
        explicit_drg=True,
        min_eigenvalue=min_eigenval,
        ranktol=1e-5,
        ignoreconc_met_idx=prep.ignoreconc_met_idx,
    )
 
    print(f'sampling lnconc samples at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')
    log_conc_samples = tfs._sample_conditional_mvn(
            basis_conc.to_observables_transform[1],
            basis_conc.to_observables_transform[0],
            result.samples.to_numpy().T,
            basis_conc.observables_ranges["log_conc"],
            min_eigenval,
        )

    log_conc=pd.DataFrame(log_conc_samples.T, columns=[m.id for idx,m in enumerate(prep.model.metabolites) if idx not in prep.ignoreconc_met_idx])


    log_conc.to_hdf(f"{args.output_folder}/metconc_samples_{args.condition}_sample{args.sample_n0}.h5", key="data", mode="w", complevel=9, complib='blosc2:zstd') ##write h5 and compress to save space
    print(f'>>>>>> writing done at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')


if __name__ == "__main__":
    main()