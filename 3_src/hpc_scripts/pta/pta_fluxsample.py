'''sample fluxes using pta : need drG samples, orthants and cobramodel'''

import sys 
import os 
import argparse
import sys, json

with open(os.path.expanduser('~/202601_13C_files/PTA_SCRIPTS/sys_path.json')) as f: ## i managed to import pta when just doing python3.9 import pta on the shell so i dumped path in this file
    paths = json.load(f)

sys.path = paths
from pta.sampling import tfs
import pandas as pd 
import numpy as np
from datetime import datetime
sys.path.append(os.path.expanduser('~/'))

import TFpta
#sectors 
import ast
# from sampler import helper_sampling as hs

#folder
BASE_FOLDERNAME='202601_13C'

base_folder=f'{os.path.expanduser("~/")}/{BASE_FOLDERNAME}_files/INPUT_FILES/'
#base_folder='/home/nathan/thermo_flux/Nathan/glucose/202409_13C/'

INPUT_FOLDER=f'{base_folder}/1_data/c_synthetic_data/'
MODEL_FOLDER=f"{base_folder}/1_data/b_processed_files/"
MODEL_BOUNDS_FOLDER=f'{base_folder}/4_results/a_fluxspace/'
INITPOINTS_FOLDER=f'{base_folder}/4_results/b_ptasampling/'
FILENAME=f'{BASE_FOLDERNAME}_tf01'


cobrafile=f'{MODEL_FOLDER}/cobra_roundp3_ymodel.xml'

drGcovsqrtfile =f'{INPUT_FOLDER}/202601_13C_tf01_drG0_cov_sqrt.csv'
drG0file=f'{INPUT_FOLDER}/202601_13C_tf01_drG0mean.csv'
met_data_file=f'{MODEL_FOLDER}/yeast_met_data_formatted.csv'
rign_file=f'{INPUT_FOLDER}/tmodel_rxn_ignore_snd.csv'

    
#not all of this is useful but at least i make sure cmodel is built the same way than for drGsampling
def parse_args(args): #could put it in a module
    def int_from_exp(arg):
        """: The type argument of add_argument accept a callable that check and convert"""
        return int(float(arg))

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ("yes", "true", "t", "1"):
            return True
        elif v.lower() in ("no", "false", "f", "0"):
            return False
        else:
            raise argparse.ArgumentTypeError("Boolean value expected.")

    #common parser for shared args
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("condition", help="Condition name")
    # common_parser.add_argument("BNAME", help="Batch name")
    parser.add_argument("pta_folder", help="Folder for input")
    parser.add_argument("output_folder", help="Folder for saving output")
    parser.add_argument("--sector_n0", type=int, default=None, help="sector index")
    parser.add_argument("--num_samples", type=int_from_exp, default=1e4, help="Number of samples to generate")
    parser.add_argument("--num_inputdrG_samples", type=int_from_exp, default=1e4, help="Number of samples to generate")

    return parser.parse_args(args)

def main():
    args=parse_args(sys.argv[1:])
    #files dependent on args
    vbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_flux_bounds_{args.condition}.csv'
    lncbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_lnconc_bounds_{args.condition}.csv'
    # log_conc_mean_file= f'{INPUT_FOLDER}/EP_{args.condition}_lncmean_ad_beta1em12_relax.csv'
    # log_conc_var_file = f'{INPUT_FOLDER}/EP_{args.condition}_lncvar_ad_beta1em12_relax.csv'
    ##sector part 
    #get tfsmodel
    #prep vars
    #get tfsmodel
    prep=TFpta.preprocess.Preprocess(cobrafile,rign_file,vbounds_file,met_data_file,args.condition,lncbounds_file,drG0file,drGcovsqrtfile)

    ##sector part 
    # if args.sector_n0 is not None:
    #     sectors_ser=pd.read_csv(sector_file,index_col=0,converters={1: ast.literal_eval})
    #     sector_dict=sectors_ser.loc[args.sector_n0].values[0] 
    #     hs.set_sector_flux(prep.model,sector_dict)

    filename_IN=f'{args.pta_folder}/{args.condition}'
    print(f'>>>>>> reading {filename_IN} at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')
    drG_samples_full=pd.read_feather(f'{filename_IN}_samples.feather')
    drG_samples=drG_samples_full.loc[np.random.choice(drG_samples_full.index,args.num_inputdrG_samples,replace=False)]
    orthants_df=pd.read_csv(f'{filename_IN}_orthants.csv',index_col=0)
    print(f'>>>>>> input read at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')
    print(f'>>>>>> sampling fluxes at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')
    df=tfs.sample_fluxes_from_drg(prep.model,drG_samples,orthants_df,num_approx_samples=args.num_samples)

    time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print(f'>>>>>> done at {time} ')
    filename_OUT=f'{args.output_folder}/{args.condition}fluxsamples'
    df.to_csv(f'{filename_OUT}.csv')

main()