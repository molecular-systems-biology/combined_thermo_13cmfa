'''find initial points for PMO sampling'''

import sys 
import os 
import argparse
import sys, json
import ast #for sectors ser convert

with open(os.path.expanduser('~/202601_13C_files/PTA_SCRIPTS/sys_path.json')) as f: ## i managed to import pta when just doing python3.9 import pta on the shell so i dumped path in this file
    paths = json.load(f)

sys.path = paths

import pandas as pd 
import numpy as np
import pta
from datetime import datetime
sys.path.append(os.path.expanduser('~/'))

import TFpta


base_folder=f'{os.path.expanduser("~/")}/202601_13C_files/INPUT_FILES/'
#base_folder='/home/nathan/thermo_flux/Nathan/glucose/202409_13C/'

INPUT_FOLDER=f'{base_folder}/1_data/c_synthetic_data/'
MODEL_FOLDER=f"{base_folder}/1_data/b_processed_files/"
MODEL_BOUNDS_FOLDER=f'{base_folder}/4_results/a_fluxspace/'
INITPOINTS_FOLDER=f'{base_folder}/4_results/b_ptasampling/'


cobrafile=f'{MODEL_FOLDER}/cobra_roundp3_ymodel.xml'

drGcovsqrtfile =f'{INPUT_FOLDER}/202601_13C_tf01_drG0_cov_sqrt.csv'
drG0file=f'{INPUT_FOLDER}/202601_13C_tf01_drG0mean.csv'
rign_file=f'{INPUT_FOLDER}/tmodel_rxn_ignore_snd.csv'

# sector_file=f'{INPUT_FOLDER}/1024_sectors_0212.csv'

# n_init=178 if condition=='TM6' else 156
##logging ?
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_args(args):
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
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("condition", help="Condition name")
    common_parser.add_argument("BNAME", help="Batch name")
    common_parser.add_argument("output_folder", help="Folder for saving output")
    common_parser.add_argument("--sector_n0", type=int, default=None, help="sector index")
 
    ##subparsers for initial points and sampling
    parser=argparse.ArgumentParser('pta : generate initial points or sample with an init points array')
    subparsers = parser.add_subparsers(dest="run_type", required=True, help="Run mode, choose between initial_points and sample")
    #init point mode
    p_init = subparsers.add_parser("initial_points", help="Generate initial points",parents=[common_parser])
    p_init.add_argument("n_direction", type=int, help="Direction index")

    #sample 
    p_sample = subparsers.add_parser("sample", help="Run sampling from initial points",parents=[common_parser])
    p_sample.add_argument("--num_samples", type=int_from_exp, default=5e5, help="number of samples")
    p_sample.add_argument("--num_dir_samples", type=int_from_exp, default=5e7, help="number of direction samples")
    p_sample.add_argument("--max_steps", type=int_from_exp, default=2**25, help="Maximum steps for sampling : determines thinning")
    p_sample.add_argument("--num_initial_steps", type=int_from_exp, default=1e7, help="Number of initial steps for sampling")
    p_sample.add_argument("--max_threads", type=int, default=32, help="Number of threads for sampling")
    p_sample.add_argument("--warmup", type = int, default=2,help="npoints to discard")
    ##add arguments to run single chains on single jobs
    p_sample.add_argument('--multijobs', action='store_true', help='Run sampling in parallel : each job will be one chain')
    p_sample.add_argument("--chain_index", type=int, default=None, help="Chain index (used only when --multijobs is enabled)")
    p_sample.add_argument("--findpoints", action='store_true', help="Find initial points in the sampling job")


    return parser.parse_args(args)

def main():
    time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print(f'>>>>>> started at {time} ')
    ##parsse args bc we need condiiton
    args=parse_args(sys.argv[1:])
    FILENAME=f'202601_13C_{args.BNAME}'

    vbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_flux_bounds_{args.condition}.csv'
    lncbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_lnconc_bounds_{args.condition}.csv'
    ##sector part 
    # sectors_ser=pd.read_csv(sector_file,index_col=0,converters={1: ast.literal_eval})
    # sector_dict=sectors_ser.loc[args.sector_n0].values[0] 
    # ##ln C distribution are sector dependent 
    # log_conc_mean_file= f'{INPUT_FOLDER}/logconcdist/EP_lncmean_ad_beta1em12_relax_sector{args.sector_n0}.csv'
    # log_conc_var_file = f'{INPUT_FOLDER}/logconcdist/EP_lncvar_ad_beta1em12_relax_sector{args.sector_n0}.csv'
    #or not 
    # log_conc_mean_file= f'{INPUT_FOLDER}/EP_{args.condition}_lncmean_ad_beta1em12_relax.csv'
    # log_conc_var_file = f'{INPUT_FOLDER}/EP_{args.condition}_lncvar_ad_beta1em12_relax.csv'
    met_data_file=f'{MODEL_FOLDER}/yeast_met_data_formatted.csv'
    
    #get tfsmodel
    tfs_model = TFpta.main(cobrafile,rign_file,vbounds_file,met_data_file,args.condition,lncbounds_file,drG0file,drGcovsqrtfile,
                            ranktol=1e-5, sector_dict=None)
    run_type=args.run_type

    if run_type.lower()=='initial_points':
        initial_points=tfs_model.safe_find_point(tfs_model.reversible_dirs[args.n_direction]) #(n,1) array
        time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        print(f'{run_type}>>>>>> done at {time} ')
        pd.Series(initial_points.flatten()).to_csv(f'{args.output_folder}/{args.condition}_initialpoint_ndir{args.n_direction}.csv') 

    elif run_type.lower()=='sample':
        ##read init points : if not multijobs, full array, otherwise one column corresponding to chain index
        if args.multijobs: #single chain job so multiple jobs
            if args.chain_index is None:
                raise ValueError("When using --multijobs, you must specify --chain_index.")
            else :
                points_full_arr=pd.read_csv(f'{INITPOINTS_FOLDER}/initialpoints_{args.condition}.csv',index_col=0).values
                points_arr=points_full_arr[:,args.chain_index].reshape(-1,1) #for a (n,1) array
                filename=f'{args.output_folder}/{args.condition}_chain{args.chain_index}'

        else: ## multiple chains job
            if args.chain_index is not None:
                print("--chain_index was provided but --multijobs is False; ignoring it'")

            if args.findpoints: ##in the same job, compute init points
                print(f'>>>>>> finding initial points at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')
                points_arr = tfs_model.get_initial_points(args.max_threads)
                time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")
                print(f'>>>>>> done at {time} ')

            else : 
                print(f'>>>>>> reading initial points')
                points_full_arr=pd.read_csv(f'{INITPOINTS_FOLDER}/initialpoints_{args.condition}_sectorfwd{args.sector_fwd}.csv',index_col=0).values
                ## randomly choose n=maxthreads numbers of columns in pointsarray ## if there are more columns than threads
                points_arr=points_full_arr[:,np.random.choice(points_full_arr.shape[1],args.max_threads,replace=False)]

            #finally, build filename
            filename=f'{args.output_folder}/{args.condition}'

        ##put params for sampling
        tfs_model.num_samples=args.num_samples
        tfs_model.num_dir_samples=args.num_dir_samples
        tfs_model.max_steps=args.max_steps
        tfs_model.num_initial_steps=args.num_initial_steps
        tfs_model.max_threads=args.max_threads
        tfs_model.warmup=args.warmup
        result = tfs_model.sample(points_arr)
        time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        print(f'{run_type} >>>>>> done at {time} ')
        result.orthants.to_csv(f'{filename}_orthants.csv')
        result.samples.to_feather(f'{filename}_samples.feather',compression='zstd',compression_level=9)  ##write feather and compress to save space


main()


