'''config file to avoid repeating paths and function call at the beginning of every notebook zzz
for 13c notebooks, with 13C_INCA environment'''

import sys
import os
import pandas as pd
# Add the parent directory to sys.path
sys.path.append(os.path.abspath('../3_src/'))

# import numpy as np

BASE_FOLDERNAME='202601_13C'
BASE_FOLDER=f'/home/nathan/PhD/04.Computational/'
INPUT_FOLDER = f"{BASE_FOLDER}{BASE_FOLDERNAME}/1_data/b_processed_files/"
SYNDATA_FOLDER = f"{BASE_FOLDER}{BASE_FOLDERNAME}/1_data/c_synthetic_data/"

MODEL_BOUNDS_FOLDER=f'{BASE_FOLDER}{BASE_FOLDERNAME}/4_results/a_fluxspace/'
EP_FOLDER=f"{BASE_FOLDER}{BASE_FOLDERNAME}/4_results/b_sampling/b0_EP_lnc/"
RESULTS_FOLDER=f"{BASE_FOLDER}{BASE_FOLDERNAME}/4_results/b_sampling/"
FILENAME=f'{BASE_FOLDERNAME}_tf01'
PLOTS_FOLDER=f'{BASE_FOLDER}{BASE_FOLDERNAME}/2_notebooks/Figures/'
###for 13c and conc sampling 
RES_13C_FOLDER=f"{BASE_FOLDER}{BASE_FOLDERNAME}/4_results/c_freeflux_run/"
# OUTPUT_FOLDER_PLOT=f"{OUTPUT_FOLDER}/Figures/"
METSAMPLING_FOLDER=f"{BASE_FOLDER}{BASE_FOLDERNAME}/4_results/e_metconc_sampling/"
DRGSAMPLING_FOLDER=f"{BASE_FOLDER}{BASE_FOLDERNAME}/4_results/b_sampling/b4_drGsampling/"
## load some common variables

cobrafile=f'{INPUT_FOLDER}/cobra_roundp3_ymodel.xml'

drGcovsqrtfile =f'{SYNDATA_FOLDER}/{FILENAME}_drG0_cov_sqrt.csv'
# drGcovfile=f'{SYNDATA_FOLDER}/{FILENAME}_drG0_cov.csv'
rign_file=f'{SYNDATA_FOLDER}/tmodel_rxn_ignore_snd.csv'

drG0file=f'{SYNDATA_FOLDER}/{FILENAME}_drG0mean.csv'

# define defaults at module scope
vbounds_file = None
lncbounds_file = None
log_conc_mean_file=None
log_conc_var_file=None
def config_prep(cond):
    from TFpta import preprocess ##inside to still use config module in another env than pta env 

    global vbounds_file, lncbounds_file ,log_conc_mean_file,log_conc_var_file

    vbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_flux_bounds_{cond}.csv'
    lncbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_lnconc_bounds_{cond}.csv'
    # log_conc_mean_file= f'{EP_FOLDER}EP_{cond}_lncmean_ad_beta1em12_relax.csv'
    # log_conc_var_file = f'{EP_FOLDER}EP_{cond}_lncvar_ad_beta1em12_relax.csv'
    met_data_file=f'{INPUT_FOLDER}/yeast_met_data_formatted.csv'
    # drgbounds_file=f'{MODEL_BOUNDS_FOLDER}/{FILENAME}_nocellconc_drG_bounds_{cond}.csv' # for ep only
    prep=preprocess.Preprocess(cobrafile,rign_file,vbounds_file,met_data_file,cond,lncbounds_file,drG0file,drGcovsqrtfile)
    return prep

def drawn_flux_file(cond):
    df=pd.read_csv(f'{SYNDATA_FOLDER}/2_fluxsamples_for13c/{cond}_drawn_samples.csv',index_col=0)
    return df


def full_samples_df(cond):
    folder=f'{RESULTS_FOLDER}/b3_flux_sampling/{cond}/'
    df_concat=pd.concat([pd.read_csv(folder+file,index_col=0) for file in os.listdir(folder)])
    return df_concat
