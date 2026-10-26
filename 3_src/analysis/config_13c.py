'''config file to avoid repeating paths and function call at the beginning of every notebook zzz
for 13c notebooks, with 13C_INCA environment'''

import sys
import os
sys.path.append('/home/nathan/thermo_flux')
# Add the parent directory to sys.path
sys.path.append(os.path.abspath('../3_src/'))
from freeflux.src.freeflux import Model
import analysis.an_func as an_func
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd ##

BASE_FOLDER='/home/nathan/PhD/04.Computational/202601_13C/'
DATA_FOLDER=f'{BASE_FOLDER}/1_data/'
SYNDATA_FOLDER=f'{DATA_FOLDER}c_synthetic_data/'
DRAWNSAMPLES_FOLDER=f'{SYNDATA_FOLDER}/2_fluxsamples_for13c/'
INPUT_FOLDER=f'{DATA_FOLDER}b_processed_files/'

RESULTS_FOLDER=f'{BASE_FOLDER}/4_results/'
SAMPLINGDISTRIB_FOLDER=f'{RESULTS_FOLDER}/b_sampling/'
FREEFLUXRUN_FOLDER=f'{RESULTS_FOLDER}/c_freeflux_run/'
PLOTS_FOLDER=f'{FREEFLUXRUN_FOLDER}/Figures/'
GEN_PLOTS_FOLDER=f'{BASE_FOLDER}/2_notebooks/Figures/'
SUPP_PLOTS_FOLDER=f'{PLOTS_FOLDER}/SuppFigures/'
ISOFOLDER=f'{RESULTS_FOLDER}/d_predict_isotopes/'
METSAMPLING_FOLDER=f"{BASE_FOLDER}/4_results/e_metconc_sampling/"

# OUTPUT_FOLDER_PLOT=f"{OUTPUT_FOLDER}/Figures/"

## load some common variables

model= Model('scerevisiae') 
MODEL_FILE=f'{SYNDATA_FOLDER}/reactions_processed_ATP_2601.xlsx'#reactions_processed_ATP_testhalf
model.read_from_file(MODEL_FILE)

def net_flux_df(condition):
    net_flux_df=an_func.get_netflux_df(model,flux_file(condition))
    net_flux_df.columns=[int(c) for c in range(len(net_flux_df.columns))]
    return net_flux_df
def flux_file(condition):
    flux_file=f'{DRAWNSAMPLES_FOLDER}/{condition}_drawn_samples.csv'
    return flux_file

def full_flux_df(condition):
    flux_df=pd.read_csv(flux_file(condition),index_col=0)
    flux_df.columns=[int(c) for c in range(len(flux_df.columns))]
    return flux_df

def bootstrap_flux_df(condition, samples_as_cols=False):
    """Load weighted flux df, optionally using sample acutal numbers as columns."""

    if not samples_as_cols:
        flux_df = pd.read_csv(
            f'{FREEFLUXRUN_FOLDER}/bootstrap_flux_df_{condition}.csv',
            index_col=0
        )
    else:
        raw = pd.read_csv(
            f'{FREEFLUXRUN_FOLDER}/bootstrap_flux_df_{condition}.csv',
            header=None,
            dtype=str
        )

        flux_df = raw.iloc[1:].copy()
        flux_df.columns = raw.iloc[0]
        flux_df = flux_df.set_index(flux_df.columns[0])
        flux_df.columns.name = None
        flux_df.index.name = None
        flux_df.columns=[int(c) for c in flux_df.columns    ]


        flux_df = flux_df.apply(pd.to_numeric, errors='coerce')

    return flux_df

# def chosen_samples(condition):
#     chosen_samples=pd.read_csv(f'{FREEFLUXRUN_FOLDER}/bootstrap_chosensamples{condition}.csv',index_col=0).values.T[0]
#     return chosen_samples

def obj_vals(condition):
    '''note that we remove the 5% worse fitting samples : numerical errors ? we are only interested in the best fitting samples anyways'''
    gsubg_fold=f'{FREEFLUXRUN_FOLDER}{condition}_newn/'
    obj_vals=an_func.get_objective_values(condition,gsubg_fold)
    # obj_vals_filtered=obj_vals[obj_vals<np.percentile(obj_vals,95)]
    return obj_vals

# def cutoff(condition):
#     return 53 if condition == 'TM6' else 10

def find_plateau(obj_vals):
    '''fit curve to find a plateau value which we use as threshold for accepted samples'''
    def saturation(x, a, b, c):
        return a * (1 - np.exp(-b*x)) + c  # exponential saturation
    y=obj_vals.sort_values()
    x = np.arange(len(y))
    popt, _ = curve_fit(saturation, x, y, p0=[max(y), 1, min(y)])
    plateau_value = popt[0] + popt[2]  # asymptote
    return plateau_value