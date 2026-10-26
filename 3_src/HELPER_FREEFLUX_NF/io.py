#-------------------------------------------------------------------------------
# Name:        io.py
# Purpose:     module to help load and format files needed to run freeflux.
# 
# 
# Author:      Nathan Fargier
#
# Created:     19/09/2024
#-------------------------------------------------------------------------------


import sys
sys.path.append('/home/nathan/thermo_flux/')

import pandas as pd
import ast
import os
# from freeflux.src.freeflux import Model
import json
import numpy as np

##### folders and files
BASE_FOLDER=os.path.expanduser('~/202601_13C_files/INPUT_FILES/')#'/home/nathan/thermo_flux/Nathan/glucose/202409_13C/'
FLUXES_FOLDER=f'{BASE_FOLDER}synthetic_data/1_flux_sampling/'

MODEL_FILE=f'{BASE_FOLDER}synthetic_data/reactions_processed.xlsx'
emu_FILE=f'{BASE_FOLDER}synthetic_data/target_EMU.json'

MS_FOLDER=f'{BASE_FOLDER}processed_files/'


def load_model(model_file):
    model= Model('scerevisiae')
    model.read_from_file(model_file)
    return model


#### Functions to load fluxes : 
## one function to load the complete dataframe, index = reactions and columns = samples No
## one function to format the fluxes for the model, we will remove the H and +H at the end of the rxn_id
## one function to transform fluxes into forward and backward fluxes for reversible reactions
## final function that takes a sample no and returns the dict, ultimate function before setting fluxes and running simulator
def load_fluxes(fluxes_file):
    flux_df=pd.read_csv(fluxes_file,index_col=0)
    return flux_df


### function to read the fluxes file and return a dict with rxn id and flux
## make up a function that fomrats the flux df with the coorect ids then return a dict with rxn id and flux
def format_fluxes_tofreeflux2(fluxes_df,model,fluxes=True):
    fluxes_df.reset_index(inplace=True)
    fluxes_df['index']=[f'A{i}' if i[0].isdigit() else i for i in fluxes_df['index'].values]
    fluxes_df['index']=fluxes_df['index'].str.replace('cell-wall_reaction','cellwall_reaction')
    ### replace all 0,1,2,3 with nothing if they are after a hyphen 
    for i in range(5):
        fluxes_df['index']=fluxes_df['index'].replace(f'-{i}','',regex=True)
        fluxes_df['index']=fluxes_df['index'].replace(f'\+{i}','',regex=True)
    ### Remove '-H', '+H', and 'H' only if preceded by a lowercase letter at the end of the string
    fluxes_df['index'] = fluxes_df['index'].replace(r'(\-H$|\+H$|(?<=[a-z])H$)', '', regex=True)

        ### replace ATPASE and IONPUMP by ATPHYD ## careful not to replace ionpump first then ionpumpm ..
    fluxes_df['index']=fluxes_df['index'].replace('ATPASE','ATPHYD',regex=True).replace('IONPUMPm','ATPHYD',regex=True).replace('IONPUMP','ATPHYD',regex=True)
    ##ADK1m and ADK1 are the same reaction in the 13cmodel
    fluxes_df['index']=fluxes_df['index'].replace('ADK1m','ADK1',regex=True)
    ### reaction BPNT and SO4SO3 are aggregated in BPNTSO4 in the 13c model, count the flux only once 
    fluxes_df['index']=fluxes_df['index'].replace('BPNT','BPNTSO4',regex=True)
    ### now some rxn_id are the same, concatenate their flux value
    if fluxes is True :
        new_df=fluxes_df.groupby('index').agg({'sum'}).reset_index()
        new_df.set_index('index',inplace=True)
        new_df.columns=new_df.columns.droplevel(1)

    else: #drg case : we use it on the drG va df and dont want to sum up values but take min-lb and max-ub
        new_df=fluxes_df.groupby('index').agg(lb=('lb', 'min'), ub=('ub', 'max')).reset_index()#.agg({'sum'}).reset_index()
        new_df.set_index('index',inplace=True)
    ##only keep reactions that are in the model
    new_df=new_df.loc[[i for i in new_df.index if i in model.reactions]]
    ###there s a number of ex reaction that dont go in the same direciton than c to e transport reaction, yet we use the latter as ex reaction (c->e).
        ## so take the opposite of the flux from initial df of fluxes ## we take abs anyways so no need
    # new_df.loc[['ACt','ETOHt','CO2t','Mant','PYRt','SUCCt']]=-new_df.loc[['ACt','ETOHt','CO2t','Mant','PYRt','SUCCt']]

    return new_df

##function to further format the fluxes further format fluxes : each reversible reaction must have 2 fluxes f and b , 
###if positive flux f has all value , if negative flux b has all value 
def f_b_fluxes_df(fluxes_df: pd.DataFrame, model) -> pd.DataFrame:
    '''using vectorized operations to split the fluxes into forward and backward fluxes : np.clip limits the values in an array. '''
    # Initialize an empty list to store the new rows for performance improvement
    new_rows = []

    def f_b_rows(row):
        if model.reactions_info[row.name].reversible:
            # Vectorized operations: split the row into forward and backward fluxes
            f_row = row.clip(lower=0)  # Forward flux (positive values)
            b_row = abs(row.clip(upper=0))  # Backward flux (negative values)
            f_row.name = row.name + '_f'
            b_row.name = row.name + '_b'
            new_rows.append(f_row)
            new_rows.append(b_row)
        else:
            unchanged_row=row.copy()
            new_rows.append(unchanged_row) # If not reversible, keep the row unchanged

    # Apply the optimized function to each row of the DataFrame
    fluxes_df.apply(f_b_rows, axis=1)
    
    # Concatenate all rows together into a new DataFrame
    new_df = pd.concat(new_rows, axis=1).T
    
    return new_df

def load_sample_fromdf(sample_no,df):
    return df.iloc[:,sample_no].to_dict()


###todo : change input format of exp data for a better one (not for inca)
def exp_data_format(condition,labeling_strategy):
    '''function to read the csv in processed files (built in 1a)
    returns a df with MDV (intensity and std) for each metabolite isotope (EMU)'''
    assert condition in ['TM6','WT'], 'condition must be either TM6 or WT'
    #put condition in lowercase
    condition=condition.lower()
    #load file
    exp_file=f'{MS_FOLDER}/202409_13C_MSdata_{condition}_{labeling_strategy}.csv'
    exp_data=pd.read_csv(exp_file,index_col=0,converters={'atom_mdv':ast.literal_eval, 'atom_ids':ast.literal_eval})
    df_label_data=exp_data.groupby(['ms_id', 'mass_isotope']).agg({'intensity':'mean'})
    ##also keep intensity_std_error
    df_label_data['std']=exp_data.groupby(['ms_id', 'mass_isotope']).agg({'intensity_std_error':'mean'})#.head(20)
    return df_label_data


#for each ms id, retrieve values in the intensity and return an array (also for std)
def data_to_MDV_func(df_label_data):
    '''we return a df of arrays, and will use the fit.set_measured_MDV method, not from file'''
    ##0 std gives a inf residual ?
    df_label_data['std']=np.where(df_label_data['std']==0,1e-5,df_label_data['std'])

    # Group by 'ms_id' and aggregate 'intensity' and 'std' into lists
    grouped_df = df_label_data.groupby('ms_id').agg(
        mean=('intensity', lambda x: [float(i) for i in x]),
        sd=('std', lambda x: [float(i) for i in x])
    )
    grouped_df['mean'] = grouped_df['mean'].apply(np.array)
    grouped_df['sd'] = grouped_df['sd'].apply(np.array)
    
    # Optional: Rename the 'ms_id' to match your desired output 'fragment_ID'
    grouped_df.rename(columns={'ms_id': 'fragment_ID'}, inplace=True)
    ## last step is to remove [' and '] in index 
    grouped_df.index=[i.replace("'",'').replace('[','').replace("]",'').replace("'",'').replace(',','') for i in grouped_df.index]

    #def replace aminoacids names
    def replace_aminoacids_name(grouped_df):
        replcaement_id={'ALA':'alaLc','ASP':'aspLc','GLU':'gluLc','GLY':'glyc','SER':'serLc',
                'THR':'thrLc','VAL':'valLc','LEU':'leuLc','ILE':'ileLc','PRO':'proLc',
                'MET':'metLc','PHE':'pheLc','TYR':'tyrLc','HIS3':'hisLc','LYS':'lysLc'}
        
            # Loop through the replacement dictionary and apply replacements to the index
        for old, new in replcaement_id.items():
            grouped_df.index = grouped_df.index.str.replace(old, new)

        return grouped_df

    return  replace_aminoacids_name(grouped_df)