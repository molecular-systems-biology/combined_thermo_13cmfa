import seaborn as sns
import numpy as np
import os 
import pandas as pd
import HELPER_FREEFLUX_NF.io as nf_io
from collections import defaultdict
import re
import csv

## make a function to draw N samples from the initial distributions 
def draw_subset_ini(net_flux_df,reaction,nsamples,accepted_samples2=None):
    '''Draw nsamples from the initial distribution of a reaction if no accepted samples are provided otherwise draw random among the accepted samples'''
    return np.random.choice(net_flux_df.loc[reaction,:],nsamples) if accepted_samples2 is None else np.random.choice(net_flux_df.loc[reaction,accepted_samples2],nsamples)


## func that takes condition, labeling strategy and returns series of sorted objective value# also need one func for net flux df 
def get_objective_values(condition, HPC_FOLDER):
    opt_files=[f for f in os.listdir(HPC_FOLDER) if f.endswith('.csv') and condition in f]
    #store obj values as dict with sample number as key
    obj_values={}
    #obj value is the header of the csv file
    for f in opt_files:
        firstline = str(pd.read_csv(f'{HPC_FOLDER}/{f}').columns[0]).strip('# ')
        obj_str = firstline.split('|')[0].strip()  # Get the part before '|'
        obj = float(obj_str)
        # obj=float(str(pd.read_csv(f'{HPC_FOLDER}/{f}').columns[0]).strip('# '))
        # obj_values[int(f.split('_')[3].strip('sample').strip('.csv'))]=obj # get sample n0 from name of file
        samplen0=re.findall(r'(?<=sample)\d+', f) # get sample n0 from name of file
        obj_values[int(samplen0[0])]=obj
    sorted_obj=pd.Series(obj_values).sort_values()
    return sorted_obj


def objvals_multiplesolvers(condition, HPC_FOLDER):
    opt_files = [f for f in os.listdir(HPC_FOLDER) if f.endswith('.csv') and f"_{condition}_" in f]
    obj_dict = defaultdict(dict)  # {sample_no: {method: obj_value}}

    for f in opt_files:
        # Extract method: between "<condition>_" and "_sample"
        method_match = re.search(fr'{condition}_(.*?)_sample', f)
        if not method_match:
            continue  # Skip if no method found
        method = method_match.group(1)

        # Extract sample number
        samplen0_match = re.search(r'(?<=sample)\d+', f)
        if not samplen0_match:
            continue
        samplen0 = int(samplen0_match.group())

        # Extract objective value from first column header
        try:
            header_line = str(pd.read_csv(os.path.join(HPC_FOLDER, f)).columns[0])
            obj_str = header_line.strip('# ').split('|')[0].strip()
            obj_value = float(obj_str)
        except Exception as e:
            print(f"Skipping {f}: {e}")
            continue

        obj_dict[samplen0][method] = obj_value

    # Return as DataFrame: index = sample number, columns = methods
    df = pd.DataFrame.from_dict(obj_dict, orient='index').dropna(how='all')
    return df.sort_index()

def get_netflux_df(model,flux_file):
    net_flux_df=nf_io.format_fluxes_tofreeflux2(pd.read_csv(flux_file,index_col=0,header=0),model)
    net_flux_df.columns=[int(c) for c in net_flux_df.columns]
    return net_flux_df

def vars_13c(condition,HPC_FOLDER,flux_file,model):
    '''condition: condition to be analyzed
    HPC_FOLDER: folder where the opt files are stored
    flux_file: file with the fluxes
    model: model used for the analysis
    OUPUT : sorted objective values, net flux df, accepted samples'''
    # get sorted objective values
    sorted_obj=get_objective_values(condition,HPC_FOLDER)
    # get net flux df
    net_flux_df=get_netflux_df(model,flux_file)
    #accepted samples are the ones with the top 20% of the objective values
    accepted_samples=sorted_obj[sorted_obj<=np.percentile(sorted_obj,20)].index
    return sorted_obj,net_flux_df,accepted_samples

## plot umap reduction in 2d, with colorbar (fig 5c). color points by their 13C score colormap
def calc_embedding(net_flux_df,accepted_samples=None):
    X=net_flux_df
    if accepted_samples is not None:
        X=net_flux_df.loc[:,accepted_samples]
    reducer = umap.UMAP(n_neighbors=15)
    embedding = reducer.fit_transform(X.T).T
    return embedding

def _load_u_fromfile(file_name,):
    with open(file_name, newline='') as csvfile:
        opt_u = list(csv.reader(csvfile))
    opt_u=np.array([float(i[0]) for i in opt_u[1:]])
    return opt_u
def load_u_fromsample(HPC_FOLDER,condition,sample_n0):
    ##get list of file names
    filename=f'{HPC_FOLDER}/optu_{condition}_combined_sample{sample_n0}.csv'
    return _load_u_fromfile(filename)

def hpd_interval(net_flux_df,accepted_samples,rxn,cred=0.95):
    samples=net_flux_df.loc[rxn, accepted_samples]
    sorted_samples=np.sort(samples)
    ci_idx_inc = int(np.floor(cred * len(sorted_samples)))
    n_cis = len(sorted_samples) - ci_idx_inc
    ci_width = sorted_samples[ci_idx_inc:] - sorted_samples[:n_cis]
    min_idx = np.argmin(ci_width)
    hpd_interval = (sorted_samples[min_idx], sorted_samples[min_idx + ci_idx_inc])
    return hpd_interval



def plot_diag(ax,cfg):
    ax.plot([0,1],[0,1],transform=ax.transAxes,**cfg['settings_diagonal'])


# def _one_KS_test_reaction(net_flux_df,reaction,accepted_samples1,accepted_samples2=None):
#     ## draw nsamples from the initial distribution
#     ## draw nsamples from the accepted samples distribution
#     ## perform KS test
#     # subset_ini=draw_subset_ini(net_flux_df,reaction,len(accepted_samples))
#     subset_accepted=net_flux_df.loc[reaction,accepted_samples1]
#     if accepted_samples2 is not None :
#         subset2_reduce=draw_subset_ini(net_flux_df,reaction,len(accepted_samples1),accepted_samples2)
#         return stats.ks_2samp(subset2_reduce,subset_accepted)
#     else : 
#         subset_ini=draw_subset_ini(net_flux_df,reaction,len(accepted_samples1))
#         return stats.ks_2samp(subset_ini,subset_accepted)

# def _KS_test_reaction(net_flux_df,reaction,accepted_samples,nrepeat,accepted_samples2=None):
#     pvalues=[]
#     for i in range(nrepeat):
#         if accepted_samples2 is None:
#             pvalues.append(_one_KS_test_reaction(net_flux_df, reaction, accepted_samples).pvalue)
#         else:
#             pvalues.append(_one_KS_test_reaction(net_flux_df, reaction, accepted_samples, accepted_samples2).pvalue)
#     return np.mean(pvalues)

# def KS_all_rxn(model,net_flux_df,accepted_samples,nrepeat,accepted_samples2=None):
#     '''model,net_flux_df,accepted_samples,nrepeat,accepted_samples2=None
#     in accepted samples2 provide the longest list of accepted samples'''
#     rxn_diff=[]
#     dict_pvals={}
#     for rxn in model.reactions : 
#         if accepted_samples2 is None :
#             pmean = _KS_test_reaction(net_flux_df,rxn,accepted_samples,nrepeat) 
#         else :
#             pmean = _KS_test_reaction(net_flux_df,rxn,accepted_samples,nrepeat,accepted_samples2)
#         if pmean <=0.05:
#             rxn_diff.append(rxn)
#             dict_pvals[rxn]=pmean

#     return rxn_diff,dict_pvals