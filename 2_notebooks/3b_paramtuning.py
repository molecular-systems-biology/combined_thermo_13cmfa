###older script, could update with args parsing

import sys

import pandas as pd
import seaborn as sns
import os
# Add the parent directory to sys.path
sys.path.append(os.path.abspath('../3_src/'))
# Now import Model from freeflux
import numpy as np
import matplotlib.pyplot as plt
import itertools

sys.path.append('/home/nathan/PhD/04.Computational/202601_13C/3_src/')
sys.path.append('../3_src/')
import analysis.io_func as io_func

from analysis import config_pta as cfp

# Import general configs:
palette ,cfg =io_func.import_plotconfigs()
sns.set_palette([palette[2]]+[palette[5]]+[palette[3]])
#global
condition='TM6'
labeling_strategy='U_1-13C'

#folders
BASE_FOLDER=cfp.RES_13C_FOLDER
HPC_FOLDER=f'{BASE_FOLDER}/param_tuning/'
OUTPUT_FOLDER=cfp.PLOTS_FOLDER+'../SuppFigures/'
#start
n_restarts=[4,8,16]

##function to get in the dir , load files and get obj value 
def get_objs(n_restart):
    n_folder=os.path.join(HPC_FOLDER,f'n{n_restart}/')
    opt_files=[f for f in os.listdir(n_folder) if f.endswith('.csv') ]
    #store obj values as dict with sample number as key
    obj_values={}
    #obj value is the header of the csv file
    for f in opt_files:
        obj=float(str(pd.read_csv(f'{n_folder}/{f}').columns[0]).strip('# '))
        obj_values[int(f.split('_')[3].strip('sample').strip('.csv'))]=obj

    return obj_values

## we will plot 2 things :A)  mean obj value and std of obj value against n_restarsts  and B) obj value for each sample, with a different marker for each n_restart

def mean_std(obj_values):
    return np.mean(list(obj_values.values())),np.std(list(obj_values.values()))

def plot_mean_std(n_restarts,ax):
    means=[]
    stds=[]
    for n in n_restarts:
        obj_values=get_objs(n)
        mean,std=mean_std(obj_values)
        means.append(mean)
        stds.append(std)
    ax.errorbar(n_restarts,means,yerr=stds,fmt='o')
    # print(means,stds)

def plot_objs(n_restarts,ax,markers=['^','o','s','x']):
    for n in n_restarts:
        obj_values=get_objs(n)
        ax.scatter(list(obj_values.keys()),list(obj_values.values()),label=f'n_{n}',marker=markers[n_restarts.index(n)],edgecolor='black',s=62.5)

#plot objs of 1 run against the other, 2 by 2 
def scatter_obs(n1,n2,ax):
    obj_values1=get_objs(n1)
    obj_values2=get_objs(n2)
    #only keep the common samples
    common_samples=set(obj_values1.keys()).intersection(set(obj_values2.keys()))
    obj_values1={k:obj_values1[k] for k in common_samples} #reoder
    obj_values2={k:obj_values2[k] for k in common_samples}
    # print(obj_values1.keys(),obj_values2.keys())
    ax.scatter(list(obj_values1.values()),list(obj_values2.values()),label=f'n_{n1} vs n_{n2}',**cfg['settings_scatterplots'])
    ## get x and y limits and set them to be the same
    xmin,xmax=ax.get_xlim()
    ymin,ymax=ax.get_ylim()
    ax.set_xlim([min(xmin,ymin),max(xmax,ymax)])
    ax.set_ylim([min(xmin,ymin),max(xmax,ymax)])
    #add a diagonal line
    ax.axline((0, 0), slope=1, **cfg['settings_diagonal'])


fig3,ax3=plt.subplots(figsize=(2.5,2.5))
io_func.scale_rcparams(fig3)

##2 combination of 3  with itertools
combinations = list(itertools.combinations(n_restarts, 2))
for n1,n2 in combinations:
    scatter_obs(n1,n2,ax3)

ax3.legend(edgecolor='none',bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f'{OUTPUT_FOLDER}SuppFig_paramtuning_scatter_objs.png',dpi=300,bbox_inches='tight')
plt.savefig(f'{OUTPUT_FOLDER}SuppFig_paramtuning_scatter_objs.svg',dpi=300,bbox_inches='tight')


fig,ax=plt.subplots(figsize=(2.5,2.5))
plot_objs(n_restarts,ax)
ax.legend(edgecolor='none',bbox_to_anchor=(1.05, 1), loc='upper left')
ax.set_xlabel('Ordered samples')
ax.set_ylabel('$^{13}$C-fitting residuals')
plt.savefig(f'{OUTPUT_FOLDER}SuppFig_paramtuning_obj_vs_samples.png',dpi=300,bbox_inches='tight')
plt.savefig(f'{OUTPUT_FOLDER}SuppFig_paramtuning_obj_vs_samples.svg',dpi=300,bbox_inches='tight')

fig2,ax2=plt.subplots(figsize=(2.5,2.5))
plot_mean_std(n_restarts,ax2)
# plt.show()
plt.savefig(f'{OUTPUT_FOLDER}SuppFig_paramtuning_mean_std_vs_nrestarts.png',dpi=300,bbox_inches='tight')
plt.savefig(f'{OUTPUT_FOLDER}SuppFig_paramtuning_mean_std_vs_nrestarts.svg',dpi=300,bbox_inches='tight')

