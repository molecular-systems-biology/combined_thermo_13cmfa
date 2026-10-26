#-------------------------------------------------------------------------------
# Name:        fit_netfluxes.py

# Purpose:     Here we build a 13C model using freeflux and processed files created in 202409_1a and 1b 
               # then we fit  MDVs (mass distribution vectors) for all labelled amino acids (in target_EMU.json)
#               to each sampled NET flux distribution (/1_flux_sampling). 
#            we fix the net fluxes of the fitter - and the fitter will estimate the forward and backward fluxes of the reactions 
#            The exchange fluxes have to be set to irreversible in the model 
#               The score of each sample is the residual (least squares) compared to the experimental data (/processed_files).

# 
# 
# Author:      Nathan Fargier
#
# Created:     19/09/2024
### version 1.1 : modified solving function to choose different solvers. Scipy and openopt problems can be defined and solved. 
# ##              also added possibility to start from best solution
#-------------------------------------------------------------------------------

import sys
import os
# sys.path.append('/home/nathan/thermo_flux/')
# sys.path.append('/home4/p313427/freeflux/')
# sys.path.append('/home2/p313417/freeflux/')
sys.path.append(os.path.expanduser('~/freeflux/'))
import pandas as pd
from scipy.linalg import qr

from . import io as nf_io
# from scipy.optimize import approx_fprime
import numpy as np
try : 
    from freeflux.src.freeflux.solver.nlpsolver import MFAModel
except: 
    from freeflux.solver.nlpsolver import MFAModel

from . import force_flux as ff

try:
    import jax
    HAS_JAX = True
except ImportError:
    HAS_JAX = False
##### folders and files
###global variables are target_emu and model
njobs=8 #parallel jobs for EMU estimation and MDV simulation 

#hpc variables
# BASE_FOLDER='/home4/p313427/13C_INPUT_FILES/'
# BASE_FOLDER=os.path.expanduser('~//PhD/04.Computational/202601_13C')

# BASE_FOLDER=os.path.expanduser('~//PhD/04.Computational/202601_13C')
BASE_FOLDER=os.path.expanduser('~/202601_13C_files/INPUT_FILES/')


"""here we load the fitter from freeflux, the fluxes dataframe and the labeling data.
 THen we fix the net fluxes of the fitter and run the fit """

###first format the labeling data 
def data_to_MDV_mix_func(df_label_data):
    '''we return a df of arrays, and will use the fit.set_measured_MDV method, not from file
    for mixed dataset : with data of mhlab and from mareike'''
    ##0 std gives a inf residual ?
    # df_label_data['std']=np.where(df_label_data['std']==0,1e-5,df_label_data['std'])

    # Group by 'ms_id' and aggregate 'intensity' and 'std' into lists
    grouped_df = df_label_data.groupby('ms_id').agg(
        mean=('mean', lambda x: [float(i) for i in x]),
        sd=('total_error', lambda x: [float(i) for i in x])
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
                'MET':'metLc','PHE':'pheLc','TYR':'tyrLc','HIS':'hisLc','LYS':'lysLc'}
        
            # Loop through the replacement dictionary and apply replacements to the index
        for old, new in replcaement_id.items():
            grouped_df.index = grouped_df.index.str.replace(old, new)

        return grouped_df

    return  replace_aminoacids_name(grouped_df)
##loading the measured MDVs
def helper_set_MDV(fit,df_label_data):
    dict={}
    for fragment_ID, row in df_label_data.iterrows():
        dict[fragment_ID]=[row['mean'], row['sd']]
    fit.model.measured_MDVs =dict

def helper_set_MDV_bootstrap(fit, df_label_data, seed):
    dict = {}
    for fragment_ID, row in df_label_data.iterrows():
        rng = np.random.default_rng(seed=seed)
        sampled_mean = rng.normal(row['mean'], row['sd'])
        dict[fragment_ID] = [sampled_mean, row['sd']]
    fit.model.measured_MDVs = dict
    print(f"measured MDVs were set, {fit.model.measured_MDVs},with seed {seed}")

def helper_prepare(fit,n_jobs):
    """call this first before setting fluxes"""
    if not HAS_JAX:
        raise ImportError(
            "helper_prepare() requires JAX. Install with:\n"
            "python -m pip install 'jax[cpu]>=0.4,<=0.4.18'"
        )

    fit._decompose_network(n_jobs=n_jobs)
    fit._calculate_null_space()

    fit._calculate_transform_matrix()
    fit._lambdify_matrix_As_and_Bs()
    fit._calculate_matrix_As_and_Bs_derivatives_p('ss', n_jobs)
    fit._calculate_substrate_MDVs()
    fit._calculate_substrate_MDV_derivatives_p('ss' )
    fit._calculate_measured_MDVs_inversed_covariance_matrix()


def helper_set_flux_bounds(fit,net_flux_df,sample_no,flux_tol=1e-4,reset=False):
    """given a flux_df and a sample_no, set the fluxes bounds in the fitter
    if the flux is irreversible, the bounds are set to the absolute value of the flux.
      This means that transport fluxes are set to positive even if they are neg in the sample. 
      its ok because I made sure the direction of transport matched the net flux sample dir before. i.e co2t is in the direction of co2 excretion. """
    if reset : 
        fit._unset_flux_bounds(fit.model.reactions) 
    dict_fluxes=nf_io.load_sample_fromdf(sample_no,net_flux_df)
    for fluxid, value in dict_fluxes.items():
        #first fix irreversible reactions
        if fit.model.reactions_info[fluxid].reversible==0 : 
            fit.set_flux_bounds(fluxid,[ abs(value)-flux_tol, abs(value)+flux_tol]) ## i add 1e-4 because [0,0] does not work
        else:
            fit.set_flux_bounds(fluxid,[ value-flux_tol, value+flux_tol])  # for others , set the net flux bounds to a small range around the value

def fit_set_measfluxes(fit,df_ex_flux_cond):
    for irr_rid in df_ex_flux_cond.index:
        fit.set_measured_flux(irr_rid,df_ex_flux_cond.loc[irr_rid,'flux'],df_ex_flux_cond.loc[irr_rid,'flux_std_error'])
    print(f'measured fluxes were set,{fit.model.measured_fluxes}')
    
def helper_prep_part2(fit,bound=True):
    """call this after setting fluxes - as it set the net fluxes. if bound is False, it will not set the bounds on the net fluxes from the flux data
    """
    if fit.model.measured_fluxes != {}:
        fit._calculate_measured_fluxes_inversed_covariance_matrix()
        fit._calculate_measured_fluxes_derivative_p('ss')
    fit._set_default_flux_bounds()
    if bound == True or bound == 'VA':
        for rxn_id in fit.model.reactions: 
            fit.model.net_fluxes_range[rxn_id]=fit.model.net_fluxes_bounds[rxn_id] ## instead of   # fit._estimate_fluxes_range(fit.model.unbalanced_metabolites)
    else :
        fit._estimate_fluxes_range(fit.model.unbalanced_metabolites)
    print('bound were set if True otherwise default', bound)

## function to set bounds for classical 13C : using va bounds instead of sampled net flux bounds
va_bounds_file=f'{BASE_FOLDER}/4_results/a_fluxspace/data_bounds_fluxes_all_2ndlaw.csv'
def set_VA_bounds(fit,model,condition):
    va_bounds=pd.read_csv(va_bounds_file,index_col=[0,1]).loc[condition]
    c13_va_bounds=nf_io.format_fluxes_tofreeflux2(va_bounds,model) # cool func to turn flux id of cmodel into 13c model ## but messes up transporters 
    for fluxid in c13_va_bounds.index:
        if fit.model.reactions_info[fluxid].reversible==0 :  ##if irreversible, we set bounds to abs value because we already set the proper direction in model definition 
            ##all irreversible reactions carry a positive flux
            max_bound=max(abs(c13_va_bounds.loc[fluxid,'lb']),abs(c13_va_bounds.loc[fluxid,'ub']))
            min_bound=min(abs(c13_va_bounds.loc[fluxid,'lb']),abs(c13_va_bounds.loc[fluxid,'ub']))
            fit.set_flux_bounds(fluxid,[ min_bound,max_bound+1e-5])
        else :
            fit.set_flux_bounds(fluxid,[ c13_va_bounds.loc[fluxid,'lb'],c13_va_bounds.loc[fluxid,'ub']+1e-5])
    
def build_constraints_nf_ff(optModel,fit,flux_bounds=True,force_flux=False):
        '''in my version of building constraints, turning off flux bounds allows to only set constraint that vtot =Nu >=0
         we can also set linear constraints for force fluxes'''
        A1 = optModel.N
        A2 = optModel.T@optModel.N
        A3 = -A2
        
        b1 = np.zeros(optModel.ntotalfluxes)
        vnet_lb = np.array([optModel.model.net_fluxes_range[rxn_id][0] for rxn_id in fit.model.reactions])
        vnet_ub = np.array([optModel.model.net_fluxes_range[rxn_id][1] for rxn_id in fit.model.reactions])

        b2 = vnet_lb
        b3 = -vnet_ub
        
        if force_flux==True:
            A4,A5,b4,b5=ff.build_ff_matrix(optModel,fit)
            A=np.vstack((A1,A2,A3,A4,A5))
            b=np.concatenate((b1,b2,b3,b4,b5))
        else:
            A = np.vstack((A1, A2, A3))
            b = np.concatenate((b1, b2, b3))

        
        if flux_bounds:
            # if optModel.solver == 'slsqp': ##this is from initial freeflux code, not to be confused with me choosing different solvers
            #     optModel.constrs = {'type': 'ineq', 'fun': lambda u: A@u - b}
            
            # elif optModel.solver == 'ralg':
                optModel.A = -A
                optModel.b = -b
        else : 
            # if optModel.solver == 'slsqp':
            #     optModel.constrs = {'type': 'ineq', 'fun': lambda u: A1@u }
            
            # elif optModel.solver == 'ralg':
                optModel.A = -A1
                optModel.b = np.zeros(optModel.ntotalfluxes)

def build_constraints_nf_ff_eq(optModel,fit,force_flux=False):
        '''in my version of building constraints, turning off flux bounds allows to only set constraint that vtot =Nu >=0
         we can also set linear constraints for force fluxes'''
        A1 = optModel.N
        Aeq = optModel.T@optModel.N
        
        b1 = np.zeros(optModel.ntotalfluxes)

        beq = np.array(fit.vnet)
        
        if force_flux==True:
            A4,A5,b4,b5=ff.build_ff_matrix(optModel,fit)
            A=np.vstack((A1,A4,A5))
            b=np.concatenate((b1,b4,b5))
        else:
            A = A1
            b = b1

        optModel.A = -A
        optModel.b = -b

        ## reduce redundant rows in Aeq ## we reorder matrix A to place the independent ones on top
        ###such that the n=rank A first rows are independent and form a submatrix A_red that we are interested in
        # QR on transpose to finds independent rows
        Q, R, piv = qr(Aeq.T, pivoting=True) ## piv lists the swapping of indexes 
        rank = np.linalg.matrix_rank(Aeq)
        row_idx = piv[:rank] ### get the indices of the independent rows [:rank] first indexes

        Aeq_red = Aeq[row_idx, :]
        beq_red = beq[row_idx]
        #then set matrices
        optModel.Aeq= Aeq_red
        optModel.beq= beq_red

def build_objective_nf(optModel,fit, fit_mdv=True,fit_fluxes=False,clipping_value=None):
        '''in my version of the objective, also possibility to fit measured fluxes now, and possibility to fit vnet from sampled vnet'''
        def _f(u):
            optModel.model.total_fluxes[:] = optModel.N@u
            
            MDV_diff = optModel._calculate_difference_sim_exp_MDVs()
            MDV_inv_cov = optModel.model.measured_MDVs_inv_cov

            if type(clipping_value) == float: ## can try # clip max values to avoid numerical instabilities ## see if it improvees fit 
                max_value = clipping_value 
                MDV_inv_cov = np.clip(MDV_inv_cov, None, max_value)

            obj = MDV_diff@MDV_inv_cov@MDV_diff

            return obj
        
        def f1(u): ### fitting net fluxes
            obj1 = _f(u)
            
            flux_diff = optModel._calculate_difference_sim_exp_fluxes()
            flux_inv_cov = optModel.model.measured_fluxes_inv_cov
            obj2 = flux_diff@flux_inv_cov@flux_diff
            
            return obj1 + obj2
            
        def f2(u):
            TN=optModel.T@optModel.N
            vnet=fit.vnet
            ob2= np.linalg.norm(TN @ u - vnet)**2
            return ob2

        if fit_mdv : 
            optModel.f = f1 if fit_fluxes else _f
        else:
            optModel.f = f2

def build_optModel(fit,condition=None,force_flux=False,fit_fluxes=False,sample_n0=None,eq_constr=False,drg_df=None):
    optModel=None
    optModel = MFAModel(fit.model, fit_fluxes,solver = 'ralg') ##here ralg arguements sets freeflux to use func designer methods
    build_objective_nf(optModel,fit, fit_mdv=True,clipping_value=None,fit_fluxes=fit_fluxes) ##very high values (1e10) are clipped to mean (~ 1e6)
    optModel.build_gradient() #init optmodel.df
    optModel.build_hessian() #init optmodel.ddf
    if force_flux:
            ## retrieve dg_lb array ##todo add drGdf file as attribute in init
            assert drg_df is not None, 'provide drg_df'

            ###commented out if we use drGsample df 
            ## reindex 0 to N the first level (sample numbers)
            new_first_level=[i for i in range(len(drg_df.index.get_level_values(0).unique()))] 
            new_index=pd.MultiIndex.from_product(  [new_first_level] + [drg_df.index.get_level_values(1).unique()] )
            drg_df.index=new_index
            #take sample n0 , now from 0 to N
            drg_df_cond=drg_df.loc[sample_n0].copy()
            drg_df_cond=nf_io.format_fluxes_tofreeflux2(drg_df_cond,fit.model,fluxes=False) ## format the drG to freeflux format
            drg_df_cond=drg_df_cond.reindex(fit.model.reactions) #match rxn list (same than fit.vnet)
            optModel.dg_lb=drg_df_cond.iloc[:, 0].values 
            print('check dg_lb',optModel.dg_lb)
            optModel.dg_ub=drg_df_cond.iloc[:, 1].values

            
    if eq_constr is False :
        build_constraints_nf_ff(optModel,fit,flux_bounds=True,force_flux=force_flux) ## flux bounds is always true and if we want larger bounds (e.g for classic 13c, we set larger flux ranges above)
    else : 
        build_constraints_nf_ff_eq(optModel,fit,force_flux=force_flux)
    optModel._initialize_total_fluxes()
    return optModel
