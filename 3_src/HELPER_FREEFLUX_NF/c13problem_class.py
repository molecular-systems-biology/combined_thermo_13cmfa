#-------------------------------------------------------------------------------
# Name:        c13problem_class.py

# Purpose:     Class to build (freeflux) and solve the 13C fitting problem with NLP solver of open opt.
#  The class allows for combined optimization of 2 models (with U-13C and 1-13C labeling strategies)
#  Uses functions from HELPER_FREEFLUX_NF.fit_netfluxes
# 
# 
# Author:      Nathan Fargier
#
# Created:     10/10/2024

##version 1.1 with multiple start solving
##version 1.2 with class and combined optimization implemented 
##version 1.3 moved the class to a separate file from the execution py script (run_fitting.py)s
#-------------------------------------------------------------------------------

import sys
import os 
# sys.path.append('/home4/p313427/freeflux/')
# sys.path.append('/home4/p313427/')

# sys.path.append('/home2/p313417/freeflux/')
# sys.path.append('/home2/p313417/')
sys.path.append(os.path.expanduser('~/'))
sys.path.append(os.path.expanduser('~/freeflux/'))

#sys.path.append(os.path.abspath('../'))

import pandas as pd
import ast
import os
try : 
    from freeflux import Model
except:
    from freeflux.src.freeflux import Model
    
from HELPER_FREEFLUX_NF import io as nf_io
import numpy as np
from HELPER_FREEFLUX_NF import fit_netfluxes as nf_fit
from HELPER_FREEFLUX_NF import solvers
from datetime import datetime
### deactivate warnings 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#with open(log_file_path, 'a'):
#    pass 
class c13_problem:
    '''class to build (freeflux) and solve the 13C fitting problem with NLP solver of open opt. 
     Arguments : condition, labeling_strategy, sample_no, MODEL_FILE, FLUX_FILE, MS_FOLDER, OUTPUT_FOLDER, OUTPUT_FILENAME, tol=1e-6, max_iters=200,
                    n_restart=3,n_jobs=8,flux_bounds=True,fit_mdv=True
     '''
    def __init__(self,condition:str,labeling_strategy:str,sample_no:int,MODEL_FILE:str,FLUX_FILE:str,MS_FOLDER:str,OUTPUT_FOLDER:str,
                OUTPUT_FILENAME:str,tol=1e-6,max_iters=200,scale=1e-2, dataset='mixed', solver='gsubg',FLAG_FILE=None,
                 n_restart=3,n_jobs=8,flux_bounds=True,fit_mdv=True,force_flux=True,fit_fluxes=False,eq_constr=True,
                 df_ex_flux=False,extended_data=False,seed=None):
        self.condition=condition
        self.labeling_strategy=labeling_strategy #U-13C or 1-13C or combined
        self.dataset=dataset ##mareike or MHlab or mixed
        self.sample_no=sample_no
        self.tol=tol
        self.max_iters=max_iters
        self.n_restart=n_restart
        self.n_jobs=n_jobs
        self.flux_bounds=flux_bounds
        self.fit_mdv=fit_mdv
        self.OUTPUT_FOLDER=OUTPUT_FOLDER
        self.MODEL_FILE=MODEL_FILE
        self.FLUX_FILE=FLUX_FILE
        self.MS_FOLDER=MS_FOLDER
        self.OUTPUT_FILENAME=OUTPUT_FILENAME
        self.optModelA=None
        self.optModelB=None
        self.modelA=self._set_model()
        self.force_flux=force_flux
        self.drg_df=None ## df with drG bounds for force_Flux
        self.fit_fluxes=fit_fluxes
        self.solver=solver
        self.scale=scale
        self.eq_constr=eq_constr
        ##seed for bootstrapping of 13C data
        self.seed=seed


        ##flag file 
        self.FLAG_FILE=FLAG_FILE
        # self.df_ex_flux=(pd.read_csv(os.path.expanduser('~/13C_INPUT_FILES/1_data/b_processed_files/202409_13C_ex_flux_data.csv'),index_col=[0,3]) 
        #     if df_ex_flux is None else df_ex_flux
        # )
        self.df_ex_flux=(pd.read_csv(f'{MS_FOLDER}/202409_13C_ex_flux_data.csv',index_col=[0,3]) 
            if df_ex_flux is None else df_ex_flux
        )
        self.extended_data=extended_data
    

    def __call__(self):
        '''create an instance of the solver class and solve it '''
        c13_solver=solvers.c13_Solver(optModelA=self.optModelA,
                                      vnet=self.fitA.vnet,
                                      tol=self.tol,
                                      max_iters=self.max_iters,
                                      n_restart=self.n_restart,
                                      solver=self.solver,
                                      scale=self.scale,
                                      OUTPUT_FOLDER=self.OUTPUT_FOLDER,
                                      OUTPUT_FILENAME=self.OUTPUT_FILENAME,
                                      eq_constr=self.eq_constr,
                                      optModelB=self.optModelB,
                                      FLAG_FILE=self.FLAG_FILE,)
        opt_obj, opt_u, is_success = c13_solver.solve()
        # opt_obj, opt_u, is_success = nf_fit.solve_multiplestarts_ralg(self.optModelA,self.tol,self.max_iters,10,self.n_restart,self.OUTPUT_FOLDER,
        #                                                                 self.OUTPUT_FILENAME,contol=1e-6,optModelB=self.optModelB,solver=self.solver,
        #                                                                 )

        return opt_obj, opt_u, is_success
    
    def _set_model(self):
        model= Model('scerevisiae')
        model.read_from_file(self.MODEL_FILE)
        return model

    def set_netflux_df(self):
        net_flux_df=nf_io.format_fluxes_tofreeflux2(pd.read_csv(f'{self.FLUX_FILE}',index_col=0),self.modelA) #{FLUXES_FOLDER}flux_samples_{condition}
        net_flux_df.columns=[i for i in range(len(net_flux_df.columns))]
        return net_flux_df

    @property
    def modelA(self):
        return self._modelA 

    @modelA.setter
    def modelA(self,model):
        self._modelA=model
        self.net_flux_df=self.set_netflux_df()

    def _get_label_df(self):
        assert self.dataset in ['mixed','MHLAB','MAREIKE'] 
        if self.dataset=='mixed':
            dataset='mix2labs'
        else : 
            dataset=self.dataset
        return (pd.read_csv(f'{self.MS_FOLDER}/202409_13C_MSdata_{dataset}_noclip_{self.condition}_{self.labeling_strategy}.csv',index_col=0)
                if self.extended_data is False or self.condition!='WT' or self.labeling_strategy!='U-13C' else 
                pd.read_csv(f'{self.MS_FOLDER}/202409_13C_MSdata_{dataset}_noclip_{self.condition}_{self.labeling_strategy}_EXTENDEDSHEN.csv',index_col=0)
                )


    def _set_fit(self,alternative_model=False):
        '''initialize the fitter with the given labelling strategy. This includes setting meas MDV, setting vnet and preparing fitter (EAMs decomposition, etc)
        To init the second fitter in case of combined optimization, set alternative_model=True'''
        fit=None
        fit = self.modelA.fitter('ss') if alternative_model==False else self.modelB.fitter('ss')
        if self.labeling_strategy=='1-13C':
            fit.set_labeling_strategy('glcDe', ['000001'],[1],[0.99])
        elif self.labeling_strategy =='U-13C':
            fit.set_labeling_strategy('glcDe', ['111111'],[0.2],[0.99])
        

        self.df_label_data=self._get_label_df()
        if self.seed is not None:
            nf_fit.helper_set_MDV_bootstrap (fit,nf_fit.data_to_MDV_mix_func(self.df_label_data),seed=self.seed)
        else : 
            nf_fit.helper_set_MDV(fit,nf_fit.data_to_MDV_mix_func(self.df_label_data))
        time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        print(f'>>>> preparing at : {time}')
        nf_fit.helper_prepare(fit,n_jobs=self.n_jobs) #prepare cut in 2 parts for later optimization whselfen running on a large number of samples

        ##after this is for the objective to regress to sample fluxes
        dict_fluxes=nf_io.load_sample_fromdf(self.sample_no,self.net_flux_df)
        
        ### part where we take abs for ex transport fluxes : we turned them around manually in stoichiometric matrix 
        vnet=[
            dict_fluxes[key] if self.modelA.reactions_info[key].reversible
              else abs(dict_fluxes[key]) for key in fit.model.reactions if key in dict_fluxes.keys()]
        fit.vnet=vnet
        
        ##measured fluxes part : 
        print(f'{self.fit_fluxes} fitfluxes')
        if self.fit_fluxes==True:
            df_ex_flux_cond=self.df_ex_flux.loc[self.condition]
            nf_fit.fit_set_measfluxes(fit,df_ex_flux_cond)

        #continue prepare
        if self.flux_bounds == True : #case for false not implemented here
            nf_fit.helper_set_flux_bounds(fit,self.net_flux_df,self.sample_no,1e-3,reset=False)
        elif self.flux_bounds == 'VA':
            
            nf_fit.set_VA_bounds(fit,self.modelA,self.condition) #va bounds file is hardcoded in fit_netfluxes.py
        nf_fit.helper_prep_part2(fit,bound=self.flux_bounds)
        
        time=datetime.now().strftime("%d/%m/%Y %H:%M:%S")    
        print(f"'>>>>>> fitter ready at : {time} ")

        if alternative_model==False:
            self.fitA=fit
        else :
            self.fitB=fit
        
    def _build_optmodel(self,alternative_model=False):
        """uses the build opt model function defined in fit net flux. allows to  differentiate between model A and model B"""
        return (
            nf_fit.build_optModel(self.fitA,self.condition,self.force_flux,fit_fluxes=self.fit_fluxes,sample_n0=self.sample_no,eq_constr=self.eq_constr,
                                  drg_df=self.drg_df) 
            if alternative_model==False
            else nf_fit.build_optModel(self.fitB,fit_fluxes=self.fit_fluxes,sample_n0=self.sample_no)
                )
     #we only use constraints from A so no need to build ff matrixes for optmodelB

