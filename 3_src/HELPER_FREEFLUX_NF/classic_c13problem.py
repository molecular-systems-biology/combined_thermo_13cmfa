#-------------------------------------------------------------------------------
# Name:        classic_c13problem_class.py

# Purpose:     variation of c13problem class to solve regular 13C fitting and use parameter continuation to build Confidence INtervals
# Class to build (freeflux) and solve the 13C fitting problem with NLP solver of open opt.
#  The class allows for combined optimization of 2 models (with U-13C and 1-13C labeling strategies)
#  Uses functions from HELPER_FREEFLUX_NF.fit_netfluxes
# 
# 
# Author:      Nathan Fargier
#
# Created:     06/03/2024

##version 1.1 Confidence intervals with nlp_opt_singleflux
###not updated to use solvers.py 
#-------------------------------------------------------------------------------

import sys
import os 
sys.path.append('/home4/p313427/freeflux/')
sys.path.append('/home4/p313427/')
#sys.path.append(os.path.abspath('../'))
from openopt import NLP
from scipy.optimize import approx_fprime

import pandas as pd
import ast
from freeflux.src.freeflux import Model
from HELPER_FREEFLUX_NF import io as nf_io
import numpy as np
from HELPER_FREEFLUX_NF import fit_netfluxes as nf_fit
from datetime import datetime
from HELPER_FREEFLUX_NF import c13problem_class as c13

### deactivate warnings 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


class classic_c13_problem(c13.c13_problem):
    def __init__(
                self, condition: str, labeling_strategy: str, sample_no: int, MODEL_FILE: str, FLUX_FILE: str,
                MS_FOLDER: str, OUTPUT_FOLDER: str, OUTPUT_FILENAME: str, dataset= 'mixed', tol=1e-6, max_iters=200, n_restart=1, n_jobs=8, 
                flux_bounds=False, df_ex_flux=None,fit_fluxes=True, ##need all args necesseray for fit if we want to use the same class for the fit
                optu=None, flux_idx =0, direction =1 , confidence_level=3.84, max_iterations=1000
                ):
        ###flux bounds ==false makes freeflux perform fva and use these bounds instead of the sample bounds
        super().__init__(
                        condition=condition,labeling_strategy= labeling_strategy, dataset=dataset,sample_no= sample_no, MODEL_FILE=MODEL_FILE, 
                        FLUX_FILE=FLUX_FILE, MS_FOLDER=MS_FOLDER,
                        OUTPUT_FOLDER=OUTPUT_FOLDER, OUTPUT_FILENAME=OUTPUT_FILENAME,tol= tol, max_iters=max_iters, n_restart=n_restart, n_jobs=n_jobs,
                        flux_bounds=flux_bounds,df_ex_flux=df_ex_flux,fit_fluxes=fit_fluxes,force_flux=False
                         )  # Initialize parent class

        # Additional attributes specific to the optimizer
        self.optu = optu ## has to be set after the first fit (call of classic_c13problem)
        self.confidence_level = confidence_level
        self.max_iterations = max_iterations
        self.direction = direction  
        self.flux_idx = flux_idx
        self.best_rss = None

    def __call__(self):
       #call super to fit the model
        opt_obj, opt_u, is_success = super().__call__()
        self.optu = opt_u #self.optu vs opt_u # confusing ?
        self.best_rss= opt_obj
        return opt_obj, opt_u, is_success

    def nlp_opt_singleflux(self):
            """
            Brute force non linear search : min or maximize a given flux given a non-linear constraint, which is built using the objective function of the previous fit. 
            The derivative of the constraint is approximated with scipy approx_fprime.
            returns obj_value, optimized_vars, is_success, sum_of_all_constrs
            """ 
            #init vars and get best rss
            self.best_rss = self._c13_get_rss() if self.best_rss is None else self.best_rss
            self.confidence_level=self.confidence_level
            #objective function designed to work for on flux idx per call
            def _build_f(c13_problem):
                def f(u): ##get vnet var
                    return (c13_problem.modelA.transform_matrix@c13_problem.optModelA.N@u)[self.flux_idx]
                return f
            def _build_df(f):
                def df(u):
                    return approx_fprime(u,f)
                return df
            f=_build_f(self)
            df=_build_df(f)

            ##now constraints
            c=self._build_c()
            dc=self.combined_dc(self.optModelA,self.optModelB)

            #formulate problem with NLP openopt
            p = NLP(f=f,df=df,c=c,dc=dc,x0=self.optu,
                    A=-self.optModelA.N,b=np.zeros(self.optModelA.ntotalfluxes),
                                fTol=self.tol,xtol=self.tol,#fEnough=1e-1,
                                maxIter = self.max_iterations,
                                output_folder=self.OUTPUT_FOLDER,output_filename=self.OUTPUT_FILENAME,
                                callback=self.callback_ralg, callback_freq=10,best_rss=self.best_rss,
            )

            p.goal = 'max' if self.direction == 1 else p.goal #set dir

            res = p.solve('gsubg',maxShoots=15) 

            print(res.ff, res.xf, res.istop>0)
            is_success = res.istop > 0
            optimized_vars = res.xf
            
            return res.ff, optimized_vars, is_success , p.point(res.xf).sum_of_all_active_constraints()


    def _build_c(self):
        '''nonlinear constraint : RSS(U) - RSS(opt_u) - confidence_level <= 0'''
        def c(u):
            curr_rss=self._c13_get_rss(u=u)
            return curr_rss-self.best_rss-self.confidence_level
        return c

    @staticmethod
    def combined_dc(optModel1,optModelU):
        '''derivative of c is just d(rss)/du'''
        optModel1.build_gradient()
        optModelU.build_gradient()
        def dc(u): ## here the derivative of constraint is the derivative of the preivous obj func(RSS(U))
            return optModel1.df(u)+optModelU.df(u)
        return dc

    def _c13_get_rss(self,u=None):
        if u is not None:
            self.optu=u
        sim_MDVs_A,sim_MDVs_B=self._u_to_SIMMDV()
        rss = self.calc_rss(sim_MDVs_A,self.fitA)+self.calc_rss(sim_MDVs_B,self.fitB)
        return rss

    def _u_to_SIMMDV(self):
        vtot=self.optModelA.N@self.optu
        self.fitA.model.total_fluxes[:]=vtot
        self.fitB.model.total_fluxes[:]=vtot
        sim_MDVs_A_all = self.fitA.calculator._calculate_MDVs()
        sim_MDVs_A = {emuid: sim_MDVs_A_all[emuid] for emuid in self.fitA.model.target_EMUs}

        sim_MDVs_B_all = self.fitB.calculator._calculate_MDVs()
        sim_MDVs_B = {emuid: sim_MDVs_B_all[emuid] for emuid in self.fitB.model.target_EMUs}
        return sim_MDVs_A,sim_MDVs_B

    @staticmethod
    def calc_rss(sim_MDVs,fit):
        rss=0
        def least_square_error(simulated,experimental,std_exp):
            """all inputs are np arrays of the same size
            here we calculate the absolute error otherwise too high values with squared"""
            ##replace 0 with 1 in std_exp
            std_exp=np.where(std_exp==0,1,std_exp)
            ##now we can calculate the least square error, absolute not squared
            error=np.sum((abs(simulated-experimental)/std_exp))
            return error
        for emuid in fit.model.target_EMUs:
            rss+=least_square_error(sim_MDVs[emuid],fit.model.measured_MDVs[emuid][0],fit.model.measured_MDVs[emuid][1])
        return rss
    
    @staticmethod
    def callback_ralg(p): 
        if p.iter %p.callback_freq==0:
            print('lin residual',p.point(p.xk).all_lin())
            try : #if we have p.c ,don't have to do an extra checking operation with if
                print('NLP obj constraint',p.c(p.xk).max())
            except:
                pass
            if p.point(p.xk).all_lin()> p.contol: # p.rk is current iter max residual
                print('--= non-feasible iter =--')
                print( 'itn:',  p.iter,p.contol)
                return False
            elif p.point(p.xk).all_lin() + p.c(p.xk).sum() <=p.contol and abs(p.best_rss-p.fk) >= 1000/(p.iter+1)  : ### increase chance to stop with iterations
                print('sol ok')
                np.savetxt(f'{p.output_folder}/optu_{p.output_filename}.csv',p.xk,delimiter=',',header=str(p.fk)) #same name so erased by final opt u
                return False #just save
            else: 
                print('lin constraints ok but not nlp')
                # np.savetxt(f'{p.output_folder}/belowtol_start_{p.restart_n}ite_{p.iter}_sol_{p.output_filename}.csv', p.xk, delimiter=',',header=str(p.fk))
                return False      
        return False
    

