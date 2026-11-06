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

##version 1.2 Confidence intervals with nlp_opt_singleflux
##  - squared weighted residuals (chi-square) instead of absolute
##  - std scaled by sqrt(wss_min/(n-p)) if the best fit is worse than expected
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
                optu=None, flux_idx =0, direction =1 , confidence_level=3.84, max_iterations=1000,FLAG_FILE=None
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
        self.sigma_scale = 1 ## set by _set_sigma_scale
        self.scaled = False ## best_rss already scaled or not
        self.FLAG_FILE = FLAG_FILE

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

            ##scale std and best_rss
            # self._set_sigma_scale()

            ##now constraints
            c=self._build_c()
            dc=self.combined_dc(self.optModelA,self.optModelB,scale=self.sigma_scale)

            #formulate problem with NLP openopt
            p = NLP(f=f,df=df,c=c,dc=dc,x0=self.optu,
                    A=-self.optModelA.N,b=np.zeros(self.optModelA.ntotalfluxes),
                                fTol=self.tol,xtol=self.tol,#fEnough=1e-1,
                                maxIter = self.max_iterations,
                                output_folder=self.OUTPUT_FOLDER,output_filename=self.OUTPUT_FILENAME,
                                best_rss=self.best_rss,
            )

            p.goal = 'max' if self.direction == 1 else p.goal

            cb = openoptCallback(p, callback_freq=10, FLAG_FILE=self.FLAG_FILE)
            p.callback = cb

            try:
                res = p.solve('gsubg', maxShoots=15)
            except TimeLimitExit:
                print(f'Time limit exit: flux {self.flux_idx}, '
                      f'direction {self.direction}')
                np.savetxt(
                    f'{self.OUTPUT_FOLDER}/interrupted_{self.OUTPUT_FILENAME}.csv',
                    cb.timeout_solution, delimiter=',',
                    header=str(cb.timeout_obj))
                # is_success=False: interrupted, not a converged bound
                return (cb.timeout_obj, cb.timeout_solution, False,
                        cb.timeout_constr)

            print(res.ff, res.istop, res.msg)
            is_success = res.istop > 0
            return (res.ff, res.xf, is_success,
                    p.point(res.xf).sum_of_all_active_constraints())

    # def _set_sigma_scale(self):
    #     '''scale std by sqrt(wss_min/(n-p)) over both datasets, never below 1.
    #     wss_min is the best_rss from the fit (unscaled std); best_rss is then divided by scale**2'''
    #     if self.scaled: ## only once, otherwise best_rss is divided again
    #         return
    #     wss_min = self.best_rss  # unscaled, from the fit
    #     n = sum(np.size(fit.model.measured_MDVs[emuid][0])
    #             for fit in (self.fitA, self.fitB) for emuid in fit.model.target_EMUs)
    #     p = len(self.optu)
    #     self.sigma_scale = max(np.sqrt(wss_min / (n - p)), 1)
    #     self.best_rss = wss_min / self.sigma_scale**2  # scaled
    #     self.scaled = True
    #     print(f"wss_min={wss_min}, n={n}, p={p}, scale={self.sigma_scale}")

    def _build_c(self):
        '''nonlinear constraint : RSS(U) - RSS(opt_u) - confidence_level <= 0'''
        def c(u):
            curr_rss=self._c13_get_rss(u=u)
            return curr_rss-self.best_rss-self.confidence_level
        return c

    @staticmethod
    def combined_dc(optModel1,optModelU,scale=1):
        '''derivative of c is just d(rss)/du, divided by scale**2 since std is scaled by scale'''
        optModel1.build_gradient()
        optModelU.build_gradient()
        def dc(u): ## here the derivative of constraint is the derivative of the preivous obj func(RSS(U))
            return (optModel1.df(u)+optModelU.df(u))/scale**2
        return dc

    def _c13_get_rss(self,u=None):
        if u is not None:
            self.optu=u
        sim_MDVs_A,sim_MDVs_B=self._u_to_SIMMDV()
        rss = (self.calc_rss(sim_MDVs_A,self.fitA,scale=self.sigma_scale)
               +self.calc_rss(sim_MDVs_B,self.fitB,scale=self.sigma_scale))
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
    def calc_rss(sim_MDVs, fit, method='squared', zero_std=1, scale=1):
        """Weighted residual error between simulated and measured MDVs.

        Parameters
        ----------
        sim_MDVs : dict
            Simulated MDVs, keyed by EMU id.
        fit : Fitter
            FreeFlux fitter holding the target EMUs and measured MDVs (mean, std).
        method : {'squared', 'absolute'}
            'squared'  : sum(((sim - exp) / std)**2), i.e. -2 log-likelihood for
                         Gaussian errors (use with the chi-square threshold, e.g. 3.84).
            'absolute' : sum(|sim - exp| / std) (previous behaviour).
        zero_std : float
            Value used to replace standard deviations equal to 0.
        scale : float
            Factor by which all standard deviations are multiplied.
        """
        if method not in ('squared', 'absolute'):
            raise ValueError("method must be 'squared' or 'absolute'")

        def weighted_error(simulated, experimental, std_exp):
            simulated = np.asarray(simulated, dtype=float)
            experimental = np.asarray(experimental, dtype=float)
            std_exp = np.asarray(std_exp, dtype=float)
            std_exp = np.where(std_exp == 0, zero_std, std_exp) * scale
            residuals = (simulated - experimental) / std_exp
            if method == 'squared':
                return np.sum(residuals**2)
            return np.sum(np.abs(residuals))

        rss = 0.0
        for emuid in fit.model.target_EMUs:
            measured_mean, measured_std = fit.model.measured_MDVs[emuid][0], fit.model.measured_MDVs[emuid][1]
            rss += weighted_error(sim_MDVs[emuid], measured_mean, measured_std)
        return rss
    
    # @staticmethod
    # def callback_ralg(p): 
    #     if p.iter %p.callback_freq==0:
    #         print('lin residual',p.point(p.xk).all_lin())
    #         try : #if we have p.c ,don't have to do an extra checking operation with if
    #             print('NLP obj constraint',p.c(p.xk).max())
    #         except:
    #             pass
    #         if p.point(p.xk).all_lin()> p.contol: # p.rk is current iter max residual
    #             print('--= non-feasible iter =--')
    #             print( 'itn:',  p.iter,p.contol)
    #             return False
    #         elif p.point(p.xk).all_lin() + p.c(p.xk).sum() <=p.contol and abs(p.best_rss-p.fk) >= 1000/(p.iter+1)  : ### increase chance to stop with iterations
    #             print('sol ok')
    #             np.savetxt(f'{p.output_folder}/optu_{p.output_filename}.csv',p.xk,delimiter=',',header=str(p.fk)) #same name so erased by final opt u
    #             return False #just save
    #         else: 
    #             print('lin constraints ok but not nlp')
    #             # np.savetxt(f'{p.output_folder}/belowtol_start_{p.restart_n}ite_{p.iter}_sol_{p.output_filename}.csv', p.xk, delimiter=',',header=str(p.fk))
    #             return False      
    #     return False
class openoptCallback:
    def __init__(self, p, callback_freq=10, FLAG_FILE=None):
        self.p = p
        self.callback_freq = callback_freq
        self.FLAG_FILE = FLAG_FILE
        # state at interruption
        self.timeout_obj = None
        self.timeout_solution = None
        self.timeout_constr = None

    def _flag_set(self):
        if not self.FLAG_FILE or not os.path.exists(self.FLAG_FILE):
            return False
        with open(self.FLAG_FILE) as f:
            return f.read().strip() == "TRUE"

    def __call__(self, *args, **kwargs):
        p = self.p
        # check the flag at every iteration
        if self._flag_set():
            self.timeout_obj = p.fk
            self.timeout_solution = np.copy(p.xk)
            self.timeout_constr = p.point(p.xk).sum_of_all_active_constraints()
            print('exit flag:', self.timeout_obj, self.timeout_constr)
            raise TimeLimitExit('Time limit exit flag received')

        if p.iter % self.callback_freq == 0:
            all_lin = p.point(p.xk).all_lin()
            c_val = np.max(p.c(p.xk))
            print(f'at step {p.iter}', 'lin residual', all_lin,
                  'NLP constraint', c_val, 'flux', p.fk)
            # save the current point if feasible
            if all_lin <= p.contol and c_val <= p.contol:
                np.savetxt(f'{p.output_folder}/optu_{p.output_filename}.csv',
                           p.xk, delimiter=',', header=str(p.fk))
        return False

class TimeLimitExit(Exception):
    """Custom exception to indicate a time limit exit."""
    pass