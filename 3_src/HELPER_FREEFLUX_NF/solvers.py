#-------------------------------------------------------------------------------
# Name:        solvers.py

# Purpose:     solver class that can use either openopt or scipy solvers
# 
# 
# Author:      Nathan Fargier
#
# Created:     19/09/2024
### version 1.1 : modified solving function to choose different solvers. Scipy and openopt problems can be defined and solved. 
### version 1.2 :separate solver class from the script with helper/prepare functions. Catch SLURM signal from run_fitting.py with exit_flag.py 
#-------------------------------------------------------------------------------

import sys
import os
# sys.path.append('/home/nathan/thermo_flux/')
# sys.path.append('/home4/p313427/freeflux/')
# sys.path.append('/home2/p313417/freeflux/')
sys.path.append(os.path.expanduser('~/'))
sys.path.append(os.path.expanduser('~/freeflux/'))

import numpy as np
# from freeflux.src.freeflux.solver.nlpsolver import MFAModel
from scipy.linalg import pinv
from openopt import NLP
from scipy.optimize import minimize,LinearConstraint,Bounds
# from . import exit_flag ##to save optim when close to finished


### to solve combined models
def combined_f(optModel1,optModelU):
    def f(u):
        return optModel1.f(u)+optModelU.f(u)
    return f

def combined_df(optModel1,optModelU):
    def df(u):
        return optModel1.df(u)+optModelU.df(u)
    return df

def combined_ddf(optModel1,optModelU):
    def ddf(u):
        return optModel1.ddf(u)+optModelU.ddf(u)
    return ddf

#define globals  for timeout case
timeout_obj=None
timeout_solution=None
timeout_all_lin=None
scalingfactor_f=None
class c13_Solver:
    def __init__(self,
            optModelA,vnet,OUTPUT_FOLDER,OUTPUT_FILENAME,
            tol=1e-6,max_iters=400,disp=True,eq_constr=False,n_restart=4,contol=1e-3,solver='gsubg',scale=1e-3,
            optModelB=None,maxflux=200,FLAG_FILE=None,**kwargs):
        self.optModelA = optModelA
        self.vnet=vnet
        self.optModelB = optModelB
        self.OUTPUT_FOLDER = OUTPUT_FOLDER
        self.OUTPUT_FILENAME = OUTPUT_FILENAME
        self.tol = tol
        self.max_iters = max_iters
        self.disp = disp
        self.eq_constr = eq_constr
        self.n_restarts = n_restart
        self.contol = contol
        self.solver = solver
        ##add global scale so openopt callback can access it 
        self.scale= scale
        global scalingfactor_f 
        scalingfactor_f=self.scale

        #
        self.maxflux=maxflux
        #flag file
        self.FLAG_FILE=FLAG_FILE
        
        ##assign solver method based on solver argument
        if solver in ['ralg','gsubg','mma','auglag','scipy_slsqp','sqlcp']:
            self._solve_standard = self._solve_openopt
        elif solver in ['SLSQP','trust-constr','COBYQA']:
            self._solve_standard = self._solve_scipy
        else:
            raise ValueError(f"Unknown solver: {solver}")
        
        ## then attributes for solving
        self.best_obj_value = float('inf')  # Initialize with a very large number
        self.best_solution = None
        self.best_status = None
        self.restart = 0


    def solve(self):
        ## main solving func
        while self.restart < self.n_restarts and self.restart < 10:
            self.restart += 1
            print(f"Restart {self.restart}/{self.n_restarts}")
            self._build_initial_flux_values_nf(previous_x0=None if self.restart % 2 == 1 else self.best_solution)

            try:
                obj, u, success, all_lin = self._solve_standard(eq_constr=self.eq_constr)

                if all_lin <= self.contol: ## we keep a current sol 
                    self._save_current_solution(u, obj, success)

                if obj < self.best_obj_value: ##and th ebest sol in optu file
                    self._checkbest_sol_save(u, obj, success, all_lin,save=True)

                if obj > 1e4: ## stopped at a bad sol so increase number of restarts #for SLSQP
                    self.n_restarts += 1
                    
            ###for signaled exit add exception here  
            except TimeLimitExit:
                print(f"Time limit exit at restart {self.restart}")
                self._checkbest_sol_save(timeout_solution,timeout_obj,False,timeout_all_lin)
                # print('timeoutsol',timeout_solution)
                # print('Saving best sol') ### TODO here i am  still losing the current iteration u and obj..
                # self._save_best_solution()
                return self.best_obj_value, self.best_solution, self.best_status

            except Exception as e:
                print(f"Error at restart {self.restart}: {e}")
                self.n_restarts += 1
        
        return self.best_obj_value, self.best_solution, self.best_status

    def _build_initial_flux_values_nf(self, ini_netfluxes = None,previous_x0=None,atol=1e-12):  ## absolute tolerance for pinv ### added by NF 01/1
            ''' we set the initial net fluxes to the pseudo-inverse of TN. except if we want to set initial guess as random values'''
            optModel = self.optModelA
            vnet_lb, vnet_ub = np.array(list(optModel.model.net_fluxes_range.values())).T
            vnet_ini = np.random.uniform(low = vnet_lb, high = vnet_ub)
            ### get initial guess from net fluxes
            u_ini = pinv(optModel.T@optModel.N,atol=atol)@vnet_ini
            optModel.x0 = u_ini
            print('max ini x0 coeff',optModel.x0.max())
            if ini_netfluxes =='random':
                optModel.x0 = np.random.rand(u_ini.shape[0])* 400 - 200 ##random values between -200 and 200
            if previous_x0 is not None:
                try :
                    assert previous_x0.shape == optModel.x0.shape
                    optModel.x0 = previous_x0
                    print('starting from previous best sol')
                except:
                    print('incorrect previous_x0')
                    pass

    def _save_current_solution(self,u,obj,success):
        '''save the solution in the output folder'''
        np.savetxt(f'{self.OUTPUT_FOLDER}/start_{self.restart}_sol_{self.OUTPUT_FILENAME}.csv', u, delimiter=',',header=f"{str(obj)} | success {success}")#dont need str
        print(f"Solution saved at {self.OUTPUT_FOLDER}/start_{self.restart}_sol_{self.OUTPUT_FILENAME}.csv") 
        
    def _save_best_solution(self):
        '''save the best solution in the output folder'''
        np.savetxt(f'{self.OUTPUT_FOLDER}/optu_{self.OUTPUT_FILENAME}.csv', self.best_solution, delimiter=',',header=f"{str(self.best_obj_value)} | success {self.best_status}")
        print(f"Best solution saved at {self.OUTPUT_FOLDER}/optu_{self.OUTPUT_FILENAME}.csv")

    def _checkbest_sol_save(self,u,obj,success,all_lin,save=True):
        '''given current obj for u, check if they are the best solution and save the best solution '''
        if obj < self.best_obj_value:
            self.best_obj_value = obj
            self.best_solution = u
            self.best_status = success if success is True else all_lin <= self.contol
        save=save if self.best_status is True else False
        if save:
            self._save_best_solution()

    def _scaled_funcs(self, func_type):
        """
        Returns the scaled objective or derivative function.
        `func_type` should be 'f', 'df', or 'ddf'
        """
        if func_type == 'f':
            base_func = combined_f(self.optModelA, self.optModelB) if self.optModelB else self.optModelA.f
        elif func_type == 'df':
            base_func = combined_df(self.optModelA, self.optModelB) if self.optModelB else self.optModelA.df
        elif func_type == 'ddf':
            base_func = combined_ddf(self.optModelA, self.optModelB) if self.optModelB else self.optModelA.ddf
        else:
            raise ValueError(f"Unknown func_type: {func_type}")

        return lambda u: self.scale * base_func(u)

    class ScipyCallback:
            def __init__(self, optModelA,optModelB=None,callback_freq=20,*args, **kwargs): ## add args and kwargs for trust constr which calls callback with more args
                self.step = 0
                self.optModelA = optModelA
                self.optModelB = optModelB
                self.callback_freq = callback_freq
                self.FLAG_FILE = os.path.expanduser('~/HELPER_FREEFLUX_NF/exit_flag.txt')

            def __call__(self, x, *args, **kwargs):
                global timeout_obj, timeout_solution,  timeout_all_lin
               
                if self.step % self.callback_freq == 0:
                    print(
                        f'at step {self.step}',
                        'min N@u ', (self.optModelA.N @ x).min(),
                        'lin residual', np.where((self.optModelA.A @ x - self.optModelA.b) > 0, self.optModelA.A @ x - self.optModelA.b, 0).sum(),
                        'obj', combined_f(self.optModelA,self.optModelB)(x) if self.optModelB else self.optModelA.f(x) # not scaled, to compare
                        )
                if os.path.exists(self.FLAG_FILE) and open(self.FLAG_FILE).read().strip() == "TRUE":
                    timeout_obj= combined_f(self.optModelA,self.optModelB)(x) if self.optModelB else self.optModelA.f(x)
                    timeout_solution=x
                    timeout_all_lin= np.where((self.optModelA.A @ x - self.optModelA.b) > 0, self.optModelA.A @ x - self.optModelA.b, 0).sum()
                    print(timeout_obj,timeout_all_lin)
                    raise TimeLimitExit('Time limit exit flag received') #return True is not supported by all scipy solvers
                # 
                self.step += 1

    def _solve_scipy(self,eq_constr=False):
        cstr= LinearConstraint(self.optModelA.A, -np.inf, self.optModelA.b)
        if eq_constr is True:
            cstreq={'type': 'eq', 'fun': lambda u: self.optModelA.Aeq@u - self.optModelA.beq} ##this is the constraint Aeqx = beq
        bounds=Bounds(lb=-self.maxflux*np.ones(len(self.optModelA.x0)),ub=self.maxflux*np.ones(len(self.optModelA.x0))) ##bounds for u so that vtot arent unrealistic
        options={'maxiter': self.max_iters, 'disp':self.disp,'gtol':self.contol}
        if self.solver.lower()=='slsqp':
            options['ftol'] = self.tol

        callback_func = self.ScipyCallback(self.optModelA,optModelB=self.optModelB)
        print(f'optModelB is not None {self.optModelB is not None} ')

        ##part to scale the function and gradient to avoid too big steps and numerical issues like positive directionnal linesearch
        func = self._scaled_funcs('f')
        dfunc = self._scaled_funcs('df')
        ddfunc = self._scaled_funcs('ddf')

        ##solve
        res = minimize(
                        fun = func,
                        x0 = self.optModelA.x0,
                        method = self.solver,
                        jac = dfunc, 
                        hess= ddfunc, ## can pass hessian for trust-constr
                        constraints = cstr if eq_constr is False else [cstr,cstreq],
                        bounds=bounds,
                        options = options,
                        callback=callback_func,
                        )
        lin_ineq= (self.optModelA.A@res.x - self.optModelA.b)
        all_lin=lin_ineq[lin_ineq>0].sum() ## sum of all lin constraints that are >0 , violating Ax <= b
        return res.fun, res.x, res.success, all_lin


    class openoptCallback:
        def __init__(self, p,callback_freq=20,FLAG_FILE=None,*args, ** kwargs):
            self.p = p
            self.callback_freq = callback_freq
            self.FLAG_FILE = FLAG_FILE

        def __call__(self, *args, **kwargs):
            if self.p.iter % self.callback_freq == 0:
                print(
                    f'at step {self.p.iter}',
                    'lin residual',self.p.point(self.p.xk).all_lin(),
                    'obj', self.p.fk*(1/scalingfactor_f) #scaling using attributes deefined in NLP
                    )
            if self.FLAG_FILE and os.path.exists(self.FLAG_FILE):
                with open(self.FLAG_FILE) as f:
                    flag_value = f.read().strip()
                if flag_value == "TRUE":
                    global timeout_obj, timeout_solution, timeout_all_lin
                    timeout_obj = self.p.fk * (1 / scalingfactor_f)
                    timeout_solution = self.p.xk
                    timeout_all_lin = self.p.point(self.p.xk).all_lin()
                    print(timeout_obj, timeout_all_lin)
                    raise TimeLimitExit('Time limit exit flag received')

                

            return False

    def _solve_openopt(self,eq_constr=False):
        func=self._scaled_funcs('f')
        dfunc=self._scaled_funcs('df')
        p  = NLP(
                f = func,
                x0 = self.optModelA.x0, 
                df =  dfunc ,
                A = self.optModelA.A,
                b = self.optModelA.b,
                Aeq=None if eq_constr is False else self.optModelA.Aeq,
                beq=None if eq_constr is False else self.optModelA.beq,
                lb=-self.maxflux*np.ones(len(self.optModelA.x0)), ## add bounds for u so that vtot arent too large
                ub=self.maxflux*np.ones(len(self.optModelA.x0)),
                xtol = self.tol, 
                ftol = self.tol, 
                fTol=self.tol,
                maxIter = self.max_iters, 
                iprint = 10 if self.disp else -1,
                callback_freq=50,
                contol=self.contol,
                scale=self.scale,
                output_folder=self.OUTPUT_FOLDER,output_filename=self.OUTPUT_FILENAME
            )
        callbackfunc=self.openoptCallback(p,FLAG_FILE=self.FLAG_FILE,callback_freq=50)
        p.callback=callbackfunc

        res = p.solve(self.solver,maxShoots=25)#gsubg :default maxshoots is 15 
        ## we have to unscale the obj 
        return res.ff*(1/self.scale), res.xf, res.istop > 0  , p.point(res.xf).all_lin() ##last one is lin ineq res

class TimeLimitExit(Exception):
    """Custom exception to indicate a time limit exit."""
    pass

# def solve_multiplestarts_ralg(optModelA,tol,max_iters,disp,n_restarts,OUTPUT_FOLDER,OUTPUT_FILENAME,contol=1e-6,optModelB=None,solver='gsubg',**kwargs):
#     '''first we define callback and solving function with callback
#     Then we call the solver for the number of restarts'''
#     def callback_ralg(p): 
#         if p.iter %p.callback_freq==0:
#             print('lin residual',p.point(p.xk).all_lin())
#             print('min N@u ',(optModelA.N@p.xk).min())
#             # try : #if we have p.c ,don't have to do an extra checking operation with if
#             #     print('flux f constraint',p.c(p.xk).max())
#             # except:
#             #     pass
#             # if p.point(p.xk).all_lin()> p.contol: # p.rk is current iter max residual
#             #     print('--= non-feasible iter =--')
#             #     print( 'itn:',  p.iter,p.contol)
#             #     return False
#             # elif p.point(p.xk).all_lin() <=p.contol and p.fk<=1e2: # i couldnt find a way to assess |f(k)-f(k-1)|<tol  here
#             #     print('sol ok')
#             #     np.savetxt(f'{p.output_folder}/start_{p.restart_n}ite_{p.iter}_sol_{p.output_filename}.csv', p.xk, delimiter=',',header=str(p.fk))
#             #     return True
#             # else: 
#             # print('constraints ok but sol keeps going')
#             # np.savetxt(f'{p.output_folder}/start_{p.restart_n}ite_{p.iter}_sol_{p.output_filename}.csv', p.xk, delimiter=',',header=str(p.fk))
#             return False      
#         return False
    
#     def _solve_flux_openopt_callback(optModelA,tol,max_iters,disp,optModelB=None,solver=solver,**kwargs):
#         p  = NLP(
#                 f =combined_f(optModelA,optModelB) if optModelB else optModelA.f, 
#                 x0 = optModelA.x0, 
#                 df = combined_df(optModelA,optModelB) if optModelB else optModelA.df, 
#                 A = optModelA.A, 
#                 b = optModelA.b, 
#                 xtol = tol, 
#                 ftol = tol, 
#                 fTol=tol,
#                 maxIter = max_iters, 
#                 iprint = 10 if disp else -1,
#                 callback=callback_ralg,
#                 callback_freq=50,
#                 contol=contol,
#                 output_folder=OUTPUT_FOLDER,output_filename=OUTPUT_FILENAME,**kwargs ##kwargs can be used to pass c and dc ?
#             )
        
#         res = p.solve(solver,maxShoots=25) ##changed to gsubg ##070125 ## way faster than ralg #default maxshoots is 15 
	
#         return res.ff, res.xf, res.istop > 0  , p.point(res.xf).all_lin() ##last one is lin ineq res

#     def _solve_flux_scipy_callback(optModelA,tol,max_iters,disp,optModelB=None,solver=solver,**kwargs):
#         # cstr={'type': 'ineq', 'fun': lambda u: -optModelA.A@u + optModelA.b}
#         cstr=LinearConstraint(optModelA.A, -np.inf, optModelA.b) ##this is the constraint Ax <= b , replacing the one above
#         bounds=Bounds(lb=-500*np.ones(len(optModelA.x0)),ub=500*np.ones(len(optModelA.x0))) ##bounds for u so that vtot arent too large
#         options={'maxiter': max_iters, 'disp':disp}
#         if solver=='slsqp' or solver.lower()=='slsqp':
#             options['ftol'] = tol

#         class ScipyCallback:
#             def __init__(self, optModelA,optModelB=None,callback_freq=20,*args, **kwargs):
#                 self.step = 0
#                 self.optModelA = optModelA
#                 self.optModelB = optModelB
#                 self.callback_freq = callback_freq

#             def __call__(self, x, *args, **kwargs):
#                 if self.step % self.callback_freq == 0:
#                     print(
#                         f'at step {self.step}',
#                         'min N@u ', (self.optModelA.N @ x).min(),
#                         'lin residual', np.where((self.optModelA.A @ x - self.optModelA.b) > 0, self.optModelA.A @ x - self.optModelA.b, 0).sum(),
#                         'obj', combined_f(self.optModelA,self.optModelB)(x) if self.optModelB else self.optModelA.f(x)
#                         )
#                 self.step += 1
#                 # np.savetxt(f'{OUTPUT_FOLDER}/start_{restart}_sol_{OUTPUT_FILENAME}.csv',x,  delimiter=',',header=str(optModelA.f(x)))
#         callback_func = ScipyCallback(optModelA,optModelB=optModelB)

#         ##part to scale the function and gradient to avoid too big steps and numerical issues like positive directionnal linesearch
#         scale=1e-4
#         def func(u):
#             return scale*combined_f(optModelA,optModelB)(u) if optModelB else scale*optModelA.f(u)
#         def dfunc(u):
#             return scale*combined_df(optModelA,optModelB)(u) if optModelB else scale*optModelA.df(u)
#         def ddfunc(u):
#             return scale*combined_ddf(optModelA,optModelB)(u) if optModelB else scale*optModelA.ddf(u)
        
#         ##solve
#         res = minimize(
#                         fun = func,#combined_f(optModelA,optModelB) if optModelB else optModelA.f,
#                         x0 = optModelA.x0,
#                         method = solver,
#                         jac = dfunc, # combined_df(optModelA,optModelB) if optModelB else optModelA.df,
#                         hess= ddfunc, #combined_ddf(optModelA,optModelB) if optModelB else optModelA.ddf,
#                         constraints = cstr,
#                         bounds=bounds,
#                         options = options,
#                         callback=callback_func,
#                         )
#         lin_ineq= (optModelA.A@res.x - optModelA.b)
#         all_lin=lin_ineq[lin_ineq>0].sum() ## sum of all lin constraints that are >0 , violating Ax <= b
#         return res.fun, res.x, res.success, all_lin
        
#     ####main part
#     best_obj_value = float('inf')  # Initialize with a very large number
#     best_solution = None
#     best_status = None

#     restart =0
#     while restart < n_restarts and restart <10: #otherwise infinite loop if pb 
#         restart += 1
#         print(f"Restart {restart}/{n_restarts}")
#         ###start 
#         ##initialize x0 from net fluxes except for even restarts where we use the previous best opt u
#         build_initial_flux_values_nf(optModelA,ini_netfluxes=None,previous_x0=None if restart == 1 or restart % 2==1  else best_solution) 
#         ### if restart == 0 else 'random') ## we could do first restart with x0=pinv(TN)@vnet_ini , then random, but initial guess is always better

#         try : ##try except in case of "SVD didnt converge"
#             if solver in ['ralg','gsubg','mma','auglag','sqlcp','scipy_slsqp']:
#                 opt_obj, opt_u, is_success ,all_lin = _solve_flux_openopt_callback(optModelA,tol, max_iters, disp,optModelB=optModelB,**kwargs)

#             elif solver in ['SLSQP','trust-constr','COBYQA']:
#                 opt_obj, opt_u, is_success, all_lin = _solve_flux_scipy_callback(optModelA,tol, max_iters, disp,optModelB=optModelB,**kwargs)

#             else:
#                 raise ValueError(f"Unknown solver: {solver}")
            
#             if  all_lin <= contol: ###save the restarts sol ## if success or if smallenough residual..
#                 # ##if it is the best residual it will be saved as optu anyways
#                 np.savetxt(f'{OUTPUT_FOLDER}/start_{restart}_sol_{OUTPUT_FILENAME}.csv',opt_u,  delimiter=',',header=f"{str(opt_obj)} | success {is_success}")

#         # If this restart gives a better objective value, update the best solution
#             if opt_obj < best_obj_value:
#                 best_obj_value = opt_obj
#                 best_solution = opt_u
#                 best_status = is_success if is_success else all_lin <= contol ## if hits maxiters success is given as False but the solution can be feas
         
#         except Exception as e:
#             print(f"Error in restart {restart + 1}/{n_restarts}: {e}")
#             n_restarts += 1  # more max restarts if an error occurs
#             continue
#     return best_obj_value, best_solution, best_status

