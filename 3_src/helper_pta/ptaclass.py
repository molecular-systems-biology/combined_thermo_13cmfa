#-------------------------------------------------------------------------------
# Name:        helper_pta.py

# Purpose:     Class to build pta solver to find initial point and samples thermodynamic space on cluster
# 
# 
# Author:      Nathan Fargier
#
# Created:     21/10/2025

##version 1.0
#-------------------------------------------------------------------------------

import sys
import os 
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats.distributions import chi2

# cobra
import cobra
from cobra.flux_analysis.variability import find_blocked_reactions
from cobra.util.array import create_stoichiometric_matrix
# pta
import pta
from pta.sampling import tfs
from pta.constants import (
    default_max_psrf,
    default_max_threads,
    default_min_eigenvalue_tds_basis,
    tfs_default_feasibility_cache_size,
    tfs_default_min_rel_region_length,
)
##pta utils for prep
from pta.utils import ceil,floor

##for thermodynamicspace basis
from component_contribution.linalg import LINALG
import numpy.linalg as la
import scipy.linalg
from pta.utils import covariance_square_root
from enkie.constants import R

#thermodynamic space 
from enkie.utils import qrvector, qvector
from pta.commons import Q

sys.path.append(os.path.expanduser('~/'))
import cvxpy as cp
##own cm
from helper_pta import conv_manager_mod as cm
##sectors 
from sampler import helper_sampling as hs

class preprocess:
    '''steps : 
    1. process codel and stoichiometric matrix : remove blocked reactions and not used mets. save rxn,met mappings
    2. get idx of rxn for which second law is ignored
    3. load thermo space data : v,lnc and  drg: 
        a. v is trivial
        b. a lognormal distribution is made from lncbounds 
        c. drG0 and Q(cov matrix) are loaded from equilibrator saved files '''
    def __init__(self,cobra_file,rign_file,vbound_file,
                 logconcmean_file,logconcvar_file,lncbounds_file,
                 drG0file, drGcovsqrtfile,drG0covfile=None,##drg0cov not sqrt used for EP
                 ):
        data_bounds_flux=pd.read_csv(vbound_file,index_col=0)
        self.model=cobra.io.read_sbml_model(cobra_file)
        drg0_prime_mean_init=pd.read_csv(drG0file,index_col=0).values.T[0]
        drg0_prime_cov_sqrt_init =pd.read_csv(drGcovsqrtfile,index_col=0).values
        if drG0covfile is not None:
            drg0_prime_cov_init = pd.read_csv(drG0covfile,index_col=0).values


        #step1 
        ## remove Mg metabolite
        initial_mets_id = [m.id for m in self.model.metabolites]
        self.model.metabolites.remove(self.model.metabolites.Mg_c)
        self.model.metabolites.remove(self.model.metabolites.Mg_e)
        self.model.metabolites.remove(self.model.metabolites.Mg_m)

        #rxn blocked
        ## reduce self.model bounds to the flux VA 
        # tol=1e-10 ##otherwise blocked rxn dont have 0 but 1e-1x
        # for i,row in data_bounds_flux.iterrows():
        #         lb,ub = row
        #         if abs(lb)<tol: lb=0
        #         if abs(ub)<tol: ub=0
        #         self.model.reactions.get_by_id(i).bounds = (lb,ub)
        self.add_v_margin(data_bounds_flux)
        print('added margins to flux bounds')#,np.max([r.bounds for r in self.model.reactions]),np.min([r.bounds for r in self.model.reactions]))
        # self.model.reactions.get_by_id("ADPATPtm-4-2").bounds = (0.0, 0.0)
        # print('ADPATPtm-4-2 bounds',self.model.reactions.get_by_id("ADPATPtm-4-2").bounds)


        ##same function than pta to add margins
        #check if reactions are blocked
        ##keep inti reactions
        init_reactions_id=[r.id for r in self.model.reactions]
        id_removed = find_blocked_reactions(self.model)
        rxn_removed = {old_idx:rid for old_idx,rid in enumerate(init_reactions_id) if rid in id_removed}
        ridx_notremoved = [
            oldidx for oldidx in range(len(init_reactions_id)) 
            if oldidx not in rxn_removed.keys()    
        ]
        old_new_r_mapping={old:new for new,old in enumerate(ridx_notremoved)}
        self.model.remove_reactions(find_blocked_reactions(self.model),remove_orphans=True)
        print(f'{len(id_removed)} reactions to remove',id_removed)
        ##mets
        ##we also removed metabolites, get their old indices
        mets_removed = {old_idx:mid for old_idx,mid in enumerate(initial_mets_id) if mid not in [m.id for m in self.model.metabolites]}
        met_idx_notremoved = [
            oldidx for oldidx in range(len(initial_mets_id)) 
            if  oldidx not in mets_removed.keys()
        ]
        S=create_stoichiometric_matrix(self.model)
        self.model=self.model ## ???

        ###2. index of reactions for which second law is ignored
        # snd_ignored=[181,182,183,184,185,186, 187, 188, 189,190, 191, 192, 193,194, 195, 196, 197,198, 218, 259, 280,281]
        # pos_ht_rxn=[205,206,207,208] ##posqt and HT are not ignore but fully defined by drG0
        # all_ignored=snd_ignored+pos_ht_rxn
        # snd_ignored_newidx=[old_new_r_mapping[old] for old in all_ignored if old in ridx_notremoved]

        # Reactions for which the second law is ignored (made a series with rxn id : flag for ignore_snd + Ht,Htm,POSQt,POSQtm
        ser_rxn_ign=pd.read_csv(rign_file,index_col=0)['0']
        all_ignored=ser_rxn_ign[ser_rxn_ign==True].index

        snd_ignored_newidx = [self.model.reactions.index(rid) for rid in all_ignored.intersection([r.id for r in self.model.reactions])]
        self.snd_ignored_newidx=snd_ignored_newidx

        rxn_constrained_idx=[ridx for ridx in range(len(self.model.reactions)) if ridx not in snd_ignored_newidx]
        map_constrained_ridx={cstridx:newidx for cstridx,newidx in enumerate(rxn_constrained_idx) }
        ##finally Sconstrained
        self.Sconstrained=S[:, [idx for idx in range(len(self.model.reactions)) if idx not in snd_ignored_newidx]]
        self.rid_constrained=[self.model.reactions[idx].id for idx in range(len(self.model.reactions)) if idx not in snd_ignored_newidx]

        ##3.load vars
        #lnc ~ a gaussian distribution
        log_conc_mean=pd.read_csv(logconcmean_file,index_col=0)['0'].values
        log_conc_var=pd.read_csv(logconcvar_file,index_col=0)['0'].values
        log_conc_cov=np.diag(log_conc_var)
        ##metabolites are already ignored in ep apprxo
        self.log_conc_cov=log_conc_cov#[met_idx_notremoved,:][:,met_idx_notremoved]
        self.log_conc_mean=log_conc_mean#[met_idx_notremoved]
        ### we also have to load lcbounds to see which metabolites have lnc =0 
        bounds_lnconc=pd.read_csv(lncbounds_file,index_col=0).loc[[m.id for m in self.model.metabolites]]
        ignore_met_id=bounds_lnconc[bounds_lnconc['ub']-bounds_lnconc['lb']==0].index
        ignoreconc_met_idx=[self.model.metabolites.index(met_id) for met_id in ignore_met_id]
        self.ignoreconc_met_idx=ignoreconc_met_idx

        ## load drG 
        idx_constrained_old=[idx for idx in range(len(drg0_prime_mean_init)) if idx not in all_ignored and idx in ridx_notremoved]
        self.drg0_prime_mean=drg0_prime_mean_init[idx_constrained_old]
        self.drg0_prime_cov_sqrt=drg0_prime_cov_sqrt_init[idx_constrained_old,:][:,idx_constrained_old]
        if drG0covfile is not None:
            self.drG0_prime_cov=drg0_prime_cov_init[idx_constrained_old,:][:,idx_constrained_old]

    #     ##try clip 
    #     if self.threshold_drG0 is not None:
    #         self.drg0_prime_cov_sqrt[abs(self.drg0_prime_cov_sqrt)<self.threshold_drG0]=0
    # # def numerics_threshold(self,threshold_conc,threshold_drG):
    # #     for 

    # @staticmethod
    ##dont use this -> too low probability to get close to bounds 
    # def lognormal_from_ln_bounds(lnL_mM, lnU_mM, coverage=0.99, overshoot=1):
    #     """
    #     Compute mu and sigma for a lognormal distribution given ln concentrations in mM,
    #     convert to M, and optionally expand bounds to avoid exact-bound issues.

    #     lnL_mM, lnU_mM: bounds in natural log (ln) space
    #     coverage: probability coverage of the interval
    #     overshoot: factor >1 to expand the interval proportionally
    #     """
    #     # Convert ln(mM) -> ln(M): ln(M) = ln(mM * 1e-3) = ln(mM) + ln(1e-3)
    #     lnL_M = lnL_mM + np.log(1e-3)
    #     lnU_M = lnU_mM + np.log(1e-3)

    #     # Expand bounds linearly in ln-space
    #     center = (lnL_M + lnU_M) / 2
    #     half_width = (lnU_M - lnL_M) / 2
    #     lnL_expanded = center - overshoot * half_width
    #     lnU_expanded = center + overshoot * half_width

    #     # z-score for the given coverage
    #     z = norm.ppf(0.5 + coverage/2)

    #     # Parameters of underlying normal
    #     mu = (lnL_expanded + lnU_expanded) / 2
    #     sigma = (lnU_expanded - lnL_expanded) / (2 * z)

    #     return mu, sigma

    def add_v_margin(self,data_bounds_flux,    margin: float = 1e-4,    round_to_digits: int = 6,):
        for r in self.model.reactions:
            lb = data_bounds_flux.loc[r.id, "lb"]
            ub = data_bounds_flux.loc[r.id, "ub"]

            if 0 <= lb <= margin:
                r.lower_bound = max(r.lower_bound, 0.0)
            else:
                r.lower_bound = max(r.lower_bound, floor(lb - margin, round_to_digits))

            if -margin <= ub <= 0:
                r.upper_bound = min(r.upper_bound, 0.0)
            else:
                r.upper_bound = min(r.upper_bound, ceil(ub + margin, round_to_digits))

class ThermodynamicSpace_mod:
    def __init__(
        problem,
        S_constraints: np.ndarray, ## uses rxn for which second law is considered but would like to use th whole model
        reaction_ids,
        metabolites,
        drg0_prime_mean,
        drg0_prime_cov_sqrt,
        log_conc_cov,
        log_conc_mean ):

        problem.S_constraints = S_constraints
        problem.reaction_ids = reaction_ids
        problem.metabolites = metabolites
        ##drg and lnc
        problem.drg0_prime_mean, problem.drg0_prime_cov_sqrt = drg0_prime_mean.reshape(-1, 1), drg0_prime_cov_sqrt
        problem.log_conc_cov=log_conc_cov
        problem.log_conc_mean=log_conc_mean.reshape(-1,1)

    # def _make_concentrations_distribution(
    #     self,
    #     metabolites:[],
    #     concentrations: pta.concentrations_prior.ConcentrationsPrior,
    # ) :
    #     """Construct the distribution over metabolite concentrations."""
    #     log_conc_mean = qvector(
    #         [
    #             Q(concentrations.get(m.miriam_id, m.compartment).log_mean)
    #             for m in metabolites
    #         ]
    #     )
    #     log_conc_cov = Q(
    #         np.diag(
    #             np.square(
    #                 qrvector(
    #                     [
    #                         Q(concentrations.get(m.miriam_id, m.compartment).log_std)
    #                         for m in metabolites
    #                     ]
    #                 ).m
    #             )
    #         )
    #     )

    #     return log_conc_mean, log_conc_cov
class ThermodynamicSpaceBasismod(pta.ThermodynamicSpaceBasis):
    '''new class to modify : 
    - pint unit 
    - add possibility to ignore concentrations of some metabolites in the drG = .. eq'''
    def __init__(
        self,
        thermodynamic_space: "ThermodynamicSpace",
        explicit_log_conc: bool = True,
        explicit_drg0: bool = True,
        explicit_drg: bool = True,
        min_eigenvalue: float = default_min_eigenvalue_tds_basis,
        ranktol: float =1e-5, ##rank tolerance for linalg qr 
        ignoreconc_met_idx: list = None,
    ):
        if ignoreconc_met_idx is None:
            ignoreconc_met_idx = []

        ##super.init calls compute basis, but made my own func so dont want that. just copied the rest of the init here.
        self._to_log_conc_transform = None
        self._to_drg0_transform = None
        self._to_drg_transform = None

        self._observables_ranges:dict = {
            "drg": [],
            "drg0": [],
            "log_conc": [],
        }
        ## mets  to ignore
        self.ignoreconc_met_idx = ignoreconc_met_idx

        self._compute_basis_mod(
            thermodynamic_space,
            explicit_log_conc,
            explicit_drg0,
            explicit_drg,
            min_eigenvalue,
            ranktol
        )

    def _compute_basis_mod(
        self,
        thermodynamic_space: "ThermodynamicSpace",
        explicit_log_conc,
        explicit_drg0,
        explicit_drg,
        tolerance,
        ranktol,T=303.15
    ):
        """Compute the basis of the thermodynamic space, only making explicit the
        specified variables.
        """
        log_conc_mean = thermodynamic_space.log_conc_mean#.m
        log_conc_cov = thermodynamic_space.log_conc_cov#.m
        drg0_prime_mean = thermodynamic_space.drg0_prime_mean#.m_as("kJ/mol")
        drg0_prime_cov_sqrt = thermodynamic_space.drg0_prime_cov_sqrt#.m_as("kJ/mol")
        S_con_T = thermodynamic_space.S_constraints.T
        ##icnore conc
        S_con_T_mod = S_con_T.copy()
        #actually cant just remove these rows
        S_con_T_mod[:, self.ignoreconc_met_idx] = 0
        ### we dont ignore conc by putting 0 but just not taking thcorrepsonding rows in S 
        S_con_T_mod=S_con_T[:,[idx for idx in range(S_con_T.shape[1]) if idx not in self.ignoreconc_met_idx]]

        RT = (R.m * T)#thermodynamic_space.parameters.T()).m_as("kJ/mol")
        n_mets = log_conc_mean.size
        n_rxns = drg0_prime_mean.size

        # Use the eigenvalue decomposition to find the minimal basis for the
        # log concentrations.
        print(log_conc_cov.shape)
        log_conc_cov_sqrt = covariance_square_root(log_conc_cov, tolerance)

        # Combined mean and squared root of the covariance of log-concentrations
        # and DrG'°.
        log_conc_drg0_mean = np.block([[log_conc_mean], [drg0_prime_mean]])

        log_conc_drg0_cov_sqrt = scipy.linalg.block_diag(
            log_conc_cov_sqrt, drg0_prime_cov_sqrt
        )

        # Construct the transform from log-concentrations and DrG'° to the
        # variables that must be explicitly represented ("observables").
        mapping_blocks = []
        if explicit_log_conc:
            mapping_blocks.append([np.identity(n_mets), np.zeros((n_mets, n_rxns))])
        if explicit_drg0:
            mapping_blocks.append([np.zeros((n_rxns, n_mets)), np.identity(n_rxns)])
        if explicit_drg:
            mapping_blocks.append([RT * S_con_T_mod, np.identity(n_rxns)])

        to_observables = np.block(mapping_blocks)
        #print shapes
        print(to_observables.shape, log_conc_drg0_mean.shape, log_conc_drg0_cov_sqrt.shape)
        observables_mean = to_observables @ log_conc_drg0_mean
        observables_cov_sqrt = to_observables @ log_conc_drg0_cov_sqrt

        # Find a square root for the covariance matrix with minimum
        # dimensionality using the QR decomposition.
        observables_cov_sqrt = LINALG.qr_rank_deficient(observables_cov_sqrt.T,ranktol).T ## use ranktol=1e-5 if we anyway clip them later
        observables_cov_sqrt[np.abs(observables_cov_sqrt) < 1e-5] = 0
        self._to_observables_transform = (observables_cov_sqrt, observables_mean)
        self._dimensionality = observables_cov_sqrt.shape[1]
        self._sigmas = la.norm(observables_cov_sqrt, axis=0).T
        ##add attributes for EP
        self.to_observables=to_observables
        self.log_conc_drg0_mean= log_conc_drg0_mean
        self.log_conc_drg0_cov_sqrt=log_conc_drg0_cov_sqrt

        # Finally, identify the transforms from the basis to the single
        # observables.
        current_id = 0
        if explicit_log_conc:
            self._observables_ranges["log_conc"] = list(
                range(current_id, current_id + n_mets)
            )
            self._to_log_conc_transform = (
                observables_cov_sqrt[self.observables_ranges["log_conc"], :],
                observables_mean[self.observables_ranges["log_conc"]],
            )
            current_id += n_mets
        if explicit_drg0:
            self._observables_ranges["drg0"] = list(
                range(current_id, current_id + n_rxns)
            )
            self._to_drg0_transform = (
                observables_cov_sqrt[self.observables_ranges["drg0"], :],
                observables_mean[self.observables_ranges["drg0"]],
            )
            current_id += n_rxns
        if explicit_drg:
            self._observables_ranges["drg"] = list(
                range(current_id, current_id + n_rxns)
            )
            self._to_drg_transform = (
                observables_cov_sqrt[self.observables_ranges["drg"], :],
                observables_mean[self.observables_ranges["drg"]],
            )


class SupTFSmodel(tfs.TFSModel):
    def __init__(self,model,thermodynamic_space,tbasis,problem,num_samples=5e5,num_dir_samples=5e7,num_initial_steps=None,max_steps=2**25,max_threads=32,warmup=2,**kwargs):
        super().__init__(model,thermodynamic_space,tbasis,**kwargs)
        self.model=model ## base class doesnt store model
        self.problem=problem
        self.num_samples=num_samples
        self.num_dir_samples=num_dir_samples
        self.max_steps=max_steps
        self.num_initial_steps=num_initial_steps
        self.max_threads= max_threads
        self.warmup=warmup
        #methdos
        self.find_rev_irr_dirs()
        # self.sampling_settings()
        self.directions_transform=None ##none for now but can modify 

    def find_rev_irr_dirs(self):
            # Find candidate optimization direactions.
        reaction_idxs_T = list(range(len(self.T.reaction_ids)))
        reaction_idxs_F = [self.F.reaction_ids.index(id) for id in self.T.reaction_ids]
        only_forward_ids_T = [
            i for i in reaction_idxs_T if self.F.lb[reaction_idxs_F[i]] >= 0
        ]
        only_backward_ids_T = [
            i for i in reaction_idxs_T if self.F.ub[reaction_idxs_F[i]] <= 0
        ]
        reversible_ids_T = [
            i
            for i in reaction_idxs_T
            if self.F.lb[reaction_idxs_F[i]] < 0 and self.F.ub[reaction_idxs_F[i]] > 0
        ]

        reversible_dirs = [(i, -1) for i in reversible_ids_T] + [
            (i, 1) for i in reversible_ids_T
        ]
        irreversible_dirs = [(i, -1) for i in only_backward_ids_T] + [
            (i, 1) for i in only_forward_ids_T
        ]

        self.reversible_dirs = reversible_dirs
        self.irreversible_dirs = irreversible_dirs

    
    def safe_find_point(self,direction):
        try:
            return tfs._find_point(self.problem,direction)
        except Exception as e:
            # log or print if you want to debug
            print(f"Failed on {direction}: {e}")
            return None   # or some default
        
    def sampling_settings(self,points_array):
        settings = tfs.pb.FreeEnergySamplerSettings()
        settings.truncation_multiplier = self.confidence_radius
        settings.feasibility_cache_size = tfs_default_feasibility_cache_size
        settings.drg_epsilon = self.drg_epsilon
        settings.flux_epsilon =  self.problem._epsilon_f[0][0] # ##nf use same than pmo pb constraints
        settings.min_rel_region_length = tfs_default_min_rel_region_length
        self.update_settings_function = lambda s, steps: tfs._fill_settings(
            s, self.num_samples, steps, points_array.shape[1], self.warmup, self.max_threads, self.num_dir_samples
        )
        print(f'settings : threads {self.max_threads}, num samples {self.num_samples}, max steps {self.max_steps},\
               num initial steps {self.num_initial_steps}, warmup {self.warmup}, num_dirsamples {self.num_dir_samples}')
        self.update_settings_function(settings, self.num_initial_steps)
        self.make_sampling_result_function = lambda r: tfs._make_sampling_result(
            self, r, self.num_samples
        )
        self.settings=settings

    def _conv_manager(self,directions_transform=None):
        if self.num_initial_steps == None :
            self.num_initial_steps = ( self.B.dimensionality**2 + self.num_samples)

        self._convergence_manager = cm.ConvergenceManager(self, self.num_initial_steps, True, directions_transform)

    def sample(self,points_array):
        if not type(points_array) == np.ndarray:
            points_array=np.hstack(points_array)
        self.sampling_settings(points_array)
        self._conv_manager(directions_transform=self.directions_transform)

        result = self._convergence_manager.run(
        self.settings,
        self.update_settings_function,
        self.make_sampling_result_function,
        points_array,
        self.max_steps,
        default_max_psrf,
        )
        maxpsrf=np.max(result.psrf)
        print(f"Sampling converged {maxpsrf<1.1} with max psrf {maxpsrf}")
        return result
    
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main(cobrafile,vbounds_file,
    logconcmean_file,logconcvar_file,lncbounds_file,
    drG0file,drGcovfile,ranktol=1e-5,sector_dict=None):
    '''
    Docstring for main
    
    :param cobrafile: modelfile
    :param vbounds_file: tfva bounds file
    :param logconcmean_file: lnc mean approx. from EP
    :param logconcvar_file: lnc var approx. from EP
    :param lncbounds_file: TVA lnc bounds file
    :param drG0file: drG0mean estimates from eQuilibrator
    :param drGcovfile: drG0cov
    :param ranktol: Description
    '''
    prep=preprocess(cobrafile,vbounds_file,logconcmean_file,logconcvar_file,lncbounds_file,drG0file,drGcovfile)

    if sector_dict is not None:
        hs.set_sector_flux(prep.model,sector_dict)


    thermodynamic_space = ThermodynamicSpace_mod(prep.Sconstrained,prep.rid_constrained,prep.model.metabolites,
                                            prep.drg0_prime_mean,prep.drg0_prime_cov_sqrt,prep.log_conc_cov,prep.log_conc_mean)
    
    #get basis mod
    tbasis=ThermodynamicSpaceBasismod(thermodynamic_space,False,False,True,min_eigenvalue=1e-10,ranktol=ranktol,ignoreconc_met_idx=prep.ignoreconc_met_idx)
    
    # # Create a PMO problem
    problem = pta.PmoProblem(prep.model, thermodynamic_space,tbasis,solver="GUROBI")

    tfs_model = SupTFSmodel(prep.model, thermodynamic_space,tbasis,problem,solver=cp.GUROBI)
    tfs_model.ignoreconc_met_idx=prep.ignoreconc_met_idx 

    return tfs_model


