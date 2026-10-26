"""Preprocessing utilities for preparing a COBRA model and thermodynamic inputs.

Contains:
  - Preprocess: class that loads model and external files (flux bounds,
    log-concentration priors, drG0 means/covariances), removes blocked
    reactions and unused metabolites, and produces constrained stoichiometric
    matrices and arrays used by downstream components.

The implementation mirrors the logic extracted from the original helper_pta.py
while being organized for testing and reuse.
"""

import os
import pandas as pd
import numpy as np
from cobra.flux_analysis.variability import find_blocked_reactions
from cobra.util.array import create_stoichiometric_matrix
from cobra.io import read_sbml_model

import math
from scipy.stats import norm
# from sampler import helper_sampling as hs

class Preprocess:
    """Loads files, sanitizes a COBRA model and provides prepared arrays.

    Parameters
    ----------
    cobra_file : str
        Path to SBML model file.
    vbound_file : str
        CSV with flux bounds indexed by reaction id, columns ['lb','ub'].
    logconcmean_file : str
        CSV with log-concentration mean estimates (EP output).
    logconcvar_file : str
        CSV with log-concentration variances (EP output).
    lncbounds_file : str
        CSV of TVA bounds for log-concentrations (used to detect fixed concs).
    drG0file : str
        CSV containing DrG'° means (equilibrator).
    drGcovsqrtfile : str
        CSV containing square-root of DrG'° covariance matrix.
    drG0covfile : str, optional
        Full covariance (not square-root) used optionally by EP-based parts.
    """

    def __init__(
        self,
        cobra_file,
        rign_file,
        vbound_file,
        met_data_file,condition, ##cond required for met file
        lncbounds_file,
        drG0file,
        drGcovsqrtfile,
        drG0covfile=None,
        sector_dict=None,
    ):
        data_bounds_flux = pd.read_csv(vbound_file, index_col=0)
        self.data_bounds_flux=data_bounds_flux
        self.model = read_sbml_model(cobra_file)

        # Load drG0 info from provided files
        drg0_prime_mean_init = pd.read_csv(drG0file, index_col=0).values.T[0]
        drg0_prime_cov_sqrt_init = pd.read_csv(
            drGcovsqrtfile, index_col=0
        ).values
        if drG0covfile is not None:
            drg0_prime_cov_init = pd.read_csv(drG0covfile, index_col=0).values
        else : 
            drg0_prime_cov_init = drg0_prime_cov_sqrt_init@drg0_prime_cov_sqrt_init.T

        # Record initial metabolites and remove specific ions (Mg)
        initial_mets_id = [m.id for m in self.model.metabolites]
        # remove Mg if present (safe)
        for mg_id in ("Mg_c", "Mg_e", "Mg_m"):
            if mg_id in [m.id for m in self.model.metabolites]:
                try:
                    self.model.metabolites.remove(self.model.metabolites.get_by_id(mg_id))
                except Exception:
                    pass

        # set flux bounds to tfva bounds
        self.add_v_margin(data_bounds_flux)

        ##aply sector
        if sector_dict is not None:
            hs.set_sector_flux(self.model, sector_dict)

        # Find blocked reactions and remove them (remember original indices)
        init_reactions_id = [r.id for r in self.model.reactions]
        rid_removed = find_blocked_reactions(self.model)
        print('removing',rid_removed)
        rxn_removed = {
            old_idx: rid for old_idx, rid in enumerate(init_reactions_id) if rid in rid_removed
        }
        ridx_notremoved = [
            oldidx for oldidx in range(len(init_reactions_id)) if oldidx not in rxn_removed.keys()
        ]
        self.ridx_not_removed=ridx_notremoved
        self.n_r_removed=len(rid_removed)

        # Remove blocked reactions and orphan metabolites
        self.model.remove_reactions(rid_removed, remove_orphans=True)

        # Stoichiometric matrix for the cleaned model
        S = create_stoichiometric_matrix(self.model)
        self.S = S

        # Reactions for which the second law is ignored (made a series with rxn id : flag for ignore_snd + Ht,Htm,POSQt,POSQtm
        ##why proton and charge transport ? because shoudnt be a drGerror or drG onc  component in determining their drG :entiely determined by direction based on drG0
        ser_rxn_ign=pd.read_csv(rign_file,index_col=0)['0']
        snd_ignored=ser_rxn_ign[ser_rxn_ign==True].index

        snd_ignored_newidx = [self.model.reactions.index(rid) for rid in snd_ignored.intersection([r.id for r in self.model.reactions])]
        snd_ignored_oldidx=[init_reactions_id.index(rid) for rid in snd_ignored]

        self.snd_ignored_newidx=snd_ignored_newidx

        # Constrained reaction indices and mapping
        self.Sconstrained = S[:, [idx for idx in range(len(self.model.reactions)) if idx not in snd_ignored_newidx]]
        self.rid_constrained = [self.model.reactions[idx].id for idx in range(len(self.model.reactions)) if idx not in snd_ignored_newidx]

        ## save mets and rxn in index,col of S_df
        self.S_df_constrained = pd.DataFrame(self.Sconstrained,columns=self.rid_constrained,index=[m.id for m in self.model.metabolites])

        # # ---- Load concentration priors ----
        # Detect metabolites with fixed concentration (lb == ub in TVA bounds)
        bounds_lnconc = pd.read_csv(lncbounds_file, index_col=0).loc[[m.id for m in self.model.metabolites]]
        self.bounds_lnconc=bounds_lnconc
        ignore_met_id = bounds_lnconc[bounds_lnconc["ub"] - bounds_lnconc["lb"] <1e-6].index
        self.ignore_met_id=ignore_met_id
        ignoreconc_met_idx = [self.model.metabolites.index(met_id) for met_id in ignore_met_id]
        self.ignoreconc_met_idx = ignoreconc_met_idx
        #### build log conc dist using metabolome data if exists otherwise default as in pta example
        ###default is LogNormalDistribution(-8.3371, 1.9885)
        S_df_constrained_redmets=self.S_df_constrained.drop(ignore_met_id)
        self.S_df_constrained_redmets=S_df_constrained_redmets

        ## from the list of metabolites, construct two matrices, mean and cov (mxm) with the either the expeirmental dist or default valeus
        met_data_cond=pd.read_csv(met_data_file,index_col=(0,1)).loc[condition]
        log_conc_mean=np.zeros(len(S_df_constrained_redmets.index))
        log_conc_cov=np.eye(len(S_df_constrained_redmets.index))* 1.9885**2
        for idx,mid in enumerate(S_df_constrained_redmets.index):
            mname=mid.split('_')[0]
            if mid  in ['nad_m','nadph_m']:
                print(f'{mid} ignored because we set constraint on redox ratio')
                continue
            try:
                mean=met_data_cond.loc[mname,'mean']
                log_conc_mean[idx]=np.log(mean*1e-3) ##the unit in pta is M
                log_conc_cov[idx,idx]=0.6**2 ## in PTA example, they use 0.2
            except Exception as e:
                print('no met data for',mid,e)
                log_conc_mean[idx]=-8.3371

        self.log_conc_mean=log_conc_mean
        self.log_conc_cov=log_conc_cov
        ##usig EP files
        # log_conc_mean = pd.read_csv(logconcmean_file, index_col=0)["0"].values
        # log_conc_var = pd.read_csv(logconcvar_file, index_col=0)["0"].values
        # log_conc_cov = np.diag(log_conc_var)
        ##using lnc distribution from between bounds 

        # ---- Load drG'° statistics restricted to constrained reactions ----
        idx_constrained_old = [
            idx for idx in range(len(drg0_prime_mean_init)) if idx not in snd_ignored_oldidx and idx in ridx_notremoved  ##old indexes
        ]
    
        self.idx_constrained_old=idx_constrained_old
        self.drg0_prime_mean = drg0_prime_mean_init[idx_constrained_old]
        self.drg0_prime_cov_sqrt = drg0_prime_cov_sqrt_init[idx_constrained_old, :]#[:, idx_constrained_old]
        # if drG0covfile is not None:
        self.drG0_prime_cov = drg0_prime_cov_init[idx_constrained_old, :][:, idx_constrained_old]

    def add_v_margin(self, data_bounds_flux, margin: float = 1e-4, round_to_digits: int = 6):
        """Apply small margins to bounds read from a CSV to avoid exact-tight bounds.
        """
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

def floor(value: float, decimals: int = 0) -> float:
    """Same as :code:`math.floor`, except that it rounds to a given number of decimals.

    Parameters
    ----------
    value : float
        Value to round.
    decimals : int, optional
        Number of decimals to round to, by default 0.

    Returns
    -------
    float
        The rounded value.
    """
    return math.floor(value * (10**decimals)) / (10**decimals)


def ceil(value: float, decimals: int = 0):
    """Same as :code:`math.ceil`, except that it rounds to a given number of decimals.

    Parameters
    ----------
    value : float
        Value to round.
    decimals : int, optional
        Number of decimals to round to, by default 0.

    Returns
    -------
    float
        The rounded value.
    """
    return math.ceil(value * (10**decimals)) / (10**decimals)

def lognormal_from_ln_bounds(lnL_mM, lnU_mM, coverage=0.99, overshoot=1):
    """
    Compute mu and sigma for a lognormal distribution given ln concentrations in mM,
    convert to M, and optionally expand bounds to avoid exact-bound issues.

    lnL_mM, lnU_mM: bounds in natural log (ln) space
    coverage: probability coverage of the interval
    overshoot: factor >1 to expand the interval proportionally
    """
    # Convert ln(mM) -> ln(M): ln(M) = ln(mM * 1e-3) = ln(mM) + ln(1e-3)
    lnL_M = lnL_mM + np.log(1e-3)
    lnU_M = lnU_mM + np.log(1e-3)

    # Expand bounds linearly in ln-space
    center = (lnL_M + lnU_M) / 2
    half_width = (lnU_M - lnL_M) / 2
    lnL_expanded = center - overshoot * half_width
    lnU_expanded = center + overshoot * half_width

    # z-score for the given coverage
    z = norm.ppf(0.5 + coverage/2)

    # Parameters of underlying normal
    mu = (lnL_expanded + lnU_expanded) / 2
    sigma = (lnU_expanded - lnL_expanded) / (2 * z)

    return mu, sigma