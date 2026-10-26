import thermo_flux
from thermo_flux.io import load_excel as ex
from thermo_flux.core.model import ThermoModel
from thermo_flux.core.reaction import ThermoReaction
from cobra.core import Reaction
from equilibrator_api import  Q_
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from thermo_flux.io import helper_load as hl
from math import ceil

"""File to skip a few long cells in the notebook and build the yeast model
version = 1.0
date=2024-06-10
author = Nathan Fargier
"""
INPUT_FOLDER='/home/nathan/PhD/04.Computational/202601_13C/1_data/a_raw_data/'

### global variables ###
model_excel=f"{INPUT_FOLDER}yeast_v3.211_starting_curated_withprotons.xlsx"
# model_excel="/home/nathan/thermo_flux/Nathan/glucose/excel_models/yeast_2025model.xlsx"

conc_file=f"{INPUT_FOLDER}ConcLimits.xlsx"
########################
macromol_reactions=['lipids_reaction','dna_reaction','cell-wall_reaction','protein_reaction','rna_reaction','storage_reaction']

def build_yeast_model(model_excel :str =model_excel,cycle_stoich: dict = None,bal_mg=False,round_dp=3) -> ThermoModel:
    #load excel model from Vakil (required to recalculate drG0')

    emodel = ex.create_thermo_model('yeast',model_excel=model_excel,
                                    keggids_csv=f"{INPUT_FOLDER}yeast_kegg_id.csv",update_dfg_bio=False )#edit_mets=EDIT_METS)    

    # #create thermo model object 
    tmodel = emodel

    # #protons have a positive charge
    for h in tmodel.proton_dict.values():
        h.charge = 1
        h.formula = 'H1'
        print(h,h.formula,h.charge)

    tmodel.reactions.get_by_id("biomass_EX").bounds = (0,500)
    tmodel.reactions.get_by_id("biomass_c").bounds = (0,500)
    tmodel.reactions.get_by_id("biomass_ce").bounds = (0,500)


    #basic checks for thermo model - h_ex and poscharge_EX should not be in the second law
    for rxn in tmodel.boundary:
        rxn.ignore_snd = True

        
    #set lower bounds of all exchange reactions to 0 except lactate because we mimic microfluidic settings
    for r in ['pyr_EX','glyc_EX','ac_EX','succ_EX','etoh_EX','co2tot_EX','gal-D_EX','man-D_EX','malt_EX'] :
    # for r in ['pyr_EX','glyc_EX','ac_EX','succ_EX','etoh_EX','co2tot_EX'] :
        tmodel.reactions.get_by_id(r).lower_bound=0
        # tmodel.reactions.get_by_id(r).ignore_snd = False
        print(tmodel.reactions.get_by_id(r).reaction)


    #assign biomass
    tmodel.metabolites.biomass_c.biomass = True
    tmodel.metabolites.biomass_e.biomass = True


    #define specific redox carier energies and charge #formation energies are retrieved from the GAMS model 
    tmodel.metabolites.focytc_c.redox = True
    tmodel.metabolites.focytc_c.dfG0 = Q_(-12.25 , "kJ/mol")
    tmodel.metabolites.focytc_c.dfG0prime = Q_(-12.25 , "kJ/mol")
    tmodel.metabolites.focytc_c.dfG_SE = Q_(1 , "kJ/mol")
    tmodel.metabolites.ficytc_c.redox = True
    tmodel.metabolites.ficytc_c.dfG0 = Q_(12.25 , "kJ/mol")
    tmodel.metabolites.ficytc_c.dfG0prime = Q_(12.25 , "kJ/mol") 
    tmodel.metabolites.ficytc_c.dfG_SE = Q_(1 , "kJ/mol")
    tmodel.metabolites.focytc_c.charge = 2
    tmodel.metabolites.ficytc_c.charge = 3


    #in imported model any charge metabolite represents a free cation that is transported 
    #this must be defined in the transported_charge field which represents the free charges (ions)
    #that are trnasported in a transport reaction
    #this should be done in model initialisation 
    for rxn in tmodel.reactions:
        transported_charge = {}
        transported_h={}
        for met, stoich in rxn.metabolites.items():
            if (met in tmodel.charge_dict.values()):
                transported_charge[met.compartment]= stoich
            if (met in tmodel.proton_dict.values()):
                transported_h[met.compartment]= stoich
        if len(transported_charge) != 0:
            print(rxn, transported_charge)
            rxn.transported_charge = transported_charge

    for met in tmodel.mg_dict.values():
        met.ignore_conc=True
        print(met,met.ignore_conc)## we transform dfG0 to the constant pMg value so ignore the concentration

    # # # # # ## pMg at 14 because we don't consider it
    if bal_mg is False:
        for comp in tmodel.compartments:
            tmodel.pMg[comp]=Q_(14,'dimensionless')

    else : 
        tmodel.pMg={'m':Q_(3.2),'c':Q_(2.9),'e':Q_(3)}


    #should not be run after having balanced the reactions : it will add the balanced charges to  transported_charge 

    # #ensure transported h are proterly defined for ComplexII 
    tmodel.reactions.ComplexII.transported_h = {'m': -2.0, 'c': 2.0}
    #the change in charge between the 2 compartment must be defined
    tmodel.reactions.ComplexII.transported_charge = {'m': 2.0, 'c': -2.0}

    ## then we balance the reactions and add protons in macromol and biomass
    for rxn in tmodel.reactions:
        thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=bal_mg,round_dp=round_dp)

    ##cell cycle part 
    if cycle_stoich is not None:
        print(tmodel.reactions.biomass_c.reaction)
        forster_massfraction={'proteins_c' : 0.45, 'lipids_c' : 0.029, 'polysacch_c' : 0.315, 'dna_c' : 0.004, 'rna_c' : 0.063,'storage_c' : 0.092} 
        for macromol,stoich in cycle_stoich.items():
            print(macromol,stoich)
            current_stoich=tmodel.reactions.biomass_c.metabolites[tmodel.metabolites.get_by_id(macromol)]
            tmodel.reactions.biomass_c.add_metabolites({tmodel.metabolites.get_by_id(macromol):-1*current_stoich})
            tmodel.reactions.biomass_c.add_metabolites({tmodel.metabolites.get_by_id(macromol):stoich*current_stoich/forster_massfraction[macromol]})

        print(tmodel.reactions.biomass_c.reaction)

    ###set the correct numbers of protons coming from atp hydrolysis costs of macromols polymerisation - no because we dont touch the macromol reaciton but only biomass
    # cycle_stoich={'proteins_c': 1, 'lipids_c': 1, 'polysacch_c': 1, 'dna_c': 1, 'rna_c': 1, 'storage_c': 1} if cycle_stoich is None else cycle_stoich
    print(cycle_stoich)
    ##if stoichiometry varies, protons from atp hydrolysis will vary as well

    ##adjust charge and protons in biomass and macromol rxn
    balance_mac_bio(tmodel)

    ##update dfg
    update_dfg_bio(tmodel)


    ##flux and objectives:
    objective_fluxbounds(tmodel)

    # ##metabolite bounds
    metabolite_bounds(tmodel,conc_file)

    return tmodel

def balance_mac_bio(tmodel,round_dp=3,hyd_prot=0.7403927813163482): ##hyd_prot =0.74 at pH 7
    ##these eractions need to be balanced first (protons added in mets)
    tmodel.metabolites.proteins_c.elements={'H':tmodel.reactions.protein_reaction.metabolites[tmodel.metabolites.h_c]+
                                                  hyd_prot*abs(tmodel.reactions.protein_reaction.metabolites[tmodel.metabolites.atp_c]) }
    tmodel.metabolites.lipids_c.elements={'H':tmodel.reactions.lipids_reaction.metabolites[tmodel.metabolites.h_c]+
                                             hyd_prot*abs(tmodel.reactions.lipids_reaction.metabolites[tmodel.metabolites.atp_c]) }
    tmodel.metabolites.dna_c.elements={'H':tmodel.reactions.dna_reaction.metabolites[tmodel.metabolites.h_c]+
                                             hyd_prot*abs(tmodel.reactions.dna_reaction.metabolites[tmodel.metabolites.atp_c]) }
    tmodel.metabolites.rna_c.elements={'H':tmodel.reactions.rna_reaction.metabolites[tmodel.metabolites.h_c]+
                                             hyd_prot*abs(tmodel.reactions.rna_reaction.metabolites[tmodel.metabolites.atp_c]) }
    tmodel.metabolites.polysacch_c.elements={'H':tmodel.reactions.get_by_id('cell-wall_reaction').metabolites[tmodel.metabolites.h_c]+
                                                  hyd_prot*abs(tmodel.reactions.get_by_id('cell-wall_reaction').metabolites[tmodel.metabolites.atp_c])}
    tmodel.metabolites.storage_c.elements={'H':tmodel.reactions.storage_reaction.metabolites[tmodel.metabolites.h_c]+
                                             hyd_prot*abs(tmodel.reactions.storage_reaction.metabolites[tmodel.metabolites.atp_c]) }
    #then balance reactions again
    for rxn in [tmodel.reactions.protein_reaction,tmodel.reactions.lipids_reaction,tmodel.reactions.get_by_id('cell-wall_reaction'),tmodel.reactions.dna_reaction,tmodel.reactions.rna_reaction,tmodel.reactions.storage_reaction]:
        thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False,round_dp=round_dp)

        ### compute the charge in each macromol
    for r in ['lipids_reaction','dna_reaction','cell-wall_reaction','protein_reaction','rna_reaction','storage_reaction']:#calculate protons from precursors
        if r=='cell-wall_reaction':
            name='polysacch_c'
        elif r=='protein_reaction':
            name='proteins_c'
        else:
            name=r[:-9]+'_c'
        ##all reactions that compose the bioamss split eq
        #take the macromol and look at its elements
        tmodel.metabolites.get_by_id(name).charge=-1*np.round(tmodel.reactions.get_by_id(r).net_elements()[0][tmodel.metabolites.charge_c],round_dp)
    
    # # #define the protons in biomass 
    tmodel.metabolites.get_by_id('biomass_c').elements={'H':67} ## canelas exp value paper in vivo framework ...
    tmodel.metabolites.get_by_id('biomass_e').elements={'H':67}
    ##balance to consider charge from released protons
    for rxn in [tmodel.reactions.biomass_c,tmodel.reactions.biomass_ce,tmodel.reactions.biomass_EX]:
        thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False,round_dp=round_dp)
    ##charge in bio
    tmodel.metabolites.biomass_c.charge = -1*np.round(tmodel.reactions.biomass_c.net_elements()[0][tmodel.metabolites.charge_c],round_dp)
    tmodel.metabolites.biomass_e.charge =tmodel.metabolites.biomass_c.charge 
    #then balance again biomass rxn
    for rxn in [tmodel.reactions.biomass_c,tmodel.reactions.biomass_ce,tmodel.reactions.biomass_EX]:
        thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False,round_dp=round_dp)

def update_dfg_bio(tmodel:ThermoModel,round_dp=3):
    tmodel.update_thermo_info(fit_unknown_dfG0=False,round_dp=round_dp)#
    #is in fact in j/gDW : initial untransformed value taken from Battley dfGbm of -80.27kJ/cmol *1/yeastdensity of 26.4cmol/gdW
    tmodel.metabolites.biomass_c.dfG0=Q_(-3040,'kJ/mol') # could also be dfG0 
    tmodel.metabolites.biomass_e.dfG0=Q_(-3040,'kJ/mol')#tmodel.metabolites.biomass_c.dfG0prime#-tmodel.reactions.biomass_ce.drGtransport 

    tmodel.update_thermo_info(fit_unknown_dfG0=False,round_dp=round_dp)#,improved_generative_model=np.sqrt(6.5),
    print(tmodel.reactions.biomass_c.drG0prime,tmodel.reactions.biomass_EX.drG0prime)
    print(tmodel.metabolites.biomass_c.dfG0prime,tmodel.reactions.biomass_ce.drGtransport+tmodel.reactions.biomass_ce.drG0prime)


def objective_fluxbounds(tmodel:ThermoModel) -> None :
    #set flux bounds to match GAMS

    tmodel.objective = tmodel.reactions.biomass_EX

    for rxn in tmodel.reactions:
        if rxn.lower_bound == -1000:
            rxn.lower_bound = -500
            
        if rxn.upper_bound == 1000:
            rxn.upper_bound = 500


    tmodel.reactions.get_by_id("etoh_EX").upper_bound = 500
    tmodel.reactions.get_by_id("etoh_EX").lower_bound = 0
    tmodel.reactions.get_by_id("glyc_EX").upper_bound = 500
    tmodel.reactions.get_by_id("glyc_EX").lower_bound = 0
    tmodel.reactions.get_by_id("ac_EX").upper_bound = 500
    tmodel.reactions.get_by_id("ac_EX").lower_bound = 0
    tmodel.reactions.get_by_id("pyr_EX").upper_bound = 500
    tmodel.reactions.get_by_id("pyr_EX").lower_bound = 0
    tmodel.reactions.get_by_id("succ_EX").upper_bound = 500
    tmodel.reactions.get_by_id("succ_EX").lower_bound = 0
    tmodel.reactions.get_by_id("glc-D_EX").upper_bound = 0
    tmodel.reactions.get_by_id("glc-D_EX").lower_bound = -500

def metabolite_bounds(tmodel: ThermoModel,conc_file:str, list_comp=['c','e','m']) -> None:

    # use the same met bounds as in gams model
    ConcLimits_df = pd.read_excel(conc_file, 
                    header=[0,1], 
                    index_col=[0], 
                    )

    for met in [m for m in tmodel.metabolites if m.compartment in list_comp]:
        if met.id[:-2] not in  ['charge','biomass', 'cation','Mg'] and met.id[:-2] in ConcLimits_df.index:   
            comp = met.compartment
            ub = ConcLimits_df.loc[met.id[:-2]][comp]['ub']
            lb = ConcLimits_df.loc[met.id[:-2]][comp]['lb']

            met.upper_bound = Q_(ub, 'mM')
            met.lower_bound = Q_(lb, 'mM')

    ## adjust pi_c concentration bounds


# file ='/home/nathan/thermo_flux/Nathan/glucose/202406METB/202406METB_regression/202406METB_canelascond_regression_modelv211.sol'

def add_vars(tmodel,sigmac_limit=130,file =None,gdiss_constraint=True): ## function to add one condition gurobi model
    tmodel.m = None #clear any previously build optimization models 
    tmodel.add_TFBA_variables(conds=['one'],error_type = 'covariance', qnorm ='sep_norm' , alpha = 0.95, 
                         gdiss_constraint = gdiss_constraint, sigmac_limit = sigmac_limit,obj2=None,
                        weight_growth=1,weight_obj2=-1,enthalpy=False
                        )#change entropy limit according to regression result

    #regression free direction and metabolite bounds of only-c-mets and e mets
    
    # thermo_flux.solver.gurobi.model_start(tmodel,file,ignore_vars=['all'],fix_vars=['qm'],fix='bound')###testt

    tmodel.m.update()

# def adjust_forstercoeff(tmodel):
#     macromol_rxn=['lipids_reaction','dna_reaction','cell-wall_reaction','protein_reaction','rna_reaction']
#     macromols=['lipids_c','dna_c','polysacch_c','proteins_c','rna_c','storage_c']
#     forster_massfraction={'protein_reaction' : 0.45, 'lipids_reaction' : 0.029, 'cell-wall_reaction' : 0.315, 'dna_reaction' : 0.004, 'rna_reaction' : 0.063} 
#     for rid in macromol_rxn:
#         model_r=tmodel.reactions.get_by_id(rid)
#         for met,stoich in [(m,s) for (m,s) in model_r.metabolites.items() if m.id not in macromols ]: #substrates
#             model_r.add_metabolites({met:-1*stoich}) # first remove then add the adjusted coeff
#             model_r.add_metabolites({met:stoich/forster_massfraction[rid]})
