import pandas as pd
from cobra.util.array import create_stoichiometric_matrix
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from equilibrator_api import R
from equilibrator_api import ComponentContribution, Q_



def open_solfile2(tmodel, file,vars=[] ,get_obj=False,nconds=1):
    '''open a solution file and return a pandas Multiindex dataframe
    used for getting regression results + Objective value
    vars has to be a list of strings with the variable names as they appear in the .sol file
    doesn't work for all types of variables (only the vars which are specificaly related to metabolites or reactions)
    output: (dataframe, obj.value)
    '''

#First define both dataframes and read the solution file
    sol= pd.DataFrame(columns = ['xn', 'names'])
    obj_val=None
    with open(file, 'r') as f:  
        for line in f:
            line.strip()
            if get_obj==True and 'Objective value' in line :
                obj_val=float(line.split()[-1])
                #happens only once so skip this if already done
                get_obj=False
            if len(line.split()) ==2:
                (names, xn) = line.split()
                split=[{'xn':float(xn), 'names':names}]
                df_split=pd.DataFrame.from_records(split)
                sol=pd.concat([sol,df_split])         

    sol.set_index('names', inplace=True)
    # df_sol=pd.DataFrame(index=vars)
    for var in vars: 
        df=sol[sol.index.str.startswith(var)].copy()
        if vars==['qm']:
            df.reset_index(inplace=True)
        else:
            if len(df.index)%len(tmodel.reactions)==0:
                multiindex=pd.MultiIndex.from_product([['d'+str(i) for i in range(1,nconds+1)],[rxn.id for rxn in tmodel.reactions]])
                df.index=multiindex
            if len(df.index)%len(tmodel.metabolites)==0:
                multiindex=pd.MultiIndex.from_product([['d'+str(i) for i in range(1,nconds+1)],[met.id for met in tmodel.metabolites]])
                df.index=multiindex
        # multiindex = pd.MultiIndex.from_product([df_sol.index, df.index])
        # df_sol=pd.concat([df_sol,df])
    
    if obj_val is not None:
        return df, obj_val
    else:
        return df
    

#analyse which constraint is violated and by how much 
def constraint_check(gurm,constrname: str)->None:
    """Check the constraint in the model and prints several attributes of the constraint and its variables"""
    c = gurm.getConstrByName(constrname)
    row = gurm.getRow(c)
    print(row, c.sense, c.rhs, '\n')

    total = 0 
    for var_idx in range(row.size()):
        var = row.getVar(var_idx)
        coeff = row.getCoeff(var_idx)
        start = var.start
        try:
            sol=var.x
        except:
            sol = 0
        lb = var.LB
        ub = var.UB
        print("coeff=",coeff,"var=",var,"start=",start,"lb=",lb,"ub=",ub,"coeff*start=",coeff*start,"sol=",sol)
        total+=coeff*start
    print('\n',total)


def check_metabolite_bounds(tmodel,metabolite_id:str) :
    """Check the bounds of a metabolite in the model in compartments c and m and return the results in a dictionary.
    model needs to be solved
    if the model is a multiscenario model, return compartmentalized conc and cellconc for each scenario in a dataframe, along with dict of bounds"""
    # First check in which compartments the metabolite is 
    comps = [met.compartment for met in tmodel.metabolites if metabolite_id == met.id.split('_')[0] if met.compartment in ['c','m']]
    # Create a dictionary to store the results
    results = {}
    cell_conc=0
    volume={'c':0.925,'m':0.075}
    # Get the index of the metabolite in the model
    if tmodel.m.numScenarios==0:
        for comp in comps:
            idx = tmodel.metabolites.index(tmodel.metabolites.get_by_id(metabolite_id + '_' + comp))
            
            # Get the lower and upper bounds of the metabolite in the compartment
            lb = np.exp(tmodel.mvars['ln_conc'][0][idx].lb)
            ub = np.exp(tmodel.mvars['ln_conc'][0][idx].ub)
            
            # Get the concentration value of the metabolite in the compartment
            conc = np.exp(tmodel.mvars['ln_conc'][0][idx].x)
            
            # Store the results in the dictionary
            results[comp] = {'Bounds': [lb, ub], 'Concentration': conc}
            cell_conc+=conc*volume[comp]
        results['cell_conc']=cell_conc
        # Return the dictionary of results
        return results
    
    else:## multiscenario
        numscen=tmodel.m.numScenarios
        df_conc=pd.DataFrame(index=['c','m','cellconc'],columns=[i for i in range(0,numscen)])
        #first make dict of bounds : 
        bounds={}
        for comp in comps:
            idx = tmodel.metabolites.index(tmodel.metabolites.get_by_id(metabolite_id + '_' + comp))
            lb = np.exp(tmodel.mvars['ln_conc'][0][idx].lb)
            ub = np.exp(tmodel.mvars['ln_conc'][0][idx].ub)
            bounds[comp]=[lb,ub]
        #then get the concs for each scenario
        for scen in range(numscen):
            tmodel.m.params.ScenarioNumber = scen
            tmodel.m.update()
            cell_conc=0
            for comp in comps:
                idx = tmodel.metabolites.index(tmodel.metabolites.get_by_id(metabolite_id + '_' + comp))
                conc = np.exp(tmodel.mvars['ln_conc'][0][idx].ScenNX) #scenario specific value of the var
                df_conc.loc[comp,scen]=conc
                cell_conc+=conc*volume[comp]
            df_conc.loc['cellconc',scen]=cell_conc

        return df_conc,bounds





from cobra.util.array import create_stoichiometric_matrix
def flux_of_metabolite(tmodel, fluxes_df:pd.DataFrame, MET:str)->pd.DataFrame:
    """given the stoicheometry matrix, the flux dataframe for all reactions and the metabolite we look for 
        we get the values of fluxes that produce or consume this mets, 
        returns a dataframe with the fluxes values of the reactions that produce or consume the metabolite"""
    S=create_stoichiometric_matrix(tmodel)
    S_T_df = pd.DataFrame(S.T, columns=[m.id for m in tmodel.metabolites], index=[r.id for r in tmodel.reactions])
    ##MET_rxns_stoich is a series of reactions with how much of the metabolite is produced or consumed
    MET_rxns_stoich=S_T_df[MET]
    MET_rxns_stoich=MET_rxns_stoich[MET_rxns_stoich!=0] ##drop 0 values : rxn that dont produce or consume the metabolite
    #multiply fluxes values by the stoicheometry of the metabolite in the reaction
    df_met_flux=fluxes_df.loc[MET_rxns_stoich.index]*MET_rxns_stoich
    return df_met_flux ##get the fluxes values of the reactions that produce or consume the metabolite

### i will iterate over the scenarios that are already made , and fix bounds for fluxes according to current values
# i will also remove the bounds of each of the 10 best mets
##then iterate on mets :
#to change the objective of each scenario to maximize the concentration of 1 metabolite
# put it in a dataframe and minimize the min, 
##all of the above will be one function, which i will call in a loop over the best mets and plot the results of this function
#then change the met and repeat min max
# met_names=df_best_mets_plot.index.to_list()

def VA_met_getcellconc2(tmodel,met_name_to_VA : str,fluxes : pd.DataFrame,directions :pd.DataFrame, conctype=['cellconc'],reset_met_bounds=False,fix_directions=False):
    """for the gurobi model, minimize and maximize concentration of a metabolite while keeping the fluxes constant as in the passed dataframe
    can return cellconc or compartmentalized conc"""
    numscen=tmodel.m.numScenarios
    ##initialize dataframe to store cell conc max and min across scenario
    df_cellconc=pd.DataFrame(index=[1,-1],columns=[i for i in range(numscen)])
    for opt in [1,-1]: ##1 for maximize, -1 for minimize
        ## first can remove bounds of best mets, no matter the scenario
        ## we reset bounds to the initial met bounds # for case 2  ?
        if reset_met_bounds:
            for comp in ['c','m']:
                idx=tmodel.metabolites.index(tmodel.metabolites.get_by_id(met_name_to_VA+'_'+comp))
                # tmodel.mvars['ln_conc'][0][idx].lb=np.log(min(tmodel.metabolites[idx].lower_bound.m,data_bounds_df.loc[met_name_to_VA+f'_{comp}','lb'])) # in mM 
                # tmodel.mvars['ln_conc'][0][idx].ub=np.log(max(tmodel.metabolites[idx].upper_bound.m,data_bounds_df.loc[met_name_to_VA+f'_{comp}','ub']))
                tmodel.mvars['ln_conc'][0][idx].lb=np.log(1e-4)
                tmodel.mvars['ln_conc'][0][idx].ub=np.log(1e3)
                #otherwise its just default metabolite bounds
                tmodel.m.update()
                print(tmodel.metabolites.get_by_id(met_name_to_VA+'_'+comp).id,comp, np.exp(tmodel.mvars['ln_conc'][0][idx].lb),np.exp(tmodel.mvars['ln_conc'][0][idx].ub))
        ##fix fluxes, on each current scenario flux values
        for scen in range(numscen):
            tmodel.m.params.ScenarioNumber=scen
            tmodel.m.update()

            for rxn_id,row in fluxes.iterrows():## we iterate over previously saved flux df
                rxn=tmodel.reactions.get_by_id(rxn_id)
                rxn_index=tmodel.reactions.index(rxn)
                tmodel.mvars['v'][0][rxn_index].ScenNLB=row[scen]#-1e-5
                tmodel.mvars['v'][0][rxn_index].ScenNUB=row[scen]#+1e-5

            if fix_directions:
                for rxn_id,row in directions.iterrows(): #we also fix b variable
                    rxn=tmodel.reactions.get_by_id(rxn_id)
                    rxn_index=tmodel.reactions.index(rxn)
                    tmodel.mvars['b'][0][rxn_index].ScenNLB=row[scen]
                    tmodel.mvars['b'][0][rxn_index].ScenNUB=row[scen]

        ##last is to change the objective, valid across scenarios
        #but will also have to iterate over met_names
        opt_dir=GRB.MAXIMIZE if opt==1 else GRB.MINIMIZE
        if len(conctype)==1:
            idx=tmodel.metabolites.index(tmodel.metabolites.get_by_id(met_name_to_VA+'_'+conctype))
            tmodel.m.setObjective(tmodel.mvars['ln_conc'][0][idx],opt_dir)
        else :
            idx_c=tmodel.metabolites.index(tmodel.metabolites.get_by_id(met_name_to_VA+'_c'))
            idx_m=tmodel.metabolites.index(tmodel.metabolites.get_by_id(met_name_to_VA+'_m'))
            tmodel.m.setObjective(tmodel.mvars['ln_conc'][0][idx_c]+tmodel.mvars['ln_conc'][0][idx_m],opt_dir)

        tmodel.m.params.OutputFlag=0
        tmodel.m.optimize()
        
        ##get the results with the modified function check_metabolite_bounds
        try :
            df_cellconc.loc[opt]=check_metabolite_bounds(tmodel,met_name_to_VA)[0].loc[conctype].values
        except Exception as e :
            print(e)

    return df_cellconc
# df_cellconc=VA_met_getcellconc(tmodel,met_names,'akg')

from gurobipy import GRB
### i will iterate over the scenarios that are already made , and fix bounds for fluxes according to current values
# i will also remove the bounds of each of the 10 best mets
##then iterate on mets :
#to change the objective of each scenario to maximize the concentration of 1 metabolite
# put it in a dataframe and minimize the min, 
##all of the above will be one function, which i will call in a loop over the best mets and plot the results of this function
#then change the met and repeat min max
# met_names=df_best_mets_plot.index.to_list()

def VA_fluxes_multiscenario(tmodel,flux_name_toVA : str,fluxes : pd.DataFrame,epsilon=1e-4,directions=None ,fix_directions=False):
    """for the gurobi model, minimize and maximize flux while keeping the fluxes constant as in the passed dataframe
    returns flux bounds
    use this function when a multiscenario model has already been built (after a gurscan for example)"""
    numscen=tmodel.m.numScenarios
    ##initialize dataframe to store cell conc max and min across scenario
    df_fluxes=pd.DataFrame(index=[1,-1],columns=[i for i in range(numscen)])
    for opt in [1,-1]: ##1 for maximize, -1 for minimize
        ## first can remove bounds of best mets, no matter the scenario
        ## we reset bounds to the initial met bounds # for case 2  ?
        ##fix fluxes, on each current scenario flux values
        for scen in range(numscen):
            tmodel.m.params.ScenarioNumber=scen
            tmodel.m.update()

            for rxn_id,row in fluxes.iterrows():## we iterate over previously saved flux df
                rxn=tmodel.reactions.get_by_id(rxn_id)
                rxn_index=tmodel.reactions.index(rxn)
                if rxn_id==flux_name_toVA:
                    tmodel.mvars['v'][0][rxn_index].ScenNLB=-GRB.UNDEFINED
                    tmodel.mvars['v'][0][rxn_index].ScenNUB=GRB.UNDEFINED
                else:
                    tmodel.mvars['v'][0][rxn_index].ScenNLB=row[scen]-epsilon
                    tmodel.mvars['v'][0][rxn_index].ScenNUB=row[scen]+epsilon
                

            if fix_directions:
                for rxn_id,row in directions.iterrows(): #we also fix b variable
                    rxn=tmodel.reactions.get_by_id(rxn_id)
                    rxn_index=tmodel.reactions.index(rxn)
                    tmodel.mvars['b'][0][rxn_index].ScenNLB=row[scen]
                    tmodel.mvars['b'][0][rxn_index].ScenNUB=row[scen]

        ##last is to change the objective, valid across scenarios
        #but will also have to iterate over met_names
        opt_dir=GRB.MAXIMIZE if opt==1 else GRB.MINIMIZE
        
        idx=tmodel.reactions.index(tmodel.reactions.get_by_id(flux_name_toVA))
        tmodel.m.setObjective(tmodel.mvars['v'][0][idx],opt_dir)

        tmodel.m.params.OutputFlag=0
        tmodel.m.optimize()
        
        ##get the results with the modified function check_metabolite_bounds
        try :
            df_fluxes.loc[opt]=get_scenObjVal(tmodel)
        except Exception as e :
            print(e)



    return df_fluxes
# df_cellconc=VA_met_getcellconc(tmodel,met_names,'akg')

#make a fuction to retrieve scennObjVAl for each scenario
def get_scenObjVal(tmodel):
    numscen=tmodel.m.numScenarios
    df_scenObjVal=pd.DataFrame(index=[i for i in range(numscen)],columns=['scenObjVal'])
    for scen in range(numscen):
        tmodel.m.params.ScenarioNumber=scen
        tmodel.m.update()
        df_scenObjVal.loc[scen]=tmodel.m.ScenNObjVal
    return df_scenObjVal.T.values[0]