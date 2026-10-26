#-------------------------------------------------------------------------------
# Name:        force_flux.py

# Purpose:     Helper functinos to build force flux model
###             ralg lin constraint are Ax<=b. in freeflux, we build A and b then take the opposite. So here give -A and -b to freeflux func.

# 
# 
# Author:      Nathan Fargier
#
# Created:     07/03/2025
#-------------------------------------------------------------------------------

import sys
sys.path.append('/home/nathan/thermo_flux/')
sys.path.append('/home4/p313427/freeflux/')


import numpy as np


##import pinv2


def _find_vplus_vminus(optModel,row):
    if optModel.T[row].sum() != 0:
        #return the column with nn zero value for this row
        idx,=np.where(optModel.T[row]!=0)
        return idx
    else: 
        idx_plus,=np.where(optModel.T[row]==1)
        idx_minus,=np.where(optModel.T[row]==-1)
        return idx_plus, idx_minus
    


def build_ff_matrix(optModel,fit):
    '''build the matrixes A4,A5,b4,b5 for the force flux model'''
    #get the dg bounds
    dg_lb = optModel.dg_lb
    dg_ub = optModel.dg_ub
    
    dict_vplus_vminus={row:_find_vplus_vminus(optModel,row) for row in range(optModel.T.shape[0])}

    vpos_mask = {i:vnet for i,vnet in enumerate(fit.vnet) if vnet>0}
    vneg_mask= {i:vnet for i,vnet in enumerate(fit.vnet) if vnet<0}

    #the v- that we keep are the vexch corresponding to vnets>0 within the vtot
    list_index_minus = [dict_vplus_vminus[net_index][1] for net_index in vpos_mask.keys() 
                        if isinstance(dict_vplus_vminus[net_index], tuple) and len(dict_vplus_vminus[net_index]) == 2]
    #the v+ that we keep are the vexch corresponding to vnets<0 within the vtot
    list_index_plus = [dict_vplus_vminus[net_index][0] for net_index in vneg_mask.keys() 
                        if isinstance(dict_vplus_vminus[net_index], tuple) and len(dict_vplus_vminus[net_index]) == 2]

    #func to make the diagonal mask matrix
    def make_diag_mask(n,indexes):
        '''fill the diagonal of a matrix with 1 at the indexes given in the list indexes'''
        M=np.zeros((n,n))
        for i in indexes:
            M[i,i]=1
        return M

    M_neg=make_diag_mask(optModel.T.shape[1],list_index_minus) # gives indexes of v-
    M_pos=make_diag_mask(optModel.T.shape[1],list_index_plus) # gives indexes of v+

    ##construct opposite of A3 and A4 : M@N (minus sign because the ineq is vnet/x -vtot[i]   <=0) ## but minus sign is added in the build constraint function
    A4=M_neg@optModel.N
    A5=M_pos@optModel.N

    T = 303.15  # 30°C
    R = 0.008314  # Universal gas constant
    exp_neg_factors = np.exp(-dg_lb / (R * T))  # For vnet > 0
    exp_pos_factors = np.exp(dg_ub / (R * T)) 
    #clip minimum pos value of exp factors otherwise i end up with some b ~ 1e-(toomuch) 
    # exp_neg_factors[abs(exp_neg_factors)<1e-10]=1e-10
    # exp_pos_factors[abs(exp_pos_factors)<1e-10]=1e-10
    # # iterate over vpos and vneg masks and get the corresponding or dg ub
    b4_vals = {dict_vplus_vminus[net_index][1][0]: vnet_val / (exp_neg_factors[net_index] - 1) 
            for net_index, vnet_val in vpos_mask.items() 
            if isinstance(dict_vplus_vminus[net_index], tuple) and len(dict_vplus_vminus[net_index]) == 2} #only rev reactions

    b5_vals = {dict_vplus_vminus[net_index][0][0]: -  vnet_val / (exp_pos_factors[net_index] - 1)
            for net_index, vnet_val in vneg_mask.items() 
            if isinstance(dict_vplus_vminus[net_index], tuple) and len(dict_vplus_vminus[net_index]) == 2}


    #now make the b arrays
    def make_b_mat(b_vals,n):
        b=np.zeros(n)
        for index, val in b_vals.items():
            b[index]=val
        return b
    ##construct b4 and b5 : since openopt is solving Ax<=b,  minus sign is added in the build constraint functio
    b4=make_b_mat(b4_vals,optModel.T.shape[1])
    b5=make_b_mat(b5_vals,optModel.T.shape[1])

    print('clip forceflux b matrixes')
    b4[abs(b4)<1e-10]=0
    b5[abs(b5)<1e-10]=0

    return A4,A5,b4,b5

