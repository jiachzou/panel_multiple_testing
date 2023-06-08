import pandas as pd
import numpy as np

def panel_unordered(posi_log_pval_matrix):
    '''
    # Input: a matrix of log p values. dimension J by N. J: number of features. N: number of units.
    # Output: a pd.DataFrame of sorted multiple testing evidence.
    '''
    if type(posi_log_pval_matrix)!=pd.core.frame.DataFrame:
        posi_log_pval_matrix= pd.DataFrame(posi_log_pval_matrix)

    posi_log_pval_matrix = posi_log_pval_matrix.fillna(float('inf'))
    
    time_series_Bonf_rejection_table = pd.DataFrame(False, index=posi_log_pval_matrix.index, 
                                                    columns=['rho_inv.N.p_1','rho_inv.N','N','p_1'])

    K_set = (posi_log_pval_matrix != float('inf')).sum(axis=1)
    M_set = (posi_log_pval_matrix != float('inf')).sum(axis=0)
    N_vec = []

    for i_row in range(posi_log_pval_matrix.shape[0]):
        this_pval_row = posi_log_pval_matrix.iloc[i_row,:]
        meaningful_ind = ~(this_pval_row == float('inf'))
        N_vec.append(M_set[meaningful_ind].sum())
    N_vec = np.array(N_vec)
    rho_inv = (K_set[N_vec>0] / N_vec[N_vec>0]).sum()
    rho = 1 / rho_inv
    
    for i_row in range(posi_log_pval_matrix.shape[0]):
        this_pval_row = posi_log_pval_matrix.iloc[i_row,:]
        meaningful_ind = ~(this_pval_row == float('inf'))

        if meaningful_ind.sum() > 0:
            my_df = N_vec[i_row]
            p_1 = this_pval_row[meaningful_ind].min()

            bonf_level = np.exp(p_1) * my_df * rho_inv
            time_series_Bonf_rejection_table.iloc[i_row,:] = [bonf_level,my_df*rho_inv,my_df,p_1]
        else:
            time_series_Bonf_rejection_table.iloc[i_row,:] = [np.nan,0,0,np.nan]

    time_series_Bonf_rejection_table['p_1'] = np.exp(time_series_Bonf_rejection_table['p_1'].astype(float))
    time_series_Bonf_rejection_table = time_series_Bonf_rejection_table.sort_values('rho_inv.N.p_1')
    time_series_Bonf_rejection_table['rho'] = rho

    return time_series_Bonf_rejection_table
