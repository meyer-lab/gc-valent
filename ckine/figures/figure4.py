"""
This creates Figure 4.
"""

import string
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import least_squares, fsolve
from .figureCommon import subplotLabel, getSetup
from .figureS5 import calc_dose_response
from ..imports import import_pstat, import_Rexpr, import_samples_2_15

ckineConc, cell_names_pstat, IL2_data, IL15_data, IL2_data2, IL15_data2 = import_pstat()
unkVec_2_15, scales = import_samples_2_15(N=100)  # use all rates
_, receptor_data, cell_names_receptor = import_Rexpr()

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 4), multz={0: 1, 2: 1, 4: 1, 6: 1, 8: 1, 10: 1, 12: 1, 14: 1})

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
    
    df = pd.DataFrame(columns=['Time Point', 'Cell Type', 'IL', 'Data Type', 'EC50'])
    
    x0 = [1, 2., 1000.]
    tps = ['30 mins', '1 hr', '2 hrs', '4 hrs']
    tps_num = np.array([0.5, 1., 2., 4.]) * 60.
    data_types = []
    cell_types = []
    EC50s_2 = np.zeros(len(cell_names_pstat)*len(tps)*4)
    EC50s_15 = np.zeros(len(cell_names_pstat)*len(tps)*4)
    
    for j, name in enumerate(cell_names_pstat):
        # match cells
        assert cell_names_pstat[j] == cell_names_receptor[j]
        celltype_data_2 = IL2_data[(j * 4):((j + 1) * 4)]
        celltype_data_2b = IL2_data2[(j * 4):((j + 1) * 4)]
        celltype_data_15 = IL15_data[(j * 4):((j + 1) * 4)]
        celltype_data_15b = IL15_data2[(j * 4):((j + 1) * 4)]
        data_types.append(np.tile(np.array('Predicted'), len(tps)))
        # predicted EC50
        EC50_2, EC50_15 = calculate_experimental_EC50(x0, receptor_data[j], tps_num, celltype_data_2, celltype_data_15)
        if j >=4:
            EC50_2b, EC50_15b = calculate_experimental_EC50(x0, receptor_data[j], tps_num, celltype_data_2b, celltype_data_15b)
        for l, item in enumerate(EC50_2):
            EC50s_2[(2*len(tps)*j)+l] = item
            EC50s_2[(len(tps)*len(cell_names_pstat)*2)+(2*len(tps)*j)+l] = EC50_2b[l]
            EC50s_15[(2*len(tps)*j)+l] = EC50_15[l]
            EC50s_15[(len(tps)*len(cell_names_pstat)*2)+(2*len(tps)*j)+l] = EC50_15b[l]
        # experimental EC50
        for i, _ in enumerate(tps):
            timepoint_data_2 = celltype_data_2[i]
            timepoint_data_2b = celltype_data_2b[i]
            timepoint_data_15 = celltype_data_15[i]
            timepoint_data_15b = celltype_data_15b[i]
            EC50s_2[len(tps)+(2*len(tps)*j)+i] = nllsq(x0, np.log10(ckineConc.astype(np.float)*10**4), timepoint_data_2)
            EC50s_2[(len(tps)*len(cell_names_pstat)*2)+len(tps)+(2*len(tps)*j)+i] = nllsq(x0, np.log10(ckineConc.astype(np.float)*10**4), timepoint_data_2b)
            EC50s_15[len(tps)+(2*len(tps)*j)+i] = nllsq(x0, np.log10(ckineConc.astype(np.float)*10**4), timepoint_data_15)
            EC50s_15[(len(tps)*len(cell_names_pstat)*2)+len(tps)+(2*len(tps)*j)+i] = nllsq(x0, np.log10(ckineConc.astype(np.float)*10**4), timepoint_data_15)
        data_types.append(np.tile(np.array('Experimental'), len(tps)))
        cell_types.append(np.tile(np.array(name), len(tps)*2)) # for both experimental and predicted
    
    EC50 = np.concatenate((EC50s_2, EC50s_15), axis=None)
    EC50 = EC50 - 4 # account for 10^4 multiplication
    data_types = np.tile(np.array(data_types).reshape(80,), 4) #for IL2 and IL15, 2 replicates
    cell_types = np.tile(np.array(cell_types).reshape(80,), 4)
    IL = np.concatenate((np.tile(np.array('IL2'), len(cell_names_pstat)*len(tps)*4), np.tile(np.array('IL15'), len(cell_names_pstat)*len(tps)*4)), axis=None)
    print(np.tile(np.array(tps), len(cell_names_pstat)*8).shape, IL.shape, EC50.shape)
    data = {'Time Point':np.tile(np.array(tps), len(cell_names_pstat)*8), 'IL':IL, 'Cell Type':cell_types.reshape(320,), 'Data Type':data_types.reshape(320,), 'EC-50':EC50}
    df = pd.DataFrame(data) 
    
    print(df.loc[(df['Time Point'] == '30 mins') & (df["IL"] == 'IL2'), "Data Type"])
    
    catplot_comparison(ax, df, tps)

    return f


def catplot_comparison(ax, df, tps):
    for i, tp in enumerate(tps):
        sns.catplot(x="Cell Type", y="EC-50", hue="Data Type", data=df.loc[(df['Time Point'] == tp) & (df["IL"] == 'IL2')], legend=False, ax=ax[i])
        if i == 3:
            sns.catplot(x="Cell Type", y="EC-50", hue="Data Type", data=df.loc[(df['Time Point'] == tp) & (df["IL"] == 'IL15')], legend=True, legend_out=True, ax=ax[4+i])
        else:
            sns.catplot(x="Cell Type", y="EC-50", hue="Data Type", data=df.loc[(df['Time Point'] == tp) & (df["IL"] == 'IL15')], legend=False, ax=ax[4+i])
        ax[i].set_ylabel('IL2 log[EC50] ('+tp+')')
        ax[4+i].set_ylabel('IL15 log[EC50] ('+tp+')')
        ax[i].set_xticklabels(ax[i].get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)
        ax[4+i].set_xticklabels(ax[2*i].get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)


def calculate_experimental_EC50(x0, cell_receptor_data, tps, IL2_pstat, IL15_pstat):
    IL2_activity, IL15_activity = calc_dose_response(unkVec_2_15, scales, cell_receptor_data, tps, ckineConc, IL2_pstat, IL15_pstat)
    avg_EC50_2 = np.zeros(len(tps))
    avg_EC50_15 = np.zeros(len(tps))
    for i, _ in enumerate(tps):
        # take average for all samples
        total_EC50_2 = 0
        total_EC50_15 = 0
        for j in range(IL2_activity.shape[1]):
            total_EC50_2 = total_EC50_2 + nllsq(x0, np.log10(ckineConc.astype(np.float)*10**4), IL2_activity[:, j, i])
            total_EC50_15 = total_EC50_15 + nllsq(x0, np.log10(ckineConc.astype(np.float)*10**4), IL15_activity[:, j, i])
        avg_EC50_2[i] = total_EC50_2 / IL2_activity.shape[1]
        avg_EC50_15[i] = total_EC50_15 / IL15_activity.shape[1]
    return avg_EC50_2, avg_EC50_15


def nllsq(x0, xdata, ydata):
    lsq_res = least_squares(residuals, x0, args=(xdata, ydata), bounds=([0., 0., 0.], [10., 10., 10**5.]), jac='3-point')
    y = hill_equation(xdata, lsq_res.x)
    y_halfMax = np.amax(y) / 2
    halfMaxConc = fsolve(hill_equation, 2, args=(lsq_res.x, y_halfMax))
    return halfMaxConc


def hill_equation(x, x0, solution=0):
    k = x0[0]
    n = x0[1]
    A = x0[2]
    xk = np.power(x / k, n)
    return (A * xk / (1.0 + xk)) - solution


def residuals(x0, x, y):
    return hill_equation(x, x0) - y
