"""
This creates Figure 4. Comparison of Experimental verus Predicted Activity across IL2 and IL15 concentrations.
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
    tps = np.array([0.5, 1., 2., 4.]) * 60.
    data_types = []
    cell_types = []
    EC50s_2 = np.zeros(len(cell_names_pstat) * len(tps) * 4)
    EC50s_15 = np.zeros(len(cell_names_pstat) * len(tps) * 4)

    for i, name in enumerate(cell_names_pstat):
        assert cell_names_pstat[i] == cell_names_receptor[i]
        celltype_data_2 = IL2_data[(i * 4):((i + 1) * 4)]
        celltype_data_2b = IL2_data2[(i * 4):((i + 1) * 4)]  # second replicate
        celltype_data_15 = IL15_data[(i * 4):((i + 1) * 4)]
        celltype_data_15b = IL15_data2[(i * 4):((i + 1) * 4)]  # second replicate
        data_types.append(np.tile(np.array('Predicted'), len(tps)))
        # predicted EC50
        # first replicate
        EC50_2, EC50_15 = calculate_predicted_EC50(x0, receptor_data[i], tps, celltype_data_2, celltype_data_15)
        for j, item in enumerate(EC50_2):
            EC50s_2[(2 * len(tps) * i) + j] = item
            EC50s_15[(2 * len(tps) * i) + j] = EC50_15[j]
        # second replicate: missing two cell types
        if i <= 1 or i >= 4:
            # handle cases with missing IL15 measurements
            if i == 0 or i == 1:
                celltype_data_15b[0, 10] = celltype_data_15b[0, 11] = 0  # missing data set to zero for scaling
                EC50_2b, EC50_15b = calculate_predicted_EC50(x0, receptor_data[i], tps, celltype_data_2b, celltype_data_15b)
                for k, item in enumerate(EC50_2b):
                    EC50s_2[(len(tps) * len(cell_names_pstat) * 2) + (2 * len(tps) * i) + k] = item
                    EC50s_15[(len(tps) * len(cell_names_pstat) * 2) + (2 * len(tps) * i) + k] = EC50_15b[k]
            else:
                EC50_2b, EC50_15b = calculate_predicted_EC50(x0, receptor_data[i], tps, celltype_data_2b, celltype_data_15b)
                for l, item in enumerate(EC50_2b):
                    EC50s_2[(len(tps) * len(cell_names_pstat) * 2) + (2 * len(tps) * i) + l] = item
                    EC50s_15[(len(tps) * len(cell_names_pstat) * 2) + (2 * len(tps) * i) + l] = EC50_15b[l]
        # experimental EC50
        for m, _ in enumerate(tps):
            timepoint_data_2 = celltype_data_2[m]
            timepoint_data_2b = celltype_data_2b[m]
            timepoint_data_15 = celltype_data_15[m]
            timepoint_data_15b = celltype_data_15b[m]
            # first replicate
            EC50s_2[len(tps) + (2 * len(tps) * i) + m] = nllsq_EC50(x0, np.log10(ckineConc.astype(np.float) * 10**4), timepoint_data_2)
            EC50s_15[len(tps) + (2 * len(tps) * i) + m] = nllsq_EC50(x0, np.log10(ckineConc.astype(np.float) * 10**4), timepoint_data_15)
            # second replicate: missing two cell types
            if i <= 1 or i >= 4:
                # handle cases with missing IL15 measurements
                if i == 0 or i == 1 and m == 0:
                    timepoint_data_15b_mod = np.delete(timepoint_data_15b, [10, 11])
                    ckineConc_mod = np.delete(ckineConc, [10, 11])
                    EC50s_15[(len(tps) * len(cell_names_pstat) * 2) + len(tps) + (2 * len(tps) * i) + m] = nllsq_EC50(x0, np.log10(ckineConc_mod.astype(np.float) * 10**4), timepoint_data_15b_mod)
                else:
                    EC50s_15[(len(tps) * len(cell_names_pstat) * 2) + len(tps) + (2 * len(tps) * i) + m] = nllsq_EC50(x0, np.log10(ckineConc.astype(np.float) * 10**4), timepoint_data_15b)
                EC50s_2[(len(tps) * len(cell_names_pstat) * 2) + len(tps) + (2 * len(tps) * i) + m] = nllsq_EC50(x0, np.log10(ckineConc.astype(np.float) * 10**4), timepoint_data_2b)
        data_types.append(np.tile(np.array('Experimental'), len(tps)))
        cell_types.append(np.tile(np.array(name), len(tps) * 2))  # for both experimental and predicted

    EC50 = np.concatenate((EC50s_2, EC50s_15), axis=None)
    EC50 = EC50 - 4  # account for 10^4 multiplication
    data_types = np.tile(np.array(data_types).reshape(80,), 4)  # for IL2 and IL15, 2 replicates
    cell_types = np.tile(np.array(cell_types).reshape(80,), 4)
    IL = np.concatenate((np.tile(np.array('IL2'), len(cell_names_pstat) * len(tps) * 4), np.tile(np.array('IL15'), len(cell_names_pstat) * len(tps) * 4)), axis=None)
    data = {'Time Point': np.tile(np.array(tps), len(cell_names_pstat) * 8), 'IL': IL, 'Cell Type': cell_types.reshape(320,), 'Data Type': data_types.reshape(320,), 'EC-50': EC50}
    df = pd.DataFrame(data)
    df = df.loc[df['EC-50'] != -4.]  # delete cells that did not have a valid replicate (0 - 4)

    catplot_comparison(ax, df, tps)

    return f


def catplot_comparison(ax, df, tps):
    """ Construct EC50 catplots for each time point for IL2 and IL15. """
    tps_str = ['30 mins', '1 hr', '2 hrs', '4 hrs']
    for i, tp in enumerate(tps):
        sns.catplot(x="Cell Type", y="EC-50", hue="Data Type", data=df.loc[(df['Time Point'] == tp) & (df["IL"] == 'IL2')], legend=False, ax=ax[i])
        ax[i].get_legend().set_visible(False)
        if i == 3:
            sns.catplot(x="Cell Type", y="EC-50", hue="Data Type", data=df.loc[(df['Time Point'] == tp) & (df["IL"] == 'IL15')], legend=True, legend_out=True, ax=ax[4 + i])
            ax[4 + i].legend(bbox_to_anchor=(1.05, 1), loc='right')
        else:
            sns.catplot(x="Cell Type", y="EC-50", hue="Data Type", data=df.loc[(df['Time Point'] == tp) & (df["IL"] == 'IL15')], legend=False, ax=ax[4 + i])
            ax[4 + i].get_legend().set_visible(False)
        ax[i].set(ylabel=('IL2 log[EC50] (' + tps_str[i] + ')'), ylim=(-3., 1.))
        ax[4 + i].set(ylabel=('IL15 log[EC50] (' + tps_str[i] + ')'), ylim=(-3., 1.))
        ax[i].set_xticklabels(ax[i].get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)
        ax[4 + i].set_xticklabels(ax[4 + i].get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)


def calculate_predicted_EC50(x0, cell_receptor_data, tps, IL2_pstat, IL15_pstat):
    """ Calculate average EC50 from model predictions. """
    IL2_activity, IL15_activity = calc_dose_response(unkVec_2_15, scales, cell_receptor_data, tps, ckineConc, IL2_pstat, IL15_pstat)
    avg_EC50_2 = np.zeros(len(tps))
    avg_EC50_15 = np.zeros(len(tps))
    for i, _ in enumerate(tps):
        # take average for all samples
        total_EC50_2 = 0
        total_EC50_15 = 0
        for j in range(IL2_activity.shape[1]):
            total_EC50_2 = total_EC50_2 + nllsq_EC50(x0, np.log10(ckineConc.astype(np.float) * 10**4), IL2_activity[:, j, i])
            total_EC50_15 = total_EC50_15 + nllsq_EC50(x0, np.log10(ckineConc.astype(np.float) * 10**4), IL15_activity[:, j, i])
        avg_EC50_2[i] = total_EC50_2 / IL2_activity.shape[1]
        avg_EC50_15[i] = total_EC50_15 / IL15_activity.shape[1]
    return avg_EC50_2, avg_EC50_15


def nllsq_EC50(x0, xdata, ydata):
    """ Performs nonlinear least squares on activity measurements to determine parameters of Hill equation and outputs EC50. """
    lsq_res = least_squares(residuals, x0, args=(xdata, ydata), bounds=([0., 0., 0.], [10., 10., 10**5.]), jac='3-point')
    y = hill_equation(xdata, lsq_res.x)
    y_halfMax = np.amax(y) / 2
    halfMaxConc = fsolve(hill_equation, 2, args=(lsq_res.x, y_halfMax))
    return halfMaxConc


def hill_equation(x, x0, solution=0):
    """ Calculates EC50 from Hill Equation. """
    k = x0[0]
    n = x0[1]
    A = x0[2]
    xk = np.power(x / k, n)
    return (A * xk / (1.0 + xk)) - solution


def residuals(x0, x, y):
    """ Residual function for Hill Equation. """
    return hill_equation(x, x0) - y
