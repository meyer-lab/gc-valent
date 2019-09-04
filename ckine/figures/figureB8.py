"""
This creates Figure 8. Comparison of Experimental verus Predicted Activity across IL2 and IL15 concentrations.
"""

import string
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from scipy.optimize import least_squares
from scipy.stats import pearsonr
from .figureCommon import subplotLabel, getSetup
from .figureB6 import calc_dose_response_mutein
from ..model import receptor_expression
from ..imports import import_muteins, import_Rexpr, import_samples_2_15

dataMean,_ =import_muteins()
dataMean.reset_index(inplace=True)
data, _, cellNames =import_Rexpr()
data.reset_index(inplace=True)
unkVec_2_15, scales = import_samples_2_15(N=1)  # use one rate
receptor_data, _, cell_names_receptor = import_Rexpr()
muteinC = dataMean.Concentration.unique()
tps = np.array([0.5, 1., 2., 4.]) * 60.



def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 3), multz={0:2})

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
        
    ligand_order = ['IL2-060', 'IL2-062', 'IL2-088', 'IL2-097']
    cell_order = ['NK', 'CD8+', 'T-reg', 'Naive Treg', 'Mem Treg', 'T-helper', 'Naive Th', 'Mem Th']
    
    # main routine for EC-50 analysis
    df = pd.DataFrame(columns=['Time Point', 'Cell Type', 'Mutein', 'Data Type', 'EC50'])

    x0 = [1, 2., 1000.]
    data_types = []
    cell_types = []
    mutein_types = []
    EC50s = np.zeros(len(cell_order) * len(tps) * 4 * 2) # EC50 for all cell types, tps, muteins, and expr/pred
    
    for i, cell_name in enumerate(cell_order):
                
        IL2Ra = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\alpha$'), "Count"].item()
        IL2Rb = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\beta$'), "Count"].item()
        gc = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == '$\\gamma_{c}$'), "Count"].item()
        cell_receptors = receptor_expression(np.array([IL2Ra, IL2Rb, gc]).astype(np.float), unkVec_2_15[17], unkVec_2_15[20], unkVec_2_15[19], unkVec_2_15[21])

        for j, ligand_name in enumerate(ligand_order):
            
            celltype_data = np.zeros((12, 4))
            for k, conc in enumerate(dataMean.Concentration.unique()):
                celltype_data[k, :] = np.array(dataMean.loc[(dataMean["Cells"] == cell_name) & (dataMean["Ligand"] == ligand_name) & (dataMean["Concentration"] == conc), "RFU"])
                            
            # predicted EC50
            pred_EC50 = calculate_predicted_EC50(x0, cell_receptors, tps, celltype_data)
            for l, item in enumerate(pred_EC50):
                EC50s[(8 * j) + (32 * i) + l] = item
            data_types.extend(np.tile(np.array('Predicted'), len(tps)))

            # experimental EC50
            for m, _ in enumerate(tps):
                timepoint_data = celltype_data[:, m]
                EC50s[(8 * j) + (32 * i) + len(tps) + m] = nllsq_EC50(x0, np.log10(muteinC * 10**4), timepoint_data)
            data_types.extend(np.tile(np.array('Experimental'), len(tps)))
            cell_types.extend(np.tile(np.array(cell_name), len(tps) * 2))  # for both experimental and predicted
            mutein_types.extend(np.tile(np.array(ligand_name), len(tps) * 2))

    EC50s = EC50s - 4  # account for 10^4 multiplication
    dataframe = {'Time Point':np.tile(tps, len(cell_order) * len(ligand_order) * 2), 'Mutein':mutein_types, 'Cell Type':cell_types, 'Data Type':data_types, 'EC-50':EC50s}
    df = pd.DataFrame(dataframe)

    catplot_comparison(ax[0], df)  # compare experiments to model predictions

    return f


def catplot_comparison(ax, df):
    """ Construct EC50 catplots for each time point for IL2 and IL15. """
    # set a manual color palette
    col_list = ["violet", "goldenrod"]
    col_list_palette = sns.xkcd_palette(col_list)
    sns.set_palette(col_list_palette)
    # plot predicted EC50
    sns.catplot(x="Cell Type", y="EC-50", hue="Mutein",
                data=df.loc[(df['Time Point'] == 60.) & (df["Data Type"] == 'Predicted')],
                legend=True, legend_out=False, ax=ax, marker='^')

    # plot experimental EC50
    sns.catplot(x="Cell Type", y="EC-50", hue="Mutein",
                data=df.loc[(df['Time Point'] == 60.) & (df["Data Type"] == 'Experimental')],
                legend=True, legend_out=False, ax=ax, marker='o')

    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, rotation_mode="anchor", ha="right", position=(0, 0.02))
    ax.set_xlabel("")  # remove "Cell Type" from xlabel
    ax.set_ylabel(r"EC-50 (log$_{10}$[nM])")


def plot_corrcoef(ax, tps):
    """ Plot correlation coefficients between predicted and experimental data for all cell types. """
    corr_coefs = np.zeros(2 * len(cell_names_receptor))
    for i, _ in enumerate(cell_names_receptor):
        assert cell_names_receptor[i] == cell_names_pstat[i]
        experimental_2 = IL2_data_avg[(i * 4):((i + 1) * 4)]
        experimental_15 = IL15_data_avg[(i * 4):((i + 1) * 4)]
        predicted_2, predicted_15 = calc_dose_response(unkVec_2_15, scales, receptor_data[i], tps, muteinC, experimental_2, experimental_15)
        corr_coef2 = pearsonr(experimental_2.flatten(), np.squeeze(predicted_2).T.flatten())
        corr_coef15 = pearsonr(experimental_15.flatten(), np.squeeze(predicted_15).T.flatten())
        corr_coefs[i] = corr_coef2[0]
        corr_coefs[len(cell_names_receptor) + i] = corr_coef15[0]

    x_pos = np.arange(len(cell_names_receptor))
    ax.bar(x_pos - 0.15, corr_coefs[0:len(cell_names_receptor)], width=0.3, color='darkorchid', label='IL2', tick_label=cell_names_receptor)
    ax.bar(x_pos + 0.15, corr_coefs[len(cell_names_receptor):(2 * len(cell_names_receptor))], width=0.3, color='goldenrod', label='IL15', tick_label=cell_names_receptor)
    ax.set(ylabel=("Correlation"), ylim=(0., 1.))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, rotation_mode="anchor", ha="right", position=(0, 0.02))


def global_legend(ax):
    """ Create legend for colors and markers in subplots A-C. """
    purple = mpatches.Patch(color='darkorchid', label='IL-2')
    yellow = mpatches.Patch(color='goldenrod', label='IL-15')
    circle = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=6, label='Experimental')
    triangle = mlines.Line2D([], [], color='black', marker='^', linestyle='None', markersize=6, label='Predicted')
    ax.legend(handles=[purple, yellow, circle, triangle], bbox_to_anchor=(1.02, 1), loc="upper left")


def calculate_predicted_EC50(x0, cell_receptor_data, tps, pstat):
    """ Calculate average EC50 from model predictions. """
    activity = calc_dose_response_mutein(unkVec_2_15, [1., 1., 5.], tps, muteinC, cell_receptor_data, pstat)
    EC50 = np.zeros(len(tps))
    # calculate EC50 for each timepoint
    for i, _ in enumerate(tps):
        EC50[i] = nllsq_EC50(x0, np.log10(muteinC.astype(np.float) * 10**4), activity[:, i])
    return EC50


def nllsq_EC50(x0, xdata, ydata):
    """ Performs nonlinear least squares on activity measurements to determine parameters of Hill equation and outputs EC50. """
    lsq_res = least_squares(residuals, x0, args=(xdata, ydata), bounds=([0., 0., 0.], [10., 10., 10**5.]), jac='3-point')
    return lsq_res.x[0]


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
