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
from .figureB6 import calc_dose_response_mutein, organize_expr_pred, mutein_scaling
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

mutaff = {
    "IL2-060": [1., 1., 5.],  # Wild-type, but dimer
    "IL2-062": [1., 100., 5.],  # Weaker b-g
    "IL2-088": [10., 1., 5.],  # Weaker CD25
    "IL2-097": [10., 100., 5.]  # Both
}

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 3), multz={0:2})

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    ligand_order = ['IL2-060', 'IL2-062', 'IL2-088', 'IL2-097']
    cell_order = ['NK', 'CD8+', 'T-reg', 'Naive Treg', 'Mem Treg', 'T-helper', 'Naive Th', 'Mem Th']

    df = pd.DataFrame(columns=['Cells', 'Ligand', 'Time Point', 'Concentration', 'Activity Type', 'Replicate', 'Activity'])  # make empty dataframe for all cell types

    # loop for each cell type and mutein
    for _, cell_name in enumerate(cell_order):

        IL2Ra = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\alpha$'), "Count"].item()
        IL2Rb = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\beta$'), "Count"].item()
        gc = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == '$\\gamma_{c}$'), "Count"].item()
        receptors = np.array([IL2Ra, IL2Rb, gc]).astype(np.float)

        for _, ligand_name in enumerate(ligand_order):

            # append dataframe with experimental and predicted activity
            df = organize_expr_pred(df, cell_name, ligand_name, receptors, muteinC, tps, unkVec_2_15)

    # determine scaling constants
    scales = np.squeeze(mutein_scaling(df, unkVec_2_15))
    
    EC50_df = calculate_EC50s(df, scales, cell_order, ligand_order)
    
    catplot_comparison(ax[0], EC50_df)  # compare experiments to model predictions
    
    return f


def catplot_comparison(ax, df):
    """ Construct EC50 catplots for 1 hour time point. """
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


def calculate_EC50s(df, scales, cell_order, ligand_order):
    """ Scales model predictions to experimental data for all cell types, muteins, and time points. """
    
    x0 = [1, 2., 1000.]
    data_types = []
    cell_types = []
    mutein_types = []
    EC50s = np.zeros(len(cell_order) * len(tps) * 4 * 2) # EC50 for all cell types, tps, muteins, and expr/pred
    pred_data = np.zeros((12, 4))
    expr_data = pred_data.copy()
    cell_groups = [['T-reg', 'Mem Treg', 'Naive Treg'], ['T-helper', 'Mem Th', 'Naive Th'], ['NK'], ['CD8+']]

    for i, cell_name in enumerate(cell_order):
        for j, ligand_name in enumerate(ligand_order):
            for k, conc in enumerate(df.Concentration.unique()):
                for l, tp in enumerate(tps):
                    pred_data[k, l] = df.loc[(df["Cells"] == cell_name) & (df["Ligand"] == ligand_name) & (
                            df["Activity Type"] == 'predicted') & (df["Concentration"] == conc) & (df["Time Point"] == tp), "Activity"]
                    expr_data[k, l] = df.loc[(df["Cells"] == cell_name) & (df["Ligand"] == ligand_name) & (
                            df["Activity Type"] == 'experimental') & (df["Concentration"] == conc) & (df["Time Point"] == tp), "Activity"]

            # scale predicted data
            for m, cell_names in enumerate(cell_groups):
                if cell_name in cell_names:
                    pred_data[:, :] = scales[m, 1] * pred_data[:, :] / (pred_data[:, :] + scales[m, 0])
            
            # calculate predicted and experimental EC50s for all time points
            for n, _ in enumerate(tps):
                EC50s[(8 * j) + (32 * i) + n] = nllsq_EC50(x0, np.log10(muteinC.astype(np.float) * 10**4), pred_data[:, n])
                EC50s[(8 * j) + (32 * i) + len(tps) + n] = nllsq_EC50(x0, np.log10(muteinC.astype(np.float) * 10**4), expr_data[:, n])
                
            data_types.extend(np.tile(np.array('Predicted'), len(tps)))
            data_types.extend(np.tile(np.array('Experimental'), len(tps)))
            cell_types.extend(np.tile(np.array(cell_name), len(tps) * 2))  # for both experimental and predicted
            mutein_types.extend(np.tile(np.array(ligand_name), len(tps) * 2))
            
    EC50s = EC50s - 4  # account for 10^4 multiplication
    dataframe = {'Time Point':np.tile(tps, len(cell_order) * len(ligand_order) * 2), 'Mutein':mutein_types, 'Cell Type':cell_types, 'Data Type':data_types, 'EC-50':EC50s}
    df = pd.DataFrame(dataframe)                                                           
                                                                      
    return df


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
