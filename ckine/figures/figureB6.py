"""
This creates Figure 6.
"""
import string
import numpy as np
import seaborn as sns
from .figureCommon import subplotLabel, getSetup
from .figureB1 import runIL2simple
from ..model import receptor_expression
from ..imports import import_muteins, import_Rexpr
from ..make_tensor import rxntfR


dataMean, _ = import_muteins()
dataMean.reset_index(inplace=True)
dataMean['Concentration'] = np.log10(dataMean['Concentration'])
data, _, _ = import_Rexpr()
data.reset_index(inplace=True)

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 5), empts=[18, 19])
    
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
    
    tps = np.array([0.5, 1.0, 2.0, 4.0]) * 60.0
    muteinC = dataMean.Concentration.unique()
    axis = 0
    
    first_group_ligands = ['IL2-060', 'IL2-062']
        
    for i, ligand_name in enumerate(first_group_ligands):
        for j, cell_name in enumerate(dataMean.Cells.unique()):
            print(cell_name)
            IL2Ra = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\alpha$'), "Count"]
            IL2Rb = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\beta$'), "Count"].item()
            gc = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == '$\\gamma_{c}$'), "Count"].item()
            cell_receptors = receptor_expression(np.array([IL2Ra, IL2Rb, gc]).astype(np.float), rxntfR[17], rxntfR[20], rxntfR[19], rxntfR[21])
            predicted = calc_dose_response_mutein([1., 1., 5.], tps, muteinC, cell_receptors) #TODO: Verify keeping 5x weaker endosomal assumption
            axis = i*9+j
            if axis == 17:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name) & (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, legend='full')
                ax[axis].legend(loc='lower right', title="time (hours)")
            else:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name) & (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, legend=False)
            ax[axis].set(xlabel=("[" + ligand_name + "] (log$_{10}$[nM])"), ylabel="Activity", title=cell_name)
    return f

def calc_dose_response_mutein(input_params, tps, muteinC, cell_receptors):
    """ Calculates activity for a given cell type at various mutein concentrations and timepoints. """
    
    total_activity = np.zeros((len(muteinC), len(tps)))

    # loop for each mutein concentration
    for i, conc in enumerate(muteinC):
        yOut = runIL2simple(input_params, conc, tps=tps, input_receptors=cell_receptors, adj_receptors=True) #see if can do mult. tps
        activity = np.dot(yOut, getTotalActiveSpecies().astype(np.float))
        total_activity[i, :] = np.reshape(activity1, (-1, 4))  # save the activity from this concentration for all 4 tps
    
    """
    # scale receptor/cell measurements to pSTAT activity for each sample
    for j in range(len(scales)):
        scale1, scale2 = optimize_scale(total_activity1[:, j, :], total_activity2[:, j, :], exp_data_1, exp_data_2)  # find optimal constants
        total_activity1[:, j, :] = scale2 * total_activity1[:, j, :] / (total_activity1[:, j, :] + scale1)  # adjust activity for this sample
        total_activity2[:, j, :] = scale2 * total_activity2[:, j, :] / (total_activity2[:, j, :] + scale1)  # adjust activity for this sample
    """

    return total_activity


def plot_dose_response(ax1, ax2, mutein1_activity, mutein2_activity, cell_type, tps, muteinC, mutein_name, legend=False):
    """ Plots both mutein activity in different plots where each plot has multiple timepoints and mutein concentrations. """
    colors = cm.rainbow(np.linspace(0, 1, tps.size))

    # plot the values with each time as a separate color
    for tt in range(tps.size):
        plot_conf_int(ax1, np.log10(muteinC.astype(np.float)), activity1[:, tt], colors[tt])  # never a legend for first mutein
        if legend:
            plot_conf_int(ax2, np.log10(muteinC.astype(np.float)), activity2[:, tt], colors[tt], (tps[tt] / 60.0).astype(str))
            ax2.legend(title="time (hours)")
        else:
            plot_conf_int(ax2, np.log10(muteinC.astype(np.float)), activity2[:, tt], colors[tt])

    # plots for input cell type
    ax1.set(xlabel=mutein_name[0] + "(log$_{10}$[nM])", ylabel="Activity", title=cell_type)
    ax2.set(xlabel=mutein_name[1] + "(log$_{10}$[nM])", ylabel="Activity", title=cell_type)


def optimize_scale(model_act2, model_act15, exp_act2, exp_act15):
    """ Formulates the optimal scale to minimize the residual between model activity predictions and experimental activity measurments for a given cell type. """
    exp_act2 = exp_act2.T  # transpose to match model_act
    exp_act15 = exp_act15.T

    # scaling factors are sigmoidal and linear, respectively
    guess = np.array([100.0, np.mean(exp_act2 + exp_act15) / np.mean(model_act2 + model_act15)])

    def calc_res(sc):
        """ Calculate the residuals. This is the function we minimize. """
        scaled_act2 = sc[1] * model_act2 / (model_act2 + sc[0])
        scaled_act15 = sc[1] * model_act15 / (model_act15 + sc[0])
        err = np.hstack((exp_act2 - scaled_act2, exp_act15 - scaled_act15))
        return np.reshape(err, (-1,))

    # find result of minimization where both params are >= 0
    res = least_squares(calc_res, guess, bounds=(0.0, np.inf))
    return res.x
