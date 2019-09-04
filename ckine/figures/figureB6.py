"""
This creates Figure 6.
"""
import string
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
from scipy.optimize import least_squares
from .figureCommon import subplotLabel, getSetup, plot_conf_int
from .figureB1 import runIL2simple
from ..model import receptor_expression
from ..imports import import_muteins, import_Rexpr, import_samples_2_15

dataMean, _ = import_muteins()
dataMean.reset_index(inplace=True)
data, _, _ = import_Rexpr()
data.reset_index(inplace=True)
unkVec_2_15, _ = import_samples_2_15(N=5)


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((14, 8), (4, 8))

    for ii, item in enumerate(ax):
        if ii < 26:
            subplotLabel(item, string.ascii_uppercase[ii])
        else:
            subplotLabel(item, 'A' + string.ascii_uppercase[ii - 26])

    tps = np.array([0.5, 1.0, 2.0, 4.0]) * 60.0
    muteinC = dataMean.Concentration.unique()
    dataMean["Concentration"] = np.log10(dataMean["Concentration"])  # logscale for plotting
    
    print(dataMean.loc[(dataMean["Cells"] == 'T-reg') & (dataMean["Ligand"] == 'IL2-097')])

    axis = 0
    ligand_order = ['IL2-060', 'IL2-062', 'IL2-088', 'IL2-097']
    cell_order = ['NK', 'CD8+', 'T-reg', 'Naive Treg', 'Mem Treg', 'T-helper', 'Naive Th', 'Mem Th']

    pred_data = np.zeros((len(muteinC), len(tps), unkVec_2_15.shape[1]))  # make empty array for predicted data at each concentration, tp, and parameter sample for a single cell type
    df = pd.DataFrame(columns=['Mutein', 'Activity Type', 'Cell Type', 'Time Point', 'Concentration', 'Replicate', 'Activity']) # make empty dataframe for all cell types
    
    iter = len(tps) * len(muteinC)

    # loop for each mutein and cell type
    for i, ligand_name in enumerate(ligand_order):
        for j, cell_name in enumerate(cell_order):
            
            # organize experimental pstat data
            exp_data = np.zeros((12, 4))
            mutein_conc = np.zeros(48,)
            for k, conc in enumerate(dataMean.Concentration.unique()):
                exp_data[k, :] = np.array(dataMean.loc[(dataMean["Cells"] == cell_name) & (dataMean["Ligand"] == ligand_name) & (dataMean["Concentration"] == conc), "RFU"])
                mutein_conc[(k * 4):((k + 1) * 4)] = conc
                
            df2 = pd.DataFrame({'Mutein':np.tile(np.array(ligand_name), iter), 'Activity Type':np.tile(np.array('experimental'), iter), 'Cell Type':np.tile(np.array(cell_name), iter), 'Time Point':np.tile(np.array([0.5, 1.0, 2.0, 4.0]), len(muteinC)), 'Concentration':mutein_conc, 'Replicate': np.zeros(iter), 'Activity':exp_data.reshape(iter,)})
            #print(df2)
            df = df.append(df2, ignore_index=True)

            # calculate predicted dose response
            IL2Ra = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\alpha$'), "Count"].item()
            IL2Rb = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\beta$'), "Count"].item()
            gc = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == '$\\gamma_{c}$'), "Count"].item()
            
            for l in range(unkVec_2_15.shape[1]):
                
                cell_receptors = receptor_expression(np.array([IL2Ra, IL2Rb, gc]).astype(np.float), unkVec_2_15[17, l], unkVec_2_15[20, l], unkVec_2_15[19, l], unkVec_2_15[21, l]) 
                pred_data[:, :, l] = calc_dose_response_mutein(unkVec_2_15[:, l], [1., 1., 5.], tps, muteinC, cell_receptors, exp_data)
                df3 = pd.DataFrame({"Activity Type":np.tile(np.array('predicted'), iter), 'Cell Type':np.tile(np.array(cell_name), iter), 'Replicate':np.tile(np.array(1 + l), iter), 'Activity':pred_data[:, :, l].reshape(iter,)})
                df = df.append(df3, ignore_index=True)
                            

            # plot experimental and predicted date with a legend for the last subplot
            axis = i * 8 + j
            if axis == 15:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name) &
                                                                                          (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, palette=cm.rainbow, legend='full')
                ax[axis].legend(loc='lower right', title="time (hours)")
            else:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name)
                                                                                          & (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, palette=cm.rainbow, legend=False)
            plot_dose_response(ax[axis], pred_data, tps, muteinC)
            ax[axis].set(xlabel=("[" + ligand_name + "] (log$_{10}$[nM])"), ylabel="Activity", title=cell_name)

    return f

def mutein_scaling(df):
    return 0


def calc_dose_response_mutein(unkVec, input_params, tps, muteinC, cell_receptors, exp_data):
    """ Calculates activity for a given cell type at various mutein concentrations and timepoints. """

    total_activity = np.zeros((len(muteinC), len(tps)))

    # loop for each mutein concentration
    for i, conc in enumerate(muteinC):
        active_ckine = runIL2simple(unkVec, input_params, conc, tps=tps, input_receptors=cell_receptors, adj_receptors=True)
        total_activity[i, :] = np.reshape(active_ckine, (-1, 4))  # save the activity from this concentration for all 4 tps

    # scale receptor/cell measurements to pSTAT activity for each sample
    #scale1, scale2 = optimize_scale(total_activity[:, :], exp_data)  # find optimal constants
    #total_activity[:, :] = scale2 * total_activity[:, :] / (total_activity[:, :] + scale1)  # adjust activity
    
    print('total_activity.shape: ', total_activity.shape)

    return total_activity


def plot_dose_response(ax, mutein_activity, tps, muteinC):
    """ Plots predicted activity for multiple timepoints and mutein concentrations. """
    colors = cm.rainbow(np.linspace(0, 1, tps.size))

    for tt in range(tps.size):
        plot_conf_int(ax, np.log10(muteinC.astype(np.float)), mutein_activity[:, tt, :], colors[tt])


def optimize_scale(model_act, exp_act):
    """ Formulates the optimal scale to minimize the residual between model activity predictions and experimental activity measurments for a given cell type. """

    # scaling factors are sigmoidal and linear, respectively
    guess = np.array([100.0, np.mean(exp_act) / np.mean(model_act)])

    def calc_res(sc):
        """ Calculate the residuals. This is the function we minimize. """
        scaled_act = sc[1] * model_act / (model_act + sc[0])
        err = exp_act - scaled_act
        return err.flatten()

    # find result of minimization where both params are >= 0
    res = least_squares(calc_res, guess, bounds=(0.0, np.inf))
    return res.x
