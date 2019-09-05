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

    ligand_order = ['IL2-060', 'IL2-062', 'IL2-088', 'IL2-097']
    cell_order = ['NK', 'CD8+', 'T-reg', 'Naive Treg', 'Mem Treg', 'T-helper', 'Naive Th', 'Mem Th']

    df = pd.DataFrame(columns=['Cells', 'Ligand', 'Time Point', 'Concentration', 'Activity Type', 'Replicate', 'Activity']) # make empty dataframe for all cell types
    
    num = len(tps) * len(muteinC)

    # loop for each mutein and cell type
    for _, cell_name in enumerate(cell_order):
        
        IL2Ra = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\alpha$'), "Count"].item()
        IL2Rb = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == 'IL-2R$\\beta$'), "Count"].item()
        gc = data.loc[(data["Cell Type"] == cell_name) & (data["Receptor"] == '$\\gamma_{c}$'), "Count"].item()
        receptors = np.array([IL2Ra, IL2Rb, gc]).astype(np.float)
            
        for _, ligand_name in enumerate(ligand_order):
            
            df = df.append(organize_expr_pred(df, dataMean, cell_name, ligand_name, receptors, muteinC, tps, num), ignore_index=True)

    scales = mutein_scaling(df)
        
    pred_data = np.zeros((12, 4, unkVec_2_15.shape[1]))
    for i, cell_name in enumerate(cell_order):
        for j, ligand_name in enumerate(ligand_order):
            # plot experimental data
            axis = j * 8 + i
            if axis == 31:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name)
                                                                                          & (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, palette=cm.rainbow, legend='full')
                ax[axis].legend(loc='lower right', title="time (hours)")
            else:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name) 
                                                                                          & (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, palette=cm.rainbow, legend=False)
            
            for k, conc in enumerate(dataMean.Concentration.unique()):
                for l, tp in enumerate(tps):
                    print(df.loc[(df["Cells"] == cell_name) & (df["Ligand"] == ligand_name) & (df["Activity Type"] == 'predicted') & (df["Concentration"] == conc) & df["Time Point"] == tp])
                    pred_data[k, l, :] = np.array(df.loc[(df["Cells"] == cell_name) & (df["Ligand"] == ligand_name) & (df["Activity Type"] == 'predicted') & (df["Concentration"] == conc) & df["Time Point"] == tp])

            if cell_name in ['T-reg', 'Mem Treg', 'Naive Treg']:
                plot_dose_response(ax[axis], scales[0, 1] * pred_data / (pred_data + scales[0, 0]), tps, muteinC)
            if cell_name in ['T-helper', 'Mem Th', 'Naive Th']:
                plot_dose_response(ax[axis], scales[1, 1] * pred_data / (pred_data + scales[1, 0]), tps, muteinC)
            if cell_name in ['NK', 'CD8+']:
                plot_dose_response(ax[axis], scales[2, 1] * pred_data / (pred_data + scales[2, 0]), tps, muteinC)
            ax[axis].set(xlabel=("[" + ligand_name + "] (log$_{10}$[nM])"), ylabel="Activity", title=cell_name)

    return f




def organize_expr_pred(df, pstat_df, cell_name, ligand_name, receptors, muteinC, tps, num):
    """ Creates dataframe with experimental and predicted activity for a given cell type and mutein. """
    
    # organize experimental pstat data
    exp_data = np.zeros((12, 4))
    mutein_conc = exp_data.copy()
    for i, conc in enumerate(pstat_df.Concentration.unique()):
        exp_data[i, :] = np.array(pstat_df.loc[(pstat_df["Cells"] == cell_name) & (pstat_df["Ligand"] == ligand_name) & (pstat_df["Concentration"] == conc), "RFU"])
        mutein_conc[i, :] = conc
    
    df_exp = pd.DataFrame({'Cells':np.tile(np.array(cell_name), num), 'Ligand':np.tile(np.array(ligand_name), num), 'Time Point':np.tile(np.array([0.5, 1.0, 2.0, 4.0]), 12), 'Concentration':mutein_conc.reshape(num,), 'Activity Type':np.tile(np.array('experimental'), num), 'Replicate':np.zeros((num,)), 'Activity':exp_data.reshape(num,)})
    df = df.append(df_exp, ignore_index=True)
    
    # calculate predicted dose response
    pred_data = np.zeros((12, 4, unkVec_2_15.shape[1]))
    avg_pred_data = np.zeros((12, 4))
    for j in range(unkVec_2_15.shape[1]):
        cell_receptors = receptor_expression(receptors, unkVec_2_15[17, j], unkVec_2_15[20, j], unkVec_2_15[19, j], unkVec_2_15[21, j])
        pred_data[:, :, j] = calc_dose_response_mutein(unkVec_2_15[:, j], [1., 1., 5.], tps, muteinC, cell_receptors, exp_data)
        df_pred = pd.DataFrame({'Cells':np.tile(np.array(cell_name), num), 'Ligand':np.tile(np.array(ligand_name), num), 'Time Point':np.tile(np.array([0.5, 1.0, 2.0, 4.0]), 12), 'Concentration':mutein_conc.reshape(num,), 'Activity Type':np.tile(np.array('predicted'), num), 'Replicate':np.tile(np.array(j+1), num), 'Activity':pred_data[:, :, j].reshape(num,)})
        df = df.append(df_pred, ignore_index=True)
    
    avg_pred_data = np.mean(pred_data, axis=2)
    df_avg = pd.DataFrame({'Cells':np.tile(np.array(cell_name), num), 'Ligand':np.tile(np.array(ligand_name), num), 'Time Point':np.tile(np.array([0.5, 1.0, 2.0, 4.0]), 12), 'Concentration':mutein_conc.reshape(num,), 'Activity Type':np.tile(np.array('average'), num), 'Replicate':np.tile(np.array(1000), num), 'Activity':avg_pred_data[:, :].reshape(num,)})
    df = df.append(df_avg, ignore_index=True)
    
    return df
    

def mutein_scaling(df):
    cell_groups = [['T-reg', 'Mem Treg', 'Naive Treg'], ['T-helper', 'Mem Th', 'Naive Th'], ['NK', 'CD8+']]
    #df_scaled = pd.DataFrame(columns=['Cells', 'Ligand', 'Time Point', 'Concentration', 'Activity Type', 'Replicate', 'Activity'])
    scales = np.zeros((3, 2))
    for i, cells in enumerate(cell_groups):
        subset_df = df[df['Cells'].isin(cells)]
        scales[i, :] = optimize_scale(np.array(subset_df.loc[(subset_df["Activity Type"] == 'average'), "Activity"]), np.array(subset_df.loc[(subset_df['Activity Type'] == 'experimental'), "Activity"]))
        #subset_df['Activity'] = scales[i, 0] * subset_df['Activity'] / (subset_df['Activity'] + scales[i, 1])
        #df_scaled = df_scaled.append(subset_df, ignore_index=True)
        
    return scales


def calc_dose_response_mutein(unkVec, input_params, tps, muteinC, cell_receptors, exp_data):
    """ Calculates activity for a given cell type at various mutein concentrations and timepoints. """

    total_activity = np.zeros((len(muteinC), len(tps)))

    # loop for each mutein concentration
    for i, conc in enumerate(muteinC):
        active_ckine = runIL2simple(unkVec, input_params, conc, tps=tps, input_receptors=cell_receptors, adj_receptors=True)
        total_activity[i, :] = np.reshape(active_ckine, (-1, 4))  # save the activity from this concentration for all 4 tps

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
