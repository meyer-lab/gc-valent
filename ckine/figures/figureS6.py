"""
This creates Figure S6. Full panel of measured vs simulated for IL2.
"""
import string
import numpy as np
import matplotlib.cm as cm
from scipy.optimize import minimize
from .figureCommon import subplotLabel, getSetup, plot_conf_int, plot_scaled_pstat
from ..model import runCkineUP, getTotalActiveSpecies, receptor_expression
from ..imports import import_Rexpr, import_pstat, import_samples_2_15


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 5))

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    _, receptor_data, cell_names_receptor = import_Rexpr()
    unkVec_2_15, scales = import_samples_2_15(N=100)  # use all rates
    ckineConc, cell_names_pstat, IL2_data, IL15_data = import_pstat()
    axis = 0

    for i, _ in enumerate(cell_names_pstat):
        # plot matching experimental and predictive pSTAT data for the same cell type
        for j in range(receptor_data.shape[0]):
            if cell_names_pstat[i] == cell_names_receptor[j]:
                if axis == 9:  # only plot the legend for the last entry
                    dose_response(ax[axis], unkVec_2_15, scales, cell_names_receptor[j],
                                  receptor_data[j], 0, ckineConc, IL2_data[(i * 4):((i + 1) * 4)])  # IL-2
                    dose_response(ax[axis+10], unkVec_2_15, scales, cell_names_receptor[j],
                                  receptor_data[j], 1, ckineConc, IL15_data[(i * 4):((i + 1) * 4)], legend=True)  # IL-15
                else:
                    dose_response(ax[axis], unkVec_2_15, scales, cell_names_receptor[j],
                                  receptor_data[j], 0, ckineConc, IL2_data[(i * 4):((i + 1) * 4)])  # IL-2
                    dose_response(ax[axis+10], unkVec_2_15, scales, cell_names_receptor[j],
                                  receptor_data[j], 1, ckineConc, IL15_data[(i * 4):((i + 1) * 4)])  # IL-15
                plot_scaled_pstat(ax[axis], np.log10(ckineConc.astype(np.float)), IL2_data[(i * 4):((i + 1) * 4)])
                plot_scaled_pstat(ax[axis+10], np.log10(ckineConc.astype(np.float)), IL15_data[(i * 4):((i + 1) * 4)])
                axis = axis + 1

    f.tight_layout(w_pad=0.1, h_pad=1.0)

    return f


def dose_response(ax, unkVec, scales, cell_type, cell_data, cytokIDX, cytokC, exp_data, legend=False):
    """ Shows activity for a given cell type at various cytokine concentrations and timepoints. """
    tps = np.array([0.5, 1., 2., 4.]) * 60.
    PTS = 12  # number of cytokine concentrations
    colors = cm.rainbow(np.linspace(0, 1, tps.size))

    rxntfr = unkVec.T.copy()
    split = rxntfr.shape[0]  # number of parameter sets used (& thus the number of yOut replicates)
    total_activity = np.zeros((PTS, split, tps.size))

    # loop for each IL2 concentration
    for i in range(PTS):
        for ii in range(rxntfr.shape[0]):
            rxntfr[ii, cytokIDX] = cytokC[i]
            # updates rxntfr for receptor expression for IL2Ra, IL2Rb, gc
            rxntfr[ii, 22] = receptor_expression(cell_data[0], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
            rxntfr[ii, 23] = receptor_expression(cell_data[1], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
            rxntfr[ii, 24] = receptor_expression(cell_data[2], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
        yOut, retVal = runCkineUP(tps, rxntfr)
        assert retVal >= 0  # make sure solver is working
        activity = np.dot(yOut, getTotalActiveSpecies().astype(np.float))
        for j in range(split):
            total_activity[i, j, :] = activity[(4 * j):((j + 1) * 4)]  # save the activity from this concentration for all 4 tps

    for j in range(len(scales)):
        guess = np.array([scales[j, 0], 100.])  # first value is 
        scale1, scale2 = optimize_scale(guess, total_activity[:, j, :], exp_data)
        total_activity[:, j, :] = scale2 * total_activity[:, j, :] / (total_activity[:, j, :] + scale1)  # account for pSTAT5 saturation

    # plot the values with each time as a separate color
    for tt in range(tps.size):
        if legend:
            plot_conf_int(ax, np.log10(cytokC.astype(np.float)), total_activity[:, :, tt], colors[tt], (tps[tt] / 60.).astype(str))
            ax.legend(title='time (hours)', loc='center left', borderaxespad=10.)
        else:
            plot_conf_int(ax, np.log10(cytokC.astype(np.float)), total_activity[:, :, tt], colors[tt])

    # plots for input cell type
    if cytokIDX == 0:
        ax.set(xlabel=r'[IL-2] (log$_{10}$[nM])', ylabel='Activity', title=cell_type)
    elif cytokIDX == 1:
        ax.set(xlabel=r'[IL-15] (log$_{10}$[nM])', ylabel='Activity', title=cell_type)


def optimize_scale(scale_guess, model_act, exp_act):
    """ Formulates the optimal scale to minimize the residual between model activity predictions and experimental activity measurments for a given cell type. """
    exp_act = exp_act.T  # transpose to match model_act

    def calc_res(sc):
        """ Calculate the residual.. This is the function we minimize. """
        scaled_act = sc[1] * model_act / (model_act + sc[0])
        return np.sum(np.square(exp_act - scaled_act))  # return sum of squared error (a scalar)

    res = minimize(calc_res, scale_guess, bounds=((0, None), (0, None)))
    return res.x
