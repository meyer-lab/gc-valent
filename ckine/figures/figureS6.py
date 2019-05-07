"""
This creates Figure S6. Full panel of measured vs simulated for IL2.
"""
import string
import numpy as np
import matplotlib.cm as cm
from .figureCommon import subplotLabel, getSetup, import_samples_2_15, plot_conf_int, receptor_expression, import_pstat, plot_scaled_pstat
from ..model import runCkineUP, getTotalActiveSpecies
from ..tensor_generation import import_Rexpr


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 4))

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    _, data_Visterra, cell_names_Visterra = import_Rexpr()
    unkVec_2_15, _ = import_samples_2_15()
    ckineConc, _, IL2_data, _ = import_pstat()

    # Delete and organize data to match the cells between receptor level and pSTAT5 data
    data_Visterra = np.delete(data_Visterra, [0, 8], 0)
    cell_names_Visterra.remove('CD3+CD4+')
    cell_names_Visterra.remove('NKT')
    match_cells = [5, 6, 7, 8, 9, 10, 0, 2, 3, 4]

    for i in range(data_Visterra.shape[0]):
        if i == (data_Visterra.shape[0] - 1):  # only plot the legend for the last entry
            IL2_dose_response(ax[i], unkVec_2_15, cell_names_Visterra[i], data_Visterra[i], ckineConc, legend=True)
        else:
            IL2_dose_response(ax[i], unkVec_2_15, cell_names_Visterra[i], data_Visterra[i], ckineConc)

    for axis, j in enumerate(match_cells):
        plot_scaled_pstat(ax[axis], np.log10(ckineConc.astype(np.float)), IL2_data[(j * 4):((j + 1) * 4)])

    f.tight_layout(w_pad=0.1, h_pad=1.0)

    return f


def IL2_dose_response(ax, unkVec, cell_type, cell_data, cytokC, legend=False):
    """ Shows activity for a given cell type at various IL2 concentrations """
    tps = np.array([0.5, 1., 2., 4.])
    PTS = 12  # number of cytokine concentrations
    colors = cm.rainbow(np.linspace(0, 1, tps.size))

    rxntfr = unkVec.T.copy()
    split = rxntfr.shape[0]  # number of parameter sets used (& thus the number of yOut replicates)
    total_activity = np.zeros((PTS, split, tps.size))

    # loop for each IL2 concentration
    for i in range(PTS):
        for ii in range(rxntfr.shape[0]):
            rxntfr[ii, 0] = cytokC[i]
            # updates rxntfr for receptor expression for IL2Ra, IL2Rb, gc
            rxntfr[ii, 22] = receptor_expression(cell_data[0], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
            rxntfr[ii, 23] = receptor_expression(cell_data[1], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
            rxntfr[ii, 24] = receptor_expression(cell_data[2], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
        yOut, retVal = runCkineUP(tps, rxntfr)
        assert retVal >= 0  # make sure solver is working
        activity = np.dot(yOut, getTotalActiveSpecies().astype(np.float))
        for j in range(split):
            total_activity[i, j, :] = activity[(4 * j):((j + 1) * 4)]  # save the activity from this concentration for all 4 tps

    # calculate total activity for a given cell type (across all IL2 concentrations & time points)
    avg_total_activity = np.sum(total_activity) / (split * tps.size)

    # plot the values with each time as a separate color
    for tt in range(tps.size):
        plot_conf_int(ax, np.log10(cytokC.astype(np.float)), total_activity[:, :, tt] / avg_total_activity, colors[tt], (tps[tt]).astype(str))

    # plots for input cell type
    ax.set(xlabel=r'IL-2 concentration (log$_{10}$[nM])', ylabel='Activity', title=cell_type)
    if legend is True:
        ax.legend(title='time (hours)', loc='center left', borderaxespad=10.)
