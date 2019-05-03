"""
This creates Figure S7. Full panel of measured vs simulated for IL15.
"""
import string
import numpy as np
import matplotlib.cm as cm
from .figureCommon import subplotLabel, getSetup, import_samples_2_15, plot_conf_int, import_Rexpr, receptor_expression, import_pstat, plot_scaled_pstat
from ..model import runCkineUP, getTotalActiveSpecies


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 4))

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    data_Visterra, cell_names_Visterra = import_Rexpr()
    unkVec_2_15, _ = import_samples_2_15()
    ckineConc, _, _, IL15_data = import_pstat()

        if i == (data_Visterra.shape[0] - 1): # only plot the legend for the last entry
            IL15_dose_response(ax[i], unkVec_2_15, cell_names_Visterra[i], data_Visterra[i], ckineConc, legend=True)
        else:
            IL15_dose_response(ax[i], unkVec_2_15, cell_names_Visterra[i], data_Visterra[i], ckineConc)

    match_cells = [6, 7, 12, 0, 1, 9, 2, 3, 10]

    for j, axis in enumerate(match_cells):
        plot_scaled_pstat(ax[axis], np.log10(ckineConc.astype(np.float)), IL15_data[(j*4):((j+1)*4)])

    f.tight_layout(w_pad=0.1, h_pad=1.0)

    return f

def IL15_dose_response(ax, unkVec, cell_type, cell_data, cytokC, legend=False):
    """ Shows activity for a given cell type at various IL15 concentrations """
    tps = np.array([0.5, 1., 2., 4.])
    PTS = 12 # number of cytokine concentrations
    #cytokC = np.logspace(-4.0, 2.0, PTS) # vary cytokine concentration from 1 pm to 100 nm
    colors = cm.rainbow(np.linspace(0, 1, tps.size))

    rxntfr = unkVec.T.copy()
    split = rxntfr.shape[0] # number of parameter sets used (& thus the number of yOut replicates)
    total_activity = np.zeros((PTS, split, tps.size))

     # loop for each IL15 concentration
    for i in range(PTS):
        for ii in range(rxntfr.shape[0]):
            rxntfr[ii, 0] = cytokC[i]
            # updates rxntfr for receptor expression for IL2Rb, gc, IL15Ra
            rxntfr[ii, 23] = receptor_expression(cell_data[1], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
            rxntfr[ii, 24] = receptor_expression(cell_data[2], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
            rxntfr[ii, 25] = receptor_expression(cell_data[3], rxntfr[ii, 17], rxntfr[ii, 20], rxntfr[ii, 19], rxntfr[ii, 21])
        yOut, retVal = runCkineUP(tps, rxntfr)
        assert retVal >= 0 # make sure solver is working
        activity = np.dot(yOut, getTotalActiveSpecies().astype(np.float))
        for j in range(split):
            total_activity[i, j, :] = activity[(4*j):((j+1)*4)] # save the activity from this concentration for all 4 tps

    # calculate total activity for a given cell type (across all IL15 concentrations & time points)
    avg_total_activity = np.sum(total_activity)/(split*tps.size)

    # plot the values with each time as a separate color
    for tt in range(tps.size):
        plot_conf_int(ax, np.log10(cytokC.astype(np.float)), total_activity[:, :, tt]/avg_total_activity, colors[tt], (tps[tt]).astype(str))

    # plots for input cell type
    ax.set(xlabel=r'IL-15 concentration (log$_{10}$[nM])', ylabel='Activity', title=cell_type)
    if legend is True:
        ax.legend(title='time (min)', loc='center left', borderaxespad=10.)
