"""
This creates Figure 4.
"""
import string
import numpy as np
import os
import pandas as pds
import matplotlib.cm as cm
from .figureCommon import subplotLabel, getSetup, load_cells, import_samples_2_15
from ..plot_model_prediction import pstat


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4), mults=[0], multz={0: 3})

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    data, cell_names = load_cells()
    unkVec = import_samples_2_15()
    
    all_cells(ax[0], data, cell_names, unkVec[:, 0])

    f.tight_layout()

    return f

def single_cell_act(unkVec, cytokC):
    """ Cytokine activity for all IL2 doses for single cell line. """
    pstat5 = pstat()
    act = np.zeros((cytokC.shape[0]))
    act = np.fromiter((pstat5.singleCalc(unkVec, 0, x) for x in cytokC), np.float64)
    return act / np.max(act)    # normalize to maximal activity
    
    
def all_cells(ax, cell_data, cell_names, unkVec):
    """ Loops through all cell types and calculates activities. """
    cell_data = cell_data.values    # convert to numpy array
    PTS = 60    # number of cytokine concentrations that are used
    cytokC = np.logspace(-5, 0, PTS)
    numCells = cell_data.shape[1] - 1   # first column is receptor names
    results = np.zeros((PTS, numCells))
    
    colors = cm.rainbow(np.linspace(0, 1, numCells))
    
    for ii in range(0, numCells):       # for all cell types
        unkVec[22:30] = cell_data[:, ii+1]  # place cell data into unkVec
        act = single_cell_act(unkVec, cytokC)
        ax.plot(np.log10(cytokC), act, label=cell_names[ii], c=colors[ii])
        if (act[0] > 0.1):
            print(cell_names[ii]) # tells us that proB_FrBC_BM and T_DP_Th cells respond at the lowest IL2 conc.
        
    ax.legend(loc='upper left', bbox_to_anchor=(1.6, 1.2))
