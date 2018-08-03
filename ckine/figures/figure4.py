"""
This creates Figure 4.
"""
import string
import numpy as np
import os
import pandas as pds
from .figureCommon import subplotLabel, getSetup, load_cells, import_samples_2_15
from ..plot_model_prediction import pstat


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

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
    return act
    
    
def all_cells(ax, cell_data, cell_names, unkVec):
    """ Loops through all cell types and calculates activities. """
    cell_data = cell_data.values    # convert to numpy array
    PTS = 30    # number of cytokine concentrations that are used
    cytokC = np.logspace(-3.3, 2.7, PTS)
    numCells = cell_data.shape[1] - 1   # first column is receptor names
    results = np.zeros((PTS, numCells))
    
    # for ii in range(0, numCells):       # for all cell types
    for ii in range(0, 2):      # just testing first two cell types for now
        unkVec[22:30] = cell_data[:, ii+1] # place cell data in all rows of unkVec
        ax.plot(np.log10(cytokC), single_cell_act(unkVec, cytokC), label=cell_names[ii])
        print("done with " + str(ii) + " cell type")
        
    ax.legend()
    # plot: x-axis is cytokC, y_axis is activity, different lines are cell types
    #ax.plot(cytokC, results[:, 0], title=cell_names[0])    

