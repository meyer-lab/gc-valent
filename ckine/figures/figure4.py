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
    
    all_cells(ax[0], data, cell_names, unkVec)

    f.tight_layout()

    return f

def single_cell_act(unkVec):
    """ Cytokine activity for all IL2 doses for single cell line. """
    pstat5 = pstat()
    PTS = 30
    cytokC = np.logspace(-3.3, 2.7, PTS)
    act = np.zeros((PTS, unkVec.shape[1]))
    for ii in range(unkVec.shape[1]):
        act[:, ii] = np.fromiter((pstat5.singleCalc(unkVec[:, ii], 0, x) for x in cytokC), np.float64)
    act_mean = np.mean(act, axis=1)
    return act_mean
    
    
def all_cells(ax, cell_data, cell_names, unkVec):
    """ Loops through all cell types and calculates activities. """
    cell_data = cell_data.values    # convert to numpy array
    for ii in range(1, cell_data.shape[1]):
        for n in range(unkVec.shape[1]):
            unkVec[22:30, n] = cell_data[:, ii] # place cell data in all rows of unkVec
        print(single_cell_act(unkVec))
        

