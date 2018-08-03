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

    f.tight_layout()

    return f

def IL2_vs_cells(ax, cell_data, cell_names, unkVec):
    """ Cytokine activity vs IL2 dose for all cell populations. """
    PTS = 30
    cytokC = np.logspace(-3.3, 2.7, PTS)
    act = np.fromiter((pstat.singleCalc(unkVec, 0, x) for x in cytokC), np.float64)

