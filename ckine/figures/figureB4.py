"""
This creates Figure 4.
"""
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ckine.imports import import_muteins
from ckine.tensor import perform_decomposition, find_R2X, z_score_values
from ckine.figures.figureCommon import plot_cells
from .figureCommon import subplotLabel, getSetup


def makeFigure():
    """ Build the figure. """
    # Get list of axis objects
    ax, f = getSetup((7, 6), (2, 2))

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    dataMean, dataTensor = import_muteins()
    dataMean.reset_index(inplace=True)

    concs = dataMean['Concentration'].unique()
    ts = dataMean['Time'].unique()
    cells = dataMean['Cells'].unique()
    ligs = dataMean['Ligand'].unique()

    dataTensor = z_score_values(dataTensor, 0)
    parafac = perform_decomposition(dataTensor, 2, weightFactor=3)
    logging.info(find_R2X(dataTensor, parafac))
    
    # Cells
    plot_cells(ax[0], parafac[1], 1, 2, ligs)
    
    # Ligands
    plot_cells(ax[1], parafac[1], 1, 2, ligs)
    
    # Timepoints
    ax[2].plot(ts, parafac[2]);
    
    # Concentration
    ax[3].semilogx(concs, parafac[3]);

    return f
