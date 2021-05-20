"""
This creates Figure S2, dose response of all IL-2 cytokines using binding model.
"""

import numpy as np
from .figureCommon import getSetup, plotDoseResponses
from ..MBmodel import runFullModel


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    modelDF = runFullModel(time=[0.5, 1], saveDict=False)

    ligands = modelDF.Ligand.unique()
    cells = ["Treg", "Thelper", "NK", "CD8"]
    ax, f = getSetup((10, 16), (len(ligands), 4))

    for i, lig in enumerate(ligands):
        for j, cell in enumerate(cells):
            plotDoseResponses(ax[4 * i + j], modelDF, lig, cell)

    return f
