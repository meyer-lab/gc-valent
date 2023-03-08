"""
This creates Figure S4, dose response of all IL-2 cytokines using binding model.
"""

import matplotlib.pyplot as plt
from .figureCommon import getSetup, plotDoseResponses
from ..MBmodel import runFullModel

plt.rcParams['svg.fonttype'] = 'none'


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    modelDF = runFullModel(time=[0.5, 1], saveDict=False)

    ligands = modelDF.Ligand.unique()
    cells = ["Treg", "Thelper", "NK", "NKBright", "CD8"]
    ax, f = getSetup((12.5, 16), (len(ligands), 5))

    for i, lig in enumerate(ligands):
        for j, cell in enumerate(cells):
            plotDoseResponses(ax[5 * i + j], modelDF, lig, cell)

    return f
