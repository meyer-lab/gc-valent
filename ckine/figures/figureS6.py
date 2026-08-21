"""
This creates Figure S7, assymmetric results.
"""

import os
import matplotlib.pyplot as plt
from .figureCommon import subplotLabel, getSetup, getLigDict, ligandPlot, ligand_ratio_plot
from ..MBmodel import runFullModelMeyer

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
plt.rcParams['svg.fonttype'] = 'none'


def makeFigure():
    """Get a list of the axis objects an@d create a figure"""

    ax, f = getSetup((7, 7), (3, 3))
    subplotLabel(ax)

    modelDF = runFullModelMeyer().reset_index()

    ligandPlot(modelDF, "Treg", ax[0], live_dead=False)
    ligandPlot(modelDF, "Thelper", ax[1], live_dead=False)
    ligandPlot(modelDF, "NK", ax[2], live_dead=False)
    ligandPlot(modelDF, "NKBright", ax[3], live_dead=False)
    ligandPlot(modelDF, "CD8", ax[4], live_dead=False)
    ligand_ratio_plot(modelDF, "Treg", "Thelper", ax[5], live_dead=False)
    ligand_ratio_plot(modelDF, "Treg", "NK", ax[6], live_dead=False)
    ligand_ratio_plot(modelDF, "Treg", "NKBright", ax[7], live_dead=False)
    ligand_ratio_plot(modelDF, "Treg", "CD8", ax[8], live_dead=False)

    return f
