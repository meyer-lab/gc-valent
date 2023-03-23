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
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 2), (1, 5))
    subplotLabel(ax)

    modelDF = runFullModelMeyer().reset_index()

    ligandPlot(modelDF, "Treg", ax[0], live_dead=True)
    ligandPlot(modelDF, "Thelper", ax[1], live_dead=True)
    ligandPlot(modelDF, "NK", ax[2], live_dead=True)
    ligandPlot(modelDF, "NKBright", ax[3], live_dead=True)
    ligandPlot(modelDF, "CD8", ax[4], live_dead=True)

    

    return f
