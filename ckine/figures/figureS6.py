"""
This creates Figure S7, assymmetric results.
"""

import os
import matplotlib.pyplot as plt
from .figureCommon import subplotLabel, getSetup, getLigDict, ligandPlot
from ..MBmodel import runFullModelMeyer

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
plt.rcParams['svg.fonttype'] = 'none'


def makeFigure():
    """Get a list of the axis objects an@d create a figure"""

    ax, f = getSetup((8, 2), (1, 4))
    subplotLabel(ax)

    modelDF = runFullModelMeyer().reset_index()

    ligandPlot(modelDF, "Treg", ax[0], live_dead=False)
    ligandPlot(modelDF, "Thelper", ax[1], live_dead=False)
    ligandPlot(modelDF, "NK", ax[2], live_dead=False)
    ligandPlot(modelDF, "CD8", ax[3], live_dead=False)

    return f
