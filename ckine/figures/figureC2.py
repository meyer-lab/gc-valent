"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
from copy import copy
from .figureCommon import subplotLabel, getSetup, overlayCartoon
from ..imports import import_pstat_all
from ..tensorFac import makeTensor, factorTensor, R2Xplot, plot_tFac_Ligs, plot_tFac_Time, plot_tFac_Conc, plot_tFac_Cells

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((12, 6), (2, 4), multz={0: 1})
    axlabel = copy(ax)
    del axlabel[6]
    subplotLabel(axlabel)
    ax[0].axis("off")
    ax[6].axis("off")

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all()
    respTensor = makeTensor(respDF)
    tFacAllM = factorTensor(respTensor, 3)
    tFacAllM.normalize()

    R2Xplot(ax[1], respTensor, 5)
    plot_tFac_Ligs(ax[2], tFacAllM, respDF)
    plot_tFac_Conc(ax[3], tFacAllM, respDF)
    plot_tFac_Cells(ax[4], tFacAllM, respDF)
    plot_tFac_Time(ax[5], tFacAllM, respDF)

    return f
