"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
from .figureCommon import subplotLabel, getSetup
from ..imports import import_pstat_all
from ..tensorFac import makeTensor, factorTensor, R2Xplot, plot_tFac_Ligs, plot_tFac_Time, plot_tFac_Conc, plot_tFac_Cells

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((12, 6), (2, 4), multz={0: 2})
    subplotLabel(ax)
    ax[0].axis("off")
    numComps = 3

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True)
    respTensor = makeTensor(respDF)
    tFacAllM = factorTensor(respTensor, numComps=numComps)
    tFacAllM.normalize()

    R2Xplot(ax[1], respTensor, 5)
    plot_tFac_Ligs(ax[2], tFacAllM, respDF, numComps=numComps)
    plot_tFac_Conc(ax[3], tFacAllM, respDF)
    plot_tFac_Cells(ax[4], tFacAllM, respDF, numComps=numComps)
    plot_tFac_Time(ax[5], tFacAllM, respDF)

    return f
