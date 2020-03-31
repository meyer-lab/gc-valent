"""
This creates Figure 1 for Single Cell FC data analysis. Examples of PCA loadings/scores plots and comparisons to gating.
"""

import os
from .figureCommon import subplotLabel, getSetup
from ..flow import importF
from ..PCA import EC50_PC_Scan, loadingPlot

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((12.5, 5), (2, 5))
    PCscanVecT = [-2, 2, 6]
    PCscanVecNk = [-1, 1, 3]
    loadingT = []
    loadingNk = []

    subplotLabel(ax)

    gates = [False, 'treg', 'nonTreg']
    Titles = ["Tcells", "T-regs", "T Helper"]
    Tsample, _ = importF(path_here + "/data/flow/2019-04-18 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate - NEW PBMC LOT/", "A")
    Nksample, _ = importF(path_here + "/data/flow/2019-03-15 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate/", "A")

    for i, cell in enumerate(gates):
        EC50_PC_Scan(Tsample, PCscanVecT, ax[2 * i], cell, Tcells=True, PC1=True)
        ax[2 * i].set_title(Titles[i] + " PC1 Scan")
        loadingT = EC50_PC_Scan(Tsample, PCscanVecT, ax[2 * i + 1], cell, Tcells=True, PC1=False)
        ax[2 * i + 1].set_title(Titles[i] + " PC2 Scan")
        loadingPlot(loadingT, ax=ax[6], Tcells=True)
        ax[6].set_title(Titles[i] + " Loadings")

    EC50_PC_Scan(Nksample, PCscanVecNk, ax[7], 'nk', Tcells=False, PC1=True)
    ax[7].set_title("Nk PC1 Scan")
    loadingNk = EC50_PC_Scan(Nksample, PCscanVecNk, ax[8], 'nk', Tcells=False, PC1=False)
    ax[8].set_title("Nk PC2 Scan")
    loadingPlot(loadingNk, ax=ax[9], Tcells=False)
    ax[9].set_title("Nk Loadings")

    return f
