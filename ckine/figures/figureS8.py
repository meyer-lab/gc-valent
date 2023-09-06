"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os.path import join
from copy import copy
from .figureCommon import subplotLabel, getSetup, getLigDict, get_cellTypeDict, make_EC50_DF
from ..imports import import_pstat_all, import_pstat_all_meyer


path_here = os.path.dirname(os.path.dirname(__file__))
plt.rcParams['svg.fonttype'] = 'none'
ligDict = getLigDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((7.5, 13), (6, 3), multz={0: 1, 3: 1, 6: 1, 9: 1})
    axlabel = copy(ax)
    del axlabel[1]
    subplotLabel(axlabel)
    ax[0].axis("off")
    ax[1].axis("off")
    ax[8].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1

    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]
    make_EC50_DF(respDF, 0.5)

    respDF = import_pstat_all_meyer()
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency", "pSTAT5": "Mean"}, axis=1)
    make_EC50_DF(respDF, np.nan, meyer=True)



    return f
