"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
from os.path import dirname, join
from .figureCommon import subplotLabel, getSetup
from ..imports import import_pstat_all
from ..tensorFac import makeTensor, factorTensor, R2Xplot, plot_tFac_Ligs, plot_tFac_Time, plot_tFac_Conc, plot_tFac_Cells

path_here = os.path.dirname(os.path.dirname(__file__))

def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((12, 6), (2, 4), multz={0: 2})
    subplotLabel(ax)
    ax[0].axis("off")
    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ  $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    doses = respDF.Dose.unique()
    pSTATcomp(ax[1], respDF, mutAffDF, "IL2Rα $K_{D}$ (nM)", "Treg", 1, doses[7])

    return f



def pSTATcomp(ax, respDF, mutAffDF, recAff, cell, time, conc):
    """Plots monovalent and bivalent lines for a single cell type by affinity"""
    respDF = respDF.loc[(respDF.Cell == cell) & ((respDF.Time == time) & (respDF.Dose == conc))]
    print(respDF)
    respDF = respDF.merge(mutAffDF)
    respDF = respDF.loc[respDF.Ligand != "IL15"]
    print(respDF)
    sns.scatterplot(data=respDF, x=recAff, y="Mean", hue="Ligand", style="Bivalent", ax=ax)
    ax.set(xscale="log")