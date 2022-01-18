"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
from os.path import join
import pandas as pd
import seaborn as sns
from .figureCommon import subplotLabel, getSetup, getLigDict
from ..imports import import_pstat_all
from ..tensorFac import makeTensor, factorTensor, R2Xplot, plot_tFac_Ligs, plot_tFac_Time, plot_tFac_Conc, plot_tFac_Cells, swarmPlotFac

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((15, 6), (2, 5), multz={0: 2})
    subplotLabel(ax)
    ax[0].axis("off")
    numComps = 3

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True)
    respTensor = makeTensor(respDF)
    tFacAllM = factorTensor(respTensor, numComps=numComps)
    tFacAllM.normalize()

    R2Xplot(ax[1], respTensor, 5)
    ligCompDF = plot_tFac_Ligs(ax[2], tFacAllM, respDF, numComps=numComps)
    plot_tFac_Conc(ax[3], tFacAllM, respDF)
    plot_tFac_Cells(ax[4], tFacAllM, respDF, numComps=numComps)
    plot_tFac_Time(ax[5], tFacAllM, respDF)
    swarmPlotFac(ax[7], ligCompDF)

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ $K_{D}$ (nM)"}, axis=1)
    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]
    ligCompDF = ligCompDF.loc[(ligCompDF["Lig Name"] != "IL15") & (ligCompDF["Lig Name"] != "IL2")]
    mutAffDF = mutAffDF.rename({"Ligand": "Lig Name"}, axis=1)
    affCompPlot(ax[6], ligCompDF, mutAffDF, "IL2Rα $K_{D}$ (nM)")

    return f


def affCompPlot(ax, ligCompDF, mutAffDF, affinity):
    """Plots component values for a given affinity"""
    affCompDF = ligCompDF.merge(mutAffDF)
    sns.scatterplot(data=affCompDF, x=affinity, y="Component_Val", style="Valency", hue="Component", palette=sns.color_palette(n_colors=3), ax=ax)
    if affinity == "IL2Rα $K_{D}$ (nM)":
        ax.set(xscale="log", xlim=(1e-1, 1e1), ylim=(0, 0.6), ylabel="Component Weight")
    elif affinity == "IL2Rβ $K_{D}$ (nM)":
        ax.set(xscale="log", xlim=(1e0, 30), ylim=(0, 0.6), ylabel="Component Weight")
        ax.get_legend().remove()
