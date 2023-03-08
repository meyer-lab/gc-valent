"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
from copy import copy
import seaborn as sns
import matplotlib.pyplot as plt
from .figureCommon import subplotLabel, getSetup, getLigDict
from ..imports import import_pstat_all
from ..tensorFac import makeTensor, factorTensor, R2Xplot, plot_tFac_Ligs, plot_tFac_Time, plot_tFac_Conc, plot_tFac_Cells, facScatterPlot

path_here = os.path.dirname(os.path.dirname(__file__))
plt.rcParams['svg.fonttype'] = 'none'
ligDict = getLigDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((15, 6), (2, 5), multz={0: 2})
    axlabel = copy(ax)
    del axlabel[2]
    subplotLabel(axlabel)
    ax[0].axis("off")
    ax[2].axis("off")
    numComps = 3

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True)
    respTensor = makeTensor(respDF)
    tFacAllM = factorTensor(respTensor, numComps=numComps)
    tFacAllM.normalize()

    R2Xplot(ax[1], respTensor, 5)
    ligCompDF = plot_tFac_Ligs(ax[3], tFacAllM, respDF, numComps=numComps)
    plot_tFac_Conc(ax[4], tFacAllM, respDF)
    plot_tFac_Cells(ax[5], tFacAllM, meyer_data=False, numComps=numComps)
    plot_tFac_Time(ax[6], tFacAllM, respDF)
    facScatterPlot(ax[7], ligCompDF)

    legend = ax[3].get_legend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[2].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 10})  # use this to place universal legend later
    ax[3].get_legend().remove()
    ax[5].get_legend().remove()

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
