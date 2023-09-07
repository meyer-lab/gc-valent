"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from os.path import join
from copy import copy
from .figureCommon import subplotLabel, getSetup, getLigDict, get_cellTypeDict, getLigandLegend, Wass_KL_Dist, hillRatioDosePlot
from ..imports import import_pstat_all

path_here = os.path.dirname(os.path.dirname(__file__))
plt.rcParams['svg.fonttype'] = 'none'
ligDict = getLigDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((7.5, 15.1666), (7, 3), multz={0: 1, 3: 1, 6: 1, 9: 1})
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
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    time = 4.0
    pseudo = 0.15

    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]
    ratioConc(ax[2:4], respDF, r"T$_{reg}$", "NK", time, mutAffDF, pseudo=pseudo, legend=True)
    ratioConc(ax[4:6], respDF, r"T$_{reg}$", r"NK$^{Bright}$", time, mutAffDF, pseudo=pseudo, legend=True)
    ratioConc(ax[6:8], respDF, r"T$_{reg}$", r"CD8$^{+}$", time, mutAffDF, pseudo=pseudo, legend=True)
    ratioConc(ax[8:10], respDF, r"T$_{reg}$", r"T$_{helper}$", time, mutAffDF, pseudo=pseudo, legend=True)

    legend = getLigandLegend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[1].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 10})  # use this to place universal legend later
    cellTarget = "Treg"
    Wass_KL_Dist(ax[9:11], cellTarget, 10)
    Wass_KL_Dist(ax[11:13], cellTarget, 10, RNA=True)
    ax[9].set(xlim=(1e-1, 1e1))
    #CITE_RIDGE(ax[13], cellTarget)

    return f


def ratioConc(ax, respDF, cell1, cell2, time, mutAffDF, pseudo=0.1, legend=False):
    """Plots Ratio of cell 1 to cell 2 over a range of concentrations"""
    hillDF = hillRatioDosePlot(ax[0], respDF, time, cell1, cell2, pseudo=pseudo)
    fitDF = pd.DataFrame()
    for ligand in hillDF.Ligand.unique():
        for valency in hillDF.loc[hillDF.Ligand == ligand].Valency.unique():
            isoData = hillDF.loc[(hillDF.Ligand == ligand) & (hillDF.Valency == valency)]
            fitDF = pd.concat([fitDF, pd.DataFrame({"Ligand": [ligand], "Valency": valency, cell2 + " Max": isoData.Ratio.max(),
                                                    cell2 + " Dose": isoData.loc[isoData.Ratio == isoData.Ratio.max()].Dose.values})])

    ax[0].set(title="Ratio of " + cell1 + " to " + cell2)

    if legend:
        h, l = ax[0].get_legend_handles_labels()
        ax[0].legend(h[-3:], l[-3:])

    fitDF = fitDF.merge(mutAffDF)

    maxLineDF = pd.DataFrame()
    doseLineDF = pd.DataFrame()
    affs = np.linspace(-1, 1, 100).flatten()
    for valency in fitDF.Valency.unique():
        valData = fitDF.loc[fitDF.Valency == valency]
        mMax, bMax = np.polyfit(np.log10(valData["IL2Rα $K_{D}$ (nM)"].values), valData[cell2 + " Max"], 1)
        mDose, bDose = np.polyfit(np.log10(valData["IL2Rα $K_{D}$ (nM)"].values), np.log10(valData[cell2 + " Dose"]), 1)
        maxLineDF = pd.concat([maxLineDF, pd.DataFrame({"Valency": valency, "IL2Rα $K_{D}$ (nM)": np.power(10, affs), cell2 + " Max": mMax * affs + bMax})])
        doseLineDF = pd.concat([doseLineDF, pd.DataFrame({"Valency": valency, "IL2Rα $K_{D}$ (nM)": np.power(10, affs), cell2 + " Dose": np.power(10, mDose * affs + bDose)})])

    maxLineDF, doseLineDF = maxLineDF.reset_index(), doseLineDF.reset_index()
    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", hue="Ligand", style="Valency", ax=ax[1], palette=ligDict, legend=False)
    sns.lineplot(data=maxLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", style="Valency", ax=ax[1], color="k", linewidth=1., legend=False)
    ax[1].set(xscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(0, None), ylabel=cell1 + "/" + cell2 + " Max Magnitude")
