"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from os.path import join
from copy import copy
from .figureCommon import subplotLabel, getSetup, getLigDict, get_cellTypeDict, hillRatioDosePlot
from ..imports import import_pstat_all


path_here = os.path.dirname(os.path.dirname(__file__))
plt.rcParams['svg.fonttype'] = 'none'
ligDict = getLigDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((11.25, 2.5), (1, 5))
    axlabel = copy(ax)
    del axlabel[0]
    subplotLabel(axlabel)
    ax[0].axis("off")
    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    time = 4.0
    pseudo = 0.15

    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]

    maxRatConc(ax[1], respDF, r"T$_{reg}$", "NK", time, mutAffDF, pseudo=pseudo, legend=True)
    maxRatConc(ax[1], respDF, r"T$_{reg}$", r"NK$^{Bright}$", time, mutAffDF, pseudo=pseudo, legend=True)
    maxRatConc(ax[2], respDF, r"T$_{reg}$", r"CD8$^{+}$", time, mutAffDF, pseudo=pseudo, legend=True)
    maxRatConc(ax[3], respDF, r"T$_{reg}$", r"T$_{helper}$", time, mutAffDF, pseudo=pseudo, legend=True)

    return f


def maxRatConc(ax, respDF, cell1, cell2, time, mutAffDF, pseudo=0.1, legend=False):
    """Plots Ratio of cell 1 to cell 2 over a range of concentrations"""
    hillDF = hillRatioDosePlot(ax, respDF, time, cell1, cell2, pseudo=pseudo, plot=False)
    fitDF = pd.DataFrame()
    for ligand in hillDF.Ligand.unique():
        for valency in hillDF.loc[hillDF.Ligand == ligand].Valency.unique():
            isoData = hillDF.loc[(hillDF.Ligand == ligand) & (hillDF.Valency == valency)]
            fitDF = pd.concat([fitDF, pd.DataFrame({"Ligand": [ligand], "Valency": valency, cell2 + " Max": isoData.Ratio.max(),
                                                    cell2 + " Dose": isoData.loc[isoData.Ratio == isoData.Ratio.max()].Dose.values})])

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

    fitDF, doseLineDF = fitDF.reset_index(), doseLineDF.reset_index()
    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", hue="Ligand", style="Valency", ax=ax, palette=ligDict, legend=True)
    sns.lineplot(data=doseLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", style="Valency", ax=ax, color="k", linewidth=1., legend=True)
    ax.set(xscale="log", yscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(1e-2, 1e2), ylabel=cell1 + "/" + cell2 + " Max Dose (nM)")
