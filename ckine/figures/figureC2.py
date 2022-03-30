"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.optimize import least_squares
from os.path import join
from copy import copy
from .figureCommon import subplotLabel, getSetup, getLigDict, get_cellTypeDict, getLigandLegend, Wass_KL_Dist, CITE_RIDGE, CITE_SVM
from ..imports import import_pstat_all


path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((14, 18), (6, 4), multz={0: 2, 4: 1, 8: 1, 12: 1})
    axlabel = copy(ax)
    del axlabel[1]
    del axlabel[11]
    del axlabel[13]
    del axlabel[13]
    subplotLabel(axlabel)
    ax[0].axis("off")
    ax[1].axis("off")
    ax[11].axis("off")
    ax[12].axis("off")
    ax[15].axis("off")
    ax[16].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    respDF = respDF.loc[(respDF.Ligand != "IL15") & (respDF.Ligand != "IL2")]
    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]

    pseudo = 20
    ratioConc(ax[2:5], respDF, r"T$_{reg}$", "NK", 4, mutAffDF, pseudo=pseudo, cutoff=200, legend=True)
    ratioConc(ax[5:8], respDF, r"T$_{reg}$", r"CD8$^{+}$", 4, mutAffDF, pseudo=pseudo, cutoff=250, legend=True)
    ratioConc(ax[8:11], respDF, r"T$_{reg}$", r"T$_{helper}$", 4, mutAffDF, pseudo=pseudo, cutoff=150, legend=True)

    legend = getLigandLegend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[1].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 10})  # use this to place universal legend later
    cellTarget = "Treg"
    Wass_KL_Dist(ax[13:15], cellTarget, 10)
    CITE_RIDGE(ax[17], cellTarget)
    CITE_SVM(ax[18], cellTarget, sampleFrac=0.2)

    return f


def pSTATcomp(ax, respDF, mutAffDF, recAff, cell, time, conc, legend=False):
    """Plots monovalent and bivalent lines for a single cell type by affinity"""
    respDF = respDF.loc[(respDF.Cell == cell) & ((respDF.Time == time) & (respDF.Dose == conc))]
    respDF = respDF.groupby(["Ligand", "Valency"]).Mean.mean().reset_index()
    respDF = respDF.merge(mutAffDF)
    respDF = respDF.loc[(respDF.Ligand != "IL15") & (respDF.Ligand != "IL2")]
    sns.scatterplot(data=respDF, x=recAff, y="Mean", hue="Ligand", style="Valency", ax=ax, legend=legend, palette=ligDict)
    ax.set(xscale="log", title=cell + " at " + str(conc) + " nM  and " + str(time) + " hours")


def gaussian(x, concs):
    """Returns Gaussian for given concs"""
    amplitude, mean, stddev = x[0], x[1], x[2]
    return amplitude * np.exp(-((concs - mean) / 4 / stddev)**2)


def gaussian_residuals(x, concs, ratios):
    """ Residual function for Gaussian. """
    return gaussian(x, concs) - ratios


def ratioConc(ax, respDF, cell1, cell2, time, mutAffDF, pseudo, cutoff, legend=False):
    """Plots Ratio of cell 1 to cell 2 over a range of concentrations"""
    respDF1 = respDF.loc[(respDF.Cell == cell1) & ((respDF.Time == time))]
    respDF2 = respDF.loc[(respDF.Cell == cell2) & ((respDF.Time == time))]
    respDF1 = respDF1.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    respDF2 = respDF2.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    respDF1 = respDF1.pivot(index=["Ligand", "Dose", "Valency"], columns="Cell", values="Mean").reset_index()
    respDF2 = respDF2.pivot(index=["Ligand", "Dose", "Valency"], columns="Cell", values="Mean").reset_index()
    ratioDF = respDF1.merge(respDF2)
    ratioDF["Ratio"] = ratioDF[cell1] / (ratioDF[cell2] + pseudo)
    ratioDF = ratioDF.dropna()

    ratioDF.loc[ratioDF.Dose < np.sort(ratioDF.Dose.unique())[3], "Ratio"] = 1
    ratioDF = ratioDF.loc[(ratioDF.Ratio <= cutoff) & (ratioDF.Ratio > 0) & (ratioDF[cell2] > 0)]
    ratioDF["Ratio"] = np.log10(ratioDF["Ratio"].values)

    doses = np.log10(np.logspace(np.log10(ratioDF.Dose.min()), np.log10(ratioDF.Dose.max()), 100))
    x0 = [3, 0, 2]
    gaussDF = pd.DataFrame()
    fitDF = pd.DataFrame()
    for ligand in ratioDF.Ligand.unique():
        for valency in ratioDF.loc[ratioDF.Ligand == ligand].Valency.unique():
            isoData = ratioDF.loc[(ratioDF.Ligand == ligand) & (ratioDF.Valency == valency)]
            xData = np.nan_to_num(np.log10(isoData.Dose.values))
            yData = np.nan_to_num(isoData.Ratio.values)
            if xData.size > 7:
                fit = least_squares(gaussian_residuals, x0, args=(xData, yData), bounds=([0, -4, 0], [3.5, 2, 6]), jac="3-point")
                gaussDF = gaussDF.append(pd.DataFrame({"Ligand": ligand, "Valency": valency, "Dose": np.power(10, doses), "Ratio": np.power(10, gaussian(fit.x, doses))}))
                fitDF = fitDF.append(pd.DataFrame({"Ligand": [ligand], "Valency": valency, cell2 + " Max": np.power(10, fit.x[0]), cell2 + " Dose": np.power(10, fit.x[1])}))

    gaussDF = gaussDF.reset_index()
    ratioDF["Ratio"] = np.power(10, ratioDF["Ratio"].values)

    sns.scatterplot(data=ratioDF, x="Dose", y="Ratio", hue="Ligand", style="Valency", size="Valency", ax=ax[0], palette=ligDict, legend=False)
    sns.lineplot(data=gaussDF, x="Dose", y="Ratio", hue="Ligand", size="Valency", ax=ax[0], palette=ligDict, legend=legend, sizes=(1, 2.5))
    ax[0].set(xscale="log", title="Ratio of " + cell1 + " to " + cell2)

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
        maxLineDF = maxLineDF.append(pd.DataFrame({"Valency": valency, "IL2Rα $K_{D}$ (nM)": np.power(10, affs), cell2 + " Max": mMax * affs + bMax}))
        doseLineDF = doseLineDF.append(pd.DataFrame({"Valency": valency, "IL2Rα $K_{D}$ (nM)": np.power(10, affs), cell2 + " Dose": np.power(10, mDose * affs + bDose)}))

    maxLineDF, doseLineDF = maxLineDF.reset_index(), doseLineDF.reset_index()
    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", hue="Ligand", style="Valency", ax=ax[1], palette=ligDict, legend=False)
    sns.lineplot(data=maxLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", style="Valency", ax=ax[1], color="k", linewidth=1., legend=False)
    ax[1].set(xscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(0, None), ylabel=cell1 + " to " + cell2 + " Ratio Max Magnitude")

    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", hue="Ligand", style="Valency", ax=ax[2], palette=ligDict, legend=False)
    sns.lineplot(data=doseLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", style="Valency", ax=ax[2], color="k", linewidth=1., legend=False)
    ax[2].set(xscale="log", yscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(1e-2, 1e2), ylabel=cell1 + " to " + cell2 + " Ratio Max Dose")
