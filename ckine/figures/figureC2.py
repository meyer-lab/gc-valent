"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.optimize import least_squares
from os.path import dirname, join
from .figureCommon import subplotLabel, getSetup, getLigDict, get_cellTypeDict
from ..imports import import_pstat_all

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((12, 12), (4, 4))
    subplotLabel(ax)
    ax[0].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ  $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    respDF = respDF.loc[(respDF.Ligand != "IL15") & (respDF.Ligand != "IL2")]
    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]

    doses = respDF.Dose.unique()
    dose = doses[3]
    ratioConc(ax[4:7], respDF, r"T$_{reg}$", "NK", 4, mutAffDF)
    ratioConc(ax[8:11], respDF, r"T$_{reg}$", r"T$_{helper}$", 4, mutAffDF)
    ratioConc(ax[12:15], respDF, r"T$_{reg}$", r"CD8$^{+}$", 4, mutAffDF)

    pSTATcomp(ax[1], respDF, mutAffDF, "IL2Rα $K_{D}$ (nM)", r"T$_{reg}$", 4, dose, legend=True)
    pSTATcomp(ax[2], respDF, mutAffDF, "IL2Rα $K_{D}$ (nM)", "NK", 4, dose)
    pSTATcomp(ax[3], respDF, mutAffDF, "IL2Rβ  $K_{D}$ (nM)", r"CD8$^{+}$", 4, dose)
    """
    pSTATcomp(ax[4], respDF, mutAffDF, "IL2Rβ  $K_{D}$ (nM)", "NK", 1, dose)
    pSTATcomp(ax[5], respDF, mutAffDF, "IL2Rα $K_{D}$ (nM)", r"CD8$^{+}$", 1, dose)
    pSTATcomp(ax[6], respDF, mutAffDF, "IL2Rβ  $K_{D}$ (nM)", r"CD8$^{+}$", 1, dose)
    """

    legend = ax[1].get_legend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[0].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 8})  # use this to place universal legend later
    ax[1].get_legend().remove()

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


def ratioConc(ax, respDF, cell1, cell2, time, mutAffDF):
    """Plots Ratio of cell 1 to cell 2 over a range of concentrations"""
    respDF1 = respDF.loc[(respDF.Cell == cell1) & ((respDF.Time == time))]
    respDF2 = respDF.loc[(respDF.Cell == cell2) & ((respDF.Time == time))]
    respDF1 = respDF1.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    respDF2 = respDF2.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    respDF1 = respDF1.pivot(index=["Ligand", "Dose", "Valency"], columns="Cell", values="Mean").reset_index()
    respDF2 = respDF2.pivot(index=["Ligand", "Dose", "Valency"], columns="Cell", values="Mean").reset_index()
    ratioDF = respDF1.merge(respDF2)
    ratioDF["Ratio"] = ratioDF[cell1] / ratioDF[cell2]
    ratioDF = ratioDF.dropna()

    ratioDF.loc[ratioDF.Dose < np.sort(ratioDF.Dose.unique())[3], "Ratio"] = 1
    ratioDF = ratioDF.loc[(ratioDF.Ratio <= 3000) & (ratioDF.Ratio > 0)]
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
            if xData.size > 5:
                fit = least_squares(gaussian_residuals, x0, args=(xData, yData), bounds=([0, -4, 0], [3.5, 2, 6]), jac="3-point")
                gaussDF = gaussDF.append(pd.DataFrame({"Ligand": ligand, "Valency": valency, "Dose": np.power(10, doses), "Ratio": np.power(10, gaussian(fit.x, doses))}))
                fitDF = fitDF.append(pd.DataFrame({"Ligand": [ligand], "Valency": valency, cell2 + " Max": np.power(10, fit.x[0]), cell2 + " Dose": np.power(10, fit.x[1])}))

    gaussDF = gaussDF.reset_index()
    ratioDF["Ratio"] = np.power(10, ratioDF["Ratio"].values)

    sns.scatterplot(data=ratioDF, x="Dose", y="Ratio", hue="Ligand", style="Valency", ax=ax[0], palette=ligDict, legend=False)
    sns.lineplot(data=gaussDF, x="Dose", y="Ratio", hue="Ligand", style="Valency", ax=ax[0], palette=ligDict, legend=False)
    ax[0].set(xscale="log", title="Ratio of " + cell1 + " to " + cell2)
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

    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", hue="Ligand", style="Valency", ax=ax[1], palette=ligDict, legend=False)
    sns.lineplot(data=maxLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", style="Valency", ax=ax[1], color="k", legend=False)
    ax[1].set(xscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(0, None))

    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", hue="Ligand", style="Valency", ax=ax[2], palette=ligDict, legend=False)
    sns.lineplot(data=doseLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", style="Valency", ax=ax[2], color="k", legend=False)
    ax[2].set(xscale="log", yscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(1e-2, 1e2))
