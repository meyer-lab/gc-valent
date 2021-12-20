"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
from os.path import dirname, join
from .figureCommon import subplotLabel, getSetup
from ..imports import import_pstat_all
from scipy.optimize import least_squares

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((12, 6), (2, 4))
    subplotLabel(ax)
    ax[0].axis("off")
    ax[1].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ  $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    respDF = respDF.loc[respDF.Ligand != "IL15"]
    mutAffDF = mutAffDF.loc[mutAffDF.Ligand != "IL15"]

    affPlot(ax[2], respDF, mutAffDF)
    legend = ax[2].get_legend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[1].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 10})  # use this to place universal legend later
    ax[2].get_legend().remove()
    dosePlot(ax[3], respDF, 1, "Treg")
    dosePlot(ax[4], respDF, 1, "Thelper")
    dosePlot(ax[5], respDF, 1, "CD8")
    dosePlot(ax[6], respDF, 1, "NK")

    return f


def affPlot(ax, respDF, mutAffDF):
    """Plots the various affinities for IL-2 Muteins"""
    mutAffDF["Valency"] = 0
    for ligand in respDF.Ligand.unique():
        valencies = respDF.loc[respDF.Ligand == ligand].Valency.unique()
        if valencies.size == 2:
            mutAffDF.loc[mutAffDF.Ligand == ligand, "Valency"] = "Monovalent and Bivalent"
        elif valencies == 1:
            mutAffDF.loc[mutAffDF.Ligand == ligand, "Valency"] = "Monovalent"
        elif valencies == 2:
            mutAffDF.loc[mutAffDF.Ligand == ligand, "Valency"] = "Bivalent"
    sns.scatterplot(data=mutAffDF, x="IL2Rα $K_{D}$ (nM)", y="IL2Rβ  $K_{D}$ (nM)", hue="Ligand", style="Valency", ax=ax)


def dosePlot(ax, respDF, time, cell):
    """Plots the various affinities for IL-2 Muteins"""
    doses = np.log10(np.logspace(np.log10(respDF.Dose.min()), np.log10(respDF.Dose.max()), 100)) + 4
    x0 = [4, 1, 2]
    hillDF = pd.DataFrame()
    respDF = respDF.loc[(respDF.Cell == cell) & (respDF.Time == time)]
    for ligand in respDF.Ligand.unique():
        for valency in respDF.loc[respDF.Ligand == ligand].Valency.unique():
            isoData = respDF.loc[(respDF.Ligand == ligand) & (respDF.Valency == valency)]
            xData = np.nan_to_num(np.log10(isoData.Dose.values)) + 4
            yData = np.nan_to_num(isoData.Mean.values)
            fit = least_squares(hill_residuals, x0, args=(xData, yData), bounds=([0.0, 0.0, 2], [5, 10.0, 6]), jac="3-point")
            hillDF = hillDF.append(pd.DataFrame({"Ligand": ligand, "Valency": valency, "Dose": np.power(10, doses - 4), "pSTAT": hill_equation(fit.x, doses)}))
    hillDF = hillDF.groupby(["Ligand", "Valency", "Dose"]).pSTAT.mean().reset_index()
    respDF = respDF.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    sns.lineplot(data=hillDF, x="Dose", y="pSTAT", hue="Ligand", style="Valency", ax=ax)

    sns.scatterplot(data=respDF, x="Dose", y="Mean", hue="Ligand", style="Valency", ax=ax, legend=False)
    ax.set(xscale="Log")


def hill_equation(x, dose):
    """ Calculates EC50 from Hill Equation. """
    #print(x, dose)
    EMax = np.power(10, x[0])
    n = x[1]
    EC50 = x[2]
    return EMax * np.power(dose, n) / (np.power(EC50, n) + np.power(dose, n))


def hill_residuals(x, dose, y):
    """ Residual function for Hill Equation. """
    return hill_equation(x, dose) - y
