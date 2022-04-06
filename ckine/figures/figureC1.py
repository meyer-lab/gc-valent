"""
This creates Figure 1. Initial look at flow data.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
from copy import copy
from scipy.optimize import least_squares
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from os.path import dirname, join
from .figureCommon import subplotLabel, getSetup, getLigDict, get_doseLimDict, get_cellTypeDict
from ..imports import import_pstat_all

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
limDict = get_doseLimDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((15, 18), (6, 4), multz={0: 1, 4: 7})
    axlabel = copy(ax)
    del axlabel[1]
    subplotLabel(axlabel)
    ax[0].axis("off")
    ax[1].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ  $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    respDF = respDF.loc[respDF.Ligand != "IL15"]
    mutAffDF = mutAffDF.loc[mutAffDF.Ligand != "IL15"]

    affPlot(ax[2], respDF, mutAffDF)
    legend = ax[2].get_legend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[1].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 10})  # use this to place universal legend later
    ax[2].get_legend().remove()
    fullHeatMap(ax[3], respDF)
    dosePlot(ax[4], respDF, 1, r"T$_{reg}$", ligList=["IL2", "R38Q N-term", "N88D C-term", "WT N-term"], legend=True)
    dosePlot(ax[5], respDF, 1, r"T$_{helper}$", ligList=["IL2", "R38Q N-term", "N88D C-term", "WT N-term"])
    dosePlot(ax[6], respDF, 1, r"CD8$^{+}$", ligList=["IL2", "R38Q N-term", "N88D C-term", "WT N-term"])
    dosePlot(ax[7], respDF, 1, "NK", ligList=["IL2", "R38Q N-term", "N88D C-term", "WT N-term"])
    dosePlot(ax[8], respDF, 1, r"T$_{reg}$", ligList=["WT N-term", "WT C-term"], legend=True)
    dosePlot(ax[9], respDF, 1, r"CD8$^{+}$", ligList=["WT N-term", "WT C-term"])
    dosePlot(ax[10], respDF, 4, r"T$_{reg}$", ligList=["WT N-term", "WT C-term"])
    dosePlot(ax[11], respDF, 4, r"CD8$^{+}$", ligList=["WT N-term", "WT C-term"])
    dosePlot(ax[12], respDF, 1, r"T$_{reg}$ $IL2Ra^{hi}$", ligList=["F42Q N-Term", "H16N N-term"], legend=True)
    dosePlot(ax[13], respDF, 1, r"T$_{reg}$ $IL2Ra^{lo}$", ligList=["F42Q N-Term", "H16N N-term"])
    dosePlot(ax[14], respDF, 1, r"T$_{helper}$ $IL2Ra^{hi}$", ligList=["F42Q N-Term", "H16N N-term"])
    dosePlot(ax[15], respDF, 1, r"T$_{helper}$ $IL2Ra^{lo}$", ligList=["F42Q N-Term", "H16N N-term"])

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
    sns.scatterplot(data=mutAffDF, x="IL2Rα $K_{D}$ (nM)", y="IL2Rβ  $K_{D}$ (nM)", hue="Ligand", ax=ax, palette=ligDict)


def fullHeatMap(ax, respDF):
    """Plots the various affinities for IL-2 Muteins"""
    heatmapDF = pd.DataFrame()
    respDFhm = copy(respDF)
    respDFhm.loc[respDF.Valency == 1, "Ligand"] = respDFhm.loc[respDFhm.Valency == 1, "Ligand"] + " (Mon)"
    respDFhm.loc[respDF.Valency == 2, "Ligand"] = respDFhm.loc[respDFhm.Valency == 2, "Ligand"] + " (Biv)"
    respDFhm = respDFhm.groupby(["Ligand", "Cell", "Dose", "Time"]).Mean.mean().reset_index()
    for ligand in respDFhm.Ligand.unique():
        for dose in respDFhm.Dose.unique():
            row = pd.DataFrame()
            row["Ligand/Dose"] = [ligand + " - " + str(dose) + " (nM)"]
            for cell in respDF.Cell.unique():
                normMax = respDFhm.loc[(respDFhm.Ligand == ligand) & (respDFhm.Cell == cell)].Mean.max()
                for time in respDFhm.Time.unique():
                    entry = respDFhm.loc[(respDFhm.Ligand == ligand) & (respDFhm.Dose == dose) & (respDFhm.Cell == cell) & (respDFhm.Time == time)].Mean.values / normMax
                    if np.isnan(entry):
                        row[cell + " - " + str(time) + " hrs"] = 0
                    elif entry.size < 1:
                        row[cell + " - " + str(time) + " hrs"] = 0
                    else:
                        row[cell + " - " + str(time) + " hrs"] = entry
            heatmapDF = pd.concat([heatmapDF, row])
    heatmapDF = heatmapDF.set_index("Ligand/Dose")
    sns.heatmap(data=heatmapDF, vmin=0, vmax=1, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")


def dosePlot(ax, respDF, time, cell, ligList=False, legend=False):
    """Plots the various affinities for IL-2 Muteins"""
    doses = np.log10(np.logspace(np.log10(respDF.Dose.min()), np.log10(respDF.Dose.max()), 100)) + 4
    x0 = [4, 1, 2]
    hillDF = pd.DataFrame()
    if not ligList:
        Ligands = respDF.Ligand.unique()
    else:
        Ligands = ligList

    #maxobs = respDF.loc[(respDF.Dose == respDF.Dose.max()) & (respDF.Cell == cell) & (respDF.Ligand == "IL2")].Mean.max()
    respDF = respDF.loc[(respDF.Cell == cell) & (respDF.Time == time) & (respDF.Ligand.isin(Ligands + ["IL2"]))]

    for ligand in respDF.Ligand.unique():
        for valency in respDF.loc[respDF.Ligand == ligand].Valency.unique():
            isoData = respDF.loc[(respDF.Ligand == ligand) & (respDF.Valency == valency)]
            xData = np.nan_to_num(np.log10(isoData.Dose.values)) + 4
            yData = np.nan_to_num(isoData.Mean.values)
            fit = least_squares(hill_residuals, x0, args=(xData, yData), bounds=([0.0, 0.0, 2], [5, 10.0, 6]), jac="3-point")
            hillDF = pd.concat([hillDF, pd.DataFrame({"Ligand": ligand, "Valency": valency, "Dose": np.power(10, doses - 4), "pSTAT": hill_equation(fit.x, doses)})])

    maxobs = hillDF.loc[(hillDF.Ligand == "IL2")].pSTAT.max()
    respDF = respDF.loc[(respDF.Cell == cell) & (respDF.Time == time) & (respDF.Ligand.isin(Ligands))]
    hillDF = hillDF.loc[(hillDF.Ligand.isin(Ligands))]

    hillDF = hillDF.groupby(["Ligand", "Valency", "Dose"]).pSTAT.mean().reset_index()
    respDF = respDF.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    hillDF["pSTAT"] /= maxobs
    respDF["Mean"] /= maxobs
    sns.lineplot(data=hillDF, x="Dose", y="pSTAT", hue="Ligand", size="Valency", ax=ax, palette=ligDict, sizes=(1, 2.5))
    sns.scatterplot(data=respDF, x="Dose", y="Mean", hue="Ligand", size="Valency", ax=ax, legend=False, palette=ligDict)
    ax.set(xscale="log", xlim=(1e-4, 1e2), title=cell + " at " + str(time) + " hours", ylim=(0, 1.2))
    if legend:
        h, l = ax.get_legend_handles_labels()
        ax.legend(h[-3:], l[-3:])
    else:
        ax.get_legend().remove()


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


def PCAheatmap(ax, respDF):
    """Plots the various affinities for IL-2 Muteins"""
    heatmapDF = pd.DataFrame()
    respDFhm = copy(respDF)
    respDFhm.loc[respDF.Valency == 1, "Ligand"] = respDFhm.loc[respDFhm.Valency == 1, "Ligand"] + " (Mon)"
    respDFhm.loc[respDF.Valency == 2, "Ligand"] = respDFhm.loc[respDFhm.Valency == 2, "Ligand"] + " (Biv)"
    respDFhm = respDFhm.groupby(["Ligand", "Cell", "Dose", "Time"]).Mean.mean().reset_index()
    ligCol = []
    doseCol = []
    cellCol = []
    timeCol = []
    valCol = []

    for i, ligand in enumerate(respDFhm["Ligand"].unique()):
        for dose in respDFhm.Dose.unique():
            row = pd.DataFrame()
            row["Ligand/Dose"] = [ligand + " - " + str(dose) + " (nM)"]
            for cell in respDF.Cell.unique():
                normMax = respDFhm.loc[(respDFhm.Ligand == ligand) & (respDFhm.Cell == cell)].Mean.max()
                for time in respDFhm.Time.unique():
                    entry = respDFhm.loc[(respDFhm.Ligand == ligand) & (respDFhm.Dose == dose) & (respDFhm.Cell == cell) & (respDFhm.Time == time)].Mean.values / normMax
                    if np.isnan(entry):
                        row[cell + " - " + str(time) + " hrs"] = 0
                    elif entry.size < 1:
                        row[cell + " - " + str(time) + " hrs"] = 0
                    else:
                        row[cell + " - " + str(time) + " hrs"] = entry
            if len(ligand.split(" ")) == 2:
                ligentry = ligand.split(" ")[0]
            else:
                ligentry = ligand.split(" ")[0] + " " + ligand.split(" ")[1]
            if ligand.split(" ")[-1] == "(Mon)":
                valentry = 1
            else:
                valentry = 2
            ligCol.append(ligentry)
            doseCol.append(np.log10(dose))
            valCol.append(valentry)
            heatmapDF = pd.concat([heatmapDF, row])

    for cell in respDF.Cell.unique():
        for time in respDFhm.Time.unique():
            cellCol.append(cell)
            timeCol.append(time)

    heatmapDF = heatmapDF.set_index("Ligand/Dose")

    PCAmat = heatmapDF.to_numpy()
    scaler = StandardScaler()
    pSTATpca = scaler.fit_transform(PCAmat)
    pca = PCA(n_components=2)
    scores = pca.fit_transform(pSTATpca)
    varExp = pca.explained_variance_ratio_ * 100
    loadings = pca.components_
    scoresDF = pd.DataFrame({"Ligand": ligCol, "Dose (log)": doseCol, "Valency": valCol, "Component 1": scores[:, 0], "Component 2": scores[:, 1]})
    loadingsDF = pd.DataFrame({"Cell": cellCol, "Time": timeCol, "Component 1": loadings[0, :], "Component 2": loadings[1, :]})

    sns.scatterplot(data=scoresDF, x="Component 1", y="Component 2", hue="Ligand", size="Dose (log)", style="Valency", palette=ligDict, ax=ax[0])
    ax[0].set(xlim=(-13, 13), ylim=(-6, 6), xlabel="PC1 (" + str(varExp[0])[0:4] + "%)", ylabel="PC2 (" + str(varExp[1])[0:4] + "%)")
    #ax[0].legend(loc="lower left")
    sns.scatterplot(data=loadingsDF, x="Component 1", y="Component 2", hue="Cell", style="Time", ax=ax[1])
    ax[1].set(xlim=(-0.3, 0.3), ylim=(-0.4, 0.4), xlabel="PC1 (" + str(varExp[0])[0:4] + "%)", ylabel="PC2 (" + str(varExp[1])[0:4] + "%)")
