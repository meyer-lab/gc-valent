"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
from os.path import dirname, join
from .figureCommon import subplotLabel, getSetup
from ..imports import import_pstat_all

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((12, 6), (2, 4))
    subplotLabel(ax)
    ax[0].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ  $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    doses = respDF.Dose.unique()
    dose = doses[4]
    ratioConc(ax[7], respDF, "Treg", "NK", 0.5)
    pSTATcomp(ax[1], respDF, mutAffDF, "IL2Rα $K_{D}$ (nM)", "Treg", 1, dose, legend=True)
    pSTATcomp(ax[2], respDF, mutAffDF, "IL2Rα $K_{D}$ (nM)", "NK", 1, dose)
    pSTATcomp(ax[3], respDF, mutAffDF, "IL2Rβ  $K_{D}$ (nM)", "Treg", 1, dose)
    pSTATcomp(ax[4], respDF, mutAffDF, "IL2Rβ  $K_{D}$ (nM)", "NK", 1, dose)
    pSTATcomp(ax[5], respDF, mutAffDF, "IL2Rα $K_{D}$ (nM)", "CD8", 1, dose)
    pSTATcomp(ax[6], respDF, mutAffDF, "IL2Rβ  $K_{D}$ (nM)", "CD8", 1, dose)
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
    sns.scatterplot(data=respDF, x=recAff, y="Mean", hue="Ligand", style="Valency", ax=ax, legend=legend)
    ax.set(xscale="log", title=cell)


def ratioConc(ax, respDF, cell1, cell2, time):
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
    sns.scatterplot(data=ratioDF, x="Dose", y="Ratio", hue="Ligand", style="Valency", ax=ax)
    ax.set(xscale="log")
