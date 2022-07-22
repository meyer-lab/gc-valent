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
from .figureCommon import (
    subplotLabel,
    getSetup,
    getLigDict,
    get_cellTypeDict,
    getLigandLegend,
    Wass_KL_Dist,
    CITE_RIDGE,
    CITE_SVM,
)
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
    subplotLabel(axlabel)
    ax[0].axis("off")
    ax[1].axis("off")
    ax[11].axis("off")
    ax[12].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename(
        {
            "Mutein": "Ligand",
            "IL2RaKD": "IL2Rα $K_{D}$ (nM)",
            "IL2RBGKD": "IL2Rβ $K_{D}$ (nM)",
        },
        axis=1,
    )

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    time = 4.0
    pseudo = 0.15

    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]

    ratioConc(
        ax[2:5], respDF, r"T$_{reg}$", "NK", time, mutAffDF, pseudo=pseudo, legend=True
    )
    ratioConc(
        ax[5:8],
        respDF,
        r"T$_{reg}$",
        r"CD8$^{+}$",
        time,
        mutAffDF,
        pseudo=pseudo,
        legend=True,
    )
    ratioConc(
        ax[8:11],
        respDF,
        r"T$_{reg}$",
        r"T$_{helper}$",
        time,
        mutAffDF,
        pseudo=pseudo,
        legend=True,
    )

    legend = getLigandLegend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[1].legend(
        legend.legendHandles, labels, loc="upper left", prop={"size": 10}
    )  # use this to place universal legend later
    cellTarget = "Treg"
    Wass_KL_Dist(ax[13:15], cellTarget, 10)
    Wass_KL_Dist(ax[15:17], cellTarget, 10, RNA=True)
    CITE_RIDGE(ax[17], cellTarget)
    CITE_SVM(ax[18], cellTarget, sampleFrac=0.2)

    return f


def pSTATcomp(ax, respDF, mutAffDF, recAff, cell, time, conc, legend=False):
    """Plots monovalent and bivalent lines for a single cell type by affinity"""
    respDF = respDF.loc[
        (respDF.Cell == cell) & ((respDF.Time == time) & (respDF.Dose == conc))
    ]
    respDF = respDF.groupby(["Ligand", "Valency"]).Mean.mean().reset_index()
    respDF = respDF.merge(mutAffDF)
    respDF = respDF.loc[(respDF.Ligand != "IL15") & (respDF.Ligand != "IL2")]
    sns.scatterplot(
        data=respDF,
        x=recAff,
        y="Mean",
        hue="Ligand",
        style="Valency",
        ax=ax,
        legend=legend,
        palette=ligDict,
    )
    ax.set(
        xscale="log",
        title=cell + " at " + str(conc) + " nM  and " + str(time) + " hours",
    )


def gaussian(x, concs):
    """Returns Gaussian for given concs"""
    amplitude, mean, stddev = x[0], x[1], x[2]
    return amplitude * np.exp(-(((concs - mean) / 4 / stddev) ** 2))


def gaussian_residuals(x, concs, ratios):
    """Residual function for Gaussian."""
    return gaussian(x, concs) - ratios


def ratioConc(ax, respDF, cell1, cell2, time, mutAffDF, pseudo=0.1, legend=False):
    """Plots Ratio of cell 1 to cell 2 over a range of concentrations"""
    hillDF = hillRatioDosePlot(ax[0], respDF, time, cell1, cell2, pseudo=pseudo)
    fitDF = pd.DataFrame()
    for ligand in hillDF.Ligand.unique():
        for valency in hillDF.loc[hillDF.Ligand == ligand].Valency.unique():
            isoData = hillDF.loc[
                (hillDF.Ligand == ligand) & (hillDF.Valency == valency)
            ]
            fitDF = pd.concat(
                [
                    fitDF,
                    pd.DataFrame(
                        {
                            "Ligand": [ligand],
                            "Valency": valency,
                            cell2 + " Max": isoData.Ratio.max(),
                            cell2
                            + " Dose": isoData.loc[
                                isoData.Ratio == isoData.Ratio.max()
                            ].Dose.values,
                        }
                    ),
                ]
            )

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
        mMax, bMax = np.polyfit(
            np.log10(valData["IL2Rα $K_{D}$ (nM)"].values), valData[cell2 + " Max"], 1
        )
        mDose, bDose = np.polyfit(
            np.log10(valData["IL2Rα $K_{D}$ (nM)"].values),
            np.log10(valData[cell2 + " Dose"]),
            1,
        )
        maxLineDF = pd.concat(
            [
                maxLineDF,
                pd.DataFrame(
                    {
                        "Valency": valency,
                        "IL2Rα $K_{D}$ (nM)": np.power(10, affs),
                        cell2 + " Max": mMax * affs + bMax,
                    }
                ),
            ]
        )
        doseLineDF = pd.concat(
            [
                doseLineDF,
                pd.DataFrame(
                    {
                        "Valency": valency,
                        "IL2Rα $K_{D}$ (nM)": np.power(10, affs),
                        cell2 + " Dose": np.power(10, mDose * affs + bDose),
                    }
                ),
            ]
        )

    maxLineDF, doseLineDF = maxLineDF.reset_index(), doseLineDF.reset_index()
    sns.scatterplot(
        data=fitDF,
        x="IL2Rα $K_{D}$ (nM)",
        y=cell2 + " Max",
        hue="Ligand",
        style="Valency",
        ax=ax[1],
        palette=ligDict,
        legend=False,
    )
    sns.lineplot(
        data=maxLineDF,
        x="IL2Rα $K_{D}$ (nM)",
        y=cell2 + " Max",
        style="Valency",
        ax=ax[1],
        color="k",
        linewidth=1.0,
        legend=False,
    )
    ax[1].set(
        xscale="log",
        title="Ratio of " + cell1 + " to " + cell2,
        xlim=(1e-1, 1e1),
        ylim=(0, None),
        ylabel=cell1 + " to " + cell2 + " Ratio Max Magnitude",
    )

    sns.scatterplot(
        data=fitDF,
        x="IL2Rα $K_{D}$ (nM)",
        y=cell2 + " Dose",
        hue="Ligand",
        style="Valency",
        ax=ax[2],
        palette=ligDict,
        legend=False,
    )
    sns.lineplot(
        data=doseLineDF,
        x="IL2Rα $K_{D}$ (nM)",
        y=cell2 + " Dose",
        style="Valency",
        ax=ax[2],
        color="k",
        linewidth=1.0,
        legend=False,
    )
    ax[2].set(
        xscale="log",
        yscale="log",
        title="Ratio of " + cell1 + " to " + cell2,
        xlim=(1e-1, 1e1),
        ylim=(1e-2, 1e2),
        ylabel=cell1 + " to " + cell2 + " Ratio Max Dose",
    )


def hillRatioDosePlot(ax, respDF, time, targCell, offTargCell, pseudo=0.2):
    """Plots the various affinities for IL-2 Muteins"""
    doses = (
        np.log10(
            np.logspace(np.log10(respDF.Dose.min()), np.log10(respDF.Dose.max()), 100)
        )
        + 4
    )
    x0 = [4, 1, 2]
    hillDF = pd.DataFrame()
    Ligands = respDF.Ligand.unique()
    respDF = respDF.loc[(respDF.Time == time)]

    for ligand in respDF.Ligand.unique():
        for valency in respDF.loc[respDF.Ligand == ligand].Valency.unique():
            targIsoData = respDF.loc[
                (respDF.Ligand == ligand)
                & (respDF.Valency == valency)
                & (respDF.Cell == targCell)
            ]
            targXData = np.nan_to_num(np.log10(targIsoData.Dose.values)) + 4
            targYData = np.nan_to_num(targIsoData.Mean.values)
            targFit = least_squares(
                hill_residuals,
                x0,
                args=(targXData, targYData),
                bounds=([0.0, 0.0, 2], [5, 10.0, 6]),
                jac="3-point",
            )

            offTIsoData = respDF.loc[
                (respDF.Ligand == ligand)
                & (respDF.Valency == valency)
                & (respDF.Cell == offTargCell)
            ]
            offTargXData = np.nan_to_num(np.log10(offTIsoData.Dose.values)) + 4
            offTargYData = np.nan_to_num(offTIsoData.Mean.values)
            offTargFit = least_squares(
                hill_residuals,
                x0,
                args=(offTargXData, offTargYData),
                bounds=([0.0, 0.0, 2], [5, 10.0, 6]),
                jac="3-point",
            )
            hillDF = pd.concat(
                [
                    hillDF,
                    pd.DataFrame(
                        {
                            "Ligand": ligand,
                            "Valency": valency,
                            "Cell": targCell,
                            "Dose": np.power(10, doses - 4),
                            targCell: hill_equation(targFit.x, doses),
                            offTargCell: hill_equation(offTargFit.x, doses),
                        }
                    ),
                ]
            )

    for cell in [targCell, offTargCell]:
        maxobs = hillDF.loc[(hillDF.Ligand == "IL2")][cell].max()
        hillDF[cell] /= maxobs
    hillDF["Ratio"] = hillDF[targCell] / (pseudo + hillDF[offTargCell])

    hillDF = hillDF.groupby(["Ligand", "Valency", "Dose"]).Ratio.mean().reset_index()
    hillDF = hillDF.loc[(hillDF.Ligand != "IL15") & (hillDF.Ligand != "IL2")]
    sns.lineplot(
        data=hillDF,
        x="Dose",
        y="Ratio",
        hue="Ligand",
        size="Valency",
        ax=ax,
        palette=ligDict,
        sizes=(1, 2.5),
    )
    ax.set(xscale="log", xlim=(1e-4, 1e2), ylim=(0, 5))
    return hillDF


def hill_equation(x, dose):
    """Calculates EC50 from Hill Equation."""
    # print(x, dose)
    EMax = np.power(10, x[0])
    n = x[1]
    EC50 = x[2]
    return EMax * np.power(dose, n) / (np.power(EC50, n) + np.power(dose, n))


def hill_residuals(x, dose, y):
    """Residual function for Hill Equation."""
    return hill_equation(x, dose) - y
