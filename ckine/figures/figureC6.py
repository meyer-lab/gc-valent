"""
This creates Figure 4, fitting of multivalent binding model to Gc Data.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
from copy import copy
from .figureCommon import (
    subplotLabel,
    getSetup,
    getLigDict,
    ligandPlot,
    ligand_ratio_plot,
)
from ..MBmodel import getKxStar, polyc, runFullModelMeyer
from ..flow_meyer import make_flow_df
from ..imports import importReceptors

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((12, 9), (3, 4), multz={2: 1})
    subplotLabel(ax)

    # make_flow_df()
    modelDF = runFullModelMeyer().reset_index()

    assymetry_Plot(ax[3])
    ligandPlot(modelDF, "Treg", ax[4], live_dead=True)
    ligandPlot(modelDF, "Thelper", ax[5], live_dead=True)
    ligandPlot(modelDF, "NK", ax[6], live_dead=True)
    ligandPlot(modelDF, "CD8", ax[7], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "Thelper", ax[8], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "NK", ax[9], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "CD8", ax[10], live_dead=True)

    return f


affsDict = {
    "IL2": np.array([[1.00000000e02, 7.51879699e09], [1.00000000e08, 1.00000000e02]]),
    "R38Q/H16N": np.array(
        [[1.00000000e02, 4.47427293e07], [1.40845070e09, 1.00000000e02]]
    ),
    "Live/Dead": np.array(
        [
            [1.00000000e02, 4.47427293e07],
            [1.40845070e09, 1.00000000e02],
            [1.66666667e08, 1.00000000e02],
        ]
    ),
}


def assymetry_Plot(ax):
    """Plots theoretical selectivity vs. potency for Tregs for three constructs - WT, R38Q/H16N and Live/Dead"""
    IL2RBaff = np.logspace(6, 10, num=100)
    targCell = "Treg"
    cellTypes = ["Treg", "Thelper", "CD8", "NK"]
    ligs = ["IL2", "R38Q/H16N", "Live/Dead"]
    recDF = importReceptors()
    specDF = pd.DataFrame()
    pointDF = pd.DataFrame()

    for aff in IL2RBaff:
        for lig in ligs:
            targ = 0
            offTarg = 0
            for cell in cellTypes:
                recCount = np.ravel(
                    [
                        recDF.loc[
                            (recDF.Receptor == "IL2Ra") & (recDF["Cell Type"] == cell)
                        ].Mean.values,
                        recDF.loc[
                            (recDF.Receptor == "IL2Rb") & (recDF["Cell Type"] == cell)
                        ].Mean.values,
                    ]
                )
                if lig == "Live/Dead":
                    Affs = copy(affsDict[lig])
                    Affs[0, 1] = aff
                    bound = polyc(
                        1e-9, getKxStar(), recCount, [[2, 2, 2]], [1.0], Affs
                    )[0][1]
                else:
                    Affs = copy(affsDict[lig])
                    Affs[0, 1] = aff
                    bound = polyc(1e-9, getKxStar(), recCount, [[4, 4]], [1.0], Affs)[
                        0
                    ][1]
                if cell == targCell:
                    targ += bound
                else:
                    offTarg += bound
            specDF = pd.concat(
                [
                    specDF,
                    pd.DataFrame(
                        {
                            "IL2RB Affinity": [aff],
                            "Ligand": lig,
                            "Potency": targ,
                            "Selectivity": targ / offTarg,
                        }
                    ),
                ]
            )

    for lig in ligs:
        targ = 0
        offTarg = 0
        for cell in cellTypes:
            recCount = np.ravel(
                [
                    recDF.loc[
                        (recDF.Receptor == "IL2Ra") & (recDF["Cell Type"] == cell)
                    ].Mean.values,
                    recDF.loc[
                        (recDF.Receptor == "IL2Rb") & (recDF["Cell Type"] == cell)
                    ].Mean.values,
                ]
            )
            if lig == "Live/Dead":
                Affs = copy(affsDict[lig])
                bound = polyc(1e-9, getKxStar(), recCount, [[2, 2, 2]], [1.0], Affs)[0][
                    1
                ]
            else:
                Affs = copy(affsDict[lig])
                bound = polyc(1e-9, getKxStar(), recCount, [[4, 4]], [1.0], Affs)[0][1]
            if cell == targCell:
                targ += bound
            else:
                offTarg += bound
        pointDF = pd.concat(
            [
                pointDF,
                pd.DataFrame(
                    {"Ligand": [lig], "Potency": targ, "Selectivity": targ / offTarg}
                ),
            ]
        )

    pointDF.Potency /= specDF.Potency.max()
    pointDF.Selectivity /= specDF.Selectivity.max()
    specDF.Potency /= specDF.Potency.max()
    specDF.Selectivity /= specDF.Selectivity.max()
    specDF = specDF.reset_index().drop("index", axis=1)
    pointDF = pointDF.reset_index().drop("index", axis=1)

    sns.scatterplot(
        data=pointDF, x="Potency", y="Selectivity", hue="Ligand", ax=ax, palette=ligDict
    )
    sns.lineplot(
        data=specDF, x="Potency", y="Selectivity", hue="Ligand", ax=ax, palette=ligDict
    )


def YT1_Plot(ax, estRec):
    """Tracks interactions of Valency and ligand affinity on Treg activation"""
    expData = pd.read_csv(join(path_here, "data/YT1TetraSTAT.csv"))
    expData.Mean = expData.Mean.clip(lower=0)
    concs = np.logspace(
        start=np.log10(expData.Concentration.min()),
        stop=np.log10(expData.Concentration.max()),
        num=101,
        endpoint=True,
    )
    outputDF = pd.DataFrame()
    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    recCount = np.ravel([0, estRec])
    recBoundDF = pd.DataFrame()

    for rowTup in expData.iterrows():
        row = rowTup[1]
        Affs = mutAffDF.loc[(mutAffDF.Mutein == row.Ligand)]
        Affs = np.power(
            np.array([Affs["IL2RaKD"].values, Affs["IL2RBGKD"].values]) / 1e9, -1
        )
        Affs = np.reshape(Affs, (1, -1))
        Affs = np.repeat(Affs, 2, axis=0)
        np.fill_diagonal(Affs, 1e2)  # Each cytokine can only bind one a and one b
        val = row.Valency
        dose = row.Concentration
        predVal = polyc(dose / 1e9, getKxStar(), recCount, [[val, val]], [1.0], Affs)[
            0
        ][1]
        recBoundDF = pd.concat(
            [
                recBoundDF,
                pd.DataFrame(
                    {
                        "Concentration": [dose],
                        "Valency": [val],
                        "Ligand": row.Ligand,
                        "Active Binding Complexes": predVal,
                    }
                ),
            ]
        )

    convFact = expData.Mean.mean() / recBoundDF["Active Binding Complexes"].mean()

    for mut in expData.Ligand.unique():
        for valency in expData.loc[expData.Ligand == mut].Valency.unique():
            for dose in concs:
                Affs = mutAffDF.loc[(mutAffDF.Mutein == mut)]
                Affs = np.power(
                    np.array([Affs["IL2RaKD"].values, Affs["IL2RBGKD"].values]) / 1e9,
                    -1,
                )
                Affs = np.reshape(Affs, (1, -1))
                Affs = np.repeat(Affs, 2, axis=0)
                np.fill_diagonal(
                    Affs, 1e2
                )  # Each cytokine can only bind one a and one b
                predVal = (
                    polyc(
                        dose / 1e9,
                        getKxStar(),
                        recCount,
                        [[valency, valency]],
                        [1.0],
                        Affs,
                    )[0][1]
                    * convFact
                )
                outputDF = pd.concat(
                    [
                        outputDF,
                        pd.DataFrame(
                            {
                                "Concentration": [dose],
                                "Valency": [valency],
                                "Ligand": mut,
                                "pSTAT5": predVal,
                                "Type": "Predicted",
                            }
                        ),
                    ]
                )
    outputDF = outputDF.reset_index()
    sns.lineplot(
        data=outputDF,
        x="Concentration",
        y="pSTAT5",
        style="Valency",
        hue="Ligand",
        palette=ligDict,
        ax=ax,
    )
    sns.scatterplot(
        data=expData,
        x="Concentration",
        y="Mean",
        style="Valency",
        hue="Ligand",
        palette=ligDict,
        ax=ax,
    )
    ax.set(xscale="log", ylabel="pSTAT5")
