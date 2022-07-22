"""
This creates Figure 6, assymmetric ligand design and results.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
from .figureCommon import subplotLabel, getSetup, getLigDict, ligandPlot, ligand_ratio_plot
from ..MBmodel import getKxStar, polyc, runFullModelMeyer
from ..flow_meyer import make_flow_df

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 9), (4, 4))
    subplotLabel(ax)

    # make_flow_df()
    modelDF = runFullModelMeyer().reset_index()

    ligandPlot(modelDF, "Treg", ax[3], live_dead=True)
    ligandPlot(modelDF, "Thelper", ax[4], live_dead=True)
    ligandPlot(modelDF, "NK", ax[5], live_dead=True)
    ligandPlot(modelDF, "CD8", ax[6], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "Thelper", ax[7], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "NK", ax[8], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "CD8", ax[9], live_dead=True)

    return f


def YT1_Plot(ax, estRec):
    """Tracks interactions of Valency and ligand affinity on Treg activation"""
    expData = pd.read_csv(join(path_here, "data/YT1TetraSTAT.csv"))
    expData.Mean = expData.Mean.clip(lower=0)
    concs = np.logspace(start=np.log10(expData.Concentration.min()), stop=np.log10(expData.Concentration.max()), num=101, endpoint=True)
    outputDF = pd.DataFrame()
    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    recCount = np.ravel([0, estRec])
    recBoundDF = pd.DataFrame()

    for rowTup in expData.iterrows():
        row = rowTup[1]
        Affs = mutAffDF.loc[(mutAffDF.Mutein == row.Ligand)]
        Affs = np.power(np.array([Affs["IL2RaKD"].values, Affs["IL2RBGKD"].values]) / 1e9, -1)
        Affs = np.reshape(Affs, (1, -1))
        Affs = np.repeat(Affs, 2, axis=0)
        np.fill_diagonal(Affs, 1e2)  # Each cytokine can only bind one a and one b
        val = row.Valency
        dose = row.Concentration
        predVal = polyc(dose / 1e9, getKxStar(), recCount, [[val, val]], [1.0], Affs)[0][1]
        recBoundDF = pd.concat([recBoundDF, pd.DataFrame({"Concentration": [dose], "Valency": [val], "Ligand": row.Ligand, "Active Binding Complexes": predVal})])

    convFact = expData.Mean.mean() / recBoundDF["Active Binding Complexes"].mean()

    for mut in expData.Ligand.unique():
        for valency in expData.loc[expData.Ligand == mut].Valency.unique():
            for dose in concs:
                Affs = mutAffDF.loc[(mutAffDF.Mutein == mut)]
                Affs = np.power(np.array([Affs["IL2RaKD"].values, Affs["IL2RBGKD"].values]) / 1e9, -1)
                Affs = np.reshape(Affs, (1, -1))
                Affs = np.repeat(Affs, 2, axis=0)
                np.fill_diagonal(Affs, 1e2)  # Each cytokine can only bind one a and one b
                predVal = polyc(dose / 1e9, getKxStar(), recCount, [[valency, valency]], [1.0], Affs)[0][1] * convFact
                outputDF = pd.concat([outputDF, pd.DataFrame({"Concentration": [dose], "Valency": [valency], "Ligand": mut, "pSTAT5": predVal, "Type": "Predicted"})])
    outputDF = outputDF.reset_index()
    sns.lineplot(data=outputDF, x="Concentration", y="pSTAT5", style="Valency", hue="Ligand", palette=ligDict, ax=ax)
    sns.scatterplot(data=expData, x="Concentration", y="Mean", style="Valency", hue="Ligand", palette=ligDict, ax=ax)
    ax.set(xscale="log", ylabel="pSTAT5")
