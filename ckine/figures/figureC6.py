"""
This creates Figure 4, fitting of multivalent binding model to Gc Data.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
from copy import copy
from matplotlib import rcParams
from .figureCommon import subplotLabel, getSetup, get_cellTypeDict, get_doseLimDict, getLigDict
from ..MBmodel import getKxStar, polyc

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
rcParams['svg.fonttype'] = 'none'


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((3, 3), (1, 1))
    subplotLabel(ax)

    # minSolved = minimize(runFullModel, x0=-12.0, args=([0.5, 1], False, True))
    # print(minSolved)

    YT1_Plot(ax[0], estRec=500)
    #YT1_Plot(ax[1], estRec=2500)
    #YT1_Plot(ax[2], estRec=5000)

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
        recBoundDF = recBoundDF.append(pd.DataFrame({"Concentration": [dose], "Valency": [val], "Ligand": row.Ligand, "Active Binding Complexes": predVal}))

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
                outputDF = outputDF.append(pd.DataFrame({"Concentration": [dose], "Valency": [valency], "Ligand": mut, "pSTAT5": predVal, "Type": "Predicted"}))
    print(outputDF)
    outputDF = outputDF.reset_index()
    #sns.lineplot(data=outputDF, x="Concentration", y="pSTAT5", style="Valency", hue="Ligand", palette=ligDict, ax=ax)
    sns.scatterplot(data=expData, x="Concentration", y="Mean", style="Valency", hue="Ligand", palette=ligDict, ax=ax)
    ax.set(xscale="log", xlim=(0.01, 1000), ylabel="pSTAT5")
