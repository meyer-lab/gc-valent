"""
This creates Figure 4, fitting of multivalent binding model to Gc Data.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import join
from .figureCommon import subplotLabel, getSetup, get_cellTypeDict, get_doseLimDict, getLigDict
from ..MBmodel import getKxStar, polyc
from ..flow_meyer import make_flow_df

path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 9), (4, 4))
    subplotLabel(ax)

    # minSolved = minimize(runFullModel, x0=-12.0, args=([0.5, 1], False, True))
    # print(minSolved)

    YT1_Plot(ax[0], estRec=500)
    YT1_Plot(ax[1], estRec=2500)
    YT1_Plot(ax[2], estRec=5000)
    # make_flow_df()

    MeyerDF = pd.read_csv(join(path_here, "data/Meyer_Flow.csv"))
    ligandPlot(MeyerDF, "Treg", ax[3], live_dead=False)
    ligandPlot(MeyerDF, "Thelper", ax[4], live_dead=False)
    ligandPlot(MeyerDF, "NK", ax[5], live_dead=False)
    ligandPlot(MeyerDF, "CD8+", ax[6], live_dead=False)
    ligand_ratio_plot(MeyerDF, "Treg", "Thelper", ax[7], live_dead=False)
    ligand_ratio_plot(MeyerDF, "Treg", "NK", ax[8], live_dead=False)
    ligand_ratio_plot(MeyerDF, "Treg", "CD8+", ax[9], live_dead=False)

    return f


palette_dict = {"R38Q/H16N N-term": "darkorchid",
                "Live/Dead": "Orange"}


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


def ligandPlot(DF, cell, ax, live_dead=False):
    """Plots a cell type response"""
    if live_dead:
        DF = DF.loc[DF.Cell == cell]
    else:
        DF = DF.loc[(DF.Cell == cell) & (DF.Ligand != "Live/Dead")]
    sns.lineplot(data=DF, x="Dose", y="pSTAT5", hue="Ligand", style="Valency", palette=palette_dict, ax=ax)
    ax.set(xscale="log", xlabel="Dose (nM)", ylabel="pSTAT5 (MFI)", title=cell)


def ligand_ratio_plot(DF, cell1, cell2, ax, live_dead=False):
    """Plots a cell type response"""
    ratioDF = pd.DataFrame()
    if live_dead:
        for dose in DF.Dose.unique():
            for ligand in DF.Ligand.unique():
                for valency in DF.loc[DF.Ligand == ligand].Valency.unique():
                    ratio = DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell1)].pSTAT5.values / \
                        DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell2)].pSTAT5.values
                    ratioDF = pd.concat([ratioDF, pd.DataFrame({"Dose": dose, "Ligand": ligand, "Valency": valency, "Ratio": ratio})])

    else:
        DF = DF.loc[(DF.Ligand != "Live/Dead")]
        for dose in DF.Dose.unique():
            for ligand in DF.Ligand.unique():
                for valency in DF.loc[DF.Ligand == ligand].Valency.unique():
                    ratio = DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell1)].pSTAT5.values / \
                        DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell2)].pSTAT5.values
                    ratioDF = pd.concat([ratioDF, pd.DataFrame({"Dose": dose, "Ligand": ligand, "Valency": valency, "Ratio": ratio})])

    ratioDF = ratioDF.reset_index()
    sns.lineplot(data=ratioDF, x="Dose", y="Ratio", hue="Ligand", style="Valency", palette=palette_dict, ax=ax)
    ax.set(xscale="log", xlabel="Dose (nM)", ylabel="Ratio", title=cell1 + " to " + cell2 + " Ratio")
