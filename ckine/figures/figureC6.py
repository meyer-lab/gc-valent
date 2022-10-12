"""
This creates Figure 6, assymmetric ligand design and results.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from copy import copy
from os.path import join
from .figureCommon import subplotLabel, getSetup, getLigDict, ligand_ratio_plot
from ..MBmodel import getKxStar, polyc, runFullModelMeyer
from ..flow_meyer import make_flow_df
from ..imports import importReceptors, import_pstat_all_meyer
from ..tensorFac import makeTensor, factorTensor, R2Xplot, plot_tFac_Ligs, plot_tFac_Conc, plot_tFac_Cells


path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
plt.rcParams['svg.fonttype'] = 'none'


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((9, 6.75), (3, 4))
    axlabel = copy(ax)
    del axlabel[3]
    del axlabel[-1]
    subplotLabel(axlabel)
    ax[2].axis("off")
    ax[3].axis("off")

    # make_flow_df()
    modelDF = runFullModelMeyer().reset_index()

    assymetry_Plot(ax[4])
    ligand_ratio_plot(modelDF, "Treg", "Thelper", ax[5], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "NK", ax[6], live_dead=True)
    ligand_ratio_plot(modelDF, "Treg", "CD8", ax[7], live_dead=True)

    respDF = import_pstat_all_meyer()
    respDF.loc[(respDF.Valency == 1), "Ligand"] = (respDF.loc[(respDF.Valency == 1)].Ligand + " (Mono)").values
    respDF.loc[(respDF.Valency == 2), "Ligand"] = (respDF.loc[(respDF.Valency == 2)].Ligand + " (Biv)").values
    respDF.loc[(respDF.Valency == 4), "Ligand"] = (respDF.loc[(respDF.Valency == 4)].Ligand + " (Tetra)").values

    respDF = respDF.rename({"pSTAT5": "Mean"}, axis=1)
    respDF["Time"] = 1
    numComps = 3
    respTensor = makeTensor(respDF, meyer_data=True)
    tFacAllM = factorTensor(respTensor, numComps=numComps)
    tFacAllM.normalize()

    R2Xplot(ax[8], respTensor, 5)
    plot_tFac_Ligs(ax[9], tFacAllM, respDF, numComps=numComps)
    plot_tFac_Cells(ax[10], tFacAllM, meyer_data=True, numComps=numComps)
    plot_tFac_Conc(ax[11], tFacAllM, respDF, meyer_data=True)

    return f


affsDict = {"IL2": np.array([[1.00000000e+02, 7.51879699e+09], [1.00000000e+08, 1.00000000e+02]]),
            "R38Q/H16N": np.array([[1.00000000e+02, 4.47427293e+07], [1.40845070e+09, 1.00000000e+02]]),
            "Live/Dead": np.array([[1.00000000e+02, 4.47427293e+07], [1.40845070e+09, 1.00000000e+02], [1.66666667e+08, 1.00000000e+02]])}


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
                recCount = np.ravel([recDF.loc[(recDF.Receptor == "IL2Ra") & (recDF["Cell Type"] == cell)].Mean.values,
                                     recDF.loc[(recDF.Receptor == "IL2Rb") & (recDF["Cell Type"] == cell)].Mean.values])
                if lig == "Live/Dead":
                    Affs = copy(affsDict[lig])
                    Affs[0, 1] = aff
                    bound = polyc(1e-9, getKxStar(), recCount, [[2, 2, 2]], [1.0], Affs)[0][1]
                else:
                    Affs = copy(affsDict[lig])
                    Affs[0, 1] = aff
                    bound = polyc(1e-9, getKxStar(), recCount, [[4, 4]], [1.0], Affs)[0][1]
                if cell == targCell:
                    targ += bound
                else:
                    offTarg += bound
            specDF = pd.concat([specDF, pd.DataFrame({"IL2RB Affinity": [aff], "Ligand": lig, "Potency": targ, "Selectivity": targ / offTarg})])

    for lig in ligs:
        targ = 0
        offTarg = 0
        for cell in cellTypes:
            recCount = np.ravel([recDF.loc[(recDF.Receptor == "IL2Ra") & (recDF["Cell Type"] == cell)].Mean.values,
                                 recDF.loc[(recDF.Receptor == "IL2Rb") & (recDF["Cell Type"] == cell)].Mean.values])
            if lig == "Live/Dead":
                Affs = copy(affsDict[lig])
                bound = polyc(1e-9, getKxStar(), recCount, [[2, 2, 2]], [1.0], Affs)[0][1]
            else:
                Affs = copy(affsDict[lig])
                bound = polyc(1e-9, getKxStar(), recCount, [[4, 4]], [1.0], Affs)[0][1]
            if cell == targCell:
                targ += bound
            else:
                offTarg += bound
        pointDF = pd.concat([pointDF, pd.DataFrame({"Ligand": [lig], "Potency": targ, "Selectivity": targ / offTarg})])

    pointDF.Potency /= specDF.Potency.max()
    pointDF.Selectivity /= specDF.Selectivity.max()
    specDF.Potency /= specDF.Potency.max()
    specDF.Selectivity /= specDF.Selectivity.max()
    specDF = specDF.reset_index().drop("index", axis=1)
    pointDF = pointDF.reset_index().drop("index", axis=1)

    sns.scatterplot(data=pointDF, x="Potency", y="Selectivity", hue="Ligand", ax=ax, palette=ligDict)
    sns.lineplot(data=specDF, x="Potency", y="Selectivity", hue="Ligand", ax=ax, palette=ligDict)


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
