"""
Figure 5. Optimization of ligands and experimental tetravalent reads
"""
from os.path import dirname
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
from scipy.optimize import minimize, Bounds, NonlinearConstraint
from .figureCommon import subplotLabel, getSetup, ligandPlot, ligand_ratio_plot
from ..MBmodel import polyc, getKxStar, runFullModelMeyer
from ..imports import getBindDict, importReceptors

path_here = dirname(dirname(__file__))
rcParams['svg.fonttype'] = 'none'


def makeFigure():
    """ Make figure 6. """
    # Get list of axis objects
    ax, f = getSetup((8, 12), (5, 3))
    subplotLabel(ax)
    #optimizeDesign([ax[0], ax[3]], ["Treg"], ["Thelper", "NK", "CD8"])
    #optimizeDesign([ax[1], ax[4]], ["NK"], ["Thelper", "Treg", "CD8"], legend=False)
    #optimizeDesign([ax[2], ax[5]], ["Thelper"], ["Treg", "NK", "CD8"], IL7=True, legend=False)

    modelDF = runFullModelMeyer().reset_index()
    modelDF = modelDF.loc[modelDF.Date == "7/22/22"]

    ligandPlot(modelDF, "Treg", ax[6], live_dead=False)
    ligandPlot(modelDF, "Thelper", ax[7], live_dead=False)
    ligandPlot(modelDF, "NK", ax[8], live_dead=False)
    ligandPlot(modelDF, "CD8", ax[9], live_dead=False)
    ligand_ratio_plot(modelDF, "Treg", "Thelper", ax[10], live_dead=False)
    ligand_ratio_plot(modelDF, "Treg", "NK", ax[11], live_dead=False)
    ligand_ratio_plot(modelDF, "Treg", "CD8", ax[12], live_dead=False)

    return f


cellTypeDict = {"Treg": r"T$_{reg}$",
                "Thelper": r"T$_{helper}$",
                "NK": "NK",
                "CD8": r"CD8$^{+}$"}


def cytBindingModelOpt(x, val, cellType, IL7=False):
    """Runs binding model for a given mutein, valency, dose, and cell type. """
    recQuantDF = importReceptors()
    convDict = getBindDict()
    Kx = getKxStar()

    if IL7:
        affs = [[np.power(10, x[0])]]
        recCount = np.ravel([recQuantDF.loc[(recQuantDF.Receptor == "IL7Ra") & (recQuantDF["Cell Type"] == cellType)].Mean.values])
        output = polyc(1e-9 / val, Kx, recCount, [[val]], [1.0], affs)[0][0]  # IL7Ra binding only
    else:
        affs = [[np.power(10, x[0]), 1e2], [1e2, np.power(10, x[1])]]
        recCount = np.ravel([recQuantDF.loc[(recQuantDF.Receptor == "IL2Ra") & (recQuantDF["Cell Type"] == cellType)].Mean.values,
                             recQuantDF.loc[(recQuantDF.Receptor == "IL2Rb") & (recQuantDF["Cell Type"] == cellType)].Mean.values])
        output = polyc(1e-9 / val, Kx, recCount, [[val, val]], [1.0], affs)[0][1]  # IL2RB binding only

    if not IL7:
        output *= np.mean(convDict.loc[(convDict.Cell == cellType)].Scale.values)

    return output


def minSelecFunc(x, val, targCell, offTCells, IL7=False):
    """Provides the function to be minimized to get optimal selectivity"""
    offTargetBound = 0

    targetBound = cytBindingModelOpt(x, val, targCell[0], IL7)
    for cellT in offTCells:
        offTargetBound += cytBindingModelOpt(x, val, cellT, IL7)

    return (offTargetBound) / (targetBound)


def optimizeDesign(ax, targCell, offTcells, IL7=False, legend=True):
    """ A more general purpose optimizer """
    vals = np.arange(1.01, 10, step=0.15)
    sigDF = pd.DataFrame()

    if IL7:
        optDF = pd.DataFrame(columns={"Valency", "Selectivity", "IL7Rα"})
        X0 = [8]  # Ka IL7
    else:
        optDF = pd.DataFrame(columns={"Valency", "Selectivity", "IL2Rα", r"IL-2Rβ/γ$_c$"})
        if targCell[0] == "NK":
            X0 = [6.0, 8]  # IL2Ra, IL2Rb
        else:
            X0 = [9.0, 6.0]  # IL2Ra, IL2Rb

    optBnds = Bounds(np.full_like(X0, 6.0), np.full_like(X0, 9.0))

    for i, val in enumerate(vals):
        if i == 0:
            optimized = minimize(minSelecFunc, X0, bounds=optBnds, args=(val, targCell, offTcells, IL7), jac="3-point")

            targLB = cytBindingModelOpt(optimized.x, val, targCell[0], IL7) / 1.01
            bindConst = NonlinearConstraint(lambda x: cytBindingModelOpt(x, val, targCell[0], IL7), targLB, np.inf)
        else:
            optimized = minimize(minSelecFunc, X0, bounds=optBnds, args=(val, targCell, offTcells, IL7), jac="3-point", constraints=bindConst)

        fitX = 1.0e9 / np.power(10.0, optimized.x)

        if IL7:
            optDF = pd.concat([optDF, pd.DataFrame({"Valency": [val], "Selectivity": [len(offTcells) / optimized.fun], "IL7Rα": fitX[0]})])
            sigDF = pd.concat([sigDF, pd.DataFrame({"Cell Type": [targCell[0]], "Target": ["Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, targCell[0], IL7)]})])
            for cell in offTcells:
                sigDF = pd.concat([sigDF, pd.DataFrame({"Cell Type": [cell], "Target": ["Off-Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, cell, IL7)]})])
        else:
            optDF = pd.concat([optDF, pd.DataFrame({"Valency": [val], "Selectivity": [len(offTcells) / optimized.fun], "IL2Rα": fitX[0], r"IL-2Rβ/γ$_c$": fitX[1]})])
            sigDF = pd.concat([sigDF, pd.DataFrame({"Cell Type": [targCell[0]], "Target": ["Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, targCell[0], IL7)]})])
            for cell in offTcells:
                sigDF = pd.concat([sigDF, pd.DataFrame({"Cell Type": [cell], "Target": ["Off-Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, cell, IL7)]})])
    # Normalize to valency 1
    for cell in targCell + offTcells:
        sigDF.loc[sigDF["Cell Type"] == cell, "pSTAT"] = sigDF.loc[sigDF["Cell Type"] == cell, "pSTAT"].div(sigDF.loc[(sigDF["Cell Type"] == cell) & (sigDF.Valency == vals[0])].pSTAT.values[0])

    sigDF = sigDF.replace(cellTypeDict)
    sigDF = sigDF.reset_index()
    optDF = optDF.reset_index()

    if IL7:
        sns.lineplot(x="Valency", y="pSTAT", hue="Cell Type", style="Target", data=sigDF, ax=ax[0], palette="husl", hue_order=cellTypeDict.values())
        ax[0].set_title(cellTypeDict[targCell[0]] + " selectivity with IL-7 mutein", fontsize=7)

        sns.lineplot(x="Valency", y="IL7Rα", data=optDF, ax=ax[1], palette="crest")
        ax[1].set(yscale="log", ylabel=r"IL7·7Rα $K_D$ (nM)")

    else:
        sns.lineplot(x="Valency", y="pSTAT", hue="Cell Type", style="Target", data=sigDF, ax=ax[0], palette="husl", hue_order=cellTypeDict.values())
        ax[0].set_title(cellTypeDict[targCell[0]] + " selectivity with IL-2 mutein", fontsize=7)

        if targCell[0] == "NK":
            affDF = pd.melt(optDF, id_vars=['Valency'], value_vars=[r"IL-2Rβ/γ$_c$"])
            sns.lineplot(x="Valency", y="value", data=affDF, ax=ax[1])
            ax[1].set(yscale="log", ylabel=r"IL2·β/γ$_c$ K$_D$ (nM)")
        else:
            affDF = pd.melt(optDF, id_vars=['Valency'], value_vars=['IL2Rα', r"IL-2Rβ/γ$_c$"])
            affDF = affDF.rename(columns={"variable": "Receptor"})
            sns.lineplot(x="Valency", y="value", hue="Receptor", data=affDF, ax=ax[1])
            ax[1].set(yscale="log", ylabel=r"IL2· $K_D$ (nM)")

    ax[0].set_ylim(bottom=0.0, top=3)
    ax[1].set_ylim(bottom=0.1, top=2000)

    ax[0].set_xticks(np.arange(1, 11))
    ax[1].set_xticks(np.arange(1, 11))
    if not legend:
        ax[0].get_legend().remove()
