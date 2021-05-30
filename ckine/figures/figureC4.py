"""
Figure 6. Optimization of Ligands
"""
from os.path import dirname, join
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import minimize
from .figureCommon import subplotLabel, getSetup
from ..MBmodel import polyc, getKxStar
from ..imports import getBindDict, importReceptors

path_here = dirname(dirname(__file__))


def makeFigure():
    """ Make figure 6. """
    # Get list of axis objects
    ax, f = getSetup((6, 4.5), (2, 3))
    subplotLabel(ax)
    optimizeDesign([ax[0], ax[3]], ["Treg"], ["Thelper", "NK", "CD8"])
    optimizeDesign([ax[1], ax[4]], ["NK"], ["Thelper", "Treg", "CD8"], legend=True)
    optimizeDesign([ax[2], ax[5]], ["Thelper"], ["Treg", "NK", "CD8"], IL7=True, legend=True)

    return f


cellTypeDict = {"Treg": r"$T_{reg}$",
                "Thelper": r"$T_{helper}$",
                "NK": "NK",
                "CD8": r"$CD8^{+}$"}


cell_order = ["$T_{reg}$", "$T_{helper}$", "NK", "$CD8^{+}$"]


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
    vals = np.arange(1, 10, step=0.25)
    sigDF = pd.DataFrame()

    if IL7:
        vals = np.arange(1.1, 10.1, step=0.25)
        optDF = pd.DataFrame(columns={"Valency", "Selectivity", "IL7Rα"})
        X0 = [8]
        optBnds = [(6, 11)]  # Ka IL7, Kx
    else:
        optDF = pd.DataFrame(columns={"Valency", "Selectivity", "IL2Rα", r"$IL-2Rβ/γ_c$"})
        X0 = [8, 8]
        optBnds = [(6, 11),  # Ka IL2Ra
                   (6, 11)]  # Ka IL2Rb

    for i, val in enumerate(vals):
        if i == 0:
            optimized = minimize(minSelecFunc, X0, bounds=optBnds, args=(val, targCell, offTcells, IL7), jac="3-point")
            initialTargBind = cytBindingModelOpt(optimized.x, val, targCell[0], IL7)

            def bindConstraintFun(x):
                return cytBindingModelOpt(x, val, targCell[0], IL7) - 0.5 * initialTargBind

            bindConst = {'type': 'ineq', 'fun': bindConstraintFun}

        else:
            optimized = minimize(minSelecFunc, X0, bounds=optBnds, args=(val, targCell, offTcells, IL7), jac="3-point", constraints=bindConst)
        if IL7:
            IL7RaKD = 1e9 / np.power(10, optimized.x[0])
            optDF = optDF.append(pd.DataFrame({"Valency": [val], "Selectivity": [len(offTcells) / optimized.fun], "IL7Rα": IL7RaKD}))
            sigDF = sigDF.append(pd.DataFrame({"Cell Type": [targCell[0]], "Target": ["Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, targCell[0], IL7)]}))
            for cell in offTcells:
                sigDF = sigDF.append(pd.DataFrame({"Cell Type": [cell], "Target": ["Off-Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, cell, IL7)]}))
        else:
            IL2RaKD = 1e9 / np.power(10, optimized.x[0])
            IL2RBGKD = 1e9 / np.power(10, optimized.x[1])
            optDF = optDF.append(pd.DataFrame({"Valency": [val], "Selectivity": [len(offTcells) / optimized.fun], "IL2Rα": IL2RaKD, r"$IL-2Rβ/γ_c$": IL2RBGKD}))
            sigDF = sigDF.append(pd.DataFrame({"Cell Type": [targCell[0]], "Target": ["Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, targCell[0], IL7)]}))
            for cell in offTcells:
                sigDF = sigDF.append(pd.DataFrame({"Cell Type": [cell], "Target": ["Off-Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, cell, IL7)]}))
    # Normalize to valency 1
    for cell in targCell + offTcells:
        sigDF.loc[sigDF["Cell Type"] == cell, "pSTAT"] = sigDF.loc[sigDF["Cell Type"] == cell, "pSTAT"].div(sigDF.loc[(sigDF["Cell Type"] == cell) & (sigDF.Valency == vals[0])].pSTAT.values[0])

    sigDF = sigDF.replace(cellTypeDict)

    if IL7:
        sns.lineplot(x="Valency", y="pSTAT", hue="Cell Type", style="Target", data=sigDF, ax=ax[0], palette="husl", hue_order=cell_order)
        ax[0].set_title(cellTypeDict[targCell[0]] + " selectivity with IL-7 mutein", fontsize=7)
        ax[0].set_ylim(bottom=0.0, top=4)

        sns.lineplot(x="Valency", y="IL7Rα", data=optDF, ax=ax[1], palette="crest")
        ax[1].set(yscale="log", ylabel=r"IL7·7Rα $K_D$ (nM)")
        ax[1].set_yticks([10, np.power(10, 2), np.power(10, 3)])

    else:
        sns.lineplot(x="Valency", y="pSTAT", hue="Cell Type", style="Target", data=sigDF, ax=ax[0], palette="husl", hue_order=cell_order)
        ax[0].set_title(cellTypeDict[targCell[0]] + " selectivity with IL-2 mutein", fontsize=7)
        ax[0].set_ylim(bottom=0.0, top=3)

        if targCell == ["NK"]:
            affDF = pd.melt(optDF, id_vars=['Valency'], value_vars=[r"$IL-2Rβ/γ_c$"])
            affDF = affDF.rename(columns={"variable": "Receptor"})
            sns.lineplot(x="Valency", y="value", hue="Receptor", data=affDF, ax=ax[1])
            ax[1].set(yscale="log", ylabel=r"$IL2·β/γ_c$ $K_D$ (nM)")

        else:
            affDF = pd.melt(optDF, id_vars=['Valency'], value_vars=['IL2Rα', r"$IL-2Rβ/γ_c$"])
            affDF = affDF.rename(columns={"variable": "Receptor"})
            sns.lineplot(x="Valency", y="value", hue="Receptor", data=affDF, ax=ax[1])
            ax[1].set(yscale="log", ylabel=r"IL2· $K_D$ (nM)")

    ax[0].set_xticks(np.arange(1, 11))
    ax[1].set_xticks(np.arange(1, 11))
    if not legend:
        ax[0].get_legend().remove()

    return optimized
