"""
Figure 6. Optimization of Ligands
"""
from os.path import dirname, join
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import minimize
from .figureCommon import subplotLabel, getSetup
from ..MBmodel import polyc

path_here = dirname(dirname(__file__))


def makeFigure():
    """ Make figure 6. """
    # Get list of axis objects
    ax, f = getSetup((5, 6), (3, 2))
    subplotLabel(ax)
    optimizeDesign(ax[0:2], ["Treg"], ["Thelper", "NK", "CD8"])
    optimizeDesign(ax[2:4], ["NK"], ["Thelper", "Treg", "CD8"])
    optimizeDesign(ax[4:6], ["Thelper"], ["Treg", "NK", "CD8"], IL7=True)

    return f


def cytBindingModelOpt(x, val, cellType, IL7=False):
    """Runs binding model for a given mutein, valency, dose, and cell type. """
    recQuantDF = pd.read_csv(join(path_here, "data/RecQuantitation.csv"))
    convDict = pd.read_csv(join(path_here, "data/BindingConvDict.csv"))
    recCount = recQuantDF[["Receptor", cellType]]
    Kx = np.power(10, -12.0)  # TODO: Place fitting Kx

    if IL7:
        affs = [[np.power(10, x[0])]]
        recCount = [recCount.loc[(recCount.Receptor == "IL7Ra")][cellType].values]
        recCount = np.ravel(np.power(10, recCount))
        output = polyc(1e-9 / val, Kx, recCount, [[val]], [1.0], affs)[1][0][0]  # IL7Ra binding only
    else:
        affs = [[np.power(10, x[0]), 1e2], [1e2, np.power(10, x[1])]]
        recCount = [recCount.loc[(recCount.Receptor == "IL2Ra")][cellType].values, recCount.loc[(recCount.Receptor == "IL2Rb")][cellType].values]
        recCount = np.ravel(np.power(10, recCount))
        output = polyc(1e-9 / val, Kx, recCount, [[val, val]], [1.0], affs)[1][0][1]  # IL2RB binding only

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


def optimizeDesign(ax, targCell, offTcells, IL7=False):
    """ A more general purpose optimizer """
    vals = np.arange(1, 9, step=1)
    sigDF = pd.DataFrame()

    if IL7:
        optDF = pd.DataFrame(columns={"Valency", "Selectivity", "IL7Ra"})
        X0 = [8]
        optBnds = [(6, 11)]  # Ka IL7, Kx
    else:
        optDF = pd.DataFrame(columns={"Valency", "Selectivity", "IL2Ra", "IL2RBG"})
        X0 = [8, 8]
        optBnds = [(6, 11),  # Ka IL2Ra
                   (6, 11)]  # Ka IL2Rb

    for val in vals:
        optimized = minimize(minSelecFunc, X0, bounds=optBnds, args=(val, targCell, offTcells, IL7), jac="3-point")
        print(val)
        print(optimized.fun)
        if IL7:
            IL7RaKD = 1e9 / np.power(10, optimized.x[0])
            optDF = optDF.append(pd.DataFrame({"Valency": [val], "Selectivity": [len(offTcells) / optimized.fun], "IL7Ra": IL7RaKD}))
            sigDF = sigDF.append(pd.DataFrame({"Cell Type": [targCell[0]], "Target": ["Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, targCell[0], IL7)]}))
            for cell in offTcells:
                sigDF = sigDF.append(pd.DataFrame({"Cell Type": [cell], "Target": ["Off-Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, cell, IL7)]}))
        else:
            IL2RaKD = 1e9 / np.power(10, optimized.x[0])
            IL2RBGKD = 1e9 / np.power(10, optimized.x[1])
            optDF = optDF.append(pd.DataFrame({"Valency": [val], "Selectivity": [len(offTcells) / optimized.fun], "IL2Ra": IL2RaKD, "IL2RBG": IL2RBGKD}))
            sigDF = sigDF.append(pd.DataFrame({"Cell Type": [targCell[0]], "Target": ["Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, targCell[0], IL7)]}))
            for cell in offTcells:
                sigDF = sigDF.append(pd.DataFrame({"Cell Type": [cell], "Target": ["Off-Target"], "Valency": [val], "pSTAT": [cytBindingModelOpt(optimized.x, val, cell, IL7)]}))

    # Normalize to valency 1
    for cell in targCell + offTcells:
        sigDF.loc[sigDF["Cell Type"] == cell, "pSTAT"] = sigDF.loc[sigDF["Cell Type"] == cell, "pSTAT"].div(sigDF.loc[(sigDF["Cell Type"] == cell) & (sigDF.Valency == 1)].pSTAT.values[0])

    if IL7:
        sns.lineplot(x="Valency", y="pSTAT", hue="Cell Type", style="Target", data=sigDF, ax=ax[0], palette="husl")
        ax[0].set(title=targCell[0] + " Selectivity with IL-7 mutein")
        ax[0].set_ylim(bottom=0.0)

        sns.lineplot(x="Valency", y="IL7Ra", data=optDF, ax=ax[1], palette="crest")
        ax[1].set(yscale="log", ylabel=r"IL7·7Rα $K_D$ (nM)")

    else:
        sns.lineplot(x="Valency", y="pSTAT", hue="Cell Type", style="Target", data=sigDF, ax=ax[0], palette="husl")
        ax[0].set(title=targCell[0] + " Selectivity with IL-2 mutein")
        ax[0].set_ylim(bottom=0.0)

        if targCell == ["NK"]:
            affDF = pd.melt(optDF, id_vars=['Valency'], value_vars=['IL2RBG'])
            affDF = affDF.rename(columns={"variable": "Receptor"})
            sns.lineplot(x="Valency", y="value", hue="Receptor", data=affDF, ax=ax[1])
            ax[1].set(yscale="log", ylabel=r"IL2· $K_D$ (nM)")
        else:
            affDF = pd.melt(optDF, id_vars=['Valency'], value_vars=['IL2Ra', 'IL2RBG'])
            affDF = affDF.rename(columns={"variable": "Receptor"})
            sns.lineplot(x="Valency", y="value", hue="Receptor", data=affDF, ax=ax[1])
            ax[1].set(yscale="log", ylabel=r"IL2· $K_D$ (nM)")

    return optimized
