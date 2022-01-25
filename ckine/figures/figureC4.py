"""
This creates Figure 4, fitting of multivalent binding model to Gc Data.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import minimize
from copy import copy
from .figureCommon import subplotLabel, getSetup
from ..PCA import nllsq_EC50
from ..MBmodel import getKxStar, runFullModel, cytBindingModel, polyc
from ..imports import getBindDict, importReceptors

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 9), (4, 4), multz={9: 1})
    axlabel = copy(ax)
    del axlabel[1]
    del axlabel[1]
    del axlabel[1]
    del axlabel[2]
    subplotLabel(axlabel)

    ax[0].axis("off")
    ax[1].axis("off")
    ax[2].axis("off")
    ax[3].axis("off")
    ax[5].axis("off")

    minSolved = minimize(runFullModel, x0=-12.0, args=([0.5, 1], False, True))
    print(minSolved)

    modelDF = runFullModel(time=[0.5, 1.0], saveDict=False, singleCell=True)  # Change to save

    print(r2_score(modelDF.Experimental.values, modelDF.Predicted.values))
    Pred_Exp_plot(ax[4], modelDF)
    legend = ax[4].get_legend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[5].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 8.5})  # use this to place universal legend later
    ax[4].get_legend().remove()

    R2_Plot_Cells(ax[6], modelDF)
    R2_Plot_Ligs(ax[7], modelDF)
    MonVsBivalent(ax[8], modelDF, ligs=True)

    R2_Plot_Conc(ax[9], modelDF)
    timePlot(ax[10])

    IL2RaEffPlot(ax[11], modelDF, "Treg", IL2RBaff=2e8, IL2Ra_affs=np.array([1e8, 1e9, 1e10]))
    IL2RaEffPlot(ax[12], modelDF, "NK", IL2RBaff=2e8, IL2Ra_affs=np.array([1e8, 1e9, 1e10]))
    recSigPlot(ax[13], modelDF, IL2RBrec=1000, IL2Rarecs=[100, 1000, 10000], IL2RBaff=2e8, IL2Ra_aff=1e8)
    recSigPlot(ax[14], modelDF, IL2RBrec=1000, IL2Rarecs=[100, 1000, 10000], IL2RBaff=2e8, IL2Ra_aff=1e10)
    return f


def Pred_Exp_plot(ax, df):
    """Plots all experimental vs. Predicted Values"""
    sns.scatterplot(x="Experimental", y="Predicted", hue="Cell", style="Valency", data=df, ax=ax, alpha=0.35)
    ax.set(xlim=(0, 60000), ylim=(0, 60000))


def R2_Plot_Cells(ax, df):
    """Plots all accuracies per cell"""
    accDF = pd.DataFrame(columns={"Cell Type", "Valency", "Accuracy"})
    cellTypes = ['Treg $IL2Ra^{hi}$', 'Treg', 'Treg $IL2Ra^{lo}$', 'Thelper $IL2Ra^{hi}$', 'Thelper', 'Thelper $IL2Ra^{lo}$', 'CD8', 'NK']
    for cell in cellTypes:
        for val in df.Valency.unique():
            preds = df.loc[(df.Cell == cell) & (df.Valency == val)].Predicted.values
            exps = df.loc[(df.Cell == cell) & (df.Valency == val)].Experimental.values
            r2 = r2_score(exps, preds)
            accDF = accDF.append(pd.DataFrame({"Cell Type": [cell], "Valency": [val], "Accuracy": [r2]}))

    sns.barplot(x="Cell Type", y="Accuracy", hue="Valency", data=accDF, ax=ax)
    ax.set(ylim=(0, 1), ylabel=r"Accuracy ($R^2$)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")


def R2_Plot_Ligs(ax, df):
    """Plots all accuracies per ligand"""
    accDF = pd.DataFrame(columns={"Ligand", "Valency", "Accuracy"})
    for ligand in df.Ligand.unique():
        for val in df.loc[df.Ligand == ligand].Valency.unique():
            preds = df.loc[(df.Ligand == ligand) & (df.Valency == val)].Predicted.values
            exps = df.loc[(df.Ligand == ligand) & (df.Valency == val)].Experimental.values
            r2 = r2_score(exps, preds)
            if val == 1:
                accDF = accDF.append(pd.DataFrame({"Ligand": [ligand + " (Mono)"], "Valency": [val], "Accuracy": [r2]}))
            else:
                accDF = accDF.append(pd.DataFrame({"Ligand": [ligand + " (Biv)"], "Valency": [val], "Accuracy": [r2]}))
    sns.barplot(x="Ligand", y="Accuracy", data=accDF, ax=ax)
    ax.set(ylim=(0, 1), ylabel=r"Accuracy ($R^2$)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")


def R2_Plot_Conc(ax, df):
    """Plots all accuracies per concentration"""
    accDF = pd.DataFrame(columns={"Concentration", "Valency", "Accuracy"})
    for conc in df.Dose.unique():
        for val in df.loc[(df.Dose == conc)].Valency.unique():
            preds = df.loc[(df.Dose == conc) & (df.Valency == val)].Predicted.values
            exps = df.loc[(df.Dose == conc) & (df.Valency == val)].Experimental.values
            r2 = r2_score(exps, preds)
            accDF = accDF.append(pd.DataFrame({"Concentration": [conc], "Valency": [val], "Accuracy": [r2]}))
    sns.scatterplot(x="Concentration", y="Accuracy", hue="Valency", data=accDF, ax=ax)
    ax.set(ylim=(0, 1), ylabel=r"Accuracy ($R^2$)", xlabel="Dose (nM)", xscale="log")


def MonVsBivalent(ax, dfAll, ligs=True):
    """Compares accuracy of monovalent vs bivalent predictions"""
    df = dfAll.loc[(dfAll.Valency == 2)].copy()
    dates = df.Date.unique()
    df["MonPredict"] = np.zeros(df.shape[0])

    for date in dates:
        dfDate = df.loc[(df.Date == date)]
        ligands = dfDate.Ligand.unique()
        concs = dfDate.Dose.unique()
        cellTypes = dfDate.Cell.unique()
        times = dfDate.Time.unique()

        for lig in ligands:
            for conc in concs:
                for cell in cellTypes:
                    predVal = cytBindingModel(lig, 1, conc * 2, cell)
                    for time in times:
                        df.loc[(df.Date == date) & (df.Ligand == lig) & (df.Dose == conc) & (df.Cell == cell) & (df.Time == time), "MonPredict"] = predVal
    dfCopy = copy(df)

    for date in dates:
        for cell in cellTypes:
            if cell[-1] == "$":  # if it is a binned pop, use ave fit
                predVecBin = df.loc[(df.Date == date) & (df.Cell == cell)].MonPredict.values
                expVec = dfCopy.loc[(dfCopy.Date == date) & (dfCopy.Cell == cell[0:-13])].Experimental.values
                predVec = dfCopy.loc[(dfCopy.Date == date) & (dfCopy.Cell == cell[0:-13])].MonPredict.values
                slope = np.linalg.lstsq(np.reshape(predVec, (-1, 1)), np.reshape(expVec, (-1, 1)), rcond=None)[0][0]
                df.loc[(df.Date == date) & (df.Cell == cell), "MonPredict"] = predVecBin * slope
            else:
                expVec = df.loc[(df.Date == date) & (df.Cell == cell)].Experimental.values
                predVec = df.loc[(df.Date == date) & (df.Cell == cell)].MonPredict.values
                slope = np.linalg.lstsq(np.reshape(predVec, (-1, 1)), np.reshape(expVec, (-1, 1)), rcond=None)[0][0]
                df.loc[(df.Date == date) & (df.Cell == cell), "MonPredict"] = predVec * slope

    if ligs:
        accDF = pd.DataFrame(columns={"Ligand", "Prediction Label", "Accuracy"})
        for ligand in df.Ligand.unique():
            BivPreds = df.loc[(df.Ligand == ligand)].Predicted.values
            MonPreds = df.loc[(df.Ligand == ligand)].MonPredict.values
            exps = df.loc[(df.Ligand == ligand)].Experimental.values
            r2Biv = r2_score(exps, BivPreds)
            r2Mon = r2_score(exps, MonPreds)
            accDF = accDF.append(pd.DataFrame({"Ligand": [ligand], "Prediction Label": "Monovalent (Incorrect)", "Accuracy": [r2Mon]}))
            accDF = accDF.append(pd.DataFrame({"Ligand": [ligand], "Prediction Label": "Bivalent (Correct)", "Accuracy": [r2Biv]}))
        sns.barplot(x="Ligand", y="Accuracy", hue="Prediction Label", data=accDF, ax=ax)
        ax.set(ylim=(0, 1), ylabel=r"Accuracy ($R^2$)")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    else:
        accDF = pd.DataFrame(columns={"Cell Type", "Prediction Label", "Accuracy"})
        for cellType in df.Cell.unique():
            BivPreds = df.loc[(df.Cell == cellType)].Predicted.values
            MonPreds = df.loc[(df.Cell == cellType)].MonPredict.values
            exps = df.loc[(df.Cell == cellType)].Experimental.values
            r2Biv = r2_score(exps, BivPreds)
            r2Mon = r2_score(exps, MonPreds)
            accDF = accDF.append(pd.DataFrame({"Cell Type": [cellType], "Prediction Label": [1], "Accuracy": [r2Mon]}))
            accDF = accDF.append(pd.DataFrame({"Cell Type": [cellType], "Prediction Label": [2], "Accuracy": [r2Biv]}))
        sns.barplot(x="Cell Type", y="Accuracy", hue="Prediction Label", data=accDF, ax=ax)
        ax.set(ylim=(0, 1), ylabel=r"Accuracy ($R^2$)")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

    return df


def EC50comp(ax, dfAll, time):
    """Predicts EC50s for each drug for 4 cell types, and plots"""
    x0exp = [4, 2.0, 1000.0]
    x0pred = [4, 2.0, 1000.0]
    df = dfAll.loc[(dfAll.Time == time)].copy()
    Ligands = df.Ligand.unique()
    Cells = df.Cell.unique()

    EC50df = pd.DataFrame(columns=["Cell Type", "Ligand", "EC50", "Exp/Pred"])

    for ligand in Ligands:
        Valencies = df.loc[(df.Ligand == ligand), "Valency"].unique()
        for valency in Valencies:
            for cell in Cells:
                #dates = df.loc[(df.Ligand == ligand) & (df.Cell == cell) & (df.Valency == valency)].Date.values
                dosesExp = df.loc[(df.Ligand == ligand) & (df.Cell == cell) & (df.Valency == valency)].Dose.values
                doseMax, doseMin = np.log10(np.amax(dosesExp)) + 4, np.log10(np.amin(dosesExp))
                dosesPredMB = np.logspace(doseMin, doseMax, 40)
                dosesPred = np.log10(dosesPredMB) + 4
                dosesExp = np.log10(dosesExp) + 4

                expVals = df.loc[(df.Ligand == ligand) & (df.Cell == cell) & (df.Valency == valency)].Experimental.values
                predVals = cytBindingModel(ligand, valency, dosesPredMB, cell)
                EC50exp = nllsq_EC50(x0exp, dosesExp, expVals) - 4
                EC50pred = nllsq_EC50(x0pred, dosesPred, predVals) - 4

                if valency == 1:
                    EC50df = EC50df.append(pd.DataFrame({"Cell Type": [cell], "Ligand": [ligand + " (Mono)"], "EC50": [EC50exp], "Exp/Pred": ["Experimental"]}))
                    EC50df = EC50df.append(pd.DataFrame({"Cell Type": [cell], "Ligand": [ligand + " (Mono)"], "EC50": [EC50pred], "Exp/Pred": ["Predicted"]}))
                else:
                    EC50df = EC50df.append(pd.DataFrame({"Cell Type": [cell], "Ligand": [ligand + " (Biv)"], "EC50": [EC50exp], "Exp/Pred": ["Experimental"]}))
                    EC50df = EC50df.append(pd.DataFrame({"Cell Type": [cell], "Ligand": [ligand + " (Biv)"], "EC50": [EC50pred], "Exp/Pred": ["Predicted"]}))

    EC50df = EC50df.loc[(EC50df["Cell Type"].isin(["Treg", "Thelper"]))]
    sns.scatterplot(x="Ligand", y="EC50", hue="Cell Type", style="Exp/Pred", data=EC50df, ax=ax)
    ax.set(ylabel=r"log$_{10}$EC50 (nM)", ylim=(-2, 6))
    ax.set_xticklabels(EC50df.Ligand.unique(), rotation=45)


def timePlot(ax):
    """Plots all experimental vs. Predicted Values"""
    times = [[0.5], [1.], [2.], [4.]]
    accDF = pd.DataFrame(columns={"Time", "Valency", "Accuracy"})
    for time in times:
        df = runFullModel(time=time, saveDict=False)
        for val in df.Valency.unique():
            preds = df.loc[(df.Time == time[0]) & (df.Valency == val)].Predicted.values
            exps = df.loc[(df.Time == time[0]) & (df.Valency == val)].Experimental.values
            r2 = r2_score(exps, preds)
            accDF = accDF.append(pd.DataFrame({"Time": time, "Valency": [val], "Accuracy": [r2]}))
    sns.barplot(x="Time", y="Accuracy", hue="Valency", data=accDF, ax=ax)
    ax.set(ylim=(0, 1), ylabel=r"Accuracy ($R^2$)")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)


def IL2RaEffPlot(ax, modelDF, cell, IL2RBaff, IL2Ra_affs):
    """Tracks interactions of Valency and ligand affinity on Treg activation"""
    concs = np.logspace(start=np.log10(modelDF.Dose.min()), stop=np.log10(modelDF.Dose.max()), num=101, endpoint=True)
    bindDict = getBindDict().groupby(["Cell"])['Scale'].mean().reset_index()
    convFact = bindDict.loc[bindDict.Cell == cell].Scale.values
    recDF = importReceptors()
    recCount = np.ravel([recDF.loc[(recDF.Receptor == "IL2Ra") & (recDF["Cell Type"] == cell)].Mean.values,
                         recDF.loc[(recDF.Receptor == "IL2Rb") & (recDF["Cell Type"] == cell)].Mean.values])
    valencies = [1, 2]
    outputDF = pd.DataFrame()

    for val in valencies:
        for alphAff in IL2Ra_affs:
            for dose in concs:
                Affs = np.array([alphAff, IL2RBaff])
                Affs = np.reshape(Affs, (1, -1))
                Affs = np.repeat(Affs, 2, axis=0)
                np.fill_diagonal(Affs, 1e2)  # Each cytokine can only bind one a and one b
                predVal = polyc(dose / 1e9 / val, getKxStar(), recCount, [[val, val]], [1.0], Affs)[0][1] * convFact
                alphAffKDnM = (1 / alphAff) / 1e-9
                outputDF = outputDF.append(pd.DataFrame({r"IL2Rα $K_D$ (nM)": [alphAffKDnM], "Concentration": [dose], "Valency": [val], "Predicted pSTAT5 MFI": predVal}))

    outputDF = outputDF.reset_index()
    sns.lineplot(data=outputDF, x="Concentration", y="Predicted pSTAT5 MFI", style="Valency", hue=r"IL2Rα $K_D$ (nM)", ax=ax)
    ax.set(xscale="log")


def recSigPlot(ax, modelDF, IL2RBrec, IL2Rarecs, IL2RBaff, IL2Ra_aff):
    """Tracks interactions of Valency and ligand affinity on Treg activation"""
    concs = np.logspace(start=np.log10(modelDF.Dose.min()), stop=np.log10(modelDF.Dose.max()), num=101, endpoint=True)
    valencies = [1, 2, 4]
    outputDF = pd.DataFrame()
    Affs = np.array([IL2Ra_aff, IL2RBaff])
    Affs = np.reshape(Affs, (1, -1))
    Affs = np.repeat(Affs, 2, axis=0)
    np.fill_diagonal(Affs, 1e2)  # Each cytokine can only bind one a and one b

    for val in valencies:
        for alphaRec in IL2Rarecs:
            for dose in concs:
                recCount = np.ravel([alphaRec,
                                    IL2RBrec])
                predVal = polyc(dose / 1e9 / val, getKxStar(), recCount, [[val, val]], [1.0], Affs)[0][1]
                outputDF = outputDF.append(pd.DataFrame({r"IL2Rα Abundance": [alphaRec], "Concentration": [dose], "Valency": [val], "Active Binding Complexes": predVal}))

    outputDF = outputDF.reset_index()
    print(outputDF)
    sns.lineplot(data=outputDF, x="Concentration", y="Active Binding Complexes", style="Valency", hue=r"IL2Rα Abundance", ax=ax)
    ax.set(xscale="log")
