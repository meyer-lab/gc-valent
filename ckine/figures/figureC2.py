"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.optimize import least_squares
from os.path import join
from sklearn.preprocessing import StandardScaler, LabelEncoder, LabelBinarizer
from sklearn.linear_model import RidgeClassifierCV
from scipy import stats
from sklearn.neighbors import KernelDensity
from sklearn.svm import SVC
from copy import copy
from .figureCommon import subplotLabel, getSetup, getLigDict, get_cellTypeDict, getLigandLegend
from ..imports import import_pstat_all, importCITE


path_here = os.path.dirname(os.path.dirname(__file__))
ligDict = getLigDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((14, 15), (5, 4), multz={0: 2, 4: 1, 8: 1, 12: 1})
    subplotLabel(ax)
    ax[0].axis("off")
    ax[1].axis("off")

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    mutAffDF = mutAffDF.rename({"Mutein": "Ligand", "IL2RaKD": "IL2Rα $K_{D}$ (nM)", "IL2RBGKD": "IL2Rβ $K_{D}$ (nM)"}, axis=1)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all(True, False)
    respDF = respDF.replace(cellDict)
    respDF = respDF.rename({"Bivalent": "Valency"}, axis=1)
    respDF["Valency"] = respDF["Valency"] + 1
    respDF = respDF.loc[(respDF.Ligand != "IL15") & (respDF.Ligand != "IL2")]
    mutAffDF = mutAffDF.loc[(mutAffDF.Ligand != "IL15") & (mutAffDF.Ligand != "IL2")]

    doses = respDF.Dose.unique()
    ratioConc(ax[2:5], respDF, r"T$_{reg}$", "NK", 4, mutAffDF, legend=True)
    ratioConc(ax[5:8], respDF, r"T$_{reg}$", r"CD8$^{+}$", 4, mutAffDF, legend=True)
    ratioConc(ax[8:11], respDF, r"T$_{reg}$", r"T$_{helper}$", 4, mutAffDF, legend=True)

    legend = getLigandLegend()
    labels = (x.get_text() for x in legend.get_texts())
    ax[1].legend(legend.legendHandles, labels, loc="upper left", prop={"size": 10})  # use this to place universal legend later
    cellTarget = "Treg"
    Wass_KL_Dist(ax[11:13], cellTarget, 10)
    CITE_RIDGE(ax[13], cellTarget)
    CITE_SVM(ax[14], cellTarget, sampleFrac=0.2)

    return f


def pSTATcomp(ax, respDF, mutAffDF, recAff, cell, time, conc, legend=False):
    """Plots monovalent and bivalent lines for a single cell type by affinity"""
    respDF = respDF.loc[(respDF.Cell == cell) & ((respDF.Time == time) & (respDF.Dose == conc))]
    respDF = respDF.groupby(["Ligand", "Valency"]).Mean.mean().reset_index()
    respDF = respDF.merge(mutAffDF)
    respDF = respDF.loc[(respDF.Ligand != "IL15") & (respDF.Ligand != "IL2")]
    sns.scatterplot(data=respDF, x=recAff, y="Mean", hue="Ligand", style="Valency", ax=ax, legend=legend, palette=ligDict)
    ax.set(xscale="log", title=cell + " at " + str(conc) + " nM  and " + str(time) + " hours")


def gaussian(x, concs):
    """Returns Gaussian for given concs"""
    amplitude, mean, stddev = x[0], x[1], x[2]
    return amplitude * np.exp(-((concs - mean) / 4 / stddev)**2)


def gaussian_residuals(x, concs, ratios):
    """ Residual function for Gaussian. """
    return gaussian(x, concs) - ratios


def ratioConc(ax, respDF, cell1, cell2, time, mutAffDF, legend=False):
    """Plots Ratio of cell 1 to cell 2 over a range of concentrations"""
    respDF1 = respDF.loc[(respDF.Cell == cell1) & ((respDF.Time == time))]
    respDF2 = respDF.loc[(respDF.Cell == cell2) & ((respDF.Time == time))]
    respDF1 = respDF1.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    respDF2 = respDF2.groupby(["Ligand", "Valency", "Cell", "Dose"]).Mean.mean().reset_index()
    respDF1 = respDF1.pivot(index=["Ligand", "Dose", "Valency"], columns="Cell", values="Mean").reset_index()
    respDF2 = respDF2.pivot(index=["Ligand", "Dose", "Valency"], columns="Cell", values="Mean").reset_index()
    ratioDF = respDF1.merge(respDF2)
    ratioDF["Ratio"] = ratioDF[cell1] / ratioDF[cell2]
    ratioDF = ratioDF.dropna()

    ratioDF.loc[ratioDF.Dose < np.sort(ratioDF.Dose.unique())[3], "Ratio"] = 1
    ratioDF = ratioDF.loc[(ratioDF.Ratio <= 3000) & (ratioDF.Ratio > 0)]
    ratioDF["Ratio"] = np.log10(ratioDF["Ratio"].values)

    doses = np.log10(np.logspace(np.log10(ratioDF.Dose.min()), np.log10(ratioDF.Dose.max()), 100))
    x0 = [3, 0, 2]
    gaussDF = pd.DataFrame()
    fitDF = pd.DataFrame()
    for ligand in ratioDF.Ligand.unique():
        for valency in ratioDF.loc[ratioDF.Ligand == ligand].Valency.unique():
            isoData = ratioDF.loc[(ratioDF.Ligand == ligand) & (ratioDF.Valency == valency)]
            xData = np.nan_to_num(np.log10(isoData.Dose.values))
            yData = np.nan_to_num(isoData.Ratio.values)
            if xData.size > 5:
                fit = least_squares(gaussian_residuals, x0, args=(xData, yData), bounds=([0, -4, 0], [3.5, 2, 6]), jac="3-point")
                gaussDF = gaussDF.append(pd.DataFrame({"Ligand": ligand, "Valency": valency, "Dose": np.power(10, doses), "Ratio": np.power(10, gaussian(fit.x, doses))}))
                fitDF = fitDF.append(pd.DataFrame({"Ligand": [ligand], "Valency": valency, cell2 + " Max": np.power(10, fit.x[0]), cell2 + " Dose": np.power(10, fit.x[1])}))

    gaussDF = gaussDF.reset_index()
    ratioDF["Ratio"] = np.power(10, ratioDF["Ratio"].values)

    sns.scatterplot(data=ratioDF, x="Dose", y="Ratio", hue="Ligand", style="Valency", ax=ax[0], palette=ligDict, legend=False)
    sns.lineplot(data=gaussDF, x="Dose", y="Ratio", hue="Ligand", style="Valency", ax=ax[0], palette=ligDict, legend=legend)
    ax[0].set(xscale="log", title="Ratio of " + cell1 + " to " + cell2)

    if legend:
        h, l = ax[0].get_legend_handles_labels()
        ax[0].legend(h[-3:], l[-3:])

    fitDF = fitDF.merge(mutAffDF)

    maxLineDF = pd.DataFrame()
    doseLineDF = pd.DataFrame()
    affs = np.linspace(-1, 1, 100).flatten()
    for valency in fitDF.Valency.unique():
        valData = fitDF.loc[fitDF.Valency == valency]
        mMax, bMax = np.polyfit(np.log10(valData["IL2Rα $K_{D}$ (nM)"].values), valData[cell2 + " Max"], 1)
        mDose, bDose = np.polyfit(np.log10(valData["IL2Rα $K_{D}$ (nM)"].values), np.log10(valData[cell2 + " Dose"]), 1)
        maxLineDF = maxLineDF.append(pd.DataFrame({"Valency": valency, "IL2Rα $K_{D}$ (nM)": np.power(10, affs), cell2 + " Max": mMax * affs + bMax}))
        doseLineDF = doseLineDF.append(pd.DataFrame({"Valency": valency, "IL2Rα $K_{D}$ (nM)": np.power(10, affs), cell2 + " Dose": np.power(10, mDose * affs + bDose)}))

    maxLineDF, doseLineDF = maxLineDF.reset_index(), doseLineDF.reset_index()
    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", hue="Ligand", style="Valency", ax=ax[1], palette=ligDict, legend=False)
    sns.lineplot(data=maxLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Max", style="Valency", ax=ax[1], color="k", linewidth=1., legend=False)
    ax[1].set(xscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(0, None))

    sns.scatterplot(data=fitDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", hue="Ligand", style="Valency", ax=ax[2], palette=ligDict, legend=False)
    sns.lineplot(data=doseLineDF, x="IL2Rα $K_{D}$ (nM)", y=cell2 + " Dose", style="Valency", ax=ax[2], color="k", linewidth=1., legend=False)
    ax[2].set(xscale="log", yscale="log", title="Ratio of " + cell1 + " to " + cell2, xlim=(1e-1, 1e1), ylim=(1e-2, 1e2))


def Wass_KL_Dist(ax, targCell, numFactors):
    """Finds markers which have average greatest difference from other cells"""
    CITE_DF = importCITE()

    markerDF = pd.DataFrame(columns=["Marker", "Cell Type", "Amount"])
    for marker in CITE_DF.loc[:, ((CITE_DF.columns != 'CellType1') & (CITE_DF.columns != 'CellType2') & (CITE_DF.columns != 'CellType3') & (CITE_DF.columns != 'Cell'))].columns:
        markAvg = np.mean(CITE_DF[marker].values)
        targCellMark = CITE_DF.loc[CITE_DF["CellType2"] == targCell][marker].values / markAvg
        offTargCellMark = CITE_DF.loc[CITE_DF["CellType2"] != targCell][marker].values / markAvg
        if np.mean(targCellMark) > np.mean(offTargCellMark):
            kdeTarg = KernelDensity(kernel='gaussian').fit(targCellMark.reshape(-1, 1))
            kdeOffTarg = KernelDensity(kernel='gaussian').fit(offTargCellMark.reshape(-1, 1))
            minVal = np.minimum(targCellMark.min(), offTargCellMark.min()) - 10
            maxVal = np.maximum(targCellMark.max(), offTargCellMark.max()) + 10
            outcomes = np.arange(minVal, maxVal + 1).reshape(-1, 1)
            distTarg = np.exp(kdeTarg.score_samples(outcomes))
            distOffTarg = np.exp(kdeOffTarg.score_samples(outcomes))
            KL_div = stats.entropy(distOffTarg.flatten() + 1e-200, distTarg.flatten() + 1e-200, base=2)
            markerDF = markerDF.append(pd.DataFrame({"Marker": [marker], "Wasserstein Distance": stats.wasserstein_distance(targCellMark, offTargCellMark), "KL Divergence": KL_div}))

    for i, distance in enumerate(["Wasserstein Distance", "KL Divergence"]):
        ratioDF = markerDF.sort_values(by=distance)
        posCorrs = ratioDF.tail(numFactors).Marker.values

        markerDF = markerDF.loc[markerDF["Marker"].isin(posCorrs)]

        sns.barplot(data=ratioDF.tail(numFactors), x="Marker", y=distance, ax=ax[i], color='k')
        ax[i].set(yscale="log")
        ax[i].set_xticklabels(ax[i].get_xticklabels(), rotation=45)


def CITE_RIDGE(ax, targCell, numFactors=10):
    """Fits a ridge classifier to the CITE data and plots those most highly correlated with T reg"""
    ridgeMod = RidgeClassifierCV()
    RIDGE_DF = importCITE()
    cellToI = RIDGE_DF.CellType2.unique()
    RIDGE_DF = RIDGE_DF.loc[(RIDGE_DF["CellType2"].isin(cellToI)), :]
    cellTypeCol = RIDGE_DF.CellType2.values
    RIDGE_DF = RIDGE_DF.loc[:, ((RIDGE_DF.columns != 'CellType1') & (RIDGE_DF.columns != 'CellType2') & (RIDGE_DF.columns != 'CellType3') & (RIDGE_DF.columns != 'Cell'))]
    factors = RIDGE_DF.columns
    X = RIDGE_DF.values
    X = StandardScaler().fit_transform(X)

    le = LabelEncoder()
    le.fit(cellTypeCol)
    y = le.transform(cellTypeCol)

    ridgeMod = RidgeClassifierCV(cv=5)
    ridgeMod.fit(X, y)
    TargCoefs = ridgeMod.coef_[np.where(le.classes_ == targCell), :].ravel()
    TargCoefsDF = pd.DataFrame({"Marker": factors, "Coefficient": TargCoefs}).sort_values(by="Coefficient")
    TargCoefsDF = TargCoefsDF.tail(numFactors)
    sns.barplot(data=TargCoefsDF, x="Marker", y="Coefficient", ax=ax, color='k')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)


def CITE_SVM(ax, targCell, numFactors=10, sampleFrac=0.5):
    """Fits a ridge classifier to the CITE data and plots those most highly correlated with T reg"""
    SVMmod = SVC()
    SVC_DF = importCITE()
    cellToI = SVC_DF.CellType2.unique()
    SVC_DF = SVC_DF.loc[(SVC_DF["CellType2"].isin(cellToI)), :]
    SVC_DF = SVC_DF.sample(frac=sampleFrac, random_state=1)
    cellTypeCol = SVC_DF.CellType2.values
    SVC_DF = SVC_DF.loc[:, ((SVC_DF.columns != 'CellType1') & (SVC_DF.columns != 'CellType2') & (SVC_DF.columns != 'CellType3') & (SVC_DF.columns != 'Cell'))]
    factors = SVC_DF.columns
    X = SVC_DF.values
    X = StandardScaler().fit_transform(X)
    CD122col = X[:, np.where(factors == "CD122")].reshape(-1, 1)

    enc = LabelBinarizer()
    y = enc.fit_transform(cellTypeCol)
    TregY = y[:, np.where(enc.classes_ == targCell)].ravel()

    AccDF = pd.DataFrame(columns=["Markers", "Accuracy"])
    baselineAcc = SVMmod.fit(CD122col, TregY).score(CD122col, TregY)
    print(baselineAcc)
    print(np.where((factors == "CD122")))
    for marker in factors:
        SVMmod = SVC()
        print(marker)
        markerCol = X[:, np.where(factors == marker)]
        CD122MarkX = np.hstack((CD122col, markerCol.reshape(-1, 1)))
        markAcc = SVMmod.fit(CD122MarkX, TregY).score(CD122MarkX, TregY)
        print(markAcc)
        AccDF = AccDF.append(pd.DataFrame({"Markers": [marker], "Accuracy": [markAcc]}))

    AccDF = AccDF.sort_values(by="Accuracy")
    markers = copy(AccDF.tail(numFactors).Markers.values)  # Here
    AccDF.Markers = "CD122 + " + AccDF.Markers

    plot_DF = AccDF.tail(numFactors).append(pd.DataFrame({"Markers": ["CD122 only"], "Accuracy": [baselineAcc]}))
    sns.barplot(data=plot_DF, x="Markers", y="Accuracy", ax=ax, color='k')
    ax.set(ylim=(0.95, 1))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
