"""
This file contains functions that are used in multiple figures.
"""
import sys
import logging
import time
from string import ascii_lowercase
import matplotlib
import seaborn as sns
import numpy as np
import svgutils.transform as st
import pandas as pd
from sklearn.preprocessing import StandardScaler, LabelEncoder, LabelBinarizer
from sklearn.linear_model import RidgeClassifierCV
from sklearn.metrics import balanced_accuracy_score
from scipy import stats
from sklearn.neighbors import KernelDensity
from sklearn.svm import SVC
from scipy.optimize import least_squares
from copy import copy
from matplotlib import gridspec, pyplot as plt
from ..imports import import_pstat_all, importCITE, importRNACITE

matplotlib.use('AGG')
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


matplotlib.rcParams["legend.labelspacing"] = 0.2
matplotlib.rcParams["legend.fontsize"] = 8
matplotlib.rcParams["xtick.major.pad"] = 2
matplotlib.rcParams["ytick.major.pad"] = 2
matplotlib.rcParams["xtick.minor.pad"] = 1.9
matplotlib.rcParams["ytick.minor.pad"] = 1.9
matplotlib.rcParams["legend.handletextpad"] = 0.5
matplotlib.rcParams["legend.handlelength"] = 0.5
matplotlib.rcParams["legend.framealpha"] = 0.5
matplotlib.rcParams["legend.markerscale"] = 0.7
matplotlib.rcParams["legend.borderpad"] = 0.35


dosemat = np.array([84, 28, 9.333333, 3.111, 1.037037, 0.345679, 0.115226, 0.038409, 0.012803, 0.004268, 0.001423, 0.000474])


def getSetup(figsize, gridd, multz=None, empts=None, constrained=True):
    """ Establish figure set-up with subplots. """
    sns.set(style="whitegrid", font_scale=0.7, color_codes=True, palette="colorblind", rc={"grid.linestyle": "dotted", "axes.linewidth": 0.6})

    # create empty list if empts isn't specified
    if empts is None:
        empts = []

    if multz is None:
        multz = dict()

    # Setup plotting space and grid
    f = plt.figure(figsize=figsize, constrained_layout=constrained)
    gs1 = gridspec.GridSpec(*gridd, figure=f)

    # Get list of axis objects
    x = 0
    ax = list()
    while x < gridd[0] * gridd[1]:
        if x not in empts and x not in multz.keys():  # If this is just a normal subplot
            ax.append(f.add_subplot(gs1[x]))
        elif x in multz.keys():  # If this is a subplot that spans grid elements
            ax.append(f.add_subplot(gs1[x: x + multz[x] + 1]))
            x += multz[x]
        x += 1

    return (ax, f)


def subplotLabel(axs):
    """ Place subplot labels on figure. """
    for ii, ax in enumerate(axs):
        ax.text(-0.2, 1.25, ascii_lowercase[ii], transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")


def overlayCartoon(figFile, cartoonFile, x, y, scalee=1, scale_x=1, scale_y=1):
    """ Add cartoon to a figure file. """

    # Overlay Figure cartoons
    template = st.fromfile(figFile)
    cartoon = st.fromfile(cartoonFile).getroot()

    cartoon.moveto(x, y)
    cartoon.scale(scalee)

    template.append(cartoon)
    template.save(figFile)


def genFigure():
    """ Main figure generation function. """
    fdir = './output/'
    start = time.time()
    nameOut = 'figure' + sys.argv[1]

    exec('from ckine.figures import ' + nameOut)
    ff = eval(nameOut + '.makeFigure()')
    ff.savefig(fdir + nameOut + '.svg', dpi=ff.dpi, bbox_inches='tight', pad_inches=0)

    # if sys.argv[1] == 'C1':
    # Overlay Figure 1 cartoon
    # overlayCartoon(fdir + 'figureC1.svg',
    #               './ckine/graphics/muteinsCartoon.svg', 1200, 350, scalee=0.040)

    if sys.argv[1] == 'C2':
        # Overlay Figure 2 cartoon
        overlayCartoon(fdir + 'figureC2.svg',
                       './ckine/graphics/selectivityCartoon.svg', 1200, 350, scalee=0.02)

        # overlayCartoon(fdir + 'figureC2.svg',
        #               './ckine/graphics/citeCartoon.svg', 2300, 20500, scalee=0.043)

    if sys.argv[1] == 'C3':
        # Overlay Figure 3 cartoon
        overlayCartoon(fdir + 'figureC3.svg',
                       './ckine/graphics/tensor4D.svg', 20, 6, scalee=1.62)

    if sys.argv[1] == 'C4':
        # Overlay Figure 4 cartoon
        overlayCartoon(fdir + 'figureC4.svg',
                       './ckine/graphics/ModelCartoon.svg', 1450, 0, scalee=0.023)

    logging.info('%s is done after %s seconds.', nameOut, time.time() - start)


cellSTATlimDict = {"Treg": (0, 60000),
                   "Thelper": (0, 40000),
                   "CD8": (0, 30000),
                   "NK": (0, 8000),
                   "NKBright": (0, 8000)}


def plotDoseResponses(ax, df, mut, cellType, val=False):
    """Plots all experimental vs. Predicted Values"""
    #df = df.groupby(["Cell", "Valency", "Ligand", "Dose", "Time"])["Experimental"].mean().reset_index()
    palette = getLigDict()[mut]
    df = df.rename(cellTypeDict)
    if isinstance(cellType, str):
        if val:
            expData = df.loc[(df.Ligand == mut) & (df.Valency == val) & (df.Cell == cellType)]
        else:
            expData = df.loc[(df.Ligand == mut) & (df.Cell == cellType)]
        valList = expData.Valency.unique()
        predData = expData.groupby(["Cell", "Valency", "Ligand", "Dose"])["Predicted"].mean().reset_index()

        if val:
            if val == 1:
                sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", size="Valency", ax=ax, color=palette)
                sns.lineplot(x="Dose", y="Predicted", data=predData, label="Predicted", size="Valency", ax=ax, color=palette)
                ax.set(title=cellTypeDict[cellType], xlabel=r"$log_{10}$ Monomeric " + mut + " (nM)", ylabel="pSTAT5", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
            if val == 2:
                sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", size="Valency", ax=ax, color=palette)
                sns.lineplot(x="Dose", y="Predicted", data=predData, label="Predicted", size="Valency", ax=ax, color=palette)
                ax.set(title=cellTypeDict[cellType], xlabel=r"$log_{10}$ Dimeric " + mut + " (nM)", ylabel="pSTAT5", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
        else:
            if len(valList) > 1:
                sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", size="Valency", ax=ax, legend="brief", color=palette)
                sns.lineplot(x="Dose", y="Predicted", data=predData, label="Predicted", size="Valency", ax=ax, legend="brief", color=palette)
                ax.set(title=cellTypeDict[cellType], xlabel=r"$log_{10}$ " + mut + " (nM)", ylabel="pSTAT5", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
                handles, labels = ax.get_legend_handles_labels()
                ax.legend([handles[0]] + handles[4::], [labels[0]] + labels[4::])
            else:
                if valList[0] == 1:
                    sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", size="Valency", ax=ax, color=palette)
                    sns.lineplot(x="Dose", y="Predicted", data=predData, label="Predicted", size="Valency", ax=ax, color=palette)
                    ax.set(title=cellTypeDict[cellType], xlabel=r"$log_{10}$ Monomeric " + mut + " (nM)", ylabel="pSTAT5", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[0:2] + handles[4::], labels[0:2] + labels[4::])
                if valList[0] == 2:
                    sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", size="Valency", ax=ax, color=palette)
                    sns.lineplot(x="Dose", y="Predicted", data=predData, label="Predicted", size="Valency", ax=ax, color=palette)
                    ax.set(title=cellTypeDict[cellType], xlabel=r"$log_{10}$ Dimeric " + mut + " (nM)", ylabel="pSTAT5", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[0:2] + handles[4::], labels[0:2] + labels[4::])
    else:
        expData = df.loc[(df.Ligand == mut) & (df.Valency == val) & (df.Cell.isin(cellType))]
        predData = expData.groupby(["Cell", "Valency", "Ligand", "Dose"])["Predicted"].mean().reset_index()
        sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", hue="Cell", ax=ax)
        sns.lineplot(x="Dose", y="Predicted", data=predData, label="Predicted", hue="Cell", ax=ax)
        if val == 1:
            ax.set(title=cellTypeDict[cellType[0]] + "s", xlabel=r"$log_{10}$ Monomeric " + mut + " (nM)", ylabel="pSTAT5", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType[0]])
        if val == 2:
            ax.set(title=cellTypeDict[cellType[0]] + "s", xlabel=r"$log_{10}$ Dimeric " + mut + " (nM)", ylabel="pSTAT5", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType[0]])


def getLigDict():
    """Gives hue dict for ligands - consistency across """
    pSTATDF = import_pstat_all(True, False)
    ligands = pSTATDF.Ligand.unique()
    #palette = sns.color_palette("Spectral", ligands.size)
    palette = sns.color_palette(["royalblue", "#1f78b4", "limegreen", "#33a02c", "#fb9a99", "#e31a1c", "orange", "#ff7f00", "deeppink", "#6a3d9a", "orange"])

    ligDict = {}
    for i, ligand in enumerate(ligands):
        ligDict[ligand] = palette[i]
    ligDict["Live/Dead"] = palette[-1]
    return ligDict


cellTypeDict = {"Treg": r"T$_{reg}$",
                "Treg $IL2Ra^{hi}$": r"T$_{reg}$ $IL2Ra^{hi}$",
                "Treg $IL2Ra^{lo}$": r"T$_{reg}$ $IL2Ra^{lo}$",
                "Thelper $IL2Ra^{hi}$": r"T$_{helper}$ $IL2Ra^{hi}$",
                "Thelper $IL2Ra^{lo}$": r"T$_{helper}$ $IL2Ra^{lo}$",
                "Thelper": r"T$_{helper}$",
                "NK": "NK",
                "CD8": r"CD8$^{+}$",
                "NKBright": r"NK$^{Bright}$"}

doseLimDict = {r"T$_{reg}$": (0, 50000),
               r"T$_{reg}$ $IL2Ra^{hi}$": (0, 50000),
               r"T$_{reg}$ $IL2Ra^{lo}$": (0, 50000),
               r"T$_{helper}$": (0, 25000),
               r"T$_{helper}$ $IL2Ra^{hi}$": (0, 25000),
               r"T$_{helper}$ $IL2Ra^{lo}$": (0, 25000),
               "NK": (0, 5000),
               r"NK$^{Bright}$": (0, 5000),
               r"CD8$^{+}$": (0, 8000)}


def get_cellTypeDict():
    """Returns dict for beautifying cell type names"""
    return cellTypeDict


def get_doseLimDict():
    """Returns dict for dose response limits"""
    return doseLimDict


def getLigandLegend():
    """Creates dummy plot and returns handles for ligands"""
    f, ax = plt.subplots()
    ligDict = getLigDict()
    respDF = import_pstat_all(True, False)
    respDF = respDF.groupby(["Ligand", "Dose"]).Mean.mean().reset_index()
    sns.scatterplot(data=respDF, x="Dose", y="Mean", hue="Ligand", legend=True, palette=ligDict, ax=ax)
    return ax.get_legend()


def Wass_KL_Dist(ax, targCell, numFactors, RNA=False):
    """Finds markers which have average greatest difference from other cells"""
    if RNA:
        CITE_DF = importRNACITE()
    else:
        CITE_DF = importCITE()

    markerDF = pd.DataFrame(columns=["Marker", "Cell Type", "Amount"])
    for marker in CITE_DF.loc[:, ((CITE_DF.columns != 'CellType1') & (CITE_DF.columns != 'CellType2') & (CITE_DF.columns != 'CellType3') & (CITE_DF.columns != 'Cell'))].columns:
        marker_values = CITE_DF[marker].to_numpy()
        markAvg = np.mean(marker_values)
        cell_type_idxs = (CITE_DF["CellType2"] == targCell).to_numpy()

        if markAvg > 0.0001:
            targCellMark = marker_values[cell_type_idxs] / markAvg
            offTargCellMark = marker_values[np.logical_not(cell_type_idxs)] / markAvg
            if np.mean(targCellMark) > np.mean(offTargCellMark):
                kdeTarg = KernelDensity(kernel='gaussian').fit(targCellMark.reshape(-1, 1))
                kdeOffTarg = KernelDensity(kernel='gaussian').fit(offTargCellMark.reshape(-1, 1))
                minVal = np.minimum(targCellMark.min(), offTargCellMark.min()) - 10
                maxVal = np.maximum(targCellMark.max(), offTargCellMark.max()) + 10
                outcomes = np.arange(minVal, maxVal + 1).reshape(-1, 1)
                distTarg = np.exp(kdeTarg.score_samples(outcomes))
                distOffTarg = np.exp(kdeOffTarg.score_samples(outcomes))
                KL_div = stats.entropy(distOffTarg.flatten() + 1e-200, distTarg.flatten() + 1e-200, base=2)
                markerDF = pd.concat([markerDF, pd.DataFrame({"Marker": [marker], "Wasserstein Distance": stats.wasserstein_distance(targCellMark, offTargCellMark), "KL Divergence": KL_div})])

    corrsDF = pd.DataFrame()
    for i, distance in enumerate(["Wasserstein Distance", "KL Divergence"]):
        ratioDF = markerDF.sort_values(by=distance)
        posCorrs = ratioDF.tail(numFactors).Marker.values
        corrsDF = pd.concat([corrsDF, pd.DataFrame({"Distance": distance, "Marker": posCorrs})])
        markerDF = markerDF.loc[markerDF["Marker"].isin(posCorrs)]
        sns.barplot(data=ratioDF.tail(numFactors), y="Marker", x=distance, ax=ax[i], color='k')
        ax[i].set(xscale="log")
        #ax[i].set_xticklabels(ax[i].get_xticklabels(), rotation=45)
    if RNA:
        ax[0].set(title="Wasserstein Distance - RNA")
        ax[1].set(title="KL Divergence - RNA")
    else:
        ax[0].set(title="Wasserstein Distance - Surface Markers")
        ax[1].set(title="KL Divergence - Surface Markers")
    return corrsDF


def CITE_RIDGE(ax, targCell, numFactors=10, RNA=False):
    """Fits a ridge classifier to the CITE data and plots those most highly correlated with T reg"""
    ridgeMod = RidgeClassifierCV()
    if RNA:
        RIDGE_DF = importRNACITE()
    else:
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
    if RNA:
        ax.set(title="RIDGE Coefficients - RNA")
    else:
        ax.set(title="RIDGE Coefficients - Surface Markers")
    return TargCoefsDF


def CITE_SVM(ax, targCell, numFactors=10, sampleFrac=0.5, RNA=False):
    """Fits a ridge classifier to the CITE data and plots those most highly correlated with T reg"""
    SVMmod = SVC()
    if RNA:
        SVC_DF = importRNACITE()
        IL2RB = "IL2RB"
    else:
        SVC_DF = importCITE()
        IL2RB = "CD122"
    cellToI = SVC_DF.CellType2.unique()
    SVC_DF = SVC_DF.loc[(SVC_DF["CellType2"].isin(cellToI)), :]
    SVC_DF = SVC_DF.sample(frac=sampleFrac, random_state=1)
    cellTypeCol = SVC_DF.CellType2.values
    SVC_DF = SVC_DF.loc[:, ((SVC_DF.columns != 'CellType1') & (SVC_DF.columns != 'CellType2') & (SVC_DF.columns != 'CellType3') & (SVC_DF.columns != 'Cell'))]
    factors = SVC_DF.columns
    X = SVC_DF.values
    X = StandardScaler().fit_transform(X)
    CD122col = X[:, np.where(factors == IL2RB)].reshape(-1, 1)

    enc = LabelBinarizer()
    y = enc.fit_transform(cellTypeCol)
    TregY = y[:, np.where(enc.classes_ == targCell)].ravel()

    AccDF = pd.DataFrame(columns=["Markers", "Accuracy"])
    baselineAccMod = SVMmod.fit(CD122col, TregY)
    baselineAcc = balanced_accuracy_score(TregY, baselineAccMod.predict(CD122col))
    for marker in factors:
        SVMmod = SVC()
        markerCol = X[:, np.where(factors == marker)]
        CD122MarkX = np.hstack((CD122col, markerCol.reshape(-1, 1)))
        fitMod = SVMmod.fit(CD122MarkX, TregY)
        markAcc = balanced_accuracy_score(TregY, fitMod.predict(CD122MarkX))
        AccDF = pd.concat([AccDF, pd.DataFrame({"Markers": [marker], "Accuracy": [markAcc]})])

    AccDF = AccDF.sort_values(by="Accuracy")
    markers = copy(AccDF.tail(numFactors).Markers.values)  # Here
    AccDF.Markers = "IL2RB" + " + " + AccDF.Markers

    plot_DF = pd.concat([AccDF.tail(numFactors), pd.DataFrame({"Markers": ["IL2RB only"], "Accuracy": [baselineAcc]})])
    sns.barplot(data=plot_DF, x="Markers", y="Accuracy", ax=ax, color='k')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=75)
    ax.set(ylabel="Balanced Accuracy")
    if RNA:
        ax.set(title="SVM Accuracy - RNA", ylim=(0.4, 0.6))
    else:
        ax.set(title="SVM Accuracy - Surface Markers", ylim=(0.4, 0.8))
    return markers


palette_dict = {"R38Q/H16N": "darkorchid",
                "Live/Dead": "Orange"}
valency_dict = {1: "royalblue", 2: "goldenrod", 4: "limegreen"}


def get_valency_dict():
    return valency_dict


def ligandPlot(DF, cell, ax, live_dead=False):
    """Plots a cell type response"""
    if live_dead:
        DF = DF.loc[(DF.Cell == cell) & (DF.Ligand != "R38Q/H16N")]
    else:
        DF = DF.loc[(DF.Cell == cell) & (DF.Ligand != "Live/Dead")]

    plotDF = DF[["Dose", "Valency", "Ligand", "Date", "Experimental", "Predicted"]].copy()
    for date in DF.Date.unique():
        maxobs = DF.loc[(DF.Date == date), "Experimental"].max()
        plotDF.loc[plotDF.Date == date, "Experimental"] = DF.loc[(DF.Date == date)]["Experimental"].values / maxobs
        plotDF.loc[plotDF.Date == date, "Predicted"] = DF.loc[(DF.Date == date)]["Predicted"].values / maxobs
    plotDF = plotDF.replace(cellTypeDict)
    sns.lineplot(data=plotDF, x="Dose", y="Predicted", hue="Valency", palette=valency_dict, ax=ax)
    sns.lineplot(
        data=plotDF,
        x="Dose",
        y="Experimental",
        hue="Valency",
        palette=valency_dict,
        linewidth=0,
        markers=True,
        err_style="bars",
        err_kws={
            "capsize": 3,
            "elinewidth": 1},
        legend=False,
        ax=ax)
    plotDF = plotDF.groupby(["Ligand", "Valency", "Dose"]).Experimental.mean().reset_index()
    sns.scatterplot(data=plotDF, x="Dose", y="Experimental", hue="Valency", s=12, palette=valency_dict, ax=ax)
    if live_dead:
        ax.set(xscale="log", xlabel="Dose Live/Dead (nM)", ylabel="pSTAT5", title=cellTypeDict[cell], xticks=[0.0001, 0.01, 1, 100], yticks=[0, 0.5, 1], ylim=[0, 1.25])
    else:
        ax.set(xscale="log", xlabel="Dose R38Q/H16N (nM)", ylabel="pSTAT5", title=cellTypeDict[cell], xticks=[0.0001, 0.01, 1, 100], yticks=[0, 0.5, 1], ylim=[0, 1.25])


def ligand_ratio_plot(DF, cell1, cell2, ax, live_dead=False):
    """Plots a cell type response"""
    expRatioDF = pd.DataFrame()
    predRatioDF = pd.DataFrame()
    if live_dead:
        DF = DF.loc[(DF.Valency != 1) & (DF.Ligand != "R38Q/H16N")]
    else:
        DF = DF.loc[(DF.Ligand != "Live/Dead")]

    for date in DF.Date.unique():
        max1 = DF.loc[(DF.Date == date) & (DF.Cell == cell1), "Experimental"].max()
        max2 = DF.loc[(DF.Date == date) & (DF.Cell == cell2), "Experimental"].max()
        DF.loc[(DF.Date == date) & (DF.Cell == cell1), "Experimental"] /= max1
        DF.loc[(DF.Date == date) & (DF.Cell == cell1), "Predicted"] /= max1
        DF.loc[(DF.Date == date) & (DF.Cell == cell2), "Experimental"] /= max2
        DF.loc[(DF.Date == date) & (DF.Cell == cell2), "Predicted"] /= max2

    for dose in DF.Dose.unique():
        for ligand in DF.Ligand.unique():
            for valency in DF.loc[DF.Ligand == ligand].Valency.unique():
                for date in DF.Date.unique():
                    expRatio = DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell1) & (DF.Date == date)].Experimental.values / \
                        (DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell2) & (DF.Date == date)].Experimental.values + 0.1)
                    predRatio = DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell1) & (DF.Date == date)].Predicted.values / \
                        (DF.loc[(DF.Ligand == ligand) & (DF.Dose == dose) & (DF.Valency == valency) & (DF.Cell == cell2) & (DF.Date == date)].Predicted.values + 0.1)
                    expRatioDF = pd.concat([expRatioDF, pd.DataFrame({"Dose": dose, "Ligand": ligand, "Valency": valency, "Ratio": expRatio, "Date": date})])
                    predRatioDF = pd.concat([predRatioDF, pd.DataFrame({"Dose": dose, "Ligand": ligand, "Valency": valency, "Ratio": predRatio, "Date": date})])

    expRatioDF = expRatioDF.reset_index()
    predRatioDF = predRatioDF.reset_index()
    sns.lineplot(data=expRatioDF, x="Dose", y="Ratio", hue="Valency", style="Ligand", palette=valency_dict, ax=ax, markers=True, err_style="bars", err_kws={"capsize": 2, "elinewidth": 1},  legend=False)
    expRatioDF = expRatioDF.groupby(["Ligand", "Valency", "Dose"]).Ratio.mean().reset_index()
    sns.scatterplot(data=expRatioDF, x="Dose", y="Ratio", hue="Valency", style="Ligand", s=15, palette=valency_dict, ax=ax, legend=False)
    if live_dead:
        ax.set(xscale="log", xlabel="Dose Live/Dead (nM)", ylabel="Ratio", title=cellTypeDict[cell1] + " to " + cellTypeDict[cell2] + " Ratio")
    else:
        ax.set(xscale="log", xlabel="R38Q/H16N (nM)", ylabel="Ratio", title=cellTypeDict[cell1] + " to " + cellTypeDict[cell2] + " Ratio")
    ax.set(xticks=[0.0001, 0.01, 1, 100], yticks=[0, 2, 4, 6, 8])


ligDict = getLigDict()


def hillRatioDosePlot(ax, respDF, time, targCell, offTargCell, pseudo=0.2, plot=True):
    """Plots the various affinities for IL-2 Muteins"""
    doses = np.log10(np.logspace(np.log10(respDF.Dose.min()), np.log10(respDF.Dose.max()), 100)) + 4
    x0 = [4, 1, 2]
    hillDF = pd.DataFrame()
    Ligands = respDF.Ligand.unique()
    respDF = respDF.loc[(respDF.Time == time)]

    for ligand in respDF.Ligand.unique():
        for valency in respDF.loc[respDF.Ligand == ligand].Valency.unique():
            targIsoData = respDF.loc[(respDF.Ligand == ligand) & (respDF.Valency == valency) & (respDF.Cell == targCell)]
            targXData = np.nan_to_num(np.log10(targIsoData.Dose.values)) + 4
            targYData = np.nan_to_num(targIsoData.Mean.values)
            targFit = least_squares(hill_residuals, x0, method="dogbox", args=(targXData, targYData), bounds=([0.0, 0.0, 2], [5, 10.0, 6]), jac="3-point")

            offTIsoData = respDF.loc[(respDF.Ligand == ligand) & (respDF.Valency == valency) & (respDF.Cell == offTargCell)]
            offTargXData = np.nan_to_num(np.log10(offTIsoData.Dose.values)) + 4
            offTargYData = np.nan_to_num(offTIsoData.Mean.values)
            offTargFit = least_squares(hill_residuals, x0, method="dogbox", args=(offTargXData, offTargYData), bounds=([0.0, 0.0, 2], [5, 10.0, 6]), jac="3-point")
            hillDF = pd.concat([hillDF, pd.DataFrame({"Ligand": ligand, "Valency": valency, "Cell": targCell, "Dose": np.power(
                10, doses - 4), targCell: hill_equation(targFit.x, doses), offTargCell: hill_equation(offTargFit.x, doses)})])

    for cell in [targCell, offTargCell]:
        maxobs = hillDF.loc[(hillDF.Ligand == "IL2")][cell].max()
        hillDF[cell] /= maxobs
    hillDF["Ratio"] = hillDF[targCell] / (pseudo + hillDF[offTargCell])

    hillDF = hillDF.groupby(["Ligand", "Valency", "Dose"]).Ratio.mean().reset_index()
    hillDF = hillDF.loc[(hillDF.Ligand != "IL15") & (hillDF.Ligand != "IL2")]
    if plot:
        sns.lineplot(data=hillDF, x="Dose", y="Ratio", hue="Ligand", size="Valency", ax=ax, palette=ligDict, sizes=(1, 2.5))
        ax.set(xscale="log", xlim=(1e-4, 1e2), ylim=(0, 5))
    return hillDF


def hill_equation(x, dose):
    """ Calculates EC50 from Hill Equation. """
    #print(x, dose)
    EMax = np.power(10, x[0])
    n = x[1]
    EC50 = x[2]
    return EMax * np.power(dose, n) / (np.power(EC50, n) + np.power(dose, n))


def hill_residuals(x, dose, y):
    """ Residual function for Hill Equation. """
    return hill_equation(x, dose) - y


def make_EC50_DF(respDF, time, meyer=False):
    """Plots the various affinities for IL-2 Muteins"""
    doses = np.log10(np.logspace(np.log10(respDF.Dose.min()), np.log10(respDF.Dose.max()), 100)) + 4
    x0 = [4, 1, 4]
    EC50_DF = pd.DataFrame()
    if not meyer:
        respDF = respDF.loc[(respDF.Time == time)]
    else:
        respDF = respDF.rename({"pSTAT5": "Mean"}, axis=1)
    cells = respDF.Cell.unique()

    for cell in cells:
        for ligand in respDF.Ligand.unique():
            for valency in respDF.loc[respDF.Ligand == ligand].Valency.unique():
                targIsoData = respDF.loc[(respDF.Ligand == ligand) & (respDF.Valency == valency) & (respDF.Cell == cell)]
                targXData = np.nan_to_num(np.log10(targIsoData.Dose.values)) + 4
                targYData = np.nan_to_num(targIsoData.Mean.values)
                targFit = least_squares(hill_residuals, x0, method="dogbox", args=(targXData, targYData), bounds=([0.0, 0.0, 2], [5, 10.0, 6]), jac="3-point")
                EC50_DF = pd.concat([EC50_DF, pd.DataFrame({"Cell": cell, "Ligand": ligand, "Valency": valency, "EC50": [np.power(10, targFit.x[2] - 4)]})])

    EC50_DF.loc[(EC50_DF.Valency == 1), "Ligand"] = (EC50_DF.loc[(EC50_DF.Valency == 1)].Ligand + " (Mono)").values
    EC50_DF.loc[(EC50_DF.Valency == 2), "Ligand"] = (EC50_DF.loc[(EC50_DF.Valency == 2)].Ligand + " (Biv)").values
    EC50_DF.loc[(EC50_DF.Valency == 4), "Ligand"] = (EC50_DF.loc[(EC50_DF.Valency == 4)].Ligand + " (Tetra)").values

    if meyer:
        EC50_DF = EC50_DF.pivot(index='Ligand', columns='Cell', values='EC50')
        EC50_DF = EC50_DF.where(EC50_DF < 99, "N/A")
        EC50_DF.to_csv("ckine/data/EC50_DF_Meyer.csv")
    else:
        EC50_DF = EC50_DF.pivot(index='Ligand', columns='Cell', values='EC50')
        EC50_DF = EC50_DF.where(EC50_DF < 99, "N/A")
        EC50_DF.to_csv("ckine/data/EC50_DF.csv")
