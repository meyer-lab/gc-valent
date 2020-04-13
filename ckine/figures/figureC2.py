"""
This creates Figure 1 for Single Cell FC data analysis. Examples of PCA loadings/scores plots and comparisons to gating.
"""

import os
import matplotlib.lines as mlines
import pandas as pds
import numpy as np
from scipy import stats
from .figureCommon import subplotLabel, getSetup
from ..flow import importF
from ..PCA import StatGini, sampleT, sampleNK
from ..flow import gating, count_data

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((10, 10), (4, 4))

    subplotLabel(ax)

    StatMV()

    # global_legend(ax[7])

    return f


def global_legend(ax):
    """ Create legend for Inverse and Standard Gini """
    blue = mlines.Line2D([], [], color='navy', marker='o', linestyle='None', markersize=6, label='Gini Coeff')
    orange = mlines.Line2D([], [], color='darkorange', marker='o', linestyle='None', markersize=6, label='Inverse Gini Coeff')
    ax.legend(handles=[orange, blue], bbox_to_anchor=(0, 1), loc="upper left")


def StatMV():
    """
    Calculate mean and variance of a sample in a pandas dataframe, and plot.
    """
    
    dataFiles = ["/data/flow/2019-03-19 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate.zip", "/data/flow/2019-03-27 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate.zip", "/data/flow/2019-04-18 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate - NEW PBMC LOT/", "/data/flow/2019-03-15 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate.zip", "/data/flow/2019-03-27 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate.zip", "/data/flow/2019-04-18 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate - NEW PBMC LOT.zip"]
    dataFiles = ["/home/brianoj/Tplate15", "/home/brianoj/Tplate27", "/home/brianoj/Tplate418", "/home/brianoj/Nkplate15", "/home/brianoj/Nkplate27", "/home/brianoj/Nkplate418"]
    dates = ["2019-3-19", "2019-3-27", "2019-4-18", "2019-3-15", "2019-3-27", "2019-4-18"]
    rows = ['A', 'B', 'C', 'D']
    cellTypesT = [False, 'treg', 'tregMem', 'tregNaive', 'nonTreg', 'THelpMem', 'THelpN']
    cellTypesNK = [False, "nk", "nkt", "bnk", "cd"]
    TitlesT = ["Tcells", "T-regs", "Mem Treg", "Naive Treg", "T-helper", "Mem Th", "Naive Th"]
    TitlesNK = ['NKcells' "NK", "NKT", "BNK", "CD8+"]
    masterMVdf = pds.DataFrame(columns={"Date", "Time", "Cell", "Dose", "Mean", "Variance", "Skew", "Kurtosis"})
    MVdf = pds.DataFrame(columns={"Date", "Time", "Dose", "Mean", "Variance", "Skew", "Kurtosis"})
    alldata = []
    dosemat = np.array([[84, 28, 9.333333, 3.111, 1.037037, 0.345679, 0.115226, 0.038409, 0.012803, 0.004268, 0.001423, 0.000474]])

    for i, filename in enumerate(dataFiles):
        if i < 3:
            Tcells = True
        else:
            Tcells = False
        for j, row in enumerate(rows):
            sample, _ = importF(filename, row)
            if Tcells:
                statcol = "RL1-H"
                for k, cell_type in enumerate(cellTypesT):
                    if cell_type:
                        gates = gating(cell_type)
                        _, alldata = count_data(sample, gates, Tcells)
                    else:
                        for jj, samplejj in enumerate(sample):
                            _, pstat, _ = sampleT(samplejj)
                            alldata.append(pstat)

                    for ii, sampleii in enumerate(sample):  # get pstat data and put it into list form
                        dat_array = alldata[ii]
                        stat_array = dat_array[[statcol]]
                        stat_array = stat_array.to_numpy()
                        stat_array = stat_array.clip(min=1)  # remove small percentage of negative pstat values
                        if stat_array.size == 0:
                            MVdf = pds.DataFrame.from_dict({"Date": dates[i], "Time": timeFunc(row), "Cell":TitlesT[k], "Dose": dosemat[0, ii], "Mean": [0], "Variance": [0], "Skew": [0], "Kurtosis": [0]})
                        else:
                            MVdf = MVdf.append(pds.DataFrame.from_dict({"Date": dates[i], "Time": timeFunc(row), "Cell":TitlesT[k], "Dose": dosemat[0, ii], "Mean": np.mean(stat_array), "Variance": np.var(
                            stat_array), "Skew": stats.skew(stat_array), "Kurtosis": stats.kurtosis(stat_array)}))

                    MVdf['Mean'] = MVdf['Mean'] - MVdf['Mean'].min()
                    masterMVdf = masterMVdf.append(MVdf)
                    MVdf = pds.DataFrame(columns={"Date", "Time", "Dose", "Mean", "Variance", "Skew", "Kurtosis"})
            else:
                statcol = "BL2-H"
                for k, cell_type in enumerate(cellTypesNK):
                    if cell_type:
                        gates = gating(cell_type)
                        _, alldata = count_data(sample, gates, Tcells) 
                    else:
                        for jj, samplejj in enumerate(sample):
                            _, pstat, _ = sampleNK(samplejj)
                            alldata.append(pstat)

                    for ii, sampleii in enumerate(sample):  # get pstat data and put it into list form
                        dat_array = alldata[ii]
                        stat_array = dat_array[[statcol]]
                        stat_array = stat_array.to_numpy()
                        stat_array = stat_array.clip(min=1)  # remove small percentage of negative pstat values
                        if stat_array.size == 0:
                            MVdf = MVdf.append(pds.DataFrame.from_dict({"Date": dates[i], "Time": timeFunc(row), "Cell":TitlesT[k], "Dose": dosemat[0, ii], "Mean": [0], "Variance": [0], "Skew": [0], "Kurtosis": [0]}))
                        else:
                            MVdf = MVdf.append(pds.DataFrame.from_dict({"Date": dates[i], "Time": timeFunc(row), "Cell":TitlesT[k], "Dose": dosemat[0, ii], "Mean": np.mean(stat_array), "Variance": np.var(
                            stat_array), "Skew": stats.skew(stat_array), "Kurtosis": stats.kurtosis(stat_array)}))

                    MVdf['Mean'] = MVdf['Mean'] - MVdf['Mean'].min()
                    masterMVdf = masterMVdf.append(MVdf)
                    MVdf = pds.DataFrame(columns={"Date", "Time", "Dose", "Mean", "Variance", "Skew", "Kurtosis"})
    masterMVdf.to_csv("VarianceData.csv", index = False)
    #masterMVdf.to_csv("/home/brianoj/VarData/")

    return MVdf


def timeFunc(letter):
    if letter == "A":
        return 4.0
    elif letter == "B":
        return 2.0
    elif letter == "C":
        return 1.0
    elif letter == "D":
        return 0.5