"""Stores Single Cell marker levels for monomeric signaling experiments"""
import os
import pandas as pds
import numpy as np
from scipy import stats

# import feather
from .figureCommon import subplotLabel, getSetup
from ..flow import importF
from ..FCimports import getChannelDict
from ..PCA import sampleT, sampleNK
from ..flow import gating, count_data

from ..FCimports import compMatrix, applyMatrix

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((10, 10), (4, 4))

    subplotLabel(ax)

    makeDimericGMM()

    return f


def makeDimericGMM():
    """
    Calculate mean and variance of a sample in a pandas dataframe, and plot.
    """

    dataFiles = [
        "/data/flow/2019-03-19 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate.zip",
        "/data/flow/2019-03-27 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate.zip",
        "/data/flow/2019-04-18 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - Treg plate - NEW PBMC LOT/",
        "/data/flow/2019-03-15 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate.zip",
        "/data/flow/2019-03-27 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate.zip",
        "/data/flow/2019-04-18 IL-2 and IL-15 treated pSTAT5 assay - Lymphocyte gated - NK plate - NEW PBMC LOT.zip",
    ]
    dataFiles = [
        "/home/brianoj/Tplate15",
        "/home/brianoj/Tplate27",
        "/home/brianoj/Tplate418",
        "/home/brianoj/Nkplate15",
        "/home/brianoj/Nkplate27",
        "/home/brianoj/Nkplate418",
    ]
    dates = [
        "3/15/2019",
        "3/27/2019",
        "4/18/2019",
        "3/15/2019",
        "3/27/2019",
        "4/18/2019",
    ]
    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    masterMVdf = pds.DataFrame(columns={"Date", "Time", "Ligand", "Dose"})
    MVdf = pds.DataFrame(columns={"Date", "Time", "Ligand", "Dose"})
    dosemat = np.array(
        [
            [
                84,
                28,
                9.333333,
                3.111,
                1.037037,
                0.345679,
                0.115226,
                0.038409,
                0.012803,
                0.004268,
                0.001423,
                0.000474,
            ]
        ]
    )
    cellTypesT = ["treg", "nonTreg"]
    cellTypesNK = ["nk", "cd"]
    TitlesT = ["Treg", "Thelper"]
    TitlesNK = ["NK", "CD8"]
    repList = [0, 0, 0, 0, 0, 0]

    T_matrix = compMatrix("2019-11-08", "1", "A")  # Create matrix 1
    Cd8_NKmatrix = compMatrix("2019-11-08", "1", "B")  # Create matrix 2

    for i, filename in enumerate(dataFiles):
        if i < 3:
            Tcells = True
        else:
            Tcells = False
        channelDict = getChannelDict(
            Tcell=Tcells, LiveCell=False, cellType=False, IL7=False
        )
        if Tcells:
            for j, row in enumerate(rows):
                print(filename)
                sample, _ = importF(filename, row)
                for jj, subSample in enumerate(sample):
                    sample[jj] = applyMatrix(subSample, T_matrix)
                cellTypeDict = {}

                for k, cell_type in enumerate(cellTypesT):
                    gates = gating(cell_type, dates[i], True, repList[i])
                    _, alldata, indices = count_data(sample, gates, Tcells, True)
                    cellTypeDict[TitlesT[k]] = indices

                for ii, sampleii in enumerate(
                    sample
                ):  # get pstat data and put it into list form
                    dat_array = sample[ii].get_data()
                    dat_array = dat_array[channelDict.keys()].rename(
                        columns=channelDict
                    )
                    dat_array["Date"] = dates[i]
                    dat_array["Time"] = timeFunc(row)
                    dat_array["Ligand"] = cytFunc(row)
                    dat_array["Dose"] = dosemat[0, ii]
                    dat_array["Valency"] = 1
                    dat_array["Cell Type"] = "None"
                    for cell in TitlesT:
                        dat_array.loc[cellTypeDict[cell][ii].values, "Cell Type"] = cell
                    MVdf = MVdf.append(dat_array)

                if j == 3 or j == 7:
                    MVdf["pSTAT5"] = (
                        MVdf["pSTAT5"] - MVdf.loc[(MVdf.Dose <= 0.001423)].pSTAT5.mean()
                    )
                    masterMVdf = masterMVdf.append(MVdf)
                    MVdf = pds.DataFrame()
        else:
            for j, row in enumerate(rows):
                print(filename)
                sample, _ = importF(filename, row)
                if (row == "H" and i == 4) is False:
                    for jj, subSample in enumerate(sample):
                        sample[jj] = applyMatrix(subSample, Cd8_NKmatrix)
                    cellTypeDict = {}

                    for k, cell_type in enumerate(cellTypesNK):
                        gates = gating(cell_type, dates[i], True, repList[i])
                        _, alldata, indices = count_data(sample, gates, Tcells, True)
                        cellTypeDict[TitlesNK[k]] = indices

                    for ii, sampleii in enumerate(
                        sample
                    ):  # get pstat data and put it into list form
                        dat_array = sample[ii].get_data()
                        dat_array = dat_array[channelDict.keys()].rename(
                            columns=channelDict
                        )
                        dat_array["Date"] = dates[i]
                        dat_array["Time"] = timeFunc(row)
                        dat_array["Ligand"] = cytFunc(row)
                        dat_array["Dose"] = dosemat[0, ii]
                        dat_array["Valency"] = 1
                        dat_array["Cell Type"] = "None"
                        for cell in TitlesNK:
                            dat_array.loc[
                                cellTypeDict[cell][ii].values, "Cell Type"
                            ] = cell
                        MVdf = MVdf.append(dat_array)
                if j == 3 or j == 7:
                    MVdf["pSTAT5"] = (
                        MVdf["pSTAT5"] - MVdf.loc[(MVdf.Dose <= 0.001423)].pSTAT5.mean()
                    )
                    masterMVdf = masterMVdf.append(MVdf)
                    MVdf = pds.DataFrame()

    dataFiles = [
        "/home/brianoj/Muteins 060-062 T/2019-04-19 IL2-060 IL2-062 Treg plate",
        "/home/brianoj/Muteins 088-097 T/2019-04-19 IL2-088 IL2-097 Treg plate",
        "/home/brianoj/Muteins 060-088 T/2019-05-02 IL2-060 IL2-088 Treg plate",
        "/home/brianoj/Muteins 062-097 T/2019-05-02 IL2-062 IL2-097 Treg plate",
        "/home/brianoj/Muteins 060-062 Nk/2019-04-19 IL2-060 IL2-062 NK plate",
        "/home/brianoj/Muteins 088-097 Nk/2019-04-19 IL2-088 IL2-097 NK plate",
        "/home/brianoj/Muteins 060-088 Nk/2019-05-02 IL2-060 IL2-088 NK plate",
        "/home/brianoj/Muteins 062-097 Nk/2019-05-02 IL2-062 IL2-097 NK plate",
    ]
    dates = [
        "4/19/2019",
        "4/19/2019",
        "5/2/2019",
        "5/2/2019",
        "4/19/2019",
        "4/19/2019",
        "5/2/2019",
        "5/2/2019",
    ]
    repList = [0, 1, 0, 1, 0, 1, 0, 1]

    print("Starting Muteins")

    for i, filename in enumerate(dataFiles):
        if i < 4:
            Tcells = True
        else:
            Tcells = False
        if Tcells:
            channelDict = getChannelDict(
                Tcell=Tcells, LiveCell=False, cellType=False, IL7=False
            )
            for j, row in enumerate(rows):
                print(filename)
                sample, _ = importF(filename, row)
                for jj, subSample in enumerate(sample):
                    sample[jj] = applyMatrix(subSample, T_matrix)
                cellTypeDict = {}

                for k, cell_type in enumerate(cellTypesT):
                    gates = gating(cell_type, dates[i], True, repList[i])
                    _, alldata, indices = count_data(sample, gates, Tcells, True)
                    cellTypeDict[TitlesT[k]] = indices

                for ii, sampleii in enumerate(
                    sample
                ):  # get pstat data and put it into list form
                    dat_array = sample[ii].get_data()
                    dat_array = dat_array[channelDict.keys()].rename(
                        columns=channelDict
                    )
                    timelig = mutFunc(row, filename)
                    dat_array["Date"] = dates[i]
                    dat_array["Time"] = timelig[0]
                    dat_array["Ligand"] = timelig[1]
                    dat_array["Dose"] = dosemat[0, ii]
                    dat_array["Valency"] = timelig[2] + 1
                    dat_array["Cell Type"] = "None"
                    for cell in TitlesT:
                        dat_array.loc[cellTypeDict[cell][ii].values, "Cell Type"] = cell
                    MVdf = MVdf.append(dat_array)

                if j == 3 or j == 7:
                    MVdf["pSTAT5"] = (
                        MVdf["pSTAT5"] - MVdf.loc[(MVdf.Dose <= 0.001423)].pSTAT5.mean()
                    )
                    masterMVdf = masterMVdf.append(MVdf)
                    MVdf = pds.DataFrame()
        else:
            for j, row in enumerate(rows):
                print(filename)
                sample, _ = importF(filename, row)
                if (row == "H" and i == 4) is False:
                    for jj, subSample in enumerate(sample):
                        sample[jj] = applyMatrix(subSample, Cd8_NKmatrix)
                    cellTypeDict = {}

                    for k, cell_type in enumerate(cellTypesNK):
                        gates = gating(cell_type, dates[i], True, repList[i])
                        _, alldata, indices = count_data(sample, gates, Tcells, True)
                        cellTypeDict[TitlesNK[k]] = indices

                    for ii, sampleii in enumerate(
                        sample
                    ):  # get pstat data and put it into list form
                        dat_array = sample[ii].get_data()
                        dat_array = dat_array[channelDict.keys()].rename(
                            columns=channelDict
                        )
                        timelig = mutFunc(row, filename)
                        dat_array["Date"] = dates[i]
                        dat_array["Time"] = timelig[0]
                        dat_array["Ligand"] = timelig[1]
                        dat_array["Dose"] = dosemat[0, ii]
                        dat_array["Valency"] = timelig[2] + 1
                        dat_array["Cell Type"] = "None"
                        for cell in TitlesNK:
                            dat_array.loc[
                                cellTypeDict[cell][ii].values, "Cell Type"
                            ] = cell
                        MVdf = MVdf.append(dat_array)

                if j == 3 or j == 7:
                    MVdf["pSTAT5"] = (
                        MVdf["pSTAT5"] - MVdf.loc[(MVdf.Dose <= 0.001423)].pSTAT5.mean()
                    )
                    masterMVdf = masterMVdf.append(MVdf)
                    MVdf = pds.DataFrame()

    masterMVdf.pSTAT5 = masterMVdf.pSTAT5.clip(lower=0)
    masterMVdf = masterMVdf.loc[masterMVdf.Ligand != "H16L N-term"]
    masterMVdf.reset_index().to_feather("FlowDataGMM_DimWT.ftr")

    return MVdf


def timeFunc(letter):
    if letter == "A" or letter == "E":
        return 4.0
    elif letter == "B" or letter == "F":
        return 2.0
    elif letter == "C" or letter == "G":
        return 1.0
    elif letter == "D" or letter == "H":
        return 0.5


def cytFunc(letter):
    if letter == "A" or letter == "B" or letter == "C" or letter == "D":
        return "IL2"
    elif letter == "E" or letter == "F" or letter == "G" or letter == "H":
        return "IL15"


def mutFunc(letter, datafile):
    if (
        datafile
        == "/home/brianoj/Muteins 060-062 T/2019-04-19 IL2-060 IL2-062 Treg plate"
        or datafile
        == "/home/brianoj/Muteins 060-062 Nk/2019-04-19 IL2-060 IL2-062 NK plate"
    ):
        if letter == "A":
            return [4.0, "WT N-term", 1]
        elif letter == "B":
            return [4.0, "WT N-term", 1]
        elif letter == "C":
            return [4.0, "WT N-term", 1]
        elif letter == "D":
            return [0.5, "WT N-term", 1]
        elif letter == "E":
            return [4.0, "H16N N-term", 1]
        elif letter == "F":
            return [2.0, "H16N N-term", 1]
        elif letter == "G":
            return [1.0, "H16N N-term", 1]
        elif letter == "H":
            return [2.0, "WT N-term", 1]

    elif (
        datafile
        == "/home/brianoj/Muteins 088-097 T/2019-04-19 IL2-088 IL2-097 Treg plate"
        or datafile
        == "/home/brianoj/Muteins 088-097 Nk/2019-04-19 IL2-088 IL2-097 NK plate"
    ):
        if letter == "A":
            return [4.0, "R38Q N-term", 1]
        elif letter == "B":
            return [2.0, "R38Q N-term", 1]
        elif letter == "C":
            return [1.0, "R38Q N-term", 1]
        elif letter == "D":
            return [1.0, "WT N-term", 1]
        elif letter == "E":
            return [4.0, "R38Q/H16N", 1]
        elif letter == "F":
            return [2.0, "R38Q/H16N", 1]
        elif letter == "G":
            return [1.0, "R38Q/H16N", 1]
        elif letter == "H":
            return [0.5, "R38Q/H16N", 1]

    elif (
        datafile
        == "/home/brianoj/Muteins 060-088 T/2019-05-02 IL2-060 IL2-088 Treg plate"
        or datafile
        == "/home/brianoj/Muteins 060-088 Nk/2019-05-02 IL2-060 IL2-088 NK plate"
    ):
        if letter == "A":
            return [4.0, "WT N-term", 1]
        elif letter == "B":
            return [4.0, "WT N-term", 1]
        elif letter == "C":
            return [4.0, "WT N-term", 1]
        elif letter == "D":
            return [0.5, "WT N-term", 1]
        elif letter == "E":
            return [4.0, "R38Q N-term", 1]
        elif letter == "F":
            return [2.0, "R38Q N-term", 1]
        elif letter == "G":
            return [1.0, "R38Q N-term", 1]
        elif letter == "H":
            return [2.0, "R38Q N-term", 1]

    elif (
        datafile
        == "/home/brianoj/Muteins 062-097 T/2019-05-02 IL2-062 IL2-097 Treg plate"
        or datafile
        == "/home/brianoj/Muteins 062-097 Nk/2019-05-02 IL2-062 IL2-097 NK plate"
    ):
        if letter == "A":
            return [4.0, "H16N N-term", 1]
        elif letter == "B":
            return [2.0, "H16N N-term", 1]
        elif letter == "C":
            return [1.0, "H16N N-term", 1]
        elif letter == "D":
            return [1.0, "H16N N-term", 1]
        elif letter == "E":
            return [4.0, "R38Q/H16N", 1]
        elif letter == "F":
            return [2.0, "R38Q/H16N", 1]
        elif letter == "G":
            return [1.0, "R38Q/H16N", 1]
        elif letter == "H":
            return [0.5, "R38Q/H16N", 1]
