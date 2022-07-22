"""Stores Single Cell marker levels for monomeric signaling experiments"""
import os
import pandas as pds
import numpy as np
from scipy import stats
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
        "ckine/data/2019-11-08 monomer IL-2 Fc signaling/CD4 T cells - IL2-060 mono, IL2-060 dimeric",
        "ckine/data/2019-11-08 monomer IL-2 Fc signaling/CD4 T cells - IL2-062 mono, IL2-118 mono",
        "ckine/data/2019-11-27 monomer IL-2 Fc signaling/CD4 T cells - C-term IL2-060 mono, C-term V91K mono",
        "ckine/data/2019-12-05 monomer IL-2 Fc signaling/CD4 T cells - IL2-109 mono, IL2-118 mono",
        "ckine/data/2019-12-05 monomer IL-2 Fc signaling/CD4 T cells - IL2-110 mono, C-term N88D mono",
        "ckine/data/2019-11-08 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-060 mono, IL2-060 dimeric",
        "ckine/data/2019-11-08 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-062 mono, IL2-118 mono",
        "ckine/data/2019-11-27 monomer IL-2 Fc signaling/NK CD8 T cells - C-term IL2-060 mono, C-term V91K mono",
        "ckine/data/2019-12-05 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-109 mono, IL2-118 mono",
        "ckine/data/2019-12-05 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-110 mono, C-term N88D mono",
    ]
    dates = [
        "11/8/2019",
        "11/8/2019",
        "11/27/2019",
        "12/5/2019",
        "12/5/2019",
        "11/8/2019",
        "11/8/2019",
        "11/27/2019",
        "12/5/2019",
        "12/5/2019",
    ]
    repList = [0, 1, 0, 0, 1, 0, 1, 0, 0, 1]
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
    print("Starting Muteins")

    for i, filename in enumerate(dataFiles):
        if i < 5:
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
                        dat_array.loc[cellTypeDict[cell][ii].values, "Cell Type"] = cell
                    MVdf = MVdf.append(dat_array)

                if j == 3 or j == 7:
                    MVdf["pSTAT5"] = (
                        MVdf["pSTAT5"] - MVdf.loc[(MVdf.Dose <= 0.001423)].pSTAT5.mean()
                    )
                    masterMVdf = masterMVdf.append(MVdf)
                    MVdf = pds.DataFrame()

    masterMVdf.pSTAT5 = masterMVdf.pSTAT5.clip(lower=0)
    masterMVdf = masterMVdf.loc[masterMVdf.Ligand != "H16L N-term"]
    masterMVdf.reset_index().to_feather("FlowDataGMM_Mon.ftr")

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


# done done


def mutFunc(letter, datafile):
    if (
        datafile
        == "ckine/data/2019-11-08 monomer IL-2 Fc signaling/CD4 T cells - IL2-060 mono, IL2-060 dimeric"
        or datafile
        == "ckine/data/2019-11-08 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-060 mono, IL2-060 dimeric"
    ):
        if letter == "A":
            return [4.0, "WT N-term", 0]
        elif letter == "B":
            return [2.0, "WT N-term", 0]
        elif letter == "C":
            return [1.0, "WT N-term", 0]
        elif letter == "D":
            return [0.5, "WT N-term", 0]
        elif letter == "E":
            return [4.0, "WT N-term", 1]
        elif letter == "F":
            return [2.0, "WT N-term", 1]
        elif letter == "G":
            return [1.0, "WT N-term", 1]
        elif letter == "H":
            return [0.5, "WT N-term", 1]
    # done done
    elif (
        datafile
        == "ckine/data/2019-11-08 monomer IL-2 Fc signaling/CD4 T cells - IL2-062 mono, IL2-118 mono"
        or datafile
        == "ckine/data/2019-11-08 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-062 mono, IL2-118 mono"
    ):
        if letter == "A":
            return [4.0, "H16N N-term", 0]
        elif letter == "B":
            return [2.0, "H16N N-term", 0]
        elif letter == "C":
            return [1.0, "H16N N-term", 0]
        elif letter == "D":
            return [0.5, "H16N N-term", 0]
        elif letter == "E":
            return [4.0, "H16L N-term", 0]
        elif letter == "F":
            return [2.0, "H16L N-term", 0]
        elif letter == "G":
            return [1.0, "H16L N-term", 0]
        elif letter == "H":
            return [0.5, "H16L N-term", 0]
    # done done
    elif (
        datafile
        == "ckine/data/2019-11-27 monomer IL-2 Fc signaling/CD4 T cells - C-term IL2-060 mono, C-term V91K mono"
        or datafile
        == "ckine/data/2019-11-27 monomer IL-2 Fc signaling/NK CD8 T cells - C-term IL2-060 mono, C-term V91K mono"
    ):
        if letter == "A":
            return [4.0, "WT C-term", 0]
        elif letter == "B":
            return [2.0, "WT C-term", 0]
        elif letter == "C":
            return [1.0, "WT C-term", 0]
        elif letter == "D":
            return [0.5, "WT C-term", 0]
        elif letter == "E":
            return [4.0, "V91K C-term", 0]
        elif letter == "F":
            return [2.0, "V91K C-term", 0]
        elif letter == "G":
            return [1.0, "V91K C-term", 0]
        elif letter == "H":
            return [0.5, "V91K C-term", 0]
    # done not done
    elif (
        datafile
        == "ckine/data/2019-12-05 monomer IL-2 Fc signaling/CD4 T cells - IL2-109 mono, IL2-118 mono"
        or datafile
        == "ckine/data/2019-12-05 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-109 mono, IL2-118 mono"
    ):
        if letter == "A":
            return [4.0, "R38Q N-term", 0]
        elif letter == "B":
            return [2.0, "R38Q N-term", 0]
        elif letter == "C":
            return [1.0, "R38Q N-term", 0]
        elif letter == "D":
            return [0.5, "R38Q N-term", 0]
        elif letter == "E":
            return [4.0, "H16L N-term", 0]
        elif letter == "F":
            return [2.0, "H16L N-term", 0]
        elif letter == "G":
            return [1.0, "H16L N-term", 0]
        elif letter == "H":
            return [0.5, "H16L N-term", 0]
    # done
    elif (
        datafile
        == "ckine/data/2019-12-05 monomer IL-2 Fc signaling/CD4 T cells - IL2-110 mono, C-term N88D mono"
        or datafile
        == "ckine/data/2019-12-05 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-110 mono, C-term N88D mono"
    ):
        if letter == "A":
            return [4.0, "F42Q N-Term", 0]
        elif letter == "B":
            return [2.0, "F42Q N-Term", 0]
        elif letter == "C":
            return [1.0, "F42Q N-Term", 0]
        elif letter == "D":
            return [0.5, "F42Q N-Term", 0]
        elif letter == "E":
            return [4.0, "N88D C-term", 0]
        elif letter == "F":
            return [2.0, "N88D C-term", 0]
        elif letter == "G":
            return [1.0, "N88D C-term", 0]
        elif letter == "H":
            return [0.5, "N88D C-term", 0]
