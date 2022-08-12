"""File that deals with everything about importing and sampling."""
import os
from functools import lru_cache
from os.path import join
import numpy as np
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import coo_matrix
from zipfile import ZipFile

path_here = os.path.dirname(os.path.dirname(__file__))


@lru_cache(maxsize=None)
def import_pstat(combine_samples=True):
    """ Loads CSV file containing pSTAT5 levels from Visterra data. Incorporates only Replicate 1 since data missing in Replicate 2. """
    path = os.path.dirname(os.path.dirname(__file__))
    data = np.array(pd.read_csv(join(path, "ckine/data/pSTAT_data.csv"), encoding="latin1"))
    ckineConc = data[4, 2:14]
    tps = np.array([0.5, 1.0, 2.0, 4.0]) * 60.0
    # 4 time points, 10 cell types, 12 concentrations, 2 replicates
    IL2_data = np.zeros((40, 12))
    IL2_data2 = IL2_data.copy()
    IL15_data = IL2_data.copy()
    IL15_data2 = IL2_data.copy()
    cell_names = list()
    for i in range(10):
        cell_names.append(data[12 * i + 3, 1])
        # Subtract the zero treatment plates before assigning to returned arrays
        if i <= 4:
            zero_treatment = data[12 * (i + 1), 13]
            zero_treatment2 = data[8 + (12 * i), 30]
        else:
            zero_treatment = data[8 + (12 * i), 13]
            zero_treatment2 = data[8 + (12 * i), 30]
        # order of increasing time by cell type
        IL2_data[4 * i: 4 * (i + 1), :] = np.flip(data[6 + (12 * i): 10 + (12 * i), 2:14].astype(float) - zero_treatment, 0)
        IL2_data2[4 * i: 4 * (i + 1), :] = np.flip(data[6 + (12 * i): 10 + (12 * i), 19:31].astype(float) - zero_treatment2, 0)
        IL15_data[4 * i: 4 * (i + 1), :] = np.flip(data[10 + (12 * i): 14 + (12 * i), 2:14].astype(float) - zero_treatment, 0)
        IL15_data2[4 * i: 4 * (i + 1), :] = np.flip(data[10 + (12 * i): 14 + (12 * i), 19:31].astype(float) - zero_treatment2, 0)

    if combine_samples is False:
        return ckineConc, cell_names, IL2_data, IL2_data2, IL15_data, IL15_data2

    for i in range(IL2_data.shape[0]):
        for j in range(IL2_data.shape[1]):
            # take average of both replicates if specific entry isn't nan
            IL2_data[i, j] = np.nanmean(np.array([IL2_data[i, j], IL2_data2[i, j]]))
            IL15_data[i, j] = np.nanmean(np.array([IL15_data[i, j], IL15_data2[i, j]]))

    dataMean = pd.DataFrame(
        {
            "Cells": np.tile(np.repeat(cell_names, 48), 2),
            "Ligand": np.concatenate((np.tile(np.array("IL2"), 480), np.tile(np.array("IL15"), 480))),
            "Time": np.tile(np.repeat(tps, 12), 20),
            "Concentration": np.tile(ckineConc, 80),
            "RFU": np.concatenate((IL2_data.reshape(480), IL15_data.reshape(480))),
        }
    )

    return ckineConc, cell_names, IL2_data, IL15_data, dataMean


# Receptor Quant - Beads (4/23 & 4/26)


channels = {}
channels["A"] = ["VL1-H", "BL5-H", "RL1-H", "RL1-H", "RL1-H", "Width"]
channels["C"] = ["VL4-H", "VL6-H", "BL1-H", "BL3-H"]
channels["D"] = ["VL1-H", "VL1-H", "VL1-H", "VL1-H", "VL1-H"]
channels["E"] = ["VL6-H", "BL3-H", "BL5-H", "BL5-H", "BL5-H", "BL5-H", "BL5-H"]
channels["F"] = channels["G"] = channels["H"] = ["RL1-H", "RL1-H", "RL1-H", "RL1-H", "RL1-H"]
channels["I"] = ["BL1-H", "BL1-H", "BL1-H", "BL1-H", "BL1-H"]

receptors = {}
receptors["A"] = ["CD25", "CD122", "CD132", "IL15(1)", "IL15(2)", " "]
receptors["C"] = ["CD3", "CD4", "CD127", "CD45RA"]
receptors["D"] = ["CD25", "CD25", "CD25", "CD25", "CD25"]
receptors["E"] = ["CD8", "CD56", "CD122", "CD122", "CD122", "CD122", "CD122"]
receptors["F"] = ["CD132", "CD132", "CD132", "CD132", "CD132"]
receptors["G"] = ["IL15(1)", "IL15(1)", "IL15(1)", "IL15(1)", "IL15(1)"]
receptors["H"] = ["IL15(2)", "IL15(2)", "IL15(2)", "IL15(2)", "IL15(2)"]
receptors["I"] = ["CD127", "CD127", "CD127", "CD127", "CD127"]


@lru_cache(maxsize=None)
def import_pstat_all(singleCell=False, updateLigs=True):
    """ Loads CSV file containing all WT and Mutein pSTAT responses and moments"""
    WTbivDF = pd.read_csv(join(path_here, "ckine/data/WTDimericMutSingleCellData.csv"), encoding="latin1")
    monDF = pd.read_csv(join(path_here, "ckine/data/MonomericMutSingleCellData.csv"), encoding="latin1")
    respDF = pd.concat([WTbivDF, monDF])
    if singleCell:
        WTbivDFbin = pd.read_csv(join(path_here, "ckine/data/WTDimericMutSingleCellDataBin.csv"), encoding="latin1")
        monDFbin = pd.read_csv(join(path_here, "ckine/data/MonomericMutSingleCellDataBin.csv"), encoding="latin1")
        respDFbin = pd.concat([WTbivDFbin, monDFbin])
        respDFbin = respDFbin.loc[respDFbin["Bin"].isin([1, 3])]
        respDFbin.loc[respDFbin["Bin"] == 1, "Cell"] += r" $IL2Ra^{lo}$"
        respDFbin.loc[respDFbin["Bin"] == 3, "Cell"] += r" $IL2Ra^{hi}$"
        respDF = pd.concat([respDF, respDFbin])

    if updateLigs:
        respDF.loc[(respDF.Bivalent == 0), "Ligand"] = (respDF.loc[(respDF.Bivalent == 0)].Ligand + " (Mono)").values
        respDF.loc[(respDF.Bivalent == 1), "Ligand"] = (respDF.loc[(respDF.Bivalent == 1)].Ligand + " (Biv)").values

    return respDF


@lru_cache(maxsize=None)
def import_pstat_all_meyer():
    """ Loads CSV file containing all WT and Mutein pSTAT responses and moments"""
    respDF = pd.read_csv(join(path_here, "ckine/data/Meyer_Flow.csv"), encoding="latin1")
    respDF = respDF.loc[respDF.Date != "7/22/22"]
    return respDF


@lru_cache(maxsize=None)
def getBindDict():
    """Gets binding to pSTAT fluorescent conversion dictionary"""
    path = os.path.dirname(os.path.dirname(__file__))
    bindingDF = pd.read_csv(join(path, "ckine/data/BindingConvDict.csv"), encoding="latin1")
    return bindingDF


@lru_cache(maxsize=None)
def importReceptors():
    """Makes Complete receptor expression Dict"""
    path = os.path.dirname(os.path.dirname(__file__))
    recDF = pd.read_csv(join(path_here, "ckine/data/RecQuantitation.csv"))
    recDFbin = pd.read_csv(join(path_here, "ckine/data/BinnedReceptorData.csv"))
    recDFbin = recDFbin.loc[recDFbin["Bin"].isin([1, 3])]
    recDFbin.loc[recDFbin["Bin"] == 1, "Cell Type"] += r" $IL2Ra^{lo}$"
    recDFbin.loc[recDFbin["Bin"] == 3, "Cell Type"] += r" $IL2Ra^{hi}$"
    recDF = pd.concat([recDF, recDFbin])
    return recDF


def importCITE():
    """Downloads all surface markers and cell types"""
    CITEmarkerDF = pd.read_csv(join(path_here, "ckine/data/CITEdata_SurfMarkers.zip"))
    return CITEmarkerDF


def makeRNAseqDF():
    """Makes surface RNAseq DF"""
    matrix = mmread(join(path_here, "ckine/data/GSM5008737_RNA_3P-matrix.mtx.gz"))
    surfGenes = pd.read_csv("ckine/data/SurfaceGenes.csv")
    featuresGenes = pd.read_csv("ckine/data/RNAfeatures.csv")
    surfaceList = surfGenes["Gene"].values
    allFeatures = featuresGenes["Gene"].values
    surfInd = np.isin(allFeatures, surfaceList)
    cols = [i for i, x in enumerate(surfInd) if x]
    dataCoords = np.isin(matrix.row, cols)
    locList = np.where(surfInd == True)[0].tolist()
    newCoords = np.arange(0, len(locList)).tolist()

    colDict = {}
    for key in locList:
        for value in newCoords:
            colDict[key] = value
            newCoords.remove(value)
            break

    def vec_translate(a, my_dict):
        return np.vectorize(my_dict.__getitem__)(a)
    matrix2 = coo_matrix((matrix.data[dataCoords], (matrix.col[dataCoords], vec_translate(matrix.row[dataCoords], colDict))), shape=(matrix.shape[1], np.count_nonzero(surfInd)))
    geneCols = allFeatures[surfInd]
    GeneDF = pd.DataFrame(data=matrix2.toarray(), columns=geneCols)
    cellTypesDF = pd.read_csv(join(path_here, "ckine/data/CITEcellTypes.csv"))
    GeneDF = pd.concat([GeneDF, cellTypesDF], axis=1)
    GeneDF.to_csv(join(path_here, "ckine/data/RNAseqSurface.csv.zip"))


def importRNACITE():
    """Downloads all surface markers and cell types"""
    RNAsurfDF = pd.read_csv(ZipFile(join(path_here, "ckine/data/RNAseqSurface.csv.zip")).open("RNAseqSurface.csv"))
    return RNAsurfDF
