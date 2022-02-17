"""
This creates Figure S6, comparing optimum selectivity for different IL2 formats.
"""
from email.mime import base
from os.path import dirname, join
from .figureCommon import getSetup
from ..imports import importCITE, importReceptors, getBindDict
from ..MBmodel import polyc, getKxStar
from scipy.optimize import minimize, Bounds, NonlinearConstraint
from copy import copy
import pandas as pd
import seaborn as sns
import numpy as np

path_here = dirname(dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((12, 10), (1, 1))

    epitopesDF = pd.DataFrame(columns={"Epitope", "Selectivity"})

    epitopesDF = pd.read_csv(join(path_here, "data/epitopeList.csv"))

    #epitopesDF = ['CD25','CD122']

    CITE_DF = importCITE()

    # Get conv factors, average them
    convFact = convFactCalc(ax[0])
    meanConv = convFact.Weight.mean()

    # weighting idea: take sample of everything of ~3 times and then average each types amount and use that as the size

    cellList = CITE_DF["CellType2"].unique().tolist()
    print(cellList)

    sampleSizes = []
    for cellType in cellList:
        cellSample = []
        for i in np.arange(3):
            sampleDF = CITE_DF.sample(1000)
            sampleSize = int(len(sampleDF.loc[sampleDF["CellType2"] == cellType]))
            cellSample.append(sampleSize)
        meanSize = np.mean(cellSample)
        sampleSizes.append(int(meanSize))

    offTCells = cellList.copy()
    offTCells.remove('Treg')

    # For each  cellType in list
    for i, cellType in enumerate(cellList):

        # Generate sample size
        sampleSize = sampleSizes[i]

        cellDF = CITE_DF.loc[CITE_DF["CellType2"] == cellType].sample(sampleSize)

        cellType_abdundances = []
        # For each epitope (being done on per cell basis)
        for e in epitopesDF.Epitope:
            # calculate abundance based on converstion factor
            if e == 'CD25':
                convFact = 77.136987
            elif e == 'CD122':
                convFact = 332.680090
            elif e == "CD127":
                convFact = 594.379215
            else:
                convFact = meanConv

            citeVal = cellDF[e].to_numpy()
            abundance = citeVal * convFact
            cellType_abdundances.append(abundance)
            # add column with this name to epitopesDF and abundances list

        epitopesDF[cellType] = cellType_abdundances

    # EpitopeDF now contains a data of single cell abundances for each cell type for each epitope
    epitopesDF['Selectivity'] = -1
    # New column which will hold selectivity per epitope

    targCell = 'Treg'

    standardDF = epitopesDF.loc[(epitopesDF.Epitope == 'CD25')].sample()
    standard2DF = epitopesDF.loc[(epitopesDF.Epitope == 'CD122')].sample()
    standardDF = standardDF.append(standard2DF)
    standardDF['Type'] = 'Standard'

    print(standardDF.columns())
    # For each epitope


        # New form
    optSelectivity = 1 / (optimizeDesign(targCell, offTCells, selectedDF, epitope))

    epitopesDF.loc[epitopesDF['Epitope'] == epitope, 'Selectivity'] = optSelectivity  # Store selectivity in DF to be used for plots

    baseSelectivity = 1 / (selecCalc(standardDF, targCell, offTCells))


    return f


def cytBindingModel(counts, x=False, date=False):
    """Runs binding model for a given mutein, valency, dose, and cell type."""
    mut = 'IL2'
    val = 1
    doseVec = np.array([0.1])
    recCount = np.ravel(counts)

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    Affs = mutAffDF.loc[(mutAffDF.Mutein == mut)]
    Affs = np.power(np.array([Affs["IL2RaKD"].values, Affs["IL2RBGKD"].values]) / 1e9, -1)
    Affs = np.reshape(Affs, (1, -1))
    Affs = np.repeat(Affs, 2, axis=0)
    np.fill_diagonal(Affs, 1e2)  # Each cytokine can only bind one a and one b

    if doseVec.size == 1:
        doseVec = np.array([doseVec])
    output = np.zeros(doseVec.size)

    for i, dose in enumerate(doseVec):
        if x:
            output[i] = polyc(dose / 1e9, np.power(10, x[0]), recCount, [[val, val]], [1.0], Affs)[0][1]
        else:
            output[i] = polyc(dose / 1e9, getKxStar(), recCount, [[val, val]], [1.0], Affs)[0][1]  # IL2RB binding only

    return output


def cytBindingModel_bispecOpt(counts, recXaff, x=False):
    """Runs binding model for a given mutein, valency, dose, and cell type."""

    mut = 'IL2'
    val = 1
    doseVec = np.array([0.1])

    recXaff = np.power(10, recXaff)

    recCount = np.ravel(counts)

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    Affs = mutAffDF.loc[(mutAffDF.Mutein == mut)]
    Affs = np.power(np.array([Affs["IL2RaKD"].values, Affs["IL2RBGKD"].values]) / 1e9, -1)
    Affs = np.reshape(Affs, (1, -1))
    Affs = np.append(Affs, recXaff)
    holder = np.full((3, 3), 1e2)
    np.fill_diagonal(holder, Affs)
    Affs = holder

    if doseVec.size == 1:
        doseVec = np.array([doseVec])
    output = np.zeros(doseVec.size)

    for i, dose in enumerate(doseVec):
        if x:
            output[i] = polyc(dose / (val * 1e9), np.power(10, x[0]), recCount, [[val, val, val]], [1.0], Affs)[0][1]
        else:
            output[i] = polyc(dose / (val * 1e9), getKxStar(), recCount, [[val, val, val]], [1.0], Affs)[0][1]  # IL2RB binding only

    return output


def selecCalc(df, targCell, offTCells):
    """Calculates selectivity for no additional epitope"""
    targetBound = 0
    offTargetBound = 0

    cd25DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD25')]
    cd122DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD122')]

    for i, cd25Count in enumerate(cd25DF[targCell].item()):
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count]
        targetBound += cytBindingModel(counts)
    for cellT in offTCells:
        for i, cd25Count in enumerate(cd25DF[cellT].item()):
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count]
            offTargetBound += cytBindingModel(counts)

    return (offTargetBound) / (targetBound)


def minSelecFunc(x, selectedDF, targCell, offTCells, epitope):
    """Provides the function to be minimized to get optimal selectivity"""
    targetBound = 0
    offTargetBound = 0

    recXaff = x

    epitopeDF = selectedDF.loc[(selectedDF.Type == 'Epitope')]
    cd25DF = selectedDF.loc[(selectedDF.Type == 'Standard') & (selectedDF.Epitope == 'CD25')]
    cd122DF = selectedDF.loc[(selectedDF.Type == 'Standard') & (selectedDF.Epitope == 'CD122')]

    for i, epCount in enumerate(epitopeDF[targCell].item()):
        cd25Count = cd25DF[targCell].item()[i]
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count, epCount]
        targetBound += cytBindingModel_bispecOpt(counts, recXaff)
    for cellT in offTCells:
        for i, epCount in enumerate(epitopeDF[cellT].item()):
            cd25Count = cd25DF[cellT].item()[i]
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count, epCount]

            offTargetBound += cytBindingModel_bispecOpt(counts, recXaff)

    return (offTargetBound) / (targetBound)


def optimizeDesign(targCell, offTcells, selectedDF, epitope):
    """ A more general purpose optimizer """

    if targCell == "NK":
        X0 = [6.0, 8]
    else:
        X0 = [7.0]

    optBnds = Bounds(np.full_like(X0, 6.0), np.full_like(X0, 9.0))

    optimized = minimize(minSelecFunc, X0, bounds=optBnds, args=(selectedDF, targCell, offTcells, epitope), jac="3-point")
    optSelectivity = optimized.fun[0]

    return optSelectivity


cellDict = {"CD4 Naive": "Thelper",
            "CD4 CTL": "Thelper",
            "CD4 TCM": "Thelper",
            "CD4 TEM": "Thelper",
            "NK": "NK",
            "CD8 Naive": "CD8",
            "CD8 TCM": "CD8",
            "CD8 TEM": "CD8",
            "Treg": "Treg"}


markDict = {"CD25": "IL2Ra",
            "CD122": "IL2Rb",
            "CD127": "IL7Ra",
            "CD132": "gc"}


def convFactCalc(ax):
    """Fits a ridge classifier to the CITE data and plots those most highly correlated with T reg"""
    CITE_DF = importCITE()
    cellToI = ["CD4 TCM", "CD8 Naive", "NK", "CD8 TEM", "CD4 Naive", "CD4 CTL", "CD8 TCM", "Treg", "CD4 TEM"]
    markers = ["CD122", "CD127", "CD25"]
    markerDF = pd.DataFrame(columns=["Marker", "Cell Type", "Amount", "Number"])
    for marker in markers:
        for cell in cellToI:
            cellTDF = CITE_DF.loc[CITE_DF["CellType2"] == cell][marker]
            markerDF = markerDF.append(pd.DataFrame({"Marker": [marker], "Cell Type": cell, "Amount": cellTDF.mean(), "Number": cellTDF.size}))

    markerDF = markerDF.replace({"Marker": markDict, "Cell Type": cellDict})
    markerDFw = pd.DataFrame(columns=["Marker", "Cell Type", "Average"])
    for marker in markerDF.Marker.unique():
        for cell in markerDF["Cell Type"].unique():
            subDF = markerDF.loc[(markerDF["Cell Type"] == cell) & (markerDF["Marker"] == marker)]
            wAvg = np.sum(subDF.Amount.values * subDF.Number.values) / np.sum(subDF.Number.values)
            markerDFw = markerDFw.append(pd.DataFrame({"Marker": [marker], "Cell Type": cell, "Average": wAvg}))

    recDF = importReceptors()
    weightDF = pd.DataFrame(columns=["Receptor", "Weight"])

    for rec in markerDFw.Marker.unique():
        CITEval = np.array([])
        Quantval = np.array([])
        for cell in markerDF["Cell Type"].unique():
            CITEval = np.concatenate((CITEval, markerDFw.loc[(markerDFw["Cell Type"] == cell) & (markerDFw["Marker"] == rec)].Average.values))
            Quantval = np.concatenate((Quantval, recDF.loc[(recDF["Cell Type"] == cell) & (recDF["Receptor"] == rec)].Mean.values))
        weightDF = weightDF.append(pd.DataFrame({"Receptor": [rec], "Weight": np.linalg.lstsq(np.reshape(CITEval, (-1, 1)), Quantval, rcond=None)[0]}))

    return weightDF
