"""
This creates Figure S6, plotting Treg to off target signaling for vaying mutein dose for different IL2 formats
"""
from email.mime import base
from os.path import dirname, join

from matplotlib.pyplot import xlim,ylim
from .figureCommon import getSetup
from ..imports import importCITE, importReceptors
from ..MBmodel import polyc, getKxStar
from copy import copy
import pandas as pd
import seaborn as sns
import numpy as np

path_here = dirname(dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((26, 8), (1, 3))

    receoptors = {'Epitope':['CD25','CD122']}
    epitopesDF = pd.DataFrame(receoptors)


    CITE_DF = importCITE()

    # Get conv factors, average them
    convFact = convFactCalc(ax[0])
    meanConv = convFact.Weight.mean()

    # weighting idea: take sample of everything of ~3 times and then average each types amount and use that as the size

    cellList = CITE_DF["CellType2"].unique().tolist()

    sampleSizes = []
    for cellType in cellList:
        cellSample = []
        for i in np.arange(3):
            sampleDF = CITE_DF.sample(1000)
            sampleSize = int(len(sampleDF.loc[sampleDF["CellType2"] == cellType]))
            cellSample.append(sampleSize)
        meanSize = np.mean(cellSample)
        sampleSizes.append(int(meanSize))

    #offTCells = cellList.copy()
    #offTCells.remove('Treg')
    

    offTCells = ['CD8 Naive','NK', 'CD8 TEM','CD8 Proliferating','NK Proliferating','NK_CD56bright']

    print(offTCells)

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
            else:
                assert(False)
                convFact = meanConv

            citeVal = cellDF[e].to_numpy()
            abundance = citeVal * convFact
            cellType_abdundances.append(abundance)
            # add column with this name to epitopesDF and abundances list

        epitopesDF[cellType] = cellType_abdundances

    # New column which will hold selectivity per epitope

    targCell = 'Treg'

    standardDF = epitopesDF.loc[(epitopesDF.Epitope == 'CD25')].sample()
    standard2DF = epitopesDF.loc[(epitopesDF.Epitope == 'CD122')].sample()
    standardDF = standardDF.append(standard2DF)
    standardDF['Type'] = 'Standard'


    #range from 0.01 <-> 100
    doseVec = np.logspace(-12,-6,10)

    treg_sigs = np.zeros((8,20))
    offTarg_sigs = np.zeros((8,20))

    #0-2 IL2 WT
    #3-5 R38Q
    #6-7 Live/Dead
    muts = ['IL2', 'R38Q/H16N']
    vals = [1,2,4]

    for i, aff in enumerate(doseVec):
        print(aff)
        for j, mut in enumerate(muts):
            for k, val in enumerate(vals):
                n = (3*j)+k
                treg_sig, offTarg_sig = bindingCalc(standardDF, targCell, offTCells, aff, val, mut)
                treg_sigs[n,i] = treg_sig
                offTarg_sigs[n,i] = offTarg_sig

        treg_sig_bi, offTarg_sig_bi = bindingCalc_bispec(standardDF, targCell, offTCells, aff, 1)
        treg_sigs[6,i] = treg_sig_bi
        offTarg_sigs[6,i] = offTarg_sig_bi

        treg_sig_bi, offTarg_sig_bi = bindingCalc_bispec(standardDF, targCell, offTCells, aff, 2)
        treg_sigs[7,i] = treg_sig_bi
        offTarg_sigs[7,i] = offTarg_sig_bi

        #print(treg_sigs)

    # Normalizes data to 1
    def norm(data):
        return data/max(data)

    ##print(y_ticks)
    def plotSignals(types,ax):
        #Add standard colors/line types
        if 'WT' in types:
            ax.plot(norm(treg_sigs[0]),norm(offTarg_sigs[0]),label='WT', c='blue')
            ax.plot(norm(treg_sigs[1]),norm(offTarg_sigs[1]),label='WT Bival', c='green')
            ax.plot(norm(treg_sigs[2]),norm(offTarg_sigs[2]),label='WT Tetraval', c='c')
        if 'R38Q/H16N' in types:
            ax.plot(norm(treg_sigs[3]),norm(offTarg_sigs[3]),label='R38Q/H16N', c='red')
            ax.plot(norm(treg_sigs[4]),norm(offTarg_sigs[4]),label='R38Q/H16N Bival', c='y')
            ax.plot(norm(treg_sigs[5]),norm(offTarg_sigs[5]),label='R38Q/H16N Tetraval', c='orange')
        if 'Live/Dead' in types:
            ax.plot(norm(treg_sigs[6]),norm(offTarg_sigs[6]),label='CD25 Live/Dead',c='indigo')
            ax.plot(norm(treg_sigs[7]),norm(offTarg_sigs[7]),label='CD25 Bivalent Live/Dead',c='magenta')
        
        ax.set(xlabel='Treg Signaling',ylabel='Off Target Signaling')
        ax.legend()

    
    plotSignals(['WT','R38Q/H16N'],ax[0])
    plotSignals(['WT','Live/Dead'],ax[1])
    plotSignals(['R38Q/H16N','Live/Dead'],ax[2])

    return f


def cytBindingModel(counts, doseVec, val, mut, x=False, date=False):
    """Runs binding model for a given mutein, valency, dose, and cell type."""
    #mut = mut
    #doseVec = np.array([0.1])
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
            output[i] = polyc(dose / 1e9, getKxStar(), recCount, [[val, val]], [1.0], Affs)[0][1] 

    return output

def bindingCalc(df, targCell, offTCells, doseVec, val, mut):
    """Calculates selectivity for no additional epitope"""
    targetBound = 0
    offTargetBound = 0

    cd25DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD25')]
    cd122DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD122')]

    for i, cd25Count in enumerate(cd25DF[targCell].item()):
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count]
        targetBound += cytBindingModel(counts, doseVec, val, mut)

    for cellT in offTCells:
        for i, cd25Count in enumerate(cd25DF[cellT].item()):
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count]
            offTargetBound += cytBindingModel(counts, doseVec, val, mut)

    return targetBound, offTargetBound


def cytBindingModel_bispec(counts, doseVec, recXaff, val, x=False):
    """Runs binding model for a given mutein, valency, dose, and cell type."""

    mut = 'R38Q/H16N' 
    #doseVec = np.array([0.1])

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
            output[i] = polyc(dose / (val * 1e9), getKxStar(), recCount, [[val, val, val]], [1.0], Affs)[0][1] 

    return output


def bindingCalc_bispec(df, targCell, offTCells,doseVec,val):
    """Calculates selectivity for no additional epitope"""
    targetBound = 0
    offTargetBound = 0

    cd25DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD25')]
    cd122DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD122')]

    for i, cd25Count in enumerate(cd25DF[targCell].item()):
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count,cd25Count]
        targetBound += cytBindingModel_bispec(counts, doseVec, 9.0, val)

    for cellT in offTCells:
        for i, cd25Count in enumerate(cd25DF[cellT].item()):
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count,cd25Count]
            offTargetBound += cytBindingModel_bispec(counts, doseVec, 9.0, val)

    return targetBound, offTargetBound

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
