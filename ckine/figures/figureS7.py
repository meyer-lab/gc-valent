"""
This creates Figure S6, plotting Treg to off target signaling for vaying mutein dose for different IL2 formats
"""
from email.mime import base
from os.path import dirname, join
from matplotlib.pyplot import xlim, ylim
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
    ax, f = getSetup((13, 4), (1, 3))

    receptors = {'Epitope': ['CD25', 'CD122']}
    epitopesDF = pd.DataFrame(receptors)

    CITE_DF = importCITE()
    cellList = ['CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 'CD8 TEM', 'CD4 Naive', 'CD4 TEM', 'CD4 TCM', 'NK', 'NK Proliferating', 'Treg']

    sampleSizes = []
    for cellType in cellList:
        cellSample = []
        for i in np.arange(10):
            sampleDF = CITE_DF.sample(1000)
            sampleSize = int(len(sampleDF.loc[sampleDF["CellType2"] == cellType]))
            cellSample.append(sampleSize)
        meanSize = np.mean(cellSample)
        sampleSizes.append(int(meanSize))

    offTCells = ['CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 'CD8 TEM', 'CD4 Naive', 'CD4 TEM', 'CD4 TCM', 'NK', 'NK Proliferating']

    # For each  cellType in list
    for i, cellType in enumerate(cellList):

        # Generate sample size
        sampleSize = sampleSizes[i]

        cellDF = CITE_DF.loc[CITE_DF["CellType2"] == cellType].sample(sampleSize, random_state=45)  # 45 is okay

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
    print(epitopesDF)

    targCell = 'Treg'

    standardDF = epitopesDF.loc[(epitopesDF.Epitope == 'CD25')].sample()
    standard2DF = epitopesDF.loc[(epitopesDF.Epitope == 'CD122')].sample()
    standardDF = standardDF.append(standard2DF)
    standardDF['Type'] = 'Standard'

    # range from pico <-> micromolar
    doseVec = np.logspace(-18, -6, 60)

    treg_sigs = np.zeros((8, 60))
    offTarg_sigs = np.zeros((8, 60))

    # 0-2 IL2 WT
    # 3-5 R38Q
    # 6-7 Live/Dead
    muts = ['IL2', 'R38Q/H16N']
    vals = [1, 2, 4]

    for i, dose in enumerate(doseVec):
        print(dose)
        for j, mut in enumerate(muts):
            for k, val in enumerate(vals):
                n = (3 * j) + k
                treg_sig, offTarg_sig = bindingCalc(standardDF, targCell, offTCells, dose, val, mut)
                treg_sigs[n, i] = treg_sig
                offTarg_sigs[n, i] = offTarg_sig

        treg_sig_bi, offTarg_sig_bi = bindingCalc_bispec(standardDF, targCell, offTCells, dose, 1)
        treg_sigs[6, i] = treg_sig_bi
        offTarg_sigs[6, i] = offTarg_sig_bi

        treg_sig_bi, offTarg_sig_bi = bindingCalc_bispec(standardDF, targCell, offTCells, dose, 2)
        treg_sigs[7, i] = treg_sig_bi
        offTarg_sigs[7, i] = offTarg_sig_bi

    tregMax = np.amax(treg_sigs)
    offMax = np.amax(offTarg_sigs)
    print(tregMax)
    print(offMax)

    print(treg_sigs)
    print(offTarg_sigs)

    # Normalizes data to 1
    def norm(data, maxVal):
        copData = copy(data)
        return copData / maxVal

    def plotSignals(types, ax):
        # Add standard colors/line types
        if 'WT' in types:
            ax.plot(norm(offTarg_sigs[0], offMax), norm(treg_sigs[0], tregMax), label='WT', c='blue')
            ax.plot(norm(offTarg_sigs[1], offMax), norm(treg_sigs[1], tregMax), label='WT Bival', c='green')
            ax.plot(norm(offTarg_sigs[2], offMax), norm(treg_sigs[2], tregMax), label='WT Tetraval', c='c')
        if 'R38Q/H16N' in types:
            ax.plot(norm(offTarg_sigs[3], offMax), norm(treg_sigs[3], tregMax), '--', label='R38Q/H16N', c='red')
            ax.plot(norm(offTarg_sigs[4], offMax), norm(treg_sigs[4], tregMax), '--', label='R38Q/H16N Bival', c='y')
            ax.plot(norm(offTarg_sigs[5], offMax), norm(treg_sigs[5], tregMax), '--', label='R38Q/H16N Tetraval', c='orange')
        if 'Live/Dead' in types:
            ax.plot(norm(offTarg_sigs[6], offMax), norm(treg_sigs[6], tregMax), '-.', label='CD25 Live/Dead', c='indigo')
            ax.plot(norm(offTarg_sigs[7], offMax), norm(treg_sigs[7], tregMax), '-.', label='CD25 Bivalent Live/Dead', c='magenta')

        #ax.set(xlabel='Treg Signaling',ylabel='Off Target Signaling')
        ax.set_xlabel('Off Target Signaling', fontsize=12)
        ax.set_ylabel('Treg Signaling', fontsize=12)
        ax.legend()

    plotSignals(['WT', 'R38Q/H16N'], ax[0])
    plotSignals(['WT', 'Live/Dead'], ax[1])
    plotSignals(['R38Q/H16N', 'Live/Dead'], ax[2])
    f.suptitle('Treg vs. Off Target Signaling Varing Dose Concentration', fontsize=18)

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
            output[i] = polyc(dose / (val), np.power(10, x[0]), recCount, [[val, val]], [1.0], Affs)[0][1]
        else:
            output[i] = polyc(dose / (val), getKxStar(), recCount, [[val, val]], [1.0], Affs)[0][1]

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
    Affs = np.power(np.array([Affs["IL2RaKD"].values, [1]]) / 1e9, -1)
    Affs = np.reshape(Affs, (1, -1))
    Affs = np.append(Affs, recXaff)
    holder = np.full((3, 2), 1e2)
    holder[0, 0] = Affs[0]
    holder[1, 1] = Affs[1]
    holder[2, 0] = Affs[2]
    Affs = holder

    if doseVec.size == 1:
        doseVec = np.array([doseVec])
    output = np.zeros(doseVec.size)

    for i, dose in enumerate(doseVec):
        if x:
            output[i] = polyc(dose / (val), np.power(10, x[0]), recCount, [[val, val, val]], [1.0], Affs)[0][1]
        else:
            output[i] = polyc(dose / (val), getKxStar(), recCount, [[val, val, val]], [1.0], Affs)[0][1]

    return output


def bindingCalc_bispec(df, targCell, offTCells, doseVec, val):
    """Calculates selectivity for no additional epitope"""
    targetBound = 0
    offTargetBound = 0

    cd25DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD25')]
    cd122DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD122')]

    for i, cd25Count in enumerate(cd25DF[targCell].item()):
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count]
        targetBound += cytBindingModel_bispec(counts, doseVec, 9, val)

    for cellT in offTCells:
        for i, cd25Count in enumerate(cd25DF[cellT].item()):
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count]
            offTargetBound += cytBindingModel_bispec(counts, doseVec, 9, val)

    return targetBound, offTargetBound
