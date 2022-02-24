"""
This creates Figure S6, plotting Treg to off target signaling for vaying IL2Rb affinity for different IL2 formats
"""
from email.mime import base
from os.path import dirname, join
from scipy.optimize import minimize, Bounds
from .figureCommon import getSetup
from ..imports import importCITE
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

    

    #offTCells = ['CD8 Naive','NK', 'CD8 TEM','CD8 TCM','CD8 Proliferating','NK Proliferating','NK_CD56bright']
    offTCells = ['CD8 Naive','NK', 'CD8 TEM','CD8 TCM','CD8 ']

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
    betaAffs = np.logspace(-4,2,40)

    treg_sigs = np.zeros((8,40))
    offTarg_sigs = np.zeros((8,40))

    #0-2 IL2 WT
    #3-5 R38Q
    #6-7 Live/Dead
    muts = ['IL2', 'R38Q/H16N']
    vals = [1,2,4]

    

    for i, aff in enumerate(betaAffs):
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
            ax.plot(norm(treg_sigs[3]),norm(offTarg_sigs[3]),'--',label='R38Q/H16N', c='red')
            ax.plot(norm(treg_sigs[4]),norm(offTarg_sigs[4]),'--',label='R38Q/H16N Bival', c='y')
            ax.plot(norm(treg_sigs[5]),norm(offTarg_sigs[5]),'--',label='R38Q/H16N Tetraval', c='orange')
        if 'Live/Dead' in types:
            ax.plot(norm(treg_sigs[6]),norm(offTarg_sigs[6]),'-.',label='CD25 Live/Dead',c='indigo')
            ax.plot(norm(treg_sigs[7]),norm(offTarg_sigs[7]),'-.',label='CD25 Bivalent Live/Dead',c='magenta')
        
        ax.set(xlabel='Treg Signaling',ylabel='Off Target Signaling')
        ax.legend()

    
    plotSignals(['WT','R38Q/H16N'],ax[0])
    plotSignals(['WT','Live/Dead'],ax[1])
    plotSignals(['R38Q/H16N','Live/Dead'],ax[2])
    f.suptitle('Treg vs. Off Target Signaling Varing Beta Affinity', fontsize=10)

    return f


def cytBindingModel(counts, betaAffs, val, mut, x=False, date=False):
    """Runs binding model for a given mutein, valency, dose, and cell type."""
    #mut = mut
    doseVec = np.array([0.1])
    recCount = np.ravel(counts)

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    Affs = mutAffDF.loc[(mutAffDF.Mutein == mut)]

    Affs = np.power(np.array([Affs["IL2RaKD"].values, [betaAffs]]) / 1e9, -1)
    
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

def bindingCalc(df, targCell, offTCells, betaAffs, val, mut):
    """Calculates selectivity for no additional epitope"""
    targetBound = 0
    offTargetBound = 0

    cd25DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD25')]
    cd122DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD122')]

    for i, cd25Count in enumerate(cd25DF[targCell].item()):
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count]
        targetBound += cytBindingModel(counts, betaAffs, val, mut)

    for cellT in offTCells:
        for i, cd25Count in enumerate(cd25DF[cellT].item()):
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count]
            offTargetBound += cytBindingModel(counts,betaAffs, val, mut)

    return targetBound, offTargetBound


def cytBindingModel_bispec(counts, betaAffs, recXaff, val, x=False):
    """Runs binding model for a given mutein, valency, dose, and cell type."""

    mut = 'R38Q/H16N' #Try switching this to R38Q
    doseVec = np.array([0.1])

    recXaff = np.power(10, recXaff)

    recCount = np.ravel(counts)

    mutAffDF = pd.read_csv(join(path_here, "data/WTmutAffData.csv"))
    Affs = mutAffDF.loc[(mutAffDF.Mutein == mut)]
    Affs = np.power(np.array([Affs["IL2RaKD"].values, [betaAffs]]) / 1e9, -1)
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


def bindingCalc_bispec(df, targCell, offTCells,betaAffs,val):
    """Calculates selectivity for no additional epitope"""
    targetBound = 0
    offTargetBound = 0

    cd25DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD25')]
    cd122DF = df.loc[(df.Type == 'Standard') & (df.Epitope == 'CD122')]

    for i, cd25Count in enumerate(cd25DF[targCell].item()):
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count,cd25Count]
        targetBound += cytBindingModel_bispec(counts, betaAffs, 9.0, val)

    for cellT in offTCells:
        for i, cd25Count in enumerate(cd25DF[cellT].item()):
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count,cd25Count]
            offTargetBound += cytBindingModel_bispec(counts, betaAffs, 9.0, val)

    return targetBound, offTargetBound

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
        targetBound += cytBindingModel_bispec(counts, recXaff)
    for cellT in offTCells:
        for i, epCount in enumerate(epitopeDF[cellT].item()):
            cd25Count = cd25DF[cellT].item()[i]
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count, epCount]

            offTargetBound += cytBindingModel_bispec(counts, recXaff)

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
