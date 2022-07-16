"""
This creates Figure S5, RNA ID with SVM and RNA optimized Gene ID with Model
"""
from email.mime import base
from os.path import dirname, join
from .figureCommon import getSetup, Wass_KL_Dist, CITE_RIDGE, CITE_SVM
from ..imports import importRNACITE, importReceptors
from ..MBmodel import polyc, getKxStar
from scipy.optimize import minimize, Bounds
import pandas as pd
import seaborn as sns
import numpy as np


path_here = dirname(dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    #ax, f = getSetup((4, 6), (2, 2), multz={2: 1})
    ax, f = getSetup((3, 3), (1, 1))
    GenesDF = pd.read_csv(join(path_here, "data/RNAseq_TregUnique.csv"), index_col=0)
    GenesDF = pd.concat([GenesDF, pd.DataFrame({"Gene": ["IL2RB"]})], ignore_index=True)
    CITE_DF = importRNACITE()

    cellTarget = "Treg"
    #RIDGEdf = CITE_RIDGE(ax[0], cellTarget, RNA=True)
    #SVMdf = CITE_SVM(ax[1], cellTarget, sampleFrac=0.05, RNA=True)
    cellCount = 20000

    # Get conv factors, average them
    convFactDF = convFactCalcRNA()
    meanConv = convFactDF.Weight.mean()

    # If selectivity values are stored in CSV at end of run
    saveFile = True

    # weighting idea: take sample of everything of ~3 times and then average each types amount and use that as the size

    cellList = CITE_DF["CellType2"].unique().tolist()

    sampleSizes = []
    for cellType in cellList:
        cellSample = []
        for i in np.arange(3):
            sampleDF = CITE_DF.sample(cellCount)
            sampleSize = int(len(sampleDF.loc[sampleDF["CellType2"] == cellType]))
            cellSample.append(sampleSize)
        meanSize = np.mean(cellSample)
        sampleSizes.append(int(meanSize))

    offTCells = cellList.copy()
    offTCells.remove('Treg')

    # For each  cellType in list
    for i, cellType in enumerate(cellList):
        sampleSize = sampleSizes[i]
        cellDF = CITE_DF.loc[CITE_DF["CellType2"] == cellType].sample(sampleSize)
        cellType_abdundances = []
        # For each Gene (being done on per cell basis)
        for gene in GenesDF.Gene:
            # calculate abundance based on conversion factor
            if gene == 'IL2RA':
                convFact = convFactDF.loc[convFactDF["Receptor"] == "IL2Ra"].Weight.values
            elif gene == 'IL2RB':
                convFact = convFactDF.loc[convFactDF["Receptor"] == "IL2Rb"].Weight.values
            elif gene == "IL7R":
                convFact = convFactDF.loc[convFactDF["Receptor"] == "IL7Ra"].Weight.values
            else:
                convFact = meanConv

            citeVal = cellDF[gene].to_numpy()
            abundance = citeVal * convFact
            cellType_abdundances.append(abundance)
            # add column with this name to GenesDF and abundances list
        GenesDF[cellType] = cellType_abdundances

    # GeneDF now contains a data of single cell abundances for each cell type for each Gene
    GenesDF['Selectivity'] = -1
    # New column which will hold selectivity per Gene

    targCell = 'Treg'

    standardDF = GenesDF.loc[(GenesDF.Gene == 'IL2RA')].sample()
    standard2DF = GenesDF.loc[(GenesDF.Gene == 'IL2RB')].sample()
    standardDF = pd.concat([standardDF, standard2DF])
    standardDF['Type'] = 'Standard'

    for Gene in GenesDF['Gene'].unique():
        selectedDF = GenesDF.loc[(GenesDF.Gene == Gene)].sample()
        selectedDF['Type'] = 'Gene'
        selectedDF = pd.concat([selectedDF, standardDF])
        selectedDF.reset_index()
        # New form
        optSelectivity = 1 / (optimizeDesign(targCell, offTCells, selectedDF, Gene))
        GenesDF.loc[GenesDF['Gene'] == Gene, 'Selectivity'] = optSelectivity  # Store selectivity in DF to be used for plots

    baseSelectivity = 1 / (selecCalc(standardDF, targCell, offTCells))

    if saveFile:
        GenesDF = GenesDF[["Gene", "Selectivity"]]  # drops single cell info
        GenesDF.to_csv(join(path_here, "data/GeneSelectivityListRNA.csv"), index=False)

    # bar of each Gene
    GenesDF = GenesDF.sort_values(by=['Selectivity']).tail(10)
    xvalues = GenesDF['Gene']
    yvalues = (((GenesDF['Selectivity']) / baseSelectivity) * 100) - 100
    sns.barplot(x=xvalues, y=yvalues, color='k', ax=ax[0])
    ax[0].set(ylabel="Selectivity (% increase over WT IL2)", title="RNA")
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=45, ha="right")

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


def cytBindingModel_bispecOpt(counts, recXaff, IL2RA=False, x=False):
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
    if IL2RA:
        holder[2, 0] = holder[2, 2]
        holder[2, 2] = 1e2
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
    """Calculates selectivity for no additional Gene"""
    targetBound = 0
    offTargetBound = 0
    cd25DF = df.loc[(df.Type == 'Standard') & (df.Gene == 'IL2RA')]
    cd122DF = df.loc[(df.Type == 'Standard') & (df.Gene == 'IL2RB')]

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


def minSelecFunc(x, selectedDF, targCell, offTCells, Gene):
    """Provides the function to be minimized to get optimal selectivity"""
    targetBound = 0
    offTargetBound = 0
    recXaff = x
    GeneDF = selectedDF.loc[(selectedDF.Type == 'Gene')]
    cd25DF = selectedDF.loc[(selectedDF.Type == 'Standard') & (selectedDF.Gene == 'IL2RA')]
    cd122DF = selectedDF.loc[(selectedDF.Type == 'Standard') & (selectedDF.Gene == 'IL2RB')]
    if Gene == "IL2RA":
        IL2RA = True
    else:
        IL2RA = False

    for i, epCount in enumerate(GeneDF[targCell].item()):
        cd25Count = cd25DF[targCell].item()[i]
        cd122Count = cd122DF[targCell].item()[i]
        counts = [cd25Count, cd122Count, epCount]
        targetBound += cytBindingModel_bispecOpt(counts, recXaff, IL2RA=IL2RA)
    for cellT in offTCells:
        for i, epCount in enumerate(GeneDF[cellT].item()):
            cd25Count = cd25DF[cellT].item()[i]
            cd122Count = cd122DF[cellT].item()[i]
            counts = [cd25Count, cd122Count, epCount]

            offTargetBound += cytBindingModel_bispecOpt(counts, recXaff, IL2RA=IL2RA)

    return (offTargetBound) / (targetBound)


def optimizeDesign(targCell, offTcells, selectedDF, Gene):
    """ A more general purpose optimizer """
    if targCell == "NK":
        X0 = [6.0, 8]
    else:
        X0 = [7.0]
    optBnds = Bounds(np.full_like(X0, 6.0), np.full_like(X0, 9.0))
    optimized = minimize(minSelecFunc, X0, bounds=optBnds, args=(selectedDF, targCell, offTcells, Gene), jac="3-point")
    optSelectivity = optimized.fun

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


markDict = {"IL2RA": "IL2Ra",
            "IL2RB": "IL2Rb",
            "IL7R": "IL7Ra",
            "IL2RG": "gc"}


def convFactCalcRNA():
    """Fits a ridge classifier to the CITE data and plots those most highly correlated with T reg"""
    CITE_DF = importRNACITE()
    cellToI = ["CD4 TCM", "CD8 Naive", "NK", "CD8 TEM", "CD4 Naive", "CD4 CTL", "CD8 TCM", "Treg", "CD4 TEM"]
    markers = ["IL2RB", "IL7R", "IL2RA", "IL2RG"]
    markerDF = pd.DataFrame(columns=["Marker", "Cell Type", "Amount", "Number"])
    for marker in markers:
        for cell in cellToI:
            cellTDF = CITE_DF.loc[CITE_DF["CellType2"] == cell][marker]
            markerDF = pd.concat([markerDF, pd.DataFrame({"Marker": [marker], "Cell Type": cell, "Amount": cellTDF.mean(), "Number": cellTDF.size})])

    markerDF = markerDF.replace({"Marker": markDict, "Cell Type": cellDict})
    markerDFw = pd.DataFrame(columns=["Marker", "Cell Type", "Average"])
    for marker in markerDF.Marker.unique():
        for cell in markerDF["Cell Type"].unique():
            subDF = markerDF.loc[(markerDF["Cell Type"] == cell) & (markerDF["Marker"] == marker)]
            wAvg = np.sum(subDF.Amount.values * subDF.Number.values) / np.sum(subDF.Number.values)
            markerDFw = pd.concat([markerDFw, pd.DataFrame({"Marker": [marker], "Cell Type": cell, "Average": wAvg})])

    recDF = importReceptors()
    weightDF = pd.DataFrame(columns=["Receptor", "Weight"])

    for rec in markerDFw.Marker.unique():
        CITEval = np.array([])
        Quantval = np.array([])
        for cell in markerDF["Cell Type"].unique():
            CITEval = np.concatenate((CITEval, markerDFw.loc[(markerDFw["Cell Type"] == cell) & (markerDFw["Marker"] == rec)].Average.values))
            Quantval = np.concatenate((Quantval, recDF.loc[(recDF["Cell Type"] == cell) & (recDF["Receptor"] == rec)].Mean.values))
        weightDF = pd.concat([weightDF, pd.DataFrame({"Receptor": [rec], "Weight": np.linalg.lstsq(np.reshape(CITEval, (-1, 1)).astype(np.float), Quantval.astype(np.float), rcond=None)[0]})])

    return weightDF
