"""
This creates Figure S3, single cell analysis."""

import numpy as np
from .figureCommon import getSetup, plotDoseResponses
from ..MBmodel import runFullModel


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    modelDF = runFullModel(time=[0.5, 1], saveDict=False, singleCell=True)

    ligands = modelDF.Ligand.unique()
    cells = modelDF.Cell.unique()
    ligList = np.array([])
    valList = np.array([])

    for lig in ligands:
        uniqueVals = modelDF.loc[(modelDF.Ligand == lig)].Valency.unique()
        valList = np.append(valList, uniqueVals)
        ligList = np.append(ligList, np.tile(lig, uniqueVals.size))

    ax, f = getSetup((10, 25), (ligList.size, 2))

    for i, lig in enumerate(ligList):
        for j, cell in enumerate(cells):
            plotDoseResponses(ax[2 * i + j], modelDF, lig, cell, val=valList[i], singleCell=True)
            handles, labels = ax[2 * i + j].get_legend_handles_labels()
            ax[2 * i + j].legend([handles[0]] + handles[6::], [labels[0]] + labels[6::])

    return f
