"""
This creates Figure S3, single cell analysis."""

import numpy as np
from ..imports import import_pstat_all
from .figureCommon import getSetup, plotDoseResponses


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    df = import_pstat_all(singleCell=True)

    ligands = df.Ligand.unique()
    cells = df.Cell.unique()
    ligList = np.array([])
    valList = np.array([])

    for lig in ligands:
        uniqueVals = df.loc[(df.Ligand == lig)].Bivalent.unique()
        valList = np.append(valList, uniqueVals)
        ligList = np.append(ligList, np.tile(lig, uniqueVals.size))

    ax, f = getSetup((6, 25), (ligList.size, 2))

    for i, lig in enumerate(ligList):
        for j, cell in enumerate(cells):
            plotDoseResponses(ax[2 * i + j], df, lig, valList[i], cell, singleCell=True)

    return f
