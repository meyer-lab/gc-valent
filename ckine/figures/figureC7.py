"""
This creates Figure 7, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
from os.path import join, dirname
import numpy as np
import pandas as pd
import tensorly as tl
from tensorly.decomposition._cp import initialize_cp
from tensorly.cp_tensor import CPTensor
from tensorly.decomposition import non_negative_parafac
from .figureCommon import subplotLabel, getSetup
from ..imports import import_pstat_all

path_here = os.path.dirname(os.path.dirname(__file__))

def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)

    # Imports receptor levels from .csv created by figC5
    respDF = import_pstat_all()
    respTensor = makeTensor(respDF)
    nnTens, maskTens = getMaskTens()
    R2Xplot(ax[0], respTensor, 5)
    
    return f


def makeTensor(sigDF):
    """Makes tensor of data with dimensions mutein x valency x time point x concentration x cell type"""
    ligands = sigDF.Ligand.unique()
    tps = sigDF.Time.unique()
    concs = sigDF.Dose.unique()
    cellTypes = sigDF.Cell.unique()
    tensor = np.empty((len(ligands), len(tps), len(concs), len(cellTypes)))
    tensor[:] = np.nan
    for i, lig in enumerate(ligands):
        for j, tp in enumerate(tps):
            for k, conc in enumerate(concs):
                for ii, cell in enumerate(cellTypes):
                    entry = sigDF.loc[(sigDF.Ligand == lig) & (sigDF.Time == tp) & (sigDF.Dose == conc) & (sigDF.Cell == cell)].Mean.values
                    if len(entry) >= 1:
                        tensor[i, j, k, ii] = np.mean(entry)

    return tensor


def R2Xplot(ax, statTens, compNum):
    """Creates R2X plot for non-neg CP tensor decomposition"""
    varHold = np.zeros(compNum)
    for i in range(1, compNum+1):
        tFac = non_negative_parafac(np.nan_to_num(statTens, nan=0), rank=i, n_iter_max=1000)
        varHold[i-1] = calcR2X(statTens, tFac)

    ax.scatter(np.arange(1, compNum+1), varHold, c='k', s=20.)
    ax.set(ylabel="Variance Explained", xlabel="Number of Components", ylim=(0, 1), xlim=(0, compNum+1), xticks=np.arange(0, compNum+1))


def calcR2X(tensorIn, tensorFac):
    """ Calculate R2X. """
    tErr = np.nanvar(tl.cp_to_tensor(tensorFac) - tensorIn)
    return 1.0 - (tErr) / (np.nanvar(tensorIn))


def getMaskTens(tensor):
    """Returns binary mask of tensor marking nan locations, and a tensor copy with NaNs as zeros"""
    masktensor = tensor.copy()
    tensorNoNan = tensor.copy()

    masktensor[np.isnan(masktensor)] = 0
    masktensor[np.invert(np.isnan(masktensor))] = 1
    tensorNoNan[np.isnan(masktensor)] = 0

    return tensorNoNan, masktensor
