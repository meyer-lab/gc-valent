"""
This creates Figure 1 for Single Cell FC datavAnalysis. Examples of PCA plots and gating and .
"""

import string
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import least_squares
from .figureCommon import subplotLabel, getSetup
import synapseclient
from ..flow import importF, PCAall, pcaPlt, appPCA, fitPCA, sampleT, sampleTcolor, sampleNK, sampleNkcolor count_data, rawData, cellCount, pcaPltColor, pcaAllCellType, loadingPlot
from pathlib import Path
from matplotlib import pyplot as plt
from FlowCytometryTools import FCMeasurement
from sklearn import preprocessing
from sklearn.decomposition import PCA

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((10, 7.5), (3, 4))
    ax[10].axis('off')
    ax[11].axis('off')

    for i, item in enumerate(ax):
        if ii < 9:
            subplotLabel(item, string.ascii_uppercase[ii])

    Tcells = True
    plateRow  = "C"
    dose_ind = np.array([0., 6., 11.])
    Tplate = "/home/brianoj/XXXXXX"
    Tsample, _ = importF(Tplate, PlateRow)
    _, pstat_arrayT, xf_arrayT, loadingT = pcaAll(Tsample, Tcells=True) #take out titles req

    Nkplate = "/home/brianoj/XXXXXX"
    Nksample, _ = importF(Nkplate, PlateRow)
    _, pstat_arrayNk, xf_arrayNk, loadingNk = pcaAll(Nksample, Tcells=False) #take out titles req
    
    for i, col in enumerate (dose_ind)
        axT, axNk = i, i + 4
        pcaPlt(xf_arrayT[col], pstat_arrayT[col], axT, tplate=True) #put ax in here # take out title
        pcaPlt(xf_arrayNk[col], pstat_arrayT[col], axNk, tplate=False)

    loadingPlot(loadingT, ax = 3, Tcells=True)
    loadingPlot(loadingNk, ax = 7, Tcells=False)

    _, _, _, _, colormatT = pcaAllCellType(Tsample, Tcells=True)
    _, _, _, _, colormatNK = pcaAllCellType(Nksample, Tcells=False) #take out titles req

    pcaPltColor(xfarrayT[4], colormatT[4], ax=8 Tcells=True)
    pcaPltColor(xfarrayNk[4], colormatNk[4], ax=9, Tcells=False)

"""
    syn = synapseclient.Synapse()
    syn.login('borcuttjahns','XXXXXXXX')
    syn20506204 = syn.get(entity='syn20506204')
    directory = "/home/brianoj/type-I-ckine-model/ckine/data/XXXXXX"

    with zipfile.ZipFile(filepath, 'r') as zip_ref:
    zip_ref.extractall(directory)
    # Get the path to the local copy of the data file
    filepath = syn20506204.path 
"""

    return f


