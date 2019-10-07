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
from ..flow import importF, pcaAll, pcaPlt, appPCA, fitPCA, sampleT, sampleTcolor, sampleNK, sampleNKcolor, count_data, rawData, cellCount, pcaPltColor, pcaAllCellType, loadingPlot
from pathlib import Path
from matplotlib import pyplot as plt
from FlowCytometryTools import FCMeasurement
from sklearn import preprocessing
from sklearn.decomposition import PCA

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((13, 10), (3, 4))
    ax[10].axis('off')
    ax[11].axis('off')

    for i, item in enumerate(ax):
        if i < 10:
            subplotLabel(item, string.ascii_uppercase[i])

    Tcells = True
    plateRow  = "C"
    dose_ind = np.array([0., 6., 11.])
    Tplate = "/home/brianoj/Tplate418"
    Tsample, _ = importF(Tplate, plateRow)
    _, pstat_arrayT, xf_arrayT, loadingT = pcaAll(Tsample, Tcells=True) #take out titles req
    dataT, _, _ = sampleT(Tsample[0])
    PCAobjT, _ = fitPCA(dataT, Tcells=True)
    
    Nkplate = "/home/brianoj/Nkplate418"
    Nksample, _ = importF(Nkplate, plateRow)
    _, pstat_arrayNk, xf_arrayNk, loadingNk = pcaAll(Nksample, Tcells=False) #take out titles req
    dataNk, _, _ = sampleNK(Nksample[0])
    PCAobjNk, _ = fitPCA(dataNk, Tcells=False)
    
    for i, col in enumerate (dose_ind):
        axT, axNk = ax[i], ax[i + 4]
        col = int(col)
        
        dataT, _, _ = sampleT(Tsample[col])
        xfT = appPCA(dataT, PCAobjT, Tcells=True)
        pcaPlt(xfT, pstat_arrayT[col], axT, Tcells=True) 
        
        dataNk, _, _ = sampleNK(Nksample[col])
        xfNk = appPCA(dataNk, PCAobjNk, Tcells=False)
        pcaPlt(xfNk, pstat_arrayNk[col], axNk, Tcells=False)

    loadingPlot(loadingT, ax = ax[3], Tcells=True)
    loadingPlot(loadingNk, ax = ax[7], Tcells=False)

    _, _, _, _, colormatT = pcaAllCellType(Tsample, Tcells=True)
    _, _, _, _, colormatNK = pcaAllCellType(Nksample, Tcells=False) #take out titles req

    pcaPltColor(xf_arrayT[4], colormatT[4], ax=ax[8], Tcells=True)
    pcaPltColor(xf_arrayNk[4], colormatNK[4], ax=ax[9], Tcells=False)


    return f


