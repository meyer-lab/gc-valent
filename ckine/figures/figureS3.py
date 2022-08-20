"""
This creates Figure 3, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from os.path import join
from copy import copy
from .figureCommon import subplotLabel, getSetup, getLigDict, get_cellTypeDict, CITE_RIDGE, CITE_SVM


path_here = os.path.dirname(os.path.dirname(__file__))
plt.rcParams['svg.fonttype'] = 'none'
ligDict = getLigDict()
cellDict = get_cellTypeDict()


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((7, 2.5), (1, 3))
    subplotLabel(ax)

    cellTarget = "Treg"
    CITE_RIDGE(ax[0], cellTarget, RNA=True)
    CITE_SVM(ax[1], cellTarget, sampleFrac=0.2)
    CITE_SVM(ax[2], cellTarget, sampleFrac=0.2, RNA=True)

    return f
