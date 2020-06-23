"""
This creates Figure 6 for IL2Ra correlatoin data analysis.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from .figureCommon import subplotLabel, getSetup
from ..imports import channels
#from ..flow import importF, bead_regression
from ..FCimports import import_gates, importF

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    
    #plate 1 is Tcells - plate 2 is NK

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)
    doses = [None] * 12
    
    for i in range(0,12):
        n = 84/(3**i)
        doses[i] = n
    
    wellRows = ['A','B','C','D']
    for x in range(1,5):
        T_rVals, T_pVals, T_title = get_p_and_r_vals('4-18','1',wellRows[x-1],x)
        plotCorr(ax[x-1],doses,T_rVals,T_pVals,T_title)
        
        NK_rVals, NK_pVals, NK_title = get_p_and_r_vals('4-18','2',wellRows[x-1],x)
        plotCorr(ax[x+3],doses,NK_rVals,NK_pVals,NK_title)
        
    
    
    return f

def get_p_and_r_vals(date, plate, wellRow, panel):
    if plate == '1':
        cellType = "Treg"
    else:
        cellType = "NK"
    times = ['4hr','2hr','1hr','30min']
    time = times[panel-1]
    
    title = cellType + " " + time
    rVals = [None]*12
    pVals = [None]*12
    
    for i in range(0,12):
        sample, _ = importF(date,plate,wellRow,str(panel),i, comp=False)
        IL2Ra_vals = sample.data['VL1-H']
        PSTAT5_vals = sample.data['RL1-H']
        
        r, p = sp.stats.pearsonr(IL2Ra_vals, PSTAT5_vals)
        rVals[i] = r
        pVals[i] = p
        
    return rVals, pVals, title
    


def plotCorr(ax,wellNums,rVals,pVals, title):
    ax.scatter(wellNums,rVals)
    i=0
    for p in pVals:
        ax.annotate(format(p, '.2g'), (wellNums[i],rVals[i]), fontsize=6)
        i+=1
    
    ax.set_ylabel('Correlation Coefficient')
    ax.set_xlabel('Dose')
    #ax.set_xticks(range(1,13))
    ax.set_xscale('log')
    ax.set_title(title)
    return