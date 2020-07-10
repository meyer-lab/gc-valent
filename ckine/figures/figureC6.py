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

    ax, f = getSetup((10,5), (2, 4))
    subplotLabel(ax)
    
    doses = [None] * 12
    
    for i in range(0,12):
        n = 84/(3**i)
        doses[i] = n
    
    #Plots pSTAT - CD25 Correlation (plots A-D) for each exposure time
    wellRows = ['A','B','C','D']
    for x in range(1,5):
        T_rVals, T_pVals, T_title = get_pSTATcorrVals('4-18','1',wellRows[x-1],x)
        plot_pSTATCorr(ax[x-1],doses,T_rVals,T_pVals, title=T_title)
        
    path_here = os.path.dirname(os.path.dirname(__file__))
    
    #Imports receptor levels from .csv created by figC5
    receptor_levels = pd.read_csv(path_here + "/data/receptor_levels.csv", header=0)
    
    cell_types = ['T-reg','T-helper','NK','CD8+']
    
    #Calculates and plots CD25-CD122 and CD25-CD132 (plot E)
    plotCorrValBars(ax[4], receptor_levels, cell_types)
    
    
    return f

def get_pSTATcorrVals(date, plate, wellRow, panel):
    cellType = "Treg"
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

def plot_pSTATCorr(ax,wellNums,rVals,pVals, title=None):
    ax.scatter(wellNums,rVals)
    i=0
    for p in pVals:
        ax.annotate(format(p, '.2g'), (wellNums[i],rVals[i]), fontsize=6)
        i+=1
    
    ax.set_ylabel('Correlation Coefficient')
    ax.set_xlabel('Dose')
    ax.set_xscale('log')
    ax.set_title(title)
    return

def getReceptorCorrVals(R1, R2):
    """Calculates correlation between two receptors"""
    R1Count = R1['Count'].reset_index(drop=True)
    R2Count = R2['Count'].reset_index(drop=True)
    d = {'R1':R1Count,'R2':R2Count}
    rCounts = pd.DataFrame(data=d)
    rCounts = rCounts.dropna()
    #print(rCounts)
    r, p = sp.stats.pearsonr(rCounts['R1'], rCounts['R2'])
        
    return r, p
    

def plotCorrValBars(ax, receptorVals, cellTypes):
    """When given all imported receptor values, creates double bar graph with bars for each cell type"""
    df_alphaCorrelation = pd.DataFrame(columns=['Receptor','Cell Type','r-val','p-val'])
    receptors = ['CD122','CD132']
    
    #get r vals and p vals for each type and receptor pairing and put into dataframe
    for cell_type in cellTypes:
        alphaLevels = receptorVals.loc[(receptorVals['Cell Type']==cell_type) & (receptorVals['Receptor']=='CD25')]
        for receptor in receptors:
            altReceptorLevels = receptorVals.loc[(receptorVals['Cell Type']==cell_type) & (receptorVals['Receptor']==receptor)]
            rVal, pVal = getReceptorCorrVals(alphaLevels, altReceptorLevels)
            df_add = pd.DataFrame({'Receptor': receptor, 'Cell Type': cell_type, 'r-val': rVal, 'p-val': pVal}, index=[0])
            df_alphaCorrelation = df_alphaCorrelation.append(df_add, ignore_index=True)
    
    print(df_alphaCorrelation)
    
    
    #use dataframe values to create bar graph
    x_pos = np.arange(len(cellTypes))
    ax.bar(x_pos - 0.15, df_alphaCorrelation.loc[df_alphaCorrelation['Receptor']=='CD122','r-val'], width=0.3, color="darkorchid", label="CD122", tick_label=cellTypes)
    ax.bar(x_pos + 0.15, df_alphaCorrelation.loc[df_alphaCorrelation['Receptor']=='CD132','r-val'], width=0.3, color="goldenrod", label="CD132", tick_label=cellTypes)
    ax.set(ylabel=("Correlation"), ylim=(0.0, 1.0))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, fontsize=6.8, rotation_mode="anchor", ha="right")
    
    return
