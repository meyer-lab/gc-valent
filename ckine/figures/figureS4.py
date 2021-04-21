"""
This creates Figure 5 for IL2Ra correlation data analysis.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from .figureCommon import subplotLabel, getSetup
from ..imports import channels
from ..flow import importF, bead_regression
from ..FCimports import import_gates, apply_gates

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)
    print(ax[0])

    
    #cell_types = ['T-reg', 'T-helper']
    receptor_levels = getReceptors()
    
    #alphaLevels_Treg = receptor_levels.loc[(receptor_levels['Cell Type'] == 'T-reg') & (receptor_levels['Receptor'] == 'CD25')]
    #alphaLevels_Thelper = receptor_levels.loc[(receptor_levels['Cell Type'] == 'T-helper') & (receptor_levels['Receptor'] == 'CD25')]
    #alphaLevels = alphaLevels_Treg.append(alphaLevels_Thelper)
    

    plotAlphaHistogram(ax[0],'T-reg',receptor_levels)
    plotAlphaHistogram(ax[1],'T-helper', receptor_levels)
        
        #fullData = fullData.append(alphaCounts)

    #fullData.to_csv(r'ckine/output/S4data.csv', index = False)

    return f

def plotAlphaHistogram(ax, cell_type, receptor_levels):
    alphaLevels = receptor_levels.loc[(receptor_levels['Cell Type'] == cell_type) & (receptor_levels['Receptor'] == 'CD25')]
        
    alphaCounts = alphaLevels['Count'].reset_index(drop=True)
    
    print(alphaCounts.describe())
    min_ = alphaCounts.min()
    max_ = alphaCounts.max()

    logbins = np.logspace(np.log10(min_),np.log10(max_),8)
    print(logbins)
    
    #print(alphaCounts.values)
    medians, _ , _ = stats.binned_statistic(alphaCounts.values, alphaCounts.values, statistic ='median', bins = logbins)
    print(medians)

    
    ax.hist(alphaCounts, bins=logbins)
    title = cell_type + " Histogram"
    ax.set(xlabel='IL2Ra Proteins/Cell', ylabel='Number of Cells',title=title)
    ax.set_yscale('log')
    ax.set_xscale('log')
    #plt.xticks(medians, medians) 


def getReceptors():
    # import bead data and run regression to get equations
    lsq_cd25, lsq_cd122, lsq_cd132 = run_regression()

    # create dataframe with gated samples (all replicates)
    df_gates = import_gates()
    df_signal = apply_gates("4-23", "1", df_gates)
    df_signal = df_signal.append(apply_gates("4-23", "2", df_gates))
    df_signal = df_signal.append(apply_gates("4-26", "1", df_gates))
    df_signal = df_signal.append(apply_gates("4-26", "2", df_gates))

    # make new dataframe for receptor counts
    df_rec = pd.DataFrame(columns=["Cell Type", "Receptor", "Count", "Date", "Plate"])
    cell_names = ["T-reg", "T-helper", "NK", "CD8+"]
    receptors_ = ["CD25", "CD122", "CD132"]
    channels_ = ["VL1-H", "BL5-H", "RL1-H"]
    lsq_params = [lsq_cd25, lsq_cd122, lsq_cd132]
    dates = ["4-23", "4-26"]
    plates = ["1", "2"]

    # calculate receptor counts
    for _, cell in enumerate(cell_names):
        for j, receptor in enumerate(receptors_):
            for _, date in enumerate(dates):
                for _, plate in enumerate(plates):
                    data = df_signal.loc[(df_signal["Cell Type"] == cell) & (df_signal["Receptor"] == receptor) & (df_signal["Date"] == date) & (df_signal["Plate"] == plate)][channels_[j]]
                    rec_counts = np.zeros(len(data))
                    for k, signal in enumerate(data):
                        A, B, C, D = lsq_params[j]
                        rec_counts[k] = C * (((A - D) / (signal - D)) - 1)**(1 / B)
                    df_add = pd.DataFrame({"Cell Type": np.tile(cell, len(data)), "Receptor": np.tile(receptor, len(data)),
                                           "Count": rec_counts, "Date": np.tile(date, len(data)), "Plate": np.tile(plate, len(data))})
                    df_rec = df_rec.append(df_add)

    return df_rec


def run_regression():
    """ Imports bead data and runs regression to get least squares parameters for conversion of signal to receptor count. """
    sampleD, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads", "D")
    sampleE, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads/", "E")
    sampleF, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads/", "F")

    recQuant1 = np.array([0., 4407, 59840, 179953, 625180])  # CD25, CD122
    recQuant2 = np.array([0., 7311, 44263, 161876, 269561])  # CD132

    _, lsq_cd25 = bead_regression(sampleD, channels['D'], recQuant1)
    _, lsq_cd122 = bead_regression(sampleE, channels['E'], recQuant2, 2, True)
    _, lsq_cd132 = bead_regression(sampleF, channels['F'], recQuant1)

    return lsq_cd25, lsq_cd122, lsq_cd132
