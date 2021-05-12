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

    receptor_levels = getReceptors()
    df_test = receptor_levels.loc[(receptor_levels['Cell Type'] == 'T-helper') & (receptor_levels['Receptor'] == 'CD122')]
    print(df_test)
    print("Nan:", df_test.isna().sum())

    binned_tregCounts = plotAlphaHistogram(ax[0], ax[2], 'T-reg', receptor_levels, 3)
    binned_thelperCounts = plotAlphaHistogram(ax[1], ax[3], 'T-helper', receptor_levels, 3)

    binnedCounts = binned_tregCounts.append(binned_thelperCounts, ignore_index=True)

    path = path_here + '/data/BinnedReceptorData.csv'
    binnedCounts.to_csv(path, index=False)

    return f


cellTypeDict = {"T-reg": "Treg",
                "T-helper": "Thelper"}


def plotAlphaHistogram(ax1, ax2, cell_type, receptor_levels, numBins):
    """Makes two plots: histogram of IL2Ra values, and IL2Rb means in those bins"""
    alphaLevels = receptor_levels.loc[(receptor_levels['Cell Type'] == cell_type) & (receptor_levels['Receptor'] == 'CD25')]
    betaLevels = receptor_levels.loc[(receptor_levels['Cell Type'] == cell_type) & (receptor_levels['Receptor'] == 'CD122')]

    alphaCounts = alphaLevels['Count'].reset_index(drop=True)
    betaCounts = betaLevels['Count'].reset_index(drop=True)
    print("BetaMean:",betaCounts.mean())
    d = {'alpha': alphaCounts, 'beta': betaCounts}
    recepCounts = pd.DataFrame(data=d)
    print(recepCounts)
    print("Nan:", recepCounts.isna().sum())
    #recepCounts = recepCounts.dropna()
    recepCounts = recepCounts[(recepCounts[['alpha', 'beta']] != 0).all(axis=1)]

    min_ = alphaCounts.quantile(0.05)
    max_ = alphaCounts.quantile(0.95)

    logbins = np.logspace(np.log10(min_), np.log10(max_), numBins + 1)

    alphaMeans, _, binnumber = stats.binned_statistic(alphaCounts, alphaCounts, statistic='mean', bins=logbins)

    ax1.hist(alphaCounts, bins=logbins)
    title = cell_type + " Histogram"
    ax1.set(xlabel='IL2Ra Proteins/Cell', ylabel='Number of Cells', title=title)
    ax1.set_yscale('log')
    ax1.set_xscale('log')

    data = {'Count': alphaCounts, 'Bin': binnumber}
    binnedAlphaCounts = pd.DataFrame(data)
    binnedAlphaCounts['Cell Type'] = cell_type

    alphaMeans, _, _ = stats.binned_statistic(recepCounts['alpha'], recepCounts['alpha'], statistic='mean', bins=logbins)
    betaMeans, _, _ = stats.binned_statistic(recepCounts['alpha'].values, recepCounts['beta'].values, statistic='mean', bins=logbins)

    binnedRecCounts = pd.DataFrame({"Receptor": "IL2Ra", "Bin": np.arange(1, numBins + 1), "Mean": alphaMeans, "Cell Type": cellTypeDict[cell_type]})
    binnedRecCounts = binnedRecCounts.append(pd.DataFrame({"Receptor": "IL2Rb", "Bin": np.arange(1, numBins + 1), "Mean": betaMeans, "Cell Type": cellTypeDict[cell_type]}))

    ax2.scatter(alphaMeans, betaMeans)
    title = cell_type + " alpha/beta mean correlation"
    ax2.set(xlabel='IL2Ra Mean/Bin', ylabel='IL2Rb Mean/Bin', title=title)
    ax2.set_yscale('log')
    ax2.set_xscale('log')

    return binnedRecCounts


def getReceptors():
    # import bead data and run regression to get equations
    lsq_cd25, lsq_cd122, lsq_cd132 = run_regression()

    df_gates = import_gates()
    df_signal = apply_gates("4-23", "1", df_gates, correlation="CD122")
    df_signal = df_signal.append(apply_gates("4-23", "2", df_gates, correlation="CD122"))
    df_signal = df_signal.append(apply_gates("4-26", "1", df_gates, correlation="CD122"))
    df_signal = df_signal.append(apply_gates("4-26", "2", df_gates, correlation="CD122"))

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
                    data = data[data >= 0]
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
    _, lsq_cd122 = bead_regression(sampleE, channels['E'], recQuant1, 2, True)
    _, lsq_cd132 = bead_regression(sampleF, channels['F'], recQuant2)

    return lsq_cd25, lsq_cd122, lsq_cd132
