"""
This creates Figure 5 for Single Cell data analysis. Plots of mean, variance, and skew by cell type.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from .figureCommon import subplotLabel, getSetup
from ..imports import channels
from ..flow import importF, bead_regression
from ..FCimports import import_gates, apply_gates

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((10, 5), (3, 1))

    subplotLabel(ax)

    # import bead data and run regression to get equations
    lsq_cd25, lsq_cd122, lsq_cd132 = run_regression()

    # create dataframe with gated samples (all replicates)
    df_gates = import_gates()
    df_signal = apply_gates("4-23", "1", df_gates)
    df_signal = df_signal.append(apply_gates("4-23", "2", df_gates))
    df_signal = df_signal.append(apply_gates("4-26", "1", df_gates))
    df_signal = df_signal.append(apply_gates("4-26", "2", df_gates))
    print(df_signal)

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
                    data = df_signal.loc[(df_signal["Cell Type"] == cell) & (df_signal["Date"] == date) & (df_signal["Plate"] == plate)][channels_[j]]
                    rec_counts = np.zeros(len(data))
                    for k, signal in enumerate(data):
                        A, B, C, D = lsq_params[j]
                        rec_counts[k] = C * (((A - D) / (signal - D)) - 1)**(1 / B)
                    df_add = pd.DataFrame({"Cell Type": np.tile(cell, len(data)), "Receptor": np.tile(receptor, len(data)),
                                           "Count": rec_counts, "Date": np.tile(date, len(data)), "Plate": np.tile(plate, len(data))})
                    df_rec = df_rec.append(df_add)

    # write to csv
    update_path = path_here + "/data/receptor_levels.csv"
    df_rec.to_csv(str(update_path), index=False, header=True)

    # calculate mean, variance, and skew for each replicate
    df_stats = pd.DataFrame(columns=["Cell Type", "Receptor", "Mean", "Variance", "Skew", "Date", "Plate"])
    for _, cell in enumerate(cell_names):
        for j, receptor in enumerate(receptors_):
            for _, date in enumerate(dates):
                for _, plate in enumerate(plates):
                    df_subset = df_rec.loc[(df_rec["Cell Type"] == cell) & (df_rec["Receptor"] == receptor) & (df_rec["Date"] == date) & (df_rec["Plate"] == plate)]["Count"]
                    mean_ = df_subset.mean()
                    var_ = df_subset.var()
                    skew_ = df_subset.skew()
                    df_new = pd.DataFrame(columns=["Cell Type", "Receptor", "Mean", "Variance", "Skew", "Date", "Plate"])
                    df_new.loc[0] = [cell, receptor, mean_, var_, skew_, date, plate]
                    df_stats = df_stats.append(df_new)

    # plot log10 of mean, variance, and skew
    df_stats["Mean"] = np.log10(df_stats["Mean"].astype(np.float))
    df_stats["Variance"] = np.log10(df_stats["Variance"].astype(np.float))
    df_stats["Skew"] = np.log10(df_stats["Skew"].astype(np.float))
    celltype_pointplot(ax[0], df_stats, "Mean")
    celltype_pointplot(ax[1], df_stats, "Variance")
    celltype_pointplot(ax[2], df_stats, "Skew")

    return f


def celltype_pointplot(ax, df, moment):
    """ Plots a given distribution moment with SD among replicates for all cell types and receptors. """
    sns.pointplot(x="Cell Type", y=moment, hue="Receptor", data=df, ci='sd', join=False, dodge=True, ax=ax)
    ax.set_ylabel("log(" + moment + ")")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)


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

def compMatrix(date, plate, panel):
    """Applies compensation matrix given parameters date in mm-dd, plate number and panel A, B, or C."""
    path = path_here + "/data/compensation/"+date+"/Plate "+plate+"/Plate "+plate+" - "+panel+".csv"
    header_names = ['Channel1', 'Channel2', 'Comp']
    df_comp = pd.read_csv(path, header=None, skiprows=1, names=header_names)
    #type-I-ckine-model/ckine/data/compensation/04-23/Plate 1/Plate 1 - A.csv
    # Add diangonal values of 100 to compensation values
    addedChannels = []
    for i in df_comp.index:
        channelName = df_comp.iloc[i]['Channel1']
        if channelName not in addedChannels:
            addedChannels.append(channelName)
            df2 = pd.DataFrame([[channelName, channelName, 100]], columns=['Channel1','Channel2','Comp'])
            df_comp = df_comp.append(df2, ignore_index=True)
        
    df_matrix = pd.DataFrame(index=addedChannels, columns=addedChannels)
    for i in df_matrix.index:
        for c in df_matrix.columns:
            df_matrix.at[i, c] = df_comp.loc[(df_comp['Channel1'] == i) & (df_comp['Channel2'] == c), 'Comp'].iloc[0]
            #switch i and c to transpose
    
    #df_matrix now has all values in square matrix form
    return df_matrix
        