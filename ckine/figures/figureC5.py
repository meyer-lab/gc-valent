"""
This creates Figure 5 for Single Cell data analysis. Plots of mean, variance, and skew by cell type.
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import fsolve
from FlowCytometryTools import ThresholdGate, PolyGate
from .figureCommon import subplotLabel, getSetup
from ..imports import channels, receptors
from ..flow import importF, bead_regression, exp_dec

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((10, 5), (3, 1))

    subplotLabel(ax)

    print(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads/")

    # import bead data and run regression to get equations
    sampleD, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads", "D")
    sampleE, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads/", "E")
    sampleF, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads/", "F")
    sampleG, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads/", "G")
    sampleH, _ = importF(path_here + "/data/flow/2019-04-23 Receptor Quant - Beads/", "H")

    recQuant1 = np.array([0., 4407, 59840, 179953, 625180])  # CD25, CD122, IL15
    recQuant2 = np.array([0., 7311, 44263, 161876, 269561])  # CD132

    _, lsq_cd25 = bead_regression(sampleD, channels['D'], receptors['D'], recQuant1)
    _, lsq_cd122 = bead_regression(sampleE, channels['E'], receptors['E'], recQuant2, 2, True)
    _, lsq_cd132 = bead_regression(sampleF, channels['F'], receptors['F'], recQuant1)
    _, lsq_il15_1 = bead_regression(sampleG, channels['G'], receptors['G'], recQuant1)
    _, lsq_il15_2 = bead_regression(sampleH, channels['H'], receptors['H'], recQuant1)
    
    # import data
    sample1A, _ = importF(path_here + "/data/PBMC-rec-quant/04-23/Plate 1/Plate 1 - Panel 1 IL2R/", "A")
    sample2B, _ = importF(path_here + "/data/PBMC-rec-quant/04-23/Plate 1/Plate 1 - Panel 2 IL2R/", "B")
    sample3C, _ = importF(path_here + "/data/PBMC-rec-quant/04-23/Plate 1/Plate 1 - Panel 3 IL2R/", "C")
    
    panel1 = combineWells(sample1A, ['VL6-H', 'VL4-H', 'BL1-H', 'VL1-H', 'BL3-H', 'BL5-H', 'RL1-H'])
    panel2 = combineWells(sample2B, ['VL4-H', 'BL3-H', 'FSC-H', 'VL1-H', 'BL5-H', 'RL1-H'])
    panel3 = combineWells(sample3C, ['VL6-H', 'VL4-H', 'BL3-H', 'VL1-H', 'BL5-H', 'RL1-H'])
    
    # gate data to sort cell types
    df_signal = implement_gating(panel1, panel2, panel3)
    
    # make new dataframe for receptor counts
    df_rec = pd.DataFrame(columns=["Cell Type", "Receptor", "Count"])
    cell_names = ["T-reg", "Naive Treg", "Mem Treg", "T-helper", "Naive Th", "Mem Th", "NK", "NKT", "CD8+", "Naive CD8+", "Mem CD8+"]
    receptors_ = ["CD25", "CD122", "CD132", "IL15(1)", "IL15(2)"]
    channels_ = ["VL1-H", "BL5-H", "RL1-H", "RL1-H", "RL1-H"]
    lsq_params = [lsq_cd25, lsq_cd122, lsq_cd132, lsq_il15_1, lsq_il15_2]
    
    for i, cell in enumerate(cell_names):
        for j, receptor in enumerate(receptors_):
            data = df_signal.loc[(df_signal["Cell Type"] == cell)][channels_[j]]
            print(data)
            rec_counts = np.zeros(len(data))
            for k, signal in enumerate(data):
                print(lsq_params[j])
                print(signal)
                A, B, C, D = lsq_params[j]
                rec_counts[k] = C * (((A-D)/(signal-D))-1)**(1/B)
            df_add = pd.DataFrame({"Cell Type": np.tile(cell, len(data)), "Receptor": np.tile(receptor, len(data)), "Count": rec_counts})
            df_rec = df_rec.append(df_add)
    
    # write to csv
    print(df_rec)
    update_path = path_here + "/data/receptor_levels.csv"
    print(update_path)
    df_rec.to_csv(r'update_path', index = False, header=True)
    
    df_stats = pd.DataFrame(columns=["Cell Type", "Receptor", "Mean", "Variance", "Skew"])
    for i, cell in enumerate(cell_names):
        for j, receptor in enumerate(receptors_):
            df_subset = df_rec.loc[(df_rec["Cell Type"] == cell) & (df_rec["Receptor"] == receptor)]["Count"]
            mean_ = df_subset.mean()
            var_ = df_subset.var()
            skew_ = df_subset.skew()
            print(mean_, var_, skew_)
            df_new = pd.DataFrame(columns=["Cell Type", "Receptor", "Mean", "Variance", "Skew"])
            df_new.loc[0] = [cell, receptor, mean_, var_,  skew_]
            print(df_new)
            df_stats = df_stats.append(df_new)
            
    print(df_stats)
    
    # mean
    sns.pointplot(x="Cell Type", y="Mean", hue="Receptor", data=df_stats, join=False, dodge=True, ax=ax[0])
    ax[0].set_ylabel("Mean")
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)
    
    # variance
    sns.pointplot(x="Cell Type", y="Variance", hue="Receptor", data=df_stats, join=False, dodge=True, ax=ax[1])
    ax[1].set_ylabel("Variance")
    ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)
    
    # skew
    sns.pointplot(x="Cell Type", y="Skew", hue="Receptor", data=df_stats, join=False, dodge=True, ax=ax[2])
    ax[2].set_ylabel("Skew")
    ax[2].set_xticklabels(ax[2].get_xticklabels(), rotation=25, rotation_mode="anchor", ha="right", position=(0, 0.02), fontsize=7.5)
    

    return f


def combineWells(samples, channels_):
    """Accepts sample array returned from importF, and array of channels, returns transformed combined well data"""
    combinedSamples = samples[0]
    for sample in samples[1:]:
        combinedSamples.data = combinedSamples.data.append(sample.data)
    t_combinedSamples = combinedSamples.transform('tlog', channels=channels_)
    return t_combinedSamples

def implement_gating(panel_1, panel_2, panel_3):
    """ Gates all cell types and returns dataframe with relevant channel data. """
    df = pd.DataFrame(columns=["Cell Type", "VL1-H", "BL5-H", "RL1-H"])
    
    cd3cd4_gate = PolyGate([(5.6e03, 5.5e03), (7.9e03, 6.0e03), (7.9e03, 7.2e03), (5.5e03, 7.2e03)], ('VL4-H', 'VL6-H'), region='in', name='cd3cd4')
    samplecd3cd4 = panel_1.gate(cd3cd4_gate)
    thelp_gate = PolyGate([(3.1e03, 4.5e03), (3.5e03, 4.5e03), (6.4e03, 6.2e03), (6.4e03, 6.8e03), (3.1e03, 6.8e03)], ('VL1-H', 'BL1-H'), region='in', name='thelp')
    treg_gate = PolyGate([(4.9e03, 4e03), (5.6e03, 5.6e03), (6.8e03, 6.3e03), (8e03, 3.8e03), (8e03, 2.0e3), (4.9e03, 2.0e3)], ('VL1-H', 'BL1-H'), region='in', name='treg')
    samplethelp = samplecd3cd4.gate(thelp_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("T-helper", samplethelp.counts), "VL1-H": samplethelp.data[['VL1-H']].values.reshape((samplethelp.counts,)), "BL5-H": samplethelp.data[['BL5-H']].values.reshape((samplethelp.counts,)), "RL1-H": samplethelp.data[['RL1-H']].values.reshape((samplethelp.counts,))})
    df = df.append(df_add)
    sampletreg = samplecd3cd4.gate(treg_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("T-reg", sampletreg.counts), "VL1-H": sampletreg.data[['VL1-H']].values.reshape((sampletreg.counts,)), "BL5-H": sampletreg.data[['BL5-H']].values.reshape((sampletreg.counts,)), "RL1-H": sampletreg.data[['RL1-H']].values.reshape((sampletreg.counts,))})
    df = df.append(df_add)
    memthelp_gate = ThresholdGate(6750.0, ['BL3-H'], region='below') & ThresholdGate(3800.0, ['BL3-H'], region='above')
    naivethelp_gate = ThresholdGate(7250.0, ['BL3-H'], region='above') & ThresholdGate(8700.0, ['BL3-H'], region='below')
    samplememthelp = samplethelp.gate(memthelp_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("Mem Th", samplememthelp.counts), "VL1-H": samplememthelp.data[['VL1-H']].values.reshape((samplememthelp.counts,)), "BL5-H": samplememthelp.data[['BL5-H']].values.reshape((samplememthelp.counts,)), "RL1-H": samplememthelp.data[['RL1-H']].values.reshape((samplememthelp.counts,))})
    df = df.append(df_add)
    samplenaivethelp = samplethelp.gate(naivethelp_gate)  
    df_add = pd.DataFrame({"Cell Type": np.tile("Naive Th", samplenaivethelp.counts), "VL1-H": samplenaivethelp.data[['VL1-H']].values.reshape((samplenaivethelp.counts,)), "BL5-H": samplenaivethelp.data[['BL5-H']].values.reshape((samplenaivethelp.counts,)), "RL1-H": samplenaivethelp.data[['RL1-H']].values.reshape((samplenaivethelp.counts,))})
    df = df.append(df_add)
    memtreg_gate = ThresholdGate(6050.0, ['BL3-H'], region='below') & ThresholdGate(3600.0, ['BL3-H'], region='above')
    naivetreg_gate = ThresholdGate(6700.0, ['BL3-H'], region='above') & ThresholdGate(8600.0, ['BL3-H'], region='below')
    samplememtreg = samplethelp.gate(memtreg_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("Mem Treg", samplememtreg.counts), "VL1-H": samplememtreg.data[['VL1-H']].values.reshape((samplememtreg.counts,)), "BL5-H": samplememtreg.data[['BL5-H']].values.reshape((samplememtreg.counts,)), "RL1-H": samplememtreg.data[['RL1-H']].values.reshape((samplememtreg.counts,))})
    df = df.append(df_add)
    samplenaivetreg = samplethelp.gate(naivetreg_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("Naive Treg", samplenaivetreg.counts), "VL1-H": samplenaivetreg.data[['VL1-H']].values.reshape((samplenaivetreg.counts,)), "BL5-H": samplenaivetreg.data[['BL5-H']].values.reshape((samplenaivetreg.counts,)), "RL1-H": samplenaivetreg.data[['RL1-H']].values.reshape((samplenaivetreg.counts,))})
    df = df.append(df_add)
    
    nk_gate = PolyGate([(3800, 5700), (6700, 8200), (7200, 8100), (5200, 5300)], ('VL4-H', 'BL3-H'), region='in', name='nk')
    nkt_gate = PolyGate([(5500, 5500), (6000, 5000), (8000, 6500), (6900, 7400)], ('VL4-H', 'BL3-H'), region='in', name='nkt')
    samplenk = panel_2.gate(nk_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("NK", samplenk.counts), "VL1-H": samplenk.data[['VL1-H']].values.reshape((samplenk.counts,)), "BL5-H": samplenk.data[['BL5-H']].values.reshape((samplenk.counts,)), "RL1-H": samplenk.data[['RL1-H']].values.reshape((samplenk.counts,))})
    df = df.append(df_add)
    samplenkt = panel_2.gate(nkt_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("NKT", samplenkt.counts), "VL1-H": samplenkt.data[['VL1-H']].values.reshape((samplenkt.counts,)), "BL5-H": samplenkt.data[['BL5-H']].values.reshape((samplenkt.counts,)), "RL1-H": samplenkt.data[['RL1-H']].values.reshape((samplenkt.counts,))})
    df = df.append(df_add)
    
    cd8_gate = PolyGate([(5700, 5600), (5600, 7800), (7900, 7800), (7800, 6600)], ('VL4-H', 'VL6-H'), region='in', name='cd8')
    samplecd3cd8 = panel_3.gate(cd8_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("CD8+", samplecd3cd8.counts), "VL1-H": samplecd3cd8.data[['VL1-H']].values.reshape((samplecd3cd8.counts,)), "BL5-H": samplecd3cd8.data[['BL5-H']].values.reshape((samplecd3cd8.counts,)), "RL1-H": samplecd3cd8.data[['RL1-H']].values.reshape((samplecd3cd8.counts,))})
    df = df.append(df_add)
    memcd3cd8_gate = ThresholdGate(7100.0, ['BL3-H'], region='below') & ThresholdGate(4200.0, ['BL3-H'], region='above')
    naivecd3cd8_gate = ThresholdGate(7500.0, ['BL3-H'], region='above') & ThresholdGate(8700.0, ['BL3-H'], region='below')
    samplememcytotoxic = samplecd3cd8.gate(memcd3cd8_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("Mem CD8+", samplememcytotoxic.counts), "VL1-H": samplememcytotoxic.data[['VL1-H']].values.reshape((samplememcytotoxic.counts,)), "BL5-H": samplememcytotoxic.data[['BL5-H']].values.reshape((samplememcytotoxic.counts,)), "RL1-H": samplememcytotoxic.data[['RL1-H']].values.reshape((samplememcytotoxic.counts,))})
    df = df.append(df_add)
    samplenaivecytotoxic = samplecd3cd8.gate(naivecd3cd8_gate)
    df_add = pd.DataFrame({"Cell Type": np.tile("Naive CD8+", samplenaivecytotoxic.counts), "VL1-H": samplenaivecytotoxic.data[['VL1-H']].values.reshape((samplenaivecytotoxic.counts,)), "BL5-H": samplenaivecytotoxic.data[['BL5-H']].values.reshape((samplenaivecytotoxic.counts,)), "RL1-H": samplenaivecytotoxic.data[['RL1-H']].values.reshape((samplenaivecytotoxic.counts,))})
    df = df.append(df_add)
    
    return df
    
    