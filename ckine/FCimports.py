"""
This file includes various methods for flow cytometry analysis of fixed cells.
"""
import os
from os.path import dirname, join
from pathlib import Path
import matplotlib.cm as cm
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from FlowCytometryTools import FCMeasurement, PolyGate, ThresholdGate

path_here = dirname(dirname(__file__))

def combineWells(samples):
    """Accepts sample array returned from importF, and array of channels, returns combined well data"""
    combinedSamples = samples[0]
    for sample in samples[1:]:
        combinedSamples.data = combinedSamples.data.append(sample.data, ignore_index=True)
    return combinedSamples


def importF(date, plate, wellRow, panel, wellNum=None):
    """
    Import FCS files. Variable input: date in format mm-dd, plate #, panel #, and well letter. Output is a list of Data File Names in FCT Format
    Title/file names are returned in the array file --> later referenced in other functions as title/titles input argument
    """
    path_ = os.path.abspath("")
    
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            print(os.path.join(root, name))
        for name in dirs:
            print(os.path.join(root, name))
    
    pathname = path_ + "/ckine/data/flow/4-23_4-26_Receptor quant/" + date + "/Plate " + plate + "/Plate " + plate + " - Panel " + str(panel) + " IL2R/"
    print(pathname)
    
    # Declare arrays and int
    file = []
    sample = []
    z = 0
    # Read in user input for file path and assign to array file
    pathlist = Path(r"" + str(pathname)).glob("**/*.fcs")

    for path in pathlist:
        print(path)
        wellID = path.name.split("_")[1]
        if wellID[0] == wellRow:
            file.append(str(path))
        else:
            unstainedWell = FCMeasurement(ID="Unstained Sample", datafile=str(path))  # Stores data from unstainedWell separately
    file.sort()
    assert file != []
    # Go through each file and assign the file contents to entry in the array sample
    for entry in file:
        sample.append(FCMeasurement(ID="Test Sample" + str(z), datafile=entry))
        z += 1
    # The array sample contains data of each file in folder (one file per entry in array)
    channels_ = []
    if panel == 1:
        channels_ = ['VL6-H', 'VL4-H', 'BL1-H', 'VL1-H', 'BL3-H', 'BL5-H', 'RL1-H']
    elif panel == 2:
        channels_ = ['VL4-H', 'BL3-H', 'FSC-H', 'VL1-H', 'BL5-H', 'RL1-H']
    elif panel == 3:
        channels_ = ['VL6-H', 'VL4-H', 'BL3-H', 'VL1-H', 'BL5-H', 'RL1-H']

    if wellNum is None:
        combinedSamples = combineWells(sample)  # Combines all files from samples
        #combinedSamples = subtract_unstained_signal(combinedSamples, channels_, unstainedWell)  # Subtracts background
        return combinedSamples.transform("hlog", channels=channels_), unstainedWell

    #tsample = subtract_unstained_signal(sample[wellNum - 1], channels, unstainedWell)
    return tsample.transform("hlog", channels=channels), unstainedWell


def subtract_unstained_signal(sample, channels, unstainedWell):
    """ Subtract mean unstained signal from all input channels for a given sample. """
    meanBackground = np.mean(unstainedWell.data['RL1-H'])  # Calculates mean unstained signal
    for _, channel in enumerate(channels):
        for i, well in enumerate(sample[channel]):
            if well < meanBackground:
                print(sample[channel].iloc[i])
                sample[channel].iloc[i] = 0.
                print(sample[channel].iloc[i])
            if well >= meanBackground:
                sample[channel].iloc[i] = sample[channel].iloc[i] - meanBackground
    return sample


def import_gates():
    """ Imports dataframe with gates for all cell types and replicates. """
    data = pd.read_csv(join(path_here, "ckine/data/fc_gates.csv"))
    data.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
    return data


def apply_gates(date, plate, gates_df, subpopulations=False):
    """ Constructs dataframe with channels relevant to receptor quantification. """
    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])
    df = df.append(thelp_sample(date, plate, gates_df, mem_naive=subpopulations))
    df = df.append(treg_sample(date, plate, gates_df, mem_naive=subpopulations))
    df = df.append(nk_nkt_sample(date, plate, gates_df, nkt=subpopulations))
    df = df.append(cd8_sample(date, plate, gates_df, mem_naive=subpopulations))
    return df


def thelp_sample(date, plate, gates_df, mem_naive=False):
    """ Returns gated T-helper sample for a given date and plate. """
    panel1, unstainedWell = importF(date, plate, "A", 1)
    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])
    
    samplecd3cd4 = panel1.gate(eval(gates_df.loc[(gates_df["Name"] == 'CD3CD4') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    samplethelp = samplecd3cd4.gate(eval(gates_df.loc[(gates_df["Name"] == 'T-helper') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    samplethelp = subtract_unstained_signal(samplethelp, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
    df_add = pd.DataFrame({"Cell Type": np.tile("T-helper", samplethelp.counts), "Date": np.tile(date, samplethelp.counts), "Plate": np.tile(plate, samplethelp.counts), "VL1-H": samplethelp.data[['VL1-H']].values.reshape((samplethelp.counts,)), "BL5-H": samplethelp.data[['BL5-H']].values.reshape((samplethelp.counts,)), "RL1-H": samplethelp.data[['RL1-H']].values.reshape((samplethelp.counts,))})
    df = df.append(df_add)
    
    # separates memory and naive populations and adds to dataframe
    if mem_naive:
        samplenaive = samplethelp.gate(eval(gates_df.loc[(gates_df["Name"] == 'Naive Th') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        samplenaive = subtract_unstained_signal(samplenaive, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
        df_add = pd.DataFrame({"Cell Type": np.tile("Naive Th", samplenaive.counts), "Date": np.tile(date, samplenaive.counts), "Plate": np.tile(plate, samplenaive.counts), "VL1-H": samplenaive.data[['VL1-H']].values.reshape((samplenaive.counts,)), "BL5-H": samplenaive.data[['BL5-H']].values.reshape((samplenaive.counts,)), "RL1-H": samplenaive.data[['RL1-H']].values.reshape((samplenaive.counts,))})
        df = df.append(df_add)
        samplemem = samplethelp.gate(eval(gates_df.loc[(gates_df["Name"] == 'Mem Th') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        samplemem = subtract_unstained_signal(samplemem, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
        df_add = pd.DataFrame({"Cell Type": np.tile("Mem Th", samplemem.counts), "Date": np.tile(date, samplemem.counts), "Plate": np.tile(plate, samplemem.counts), "VL1-H": samplemem.data[['VL1-H']].values.reshape((samplemem.counts,)), "BL5-H": samplemem.data[['BL5-H']].values.reshape((samplemem.counts,)), "RL1-H": samplemem.data[['RL1-H']].values.reshape((samplemem.counts,))})
        df = df.append(df_add)
    
    return df


def treg_sample(date, plate, gates_df, mem_naive=False):
    """ Returns gated T-reg sample for a given date and plate. """
    panel1, unstainedWell = importF(date, plate, "A", 1)
    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])
    
    samplecd3cd4 = panel1.gate(eval(gates_df.loc[(gates_df["Name"] == 'CD3CD4') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    sampletreg = samplecd3cd4.gate(eval(gates_df.loc[(gates_df["Name"] == 'T-reg') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    sampletreg = subtract_unstained_signal(sampletreg, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
    df_add = pd.DataFrame({"Cell Type": np.tile("T-reg", sampletreg.counts), "Date": np.tile(date, sampletreg.counts), "Plate": np.tile(plate, sampletreg.counts), "VL1-H": sampletreg.data[['VL1-H']].values.reshape((sampletreg.counts,)), "BL5-H": sampletreg.data[['BL5-H']].values.reshape((sampletreg.counts,)), "RL1-H": sampletreg.data[['RL1-H']].values.reshape((sampletreg.counts,))})
    df = df.append(df_add)
    
    # separates memory and naive populations and adds to dataframe
    if mem_naive:
        samplenaive = sampletreg.gate(eval(gates_df.loc[(gates_df["Name"] == 'Naive Treg') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        samplenaive = subtract_unstained_signal(samplenaive, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
        df_add = pd.DataFrame({"Cell Type": np.tile("Naive Treg", samplenaive.counts), "Date": np.tile(date, samplenaive.counts), "Plate": np.tile(plate, samplenaive.counts), "VL1-H": samplenaive.data[['VL1-H']].values.reshape((samplenaive.counts,)), "BL5-H": samplenaive.data[['BL5-H']].values.reshape((samplenaive.counts,)), "RL1-H": samplenaive.data[['RL1-H']].values.reshape((samplenaive.counts,))})
        df = df.append(df_add)
        samplemem = sampletreg.gate(eval(gates_df.loc[(gates_df["Name"] == 'Mem Th') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        samplemem = subtract_unstained_signal(samplemem, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
        df_add = pd.DataFrame({"Cell Type": np.tile("Mem Treg", samplemem.counts), "Date": np.tile(date, samplemem.counts), "Plate": np.tile(plate, samplemem.counts), "VL1-H": samplemem.data[['VL1-H']].values.reshape((samplemem.counts,)), "BL5-H": samplemem.data[['BL5-H']].values.reshape((samplemem.counts,)), "RL1-H": samplemem.data[['RL1-H']].values.reshape((samplemem.counts,))})
        df = df.append(df_add)
    
    return df


def nk_nkt_sample(date, plate, gates_df, nkt=False):
    """ Returns gated NK sample for a given date and plate. """
    panel2, unstainedWell = importF(date, plate, "B", 2)    
    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])
    
    samplenk = panel2.gate(eval(gates_df.loc[(gates_df["Name"] == 'NK') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    samplenk = subtract_unstained_signal(samplenk, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
    df_add = pd.DataFrame({"Cell Type": np.tile("NK", samplenk.counts), "Date": np.tile(date, samplenk.counts), "Plate": np.tile(plate, samplenk.counts), "VL1-H": samplenk.data[['VL1-H']].values.reshape((samplenk.counts,)), "BL5-H": samplenk.data[['BL5-H']].values.reshape((samplenk.counts,)), "RL1-H": samplenk.data[['RL1-H']].values.reshape((samplenk.counts,))})
    df = df.append(df_add)
    
    # gates NKT population and adds to dataframe
    if nkt:
        samplenkt = panel2.gate(eval(gates_df.loc[(gates_df["Name"] == 'NKT') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        samplenkt = subtract_unstained_signal(samplenkt, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
        df_add = pd.DataFrame({"Cell Type": np.tile("NKT", samplenkt.counts), "Date": np.tile(date, samplenkt.counts), "Plate": np.tile(plate, samplenkt.counts), "VL1-H": samplenkt.data[['VL1-H']].values.reshape((samplenkt.counts,)), "BL5-H": samplenkt.data[['BL5-H']].values.reshape((samplenkt.counts,)), "RL1-H": samplenkt.data[['RL1-H']].values.reshape((samplenkt.counts,))})
        df = df.append(df_add)
    
    return df


def cd8_sample(date, plate, gates_df, mem_naive=False):
    """ Returns gated CD8+ sample for a given date and plate. """
    panel3, unstainedWell = importF(date, plate, "C", 3)
    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])
    
    samplecd8 = panel3.gate(eval(gates_df.loc[(gates_df["Name"] == 'CD8+') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    samplecd8 = subtract_unstained_signal(samplecd8, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
    df_add = pd.DataFrame({"Cell Type": np.tile("CD8+", samplecd8.counts), "Date": np.tile(date, samplecd8.counts), "Plate": np.tile(plate, samplecd8.counts), "VL1-H": samplecd8.data[['VL1-H']].values.reshape((samplecd8.counts,)), "BL5-H": samplecd8.data[['BL5-H']].values.reshape((samplecd8.counts,)), "RL1-H": samplecd8.data[['RL1-H']].values.reshape((samplecd8.counts,))})
    df = df.append(df_add)
    
    # separates memory and naive populations and adds to dataframe
    if mem_naive:
        samplenaive = samplecd8.gate(eval(gates_df.loc[(gates_df["Name"] == 'Naive CD8+') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        samplenaive = subtract_unstained_signal(samplenaive, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
        df_add = pd.DataFrame({"Cell Type": np.tile("Naive CD8+", samplenaive.counts), "Date": np.tile(date, samplenaive.counts), "Plate": np.tile(plate, samplenaive.counts), "VL1-H": samplenaive.data[['VL1-H']].values.reshape((samplenaive.counts,)), "BL5-H": samplenaive.data[['BL5-H']].values.reshape((samplenaive.counts,)), "RL1-H": samplenaive.data[['RL1-H']].values.reshape((samplenaive.counts,))})
        df = df.append(df_add)
        samplemem = samplecd8.gate(eval(gates_df.loc[(gates_df["Name"] == 'Mem CD8+') & (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        samplemem = subtract_unstained_signal(samplemem, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
        df_add = pd.DataFrame({"Cell Type": np.tile("Mem CD8+", samplemem.counts), "Date": np.tile(date, samplemem.counts), "Plate": np.tile(plate, samplemem.counts), "VL1-H": samplemem.data[['VL1-H']].values.reshape((samplemem.counts,)), "BL5-H": samplemem.data[['BL5-H']].values.reshape((samplemem.counts,)), "RL1-H": samplemem.data[['RL1-H']].values.reshape((samplemem.counts,))})
        df = df.append(df_add)
    
    return df
