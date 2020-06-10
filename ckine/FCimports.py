"""
This file includes various methods for flow cytometry analysis of fixed cells.
"""
import ast
import os
from os.path import dirname, join
from pathlib import Path
import pandas as pd
import numpy as np
from FlowCytometryTools import FCMeasurement

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

    pathname = path_ + "/ckine/data/flow/" + date + " Live PBMC Receptor Data/Plate " + plate + "/Plate " + plate + " - Panel " + str(panel) + " IL2R/"

    # Declare arrays and int
    file = []
    sample = []
    z = 0
    # Read in user input for file path and assign to array file
    pathlist = Path(r"" + str(pathname)).glob("**/*.fcs")

    for path in pathlist:
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

    if wellNum is None:
        combinedSamples = combineWells(sample)  # Combines all files from samples
        compSample = applyMatrix(combinedSamples, compMatrix(date, plate, wellRow))  # Applies compensation matrix
        return compSample, unstainedWell

    compSample = applyMatrix(sample, compMatrix(date, plate, wellRow))
    return compSample, unstainedWell


def subtract_unstained_signal(sample, channels, unstainedWell):
    """ Subtract mean unstained signal from all input channels for a given sample. """
    for _, channel in enumerate(channels):
        meanBackground = np.mean(unstainedWell.data[channel])  # Calculates mean unstained signal for given channel
        sample[channel] = np.maximum(sample[channel] - meanBackground, 0.0)

    return sample


def compMatrix(date, plate, panel, invert=True):
    """Creates compensation matrix given parameters date in mm-dd, plate number and panel A, B, or C."""
    path = path_here + "/ckine/data/compensation/0" + date + "/Plate " + plate + "/Plate " + plate + " - " + panel + ".csv"
    # imports csv file with comp values as a dataframe
    header_names = ['Channel1', 'Channel2', 'Comp']
    df_comp = pd.read_csv(path, header=None, skiprows=1, names=header_names)
    # Add diangonal values of 100 to compensation values
    addedChannels = []
    for i in df_comp.index:
        channelName = df_comp.iloc[i]['Channel1']
        if channelName not in addedChannels:  # Ensures a diagonal value is only added once for each channel
            addedChannels.append(channelName)
            df2 = pd.DataFrame([[channelName, channelName, 100]], columns=['Channel1', 'Channel2', 'Comp'])  # Creates new row for dataframe
            df_comp = df_comp.append(df2, ignore_index=True)  # Adds row
    # create square matrix from compensation values
    df_matrix = pd.DataFrame(index=addedChannels, columns=addedChannels)  # df_matrix is now a square and has exactly one row and one column for each channel
    for i in df_matrix.index:
        for c in df_matrix.columns:
            df_matrix.at[i, c] = df_comp.loc[(df_comp['Channel1'] == c) & (df_comp['Channel2'] == i), 'Comp'].iloc[0]  # Fills in square matrix by finding corresponding comp value from csv
    # df_matrix now has all values in square matrix form
    df_matrix = df_matrix.div(100)
    if invert:  # true by default, inverts matrix before returning it
        a = np.matrix(df_matrix.values, dtype=float)  # Convert to np to allow for linalg usage
        df_matrix = pd.DataFrame(np.linalg.inv(a), df_matrix.columns, df_matrix.index)  # Calculate inverse and put pack as dataframe
    return df_matrix


def applyMatrix(sample, matrix):
    """Multiples two matrices together in the order sample dot matrix"""
    holder = pd.DataFrame()  # Will hold columns not being compensated
    for c in sample.data.columns:
        if c not in matrix:  # If sample channel column is not found in matrix
            holder = holder.join(sample.data[[c]], how='right')  # Store for after calculation
            sample.data = sample.data.drop([c], axis=1)  # Removed column to allow for matrix multiplication

    cols = sample.data.columns
    matrix = matrix[cols]
    matrix = matrix.reindex(cols)
    sample.data = sample.data.dot(matrix)  # Use matrix multiplication to compensate the relevant data
    sample.data = sample.data.join(holder)  # Restore uncompensated channels to sample
    return sample


def import_gates():
    """ Imports dataframe with gates for all cell types and replicates. """
    data = pd.read_csv(join(path_here, "ckine/data/fc_gates.csv"))
    data.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
    return data


def apply_gates(date, plate, gates_df, subpopulations=False):
    """ Constructs dataframe with channels relevant to receptor quantification. """
    df, unstainedWell = thelp_sample(date, plate, gates_df, mem_naive=subpopulations)
    df = df.append(treg_sample(date, plate, gates_df, mem_naive=subpopulations))
    df = df.append(nk_nkt_sample(date, plate, gates_df, nkt=subpopulations))
    df = df.append(cd8_sample(date, plate, gates_df, mem_naive=subpopulations))
    # All samples for data and plate processed combined
    df = subtract_unstained_signal(df, ["VL1-H", "BL5-H", "RL1-H"], unstainedWell)
    # Background signal substracted
    return df


def thelp_sample(date, plate, gates_df, mem_naive=False):
    """ Returns gated T-helper sample for a given date and plate. """
    # import data and create transformed df for gating
    panel1, unstainedWell = importF(date, plate, "A", 1)
    panel1_t = panel1.transform("tlog", channels=['VL6-H', 'VL4-H', 'BL1-H', 'VL1-H', 'BL3-H'])  # Creates copy of panel1 to transform and gate

    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])  # initialize dataframe for receptor quant channels

    # implement gating, revert tlog, and add to dataframe
    samplecd3cd4 = panel1_t.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'CD3CD4') &
                                                               (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    samplethelp = samplecd3cd4.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'T-helper') &
                                                                  (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    # Gated signals based on gating values from csv
    gated_idx = np.array(samplethelp.data.index)
    panel1.set_data(panel1.data.loc[gated_idx])  # Selects only the corresponding data points from panel1(untransformed) based on gated points from panel1_t

    df_add = pd.DataFrame({"Cell Type": np.tile("T-helper", panel1.counts), "Date": np.tile(date, panel1.counts), "Plate": np.tile(plate, panel1.counts),
                           "VL1-H": panel1.data[['VL1-H']].values.reshape((panel1.counts,)), "BL5-H": panel1.data[['BL5-H']].values.reshape((panel1.counts,)),
                           "RL1-H": panel1.data[['RL1-H']].values.reshape((panel1.counts,))})
    df = df.append(df_add)

    # separates memory and naive populations and adds to dataframe
    if mem_naive:
        panel1_n = panel1.copy()
        samplenaive = samplethelp.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'Naive Th') &
                                                                     (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        gated_idx = np.array(samplenaive.data.index)
        panel1_n.set_data(panel1.data.loc[gated_idx])
        df_add = pd.DataFrame({"Cell Type": np.tile("Naive Th", samplenaive.counts), "Date": np.tile(date, samplenaive.counts), "Plate": np.tile(plate, samplenaive.counts),
                               "VL1-H": panel1_n.data[['VL1-H']].values.reshape((samplenaive.counts,)), "BL5-H": panel1_n.data[['BL5-H']].values.reshape((samplenaive.counts,)),
                               "RL1-H": panel1_n.data[['RL1-H']].values.reshape((samplenaive.counts,))})
        df = df.append(df_add)
        panel1_m = panel1.copy()
        samplemem = samplethelp.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'Mem Th') &
                                                                   (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        gated_idx = np.array(samplemem.data.index)
        panel1_m.set_data(panel1.data.loc[gated_idx])
        df_add = pd.DataFrame({"Cell Type": np.tile("Mem Th", samplemem.counts), "Date": np.tile(date, samplemem.counts), "Plate": np.tile(plate, samplemem.counts),
                               "VL1-H": panel1_m.data[['VL1-H']].values.reshape((samplemem.counts,)), "BL5-H": panel1_m.data[['BL5-H']].values.reshape((samplemem.counts,)),
                               "RL1-H": panel1_m.data[['RL1-H']].values.reshape((samplemem.counts,))})
        df = df.append(df_add)

    return df, unstainedWell


def treg_sample(date, plate, gates_df, mem_naive=False):
    """ Returns gated T-reg sample for a given date and plate. """
    # import data and create transformed df for gating
    panel1, _ = importF(date, plate, "A", 1)
    panel1_t = panel1.transform("tlog", channels=['VL6-H', 'VL4-H', 'BL1-H', 'VL1-H', 'BL3-H'])  # Creates copy of panel1 to transform and gate

    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])

    # implement gating, revert tlog, and add to dataframe
    samplecd3cd4 = panel1_t.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'CD3CD4') &
                                                               (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    sampletreg = samplecd3cd4.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'T-reg') &
                                                                 (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    # Gated signals based on gating values from csv
    gated_idx = np.array(sampletreg.data.index)
    panel1.set_data(panel1.data.loc[gated_idx])  # Selects only the corresponding data points from panel1(untransformed) based on gated points from panel1_t

    df_add = pd.DataFrame({"Cell Type": np.tile("T-reg", panel1.counts), "Date": np.tile(date, panel1.counts), "Plate": np.tile(plate, panel1.counts),
                           "VL1-H": panel1.data[['VL1-H']].values.reshape((panel1.counts,)), "BL5-H": panel1.data[['BL5-H']].values.reshape((panel1.counts,)),
                           "RL1-H": panel1.data[['RL1-H']].values.reshape((panel1.counts,))})
    df = df.append(df_add)

    # separates memory and naive populations and adds to dataframe
    if mem_naive:
        panel1_n = panel1.copy()
        samplenaive = sampletreg.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'Naive Treg') &
                                                                    (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        gated_idx = np.array(samplenaive.data.index)
        panel1_n.set_data(panel1.data.loc[gated_idx])
        df_add = pd.DataFrame({"Cell Type": np.tile("Naive Treg", samplenaive.counts), "Date": np.tile(date, samplenaive.counts), "Plate": np.tile(plate, samplenaive.counts),
                               "VL1-H": panel1_n.data[['VL1-H']].values.reshape((samplenaive.counts,)), "BL5-H": panel1_n.data[['BL5-H']].values.reshape((samplenaive.counts,)),
                               "RL1-H": panel1_n.data[['RL1-H']].values.reshape((samplenaive.counts,))})
        df = df.append(df_add)
        panel1_m = panel1.copy()
        samplemem = sampletreg.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'Mem Th') &
                                                                  (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        gated_idx = np.array(samplemem.data.index)  # Selects only the corresponding data points from panel1(untransformed) based on gated points from panel1_t
        panel1_m.set_data(panel1.data.loc[gated_idx])
        df_add = pd.DataFrame({"Cell Type": np.tile("Mem Treg", samplemem.counts), "Date": np.tile(date, samplemem.counts), "Plate": np.tile(plate, samplemem.counts),
                               "VL1-H": panel1_m.data[['VL1-H']].values.reshape((samplemem.counts,)), "BL5-H": panel1_m.data[['BL5-H']].values.reshape((samplemem.counts,)),
                               "RL1-H": panel1_m.data[['RL1-H']].values.reshape((samplemem.counts,))})
        df = df.append(df_add)

    return df


def nk_nkt_sample(date, plate, gates_df, nkt=False):
    """ Returns gated NK sample for a given date and plate. """
    # import data and create transformed df for gating
    panel2, _ = importF(date, plate, "B", 2)
    panel2_t = panel2.transform("tlog", channels=['VL4-H', 'BL3-H'])  # Creates copy of panel1 to transform and gate

    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])

    # implement gating, revert tlog, and add to dataframe
    samplenk = panel2_t.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'NK') &
                                                           (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    # Gated signals based on gating values from csv
    panel2_1 = panel2.copy()  # Not sure why this copy was created and used(maybe debugging) have to ask Zoe
    gated_idx = np.array(samplenk.data.index)  # Selects only the corresponding data points from panel1(untransformed) based on gated points from panel1_t
    panel2_1.set_data(panel2.data.loc[gated_idx])
    df_add = pd.DataFrame({"Cell Type": np.tile("NK", samplenk.counts), "Date": np.tile(date, samplenk.counts), "Plate": np.tile(plate, samplenk.counts),
                           "VL1-H": panel2_1.data[['VL1-H']].values.reshape((samplenk.counts,)), "BL5-H": panel2_1.data[['BL5-H']].values.reshape((samplenk.counts,)),
                           "RL1-H": panel2_1.data[['RL1-H']].values.reshape((samplenk.counts,))})
    df = df.append(df_add)

    # gates NKT population and adds to dataframe
    if nkt:
        samplenkt = panel2_t.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'NKT') &
                                                                (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        panel2_2 = panel2.copy()
        gated_idx = np.array(samplenkt.data.index)
        panel2_2.set_data(panel2.data.loc[gated_idx])
        df_add = pd.DataFrame({"Cell Type": np.tile("NKT", samplenkt.counts), "Date": np.tile(date, samplenkt.counts), "Plate": np.tile(plate, samplenkt.counts),
                               "VL1-H": panel2_2.data[['VL1-H']].values.reshape((samplenkt.counts,)), "BL5-H": panel2_2.data[['BL5-H']].values.reshape((samplenkt.counts,)),
                               "RL1-H": panel2_2.data[['RL1-H']].values.reshape((samplenkt.counts,))})
        df = df.append(df_add)

    return df


def cd8_sample(date, plate, gates_df, mem_naive=False):
    """ Returns gated CD8+ sample for a given date and plate. """
    # import data and create transformed df for gating
    panel3, _ = importF(date, plate, "C", 3)
    panel3_t = panel3.transform("tlog", channels=['VL4-H', 'VL6-H', 'BL3-H'])  # Creates copy of panel1 to transform and gate

    df = pd.DataFrame(columns=["Cell Type", "Date", "Plate", "VL1-H", "BL5-H", "RL1-H"])

    # implement gating, revert tlog, and add to dataframe
    samplecd8 = panel3_t.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'CD8+') &
                                                            (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
    # Gated signals based on gating values from csv
    gated_idx = np.array(samplecd8.data.index)
    panel3.set_data(panel3.data.loc[gated_idx])  # Selects only the corresponding data points from panel1(untransformed) based on gated points from panel1_t
    df_add = pd.DataFrame({"Cell Type": np.tile("CD8+", samplecd8.counts), "Date": np.tile(date, samplecd8.counts), "Plate": np.tile(plate, samplecd8.counts),
                           "VL1-H": panel3.data[['VL1-H']].values.reshape((samplecd8.counts,)), "BL5-H": panel3.data[['BL5-H']].values.reshape((samplecd8.counts,)),
                           "RL1-H": panel3.data[['RL1-H']].values.reshape((samplecd8.counts,))})
    df = df.append(df_add)

    # separates memory and naive populations and adds to dataframe
    if mem_naive:
        panel3_n = panel3.copy()
        samplenaive = samplecd8.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'Naive CD8+') &
                                                                   (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        gated_idx = np.array(samplenaive.data.index)
        panel3_n.set_data(panel3.data.loc[gated_idx])
        df_add = pd.DataFrame({"Cell Type": np.tile("Naive CD8+", samplenaive.counts), "Date": np.tile(date, samplenaive.counts), "Plate": np.tile(plate, samplenaive.counts),
                               "VL1-H": panel3_n.data[['VL1-H']].values.reshape((samplenaive.counts,)), "BL5-H": panel3_n.data[['BL5-H']].values.reshape((samplenaive.counts,)),
                               "RL1-H": panel3_n.data[['RL1-H']].values.reshape((samplenaive.counts,))})
        df = df.append(df_add)
        panel3_m = panel3.copy()
        samplemem = samplecd8.gate(ast.literal_eval(gates_df.loc[(gates_df["Name"] == 'Mem CD8+') &
                                                                 (gates_df["Date"] == date) & (gates_df["Plate"] == float(plate))]["Gate"].values[0]))
        gated_idx = np.array(samplemem.data.index)
        panel3_m.set_data(panel3.data.loc[gated_idx])
        df_add = pd.DataFrame({"Cell Type": np.tile("Mem CD8+", samplemem.counts), "Date": np.tile(date, samplemem.counts), "Plate": np.tile(plate, samplemem.counts),
                               "VL1-H": panel3_m.data[['VL1-H']].values.reshape((samplemem.counts,)), "BL5-H": panel3_m.data[['BL5-H']].values.reshape((samplemem.counts,)),
                               "RL1-H": panel3_m.data[['RL1-H']].values.reshape((samplemem.counts,))})
        df = df.append(df_add)

    return df
