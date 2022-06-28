"""
This file includes various methods for flow cytometry analysis of fixed cells.
"""
import os
from os.path import dirname, join
from pathlib import Path
import ast
import textwrap
import pandas as pd
import numpy as np
import warnings
from FlowCytometryTools import PolyGate, FCMeasurement

path_here = os.path.dirname(os.path.dirname(__file__))


warnings.filterwarnings("ignore")
gate_df = pd.DataFrame()


flow_dict = {"Pacific Blue-A": "FoxP3",
             "FITC-A": "CD8",
             "PE-A": "CD4",
             "BV605-A": "CD56",
             "Alexa Fluor 700-A": "pSTAT5",
             "APC-Cy7-A": "CD25",
             "Alexa Fluor 647-A": "CD3"}


def compile_untreated(cellFrac):
    """Adds all data from a single patient to an FC file"""
    pathname = join(path_here, "ckine/data/Flow_Data_Meyer/Untreated/")
    pathlist = Path(r"" + str(pathname)).glob("**/*.fcs")
    FCfiles = []
    for path in pathlist:
        FCfiles.append(FCMeasurement(ID="All Data", datafile=path))
    return combineWells(FCfiles, cellFrac)


def combineWells(samples, cellFrac):
    """Accepts sample array returned from importF, and array of channels, returns combined well data"""
    comp_mat = compMatrix()
    markers = np.array(["Pacific Blue-A", "FITC-A", "PE-A", "BV605-A", "Alexa Fluor 700-A", "APC-Cy7-A", "Alexa Fluor 647-A"])
    log_markers = markers[np.isin(markers, samples[0].data.columns)]
    samples[0] = applyMatrix(samples[0], comp_mat)
    samples[0] = samples[0].transform("tlog", channels=log_markers)
    combinedSamples = samples[0]
    for sample in samples[1:]:
        log_markers = markers[np.isin(markers, sample.data.columns)]
        sample = applyMatrix(sample, comp_mat)
        sample = sample.transform("tlog", channels=log_markers)
        combinedSamples.data = pd.concat([combinedSamples.data, sample.data.sample(frac=cellFrac)])
    combinedSamples.data = combinedSamples.data.rename(flow_dict, axis=1)
    combinedSamples.data = combinedSamples.data.astype(int)
    combinedSamples.data = combinedSamples.data.clip(lower=0, upper=5e5)
    return combinedSamples


def compMatrix(invert=True):
    """Creates compensation matrix given parameters date in mm-dd, plate number and panel A, B, or C."""
    path = join(path_here, "ckine/data/Flow_Data_Meyer/Compensation.csv")
    # imports csv file with comp values as a dataframe
    header_names = ['Channel1', 'Channel2', 'Comp']
    df_comp = pd.read_csv(path, header=None, index_col=False, skiprows=0, names=header_names)
    # Add diangonal values of 100 to compensation values
    addedChannels = []
    for i in df_comp.index:
        channelName = df_comp.iloc[i]['Channel1']
        if channelName not in addedChannels:  # Ensures a diagonal value is only added once for each channel
            addedChannels.append(channelName)
            df2 = pd.DataFrame([[channelName, channelName, 100]], columns=['Channel1', 'Channel2', 'Comp'])  # Creates new row for dataframe
            df_comp = pd.concat([df_comp, df2], ignore_index=True)  # Adds row
    # Create square matrix from compensation values
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


def process_sample(sample):
    """Relabels and logs a sample"""
    comp_mat = compMatrix()
    sample = applyMatrix(sample, comp_mat)
    markers = np.array(["Pacific Blue-A", "FITC-A", "PE-A", "BV605-A", "APC-Cy7-A", "Alexa Fluor 647-A"])
    log_markers = markers[np.isin(markers, sample.data.columns)]
    sample = sample.transform("tlog", channels=log_markers)
    sample.data = sample.data.rename(flow_dict, axis=1)
    return sample, log_markers


def makeGate(lowerCorner, upperCorner, channels, name):
    """Returns square gate using upper and lower corners"""
    return PolyGate([(lowerCorner[0], lowerCorner[1]), (upperCorner[0], lowerCorner[1]), (upperCorner[0], upperCorner[1]), (lowerCorner[0], upperCorner[1])], channels, region='in', name=name)


cellTypes = ["NK", "CD8+", "Thelper", "Treg"]

gate_dict = {"NK": ["NK Gate"],
             "CD8+": ["CD8 Gate"],
             "Treg": ["CD4 Gate", "Treg Gate"],
             "Thelper": ["CD4 Gate", "Thelper Gate"]}


def form_gate(gate):
    """Deconvolutes string flow gate object"""
    vertices = ast.literal_eval(textwrap.dedent(str(gate).split("Vertices: ")[1].split("Channel")[0]))
    channels = ast.literal_eval(textwrap.dedent(str(gate).split("Channel(s): ")[1].split("Name")[0]))
    return PolyGate(vertices, channels)


def live_PBMC_gate(sample, gateDF):
    """Returns singlet lymphocyte live PBMCs for a patient"""
    gates = ["PBMC Gate", "Live Gate"]
    for gate_name in gates:
        gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == gate_name)].Gate.values[0])
        sample = sample.gate(gate)
    return sample


def pop_gate(sample, cell_type, gateDF):
    """Extracts cell population sample"""
    gates = gate_dict[cell_type]
    for gate_name in gates:
        gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == gate_name)].Gate.values[0])
        pop_sample = sample.gate(gate)
    return pop_sample


lig_dict = {"01": ["R38Q/H16N N-term", 1],
            "02": ["R38Q/H16N N-term", 1],
            "03": ["R38Q/H16N N-term", 2],
            "04": ["R38Q/H16N N-term", 2],
            "05": ["R38Q/H16N N-term", 4],
            "06": ["R38Q/H16N N-term", 4],
            "07": ["Live/Dead", 2],
            "08": ["Live/Dead", 2],
            "09": ["Live/Dead", 4],
            "10": ["Live/Dead", 4]}

dose_dict = {"A": 50, "B": 12.954, "C": 3.355, "D": 0.869, "E": 0.2252, "F": 0.058, "G": 0.01512, "H": 0.00392, "I": 0.001}


def make_flow_df(subtract=True):
    """Compiles data for all populations for all patients into .csv"""
    columns = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    rows = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]
    rowminus = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]

    cell_types = ["Treg", "Thelper", "CD8+", "NK"]
    gateDF = pd.read_csv(join(path_here, "ckine/data/Meyer_Flow_Gates.csv"))
    MeyerDF = pd.DataFrame([])

    for column in columns:
        for i, row in enumerate(rows):
            sample = FCMeasurement(ID=column + row, datafile=(join(path_here, "ckine/data/Flow_Data_Meyer/Specimen_001_" + column + rowminus[i] + "_" + column + row + ".fcs")))
            sample, _ = process_sample(sample)
            sample = live_PBMC_gate(sample, gateDF)
            for cell_type in cell_types:
                pop_sample = pop_gate(sample, cell_type, gateDF)
                mean = pop_sample.data["pSTAT5"]
                mean = np.mean(mean.values[mean.values < np.quantile(mean.values, 0.995)])
                MeyerDF = pd.concat([MeyerDF, pd.DataFrame({"Ligand": lig_dict[row][0], "Valency": lig_dict[row][1], "Dose": dose_dict[column], "Cell": cell_type, "pSTAT5": [mean]})])
    MeyerDF = MeyerDF.groupby(["Ligand", "Valency", "Dose", "Cell"]).pSTAT5.mean().reset_index()

    untreatedDF = pd.DataFrame()
    untreated_sample = compile_untreated(cellFrac=1.0)
    untreated_sample = live_PBMC_gate(sample, gateDF)
    for cell_type in cell_types:
        pop_sample = pop_gate(untreated_sample, cell_type, gateDF)
        mean = pop_sample.data["pSTAT5"]
        mean = np.mean(mean.values[mean.values < np.quantile(mean.values, 0.995)])
        untreatedDF = pd.concat([untreatedDF, pd.DataFrame({"Cell": cell_type, "pSTAT5": [mean]})])

    for ligand in MeyerDF.Ligand.unique():
        for valency in MeyerDF.loc[MeyerDF.Ligand == ligand].Valency.unique():
            for dose in MeyerDF.Dose.unique():
                for cell in MeyerDF.Cell.unique():
                    MeyerDF.loc[(MeyerDF.Valency == valency) & (MeyerDF.Ligand == ligand) & (MeyerDF.Cell == cell) & (
                        MeyerDF.Dose == dose), "pSTAT5"] -= untreatedDF.loc[untreatedDF.Cell == cell].pSTAT5.values
    MeyerDF.pSTAT5 = MeyerDF.pSTAT5.clip(lower=0)
    MeyerDF.to_csv(join(path_here, "ckine/data/Meyer_Flow.csv"))
    return MeyerDF
