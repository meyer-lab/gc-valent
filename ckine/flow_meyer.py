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
from copy import copy
from FlowCytometryTools import PolyGate, FCMeasurement

path_here = os.path.dirname(os.path.dirname(__file__))


warnings.filterwarnings("ignore")
gate_df = pd.DataFrame()


flow_dict_1 = {"Pacific Blue-A": "FoxP3",
"FITC-A": "CD8",
"PE-A": "CD4",
"BV605-A": "CD56",
"Alexa Fluor 700-A": "pSTAT5",
"APC-Cy7-A": "CD25",
"Alexa Fluor 647-A": "CD3"}


flow_dict_2 = {"Pacific Blue-A": "FoxP3",
"FITC-A": "CD8",
"PE-A": "CD4",
"BV605-A": "CD56",
"Alexa Fluor 700-A": "pSTAT5",
"APC-Cy7-A": "CD25",
"PerCP-Cy5-5-A": "CD3"}


def compile_untreated(date, cellFrac):
    """Adds all data from a single patient to an FC file"""
    pathname = join(path_here, "ckine/data/Flow_Data_Meyer/" + date + "/Untreated/")
    pathlist = Path(r"" + str(pathname)).glob("**/*.fcs")
    FCfiles = []
    for path in pathlist:
        FCfiles.append(FCMeasurement(ID="All Data", datafile=path))
    return combineWells(FCfiles, cellFrac, date)


def combineWells(samples, cellFrac, date):
    """Accepts sample array returned from importF, and array of channels, returns combined well data"""
    comp_mat = compMatrix(date)
    markers = np.array(["Pacific Blue-A", "FITC-A", "PE-A", "BV605-A", "Alexa Fluor 700-A", "APC-Cy7-A", "Alexa Fluor 647-A", "PerCP-Cy5-5-A"])
    log_markers = markers[np.isin(markers, samples[0].data.columns)]
    samples[0] = applyMatrix(samples[0], comp_mat)
    samples[0] = samples[0].transform("tlog", channels=log_markers)
    combinedSamples = samples[0]
    for sample in samples[1:]:
        log_markers = markers[np.isin(markers, sample.data.columns)]
        sample = applyMatrix(sample, comp_mat)
        sample = sample.transform("tlog", channels=log_markers)
        combinedSamples.data = pd.concat([combinedSamples.data, sample.data.sample(frac=cellFrac)])
    if date == "3_17_23":
        combinedSamples.data = combinedSamples.data.rename(flow_dict_2, axis=1)
    else:
        combinedSamples.data = combinedSamples.data.rename(flow_dict_1, axis=1)
    combinedSamples.data = combinedSamples.data.astype(int)
    combinedSamples.data = combinedSamples.data.clip(lower=0, upper=5e5)
    return combinedSamples


def compMatrix(date, invert=True):
    """Creates compensation matrix given parameters date in mm-dd, plate number and panel A, B, or C."""
    path = join(path_here, "ckine/data/Flow_Data_Meyer/" + date + "/Compensation.csv")
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


def process_sample(sample, date):
    """Relabels and logs a sample"""
    comp_mat = compMatrix(date)
    sample = applyMatrix(sample, comp_mat)
    markers = np.array(["Pacific Blue-A", "FITC-A", "PE-A", "BV605-A", "APC-Cy7-A", "Alexa Fluor 647-A", "PerCP-Cy5-5-A"])
    log_markers = markers[np.isin(markers, sample.data.columns)]
    sample = sample.transform("tlog", channels=log_markers)
    if date == "3_17_23":
        sample.data = sample.data.rename(flow_dict_2, axis=1)
    else:
        sample.data = sample.data.rename(flow_dict_1, axis=1)
    return sample, log_markers


def makeGate(lowerCorner, upperCorner, channels, name):
    """Returns square gate using upper and lower corners"""
    return PolyGate([(lowerCorner[0], lowerCorner[1]), (upperCorner[0], lowerCorner[1]), (upperCorner[0], upperCorner[1]), (lowerCorner[0], upperCorner[1])], channels, region='in', name=name)


gate_dict = {"NK": ["NK Gate"],
             "NKBright": ["NK Bright Gate"],
             "CD8": ["CD8 Gate"],
             "Treg": ["CD4 Gate", "Treg Gate"],
             "Thelper": ["CD4 Gate", "Thelper Gate"]}


def form_gate(gate):
    """Deconvolutes string flow gate object"""
    vertices = ast.literal_eval(textwrap.dedent(str(gate).split("Vertices: ")[1].split("Channel")[0]))
    channels = ast.literal_eval(textwrap.dedent(str(gate).split("Channel(s): ")[1].split("Name")[0]))
    return PolyGate(vertices, channels)


def live_PBMC_gate(sample, gateDF, date):
    """Returns singlet lymphocyte live PBMCs for a patient"""
    gates = ["PBMC Gate", "Live Gate"]
    gateDF = gateDF.loc[gateDF.Date == date]
    for gate_name in gates:
        gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == gate_name)].Gate.values[0])
        sample = sample.gate(gate)
    return sample


def pop_gate(sample, cell_type, gateDF, date):
    """Extracts cell population sample"""
    gateDF = gateDF.loc[gateDF.Date == date]
    gates = gate_dict[cell_type]
    pop_sample = copy(sample)
    for gate_name in gates:
        gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == gate_name)].Gate.values[0])
        pop_sample = pop_sample.gate(gate)
    return pop_sample


lig_dict = {"01": ["R38Q/H16N", 1, 1],
            "02": ["R38Q/H16N", 1, 2],
            "03": ["R38Q/H16N", 2, 1],
            "04": ["R38Q/H16N", 2, 2],
            "05": ["R38Q/H16N", 4, 1],
            "06": ["R38Q/H16N", 4, 2],
            "07": ["Live/Dead", 2, 1],
            "08": ["Live/Dead", 2, 2],
            "09": ["Live/Dead", 4, 1],
            "10": ["Live/Dead", 4, 2]}

dose_dict = {"A": 50, "B": 12.954, "C": 3.355, "D": 0.869, "E": 0.2252, "F": 0.058, "G": 0.01512, "H": 0.00392, "I": 0.001}


def make_flow_df(subtract=True):
    """Compiles data for all populations for all patients into .csv"""
    columns = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    rows = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]
    rowminus = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    dates = ["6_25_22", "7_22_22", "8_10_22", "9_19_22", "3_17_23", "3_24_23_1", "3_24_23_2"]
    datesFormat = ["6/25/22", "7/22/22", "8/10/22", "9/19/22", "3/17/23", "3/24/23_1", "3/24/23_2"]

    cell_types = ["Treg", "Thelper", "CD8", "NK", "NKBright"]
    gateDF = pd.read_csv(join(path_here, "ckine/data/Meyer_Flow_Gates.csv"))
    MeyerDF = pd.DataFrame([])

    for j, date in enumerate(dates):
        print(date)
        for column in columns:
            for i, row in enumerate(rows):
                sample = FCMeasurement(ID=column + row, datafile=(join(path_here, "ckine/data/Flow_Data_Meyer/" + date + "/Specimen_001_" + column + rowminus[i] + "_" + column + row + ".fcs")))
                sample, _ = process_sample(sample, date)
                sample = live_PBMC_gate(sample, gateDF, datesFormat[j])
                for cell_type in cell_types:
                    pop_sample = pop_gate(sample, cell_type, gateDF, datesFormat[j])
                    mean = pop_sample.data["pSTAT5"]
                    if len(mean) == 0:
                        mean = 0
                    else:
                        mean = np.mean(mean.values[mean.values < np.quantile(mean.values, 0.995)])
                    MeyerDF = pd.concat([MeyerDF, pd.DataFrame({"Ligand": lig_dict[row][0], "Valency": lig_dict[row][1],
                                        "Dose": dose_dict[column], "Cell": cell_type, "pSTAT5": [mean], "Date": datesFormat[j], "Replicate": lig_dict[row][2]})])
    # Remove cell-deficient wells
    MeyerDF = MeyerDF.loc[(MeyerDF.Date != "8/10/22") | (MeyerDF.Dose != 12.954) | (MeyerDF.Replicate != 2) | (MeyerDF.Ligand != "Live/Dead") | (MeyerDF.Valency != 4)]
    MeyerDF = MeyerDF.loc[(MeyerDF.Date != "8/10/22") | (MeyerDF.Dose != 0.869) | (MeyerDF.Replicate != 2) | (MeyerDF.Ligand != "Live/Dead") | (MeyerDF.Valency != 4)]
    MeyerDF = MeyerDF.loc[(MeyerDF.Date != "8/10/22") | (MeyerDF.Dose != 0.2252) | (MeyerDF.Ligand != "Live/Dead") | (MeyerDF.Valency != 4)]
    MeyerDF = MeyerDF.loc[(MeyerDF.Date != "9/19/22") | (MeyerDF.Dose != 0.00392) | (MeyerDF.Ligand != "Live/Dead") | (MeyerDF.Valency != 4)]

    MeyerDF = MeyerDF.groupby(["Ligand", "Valency", "Dose", "Cell", "Date"]).pSTAT5.mean().reset_index()

    untreatedDF = pd.DataFrame()
    for i, date in enumerate(dates):
        untreated_sample = compile_untreated(date, cellFrac=1.0)
        untreated_sample = live_PBMC_gate(untreated_sample, gateDF, datesFormat[i])
        for cell_type in cell_types:
            pop_sample = pop_gate(untreated_sample, cell_type, gateDF, datesFormat[i])
            mean = pop_sample.data["pSTAT5"]
            if len(mean) == 0:
                mean = 0
            else:
                mean = np.mean(mean.values[mean.values < np.quantile(mean.values, 0.995)])
            untreatedDF = pd.concat([untreatedDF, pd.DataFrame({"Cell": cell_type, "pSTAT5": [mean], "Date": datesFormat[i]})])
    untreatedDF = untreatedDF.fillna(value=0)
    untreatedDF.pSTAT5 = untreatedDF.pSTAT5.clip(lower=0)

    for date in datesFormat:
        for cell in MeyerDF.Cell.unique():
            # Remove persistent baseline values
            MeyerDF.loc[(MeyerDF.Cell == cell) & (MeyerDF.Date == date), "pSTAT5"] -= MeyerDF.loc[(MeyerDF.Cell == cell) & (
                MeyerDF.Dose <= 0.00392) & (MeyerDF.Date == date)].pSTAT5.mean()
            MeyerDF.pSTAT5 = MeyerDF.pSTAT5.clip(lower=0)
    MeyerDF.to_csv(join(path_here, "ckine/data/Meyer_Flow.csv"))
    return MeyerDF


flow_dict_ILC = {"Pacific Blue-A": "FoxP3",
"FITC-A": "Lin",
"Alexa Fluor 647-A": "FoxP3",
"BV605-A": "CRTH2",
"APC-Cy7-A": "CD25",
"BV421-A": "CD127"}


def combineWells_ILC(samples, cellFrac):
    """Accepts sample array returned from importF, and array of channels, returns combined well data"""
    comp_mat = compMatrix(ILCs=True)
    markers = np.array(["Pacific Blue-A", "FITC-A", "B421-A", "BV605-A", "APC-Cy7-A", "Alexa Fluor 647-A"])
    log_markers = markers[np.isin(markers, samples[0].data.columns)]
    samples[0] = applyMatrix(samples[0], comp_mat)
    samples[0] = samples[0].transform("tlog", channels=log_markers)
    combinedSamples = samples[0]
    for sample in samples[1:]:
        log_markers = markers[np.isin(markers, sample.data.columns)]
        sample = applyMatrix(sample, comp_mat)
        sample = sample.transform("tlog", channels=log_markers)
        combinedSamples.data = pd.concat([combinedSamples.data, sample.data.sample(frac=cellFrac)])
    
    combinedSamples.data = combinedSamples.data.rename(flow_dict_ILC, axis=1)
    combinedSamples.data = combinedSamples.data.astype(int)
    combinedSamples.data = combinedSamples.data.clip(lower=0, upper=5e5)
    return combinedSamples


def compMatrix_ILCs(invert=True):
    """Creates compensation matrix given parameters date in mm-dd, plate number and panel A, B, or C."""
    path = join(path_here, "ckine/data/Flow_Data_Meyer/ILCs/Compensation.csv")
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


def process_sample_ILCs(sample):
    """Relabels and logs a sample"""
    comp_mat = compMatrix_ILCs()
    sample = applyMatrix(sample, comp_mat)
    markers = np.array(["Pacific Blue-A", "FITC-A", "B421-A", "BV605-A", "APC-Cy7-A", "Alexa Fluor 647-A"])
    log_markers = markers[np.isin(markers, sample.data.columns)]
    sample = sample.transform("tlog", channels=log_markers)
    sample.data = sample.data.rename(flow_dict_ILC, axis=1)
    return sample, log_markers


gate_dict_ILCs = {
             "Treg": ["T Gate", "Treg Gate"],
             "ILC2": ["ILC Gate", "ILC2 Gate"]}


def live_PBMC_gate_ILCs(sample, gateDF, donor):
    """Returns singlet lymphocyte live PBMCs for a patient"""
    gates = ["PBMC Gate", "Live Gate"]
    gateDF = gateDF.loc[gateDF.Donor == donor]
    for gate_name in gates:
        gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == gate_name)].Gate.values[0])
        sample = sample.gate(gate)
    return sample


def pop_gate_ILCs(sample, cell_type, gateDF, donor):
    """Extracts cell population sample"""
    gateDF = gateDF.loc[gateDF.Donor == donor]
    gates = gate_dict_ILCs[cell_type]
    pop_sample = copy(sample)
    for gate_name in gates:
        gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == gate_name)].Gate.values[0])
        pop_sample = pop_sample.gate(gate)
    return pop_sample


def make_flow_df_ILCs(subtract=True):
    """Compiles data for all populations for all patients into .csv"""
    donors = ["Donor_303", "Donor_334", "Donor_850", "Donor_X"]
    donorFormat = ["Donor 303", "Donor 334", "Donor 850", "Donor X"]

    cell_types = ["Treg", "ILC2"]
    gateDF = pd.read_csv(join(path_here, "ckine/data/Meyer_Flow_Gates_ILCs.csv"))
    MeyerDF = pd.DataFrame([])
    MeyerDF_SC = pd.DataFrame([])

    for i, donor in enumerate(donors):
        print(donor)
        sample = FCMeasurement(ID=donor, datafile=(join(path_here, "ckine/data/Flow_Data_Meyer/ILCs/" + donor +  "/Specimen_001_" + donor + ".fcs")))
        sample, _ = process_sample_ILCs(sample)
        sample = live_PBMC_gate_ILCs(sample, gateDF, donorFormat[i].split(" ")[-1])
        for cell_type in cell_types:
            pop_sample = pop_gate_ILCs(sample, cell_type, gateDF, donorFormat[i].split(" ")[-1])
            mean = pop_sample.data["CD25"]
            mean = mean.loc[mean >= 5]
            MeyerDF_SC = pd.concat([MeyerDF_SC, pd.DataFrame({ "Cell": cell_type, "CD25": mean.values, "Donor": donorFormat[i]})])
            if len(mean) == 0:
                mean = 0
            else:
                mean = np.mean(mean.values[mean.values < np.quantile(mean.values, 0.995)])
            MeyerDF = pd.concat([MeyerDF, pd.DataFrame({ "Cell": cell_type, "CD25": [mean], "Donor": donorFormat[i]})])

    MeyerDF.to_csv(join(path_here, "ckine/data/Meyer_Flow_ILCs.csv"))
    MeyerDF_SC.to_csv(join(path_here, "ckine/data/Meyer_Flow_ILCs_SC.csv"))
    return MeyerDF
