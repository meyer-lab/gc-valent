"""
This file includes various methods for flow cytometry analysis of fixed cells.
"""
import os
from os.path import dirname
from pathlib import Path
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt
from FlowCytometryTools import FCMeasurement
from FlowCytometryTools import PolyGate

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
        combinedSamples = subtract_unstained_signal(combinedSamples, channels, unstainedWell)  # Subtracts background
        return combinedSamples.transform("hlog", channels=channels)  # Transforms and returns

    tsample = subtract_unstained_signal(sample[wellNum - 1], channels, unstainedWell)
    return tsample.transform("hlog", channels=channels)


def subtract_unstained_signal(sample, channels, unstainedWell):
    """ Subtract mean unstained signal from all input channels for a given sample. """
    for _, channel in enumerate(channels):
        meanBackground = np.mean(unstainedWell.data[channel])  # Calculates mean unstained signal for given channel
        sample[channel] = np.maximum(sample[channel] - meanBackground, 0.0)

    return sample


# *********************************** Gating Fxns *******************************************

# Panel 1 gates:naive and memory T-regulatory and T-helper cells
def cd3cd4():
    """Function for gating CD3+CD4+ cells (generates T cells)"""
    cd3cd4_gate = PolyGate([(5.2e03, 5.4e03), (5.8e03, 7.9e03), (7.2e03, 7.9e03), (7.2e03, 5.7e03)], ("VL6-H", "VL4-H"), region="in", name="cd3cd4")
    return cd3cd4_gate


def thelper():
    """Fucntion creating and returning T helper gate on CD3+CD4+ cells"""
    thelp_gate = PolyGate([(3.2e03, 1e03), (3.2e03, 4.4e03), (5.2e03, 5.55e03), (6.4e03, 5.55e03), (6.4e03, 1e03)], ("BL1-H", "VL1-H"), region="in", name="thelp")
    return thelp_gate


def treg():
    """Function creating and returning the T reg gate on CD3+CD4+ cells"""
    treg_gate = PolyGate([(2.9e03, 4.5e03), (4.8e03, 5.9e03), (3.5e03, 6.3e03), (9.0e02, 6.3e03), (9.0e02, 4.5e03)], ("BL1-H", "VL1-H"), region="in", name="treg")
    return treg_gate


# Panel 2 gates: NK and CD56bright NK cells

# Panel 3 gates: naive and memory cytotoxic T cells

# Panel 4 TODO


def plot_Tcells(sample, cd3cd4gate, Thelpgate, Treggate):
    """Plotting naïve and memory T-regulatory and T-helper cells. Input transformed sample and gate functions for arguments"""

    fig, axs = plt.subplots(2, 2)
    sample.plot(["VL4-H", "VL6-H"], cmap=cm.viridis, gates=cd3cd4gate, gate_lw=2, ax=axs[0, 0])
    axs[0, 0].set(xlabel="CD3", ylabel="CD4", title="Singlet Lymphocytes")

    cd3cd4gated_sample = sample.gate(cd3cd4gate)
    cd3cd4gated_sample.plot(["VL1-H", "BL1-H"], cmap=cm.viridis, gates=(Thelpgate, Treggate), gate_lw=2, ax=axs[0, 1])
    axs[0, 1].set(xlabel="CD25", ylabel="CD127", title="CD3+CD4+ Cells")
    axs[0, 1].set_xlim(right=7000)
    axs[0, 1].set_ylim(top=7000)

    ThelpGated_sample = cd3cd4gated_sample.gate(Thelpgate)
    ThelpGated_sample.plot(["BL3-H"], color="blue", ax=axs[1, 0])
    axs[1, 0].set(xlabel="CD45Ra", title="T helper")

    TregGated_sample = cd3cd4gated_sample.gate(Treggate)
    TregGated_sample.plot(["BL3-H"], color="blue", ax=axs[1, 1])
    axs[1, 1].set(xlabel="CD45Ra", title="T reg")

    #compSample = applyMatrix(sample, compMatrix(date, plate, wellRow))
    return sample, unstainedWell
