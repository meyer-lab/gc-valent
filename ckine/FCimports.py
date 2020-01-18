from pylab import *
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import scale
from FlowCytometryTools import FCMeasurement, FCPlate
from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import QuadGate, ThresholdGate, PolyGate
from sklearn import preprocessing

def importF(pathname, WellRow):
    """
    Import FCS files. Variable input: name of path name to file. Output is a list of Data File Names in FCT Format
    Title/file names are returned in the array file --> later referenced in other functions as title/titles input argument
    """
    # Declare arrays and int
    file = []
    sample = []
    z = 0
    # Read in user input for file path and assign to array file
    pathlist = Path(r"" + str(pathname)).glob("**/*.fcs")
    for path in pathlist:
        wellID = path.name.split("_")[1]
        if wellID[0] == WellRow:
            file.append(str(path))
    file.sort()
    assert file != []
    # Go through each file and assign the file contents to entry in the array sample
    for entry in file:
        sample.append(FCMeasurement(ID="Test Sample" + str(z), datafile=entry))
        z += 1
    # Returns the array sample which contains data of each file in folder (one file per entry in array)
    return sample, file

def combineWells(samples, channels_):
    """Accepts sample array returned from importF, and array of channels, returns transformed combined well data"""
    combinedSamples = samples[0]
    for sample in samples1A[1:]:
        combinedSamples.data = combinedSamples.data.append(sample.data)
    t_combinedSamples = combinedSamples.transform('hlog', channels = channels_)
    return t_combinedSamples


# *********************************** Gating Fxns *******************************************

#Panel 1 gates:naive and memory T-regulatory and T-helper cells
def cd3cd4():
    """Function for gating CD3+CD4+ cells (generates T cells)"""
    cd3cd4 = PolyGate([(5.2e03,5.4e03),(5.8e03,7.9e03),(7.2e03,7.9e03),(7.2e03,5.7e03)], ('VL6-H','VL4-H'),region='in', name = 'cd3cd4')
    cd3cd4_gate = cd3cd4
    return cd3cd4_gate

def thelper():
    """Fucntion creating and returning T helper gate on CD3+CD4+ cells"""
    thelp = PolyGate([(3.2e03, 1e03), (3.2e03, 4.4e03), (5.2e03, 5.55e03), (6.4e03,5.55e03 ), (6.4e03, 1e03)], ('BL1-H','VL1-H'),region = 'in', name = 'thelp')
    thelp_gate = thelp 
    return thelp_gate

def treg():
    """Function creating and returning the T reg gate on CD3+CD4+ cells"""
    treg = PolyGate([(2.9e03, 4.5e03), (4.8e03, 5.9e03), (3.5e03, 6.3e03), (9.0e02,6.3e03 ), (9.0e02, 4.5e03)], ('BL1-H', 'VL1-H'), region='in', name='treg')
    treg_gate = treg
    return treg_gate

def tregMem():#switch to polygates
    """Function for creating and returning the Treg Memory gate on Treg CD3+CD4+ cells"""
    treg1 = QuadGate((1e+03, 1e+03), ('BL1-H', 'VL1-H'), region='top right', name='treg1')
    treg2 = QuadGate((0, 0), ('BL1-H', 'VL1-H'), region='bottom left', name='treg2')
    cd45 = ThresholdGate(6300, ('BL3-H'), region="below", name='cd45')
    tregMem_gate = treg1 & treg2 & cd3cd4() & cd45
    return tregMem_gate

def tregN():#switch to polygates
    """Function for creating and returning the T reg gate on CD4+ cells"""
    treg1 = QuadGate((1e+03, 1e+03), ('BL1-H', 'VL1-H'), region='top right', name='treg1')
    treg2 = QuadGate((0, 0), ('BL1-H', 'VL1-H'), region='bottom left', name='treg2')
    cd45 = ThresholdGate(6300, ('BL3-H'), region="above", name='cd45')
    tregN_gate = treg1 & treg2 & cd3cd4() & cd45
    return tregN_gate

def THelpMem():#switch to polygates
    """Function for creating and returning the non T reg gate on CD4+ cells"""
    thelp1 = QuadGate((2e+03,2e+03), ('BL1-H','VL1-H'),region='top right', name = 'thelp1')
    thelp2 = QuadGate((0,0), ('BL1-H','VL1-H'),region='bottom left', name = 'thelp2')
    cd45 = ThresholdGate(6300, ('BL3-H'), region="below", name='cd45')
    Thelp_gate = thelper() & cd3cd4() & cd45
    return Thelp_gate

def THelpN():#switch to polygates
    """Function for creating and returning the non T reg gate on CD4+ cells"""
    thelp1 = QuadGate((5.115e+03, 3.470e+02), ('BL1-H', 'VL1-H'), region="top left", name='thelp1')
    thelp2 = QuadGate((2.586e+03, 5.245e+03), ('BL1-H', 'VL1-H'), region="bottom right", name='thelp2')
    cd45 = ThresholdGate(6300, ('BL3-H'), region="above", name='cd45')
    ThelpN_gate = thelp1 & thelp2 & cd3cd4() & cd45
    return ThelpN_gate

#Panel 2 gates: NK and CD56bright NK cells

#Panel 3 gates: naive and memory cytotoxic T cells

#Panel 4 TODO

def plot_Tcells(sample, cd3cd4gate, Thelpgate, Treggate):
    """Plotting naïve and memory T-regulatory and T-helper cells. Input transformed sample and gate functions for arguments"""
    
    fig = figure(figsize=(12,9))
    fig.subplots_adjust(hspace=.25)
    ax1 = subplot(221)
    sample.plot(['VL4-H', 'VL6-H'], cmap=cm.viridis, gates=cd3cd4gate, gate_lw=2)
    title('Singlet Lymphocytes')
    ax1.set_ylabel('CD4')
    ax1.set_xlabel('CD3')
    
    cd3cd4gated_sample = sample.gate(cd3cd4gate)
    
    ax2 = subplot(222)
    cd3cd4gated_sample.plot(['VL1-H','BL1-H'], cmap=cm.viridis, gates=(Thelpgate,Treggate), gate_lw=2)
    title('CD3+CD4+ Cells')
    ax2.set_ylabel('CD127')
    ax2.set_xlabel('CD25')
    ax2.set_xlim(right=7000)
    ax2.set_ylim(top=7000)
    
    ThelpGated_sample = cd3cd4gated_sample.gate(Thelpgate)
    TregGated_sample = cd3cd4gated_sample.gate(Treggate)

    ax3 = subplot(223)
    ThelpGated_sample.plot(['BL3-H'], color='blue');
    title('T helper')
    ax3.set_xlabel('CD45Ra');

    ax4= subplot(224)
    TregGated_sample.plot(['BL3-H'], color='blue');
    title('T reg')
    ax4.set_xlabel('CD45Ra');
    
    return fig