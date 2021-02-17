import os
from os.path import dirname, join
from pathlib import Path
import numpy as np
import pandas as pd
from .figureCommon import subplotLabel, getSetup
from ..imports import channels
from ..FCimports import combineWells, compMatrix, applyMatrix, import_gates, apply_gates
from FlowCytometryTools import FCMeasurement, ThresholdGate, PolyGate, QuadGate
from matplotlib import pyplot as plt

path_here = os.path.dirname(os.path.dirname(__file__))

def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)

    Tcell_pathname = path_here + "/data/flow/2019-11-08 monomer IL-2 Fc signaling/CD4 T cells - IL2-060 mono, IL2-060 dimeric"
    NK_CD8_pathname = path_here + "/data/flow/2019-11-08 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-060 mono, IL2-060 dimeric"

    print(Tcell_pathname)

    Tcell_sample, _ = importF(Tcell_pathname, "A")
    NK_CD8_sample, _ = importF(NK_CD8_pathname, "A")

    Tcell_sample = combineWells(Tcell_sample)
    NK_CD8_sample = combineWells(NK_CD8_sample)

    Tcell_sample = applyMatrix(Tcell_sample, compMatrix('2019-11-08','1','A'))
    NK_CD8_sample = applyMatrix(NK_CD8_sample, compMatrix('2019-11-08','1','B'))

    Tcell_sample = Tcell_sample.transform("tlog", channels=['VL1-H', 'VL4-H', 'BL1-H','BL3-H']) #Tlog transformations
    NK_CD8_sample = NK_CD8_sample.transform("tlog", channels=['RL1-H', 'VL4-H', 'BL1-H', 'BL2-H']) #Tlog transformations

    ax[0]
    cd4_gate = ThresholdGate(6500.0, ['VL4-H'], region='above') & ThresholdGate(8000.0, ['VL4-H'], region='below')
    _ = Tcell_sample.plot(['VL4-H'], gates=cd4_gate) #CD4
    plt.title("Singlet Lymphocytes")
    plt.xlabel("CD4")
    plt.ylabel("Events")
    plt.grid()

    
    return f


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