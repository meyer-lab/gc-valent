import os
import numpy as np
import pandas as pd
from .figureCommon import subplotLabel, getSetup
from ..imports import channels
from ..FCimports import combineWells, compMatrix, applyMatrix, import_gates, apply_gates
from FlowCytometryTools import FCMeasurement

path_here = os.path.dirname(os.path.dirname(__file__))

def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)

    Tcell_pathname = path_here + "/ckine/data/flow/2019-11-08 monomer IL-2 Fc signaling/CD4 T cells - IL2-060 mono, IL2-060 dimeric"
    NK_CD8_pathname = path_here + "/ckine/data/flow/2019-11-08 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-060 mono, IL2-060 dimeric"

    Tcell_sample = importF(Tcell_pathname, "A")
    NK_CD8_sample = importF(NK_CD8_pathname, "A")

    Tcell_sample = combineWells(Tcell_sample)
    NK_CD8_sample = combineWells(NK_CD8_sample)

    Tcell_matrixPath = path_here + "/ckine/data/compensation/CD4+ comp matrix.csv"
    Cd8_NKmatrixPath = path_here + "/ckine/data/compensation/CD8-CD56 comp matrix.csv"

    Tcell_sample = applyMatrix(Tcell_sample, compMatrix(Tcell_matrixPath))
    NK_CD8_sample = applyMatrix(NK_CD8_sample, compMatrix(Cd8_NKmatrixPath))

    Tcell_sample = Tcell_sample.transform("tlog", channels=['VL1-H', 'VL4-H', 'BL1-H','BL3-H']) #Tlog transformations
    NK_CD8_sample = NK_CD8_sample.transform("tlog", channels=['RL1-H', 'VL4-H', 'BL1-H', 'BL2-H']) #Tlog transformations

    
    
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