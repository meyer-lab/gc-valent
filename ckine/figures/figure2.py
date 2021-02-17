import os
import numpy as np
import pandas as pd
from .figureCommon import subplotLabel, getSetup
from ..imports import channels
from ..FCimports import import_gates, apply_gates

path_here = os.path.dirname(os.path.dirname(__file__))

def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)

    Tcell_pathname = path_here + "/ckine/data/flow/2019-11-08 monomer IL-2 Fc signaling/CD4 T cells - IL2-060 mono, IL2-060 dimeric"
    NK_CD8_pathname = path_here + "/ckine/data/flow/2019-11-08 monomer IL-2 Fc signaling/NK CD8 T cells - IL2-060 mono, IL2-060 dimeric"

    Tcell_sample = importF(Tcell_pathname, "A")
    NK_CD8_sample = importF(NK_CD8_pathname, "A")

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