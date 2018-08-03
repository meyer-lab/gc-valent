"""
This creates Figure 4.
"""
import string
import numpy as np
import os
import pandas as pds
from .figureCommon import subplotLabel, getSetup, load_cells


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    data, cell_names = load_cells()

    f.tight_layout()

    return f

