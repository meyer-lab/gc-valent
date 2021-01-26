"""
This creates Figure 7, tensor factorization of mutant and WT biv and monovalent ligands.
"""

import os
import numpy as np
import pandas as pd
from .figureCommon import subplotLabel, getSetup
from ..MBmodel import getBindConvDict


path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)

    getBindConvDict()

    return f
