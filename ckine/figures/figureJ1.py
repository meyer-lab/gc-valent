"""
This creates Figure 4 for Single Cell data analysis. Plots of flow intensity versus receptor quantification.
"""
import os
import numpy as np
from .figureCommon import subplotLabel, getSetup

path_here_jl = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import Pkg
Pkg.activate(path_here_jl + "/gcSolver.jl")
from julia import gcSolver


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((10, 5), (2, 4))

    subplotLabel(ax)
    a = gcSolver.getUnkVec()
    print(a)
    print(type(a))

    return f
