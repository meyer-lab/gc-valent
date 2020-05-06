"""
This file contains functions that are used in multiple figures.
"""
from julia import gcSolver
from julia.api import Julia
jl = Julia(compiled_modules=False)


def getUnkVecPy():
    return gcSolver.getUnkVec()
