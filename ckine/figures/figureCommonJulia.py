"""
This file contains functions that are used in multiple figures.
"""
from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import gcSolver

def getUnkVecPy():
    return gcSolver.getUnkVec()
