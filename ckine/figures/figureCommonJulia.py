"""
This file contains functions that are used in multiple figures from Julia.
"""
from julia import gcSolver
from julia.api import Julia
jl = Julia(compiled_modules=False)


def getUnkVecPy():
    "Returns initial points for unknown quantities"
    return gcSolver.getUnkVec()
