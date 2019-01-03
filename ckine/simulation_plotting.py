"""File to simulate model with different drug combinations, used for plotting."""
import numpy as np
from .model import nSpecies, getTotalActiveCytokine
from .tensor_generation import ySolver

def input_IL2_15(final_conc, num):
    '''Function that creates the input for the solver. Takes in 1nM or 500nM for final_conc.'''
    Rexpr = np.array([3.8704, 0.734, 1.7147, 0.32010875, 0.0, 0.0, 0.0, 0.0]) #Il2ra, Il2rb, Il2rg, Il15ra, Il7r, Il9r, IL4Ra, IL21Ra in that order
    ligand_conc = np.zeros((num, 6)) # Cytokine stimulation concentrations in the following order: IL2, 15, 7, 9, 4, 21, and in nM.
    IL2 = np.linspace(final_conc, final_conc/num, num)
    IL15 = np.linspace(final_conc/num, final_conc, num)
    for ii in range(num):
        ligand_conc[ii, 0:2] = [IL2[ii], IL15[ii]]
    return ligand_conc, Rexpr, np.log(IL2/IL15)

def solve_IL2_IL15(final_conc, num, time, nSpecies = nSpecies):
    """Function to simulate model with IL2 and IL15 only at timepoint tps."""
    ligand_conc, Rexpr, xaxis = input_IL2_15(final_conc, num)
    tps = np.array([time])
    yOut = np.zeros((num, nSpecies()))
    active = np.zeros((num))
    for ii in range(num):
        yOut[ii] = ySolver(np.concatenate((ligand_conc[ii], Rexpr)), tps) 
        active[ii] = getTotalActiveCytokine(0, np.squeeze(yOut[ii])) + getTotalActiveCytokine(1, np.squeeze(yOut[ii]))
    return active, num, xaxis
