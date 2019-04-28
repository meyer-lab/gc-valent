"""
This creates Figure 1.
"""
import string
import numpy as np
import matplotlib.pyplot as plt
from .figureCommon import subplotLabel, getSetup
from ..model import runIL2simple

def makeFigure():
    # Get list of axis objects
    ax, f = getSetup((11, 6), (2, 2))

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
        
    subplotA(ax[0])
    subplotB(ax[1])
    subplotC(ax[2])
    subplotD(ax[3])

    f.tight_layout()

    return f

def dRespon(input, CD25=1.0):
    """ Calculate an IL2 dose response curve. """
    ILs = np.logspace(-3.0, 3.0)
    activee = np.array([runIL2simple(input, ii, CD25) for ii in ILs])
        
    return ILs, activee

def IC50global(x, y):
    """ Calculate half-maximal concentration w.r.t. wt. """
    return np.interp(20.0, y, x)

changesA = np.logspace(-2, 1, num=20)
changesB = np.array([1.0, 2.0, 5.0, 10.0, 20.0])
output = np.zeros((changesA.size, changesB.size))

for i, itemA in enumerate(changesA):
    for j, itemB in enumerate(changesB):
        ILs, BB = dRespon([itemA, itemB, 5.0])
        output[i, j] = IC50global(ILs, BB)
        
def subplotA(ax):
    for ii in range(output.shape[1]):
        ax.loglog(changesA, output[:, ii], label=str(changesB[ii]));
    ax.loglog([0.01, 10.], [0.17, 0.17], 'k-')
    ax.set(ylabel='Half-Maximal IL2 Concentration [nM]', xlabel='IL2Ra-IL2 Kd (fold wt)', ylim=(0.01, 20))
    ax.legend(title="IL2Rb Kd v wt");
    
def subplotB(ax):
    wt = dRespon([1.0, 1.0, 5.0])
    ten = dRespon([0.1, 5.0, 5.0])

    ax.semilogx(wt[0], wt[1], label="wt");
    ax.semilogx(ten[0], ten[1], 'r', label="10X higher/lower affinity IL2Ra/IL2Rb")
    ax.set(ylabel='Active Receptor Complexes (#/cell)', xlabel='IL2 [nM]')
    ax.legend();

def subplotC(ax):
    changesA = np.logspace(-1, 1.5, num=20)
    changesB = np.array([0.0, 0.1, 0.25, 0.5, 1.0])
    output = np.zeros((changesA.size, changesB.size))

    for i, itemA in enumerate(changesA):
        for j, itemB in enumerate(changesB):
            ILs, BB = dRespon([1.0, itemA, 5.0], CD25=itemB)
            output[i, j] = IC50global(ILs, BB)
    
    for ii in range(output.shape[1]):
        ax.loglog(changesA, output[:, ii], label=str(changesB[ii]));

    ax.loglog([0.1, 10.], [0.17, 0.17], 'k-')
    ax.set(ylabel='Half-Maximal IL2 Concentration [nM]',xlabel='IL2Rb-IL2 Kd (relative to wt)',ylim=(0.001, 10),xlim=(0.1, 10))
    ax.legend(title="CD25 rel expr");
    
def subplotD(ax):
    for i, itemA in enumerate(changesA):
        for j, itemB in enumerate(changesB):
            ILs, BB = dRespon([0.1, itemA, 5.0], CD25=itemB)
            output[i, j] = IC50global(ILs, BB)
    
    for ii in range(output.shape[1]):
        ax.loglog(changesA, output[:, ii], label=str(changesB[ii]));

    ax.loglog([0.1, 10.], [0.17, 0.17], 'k-')
    ax.set(ylabel='Half-Maximal IL2 Concentration [nM]',xlabel='IL2Rb-IL2 Kd (relative to wt)',ylim=(0.001, 10),xlim=(0.1, 10))
    ax.legend(title="CD25 rel expr");
