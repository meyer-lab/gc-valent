"""
This creates Figure 3.
"""

import string
import numpy as np
from .figureCommon import subplotLabel, getSetup
from .figureB1 import runIL2simple
from ..model import receptor_expression
from ..make_tensor import rxntfR
from ..imports import import_Rexpr

df, _, _ = import_Rexpr()
df.reset_index(inplace=True)

def makeFigure():
    """ Get a list of the axis objects and create a figure. """
    ax, f = getSetup((8, 11), (3, 3))

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
        
    cellReceptors = np.zeros((3,3))
    
    for i, cellName in enumerate(cellNames):
        IL2Ra = df.loc[(df["Cell Type"] == cellName) & (df["Receptor"] =='IL-2R$\\alpha$'), "Count"].item()
        IL2Rb = df.loc[(df["Cell Type"] == cellName) & (df["Receptor"] =='IL-2R$\\beta$'), "Count"].item()
        gc = df.loc[(df["Cell Type"] == cellName) & (df["Receptor"] =='$\\gamma_{c}$'), "Count"].item()
        cellReceptors[i, :] = receptor_expression(np.array([IL2Ra, IL2Rb, gc]).astype(np.float), rxntfR[17], rxntfR[20], rxntfR[19], rxntfR[21])

    for i in range(3):
        plot_lDeg_2Ra(ax[i], CD25_input[i], cellReceptors[i, :], adj_receptors=True)
        ax[i].set_title(cellNames[i])
        plot_lDeg_2Rb(ax[3 + i], CD25_input[i], cellReceptors[i, :], adj_receptors=True)
        ax[3 + i].set_title(cellNames[i])
        plot_lDeg_2Rb_HIGH(ax[6 + i], CD25_input[i], cellReceptors[i, :], adj_receptors=True)
        ax[6 + i].set_title(cellNames[i])

    ax[0].legend(title="IL2Ra Kd vs wt")
    ax[3].legend(title="IL2Rb Kd vs wt")
    ax[6].legend(title="IL2Rb Kd vs wt")

    return f


changesAff = np.logspace(-2, 2, num=7)
CD25_input = [1.0, 0.1, 0.0]
cellNames = ["T-reg", "T-helper", "NK"]


def ligandDeg_IL2(input_params, CD25, input_receptors=None, adj_receptors=False):
    """ Calculate an IL2 degradation curve. """
    ILs = np.logspace(-4.0, 5.0)
    ld = np.array([runIL2simple(rxntfR, input_params, ii, CD25=CD25, input_receptors=input_receptors, adj_receptors=adj_receptors, ligandDegradation=True) for ii in ILs])
    return ILs, ld


def plot_lDeg_2Ra(ax, CD25, input_receptors=None, adj_receptors=False):
    """ Plots IL2 degradation curves for various IL2Ra affinities given a CD25 relative expression rate. """
    for _, itemA in enumerate(changesAff):
        ILs, BB = ligandDeg_IL2([itemA, 1.0, 5.0], CD25, input_receptors, adj_receptors)
        ax.semilogx(ILs, BB, label=str(round(itemA, 2)))

    ax.set_ylabel('Rate of IL2 Degradation')
    ax.set_xlabel('IL2 [nM]')


def plot_lDeg_2Rb(ax, CD25, input_receptors=None, adj_receptors=False):
    """ Plots IL2 degradation curves for various IL2Rb affinities given a CD25 relative expression rate. """
    for _, itemB in enumerate(changesAff):
        ILs, BB = ligandDeg_IL2([1.0, itemB, 5.0], CD25, input_receptors, adj_receptors)
        ax.semilogx(ILs, BB, label=str(round(itemB, 2)))

    ax.set_ylabel('Rate of IL2 Degradation')
    ax.set_xlabel('IL2 [nM]')


def plot_lDeg_2Rb_HIGH(ax, CD25, input_receptors=None, adj_receptors=False):
    """ Plots IL2 degradation curves for various IL2Rb affinities given a CD25 relative expression rate. """
    for _, itemB in enumerate(changesAff):
        ILs, BB = ligandDeg_IL2([0.1, itemB, 5.0], CD25, input_receptors, adj_receptors)
        ax.semilogx(ILs, BB, label=str(round(itemB, 2)))

    ax.set_ylabel('Rate of IL2 Degradation')
    ax.set_xlabel('IL2 [nM]')

