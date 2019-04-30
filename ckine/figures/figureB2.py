"""
This creates Figure 2.
"""

import string
import numpy as np
import matplotlib.pyplot as plt
from .figureCommon import subplotLabel, getSetup
from ..model import runIL2simple


def makeFigure():
    # Get list of axis objects
    ax, f = getSetup((8, 11), (3, 3))
    
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
    
    for i in range(3):
        plot_dResp_2Ra(CD25[i], ax[i])
        ax[i].set_title(titles[i])
        plot_dResp_2Rb(CD25[i], ax[3+i])
        ax[3+i].set_title(titles[i])
        plot_dResp_2Rb_HIGH(CD25[i], ax[6+i])
        ax[6+i].set_title(titles[i])
    
    ax[0].legend(title="IL2Ra Kd vs wt")
    ax[3].legend(title="IL2Rb Kd vs wt")
    ax[6].legend(title="IL2Rb Kd vs wt")

    f.tight_layout()

    return f

changesAff = np.logspace(-2, 2, num=7)
CD25 = [1.0, 0.1, 0.0]
titles = ["CD25+", "10% CD25+", "CD25-"]


def dRespon_loc(input, CD25): # same as dRespon except with different ILs range
        """ Calculate an IL2 dose response curve. """
        ILs = np.logspace(-4.0, 3.0)
        activee = np.array([runIL2simple(input, ii, CD25) for ii in ILs])
        return ILs, activee

def plot_dResp_2Ra(CD25, ax):
    """ Plots dose response curves for various IL2Ra affinities given a CD25 relative expression rate. """
    for i, itemA in enumerate(changesAff):
        ILs, BB = dRespon_loc([itemA, 1.0, 5.0], CD25)
        ax.semilogx(ILs, BB, label=str(round(itemA, 2)))

    ax.set_ylabel('Active Receptor Complexes (#/cell)');
    ax.set_xlabel('IL2 [nM]');
    
def plot_dResp_2Rb(CD25, ax):
    """ Plots dose response curves for various IL2Rb affinities given a CD25 relative expression rate. """
    for i, itemB in enumerate(changesAff):
        ILs, BB = dRespon_loc([1.0, itemB, 5.0], CD25)
        ax.semilogx(ILs, BB, label=str(round(itemB, 2)))

    ax.set_ylabel('Active Receptor Complexes (#/cell)');
    ax.set_xlabel('IL2 [nM]');

def plot_dResp_2Rb_HIGH(CD25, ax):
    """ Plots dose response curves for various IL2Rb affinities given a CD25 relative expression rate. """
    for i, itemB in enumerate(changesAff):
        ILs, BB = dRespon_loc([0.1, itemB, 5.0], CD25)
        ax.semilogx(ILs, BB, label=str(round(itemB, 2)))

    ax.set_ylabel('Active Receptor Complexes (#/cell)');
    ax.set_xlabel('IL2 [nM]');
