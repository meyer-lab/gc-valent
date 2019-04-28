"""
This creates Figure 3.
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
        
    plot_lDeg_2Ra(1.0, ax[0]);
    ax[0].set_title("CD25+")
    plot_lDeg_2Ra(0.1, ax[1]);
    ax[1].set_title("10% CD25+")
    plot_lDeg_2Ra(0.0, ax[2]);
    ax[2].set_title("CD25-");
    ax[0].legend(title="IL2Ra Kd vs wt");
    
    plot_lDeg_2Rb(1.0, ax[3]);
    ax[3].set_title("CD25+")
    plot_lDeg_2Rb(0.1, ax[4]);
    ax[4].set_title("10% CD25+")
    plot_lDeg_2Rb(0.0, ax[5]);
    ax[5].set_title("CD25-");
    ax[3].legend(title="IL2Rb Kd vs wt");
    
    plot_lDeg_2Rb_HIGH(1.0, ax[6]);
    ax[6].set_title("CD25+")
    plot_lDeg_2Rb_HIGH(0.1, ax[7]);
    ax[7].set_title("10% CD25+")
    plot_lDeg_2Rb_HIGH(0.0, ax[8]);
    ax[8].set_title("CD25-");
    ax[6].legend(title="IL2Rb Kd vs wt");

    f.tight_layout()

    return f

changesAff = np.logspace(-2, 2, num=7)

def ligandDeg_IL2(input, CD25): 
        """ Calculate an IL2 degradation curve. """
        ILs = np.logspace(-4.0, 5.0)
        ld = np.array([runIL2simple(input, ii, CD25, True) for ii in ILs])
        return ILs, ld

def plot_lDeg_2Ra(CD25, ax):
    """ Plots IL2 degradation curves for various IL2Ra affinities given a CD25 relative expression rate. """
    for i, itemA in enumerate(changesAff):
        ILs, BB = ligandDeg_IL2([itemA, 1.0, 5.0], CD25)
        ax.semilogx(ILs, BB, label=str(round(itemA, 2)))

    ax.set_ylabel('Rate of IL2 Degradation');
    ax.set_xlabel('IL2 [nM]');
    
def plot_lDeg_2Rb(CD25, ax):
    """ Plots IL2 degradation curves for various IL2Rb affinities given a CD25 relative expression rate. """
    for i, itemB in enumerate(changesAff):
        ILs, BB = ligandDeg_IL2([1.0, itemB, 5.0], CD25)
        ax.semilogx(ILs, BB, label=str(round(itemB, 2)))

    ax.set_ylabel('Rate of IL2 Degradation');
    ax.set_xlabel('IL2 [nM]');

def plot_lDeg_2Rb_HIGH(CD25, ax):
    """ Plots IL2 degradation curves for various IL2Rb affinities given a CD25 relative expression rate. """
    for i, itemB in enumerate(changesAff):
        ILs, BB = ligandDeg_IL2([0.1, itemB, 5.0], CD25)
        ax.semilogx(ILs, BB, label=str(round(itemB, 2)))

    ax.set_ylabel('Rate of IL2 Degradation');
    ax.set_xlabel('IL2 [nM]');
