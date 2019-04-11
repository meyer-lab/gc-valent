"""
This creates Figure S5. Full panel of Geweke convergence tests.
"""
import string
import os
import pymc3 as pm
import matplotlib.cm as cm
import numpy as np
import os
from os.path import join
from .figureCommon import subplotLabel, getSetup
from ..fit import build_model as build_model_2_15


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 2))

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    plot_geweke(ax[0:3], True)

    f.tight_layout()

    return f

def plot_geweke(ax, traf):
    """ Uses geweke criterion to evaluate model convergence during fitting. """
    bmodel = build_model_2_15() # access pymc3 object from IL2_IL15 fitting
    path = os.path.dirname(os.path.abspath(__file__))
    if traf:
        trace = pm.backends.text.load(join(path, '../../IL2_model_results'), bmodel.M)
    else:
        trace = pm.backends.text.load(join(path, '../../IL2_15_no_traf'), bmodel.M)

    # use use trace to calculate geweke z-scores
    score = pm.diagnostics.geweke(trace, first=0.1, last=0.5, intervals=20)
    
    # plot the scores for rxn rates
    rxn_len = len(score[0]['rxn'])
    rxn_names = [r'$k_{4}$', r'$k_{5}$', r'$k_{16}$', r'$k_{17}$', r'$k_{22}$', r'$k_{23}$']
    colors = cm.rainbow(np.linspace(0, 1, rxn_len))
    for ii in range(rxn_len):
        ax[0].scatter(score[0]['rxn'][ii][:,0], score[0]['rxn'][ii][:,1], marker='o', s=25, color=colors[ii], label=rxn_names[ii]) # plot all rates within rxn
    ax[0].axhline(-1., c='r')
    ax[0].axhline(1., c='r')
    ax[0].set_ylim(-1.25,1.25)
    ax[0].set_xlim(0-10,.5*trace['rxn'].shape[0]/2+10)
    ax[0].set_title('Geweke Plot Comparing first 10% and Slices of the Last 50% of Chain\nDifference in Mean krev Z scores')
    ax[0].legend()
    ax[0].set_xlabel("Position in Chain")
    ax[0].set_ylabel("Geweke Score")

    # plot the scores for receptor expression rates
    rexpr_len = len(score[0]['IL2Raexpr'])
    rexpr_names = ['IL-2Rα', 'IL-2Rβ', r'$\gamma_{c}$', 'IL-15Rα']
    colors = cm.rainbow(np.linspace(0, 1, rexpr_len))
    for ii in range(rexpr_len):
        ax[1].scatter(score[0]['IL2Raexpr'][ii][:,0], score[0]['IL2Raexpr'][ii][:,1], marker='o', s=25, color=colors[ii], label=rexpr_names[ii])
    ax[1].axhline(-1., c='r')
    ax[1].axhline(1., c='r')
    ax[1].set_ylim(-1.25,1.25)
    ax[1].set_xlim(0-10,.5*trace['IL2Raexpr'].shape[0]/2+10)
    ax[1].set_title('Geweke Plot Comparing first 10% and Slices of the Last 50% of Chain\nDifference in Mean Rexpr Z scores')
    ax[1].legend()
    ax[1].set_xlabel("Position in Chain")
    ax[1].set_ylabel("Geweke Score")

    # plot the scores for scaling constant and kfwd
    ax[2].scatter(score[0]['scales'][:,0], score[0]['scales'][:,1], marker='o', s=25, color='g', label=r'$C_{5}$')
    ax[2].scatter(score[0]['kfwd'][:,0], score[0]['kfwd'][:,1], marker='o', s=25, color='b', label=r'$k_{fwd}$')
