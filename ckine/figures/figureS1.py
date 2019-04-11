"""
This creates Figure S1.
"""
import string
import os
import pymc3 as pm
import os
from os.path import join
from .figureCommon import subplotLabel, getSetup, import_samples_2_15, kfwd_info
from .figure1 import pstat_act, violinPlots, rateComp
from ..fit import build_model as build_model_2_15


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((8, 6), (2, 3), mults=[1], multz={1: 2}, empts=[3])

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    '''unkVec, scales = import_samples_2_15(Fig1=False)
    kfwd_avg, kfwd_std = kfwd_info(unkVec)
    print("kfwd = " + str(kfwd_avg) + " +/- " + str(kfwd_std))
    pstat_act(ax[0], unkVec, scales)
    rateComp(ax[1], unkVec)
    violinPlots(ax[2:4], unkVec, scales, Fig1=False)
    '''
    plot_geweke(ax[0], True)

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

    # use use trace to calculate geweke z-scores ... TODO: figure out proper arguments for first, last, intervals
    out = pm.diagnostics.geweke(trace, first=0.1, last=0.5, intervals=20)
    print(len(out[0]))
