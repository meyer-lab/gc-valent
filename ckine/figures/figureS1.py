"""
This creates Figure S1.
"""
import string
from .figureCommon import subplotLabel, getSetup, kfwd_info
from .figure1 import pstat_act, violinPlots, rateComp
from ..imports import import_samples_2_15

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((8, 6), (2, 3), mults=[1], multz={1: 2}, empts=[3])

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    unkVec, scales = import_samples_2_15(Traf=False, N=100)
    full_unkVec, full_scales = import_samples_2_15(Traf=False)
    kfwd_avg, kfwd_std = kfwd_info(full_unkVec)
    print("kfwd = " + str(kfwd_avg) + " +/- " + str(kfwd_std))
    pstat_act(ax[0], unkVec, scales)
    rateComp(ax[1], full_unkVec)
    violinPlots(ax[2:4], full_unkVec, full_scales, Traf=False)

    f.tight_layout()

    return f
