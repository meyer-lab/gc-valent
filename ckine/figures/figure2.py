"""
This creates Figure 2.
"""
from .figureCommon import subplotLabel, getSetup, rev_rxn_names, traf_names, Rexpr_names
from ..plot_model_prediction import pstat
from ..model import nParams
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import pymc3 as pm, os
from os.path import join
from ..fit import build_model


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    # Blank out for the cartoon
    ax[0].axis('off')

    subplotLabel(ax[0], 'A')

    f.tight_layout()

    return f

def import_samples():
    bmodel = build_model()
    n_params = nParams()

    path = os.path.dirname(os.path.abspath(__file__))
    trace = pm.backends.text.load(join(path, '../../IL4-7_model_results'), bmodel.M)
    kfwd = trace.get_values('kfwd', chains=[0])
    k27rev = trace.get_values('k27rev', chains=[0])
    k33rev = trace.get_values('k33rev', chains=[0])
    endo_activeEndo = trace.get_values('endo', chains=[0])
    sortF = trace.get_values('sortF', chains=[0])
    kRec_kDeg = trace.get_values('kRec_kDeg', chains=[0])
    scales = trace.get_values('scales', chains=[0])
    
    GCexpr = (328. * endo_activeEndo[:, 0]) / (1 + (kRec_kDeg[:, 0] / kRec_kDeg[:, 1]))
    IL7Raexpr = (2591. * endo_activeEndo[:, 0]) / (1 + (kRec_kDeg[:, 0] / kRec_kDeg[:, 1]))
    IL4Raexpr = (254. * endo_activeEndo[:, 0]) / (1 + (kRec_kDeg[:, 0] / kRec_kDeg[:, 1]))
    
    
    
    unkVec = np.zeros((n_params, 500))
    for ii in range (0, 500):
        unkVec[:, ii] = np.array([0., 0., 0., 0., 0., 0., kfwd[ii], 0., 0., 0., 0., 0., 0., k27rev[ii], 0., k33rev[ii], 0., endo_activeEndo[ii, 0], endo_activeEndo[ii, 1], sortF[ii], kRec_kDeg[ii, 0], kRec_kDeg[ii, 1], 0., 0., GCexpr[ii], 0., IL7Raexpr[ii], 0., IL4Raexpr[ii], 0.])
    
    return unkVec
