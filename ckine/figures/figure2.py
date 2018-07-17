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
    ''' Imports the csv files into a proper unkVec'''
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
    
    return unkVec, scales

def pstat_calc(ax):
    ''' This function generates plots that mimic the Gonnard Figures S3B and S3C based on the fitting results. '''
    PTS = 30
    cytokC_4 = np.linspace(5./14900., 250000./14900., num=PTS)
    cytokC_7 = np.linspace(1./17400., 100000./17400., num=PTS)

    # import function returns from model.py
    activity = getTotalActiveSpecies().astype(np.float64)
    ts = np.array([10.]) # was 10. in literature
    # percentage value that is used in scaling output
    y_max = 100
    unkVec, scales = import_samples()
    
    def singleCalc(unkVec, cytokine, conc):
        ''' This function generates the active vector for a given unkVec, cytokine, and concentration. '''
        unkVec = unkVec.copy()
        unkVec[cytokine] = conc

        returnn, retVal = runCkineU(ts, unkVec)

        assert retVal >= 0

        return np.dot(returnn, activity)
    
    assert unkVec.size == nParams()
    actVec_IL7 = np.fromiter((singleCalc(unkVec, 2, x) for x in cytokC_7), np.float64)
    actVec_IL4 = np.fromiter((singleCalc(unkVec, 4, x) for x in cytokC_4), np.float64)
    
    actVec = np.concatenate((actVec_IL4 * scales[0], actVec_IL7 * scales[1]))
    

class pstat:
    '''Generate values to match the pSTAT5 measurements used in fitting'''
    def __init__(self):
        self.PTS = 30
        self.cytokC_4 = np.linspace(5./14900., 250000./14900., num=self.PTS)
        self.cytokC_7 = np.linspace(1./17400., 100000./17400., num=self.PTS)
        
        # import function returns from model.py
        self.activity = getTotalActiveSpecies().astype(np.float64)
        
        self.ts = np.array([10.]) # was 10. in literature

        # percentage value that is used in scaling output
        self.y_max = 100

    def singleCalc(self, unkVec, cytokine, conc):
        """ Calculates the activity for one condition. """
        unkVec = unkVec.copy()
        unkVec[cytokine] = conc

        returnn, retVal = runCkineU(self.ts, unkVec)

        assert retVal >= 0

        return np.dot(returnn, self.activity)

    def calc(self, unkVec):
        '''This function uses an unkVec that has the same elements as the unkVec in fit.py'''
        assert unkVec.size == nParams()

        unkVec_IL2Raminus = unkVec.copy()
        unkVec_IL2Raminus[18] = 0.0 # set IL2Ra expression rate to 0

        # Calculate activities
        actVec_IL2 = np.fromiter((self.singleCalc(unkVec, 0, x) for x in self.cytokC), np.float64)
        actVec_IL2_IL2Raminus = np.fromiter((self.singleCalc(unkVec_IL2Raminus, 0, x) for x in self.cytokC), np.float64)
        actVec_IL15 = np.fromiter((self.singleCalc(unkVec, 1, x) for x in self.cytokC), np.float64)
        actVec_IL15_IL2Raminus = np.fromiter((self.singleCalc(unkVec_IL2Raminus, 1, x) for x in self.cytokC), np.float64)

        # Normalize to the maximal activity, put together into one vector
        actVec = np.concatenate((actVec_IL2 / np.max(actVec_IL2), actVec_IL2_IL2Raminus / np.max(actVec_IL2_IL2Raminus), actVec_IL15 / np.max(actVec_IL15), actVec_IL15_IL2Raminus / np.max(actVec_IL15_IL2Raminus)))

        return actVec

