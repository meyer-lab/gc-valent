"""
This creates Figure 2.
"""
from .figureCommon import subplotLabel, getSetup, traf_names, Rexpr_names
from ..plot_model_prediction import pstat
from ..model import nParams, getTotalActiveSpecies, runCkineU, getSurfaceGCSpecies, runCkineY0
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import pymc3 as pm, os
from os.path import join
from ..fit_others import build_model



def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    # Blank out for the cartoon
    ax[0].axis('off')

    #subplotLabel(ax[0], 'A')
    #pstat_plot(ax[1])
    #violinPlots(ax[2:6])
    #surf_gc(ax[6], 100.)
    #surf_gc(ax[7], 1000.)
    pretreat(ax[8])
    
    

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
    
    GCexpr = (328. * endo_activeEndo[:, 0]) / (1. + ((kRec_kDeg[:, 0]*(1.-sortF[:, 0])) / (kRec_kDeg[:, 1]*sortF[:, 0]))) # constant according to measured number per cell
    IL7Raexpr = (2591. * endo_activeEndo[:, 0]) / (1. + ((kRec_kDeg[:, 0]*(1.-sortF[:, 0])) / (kRec_kDeg[:, 1]*sortF[:, 0]))) # constant according to measured number per cell
    IL4Raexpr = (254. * endo_activeEndo[:, 0]) / (1. + ((kRec_kDeg[:, 0]*(1.-sortF[:, 0])) / (kRec_kDeg[:, 1]*sortF[:, 0]))) # constant according to measured number per cell

    unkVec = np.zeros((n_params, 500))
    for ii in range (0, 500):
        unkVec[:, ii] = np.array([0., 0., 0., 0., 0., 0., kfwd[ii], 1., 1., 1., 1., 1., 1., k27rev[ii], 1., k33rev[ii], 1., 
            endo_activeEndo[ii, 0], endo_activeEndo[ii, 1], sortF[ii], kRec_kDeg[ii, 0], kRec_kDeg[ii, 1], 0., 0.,
            np.squeeze(GCexpr[ii]), 0., np.squeeze(IL7Raexpr[ii]), 0., np.squeeze(IL4Raexpr[ii]), 0.])
    
    return unkVec, scales

def pstat_calc(unkVec, scales, cytokC):
    ''' This function performs the calculations necessary to produce the Gonnord Figures S3B and S3C. '''
    # import function returns from model.py
    activity = getTotalActiveSpecies().astype(np.float64)
    ts = np.array([10.]) # was 10. in literature
    
    def singleCalc(unkVec, cytokine, conc):
        ''' This function generates the active vector for a given unkVec, cytokine, and concentration. '''
        unkVec = unkVec.copy()
        unkVec[cytokine] = conc

        returnn, retVal = runCkineU(ts, unkVec)

        assert retVal >= 0

        return np.dot(returnn, activity)
    
    assert unkVec.size == nParams()
    actVec_IL7 = np.fromiter((singleCalc(unkVec, 2, x) for x in cytokC), np.float64)
    actVec_IL4 = np.fromiter((singleCalc(unkVec, 4, x) for x in cytokC), np.float64)
    
    actVec = np.concatenate((actVec_IL4 * scales[0], actVec_IL7 * scales[1]))
    return actVec
    
def pstat_plot(ax):
    ''' This function calls the pstat_calc function to re-generate Gonnord figures S3B and S3C with our own fitting data. '''
    PTS = 30
    cytokC_4 = np.array([5., 50., 500., 5000., 50000., 250000.]) / 14900. # 14.9 kDa according to sigma aldrich
    cytokC_7 = np.array([1., 10., 100., 1000., 10000., 100000.]) / 17400. # 17.4 kDa according to prospec bio
    cytokC_common = np.logspace(-3.8, 1.5, num=PTS)
    unkVec, scales = import_samples()
    path = os.path.dirname(os.path.abspath(__file__))
    dataIL4 = pd.read_csv(join(path, "../data/Gonnord_S3B.csv")).values # imports IL4 file into pandas array
    dataIL7 = pd.read_csv(join(path, "../data/Gonnord_S3C.csv")).values # imports IL7 file into pandas array
    
    def plot_structure(IL4vec, IL7vec, title, ax):
        ax.set_title(title)
        ax.plot(np.log10(cytokC_common), IL4vec, color='powderblue', alpha=0.5, label="IL4")
        ax.plot(np.log10(cytokC_common), IL7vec, color='b', alpha=0.5, label='IL7') 
        ax.set_ylabel('pSTAT activation' )
        ax.set_xlabel('cytokine concentration (nM)')
        # ax.legend()

    for ii in range(0,500):
        output = pstat_calc(unkVec[:,ii], scales[ii,:], cytokC_common)
        IL4_output = output[0:PTS]
        IL7_output = output[PTS:(PTS*2)]

        plot_structure(IL4_output, IL7_output, "PBMCs stimulated for 10 min.", ax)
        
    ax.scatter(np.log10(cytokC_4), dataIL4[:,1], color='powderblue', marker='^', edgecolors='k', zorder=100)
    ax.scatter(np.log10(cytokC_4), dataIL4[:,2], color='powderblue', marker='^', edgecolors='k', zorder=200)
    ax.scatter(np.log10(cytokC_7), dataIL7[:,1], color='b', marker='^', edgecolors='k', zorder=300)
    ax.scatter(np.log10(cytokC_7), dataIL7[:,2], color='b', marker='^', edgecolors='k', zorder=400)
        
def violinPlots(ax):
    """ Create violin plots of model posterior. """
    unkVec, scales = import_samples()
    unkVec = unkVec.transpose()
    
    rxn = np.array([unkVec[:, 6], unkVec[:, 13], unkVec[:, 15]])
    rxn = rxn.transpose()
    rxn = pd.DataFrame(rxn)
    traf = pd.DataFrame(unkVec[:, 17:22])
    Rexpr = np.array([unkVec[:, 24], unkVec[:, 26], unkVec[:, 28]])
    Rexpr = Rexpr.transpose()
    Rexpr = pd.DataFrame(Rexpr)
    scales = pd.DataFrame(scales)
    
    rxn.columns = ['kfwd', 'k27rev', 'k33rev']
    a = sns.violinplot(data=np.log10(rxn), ax=ax[0])  # creates names based on dataframe columns
    a.set_xticklabels(a.get_xticklabels(),
                       rotation=40,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))
    
    
    traf.columns = traf_names()
    b = sns.violinplot(data=traf, ax=ax[1])
    b.set_xticklabels(b.get_xticklabels(),
                       rotation=40,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))
    
    Rexpr.columns = ['GCexpr', 'IL7Raexpr', 'IL4Raexpr']
    c = sns.violinplot(data=Rexpr, ax=ax[2])
    c.set_xticklabels(c.get_xticklabels(),
                       rotation=40,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))
    
    scales.columns = ['IL4 scale', 'IL7 scale']
    sns.violinplot(data=scales, ax=ax[3])
    

def pretreat(ax):
    ''' This generates a plot that simulates pretreatment of IL4 or IL7 doses before being stimulated by a standard amount of the other cytokine. The y-axis represent % of inhibition. '''
    path = os.path.dirname(os.path.abspath(__file__))
    data = pd.read_csv(join(path, "../data/Gonnord_S3D.csv")).values 
    IL7_pretreat_conc = data[:, 0] # concentrations used for IL7 pretreatment followed by IL4 stimulation
    IL4_pretreat_conc = data[:, 5] # concentrations used for IL4 pretreatment followed by IL7 stimulation
    print('IL7_pretreat_conc: ' + str(IL7_pretreat_conc))
    IL4_stim_conc = 100. / 14900. # concentration used for IL4 stimulation
    IL7_stim_conc = 50. / 17400. # concentration used for IL7 stimulation
    unkVec, scales = import_samples()
    print('unkVec.shape: ' + str(unkVec.shape))
    print('unkVec[:, 0]: ' + str(unkVec[:, 0]))
    activity = getTotalActiveSpecies().astype(np.float64)
    
    ts = np.array([10.]) # assuming both pretreatment and stimulation were for 10 mins
    
    IL4_act_7pre = np.zeros((5, 500))
    IL7_act_4pre = IL4_act_7pre.copy()
    for x in range(5):
        unkVec2 = unkVec.copy()
        unkVec2[2, :] = IL7_pretreat_conc[x] # plug in IL7 pretreatment concentration
        print('unkVec2[2, 0]: ' + str(unkVec2[2,0]))
        unkVec3 = unkVec.copy()
        unkVec3[4, :] = IL4_stim_conc # plug in IL4 stimulation concentration
        y0 = np.zeros((62, 500))
        yOut = y0.copy()
        for ii in range(500):
            y0[:, ii], retVal = runCkineU(ts, unkVec2[:, ii])
            assert(retVal >= 0)
            yOut[:, ii], retVal = runCkineY0(y0[:, ii], ts, unkVec3[:, ii])
            assert(retVal >= 0)
            IL4_act_7pre[x, ii] = np.dot(yOut, activity)
         
        unkVec4 = unkVec.copy()
        unkVec4[4, :] = IL4_pretreat_conc[x] # plug in IL7 pretreatment concentration
        unkVec5 = unkVec.copy()
        unkVec5[2, :] = IL7_stim_conc # plug in IL4 stimulation concentration
        y0 = np.zeros((62, 500))
        yOut = y0.copy()
        for ii in range(500):
            y0[:, ii], retVal = runCkineU(ts, unkVec4[:, ii])
            assert(retVal >= 0)
            yOut[:, ii], retVal = runCkineY0(y0[:, ii], ts, unkVec5[:, ii])
            assert(retVal >= 0)
            IL7_act_4pre[x, ii] = np.dot(yOut, activity)
        
    print('IL7_act_4pre: ' + str(IL7_act_4pre.shape))
    print('IL4_act_7pre: ' + str(IL4_act_7pre.shape))
    
    
def surf_gc(ax, cytokC_pg):
    size = 40
    ts = np.linspace(0., 100., num=size)
    output = calc_surf_gc(ts, cytokC_pg)
    IL4vec = output[:, 0:size]
    IL7vec = output[:, size:(size*2)]
    
    for ii in range(500):
        ax.plot(ts, IL4vec[ii, :], color='powderblue', label='IL4', alpha=0.5, zorder=ii)
        ax.plot(ts, IL7vec[ii, :], color='b', label='IL7', alpha=0.5, zorder=ii)
    
    ax.set_title(str(cytokC_pg) + ' pg/mL of stimulation')
    ax.set_ylim(0,120)
    ax.set_ylabel("Surface gamma chain (% x 100)")
    ax.set_xlabel("Time (min)")
    
def calc_surf_gc(t, cytokC_pg):
    gc_species_IDX = getSurfaceGCSpecies()
    unkVec, scales = import_samples()
    
    def singleCalc(unkVec, cytokine, conc, t):
        """ Calculates the surface gc over time for one condition. """
        unkVec = unkVec.copy()
        unkVec[cytokine] = conc
        returnn, retVal = runCkineU(t, unkVec)
        assert retVal >= 0
        a = np.dot(returnn, gc_species_IDX)
        return a
    
    size = t.size
    result = np.zeros((500,size*2))
    for ii in range(500):
        # calculate IL4 stimulation
        a = singleCalc(unkVec[:, ii], 4, (cytokC_pg / 14900.), t)
        # calculate IL7 stimulation
        b = singleCalc(unkVec[:, ii], 2, (cytokC_pg / 17400.), t)
        result[ii, :] = np.concatenate((a, b))
        
    return (result / np.max(result)) * 100.


def data_path():
    path = os.path.dirname(os.path.abspath(__file__))
    dataIL4 = pd.read_csv(join(path, "../data/Gonnord_S3B.csv")).values # imports IL4 file into pandas array
    dataIL7 = pd.read_csv(join(path, "../data/Gonnord_S3C.csv")).values
    return (dataIL4, dataIL7)

