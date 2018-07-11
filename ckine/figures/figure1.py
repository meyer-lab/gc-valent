"""
This creates Figure 1.
"""
from .figureCommon import subplotLabel, getSetup
from ..plot_model_prediction import surf_IL2Rb, pstat, surf_gc
from ..model import nParams
import numpy as np
import matplotlib.pyplot as plt
import pymc3 as pm, os
from os.path import join
from ..fit import build_model


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    subplotLabel(ax[0], 'A')
    
    surf_perc(ax[0:4], 'IL2Rb')
    pstat_act(ax[4:6])
    surf_perc(ax[8:12], 'gc')

    f.tight_layout()

    return f

def surf_perc(ax, species):
    if (species == 'IL2Rb'):
        surf = surf_IL2Rb()
    elif (species == 'gc'):
        surf = surf_gc()
    else:
        print('not a valid species')
        return -1

    y_max = 10.
    ts = np.array([0., 2., 5., 15., 30., 60., 90.])
    size = len(ts)
    #unkVec = np.array([0., 0., 0., 0., 2.02E-05, 0.38406199, 0.047815524, 0.118719663, 0.088470489, 0.092582828, 0.107673095, 0.101291622, 0.089522673, 0.023247533, 1.954609999, 0.042469514, 0.300662763, 0.004387724, 2.166630382, 0.673719144, 2.077778334, 0.288524177, 0., 0.])
    unkVec = import_samples()

    def plot_structure(IL2vec, IL15vec, title, ax):
        ax.set_title(title)
        ax.scatter(ts, IL2vec, color='r', label='IL2', alpha=0.7)
        ax.scatter(ts, IL15vec, color='g', label='IL15', alpha=0.7)
        ax.set_ylim(0,(y_max + (0.2 * y_max)))
        ax.set_ylabel("Surface " + str(species) + " (% x " + str(y_max) + ')')
        ax.set_xlabel("Time (min)")
        ax.legend()


    output = surf.calc(unkVec[0]) * y_max
    IL2_1_plus = output[0:(size)]
    IL2_500_plus = output[(size):(size*2)]
    IL2_1_minus = output[(size*2):(size*3)]
    IL2_500_minus = output[(size*3):(size*4)]
    IL15_1_plus = output[(size*4):(size*5)]
    IL15_500_plus = output[(size*5):(size*6)]
    IL15_1_minus = output[(size*6):(size*7)]
    IL15_500_minus = output[(size*7):(size*8)]

    plot_structure(IL2_1_minus, IL15_1_minus, '1 nM and IL2Ra-', ax[0])
    plot_structure(IL2_500_minus, IL15_500_minus, "500 nM and IL2Ra-", ax[1])
    plot_structure(IL2_1_plus, IL15_1_plus, "1 nM and IL2Ra+", ax[2])
    plot_structure(IL2_500_plus, IL15_500_plus, "500 nM and IL2Ra+", ax[3])

    
def pstat_act(ax):
    pstat5 = pstat()
    PTS = 30
    cytokC = np.logspace(-3.3, 2.7, PTS)
    y_max = 100.
    #unkVec = np.array([0., 0., 0., 0., 2.02E-05, 0.38406199, 0.047815524, 0.118719663, 0.088470489, 0.092582828, 0.107673095, 0.101291622, 0.089522673, 0.023247533, 1.954609999, 0.042469514, 0.300662763, 0.004387724, 2.166630382, 0.673719144, 2.077778334, 0.288524177, 0., 0.])
    unkVec = import_samples()
    
    def plot_structure(IL2vec, IL15vec, title, ax):
        ax.set_title(title)
        ax.scatter(np.log10(cytokC), IL2vec, color='r', alpha=0.5, label="IL2")
        ax.scatter(np.log10(cytokC), IL15vec, color='g', alpha=0.5, label='IL15')
        ax.set_ylim(0,(y_max + (0.25*y_max)))
        ax.set_ylabel('Maximal p-STAT5 (% x ' + str(y_max) + ')')
        ax.set_xlabel('log10 of cytokine concentration (nM)')
        ax.legend()

    output = pstat5.calc(unkVec[0]) * y_max
    IL2_plus = output[0:PTS]
    IL2_minus = output[PTS:(PTS*2)]
    IL15_plus = output[(PTS*2):(PTS*3)]
    IL15_minus = output[(PTS*3):(PTS*4)]

    plot_structure(IL2_minus, IL15_minus, "IL2Ra- YT-1 cells", ax[0])
    plot_structure(IL2_plus, IL15_plus, "IL2Ra+ YT-1 cells", ax[1])
    
def import_samples():
    bmodel = build_model()

    path = os.path.dirname(os.path.abspath(__file__))
    trace = pm.backends.text.load(join(path, '../../IL2_model_results'), bmodel.M)
    kfwd = trace.get_values('kfwd', chains=[0])
    k4rev = trace.get_values('rxn__0', chains=[0])
    k5rev = trace.get_values('rxn__1', chains=[0])
    k16rev = trace.get_values('rxn__2', chains=[0])
    k17rev = trace.get_values('rxn__3', chains=[0])
    k22rev = trace.get_values('rxn__4', chains=[0])
    k23rev = trace.get_values('rxn__5', chains=[0])
    k27rev = trace.get_values('rxn__6', chains=[0])
    k31rev = trace.get_values('rxn__7', chains=[0])
    endo = trace.get_values('endo__0', chains=[0])
    activeEndo = trace.get_values('endo__1', chains=[0])
    sortF = trace.get_values('sortF', chains=[0])
    kRec = trace.get_values('kRec_kDeg__0', chains=[0])
    kDeg = trace.get_values('kRec_kDeg__1', chains=[0])
    IL2Ra_expr = trace.get_values('IL2Raexpr__0', chains=[0])
    IL2Rb_expr = trace.get_values('IL2Raexpr__1', chains=[0])
    gc_expr = trace.get_values('IL2Raexpr__2', chains=[0])
    IL15Ra_expr = trace.get_values('IL2Raexpr__3', chains=[0])
    
    # probably can't use numpy for this
    unkVec = np.zeros((nParams()), len(kfwd))
    for ii in range (0, len(kfwd)):
        unkVec = np.array([0., 0., 0., 0., kfwd[ii], k4rev[ii], k5rev[ii], k16rev[ii], k17rev[ii], k22rev[ii], k23rev[ii], k27rev[ii], k31rev[ii], endo[ii], activeEndo[ii], sortF[ii], kRec[ii], kDeg[ii], IL2Ra_expr[ii], IL2Rb_expr[ii], gc_expr[ii], IL15Ra_expr[ii], 0., 0.])
    
    return unkVec
        
        
    
    
    
    

