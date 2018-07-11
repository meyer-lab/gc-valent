"""
This creates Figure 1.
"""
from .figureCommon import subplotLabel, getSetup
from ..plot_model_prediction import surf_IL2Rb, pstat, surf_gc
import numpy as np
import matplotlib.pyplot as plt


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
    unkVec = np.array([0., 0., 0., 0., 2.02E-05, 0.38406199, 0.047815524, 0.118719663, 0.088470489, 0.092582828, 0.107673095, 0.101291622, 0.089522673, 0.023247533, 1.954609999, 0.042469514, 0.300662763, 0.004387724, 2.166630382, 0.673719144, 2.077778334, 0.288524177, 0., 0.])

    def plot_structure(IL2vec, IL15vec, title, ax):
        ax.set_title(title)
        ax.scatter(ts, IL2vec, color='r', label='IL2', alpha=0.7)
        ax.scatter(ts, IL15vec, color='g', label='IL15', alpha=0.7)
        ax.set_ylim(0,(y_max + (0.2 * y_max)))
        ax.set_ylabel("Surface " + str(species) + " (% x " + str(y_max) + ')')
        ax.set_xlabel("Time (min)")
        ax.legend()


    output = surf.calc(unkVec) * y_max
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
    unkVec = np.array([0., 0., 0., 0., 2.02E-05, 0.38406199, 0.047815524, 0.118719663, 0.088470489, 0.092582828, 0.107673095, 0.101291622, 0.089522673, 0.023247533, 1.954609999, 0.042469514, 0.300662763, 0.004387724, 2.166630382, 0.673719144, 2.077778334, 0.288524177, 0., 0.])
    
    def plot_structure(IL2vec, IL15vec, title, ax):
        ax.set_title(title)
        ax.scatter(np.log10(cytokC), IL2vec, color='r', alpha=0.5, label="IL2")
        ax.scatter(np.log10(cytokC), IL15vec, color='g', alpha=0.5, label='IL15')
        ax.set_ylim(0,(y_max + (0.25*y_max)))
        ax.set_ylabel('Maximal p-STAT5 (% x ' + str(y_max) + ')')
        ax.set_xlabel('log10 of cytokine concentration (nM)')
        ax.legend()

    output = pstat5.calc(unkVec) * y_max
    IL2_plus = output[0:PTS]
    IL2_minus = output[PTS:(PTS*2)]
    IL15_plus = output[(PTS*2):(PTS*3)]
    IL15_minus = output[(PTS*3):(PTS*4)]

    plot_structure(IL2_minus, IL15_minus, "IL2Ra- YT-1 cells", ax[0])
    plot_structure(IL2_plus, IL15_plus, "IL2Ra+ YT-1 cells", ax[1])


def gc_surf_perc(ax):
    surfGC = surf_gc()

    y_max = 10.
    ts = np.array([0., 2., 5., 15., 30., 60., 90.])
    size = len(ts)
    unkVec = np.array([0., 0., 0., 0., 2.02E-05, 0.38406199, 0.047815524, 0.118719663, 0.088470489, 0.092582828, 0.107673095, 0.101291622, 0.089522673, 0.023247533, 1.954609999, 0.042469514, 0.300662763, 0.004387724, 2.166630382, 0.673719144, 2.077778334, 0.288524177, 0., 0.])

    def plot_structure(IL2vec, IL15vec, title, ax):
        ax.set_title(title)
        ax.scatter(ts, IL2vec, color='r', label='IL2', alpha=0.7)
        ax.scatter(ts, IL15vec, color='g', label='IL15', alpha=0.7)
        ax.set_ylim(0,(y_max + (0.2 * y_max)))
        ax.set_ylabel("Surface gc (% x " + str(y_max) + ')')
        ax.set_xlabel("Time (min)")
        ax.legend()


    output = surfGC.calc(unkVec) * y_max
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
