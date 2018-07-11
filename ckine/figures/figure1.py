"""
This creates Figure 1.
"""
from .figureCommon import subplotLabel, getSetup
from ..plot_model_prediction import surf_IL2Rb, pstat, surf_gc
import numpy as np


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    subplotLabel(ax[0], 'A')
    
    IL2Rb_surf_perc(ax[0:4])

    f.tight_layout()

    return f

def IL2Rb_surf_perc(ax):
    surf2Rb = surf_IL2Rb()
    
    y_max = 10.
    ts = np.array([0., 2., 5., 15., 30., 60., 90.])
    unkVec = np.array([0., 0., 0., 0., 2.02E-05, 0.38406199, 0.047815524, 0.118719663, 0.088470489, 0.092582828, 0.107673095, 0.101291622, 0.089522673, 0.023247533, 1.954609999, 0.042469514, 0.300662763, 0.004387724, 2.166630382, 0.673719144, 2.077778334, 0.288524177, 0., 0.])
    
    def plot_structure(IL2vec, IL15vec, title, ax):
        plt.title(title)
        plt.scatter(ts, IL2vec, color='r', label='IL2', alpha=0.7, ax=ax)
        plt.scatter(ts, IL15vec, color='g', label='IL15', alpha=0.7, ax=ax)
        # plt.ylim(0,(y_max + (0.2 * y_max)))
        plt.ylabel("Surface IL2Rb (% x " + str(self.y_max) + ')')
        plt.xlabel("Time (min)")
        plt.legend()
        plt.show()


    output = surf2Rb.calc(unkVec) * y_max
    IL2_1_plus = output[0:7]
    IL2_500_plus = output[7:14]
    IL2_1_minus = output[14:21]
    IL2_500_minus = output[21:28]
    IL15_1_plus = output[28:35]
    IL15_500_plus = output[35:42]
    IL15_1_minus = output[42:49]
    IL15_500_minus = output[49:56]

    plot_structure(IL2_1_minus, IL15_1_minus, '1 nM and IL2Ra-', ax[0])
    plot_structure(IL2_500_minus, IL15_500_minus, "500 nM and IL2Ra-", ax[1])
    plot_structure(IL2_1_plus, IL15_1_plus, "1 nM and IL2Ra+", ax[2])
    plot_structure(IL2_500_plus, IL15_500_plus, "500 nM and IL2Ra+", ax[3])
