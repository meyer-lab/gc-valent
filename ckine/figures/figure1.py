"""
This creates Figure 1.
"""
from .figureCommon import subplotLabel, getSetup
from ..plot_model_prediction import surf_IL2Rb, pstat, surf_gc


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    subplotLabel(ax[0], 'A')
    
    IL2Rb_surf_perc(ax[1:5])

    f.tight_layout()

    return f

def IL2Rb_surf_perc(ax):
    surf2Rb = surf_IL2Rb()
    unkVec = np.array([0., 0., 0., 0., 2.02E-05, 0.38406199, 0.047815524, 0.118719663, 0.088470489, 0.092582828, 0.107673095, 0.101291622, 0.089522673, 0.023247533, 1.954609999, 0.042469514, 0.300662763, 0.004387724, 2.166630382, 0.673719144, 2.077778334, 0.288524177, 0., 0.])
    ax = surf2Rb.plot(unkVec)
