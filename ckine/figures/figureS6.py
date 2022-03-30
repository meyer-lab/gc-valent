"""
This creates Figure S6, plotting Treg to off target signaling for vaying mutein dose for different IL2 formats
"""
from os.path import dirname
from .figureCommon import getSetup, Wass_KL_Dist, CITE_RIDGE, CITE_SVM


path_here = dirname(dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((6, 6), (2, 2))
    cellTarget = "Treg"
    Wass_KL_Dist(ax[0:2], cellTarget, 10, RNA=True)
    CITE_RIDGE(ax[2], cellTarget, RNA=True)
    CITE_SVM(ax[3], cellTarget, sampleFrac=0.05, RNA=True)

    return f
