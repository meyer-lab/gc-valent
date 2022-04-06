"""
This creates Figure S6, RNA seq based epitope ID without binding model.
"""
import pandas as pd
from os.path import dirname
from os.path import join
from .figureCommon import getSetup, Wass_KL_Dist, CITE_RIDGE, CITE_SVM

path_here = dirname(dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((6, 6), (2, 2))
    cellTarget = "Treg"
    distDF = Wass_KL_Dist(ax[0:2], cellTarget, 10, RNA=True)
    RIDGEdf = CITE_RIDGE(ax[2], cellTarget, RNA=True)
    SVMdf = CITE_SVM(ax[3], cellTarget, sampleFrac=0.05, RNA=True)

    selectiveDF = pd.DataFrame()
    selectiveDF = pd.concat([selectiveDF, pd.DataFrame({"Method": "Distance", "Marker": distDF.Marker.values})])
    selectiveDF = pd.concat([selectiveDF, pd.DataFrame({"Method": "RIDGE", "Marker": RIDGEdf.Marker.values})])
    selectiveDF = pd.concat([selectiveDF, pd.DataFrame({"Method": "SVM", "Marker": SVMdf})])
    RNAuniqueDF = pd.DataFrame({"Gene": selectiveDF.Marker.unique()})

    RNAuniqueDF.to_csv(join(path_here, "data/RNAseq_TregUnique.csv"))

    return f
