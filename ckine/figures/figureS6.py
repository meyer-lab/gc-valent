"""
This creates Figure S6, RNA seq based epitope ID without binding model (makes DF)
"""
import pandas as pd
import matplotlib.pyplot as plt
from os.path import dirname
from os.path import join
from .figureCommon import getSetup, CITE_RIDGE

path_here = dirname(dirname(__file__))
plt.rcParams['svg.fonttype'] = 'none'

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    ax, f = getSetup((6, 3), (1, 1))
    cellTarget = "Treg"
    RIDGEdf = CITE_RIDGE(ax[0], cellTarget, RNA=True)

    selectiveDF = pd.DataFrame()
    selectiveDF = pd.concat([selectiveDF, pd.DataFrame({"Method": "RIDGE", "Marker": RIDGEdf.Marker.values})])
    RNAuniqueDF = pd.DataFrame({"Gene": selectiveDF.Marker.unique()})

    RNAuniqueDF.to_csv(join(path_here, "data/RNAseq_TregUnique.csv"))

    return f
