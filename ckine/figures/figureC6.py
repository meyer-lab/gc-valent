"""
This creates Figure 6 for IL2Ra correlatoin data analysis.
"""

import os
import numpy as np
import pandas as pd
from .figureCommon import subplotLabel, getSetup

path_here = os.path.dirname(os.path.dirname(__file__))


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((10, 5), (2, 4))
    subplotLabel(ax)

    path_here = os.path.dirname(os.path.dirname(__file__))

    # Imports receptor levels from .csv created by figC5
    receptor_levels = pd.read_csv(path_here + "/data/receptor_levels.csv", header=0)

    cell_types = ['T-reg', 'T-helper', 'NK', 'CD8+']

    for index, cell_type in enumerate(cell_types):
        i = 2 * index

        alphaLevels = receptor_levels.loc[(receptor_levels['Cell Type'] == cell_type) & (receptor_levels['Receptor'] == 'CD25')]
        betaLevels = receptor_levels.loc[(receptor_levels['Cell Type'] == cell_type) & (receptor_levels['Receptor'] == 'CD122')]
        gammaLevels = receptor_levels.loc[(receptor_levels['Cell Type'] == cell_type) & (receptor_levels['Receptor'] == 'CD132')]

        alphaCounts = alphaLevels['Count'].reset_index(drop=True)
        betaCounts = betaLevels['Count'].reset_index(drop=True)
        d = {'alpha': alphaCounts, 'beta': betaCounts}
        recepCounts = pd.DataFrame(data=d)
        recepCounts = recepCounts.dropna()
        recepCounts = recepCounts[(recepCounts[['alpha', 'beta']] != 0).all(axis=1)]

        hex1 = ax[i]
        hex1.hexbin(recepCounts['alpha'], recepCounts['beta'], xscale='log', yscale='log', mincnt=1, cmap='viridis')
        hex1.set_xlabel('CD25')
        hex1.set_ylabel('CD122')
        hex1.set_title(cell_type + ' Alpha-Beta correlation')

        alphaCounts = alphaLevels['Count'].reset_index(drop=True)
        gammaCounts = gammaLevels['Count'].reset_index(drop=True)
        d2 = {'alpha': alphaCounts, 'gamma': gammaCounts}
        recepCounts2 = pd.DataFrame(data=d2)
        recepCounts2 = recepCounts2.dropna()
        recepCounts2 = recepCounts2[(recepCounts2[['alpha', 'gamma']] != 0).all(axis=1)]

        hex2 = ax[i + 1]
        hex2.hexbin(recepCounts2['alpha'], recepCounts2['gamma'], xscale='log', yscale='log', mincnt=1, cmap='viridis')
        hex2.set_xlabel('CD25')
        hex2.set_ylabel('CD132')
        hex2.set_title(cell_type + ' Alpha-Gamma correlation')

    return f
