"""
This creates Figure 6.
"""
import numpy as np
import seaborn as sns
from .figureCommon import subplotLabel, getSetup, plot_pstat
from ..imports import import_muteins

dataMean, _ = import_muteins()
dataMean.reset_index(inplace=True)
dataMean['Concentration'] = np.log10(dataMean['Concentration'])

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 5), empts=[18, 19])
    
    axis = 0
    
    first_group_ligands = ['IL2-060', 'IL2-062']
    
    for i, ligand_name in enumerate(first_group_ligands):
        for j, cell_name in enumerate(dataMean.Cells.unique()):
            axis = i*9+j
            if axis == 17:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name) & (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, legend='full')
                ax[axis].legend(loc='lower right')
            else:
                sns.scatterplot(x="Concentration", y="RFU", hue="Time", data=dataMean.loc[(dataMean["Cells"] == cell_name) & (dataMean["Ligand"] == ligand_name)], ax=ax[axis], s=10, legend=False)
            ax[axis].set(xlabel=("[" + ligand_name + "] (log$_{10}$[nM])"), ylabel="Activity", title=cell_name)
            
    return f
