"""
This creates Figure 3.
"""
import string
import tensorly as tl
import numpy as np
import seaborn as sns
from .figureCommon import subplotLabel, getSetup, plot_cells, plot_ligands, plot_timepoints
from ..imports import import_Rexpr
from ..tensor import find_R2X, perform_decomposition
from ..make_tensor import make_tensor, n_lig

n_ligands = n_lig
cell_dim = 1  # For this figure, the cell dimension is along the second [python index 1].
values, _, mat, _, _ = make_tensor()
values = tl.tensor(values)

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    x, y = 3, 4
    ax, f = getSetup((12, 9), (x, y), mults=[0, 2], multz={0: 2, 2: 2}, empts=[7, 11])
    # Blank out for the cartoon
    ax[0].axis('off')

    data, _, cell_names = import_Rexpr()
    factors_activity = []
    for jj in range(len(mat) - 1):
        factors = perform_decomposition(values, jj + 1, cell_dim)
        factors_activity.append(factors)

    n_comps = 4
    factors_activ = factors_activity[n_comps - 1]

    bar_receptors(ax[1], data)
    plot_R2X(ax[2], values, factors_activity, n_comps=5, cells_dim=cell_dim)

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])  # Add subplot labels

    plot_timepoints(ax[5], factors_activ[0])  # Change final input value depending on need

    plot_cells(ax[3], factors_activ[1], 1, 2, cell_names, ax_pos=3)
    plot_cells(ax[6], factors_activ[1], 3, 4, cell_names, ax_pos=6)

    plot_ligands(ax[4], factors_activ[2], 1, 2, ax_pos=4, n_ligands=n_ligands, mesh=mat)
    plot_ligands(ax[7], factors_activ[2], 3, 4, ax_pos=7, n_ligands=n_ligands, mesh=mat)

    f.tight_layout()

    return f

def bar_receptors(ax, data):
    """Plot Bar graph for Receptor Expression Data. """
    sns.barplot(x="Cell Type", y="Count", hue="Receptor", data=data, ci=68, ax=ax, capsize=.2, errwidth=0.4)
    ax.legend(loc='best')
    ax.set_ylabel("Surface Receptor [# / cell]")
    ax.set_xticklabels(ax.get_xticklabels(),
                       rotation=40, rotation_mode="anchor", ha="right",
                       position=(0, 0.05), fontsize=6.5)


def plot_R2X(ax, tensor, factors_list, n_comps, cells_dim):
    """Function to plot R2X bar graph."""
    R2X_array = list()
    for n in range(n_comps):
        factors = factors_list[n]
        R2X = find_R2X(tensor, factors, cells_dim)
        R2X_array.append(R2X)
    ax.plot(range(1, n_comps + 1), R2X_array, 'ko', label='Overall R2X')
    ax.set_ylabel('R2X')
    ax.set_xlabel('Number of Components')
    ax.set_ylim(0, 1)
    ax.set_xticks(np.arange(1, n_comps + 1))
    ax.set_xticklabels(np.arange(1, n_comps + 1))
