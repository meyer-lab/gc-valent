"""
This creates Figure 3.
"""
import string
import tensorly as tl
import seaborn as sns; sns.set()
from .figureCommon import subplotLabel, getSetup, plot_cells, plot_ligands, plot_timepoints, plot_R2X
from ..imports import import_Rexpr
from ..tensor import perform_decomposition
from ..make_tensor import make_tensor, n_lig

cell_dim = 1  # For this figure, the cell dimension is along the second [python index 1].
values, _, mat, _, _ = make_tensor()
values = tl.tensor(values)


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    x, y = 3, 4
    ax, f = getSetup((7.5, 7), (x, y), mults=[2, 4], multz={2: 2, 4: 3})
    real_mults = [0, 4]

    # Blank out for the cartoon
    ax[3].axis('off')

    n_ligands = n_lig(mut=False)
    data, _, cell_names = import_Rexpr()
    factors_activity = []
    for jj in range(len(mat) - 1):
        factors = perform_decomposition(values, jj + 1, cell_dim)
        factors_activity.append(factors)

    n_comps = 3
    factors_activ = factors_activity[n_comps - 1]

    catplot_receptors(ax[2], data)
    plot_R2X(ax[4], values, factors_activity, n_comps=5, cells_dim=cell_dim)

    # Add subplot labels
    for ii, item in enumerate(ax):
        h = 2.5 if ii in real_mults else 1
        subplotLabel(item, string.ascii_uppercase[ii], hstretch=h)  # Add subplot labels

    #plot_timepoints(ax[3], factors_activ[0])  # Change final input value depending on need

    #plot_cells(ax[4], factors_activ[1], 1, 2, cell_names, ax_pos=4)
    #plot_cells(ax[8], factors_activ[1], 2, 3, cell_names, ax_pos=8)

    #plot_ligands(ax[5], factors_activ[2], 1, 2, ax_pos=5, n_ligands=n_ligands, mesh=mat, fig=f)
    #plot_ligands(ax[9], factors_activ[2], 2, 3, ax_pos=9, n_ligands=n_ligands, mesh=mat, fig=f)

    f.tight_layout()

    return f


def catplot_receptors(ax, data):
    """Plot Bar graph for Receptor Expression Data. """
    sns.set(style="ticks")
    sns.catplot(x="Cell Type", y="Count", hue="Receptor", data=data, ci=68, ax=ax, legend=False, legend_out=False)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.5, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), borderpad=0.4, labelspacing=0.2, handlelength=0.2, handletextpad=0.5, markerscale=0.7, fontsize=8)
    ax.set_ylabel("Surface Receptor [# / cell]")
    ax.set_xticklabels(ax.get_xticklabels(),
                       rotation=25, rotation_mode="anchor", ha="right",
                       position=(0, 0.02), fontsize=7.5)
