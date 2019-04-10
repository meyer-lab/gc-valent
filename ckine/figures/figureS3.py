"""
This creates Figure S3, which covers the Tucker factorization form.
"""
import string
import tensorly as tl
from tensorly import unfold
from .figureCommon import subplotLabel, getSetup, plot_cells, plot_ligands, plot_timepoints, values, set_bounds
from ..Tensor_analysis import perform_tucker, find_R2X_tucker
from ..tensor_generation import cell_names
import seaborn  as sns

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    x, y = 3, 3
    ax, f = getSetup((10, 8), (x, y), empts=[3])

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])  # Add subplot labels

    rank_list = [3, 4, 4]
    out = perform_tucker(values, rank_list)
    print(out[0])
    print(find_R2X_tucker(values, out))

    plot_core(ax[6:9], out[0])
    factors = out[1]
    plot_timepoints(ax[0], tl.to_numpy(factors[0]))

    for row in range(2):
        compNum = 2 * row + 1

        plot_cells(ax[row * y + 1], tl.to_numpy(factors[1]), compNum, compNum + 1, cell_names, ax_pos=row * y + 1, fig3=False)
        plot_ligands(ax[row * y + 2], tl.to_numpy(factors[2]), compNum, compNum + 1, ax_pos= row * y + 2, fig3=False)

    f.tight_layout()

    return f

def plot_core(ax, core):
    #Begin by unfolding the core tensor on its 3 faces.
    X1, X2, X3 = unfold(core, 0), unfold(core, 1), unfold(core, 2)
    sns.heatmap(X1,cmap="YlGnBu",cbar=True,ax=ax[0])
    ax[0].set_title("Time Dimension")
    sns.heatmap(X2,cmap="YlGnBu",cbar=True,ax=ax[1])
    ax[1].set_title("Cell Dimension")
    sns.heatmap(X3,cmap="YlGnBu",cbar=True,ax=ax[2])
    ax[2].set_title("Ligand Dimension")
