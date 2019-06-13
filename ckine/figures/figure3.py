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
    real_mults = [2, 4]

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

    PCA_receptor[ax[0], ax[1], cell_names, data]
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

def PCA_receptor(ax1, ax2, cell_names, data):
    """Plot PCA scores and loadings for Receptor Expression Data. """
    pca = PCA(n_components = 2)
    data = stats.zscore(data.astype(float), axis = 0)
    scores = pca.fit(data.T).transform(data.T) #34 cells by n_comp
    loadings = pca.components_ #n_comp by 8 receptors
    expVar = pca.explained_variance_ratio_

    colors = cm.rainbow(np.linspace(0, 1, len(cell_names)))
    markersCells = ['^', '*', 'D', 's', 'X', 'o', '4', 'H', 'P', '*', 'D', 's', 'X']  # 'o', 'd', '1', '2', '3', '4', 'h', 'H', 'X', 'v', '*', '+', '8', 'P', 'p', 'D', '_','D', 's', 'X', 'o'
    markersReceptors = ['^', '4', 'P', '*', 'D'] #'s', 'X' ,'o'
    labelReceptors = ['IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra'] #'IL9R', 'IL4Ra', 'IL21Ra']

    for ii in range(scores.shape[0]):
        ax1.scatter(scores[ii,0], scores[ii,1], c = colors[ii], marker = markersCells[ii], label = cell_names[ii])

    for jj in range(loadings.shape[1]):
        ax2.scatter(loadings[0,jj], loadings[1,jj], marker = markersReceptors[jj], label = labelReceptors[jj])

    x_max1 = np.max(np.absolute(np.asarray(ax1.get_xlim())))*1.1
    y_max1 = np.max(np.absolute(np.asarray(ax1.get_ylim())))*1.1

    x_max2 = np.max(np.absolute(np.asarray(ax2.get_xlim())))*1.1
    y_max2 = np.max(np.absolute(np.asarray(ax2.get_ylim())))*1.1

    ax1.set_xlim(-x_max1, x_max1)
    ax1.set_ylim(-y_max1, y_max1)
    ax1.set_xlabel('PC1 (' + str(round(expVar[0]*100, 2))+ '%)')
    ax1.set_ylabel('PC2 (' + str(round(expVar[1]*100, 2))+ '%)')
    ax1.set_title('Scores')
    ax1.legend(loc='upper left', bbox_to_anchor=(3.5, 1.735))

    ax2.set_xlim(-x_max2, x_max2)
    ax2.set_ylim(-y_max2, y_max2)
    ax2.set_xlabel('PC1 (' + str(round(expVar[0]*100, 2))+ '%)')
    ax2.set_ylabel('PC2 (' + str(round(expVar[1]*100, 2))+ '%)')
    ax2.set_title('Loadings')
    ax2.legend()
