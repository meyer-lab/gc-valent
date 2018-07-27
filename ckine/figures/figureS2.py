"""
This creates Figure S2.
"""
import string
import numpy as np
import itertools
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from .figureCommon import subplotLabel, getSetup
from ..tensor_generation import prepare_tensor
from ..Tensor_analysis import perform_decomposition, combo_low_high

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects

    x, y = 6, 4
    ssize = 3
    ax, f = getSetup((ssize*y, ssize*x), (x, y))
    lig = 5
    values, _, mat, _, cell_names = prepare_tensor(lig)
    factors = perform_decomposition(values, 9)

    for row in range(4):
        compNum = 2*row + 1
        plot_timepoint(ax[row*y], factors[0], compNum, compNum+1)
        plot_cells(ax[row*y + 1], factors[1], compNum, compNum+1, cell_names)
        plot_ligands(ax[row*y + 2], factors[2], mat, lig, compNum, compNum+1)
        plot_values(ax[row*y + 3], factors[3], compNum, compNum+1)

        # Set axes to center on the origin, and add labels
        for col in range(y):
            ax[row*y + col].set_xlabel('Component ' + str(compNum))
            ax[row*y + col].set_ylabel('Component ' + str(compNum+1))

            x_max = np.max(np.absolute(np.asarray(ax[row*y + col].get_xlim())))*1.1
            y_max = np.max(np.absolute(np.asarray(ax[row*y + col].get_ylim())))*1.1

            ax[row*y + col].set_xlim(-x_max, x_max)
            ax[row*y + col].set_ylim(-y_max, y_max)

    f.tight_layout()

    return f


def plot_values(ax, factors, component_x, component_y):
    """Plot the values decomposition factors matrix."""
    #Generate a plot for component x vs component y of the factors[3] above representing our values
    # The markers are for the following elements in order: 'IL2', 'IL15', 'IL7', 'IL9', 'IL4','IL21','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra.'

    markersLigand = itertools.cycle(('^', '*'))
    markersReceptors = itertools.cycle(('^', '4', 'P', '*'))
    ligand_activity = itertools.cycle(('IL2 Activity', 'IL15 Activity'))
    surface_amount = itertools.cycle(('Surface IL2Ra', 'Surface IL2Rb', 'Surface gc', 'Surface IL15Ra'))
    total_amount = itertools.cycle(('Total IL2Ra', 'Total IL2Rb', 'Total gc', 'Total IL15Ra'))

    #[0,1,6,7,8,9,14,15,16,17]
    #Set Active to color red. Set Surface to color blue. Set Total to color black
    for q,p in zip(factors[0:2, component_x - 1], factors[0:2, component_y - 1]):
        ax.plot(q, p, linestyle = '', c = 'r', marker = next(markersLigand), label = next(ligand_activity))

    for q,p in zip(factors[2:6, component_x - 1], factors[2:6, component_y - 1]):
        ax.plot(q, p, linestyle = '', c = 'b', marker = next(markersReceptors), label = next(surface_amount))
    
    for q,p in zip(factors[6::, component_x - 1], factors[6::, component_y - 1]):
        ax.plot(q, p, linestyle = '', c = 'k', marker = next(markersReceptors), label = next(total_amount))

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


def plot_timepoint(ax, factors, component_x, component_y):
    """Plot the timepoint decomposition in the first column of figS2."""
    ax.plot(factors[:, component_x - 1], factors[:, component_y - 1], color = 'k')
    ax.scatter(factors[-1, component_x - 1], factors[-1, component_y - 1], s = 12, color = 'b')

def plot_cells(ax, factors, component_x, component_y, cell_names):
    """This function plots the combination decomposition based on cell type."""
    colors = cm.rainbow(np.linspace(0, 1, len(cell_names)))
    
    markersCells = ['^', '*', 'D', 's', 'X', 'o', '^', '4', 'P', '*', 'D', 's', 'X' ,'o']
    
    
    for ii in range(len(factors[:, component_x - 1])):
        ax.scatter(factors[ii, component_x - 1], factors[ii, component_y - 1], c = colors[ii], marker = markersCells[ii], label = cell_names[ii])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

def plot_ligands(ax, factors, mat, lig, component_x, component_y):
    "This function is to plot the ligand combination dimension of the values tensor."
    IL2_low_high, IL15_low_high = combo_low_high(mat,lig)
    cm = plt.cm.get_cmap('YlOrRd')

    IL15_markers = ['*', 'D', 's', 'X' ,'o']
    labels = ['0.001 nm IL15', '0.018 nm IL15', '0.32 nm IL15', '5.6 nm IL15', '100 nm IL15']
    
    for jj in range(lig):
        ax.scatter(factors[IL15_low_high[jj],component_x - 1], factors[IL15_low_high[jj],component_y - 1], marker = IL15_markers[jj], c = [0,1,2,3,4], cmap = cm, label = labels[jj])
        
        
    leg = ax.legend()
    for k in range(lig):
        leg.legendHandles[k].set_color('black')
