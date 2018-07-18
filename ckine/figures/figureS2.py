"""
This creates Figure S2.
"""
import string
from .figureCommon import subplotLabel, getSetup
from ..tensor_generation import prepare_tensor
from ..Tensor_analysis import perform_decomposition
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    x, y = 3, 5
    ax, f = getSetup((10, 10), (x, y))

    values, _, mat, _, cell_names = prepare_tensor(2)
    factors = perform_decomposition(values, 6)

    for row in range(x):
        plot_timepoint(ax[row*y], factors[0], row+1, row+2)
        plot_cells(ax[row*y + 1], factors[1], mat, row+1, row+2, cell_names)
        plot_ligands(ax[row*y + 2], factors[2], row+1, row+2)
        plot_values(ax[row*y + 3], factors[3], row+1, row+2)

    f.tight_layout()

    return f


def plot_values(ax, factors, component_x, component_y):
    """Plot the values decomposition factors matrix."""
    #Generate a plot for component x vs component y of the factors[3] above representing our values
    labels = ['IL2', 'IL15', 'IL7', 'IL9', 'IL4','IL21','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra']
    #Set Active to color red. Set Surface to color blue. Set Total to color black

    ax.scatter(factors[0:6, component_x - 1], factors[0:6, component_y - 1], color = 'r', label = 'Ligand Activity')
    ax.scatter(factors[6:14, component_x - 1], factors[6:14, component_y - 1], color = 'b', label = 'Surface Receptor')
    ax.scatter(factors[14::, component_x - 1], factors[14::, component_y - 1], color = 'k', label = 'Total Receptor')

    ax.annotate(labels, xy=(factors[:, component_x - 1], factors[:, component_y - 1]), xytext = (0, 0), textcoords = 'offset points')

    set_axis_limits(ax, factors, component_x, component_y)

def set_axis_limits(ax, factor, component_x, component_y):
    """Function to set axes labels and limits of the designated factor matrix."""
    ax.set_xlabel('Component ' + str(component_x))
    ax.set_ylabel('Component ' + str(component_y))
    ax.set_xlim(-np.max(np.absolute(factor[:,component_x - 1]))*1.1, np.max(np.absolute(factor[:,component_x - 1]))*1.1)
    ax.set_ylim(-np.max(np.absolute(factor[:,component_y - 1]))*1.1, np.max(np.absolute(factor[:,component_y - 1]))*1.1)
    
def plot_timepoint(ax, factors, component_x, component_y):
    """Plot the timepoint decomposition in the first column of figS2."""
    ax.scatter(factors[:, component_x - 1], factors[:, component_y - 1], color = 'k')

    ax.annotate(str(1000), xy=(factors[999, component_x - 1], factors[999, component_y - 1]), xytext = (0, 0), textcoords = 'offset points')

    set_axis_limits(ax, factors, component_x, component_y)

def plot_cells(ax,factors, mat, component_x, component_y, cell_names):
    """This function plots the combination decomposition based on cell type."""
    colors = cm.rainbow(np.linspace(0, 1, len(cell_names)))
    ax.scatter(factors[:,component_x - 1], factors[:,component_y - 1], c=colors, label = cell_names)
    set_axis_limits(ax, factors, component_x, component_y)

def plot_ligands(ax, factors, component_x, component_y):
    "This function is to plot the ligand combination dimension of the values tensor."
    ax.scatter(factors[:,component_x - 1], factors[:,component_y - 1])
    set_axis_limits(ax, factors, component_x, component_y)
